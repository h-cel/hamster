#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MAIN FUNCTION TO VALIDATE
"""

###########################################################################
##--- MODULES
###########################################################################

import argparse
import calendar
import csv
import datetime
import fnmatch
import gzip
import imp
import math
import os
import random
import re
import struct
import sys
import time
import timeit
import warnings
from datetime import date, datetime, timedelta
from functools import reduce
from math import acos, atan, atan2, cos, floor, sin, sqrt

import h5py
import netCDF4 as nc4
import numpy as np
import pandas as pd
from dateutil.relativedelta import relativedelta

###########################################################################
##--- FUNCTIONS + COMMAND LINE ARGUMENTS
###########################################################################

## (1) LOADING FUNCTIONS
exec(open("hamsterfunctions.py").read())

## (2) COMMAND LINE ARGUMENTS
# read command line arguments (dates, thresholds and other flags)
args = read_cmdargs()
verbose = args.verbose
print(printsettings(args))

# just waiting a random number of seconds (max. 30s)
# to avoid overlap of path.exist and makedirs between parallel jobs (any better solution?)
if args.waiter:
    waiter = random.randint(0, 30)
    time.sleep(waiter)

###########################################################################
##--- PATHS
###########################################################################

## determine working directory
wpath = os.getcwd()
os.chdir(wpath)

## load input and output paths & input file name base(s)
print("Using paths from: " + wpath + "/" + args.pathfile)
content = imp.load_source("", wpath + "/" + args.pathfile)  # load like a python module
path_refp = content.path_ref_p
path_refe = content.path_ref_e
path_refh = content.path_ref_h
path_diag = content.path_diag

###########################################################################
##--- ADDITIONAL VALIDATION FUNCTIONS
###########################################################################


def contingency_table(ref, mod, thresh=0):
    # creates a contingency table based on 1D np.arrays
    # ATTN: mod is already binary data, i.e. 1 = detected; 0 = not detected
    ieventobs = np.where(ref > thresh)[0]
    ineventobs = np.where(ref <= thresh)[0]
    a = len(np.where(mod[ieventobs] >= 1)[0])  # hits
    b = len(np.where(mod[ineventobs] >= 1)[0])  # false alarms
    c = len(np.where(mod[ieventobs] < 1)[0])  # misses
    d = len(np.where(mod[ineventobs] < 1)[0])  # correct negatives
    return {"a": a, "b": b, "c": c, "d": d}


def try_div(x, y):
    try:
        return x / y
    except ZeroDivisionError:
        return 0


def calc_ctab_measures(cdict):
    # calculates common contingency table scores
    # scores following definitions from https://www.cawcr.gov.au/projects/verification/
    a = cdict["a"]  # hits
    b = cdict["b"]  # false alarms
    c = cdict["c"]  # misses
    d = cdict["d"]  # correct negatives
    # calculate scores
    acc = try_div(a + d, a + b + c + d)  # accuracy
    far = try_div(b, a + b)  # false alarm ratio
    fbias = try_div(a + b, a + c)  # frequency bias
    pod = try_div(a, a + c)  # probability of detection (hit rate)
    pofd = try_div(b, b + d)  # probability of false detection (false alarm rate)
    sr = try_div(a, a + b)  # success ratio
    ts = try_div(a, a + c + b)  # threat score (critical success index)
    a_random = try_div((a + c) * (a + b), a + b + c + d)
    ets = try_div(
        (a - a_random), (a + b + c + a_random)
    )  # equitable threat score (gilbert skill score)
    pss = pod - pofd  # peirce's skill score (true skill statistic)
    odr = try_div(a * d, c * b)  # odd's ratio
    return {
        "acc": acc,
        "far": far,
        "fbias": fbias,
        "pod": pod,
        "pofd": pofd,
        "sr": sr,
        "pss": pss,
        "odr": odr,
    }


def init_netcdf(ofile, lat, lon):
    print(" Writing output file: " + str(ofile))
    ncf = nc4.Dataset(ofile, "w", format="NETCDF4")
    ncf.createDimension("lat", lat.size)
    ncf.createDimension("lon", lon.size)
    latitudes = ncf.createVariable("lat", "f8", "lat")
    longitudes = ncf.createVariable("lon", "f8", "lon")
    ncf.title = "Validation statistics"
    ncf.description = ""
    today = datetime.datetime.now()
    ncf.history = "Created " + today.strftime("%d/%m/%Y %H:%M:%S")
    ncf.institution = (
        "Hydro-Climate Extremes Laboratory (H-CEL), Ghent University, Ghent, Belgium"
    )
    ncf.source = "Validation statistics for HAMSTER"
    latitudes.units = "degrees_north"
    longitudes.units = "degrees_east"
    latitudes[:] = lat[:]
    longitudes[:] = lon[:]
    ncf.close()


def calc_stats(mdata, mndata, rdata, thresh=0.001):
    bias = np.zeros(shape=(mdata.shape[1], mdata.shape[2]))
    acc = np.zeros(shape=(mdata.shape[1], mdata.shape[2]))
    pod = np.zeros(shape=(mdata.shape[1], mdata.shape[2]))
    pofd = np.zeros(shape=(mdata.shape[1], mdata.shape[2]))
    pss = np.zeros(shape=(mdata.shape[1], mdata.shape[2]))
    fbias = np.zeros(shape=(mdata.shape[1], mdata.shape[2]))
    odr = np.zeros(shape=(mdata.shape[1], mdata.shape[2]))
    print(" Processing...")
    for y in range(rdata.shape[1]):
        progress = 100 * y / (rdata.shape[1] - 1)
        if verbose and round(progress) in range(10, 100, 5):
            print(" ..." + str(round(progress)) + "%...")
        for x in range(rdata.shape[2]):
            myctab = calc_ctab_measures(
                contingency_table(rdata[:, y, x], mndata[:, y, x], thresh=thresh)
            )
            acc[y, x] = myctab["acc"]
            pod[y, x] = myctab["pod"]
            pofd[y, x] = myctab["pofd"]
            pss[y, x] = myctab["pss"]
            fbias[y, x] = myctab["fbias"]
            odr[y, x] = myctab["odr"]
            diff = mdata[:, y, x] - rdata[:, y, x]
            bias[y, x] = np.nanmean(diff)
    return {
        "acc": acc,
        "fbias": fbias,
        "pod": pod,
        "pofd": pofd,
        "pss": pss,
        "odr": odr,
        "diff": diff,
        "bias": bias,
    }


def write_to_netcdf(ofile, vals, var="P"):
    print(
        " Writing " + str(var) + " validation statistics to output file: " + str(ofile)
    )

    # create netCDF4 instance
    ncf = nc4.Dataset(ofile, "r+", format="NETCDF4")

    # create variables
    ncacc = ncf.createVariable(str(var) + "_acc", "f8", ("lat", "lon"))
    ncpod = ncf.createVariable(str(var) + "_pod", "f8", ("lat", "lon"))
    ncpofd = ncf.createVariable(str(var) + "_pofd", "f8", ("lat", "lon"))
    ncpss = ncf.createVariable(str(var) + "_pss", "f8", ("lat", "lon"))
    ncodr = ncf.createVariable(str(var) + "_odr", "f8", ("lat", "lon"))
    ncfbias = ncf.createVariable(str(var) + "_fbias", "f8", ("lat", "lon"))
    ncbias = ncf.createVariable(str(var) + "_bias", "f8", ("lat", "lon"))

    # set attributes
    ncpod.long_name = str(var) + ": probability of detection (hit rate)"
    ncacc.long_name = str(var) + ": accuracy"
    ncbias.long_name = str(var) + ": bias"
    ncfbias.long_name = str(var) + ": frequency bias"
    ncpofd.long_name = str(var) + ": probability of false detection (false alarm rate)"
    ncpss.long_name = str(var) + ": peirce´s skill score (true skill statistic)"
    ncodr.long_name = str(var) + ": odd`s ratio"

    # write to netcdf
    ncodr[:] = vals["odr"][:]
    ncpod[:] = vals["pod"][:]
    ncpss[:] = vals["pss"][:]
    ncpofd[:] = vals["pofd"][:]
    ncacc[:] = vals["acc"][:]
    ncpofd[:] = vals["pofd"][:]
    ncfbias[:] = vals["fbias"][:]
    ncbias[:] = vals["bias"][:]

    # close file
    ncf.close()


###########################################################################
###########################################################################
##--- MAIN
###########################################################################
###########################################################################


def main_validation(
    ryyyy,
    ayyyy,
    am,
    opath_diag,  # diagnosis (output)
    ipath_refp,
    ipath_refe,
    ipath_refh,
    opath,
    ofile_base,  # output
    fprec,
    fevap,
    fheat,
    verbose,
    veryverbose,
    fwrite_netcdf,
):

    ## construct precise input and storage paths
    ofilename = (
        str(ofile_base)
        + "_diag_r"
        + str(ryyyy)[-2:]
        + "_"
        + str(ayyyy)
        + "-"
        + str(am).zfill(2)
        + "_validation.nc"
    )
    ofile = opath + "/" + ofilename

    #### DISCLAIMER
    if verbose:
        disclaimer()
        print("\n PROCESSING: \t", ayyyy, "-", str(am).zfill(2) + "\n")
    ## Resets & consistency checks
    if verbose:
        print(" ! using input paths: \t")
        print("\t" + str(opath_diag))
        print(" ! using reference data from: \t")
        print("\t" + str(ipath_refp))
        print("\t" + str(ipath_refe))
        print("\t" + str(ipath_refh))
        print(" ! writing netcdf output: \t")
        print("\t" + str(ofile))
        print(
            "\n============================================================================================================"
        )
        print("\n")

    ##--1. load diagnosis data ####################################################
    if verbose:
        print(" * Reading diagnosis data...")

    # read concatenated data
    ifilename = (
        str(opath_diag)
        + "/"
        + str(ofile_base)
        + "_diag_r"
        + str(ryyyy)[-2:]
        + "_"
        + str(ayyyy)
        + "-"
        + str(am).zfill(2)
        + ".nc"
    )
    with nc4.Dataset(ifilename, mode="r") as f:
        idate_seq = nc4.num2date(f["time"][:], f["time"].units, f["time"].calendar)
    totlats, totlons = read_diagdata(
        opath_diag, ofile_base, ryyyy, idate_seq, var="grid"
    )
    glon, glat, garea = makegrid(resolution=abs(totlats[0] - totlats[1]))
    ftime = read_diagdata(opath_diag, ofile_base, ryyyy, idate_seq, var="time")
    fdays = np.unique(cal2date(ftime))
    fyears = np.unique(date2year(ftime))
    if fevap:
        E = read_diagdata(opath_diag, ofile_base, ryyyy, idate_seq, var="E")
        E_npart = read_diagdata(
            opath_diag, ofile_base, ryyyy, idate_seq, var="E_n_part"
        )
    if fprec:
        P = -read_diagdata(opath_diag, ofile_base, ryyyy, idate_seq, var="P")
        P_npart = read_diagdata(
            opath_diag, ofile_base, ryyyy, idate_seq, var="P_n_part"
        )
    if fheat:
        H = read_diagdata(opath_diag, ofile_base, ryyyy, idate_seq, var="H")
        H_npart = read_diagdata(
            opath_diag, ofile_base, ryyyy, idate_seq, var="H_n_part"
        )

    # make sure we use daily aggregates
    if fdays.size != ftime.size:
        if fevap:
            Etot = convert2daily(E, ftime, fagg="sum")
            Enparttot = convert2daily(E_npart, ftime, fagg="sum")
        if fprec:
            Ptot = convert2daily(P, ftime, fagg="sum")
            Pnparttot = convert2daily(P_npart, ftime, fagg="sum")
        if fheat:
            Htot = convert2daily(H, ftime, fagg="mean")
            Hnparttot = convert2daily(H_npart, ftime, fagg="sum")
    else:
        if fevap:
            Etot = E
            Enparttot = E_npart
        if fprec:
            Ptot = P
            Pnparttot = P_npart
        if fheat:
            Htot = H
            Hnparttot = H_npart

    ##--2. load reference data ####################################################
    """
    this part is STRICTLY CODED FOR (12-hourly) ERA-INTERIM only (so far),
    and HARDCODED too
    """
    if verbose:
        print(" * Reading reference data...")

    if fevap:
        Eref, reflats, reflons = eraloader_12hourly(
            var="e",
            datapath=ipath_refe + "/E_1deg_",
            maskpos=True,
            maskneg=False,
            uptake_years=fyears,
            uptake_dates=fdays,
            lats=totlats,
            lons=totlons,
        )
        gridcheck(totlats, reflats, totlons, reflons)

    if fheat:
        Href, reflats, reflons = eraloader_12hourly(
            var="sshf",
            datapath=ipath_refh + "/H_1deg_",
            maskpos=True,
            maskneg=False,
            uptake_years=fyears,
            uptake_dates=fdays,
            lats=totlats,
            lons=totlons,
        )
        gridcheck(totlats, reflats, totlons, reflons)

    if fprec:
        Pref, reflats, reflons = eraloader_12hourly(
            var="tp",
            datapath=ipath_refp + "/P_1deg_",
            maskpos=False,  # do NOT set this to True!
            maskneg=True,
            uptake_years=fyears,
            uptake_dates=fdays,
            lats=totlats,
            lons=totlons,
        )
        gridcheck(totlats, reflats, totlons, reflons)

    ##--3. validation #########################################################
    if verbose:
        print(" * Starting validation...")

    init_netcdf(ofile, totlats, totlons)

    if fprec:
        if verbose:
            print("   P")
        pstats = calc_stats(Ptot, Pnparttot, Pref, thresh=0.001)
        write_to_netcdf(ofile, pstats, var="P")
    if fevap:
        if verbose:
            print("   E")
        estats = calc_stats(Etot, Enparttot, Eref, thresh=0.001)
        write_to_netcdf(ofile, estats, var="E")
    if fheat:
        if verbose:
            print("   H")
        hstats = calc_stats(Htot, Hnparttot, Href, thresh=1)
        write_to_netcdf(ofile, hstats, var="H")


###########################################################################
##--- run main script
###########################################################################
main_validation(
    ryyyy=args.ryyyy,
    ayyyy=args.ayyyy,
    am=args.am,
    opath_diag=path_diag,
    ipath_refp=path_refp,
    ipath_refe=path_refe,
    ipath_refh=path_refh,
    opath=path_diag,
    ofile_base=args.expid,
    fprec=args.fprec,
    fevap=args.fevap,
    fheat=args.fheat,
    verbose=args.verbose,
    veryverbose=args.veryverbose,
    fwrite_netcdf=args.write_netcdf,
)
