#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Main script to establish source-receptor relationships based on backward trajectories, 
# from e.g. FLEXPART. Two source regions are established: moisture and heat. 
# 
# This file is part of HAMSTER, 
# originally created by Dominik Schumacher, Jessica Keune, Diego G. Miralles
# at the Hydro-Climate Extremes Lab, Department of Environment, Ghent University
# 
# https://github.com/h-cel/hamster
# 
# HAMSTER is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation v3.
#
# HAMSTER is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with HAMSTER. If not, see <http://www.gnu.org/licenses/>.
#

import argparse
import calendar
import csv
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
import datetime as datetime
from functools import reduce
from math import acos, atan, atan2, cos, floor, sin, sqrt

import h5py
import netCDF4 as nc4
import numpy as np
import pandas as pd
from dateutil.relativedelta import relativedelta

from hamsterfunctions import *


def main_attribution(
    ryyyy,
    ayyyy,
    am,
    ad,
    ipath,
    ifile_base,
    ipath_f2t,
    opath,
    ofile_base,
    mode,
    gres,
    maskfile,
    maskval,
    verbose,
    veryverbose,
    ctraj_len,
    # E criteria
    cevap_dqv,
    fevap_drh,
    cevap_drh,
    cevap_hgt,
    # P criteria
    cprec_dqv,
    cprec_rh,
    # H criteria
    cheat_dtemp,
    fheat_drh,
    cheat_drh,
    cheat_hgt,
    fheat_rdq,
    cheat_rdq,
    # pbl and height criteria
    cpbl_method,
    cpbl_strict,
    cpbl_factor,
    refdate,
    fwrite_netcdf,
    precision,
    ftimethis,
    fdry,
    fmemento,
    mattribution,
    crandomnit,
    randatt_forcall,
    explainp,
    fdupscale,
    fmupscale,
    fvariable_mass,
    fwritestats,
    strargs,
):

    # TODO: add missing features
    if fvariable_mass:
        raise SystemExit("---- ABORTED: no can do, not implemented!")

    #### INPUT PATHS (incl. year)
    ipath_pp = os.path.join(ipath_f2t, str(ryyyy))  # raw partposit data
    ipath_tr = os.path.join(ipath, str(ryyyy))  # trajectory data (h5 files)

    #### OUTPUT FILES
    ## main netcdf output
    ofilename = (
        str(ofile_base)
        + "_attr_r"
        + str(ryyyy)[-2:]
        + "_"
        + str(ayyyy)
        + "-"
        + str(am).zfill(2)
        + ".nc"
    )
    ofile = opath + "/" + ofilename
    ## additional statistic output files (*.csv)
    # monthly statistics
    sfilename = (
        str(ofile_base)
        + "_attr_r"
        + str(ryyyy)[-2:]
        + "_"
        + str(ayyyy)
        + "-"
        + str(am).zfill(2)
        + "_stats.csv"
    )
    statfile = opath + "/" + sfilename
    # trajectory-based precipitation statistics
    if fwritestats and mattribution == "linear":
        pfilename = (
            str(ofile_base)
            + "_attr_r"
            + str(ryyyy)[-2:]
            + "_"
            + str(ayyyy)
            + "-"
            + str(am).zfill(2)
            + "_p-linear-attribution.csv"
        )
        pattfile = opath + "/" + pfilename
        with open(pattfile, "w") as pfile:
            writer = csv.writer(pfile, delimiter="\t", lineterminator="\n",)
            writer.writerow(["DATE", "F_ATT", "F_POT", "P_DQDT"])

    #### INPUT FILES
    ## read netcdf mask
    mask, mlat, mlon = maskgrabber(maskfile)

    #### DISCLAIMER
    if verbose:
        disclaimer()
        print("\n PROCESSING: \t", ayyyy, "-", str(am).zfill(2))
        print(
            "\n============================================================================================================\n"
        )
    ## Resets & consistency checks
    if fwritestats and mattribution == "random":
        print(
            " ! Option <writestats> not yet available for attribution method random. Continuing anyhow..."
        )
    if mode == "oper" and precision == "f4":
        precision = "f8"
        print(
            " ! Single precision should only be used for testing. Reset to double-precision."
        )
    if verbose:
        print(" ! using raw partposit input path: \t", ipath_pp)
        print(" ! using trajectory data input path: \t", ipath_tr)
        print(" ! using internal timer: \t" + str(ftimethis))
        print(" ! using mode: \t" + str(mode))
        print(" ! using attribution method (P): \t" + str(mattribution))
        print(" ! using trajectory-based upscaling: \t" + str(explainp))
        if mattribution == "random":
            print(" ! using minimum iterations: \t" + str(crandomnit))
        print(" ! using daily upscaling: \t" + str(fdupscale))
        print(" ! using monthly upscaling: \t" + str(fmupscale))
        if fwrite_netcdf:
            print(" ! writing netcdf output: \t" + str(fwrite_netcdf))
            print(" \t ! with grid resolution: \t", str(gres))
            print(" \t ! output file: \t", opath + "/" + ofilename)
        print(" ! additional statistics in: \t" + str(statfile))
        if fwritestats and mattribution == "linear":
            print(" ! precipitation statistics in: \t" + str(pattfile))
        print(
            "\n============================================================================================================"
        )
        print(
            "\n============================================================================================================"
        )

    #### SETTINGS
    ## start timer
    if ftimethis:
        megatic = timeit.default_timer()

    ## grids
    glon, glat, garea = makegrid(resolution=gres)
    gridcheck(glat, mlat, glon, mlon)

    ## -- DATES
    dt = 6  # hardcoded for FLEXPART ERA-INTERIM with 6h
    timestep = datetime.timedelta(hours=dt)

    # get start date (to read trajectories) - NOTE: we begin at 06 UTC...
    sdate_bgn = (
        str(ayyyy)
        + "-"
        + str(am).zfill(2)
        + "-"
        + str(ad).zfill(2)
        + "-"
        + str(dt).zfill(2)
    )
    datetime_bgn = datetime.datetime.strptime(sdate_bgn, "%Y-%m-%d-%H")
    # get end date (to read trajectories) - NOTE: always 00 UTC of the 1st of the next month
    nayyyy, nam = nextmonth(datetime_bgn)
    sdate_end = str(nayyyy) + "-" + str(nam).zfill(2) + "-01-00"
    datetime_end = datetime.datetime.strptime(sdate_end, "%Y-%m-%d-%H")

    # file dates (arrival, 6h seq)
    datetime_seq, fdatetime_seq, ffdatetime_seq = timelord(
        datetime_bgn, datetime_end, timestep
    )
    # daily arrival dates (24h seq for netCDF writing)
    fdate_seq = timelord(
        datetime_bgn - timestep,
        datetime_end - timestep,
        datetime.timedelta(hours=24),
        ret="datetime",
    )
    fdateasdate = datetime2date(fdate_seq)

    # NOTE: better to keep these as lists to maintain consistency

    # calculate number of time steps, also aggregated to daily resolution
    ntime = len(fdatetime_seq)
    ndaytime = len(fdate_seq)

    ## TESTMODE
    if mode == "test":

        ctraj_len_orig = ctraj_len
        ctraj_len = (
            2  # NOTE: changes here must be accompanied by changes in 01_diagnosis!
        )

        ntime = 4  # NOTE: use multiples of 4 only, else output is not saved
        datetime_seq = datetime_seq[(4 * ctraj_len) : (4 * ctraj_len + ntime)]
        fdatetime_seq = fdatetime_seq[(4 * ctraj_len) : (4 * ctraj_len + ntime)]
        ndaytime = int(ntime / 4)
        fdate_seq = fdate_seq[ctraj_len : ctraj_len + ndaytime]
        fdateasdate = fdateasdate[ctraj_len : ctraj_len + ndaytime]

    ## -- WRITE NETCDF OUTPUT (empty, to be filled)
    if fwrite_netcdf:
        writeemptync4D(
            ofile, fdate_seq, np.arange(-ctraj_len, 1), glat, glon, strargs, precision
        )

    # traj max len, expressed in input data (6-hourly) steps
    tml = int(4 * ctraj_len)  # hardcoded for 6-hourly input
    # compact form of max traj len in days (used for array filling w/ shortened uptake dim)
    ctl = ctraj_len

    ### MEMENTO --- pre-loop to produce independent monthly output
    ## NOTE: this is irrelevant for E2P, but crucial for Had (& Ead)
    ## NOTE: we only need to know if some parcel makes it to the ABL, that's it!
    ## NOTE: must fill array with negative number whose abs exceeds max traj len
    if fmemento:
        pidlog = -999 * np.ones(shape=2100000).astype(int)

        if mode == "oper":  # skip if multi-counting somehow desired and/or if testing
            ###--- PRELOOP v2
            # NOTE: code further below this function call here could also be moved to
            #       first main loop iteration, so that no dim checking necessary
            ntrajstep = readtraj(
                idate=datetime_seq[0],
                ipath=ipath_tr,
                ifile_base=ifile_base,
                verbose=False,
            ).shape[0]

            # I believe this could be done for parcels of interest / pot. conflict only... but we'll leave it like this for now
            if ntrajstep < tml + 2 + 4:
                # only do this if data really isn't already 'there'
                preloop_dates = timelord(
                    fdatetime_seq[0] - (tml + 5) * timestep,
                    fdatetime_seq[0] - timestep,
                    timestep,
                    ret="fileformat",
                )
                preloop_files = [
                    ipath_pp + "/partposit_" + idfile + ".gz"
                    for idfile in preloop_dates
                ]
                extendarchive = grabmesomehpbl(filelist=preloop_files, verbose=verbose)
            else:
                if verbose:
                    print(
                        "\n=== \t INFO: no pre-loop needed, trajectories are long enough"
                    )

    ###--- MAIN LOOP

    ## prepare STATS
    # number of parcels
    tneval = tnnevala = tnnevalm = tnevalp = tnnevalp = tnevalh = tnnevalh = 0
    # precip. statistics
    psum = patt = punatt = pmiss = 0

    ## loop over time to read in files
    if verbose:
        print("\n=== \t Start main program...\n")
    for ix in range(ntime):
        if verbose:
            print(
                "--------------------------------------------------------------------------------------"
            )
        print("Processing " + str(fdatetime_seq[ix]))

        # safety exit hardcoded for FLEXPART–ERA-Interim runs
        if ryyyy == ayyyy and am == 12 and ix == range(ntime)[-1]:
            continue

        ## 1) read in all files associated with data --> ary is of dimension (ntrajlen x nparcels x nvars)
        ary = readtraj(
            idate=datetime_seq[ix],
            ipath=ipath_tr,
            ifile_base=ifile_base,
            verbose=verbose,
        )

        nparcel = ary.shape[1]
        ntrajleng = ary.shape[0]
        if verbose:
            print(
                " TOTAL: " + str(datetime_seq[ix]) + " has " + str(nparcel) + " parcels"
            )

        if mode == "test":
            ntot = range(1000)
        else:
            ntot = range(nparcel)

        # figure out where to store data (on which arriving day)
        arv_idx = np.where(
            np.asarray(fdateasdate)
            == (fdatetime_seq[ix] - relativedelta(hours=3)).date()
        )[0][0]

        # pre-allocate arrays (repeat at every 4th step)
        if ix % 4 == 0:
            ary_heat = np.zeros(shape=(ctl + 1, glat.size, glon.size))
            ary_etop = np.zeros(shape=(ctl + 1, glat.size, glon.size))
            # upscaling measures (currently has to be per day as well)
            if fdupscale:
                ipatt = ipmiss = 0

        # STATS: number of parcels per file
        neval = nnevala = nnevalm = nevalp = nnevalp = nevalh = nnevalh = 0

        # grab extended trajectory data
        ### CHECK/CLEANUP: CAN ALL OF THIS GO?
        # (i)  POM data not supported anymore
        # (ii) if trajectories too short: abort (let's not fix everything right away for the user)
        if fmemento:
            # NOTE: pom data can come with duplicate IDs; remove to avoid (some) trouble
            if not np.unique(ary[0, :, 0]).size == ary[0, :, 0].size:
                print(
                    "\t INFO: duplicates detected, original pom array shape=", ary.shape
                )
                _, ikeep = np.unique(ary[0, :, 0], return_index=True)
                ary = ary[:, ikeep, :]
                print("\t INFO: duplicates eliminated, new pom array shape=", ary.shape)
                nparcel = ary.shape[1]  # update
                ntot = range(nparcel)
            # NOTE: yet another ******* pom fix ... to be removed!
            #       (flex2traj output is not affected by this)
            IDs = ary[0, :, 0]
            thresidx = int((9 / 10) * ary.shape[1])  # this should do the trick
            # simply shift to indices > 2e6
            IDs[thresidx:][IDs[thresidx:] < 3000] += 2e6

            # extract what is needed from extendarchive if trajs 'too short'
            if ary.shape[0] < tml + 2 + 4 and ix < ctraj_len * 4:
                extendtrajs = np.empty(shape=(4, nparcel, 2))
                for pp in range(4):
                    allIDs = extendarchive[-(4 - pp + ix)][:, 0]
                    extendtrajs[pp, :, 0] = extendarchive[-(4 - pp + ix)][:, 0][
                        np.where(np.isin(allIDs, ary[0, :, 0]))
                    ]  # ID
                    extendtrajs[pp, :, 1] = extendarchive[-(4 - pp + ix)][:, 1][
                        np.where(np.isin(allIDs, ary[0, :, 0]))
                    ]  # hpbl

        ## 2) establish source–receptor relationships
        for i in ntot:

            ## - 2.0) only evaluate if the parcel is in target region (midpoint or arrival point)
            mlat_ind, mlon_ind = midpindex(ary[:2, i, :], glon=mlon, glat=mlat)
            alat_ind, alon_ind = arrpindex(ary[0, i, :], glon=mlon, glat=mlat)
            if (
                not mask[alat_ind, alon_ind] == maskval
                and not mask[mlat_ind, mlon_ind] == maskval
            ):
                nnevalm += 1
                nnevala += 1
                continue

            ## - 2.1) check how far back trajectory should be evaluated
            if mask[alat_ind, alon_ind] == maskval:
                ID = int(ary[0, i, 0])
                if fmemento and ary[0, i, 3] < np.max(ary[:4, i, 7]):
                    if ix < ctraj_len * 4:  # rely on (extended) traj data

                        # check if parcel has been inside before
                        is_inmask = whereinmask(
                            mask=mask,
                            maskval=maskval,
                            masklat=mlat,
                            masklon=mlon,
                            trajlat=ary[: (tml + 2), i, 2],
                            trajlon=ary[: (tml + 2), i, 1],
                        )

                        ### CHECK/CLEANUP: CAN ALL OF THIS GO? (if only trajs > length allowed)?
                        # check if parcel was 'in PBL'
                        hgt = ary[: (tml + 2), i, 3]  # consistent with max traj len
                        if ary.shape[0] < tml + 2 + 4:
                            longhpbl = np.concatenate(
                                (ary[: (tml + 2), i, 7], extendtrajs[:, i, 1])
                            )
                        else:
                            longhpbl = ary[: (tml + 2 + 4), i, 7]
                        is_inpbl = np.where(hgt < maxlastn(longhpbl, n=4)[:-4])[
                            0
                        ]  # omit last 4

                        # check where parcel inside and 'in PBL'
                        is_arrv = np.intersect1d(is_inmask, is_inpbl)

                        # now determine ihf_H
                        if is_arrv.size > 1:
                            ihf_H = is_arrv[is_arrv > 0].min() + 1
                        elif is_arrv == 0:
                            ihf_H = tml + 2  # use max traj len
                        else:
                            raise RuntimeError(
                                "--- FATAL ERROR: Schrödingers cat situation; is parcel inside and in PBL, or not?"
                            )

                    else:  # fully rely on log from now
                        istep = pidlog[ID]
                        ihf_H = min((ix - istep + 1), tml + 2)
                else:
                    ihf_H = tml + 2

            ## - 2.2) read only the most basic parcel information
            # NOTE: this could easily be done more efficiently
            hgt, hpbl, temp, qv, dens, pres = glanceparcel(ary[:4, i, :])

            # sorry, yet another date for writing the P date to the csv (preliminary).
            # because i wanted to have the hours in there as wel (not assign to day only)
            pdate = str(
                (fdatetime_seq[ix] - relativedelta(hours=3)).strftime("%Y%m%d%H")
            )

            ## - 2.3) diagnose fluxes

            ## (a) E2P, evaporation resulting in precipitation
            if not mask[mlat_ind, mlon_ind] == maskval:
                nnevalm += 1
            else:
                if (qv[0] - qv[1]) < cprec_dqv and (
                    (q2rh(qv[0], pres[0], temp[0]) + q2rh(qv[1], pres[1], temp[1])) / 2
                ) > cprec_rh:

                    # prec
                    prec = abs(qv[0] - qv[1])
                    # log some statistics
                    nevalp += 1
                    psum += prec

                    # read full parcel information
                    (
                        lons,
                        lats,
                        temp,
                        hgt,
                        qv,
                        hpbl,
                        dens,
                        pres,
                        pottemp,
                        epottemp,
                    ) = readparcel(ary[: tml + 2, i, :])
                    rh = q2rh(qv, pres, temp)

                    # calculate all required changes along trajectory
                    dq = trajparceldiff(qv[:], "diff")
                    # evaluate only until trajectory falls dry
                    ihf_E = tml + 2
                    if fdry and np.any(qv[1:ihf_E] <= 0.00005):
                        ihf_E = np.min(np.where(qv[1:ihf_E] <= 0.00005)[0] + 1)

                    # identify uptake locations
                    is_inpbl = pblcheck(
                        np.dstack((hgt[:ihf_E], hpbl[:ihf_E]))[0, :, :],
                        cpbl_strict,
                        minh=cheat_hgt,
                        fpbl=cpbl_factor,
                        method=cpbl_method,
                    )
                    is_drh = drhcheck(rh[:ihf_E], checkit=fevap_drh, maxdrh=cevap_drh)
                    is_uptk = dq[: ihf_E - 1] > cevap_dqv
                    evap_idx = np.where(
                        np.logical_and(is_inpbl, np.logical_and(is_drh, is_uptk))
                    )[0]

                    # log some stats if trajectory is without any uptakes (for upscaling)
                    if evap_idx.size == 0:
                        nnevalp += 1
                        pmiss += prec
                        if fdupscale:
                            ipmiss += prec

                    # ATTRIBUTION
                    if evap_idx.size > 0:
                        if mattribution == "linear":
                            etop = linear_attribution_p(
                                qv[:ihf_E], iupt=evap_idx, explainp=explainp
                            )
                        elif mattribution == "random":
                            etop = random_attribution_p(
                                qtot=qv[:ihf_E],
                                iupt=evap_idx,
                                explainp=explainp,
                                nmin=crandomnit,
                                forc_all=randatt_forcall,
                                verbose=verbose,
                                veryverbose=veryverbose,
                            )
                        # write to grid
                        for itj in evap_idx:
                            ary_etop[ctl - (itj + 3 - ix % 4) // 4, :, :] += gridder(
                                plon=lons[itj : itj + 2],
                                plat=lats[itj : itj + 2],
                                pval=etop[itj],
                                glon=glon,
                                glat=glat,
                            )

                        # write additional stats to csv-file (currently: ALWAYS explain="none"; also: why only linear?)
                        if fwritestats and mattribution == "linear":
                            if evap_idx.size == 0:
                                pattdata = [pdate, str(0), str(0), str(prec)]
                            elif evap_idx.size > 0:
                                etop = linear_attribution_p(
                                    qv[:ihf_E], iupt=evap_idx, explainp="none"
                                )
                                pattdata = [
                                    pdate,
                                    str(np.sum(etop[evap_idx] / prec)),
                                    str(1 - etop[-1] / prec),
                                    str(prec),
                                ]
                            append2csv(pattfile, pattdata)

                        # log some statistics (for upscaling)
                        patt += np.sum(etop[evap_idx])
                        punatt += prec - np.sum(etop[evap_idx])
                        if fdupscale:
                            ipatt += np.sum(etop[evap_idx])
                            ipmiss += prec - np.sum(etop[evap_idx])

            ## (b) H, surface sensible heat arriving in PBL (or nocturnal layer)
            if not mask[alat_ind, alon_ind] == maskval:
                nnevala += 1
            else:
                if ihf_H >= 2 and hgt[0] < np.max(hpbl[:4]):

                    # log some statistics
                    nevalh += 1

                    # read full parcel information
                    (
                        lons,
                        lats,
                        temp,
                        hgt,
                        qv,
                        hpbl,
                        dens,
                        pres,
                        pottemp,
                        epottemp,
                    ) = readparcel(ary[:ihf_H, i, :])
                    rh = q2rh(qv, pres, temp)

                    # calculate all required changes along trajectory
                    dTH = trajparceldiff(pottemp[:], "diff")

                    # identify sensible heat uptakes
                    is_inpbl = pblcheck(
                        np.dstack((hgt[:ihf_H], hpbl[:ihf_H]))[0, :, :],
                        cpbl_strict,
                        minh=cheat_hgt,
                        fpbl=cpbl_factor,
                        method=cpbl_method,
                    )
                    is_drh = drhcheck(rh[:ihf_H], checkit=fheat_drh, maxdrh=cheat_drh)
                    is_rdqv = rdqvcheck(
                        qv[:ihf_H], checkit=fheat_rdq, maxrdqv=cheat_rdq
                    )
                    is_uptk = dTH[: ihf_H - 1] > cheat_dtemp
                    heat_idx = np.where(
                        np.logical_and(
                            is_inpbl,
                            np.logical_and(is_drh, np.logical_and(is_rdqv, is_uptk)),
                        )
                    )[0]

                    # discount uptakes linearly
                    if heat_idx.size == 0:
                        # log some statistics
                        nnevalh += 1
                    if heat_idx.size > 0:
                        dTH_disc = linear_discounter(v=pottemp[:ihf_H], min_gain=0)

                    # loop through sensible heat uptakes
                    for itj in heat_idx:
                        # NOTE: hardcoded for writing daily data
                        ary_heat[ctl - (itj + 3 - ix % 4) // 4, :, :] += (
                            gridder(
                                plon=lons[itj : itj + 2],
                                plat=lats[itj : itj + 2],
                                pval=dTH_disc[itj],
                                glon=glon,
                                glat=glat,
                            )
                            / 4
                        )

                    # update parcel log
                    if fmemento:
                        pidlog[ID] = ix  # NOTE: double-check

        neval = len(ntot)
        if verbose:
            print(
                " STATS: Evaluated "
                + str(neval - nnevala)
                + " ({:.2f}".format(100 * (neval - nnevala) / (neval))
                + "%) arriving parcels inside mask (advection)."
            )
            if nnevalh != 0:
                print(
                    " --- ATTENTION: "
                    + str(nnevalh)
                    + "/"
                    + str(neval - nnevala)
                    + " arriving parcels are not associated with any heat uptakes..."
                )
            print(
                " STATS: Evaluated "
                + str(neval - nnevalm)
                + " ({:.2f}".format(100 * (neval - nnevalm) / (neval))
                + "%) midpoint parcels inside mask (precipitation)."
            )
            if nnevalm == neval:
                print(" STATS: Evaluated " + str(nevalp) + " precipitating parcels.")
            else:
                print(
                    " STATS: Evaluated "
                    + str(nevalp)
                    + " ({:.2f}".format(100 * (nevalp) / (neval - nnevalm))
                    + "%) precipitating parcels."
                )
            if nnevalp != 0:
                print(
                    " --- ATTENTION: "
                    + str(nnevalp)
                    + "/"
                    + str(nevalp)
                    + " precipitating parcels are not associated with any evap uptakes..."
                )

        ## SOME DAILY CALCULATIONS
        if (ix + 1) % 4 == 0:
            # DAILY UPSCALING of E2P, taking into account the missing trajectories (i.e. the ones without any uptakes)
            if fdupscale and (nnevalp != 0 or ipmiss != 0):
                if ipatt == 0:
                    print(
                        " \n--- WARNING: there were no trajectories with uptakes, so upscaling is impossible...\n "
                    )
                else:
                    upsfac = 1 + (ipmiss / ipatt)
                    ary_etop[:, :, :] = upsfac * ary_etop[:, :, :]
                    # corrections for final statistics
                    patt += -np.sum(ipatt) + np.sum(ary_etop)
                    pmiss += -np.sum(ipmiss)
                    if verbose:
                        print(" * Upscaling... (factor: {:.4f}".format(upsfac) + ")")
            # Convert units
            if verbose:
                print(" * Converting units...")
            ary_etop[:, :, :] = convertunits(ary_etop[:, :, :], garea, "E")
            ary_heat[:, :, :] = convertunits(ary_heat[:, :, :], garea, "H")

            if fwrite_netcdf:
                writenc4D(ofile, arv_idx, ary_etop, ary_heat, verbose)

        ## STATS summary
        tneval += neval
        tnnevala += nnevala
        tnnevalm += nnevalm
        tnevalp += nevalp
        tnnevalp += nnevalp
        tnevalh += nevalh
        tnnevalh += nnevalh

    # MONTHLY UPSCALING of E2P, taking into account the missing trajectories (i.e. the ones without any uptakes)
    if fmupscale and (pmiss != 0 or punatt != 0):
        if patt == 0:
            print(
                " \n--- WARNING: there were no trajectories with uptakes, so upscaling is impossible...\n"
            )
        else:
            upsfac = 1 + ((pmiss + punatt) / patt)
            # load full etop array and upscale
            fdata = nc4.Dataset(ofile, "r+")
            uns_etop = fdata.variables["E2P"][:]
            ups_etop = upsfac * uns_etop
            fdata.variables["E2P"][:] = ups_etop
            fdata.close()
            # corrections for final statistics
            patt += np.sum(pmiss) + np.sum(punatt)
            pmiss += -np.sum(pmiss)
            punatt += -np.sum(punatt)
            if verbose:
                print(
                    " * Monthly upscaling for unattributed precipitation... (factor: {:.4f}".format(
                        upsfac
                    )
                    + ")"
                )

    if ftimethis:
        megatoc = timeit.default_timer()
        if verbose:
            print(
                "\n=== \t End main program (total runtime so far: ",
                str(round(megatoc - megatic, 2)),
                "seconds) \n",
            )

    if verbose:
        if fwrite_netcdf:
            if psum != 0:
                # adding attributed fraction to netcdf file description
                append_attrfrac_netcdf(ofile, "{:.2f}".format(patt / psum))
            print("\n Successfully written: " + ofile + " !\n")

    writestats_02(
        statfile,
        tneval,
        tnnevala,
        tnevalh,
        tnnevalh,
        tnnevalm,
        tnevalp,
        tnnevalp,
        patt,
        psum,
        punatt,
        pmiss,
    )
    if verbose:
        with open(statfile, "r") as sfile:
            print(sfile.read())
