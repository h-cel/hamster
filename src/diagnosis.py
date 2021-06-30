#!/usr/bin/env python3
# -*- coding: utf-8 -*-

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

from hamsterfunctions import *


def main_diagnosis(
           ryyyy, ayyyy, am, ad,
           ipath, ifile_base,
           opath, ofile_base,
           mode,
           gres,
           verbose,
           veryverbose,
           fproc_npart,
           # E criteria
           fevap, cevap_dqv, fevap_drh, cevap_drh, cevap_hgt,
           # P criteria
           fprec, cprec_dqv, cprec_rh,
           # H criteria
           fheat, cheat_dtemp, fheat_drh, cheat_drh, cheat_hgt, fheat_rdq, cheat_rdq,
           # pbl and height criteria
           cpbl_method, cpbl_strict, cpbl_factor,
           refdate,
           fwrite_netcdf,
           precision,
           ftimethis,fvariable_mass,
           strargs):

    ## Perform consistency checks
    if mode=="oper" and precision=="f4":
        precision = "f8"
        print("Single precision should only be used for testing. Reset to double-precision.")
    if fvariable_mass and not fproc_npart:
        fproc_npart = True
        print("Have to process all parcels for variable mass...")

    ## Construct precise input and storage paths
    mainpath  = ipath+str(ryyyy)+"/"
    ofilename = str(ofile_base)+"_diag_r"+str(ryyyy)[-2:]+"_"+str(ayyyy)+"-"+str(am).zfill(2)+".nc"
    ofile     = opath+"/"+ofilename
    
    if verbose:
        disclaimer()
        print("\n PROCESSING: \t", 	ayyyy, "-", str(am).zfill(2))
        print("\n============================================================================================================")
        print(" ! using input path: \t", 	ipath)
        print(" ! using variable mass: \t" +str(fvariable_mass) )
        if fvariable_mass:
            print(" \t ! reference date for number of particles: \t" +str(refdate) )
        if fwrite_netcdf:
            print(" ! writing netcdf output: \t" +str(fwrite_netcdf) )
            print(" \t ! with grid resolution: \t", str(gres) )
            print(" \t ! output file: \t", opath+"/"+ofilename)
        print(" ! using internal timer: \t" +str(ftimethis) )
        print(" ! using mode: \t" +str(mode))
        print("\n============================================================================================================")
        print("\n============================================================================================================")

    ## Start timer
    if ftimethis:
        megatic = timeit.default_timer()

    ## Prepare grid
    glon, glat, garea = makegrid(resolution=gres)

    ## Handle dates
    date_bgn        = datetime.datetime.strptime(str(ayyyy)+"-"+str(am).zfill(2)+"-"+str(ad).zfill(2), "%Y-%m-%d")
    # get end date (always 00 UTC of the 1st of the next month)
    nayyyy          = (date_bgn + relativedelta(months=1)).strftime('%Y')
    nam             = (date_bgn + relativedelta(months=1)).strftime('%m')
    date_end        = datetime.datetime.strptime(str(nayyyy)+"-"+str(nam).zfill(2)+"-01-00",  "%Y-%m-%d-%H")
    timestep        = datetime.timedelta(hours=6)
    date_seq        = []
    fdate_seq       = []
    mfdate_seq      = []
    idate           = date_bgn + timestep
    while idate <= date_end:
        date_seq.append(idate.strftime('%Y%m%d%H'))
        fdate_seq.append(idate)
        mfdate_seq.append(idate-timestep/2) # -dt/2 for backward run
        idate   += timestep
    ntime           = len(date_seq)

    ##-- TESTMODE
    if mode == "test":
        ntime       = 12
        date_seq    = date_seq[0:ntime]
        fdate_seq   = fdate_seq[0:ntime]
        mfdate_seq  = mfdate_seq[0:ntime]

    ## Create empty netcdf file (to be filled)
    if fwrite_netcdf:
        writeemptync(ofile,mfdate_seq,glon,glat,strargs,precision,fproc_npart,fprec,fevap,fheat)

    # Read in reference distribution of parcels
    if fvariable_mass:
        print(" \n !!! WARNING !!! With this version, variable mass can only be applied to 01_diagnosis -- it cannot be used consistently for all steps yet! \n")
        ary_rnpart   = get_refnpart(refdate=refdate, ryyyy=ryyyy, glon=glon, glat=glat)
    
    ##-- LOOP THROUGH FILES
    if verbose:
        print("\n=== \t Start main program: 01_diagnosis...\n")

    for ix in range(ntime):
        
        if verbose:
            print("--------------------------------------------------------------------------------------")
            print("Processing "+str(fdate_seq[ix]))

        ## Read date related trajectories -> ary is of dimension (ntrajlen x nparticles x nvars)
        ary         = readtraj(idate        = date_seq[ix], 
                               ipath        = ipath+"/"+str(ryyyy), 
                               ifile_base   = ifile_base,
                               verbose      = verbose)
        ary         = calc_allvars(ary)
        dq          = trajparceldiff(ary[:,:,5], "diff")
        mrh         = np.apply_over_axes(np.mean, ary[:,:,10], 0)
        dTH         = trajparceldiff(ary[:,:,11], "diff")

        nparticle   = ary.shape[1]
        if verbose:
            print(" TOTAL: " + str(date_seq[ix]) + " has " + str(nparticle) + " parcels")

        ## TESTMODE: less parcels
        if mode == "test":
            ntot    = range(10000,10100)
        else:
            ntot    = range(nparticle)

        #smalltic = timeit.default_timer()

        ##-- LOOP OVER PARCELS TO DIAGNOSE P, E, H (and npart) and assign to grid
        if fproc_npart:
            # get midpoint indices on grid from ary
            imidi           = get_all_midpindices(ary, glon, glat)
            ary_npart       = gridall(imidi[:,1], imidi[:,0], np.repeat(1,nparticle), glon=glon, glat=glat)
        elif not fproc_npart:
            # currently just writing empty array ... to be changed
            ary_npart       = np.zeros(shape=(glat.size,glon.size))
        
        ## Precipitation
        if fprec:
            fdqv        = np.where(dq[0,:]<cprec_dqv)
            frh         = np.where(mrh[0,:]>cprec_rh)
            isprec      = np.intersect1d(fdqv,frh)
            p_ary       = ary[:,isprec,:]
            pmidi       = get_all_midpindices(p_ary, glon, glat)
            # grid
            ary_prec    = gridall(pmidi[:,1], pmidi[:,0], dq[:,isprec][0], glon=glon, glat=glat)
            ary_pnpart  = gridall(pmidi[:,1], pmidi[:,0], np.repeat(1,isprec.size), glon=glon, glat=glat)
            
        ## Evaporation
        if fevap:
            isevap      = filter_for_evap_parcels(ary, dq, cpbl_method, cpbl_strict, cpbl_factor, cevap_hgt, fevap_drh, cevap_drh, cevap_dqv, veryverbose)
            e_ary       = ary[:,isevap,:]
            emidi       = get_all_midpindices(e_ary, glon, glat)
            # grid
            ary_evap    = gridall(emidi[:,1], emidi[:,0], dq[:,isevap][0], glon=glon, glat=glat)
            ary_enpart  = gridall(emidi[:,1], emidi[:,0], np.repeat(1,isevap.size), glon=glon, glat=glat)

        ## Sensible heat
        if fheat:
            isheat      = filter_for_heat_parcels(ary, dTH, cpbl_method, cpbl_strict, cpbl_factor, cheat_hgt, fheat_drh, cheat_drh, cheat_dtemp, fheat_rdq, cheat_rdq, veryverbose)
            h_ary       = ary[:,isheat,:]
            hmidi       = get_all_midpindices(h_ary, glon, glat)
            # grid
            ary_heat    = gridall(hmidi[:,1], hmidi[:,0], dTH[:,isheat][0], glon=glon, glat=glat)
            ary_hnpart  = gridall(hmidi[:,1], hmidi[:,0], np.repeat(1,isheat.size), glon=glon, glat=glat)

        #smalltoc = timeit.default_timer()
        #print("=== \t All parcels: ",str(round(smalltoc-smalltic, 2)),"seconds \n")

        ## Convert units
        if verbose:
            print(" * Converting units...")
        if fprec:    
            ary_prec[:,:] = convertunits(ary_prec[:,:], garea, "P")
        if fevap:
            ary_evap[:,:] = convertunits(ary_evap[:,:], garea, "E")
        if fheat:
            ary_heat[:,:] = convertunits(ary_heat[:,:], garea, "H")

        ## Scale with parcel mass
        if fvariable_mass:
            print(" !!! WARNING !!! With this version, variable mass can only be applied to 01_diagnosis -- it cannot be used consistently for all steps yet!")
            if verbose: 
                print(" * Applying variable mass...")
            if fprec:    
                ary_prec[:,:]         = scale_mass(ary_prec[:,:], ary_npart[:,:], ary_rnpart)
            if fevap:
                ary_evap[:,:]         = scale_mass(ary_evap[:,:], ary_npart[:,:], ary_rnpart)
            if fheat:
                ary_heat[:,:]         = scale_mass(ary_heat[:,:], ary_npart[:,:], ary_rnpart)

        # write to netcdf file
        if fwrite_netcdf:
            #writenc(ofile,ix,ary_prec[:,:],ary_evap[:,:],ary_heat[:,:],ary_npart[:,:],ary_pnpart[:,:],ary_enpart[:,:],ary_hnpart[:,:])
            if fproc_npart:
                writenc(ofile,ix,ary_npart[:,:],'n_part')
            if fprec:
                writenc(ofile,ix,ary_prec[:,:],'P')
                writenc(ofile,ix,ary_pnpart[:,:],'P_n_part')
            if fevap:
                writenc(ofile,ix,ary_evap[:,:],'E')
                writenc(ofile,ix,ary_enpart[:,:],'E_n_part')
            if fheat:
                writenc(ofile,ix,ary_heat[:,:],'H')
                writenc(ofile,ix,ary_hnpart[:,:],'H_n_part')

        ## re-init. arrays
        if fproc_npart:
            ary_npart[:,:]  = 0
        if fprec:
            ary_pnpart[:,:] = 0
            ary_prec[:,:]   = 0
        if fevap:
            ary_enpart[:,:] = 0
            ary_evap[:,:]   = 0
        if fheat:    
            ary_hnpart[:,:] = 0
            ary_heat[:,:]   = 0

    if ftimethis:
        megatoc = timeit.default_timer()
        if verbose:
            print("\n=== \t End main program (total runtime so far: ",str(round(megatoc-megatic, 2)),"seconds) \n")

    if verbose:
        if fwrite_netcdf:
            print("\n Successfully written: "+ofile+" !")
