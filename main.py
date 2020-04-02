#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MAIN SCRIPT MOISTURE AND HEAT DIAGNOSIS
@author: dominik and jessica

To execute interactively: 
> exec(open("./main.py").read())

"""

###########################################################################
##--- MODULES
###########################################################################

import gzip
import pandas as pd
import numpy as np
import os, fnmatch
import timeit
import netCDF4 as nc4
import sys
import argparse
import time
import math as math
from datetime import datetime, timedelta, date
from math import sin,cos,acos,atan,atan2,sqrt,floor
from dateutil.relativedelta import relativedelta
import datetime as datetime
import imp
import warnings
from matplotlib import pyplot as plt
import cartopy.crs as ccrs
from cartopy.feature import NaturalEarthFeature

###########################################################################
##--- PATHS
###########################################################################

## determine working directory
wpath = os.getcwd()

## load input and output paths & input file name base(s)
with open(wpath+"/paths.txt") as f: 
    content = imp.load_source('','',f) # load like a python module
    ipath_REF = content.ipath_REF # input path (01_diagnosis)
    ipath_DGN = content.ipath_DGN # input path (01_diagnosis)
    ibase_DGN = content.ibase_DGN # input file name base(s)
    opath_DGN = content.opath_DGN # output path
    ipath_ATR = content.ipath_ATR # as above (for 02_attribution)
    ibase_ATR = content.ibase_ATR 
    opath_ATR = content.opath_ATR 
    opath_BIA = content.opath_BIA
    maskfile  = content.maskfile

###########################################################################
##--- MAIN
###########################################################################

os.chdir(wpath)

## (1) LOADING FUNCTIONS
exec(open("disclaimer.py").read())
exec(open("constants.py").read())
exec(open("metfunctions.py").read())
exec(open("01_diagnosis.py").read())
exec(open("02_attribution.py").read())
exec(open("03_biascorrection.py").read())
exec(open("hamsterfunctions.py").read())

## (2) get date, thresholds and flags from command line (job script) 
#      note: this is where we set the default values now. 
args    = read_cmdargs()
if args.ryyyy is None:
    args.ryyyy = args.ayyyy
verbose = args.verbose
print(printsettings(args,args.steps))

## (3) RUN main scripts with arguments
if args.steps == 1 or args.steps == 4:
    main_diagnosis(ryyyy=args.ryyyy, ayyyy=args.ayyyy, am=args.am, ad=args.ad,
              ipath=ipath_DGN,
              ifile_base=ibase_DGN, 
              opath=opath_DGN,
              ofile_base=args.expid,
              mode=args.mode,
              gres=args.gres,
              tdiagnosis=args.tdiagnosis,
              cheat_dtemp=args.cheat_dtemp,
              cheat_cc=args.cheat_cc,
              cevap_cc=args.cevap_cc,
              cevap_hgt=args.cevap_hgt,
              cheat_hgt=args.cheat_hgt,
              cprec_dqv=args.cprec_dqv,
              cprec_dtemp=args.cprec_dtemp,
              cprec_rh=args.cprec_rh,
              cpbl_strict=args.cpbl_strict,
              fjumps=args.fjumps,
              cjumps=args.cjumps,
              refdate=args.refdate,
              fwrite_netcdf=args.write_netcdf,
              precision=args.precision,
              ftimethis=args.timethis,
              fcc_advanced=args.cc_advanced,
              fvariable_mass=args.variable_mass,
              strargs=printsettings(args,1))

if args.steps == 2 or args.steps == 4:
    main_attribution(ryyyy=args.ryyyy, ayyyy=args.ayyyy, am=args.am, 
              ipath=ipath_ATR,
              ifile_base=ibase_ATR,
              opath=opath_ATR,
              ofile_base=args.expid,
              mode=args.mode,
              gres=args.gres,
              maskfile=maskfile,
              maskval=args.maskval,
              tdiagnosis=args.tdiagnosis,
              ctraj_len=args.ctraj_len,
              cheat_dtemp=args.cheat_dtemp,
              cheat_cc=args.cheat_cc, 
              cevap_cc=args.cevap_cc,
              cevap_hgt=args.cevap_hgt, 
              cheat_hgt=args.cheat_hgt,
              cprec_dqv=args.cprec_dqv, 
              cprec_dtemp=args.cprec_dtemp, 
              cprec_rh=args.cprec_rh,
              fjumps=args.fjumps,
              fjumpsfull=args.fjumpsfull,
              cjumps=args.cjumps,
              refdate=args.refdate,
              fwrite_netcdf=args.write_netcdf,
              precision=args.precision,
              ftimethis=args.timethis, 
              fdry=args.fallingdry,
              fmemento=args.memento,
              fcc_advanced=args.cc_advanced,
              fvariable_mass=args.variable_mass,
              strargs=printsettings(args,2))

if args.steps == 3 or args.steps == 4:
    main_biascorrection(ryyyy=args.ryyyy, ayyyy=args.ayyyy, am=args.am,
               opathA=opath_ATR, 
               opathD=opath_DGN, 
               ipathR=ipath_REF,
               opath=opath_BIA, 
               ofile_base=args.expid, # output
               maskfile=maskfile,
               maskval=args.maskval,
               set_negERA_to0=args.setnegzero,        # (only) makes sense for ERA-I data
               verbose=args.verbose,
               fdebug=args.debug,
               fwrite_netcdf=args.write_netcdf,
               precision=args.precision,
               strargs=printsettings(args,3))
