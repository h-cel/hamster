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
import csv
import random
import struct
import calendar
import h5py
import re

###########################################################################
##--- PATHS
###########################################################################

## determine working directory
wpath = os.getcwd()

## load input and output paths & input file name base(s)
content = imp.load_source('',wpath+"/paths.txt") # load like a python module
ipath_REF = content.ipath_REF # input path (01_diagnosis)
ipath_DGN = content.ipath_DGN # input path (01_diagnosis)
ibase_DGN = content.ibase_DGN # input file name base(s)
opath_DGN = content.opath_DGN # output path
ipath_ATR = content.ipath_ATR # as above (for 02_attribution)
ibase_ATR = content.ibase_ATR
opath_ATR = content.opath_ATR
opath_BIA = content.opath_BIA
maskfile  = content.maskfile
ibase_f2t = content.ibase_f2t
ipath_f2t = content.ipath_f2t
opath_f2t = content.opath_f2t
wpath_f2t = wpath

###########################################################################
##--- MAIN
###########################################################################

os.chdir(wpath)

## (1) LOADING FUNCTIONS
exec(open("disclaimer.py").read())
exec(open("constants.py").read())
exec(open("metfunctions.py").read())
exec(open("00_flex2traj.py").read())
exec(open("01_diagnosis.py").read())
exec(open("02_attribution.py").read())
exec(open("03_biascorrection.py").read())
exec(open("hamsterfunctions.py").read())

## (2) get date, thresholds and flags from command line (job script) 
#      note: this is where we set the default values now. 
args    = read_cmdargs()
verbose = args.verbose
print(printsettings(args,args.steps))

# just waiting a random number of seconds (max. 30s)
# to avoid overlap of path.exist and makedirs between parallel jobs (any better solution?)
waiter  = random.randint(0,30)
time.sleep(waiter)

# create output directories if they do not exist (in dependency of step)
if args.steps==0 and not os.path.exists(opath_f2t):
        os.makedirs(opath_f2t)
if args.steps==1 and not os.path.exists(opath_DGN):
        os.makedirs(opath_DGN)
if args.steps==2 and not os.path.exists(opath_ATR):
        os.makedirs(opath_ATR)
if args.steps==3 and not os.path.exists(opath_BIA):
        os.makedirs(opath_BIA)

## (3) RUN main scripts with arguments
if args.steps ==0:
    main_flex2traj(ryyyy=args.ryyyy, ayyyy=args.ayyyy, am=args.am, ad=args.ad,
                   tml=args.ctraj_len,
                   maskfile=maskfile,
                   maskval=args.maskval,
                   cpbl_strict=args.cpbl_strict,
                   idir=ipath_f2t,
                   odir=opath_f2t,
                   fout=ibase_f2t,
                   workdir=wpath_f2t)

if args.steps == 1:
    main_diagnosis(ryyyy=args.ryyyy, ayyyy=args.ayyyy, am=args.am, ad=args.ad,
              ipath=ipath_DGN,
              ifile_base=ibase_DGN, 
              opath=opath_DGN,
              ofile_base=args.expid,
              mode=args.mode,
              gres=args.gres,
              verbose=args.verbose,
              veryverbose=args.veryverbose,
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
              refdate=args.refdate,
              fwrite_netcdf=args.write_netcdf,
              precision=args.precision,
              ftimethis=args.timethis,
              fcc_advanced=args.cc_advanced,
              fvariable_mass=args.variable_mass,
              strargs=printsettings(args,1))

if args.steps == 2:
    main_attribution(ryyyy=args.ryyyy, ayyyy=args.ayyyy, am=args.am, ad=args.ad, 
              ipath=ipath_ATR,
              ifile_base=ibase_ATR,
              ipath_f2t=ipath_f2t,
              opath=opath_ATR,
              ofile_base=args.expid,
              mode=args.mode,
              gres=args.gres,
              maskfile=maskfile,
              maskval=args.maskval,
              verbose=args.verbose,
              veryverbose=args.veryverbose,
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
              cpbl_strict=args.cpbl_strict,
              refdate=args.refdate,
              fwrite_netcdf=args.write_netcdf,
              precision=args.precision,
              ftimethis=args.timethis, 
              fdry=args.fallingdry,
              fmemento=args.memento,
              mattribution=args.mattribution,
              crandomnit=args.ratt_nit,
              randatt_forcall=args.ratt_forcall,
              explainp=args.explainp,
              fdupscale=args.dupscale,
              fmupscale=args.mupscale,
              fcc_advanced=args.cc_advanced,
              fvariable_mass=args.variable_mass,
              fwritestats=args.writestats,
              strargs=printsettings(args,2))

if args.steps == 3:
    main_biascorrection(ryyyy=args.ryyyy, ayyyy=args.ayyyy, am=args.am,
               opathA=opath_ATR, 
               opathD=opath_DGN, 
               ipathR=ipath_REF,
               opath=opath_BIA, 
               ofile_base=args.expid, # output
               mode=args.mode,
               maskfile=maskfile,
               maskval=args.maskval,
               verbose=args.verbose,
               veryverbose=args.veryverbose,
               fuseattp=args.bc_useattp,
               bcscale=args.bc_time,
               faggbwtime=args.bc_aggbwtime,
               fdebug=args.debug,
               fwrite_netcdf=args.write_netcdf,
               fwrite_month=args.write_month,
               fwritestats=args.writestats,
               precision=args.precision,
               strargs=printsettings(args,3))
