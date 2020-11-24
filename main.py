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
path_ref = content.path_ref
path_orig = content.path_orig
path_diag = content.path_diag
path_attr = content.path_attr
path_bias = content.path_bias
maskfile  = content.maskfile
path_f2t_diag = content.path_f2t_diag
base_f2t_diag = content.base_f2t_diag
path_f2t_traj = content.path_f2t_traj
base_f2t_traj = content.base_f2t_traj

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
if args.waiter:
    waiter  = random.randint(0,30)
    time.sleep(waiter)

# create output directories if they do not exist (in dependency of step)
if args.steps==0 and args.ctraj_len==0 and not os.path.exists(path_f2t_diag):
        os.makedirs(path_f2t_diag)
        os.makedirs(path_f2t_diag+"/"+str(args.ryyyy))
if args.steps==0 and args.ctraj_len>0 and not os.path.exists(path_f2t_traj):
        os.makedirs(path_f2t_traj)
        os.makedirs(path_f2t_traj+"/"+str(args.ryyyy))
if args.steps==1 and not os.path.exists(path_diag):
        os.makedirs(path_diag)
if args.steps==2 and not os.path.exists(path_attr):
        os.makedirs(path_attr)
if args.steps==3 and not os.path.exists(path_bias):
        os.makedirs(path_bias)

## (3) RUN main scripts with arguments
if args.steps ==0:
    if args.ctraj_len==0:
        path_f2t=path_f2t_diag
        base_f2t=base_f2t_diag
    elif args.ctraj_len>0:
        path_f2t=path_f2t_traj
        base_f2t=base_f2t_traj
    main_flex2traj(ryyyy=args.ryyyy, ayyyy=args.ayyyy, am=args.am, ad=args.ad,
                   tml=args.ctraj_len,
                   maskfile=maskfile,
                   maskval=args.maskval,
                   idir=path_orig,
                   odir=path_f2t,
                   fout=base_f2t)

if args.steps == 1:
    main_diagnosis(ryyyy=args.ryyyy, ayyyy=args.ayyyy, am=args.am, ad=args.ad,
              ipath=path_f2t_diag,
              ifile_base=base_f2t_diag,
              opath=path_diag,
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
              ipath=path_f2t_traj,
              ifile_base=base_f2t_traj,
              ipath_f2t=path_orig,
              opath=path_attr,
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
               opathA=path_attr, 
               opathD=path_diag, 
               ipathR=path_ref,
               opath=path_bias, 
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
