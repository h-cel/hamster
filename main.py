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

###########################################################################
##--- PATHS
###########################################################################

## determine working directory
wpath = os.getcwd()

## load input and output paths & input file name base(s)
with open(wpath+"/paths.txt") as f: 
    content = imp.load_source('','',f) # load like a python module
    ipath_DGN = content.ipath_DGN # input path (01_diagnosis)
    ibase_DGN = content.ibase_DGN # input file name base(s)
    opath_DGN = content.opath_DGN # output path
    ipath_ATR = content.ipath_ATR # as above (for 02_attribution)
    ibase_ATR = content.ibase_ATR 
    opath_ATR = content.opath_ATR 

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
exec(open("hamsterfunctions.py").read())

## (2) get date, thresholds and flags from command line (job script) 
#      note: this is where we set the default values now. 
args    = read_cmdargs()
verbose = args.verbose
print(printsettings(args))

## (3) RUN main script with arguments
main_diagnosis(ryyyy=args.ryyyy, ayyyy=args.ayyyy, am=args.am,
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
          refdate=args.refdate,
          fwrite_netcdf=args.write_netcdf,
          ftimethis=args.timethis,
          fcc_advanced=args.cc_advanced,
          fvariable_mass=args.variable_mass,
          strargs=printsettings(args))

main_attribution(ryyyy=args.ryyyy, ayyyy=args.ayyyy, am=args.am, 
          ipath=ipath_ATR,
          ifile_base=ibase_ATR,
          opath=opath_ATR,
          ofile_base=args.expid,
          mode=args.mode,
          gres=args.gres,
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
          refdate=args.refdate,
          fwrite_netcdf=args.write_netcdf,
          ftimethis=args.timethis, 
          fdry=args.fallingdry,
          fmemento=args.memento,
          fcc_advanced=args.cc_advanced,
          fvariable_mass=args.variable_mass,
          strargs=printsettings(args))
