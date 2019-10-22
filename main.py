#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MAIN SCRIPT MOISTURE AND HEAT DIAGNOSIS
@author: dominik and jessica

To execute interactively: 
> exec(open("./main.py").read())

"""

###########################################################################
#############################    MODULES ##################################

import gzip
import pandas as pd
import numpy as np
import os, fnmatch
import timeit
import netCDF4 as nc4
import sys
import argparse
import time
from datetime import datetime, timedelta
from math import sin,cos,acos,atan,atan2,sqrt
from dateutil.relativedelta import relativedelta
import datetime as datetime

###############################################################################
## ------ USER SETTINGS
###############################################################################

## Paths
# work directory
wpath           = "/kyukon/data/gent/vo/000/gvo00090/vsc42383/tools/flexpart/hamster"
# path to input data
ipath           = "/scratch/gent/vo/000/gvo00090/D2D/data/FLEXPART/era_global/particle-o-matic_t0/gglobal/"
# path for output data
opath           = "/scratch/gent/vo/000/gvo00090/vsc42383/flexpart_data/hamster/01_diagnosis/"


###############################################################################
# ------ END USER SETTINGS
###############################################################################

os.chdir(wpath)

## (1) LOADING FUNCTIONS
exec(open("disclaimer.py").read())
exec(open("constants.py").read())
exec(open("metfunctions.py").read())
exec(open("01_diagnosis.py").read())

## (2) get date, thresholds and flags from command line (job script) 
#      note: this is where we set the default values now. 
args    = read_cmdargs()
verbose = args.verbose

## (3) RUN main script with arguments
main_diagnosis(ryyyy=args.ryyyy, ayyyy=args.ayyyy, am=args.am, 
          ipath=ipath, 
          opath=opath,
          mode=args.mode,
          gres=args.gres,
          sfnam_base=args.expid,
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
          fcc_advanced=args.cc_advanced,
          fvariable_mass=args.variable_mass)
