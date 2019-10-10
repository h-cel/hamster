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

## Time period
ryyyy           = 2002
ayyyy           = 2002
am              = 5
## Experiment ID (choose a letter or short name)
expID           = "FXvH_r"
## mode (test/oper)
mode            = "test"

## DIAGNOSIS SETTINGS
tdTH            = 1.0       # used for E,H,P (if cprec_dqv==None)
cheat_cc         = 0.7       # for CC criterion of H, E diagnosis (lower = more strict)
f_dqTds         = 0.7       # for H, E diagnosis (lower = more strict)
cevap_hgt          = 0         # up to which height should E be considered? max(cevap_hgt, BLh_max)
cheat_hgt          = 0         # up to which height should P be considered? max(cheat_hgt, BLh_max)
cprec_dqv        = None      # 
cprec_dtemp     = 0         #
cprec_rh         = 80        # 

# Optional flags
write_netcdf  = True      # write netcdf output
timethis      = True      # check runtime of diagnoser & gridder
scale_mass    = False     # scale mass with number of particles 
cc_advanced    = False     # use advanced Clausius-Clapeyron criteria
verbose         = True      # use as global variable

###############################################################################
# ------ END USER SETTINGS
###############################################################################

os.chdir(wpath)

## (1) LOADING FUNCTIONS
exec(open("constants.py").read())
exec(open("metfunctions.py").read())
exec(open("01_diagnosis.py").read())

## (2) RUN
readNmore(ryyyy=ryyyy, ayyyy=ayyyy, am=am, 
          ipath=ipath, 
          opath=opath,
          mode=mode,
          gres=1,
          sfnam_base=expID,
          cheat_dtemp=tdTH,
          cheat_cc=0.7, 
          cevap_cc=0.7,
          cevap_hgt=0, 
          cheat_hgt=0,
          cprec_dqv=None, 
          #cprec_dtemp=0, 
          cprec_rh=80,
          fwrite_netcdf=write_netcdf,
          ftimethis=timethis, 
          fcc_advanced=cc_advanced)
