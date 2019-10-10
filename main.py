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
tdTH            = 1.0       # used for E,H,P (if P_dq_min==None)
f_dqsdT         = 0.7       # for CC criterion of H, E diagnosis (lower = more strict)
f_dqTds         = 0.7       # for H, E diagnosis (lower = more strict)
hmax_E          = 0         # up to which height should E be considered? max(hmax_E, BLh_max)
hmax_H          = 0         # up to which height should P be considered? max(hmax_H, BLh_max)
P_dq_min        = None      # 
P_dT_thresh     = 0         #
P_RHmin         = 80        # 

# Optional flags
ffwrite_netcdf  = True      # write netcdf output
fftimethis      = True      # check runtime of diagnoser & gridder
ffscale_mass    = False     # scale mass with number of particles 
fcc_advanced    = False     # use advanced Clausius-Clapeyron criteria
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
          dTH_thresh=tdTH,
          f_dqsdT=0.7, 
          f_dTdqs=0.7,
          hmax_E=0, 
          hmax_H=0,
          P_dq_min=None, 
          P_dT_thresh=0, 
          P_RHmin=80,
          write_netcdf=ffwrite_netcdf,
          timethis=fftimethis, fcc_advanced=fcc_advanced)
