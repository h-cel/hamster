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


def main_flex2traj(ryyyy, ayyyy, am, ad, tml, maskfile, maskval,
        idir, odir, fout, verbose):

    ###--- MISC ---################################################################
    logo =""" 
        Hello, user. 
        
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       %  __ _           ____  _              _  %
      %  / _| | _____  _|___ \| |_ _ __ __ _ (_)  %  
     %  | |_| |/ _ \ \/ / __) | __| '__/ _` || |   %
     %  |  _| |  __/>  < / __/| |_| | | (_| || |   %
      % |_| |_|\___/_/\_\_____|\__|_|  \__,_|/ |  %
       %                                    |__/ %                                     
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    """
    ###############################################################################
    ###--- SETUP ---###############################################################
        
    #******************************************************************************
    ## UNCOMMENT line ---> variable not saved
    selvars=np.asarray([
    0, # pid        |        [ ALWAYS ] * 0
    1, # x          |        [ ALWAYS ] * 1
    2, # y          |        [ ALWAYS ] * 2
    3, # z          |        [ ALWAYS ] * 3
    #4, # itramem    |        [ NEVER  ]
    5, # oro        |        [ OPTIONAL ] # only needed for dry static energy   
    #6, # pv         |        [ OPTIONAL ]
    7, # qq         |        [ ALWAYS ] * 4
    8, # rho        |        [ ALWAYS ] * 5
    9, # hmix       |        [ ALWAYS ] * 6
    #10,# tropo      |        [ OPTIONAL ] * 7  # needed for droughtpropag
    11,# temp       |        [ ALWAYS ] * 8
    #12,# mass       |        [ NEVER! ]
    ])
    thevars = np.asarray(["pid","x","y","z","itramem","oro","pv",
                          "qv","rho","hmix","tropo","temp","mass"])
    #******************************************************************************
    
    # last day of month
    ed   = int(calendar.monthrange(ayyyy, am)[1])

    dt_h = 6 # hardcoded, as further edits would be necessary if this was changed!
    time_bgn = datetime.datetime(year=ayyyy, month=am, day=ad, hour=6)
    # add 6 hours to handle end of month in same way as any other period
    time_end = datetime.datetime(year=ayyyy, month=am, day=ed, hour=18) + datetime.timedelta(hours=dt_h)
    # convert trajectory length from day to dt_h (!=6); +2 needed ;)
    ntraj = tml*(24//dt_h) + 2 
    
    ###############################################################################
    ###--- MAIN ---################################################################
    
    if verbose: print(logo)

    ##---0.) pepare directories
    outdir = odir+"/"+str(ryyyy)
    if not os.path.exists(outdir): # could use isdir too
        os.makedirs(outdir)
    
    ##---1.) load netCDF mask
    if maskfile is None or maskval==-999:
        mask = mlat = mlon = None
    else:
        mask, mlat, mlon = maskgrabber(maskfile)
        
    ##---2.) create datetime object (covering arrival period + trajectory length)
    fulltime_str = f2t_timelord(ntraj_d=tml, dt_h=dt_h,
                               tbgn=time_bgn, tend=time_end)
    
    #---3.) handle first step
    if verbose: print("\n---- Reading files to begin constructing trajectories ...\n")
    data, trajs = f2t_establisher(partdir=idir+'/'+str(ryyyy), selvars=selvars,
                                 time_str=fulltime_str[:ntraj], ryyyy=ryyyy,                      
                                 mask=mask, maskval=maskval, mlat=mlat, mlon=mlon,
                                 outdir=outdir, fout=fout,
                                 verbose=verbose)
    
    ##---4.) continue with next steps
    if verbose: print("\n\n---- Adding more files ... ")
    for ii in range(1, len(fulltime_str)-ntraj+1): # CAUTION: INDEXING from 1!
        data, trajs = f2t_ascender(data=data, partdir=idir+'/'+str(ryyyy), selvars=selvars,
                                   time_str=fulltime_str[ii:ntraj+ii], ryyyy=ryyyy,
                                   mask=mask, maskval=maskval, mlat=mlat, mlon=mlon,
                                   outdir=outdir, fout=fout,
                                   verbose=verbose)
   
    ##---5.) done
    if verbose: 
        print("\n\n---- Done! \n     Files with base '"+fout+"' written to:\n    ",odir+'/'+str(ryyyy))
        print("     Dimensions: nstep x nparcel x nvar\n     Var order: ", end='')
        print(*thevars[selvars].tolist(), sep=', ')
        print("\n     All done!")
