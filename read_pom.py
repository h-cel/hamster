#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Script to read in FLEXPART trajectories (pom)

@author: jessica

To execute interactively: exec(open("./read_pom.py").read())

"""

#############################  -- MODULES ###########################################

import gzip
import os, fnmatch
import timeit
import sys
from copy import deepcopy
import datetime as datetime
#from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta
import pandas as pd
import numpy as np


#############################  -- SETTINGS ##########################################

## -- SETTINGS
ryear       = "2002"
ayear       = "2002"
amonth      = "2"
mode        = "test"
ipath       = "/scratch/gent/vo/000/gvo00090/D2D/data/FLEXPART/era_global/particle-o-matic_t0/gglobal/"
ifile_base  = ["terabox_NH_AUXTRAJ_", "terabox_SH_AUXTRAJ_"]
mainpath    = ipath+str(ryear)+"/"

# -- DATES
date_bgn        = datetime.datetime.strptime(str(ayear)+"-"+str(amonth).zfill(2)+"-01", "%Y-%m-%d")
date_end        = date_bgn + relativedelta(months=1)
timestep        = datetime.timedelta(hours=6)
date_seq        = []
idate           = date_bgn
while idate < date_end:
    date_seq.append(idate.strftime('%Y%m%d%H'))
    idate   += timestep

print(date_seq)

if mode == "test":
    date_seq = date_seq[0:1]
    print("....\n")
    print("TESTMODE: only processing "+str(date_seq))
    print("....\n")

megatic = timeit.default_timer()
for ix in range(len(date_seq)): 

    # initialize dataar empty for each date (ix)
    dataar  = None
    # loop over ifile_base, concatenating files for the same date 
    for iifile_base in ifile_base:
        # Check if file exists
        ifile   = str(ipath+str(ryear)+"/"+iifile_base+date_seq[ix]+".dat.gz")
        if not os.path.isfile(ifile):
            print(ifile + " does not exist!")
            break
        elif os.path.isfile(ifile):
            # Read file
            print("--------------------------------------------------------------------------------------")
            print("Reading " + ifile)
            ary_dim     = pd.read_table(gzip.open(ifile, 'rb'), sep="\s+", header=None, skiprows=1, nrows=1)
            nparticle   = int(ary_dim[0])
            ntrajstep   = int(ary_dim[1])
            ncolumn     = int(ary_dim[2])
            print("nparticle = ",nparticle, " |  ntrajstep=",ntrajstep,"  | ncolumn=",ncolumn)
            print("--------------------------------------------------------------------------------------") 
            ary_dat     = pd.read_table(gzip.open(ifile, 'rb'), sep="\s+", header=None, skiprows=2)
            datav       = (np.asarray(ary_dat).flatten('C'))
            if dataar is None:
                dataar      = np.reshape(datav, (ntrajstep,nparticle,ncolumn), order='F')
            else:
                dataar      = np.append(dataar, np.reshape(datav, (ntrajstep,nparticle,ncolumn), order='F'), axis=1)
            # flip time axis    (TODO: flip axis depending on forward/backward flag)
            dataar          = dataar[::-1,:,:]

            #######################################################################
            #pID   = array[0,ID,0]   # ID
            #lats  = array[:,ID,2]   # LAT
            #lons  = array[:,ID,1])  # LON
            ######################### for ERA-global ###############################
            ##if dist_on_sphere(lats[0],lons[0],lats[1],lons[1])>1620:
            ##    return(0)
            ########################################################################
            #temp  = array[:,ID,8]   # temperature (K)
            #ztra  = array[:,ID,3]   # height (m)
            ##topo  = array[:,ID,4]  # topography (m) 
            #qv    = array[:,ID,5]   # specific humidity (kg/kg)
            #hmix  = array[:,ID,7]   # ABL height (m)
            #dens  = array[:,ID,6]   # density (kg/m^3)
            #######################################################################


megatoc = timeit.default_timer()
print("\n Runtime: ",str(round(megatoc-megatic, 2)),"seconds")
