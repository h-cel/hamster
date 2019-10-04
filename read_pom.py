#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Script to read in FLEXPART trajectories (pom)

@author: jessica

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

ryear="2002"
ayear="2002"
amonth="2"
ipath="/scratch/gent/vo/000/gvo00090/D2D/data/FLEXPART/era_global/particle-o-matic_t0/gglobal/"
opath="/scratch/gent/vo/000/gvo00090/vsc42383/flexpart_data/hamster/01_diagnosis/"
sfnam_base="FXvG_r"

mainpath  = ipath+str(ryear)+"/"

# DATES
date_bgn      = datetime.datetime.strptime(str(ayear)+"-"+str(amonth).zfill(2)+"-01", "%Y-%m-%d")
date_end      = date_bgn + relativedelta(months=1)
timestep      = datetime.timedelta(hours=6)
date_seq      = []
idate         = date_bgn
while idate < date_end:
	date_seq.append(idate.strftime('%Y%m%d%H'))
	idate += timestep

print(date_seq)

for ix in range(len(date_seq)): 

	# Check if file exists
	ifile	= ipath+str(ryear)+"/"+"terabox_NH_AUXTRAJ_"+date_seq[ix]+".dat.gz"
	if not os.path.isfile(ifile)
		print(ifile + " does not exist!")
		break

	# Read file
	print("Reading " + ifile)
        ary_dim 	= pd.read_table(gzip.open(ifile, 'rb'), sep="\s+", header=None, skiprows=1, nrows=1)
	nparticle 	= int(ary_dim[0])
        ntrajstep 	= int(ary_dim[1])
        ncolumn   	= int(ary_dim[2])
        print("-----------------------------------------------------------------------------")
        print("NH: nparticle = ",nparticle, " |  ntrajstep=",ntrajstep,"  | ncolumn=",ncolumn)
        print("-----------------------------------------------------------------------------")	
        ary_dat 	= pd.read_table(gzip.open(ifile, 'rb'), sep="\s+", header=None, skiprows=2)
        ary 		= np.fliplr(np.asarray(ary_pd)) # flips axis 1
	# data 		= (np.asarray(ary_dat).flatten('C'))
	# dataar 	= np.reshape(data, (ntrajstep,nparticle,ncolumn), order='F')
	
	# Check if file exists
	ifile	= ipath+str(ryear)+"/"+"terabox_SH_AUXTRAJ_"+date_seq[ix]+".dat.gz"
	if not os.path.isfile(ifile)
		print(ifile + " does not exist!")
		break
	
	# Read file
	print("Reading " + ifile)
        ary_dim 	= pd.read_table(gzip.open(ifile, 'rb'), sep="\s+", header=None, skiprows=1, nrows=1)
	nparticle 	= int(ary_dim[0])
        ntrajstep 	= int(ary_dim[1])
        ncolumn   	= int(ary_dim[2])
        print("-----------------------------------------------------------------------------")
        print("SH: nparticle = ",nparticle, " |  ntrajstep=",ntrajstep,"  | ncolumn=",ncolumn)
        print("-----------------------------------------------------------------------------")	
        ary_dat 	= pd.read_table(gzip.open(ifile, 'rb'), sep="\s+", header=None, skiprows=2)
        ary 		= npappend( ary, np.fliplr(np.asarray(ary_pd)), axis=0) # flips axis 1 and append to ary


