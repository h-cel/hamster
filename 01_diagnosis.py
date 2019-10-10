#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MAIN FUNCTIONS FOR 01_diagnosis
"""


def readpom(idate,     # run year
            ipath,      # input data path
            ifile_base):# loop over ifile_base filenames for each date

    """
    INPUT
        - idate :       date as string [YYYYMMDDHH]
        - ipath :       path where input files are located
        - ifile_base :  base filename(s); loop over filenames possible
    ACTION
        reads trajectories into 3D array of dimension (ntrajlength x nparticles x nvars),
        flipping time axis (HARDCODED: from backward to 'forward', i.e. to 1 = now, 0 = previous)
        and concatenates all information of files, 
            - [ifile_base][idate].dat.gz
        e.g. terabox_NH_AUXTRAJ_2002080100.dat.gz and terabox_SH_AUXTRAJ_2002080100.dat.gz
        to one array (nparticle = SUM ( nparticle[*] ) for all files of ifile_base)
    RETURNS
        - dataar :      data array of dimension (ntrajlength x nparticles x nvars)
    """

    dataar  = None
    # loop over ifile_base, concatenating files for the same date 
    for iifile_base in ifile_base:
        # Check if file exists
        ifile   = str(ipath+"/"+iifile_base+idate+".dat.gz")
        if not os.path.isfile(ifile):
            print(ifile + " does not exist!")
            break
        elif os.path.isfile(ifile):
            # Read file
            if verbose:
                print("--------------------------------------------------------------------------------------")
                print("Reading " + ifile)
            ary_dim     = pd.read_table(gzip.open(ifile, 'rb'), sep="\s+", header=None, skiprows=1, nrows=1)
            nparticle   = int(ary_dim[0])
            ntrajstep   = int(ary_dim[1])
            nvars       = int(ary_dim[2])
            if verbose:
                print("nparticle = ",nparticle, " |  ntrajstep=",ntrajstep,"  | =",nvars)
                print("--------------------------------------------------------------------------------------")
            ary_dat     = pd.read_table(gzip.open(ifile, 'rb'), sep="\s+", header=None, skiprows=2)
            datav       = (np.asarray(ary_dat).flatten('C'))
            if dataar is None:
                dataar      = np.reshape(datav, (ntrajstep,nparticle,nvars), order='F')
            else:
                dataar      = np.append(dataar, np.reshape(datav, (ntrajstep,nparticle,nvars), order='F'), axis=1)
            # flip time axis    (TODO: flip axis depending on forward/backward flag)
            dataar          = dataar[::-1,:,:]

    return(dataar)

def readparcel(parray):
    
    ## parcel information
    lats    = parray[:,2]                   # latitude
    lons    = parray[:,1]                   # longitude
    lons[lons>180.0] -= 360                 # transform coordinates from [0 ... 360] to [-180 ... 180]
    temp    = parray[:,8]                   # temperature (K)
    ztra    = parray[:,3]                   # height (m)
    #topo   = parray[:,4]                   # topography (m) 
    qv      = parray[:,5]                   # specific humidity (kg kg-1)
    hpbl    = parray[:,7]                   # ABL height (m)
    dens    = parray[:,6]                   # density (kg m-3)
    pres    = calc_pres(dens,temp)          # pressure (Pa)
    pottemp = calc_pottemp(pres, qv, temp)  # potential temperature (K)
    epottemp= calc_pottemp_e(pres, qv, temp)# equivalent potential temperature (K)

    return lons, lats, temp, ztra, qv, hpbl, dens, pres, pottemp, epottemp 

def parceldiff(pvals, meval):
    # difference 
    if meval in ['diff']:
        dpval   = pvals[0] - pvals[1]
    # mean
    if meval in ['mean']:
        dpval   = np.mean(pvals)
    # mean
    if meval in ['max']:
        dpval   = np.max(pvals)
    return(dpval)

def gridder(plon, plat, pval,
            glat, glon):
    """
    INPUT
        - plon, plat: parcel longitutde and latitude
        - glon, glat: grid longitude and latitutde
        - pval      : parcel value to be assigned to grid
    ACTION
        1. calculated midpoint of two coordinates
        2. assigns val to gridcell corresponding to midpoint
    RETURN
        - array of dimension (glat.size x glon.size) with 0's and one value assigned
    """
    # 1. calculate midpoint
    lat_mid,lon_mid = midpoint_on_sphere(plat[0],plon[0],plat[1],plon[1]) # use own function to calculate midpoint position
    if (lon_mid>179.5): lon_mid -= 360    # now shift all coords that otherwise would be allocated to +180 deg to - 180
    # 2. get grid index
    ind_lat = np.argmin(np.abs(glat-lat_mid))    # index on grid # ATTN, works only for 1deg grid
    ind_lon = np.argmin(np.abs(glon-lon_mid))    # index on grid # ATTN, works only for 1deg grid
    # and assign pval to gridcell (init. with 0's)
    gval    = np.zeros(shape=(glat.size, glon.size))       # shape acc. to pre-allocated result array of dim (ntime, glat.size, glon.size)
    gval[ind_lat,ind_lon]    += pval
    return(gval)

def default_thresholds(P_dq_min):
    if P_dq_min == None:
        #if verbose:
        #    print("\n--- INFO: P_dq_min is calculated based on d(pottemp)-threshold!")
        dummy_dq = 0.2 # this choice doesn't matter too much...
        P_dq_min = -(1/(calc_pottemp_e(PREF, (5+dummy_dq)/1e3, TREF+15) - 
                       calc_pottemp_e(PREF, 5/1e3, TREF+15)))*dummy_dq/1e3
        #print("P_dq_min = ", 1e3*P_dq_min, "g/kg")
    elif P_dq_min > 0:
        raise SystemExit("------ FATAL ERROR: P_dq_min should be negative (and in kg/kg)!")
    return P_dq_min

def convertunits(ary_val, garea, var):
    """
    INPUT
        - aryval
    ACTION
        - calculates grid cell values 
    RETURN
        - returns P and E as mm
        - returns H as W m-2
    """
    if var in ['P','E']:
        return(PMASS*ary_val/(1e6*garea))
    if var in ['H']:
        return(PMASS*ary_val*CPD/(1e6*garea*6*3600))


############################################################################
#############################    SETTINGS ##################################

def readNmore(
           ryyyy, ayyyy, am,
           ipath, opath,
           mode,
           gres,
           sfnam_base,
           dTH_thresh=0., # used for E,H,P (if P_dq_min==None)
           f_dqsdT=0.7, f_dTdqs=0.7, # for H, E diagnosis (lower = more strict)
           hmax_E=0, hmax_H=0, # set min ABLh, disabled if 0 
           P_dq_min=None, P_dT_thresh=0, P_RHmin=80, # P settings
           fwrite_netcdf=True,ftimethis=True,fcc_advanced=False):

    """
    comments
    
    - more sophisticated ABL criteria as in versionD are highly beneficial to
      avoid E-overestimation over deserts; for H, it does not seem to change
      as much.
    - with the current configuration, there are only 4 parameters:
        
        dTH_thresh = 1. (Kelvin),
        f_dqdst == f_dTdqs,
        P_dT_thresh = 0. (Kelvin), # not a good idea to increase this a lot    
        P_RHmin=80 (%) 
        
        thus, the previosuly introduced dz-Parameter could return,
        with the advantage of being used both for E,H & P 
        (as of now, dz>0 is used for P anyways).
    """

    ## construct precise input and storage paths
    mainpath  = ipath+str(ryyyy)+"/"
    sfilename = str(sfnam_base)+str(ryyyy)[-2:]+"_"+str(ayyyy)+"-"+str(am).zfill(2)+".nc"

    ########### LOG W/IN PYTHON SCRIPT by redirecting output #############
    
    if verbose:
        print("\n============================================================================================================")
        print(os.system("figlet -f bubble hamster"))
        print("\n PROCESSING: \t", 	ayyyy, "-", str(am).zfill(2))
        print("\n INPUT PATH: \t", 	ipath)
        print("\n OUTPUT FILE: \t", 	opath+sfilename)
        print("\n============================================================================================================")
        print("\n============================================================================================================")
        
    ## start timer
    if ftimethis:
        megatic = timeit.default_timer()
    
    glon, glat, garea = makegrid(resolution=gres)

    ## -- DATES
    date_bgn        = datetime.datetime.strptime(str(ayyyy)+"-"+str(am).zfill(2)+"-01", "%Y-%m-%d")
    date_end        = date_bgn + relativedelta(months=1)
    timestep        = datetime.timedelta(hours=6)
    date_seq        = []
    fdate_seq       = []
    idate           = date_bgn
    while idate < date_end:
        date_seq.append(idate.strftime('%Y%m%d%H'))
        fdate_seq.append(idate)
        idate   += timestep
    ntime           = len(date_seq)

    # TESTMODE
    if mode == "test":
        ntime       = 1
        date_seq    = date_seq[0:ntime]
        fdate_seq   = fdate_seq[0:ntime]
        print("....\n")
        print("TESTMODE: only processing "+str(date_seq) + " and looping over 1000 parcels")
        print("....\n")

    ## pre-allocate arrays
    ary_heat     = np.zeros(shape=(ntime,glat.size,glon.size))
    ary_evap     = np.zeros(shape=(ntime,glat.size,glon.size))
    ary_prec     = np.zeros(shape=(ntime,glat.size,glon.size))
    ary_npart    = np.zeros(shape=(ntime,glat.size,glon.size))

    # set some default thresholds
    P_dq_min    = default_thresholds(P_dq_min) 

    for ix in range(ntime):
        print("Processing "+str(fdate_seq[ix]))
        ## 1) read in all files associated with data --> ary is of dimension (ntrajlen x nparticles x nvars)
        ary = readpom( idate    = date_seq[ix], 
                       ipath    = "/scratch/gent/vo/000/gvo00090/D2D/data/FLEXPART/era_global/particle-o-matic_t0/gglobal/"+str(ryyyy), 
                       ifile_base = ["terabox_NH_AUXTRAJ_", "terabox_SH_AUXTRAJ_"])
        nparticle   = ary.shape[1]
        print("TOTAL: " + str(date_seq[ix]) + " has " + str(nparticle) + " parcels")

        #bar = Bar('Processing', suffix='%(percent)d%%', fill="*")
        if mode == "test":
            ntot    = range(1000)
        else:
            ntot    = range(nparticle)

        ## 2) diagnose P, E, H and npart per grid cell
        for i in ntot:

            ## - 2.1) read parcel information
            lons, lats, temp, ztra, qv, hpbl, dens, pres, pottemp, epottemp = readparcel(ary[:,i,:])

            ## - 2.2) parcel changes / criteria
            dq          = parceldiff(qv, 'diff') 
            hpbl_max    = parceldiff(hpbl, 'max')
            dT          = parceldiff(temp, 'diff')
            dTH         = parceldiff(pottemp, 'diff')
            dTHe        = parceldiff(epottemp, 'diff')
            dz          = parceldiff(ztra, 'diff')

            ## - 2.3) diagnose fluxes

            ## (a) number of parcels
            ary_npart[ix,:,:] += gridder(plon=lons, plat=lats, pval=int(1), glon=glon, glat=glat)

            ## (b) precipitation
            if ( dq < P_dq_min and 
                 q2rh(qv[0], pres[0], temp[0]) > P_RHmin  and
                 q2rh(qv[1], pres[1], temp[1]) > P_RHmin ):
                ary_prec[ix,:,:] += gridder(plon=lons, plat=lats, pval=dq, glon=glon, glat=glat)

            ## (c) evaporation
            if fcc_advanced:
                if ( ztra[0] <  max(hmax_E, hpbl_max)  and
                     ztra[1] <  max(hmax_E, hpbl_max)  and
                     (dTHe - dTH) > dTH_thresh and
                     ( (dT > 0 and dT       < f_dTdqs * (dq) * dTdqs(p_hPa=pres[1]/1e2, q_kgkg=qv[1])) or
                       (dT < 0 and abs(dTH) < f_dTdqs * (dq) * dTdqs(p_hPa=pres[1]/1e2, q_kgkg=qv[1]))
                     )
                   ):
                    ary_evap[ix,:,:] += gridder(plon=lons, plat=lats, pval=dq, glon=glon, glat=glat)
            else:
                if ( ztra[0] <  max(hmax_E, hpbl_max)  and
                     ztra[1] <  max(hmax_E, hpbl_max)  and
                     (dTHe - dTH) > dTH_thresh and
                     abs(dTH) < f_dTdqs * (dq) * dTdqs(p_hPa=pres[1]/1e2, q_kgkg=qv[1]) ):
                    ary_evap[ix,:,:] += gridder(plon=lons, plat=lats, pval=dq, glon=glon, glat=glat)


            ## (d) sensible heat
            if fcc_advanced:
                if ( ztra[0] <  max(hmax_H, hpbl_max) and 
                     ztra[1] <  max(hmax_H, hpbl_max) and 
                     (dTH > dTH_thresh) and 
                     ( (dT > 0 and abs(dq) < f_dqsdT * (dT)  * dqsdT(p_hPa=pres[1]/1e2, T_degC=temp[1]-TREF)) or
                       (dT < 0 and abs(dq) < f_dqsdT * (dTH) * dqsdT(p_hPa=pres[1]/1e2, T_degC=temp[1]-TREF))
                     )
                   ):
                    ary_heat[ix,:,:] += gridder(plon=lons, plat=lats, pval=dTH, glon=glon, glat=glat) 
            else:
                if ( ztra[0] <  max(hmax_H, hpbl_max) and 
                     ztra[1] <  max(hmax_H, hpbl_max) and 
                     (dTH > dTH_thresh) and 
                     abs(dq) < f_dqsdT * (dTH) * dqsdT(p_hPa=pres[1]/1e2, T_degC=temp[1]-TREF) ):
                    ary_heat[ix,:,:] += gridder(plon=lons, plat=lats, pval=dTH, glon=glon, glat=glat) 


        # Convert units
        ary_prec[ix,:,:] = convertunits(ary_prec[ix,:,:], garea, "P")
        ary_evap[ix,:,:] = convertunits(ary_evap[ix,:,:], garea, "E")
        ary_heat[ix,:,:] = convertunits(ary_heat[ix,:,:], garea, "H")

    if ftimethis:
        megatoc = timeit.default_timer()
        print("\n=======    main loop completed, total runtime so far: ",str(round(megatoc-megatic, 2)),"seconds")
    
    ###########################################################################    
    
    if fwrite_netcdf:
            
        ### delete nc file if it is present (avoiding error message)
        try:
            os.remove(opath+sfilename)
        except OSError:
            pass
        
        ### create netCDF4 instance
        nc_f = nc4.Dataset(opath+sfilename,'w', format='NETCDF4')
        
        ### create dimensions ###
        nc_f.createDimension('time', ntime)
        nc_f.createDimension('lat', glat.size)
        nc_f.createDimension('lon', glon.size)
    
        ### create variables
        times       = nc_f.createVariable('time', 'i4', 'time')
        latitudes   = nc_f.createVariable('lat', 'f4', 'lat')
        longitudes  = nc_f.createVariable('lon', 'f4', 'lon')
        heats       = nc_f.createVariable('H', 'f4', ('time','lat','lon'))
        evaps       = nc_f.createVariable('E', 'f4', ('time','lat','lon'))
        precs       = nc_f.createVariable('P', 'f4', ('time','lat','lon'))
        nparts      = nc_f.createVariable('n_part', 'f4', ('time','lat','lon'))
    
        ### set attributes
        nc_f.description    = "FLEXPART: 01_diagnosis of upward land surface fluxes and precipitation"
        times.units         = 'hours since 1900-01-01 00:00:00'
        times.calendar      = 'Standard' # do NOT use gregorian here!
        latitudes.units     = 'degrees_north'
        longitudes.units    = 'degrees_east'
        heats.units         = 'W m-2'
        heats.long_name	    = 'surface sensible heat flux'
        evaps.units         = 'mm'
        evaps.long_name	    = 'evaporation'
        precs.units         = 'mm'
        precs.long_name	    = 'precipitation'
        nparts.units        = 'int'
        nparts.long_name    = 'number of parcels (mid pos.)'
    
        ### write data
        times[:]            = nc4.date2num(fdate_seq, times.units, times.calendar)
        longitudes[:]       = glon
        latitudes[:]        = glat
        
        heats[:]            = ary_heat[:]
        evaps[:]            = ary_evap[:]
        precs[:]            = ary_prec[:]     
        nparts[:]           = ary_npart[:]
        
        ### close file
        nc_f.close()
        
        print("\n===============================================================")
        print("\n======  Successfully written: "+opath+sfilename+ " !")
        print("\n===============================================================")
        
        
