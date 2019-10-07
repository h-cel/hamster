#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MAIN FUNCTIONS FOR 01_diagnosis
"""

def diagnoser(parray,
              dTH_thresh=tdTH, 
              f_dqsdT=1.0, 
              f_dTdqs=1.0, 
              hmax_E=0, 
              hmax_H=0, 
              P_dq_min=None, 
              P_dT_thresh=0, 
              P_RHmin=80):
    
    """
    INPUTS
    ACTION
    RETURNS
    """
        
    ###### CONSTANTS #######
    R_specific = 287.057                    # round(8.3144598 / (28.9645/1e3), 3)
    DALR       = (-9.8/1e3)                 # dry adiabatic lapse rate

    ####################################################################### 
    ###### PARTICLE INFORMATION ########
    lats    = parray[:,2]                   # latitude
    lons    = parray[:,1]                   # longitude
    # skip this parcel if unreasonable jump in data
    if dist_on_sphere(lats[0],lons[0],lats[1],lons[1])>1620:
        return(0)
    temp    = parray[:,8]                   # temperature (K)
    ztra    = parray[:,3]                   # height (m)
    #topo   = parray[:,4]                   # topography (m) 
    qv      = parray[:,5]                   # specific humidity (kg kg-1)
    hmix    = parray[:,7]                   # ABL height (m)
    dens    = parray[:,6]                   # density (kg m-3)
    ## calculate everything else that is needed at least once from here
    dq      = qv[0] - qv[1]                 # humidity change (kg kg-1)
    BLh_max = np.max(hmix[:2])              # max. boundary layer height (m)
    pres    = dens*R_specific*temp          # pressure (Pa)
    potT    = calc_theta(pres, qv, temp)    # potential temperature (K)
    eqvpotT = calc_theta_e(pres, qv, temp)  # equivalent potential temperature (K)
    dTH     = potT[0] - potT[1]             # potential temperature difference
    dTHe    = (eqvpotT[0]-eqvpotT[1])       # equiv. pot. temperature difference (K)
    dz      = ztra[0] - ztra[1]             # height difference (rise/descent) (m)
    #dz = (ztra[0]+topo[0]) - (ztra[1]+topo[1]) # WHY WOULD THIS BE CORRECT?
    ####################################################################### 
    
    counter = 0 # not elegant, but it works
    
    ########### PRECIPITATION ############
    if ( 
          dq < P_dq_min and
          q2rh(qv[0], dens[0]*R_specific*temp[0], temp[0])>P_RHmin and 
          q2rh(qv[1], dens[1]*R_specific*temp[1], temp[1])>P_RHmin 
          ):
            counter += 1
    elif ( 
        dq < P_dq_min and
        dz > 0 # not really a parameter, but rather: is the parcel ascending?
        ): 
        ## moist adiabatic lapse rate, T @ Lifting Condensation Level
        T_LCL, MALR = moist_ascender(p_Pa=pres[1], q_kgkg=qv[1], T_K=temp[1])
        ## calculate whether LCL was reached during ascent
        dz_reachLCL = (T_LCL-temp[1])/DALR
        if dz > dz_reachLCL:
            dz_rem = dz - dz_reachLCL    # remaining ascent
            T_moi = T_LCL  + dz_rem*MALR # hypoth. T if partially moist adiabatic
            T_dry = temp[1]+ dz*DALR     # hypoth. T if purely dry adibatic
            ## check if temperatures are somewhat consistent
            if ( # yes, yet another if
                # check if temperature change closer to what we expect in case
                # of condensation occurring during ascent
                abs(temp[0]-T_moi) < abs(temp[0]-T_dry) and 
                (temp[0]-T_moi) < P_dT_thresh 
                # above: if air is 'too warm', ABL 'detrainment' in large-scale
                # subsidence areas more likely (results in pseudo-precip)
                ): 
                counter += 1

    ########### EVAPORATION & SENSIBLE HEAT ############
    if (
        (ztra[0] < hmix[0]) and
        (ztra[1] < hmix[1]) and
        (dTH > dTH_thresh) and 
        ((dTHe - dTH) > dTH_thresh) and         
        abs(dq) < f_dqsdT*(dTH)*dqsdT(p_hPa=dens[1]*R_specific*temp[1]/1e2, T_degC=temp[1]-273.15) and
        abs(dTH) < f_dTdqs*(dq)*dTdqs(p_hPa=dens[1]*R_specific*temp[1]/1e2, q_kgkg=qv[1])
        ):
        ### EVAP-HEAT ###
        counter += 6
        return(counter)
        
    ########### SENSIBLE HEAT ############
    if (   
        (ztra[0] <  max(hmax_H, BLh_max)) and 
        (ztra[1] <  max(hmax_H, BLh_max)) and 
        (dTH > dTH_thresh) and 
        abs(dq) < f_dqsdT*(dTH)*dqsdT(p_hPa=dens[1]*R_specific*temp[1]/1e2, T_degC=temp[1]-273.15)
        ):
        ### HEAT ###
        counter += 4
    
    ########### EVAPORATION ############
    if ( 
        (ztra[0] <  max(hmax_E, BLh_max)) and 
        (ztra[1] <  max(hmax_E, BLh_max)) and
        ((dTHe - dTH) > dTH_thresh) and
         abs(dTH) < f_dTdqs*(dq)*dTdqs(p_hPa=dens[1]*R_specific*temp[1]/1e2, q_kgkg=qv[1])
       ):
        ### EVAP ###
        counter += 2

    return(counter)
    
def freadpom(idate,     # run year
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
            print("--------------------------------------------------------------------------------------")
            print("Reading " + ifile)
            ary_dim     = pd.read_table(gzip.open(ifile, 'rb'), sep="\s+", header=None, skiprows=1, nrows=1)
            nparticle   = int(ary_dim[0])
            ntrajstep   = int(ary_dim[1])
            nvars       = int(ary_dim[2])
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

def gridder(parray, 
            code,
            fileID, 
            glat, glon, gd_area, 
            logger, 
            heat, evap, prec):
    
    """
    INPUT 
        -
    ACTION
    RETURNS
    """

    ###### CONSTANTS #######
    R_specific = 287.057                        # round(8.3144598 / (28.9645/1e3), 3)
    cpd        = 1005.7
    mass_part  = 2548740090557.712              # 5.09198e18/1997842
    
    #######################################################################
    # Save midpoint position 
    lats        = parray[:,2]                   # lat
    lons        = parray[:,1]                   # lon
    lons[lons>180.0] -= 360                     # enforce proper coordinates for midpoint function [-180.0 .. 180.0]
    lat_mid,lon_mid = midpoint_on_sphere(lats[0],lons[0],lats[1],lons[1]) # use own function to calculate midpoint position
    if (lon_mid>179.5): lon_mid -= 360    # now shift all coords that otherwise would be allocated to +180 deg to - 180
    lat_ix = np.argmin(np.abs(glat-lat_mid))    # index on grid # ATTN, works only for 1deg grid
    lon_ix = np.argmin(np.abs(glon-lon_mid))    # index on grid # ATTN, works only for 1deg grid
    
    ###################################################################
    logger[fileID,lat_ix,lon_ix] += 1 ## store npart per pixel
    ###################################################################
    
    if code > 0:
        qv      = parray[:,5]                   # specific humidity (kg/kg)
        dq      = qv[0] - qv[1]                 # specific humidity change 
    
        if code>=4:
            temp  = parray[:,8]                 # temperature (K)
            dens  = parray[:,6]                 # density (kg/m^3)
            pres  = dens*R_specific*temp        # pressure (Pa)
            theta = calc_theta(pres, qv, temp)  # potential temperature 
            dT    = (theta[0] - theta[1])*cpd   # <<-------------------------------- assuming dry air; check this!
        if code >= 4:
            heat[fileID,lat_ix,lon_ix] += mass_part*dT/(1e6*gd_area[lat_ix]*6*3600) # Wm-2
        if code in [2,3,6,7]: # I know, so amateur-like. sorry.    
            evap[fileID,lat_ix,lon_ix] += mass_part*dq/(1e6*gd_area[lat_ix]) # mm
        if code%2-1==0:
            prec[fileID,lat_ix,lon_ix] += mass_part*dq/(1e6*gd_area[lat_ix]) # mm
            
    #-- DONE!
    return(logger, heat, evap, prec)


############################################################################
#############################    SETTINGS ##################################

def readNmore(
           ryyyy, ayyyy, am,
           ipath, opath,
           mode,
           sfnam_base,           
           dTH_thresh=tdTH,          # used for E,H,P (if P_dq_min==None)
           f_dqsdT=0.7, f_dTdqs=0.7, # for H, E diagnosis (lower = more strict)
           hmax_E=0, hmax_H=0, # set min ABLh, disabled if 0 
           P_dq_min=None, P_dT_thresh=0, P_RHmin=80, # P settings
           verbose=True,write_netcdf=True,timethis=True):

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
        
    ##########################    EXPERIMENTAL    #############################
    if P_dq_min == None:
        #if verbose:
        #    print("\n--- INFO: P_dq_min is calculated based on d(theta)-threshold!")
        dummy_dq = 0.2 # this choice doesn't matter too much...
        P_dq_min = -(1/(calc_theta_e(1013.25e2, (5+dummy_dq)/1e3, 273.15+15) - 
                       calc_theta_e(1013.25e2, 5/1e3, 273.15+15)))*dummy_dq/1e3
        #print("P_dq_min = ", 1e3*P_dq_min, "g/kg")
    elif P_dq_min > 0:
        raise SystemExit("------ FATAL ERROR: P_dq_min should be negative (and in kg/kg)!")
    ###########################################################################
        
    ###########################################################################
    #############################     SETUP      ##############################
    
    ## start timer
    if timethis:
        megatic = timeit.default_timer()
    
    ## prepare grid
    resolution = 1. # in degrees
    glat = np.arange(-90,90+resolution,resolution) # lats from -90 to + 90 deg
    glon = np.arange(-180,180,resolution)
    nlat = glat.size
    nlon = glon.size

    ###########################################################################
    ###########################################################################
    
    ## grab gridded areas prior to loop to save some CPU time
    gd_area = gridded_area_exact(glat, res=1.0, R=6371)
    
    ###########################################################################

    # -- DATES
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
    ary_heat     = np.zeros(shape=(ntime,nlat,nlon))
    ary_evap     = np.zeros(shape=(ntime,nlat,nlon))
    ary_prec     = np.zeros(shape=(ntime,nlat,nlon))
    npart_log    = np.zeros(shape=(ntime,nlat,nlon), dtype=int)

    for ix in range(ntime):
        print("Processing "+str(fdate_seq[ix]))
        ## 1) read in all files associated with data --> ary is of dimension (ntrajlen x nparticles x nvars)
        ary = freadpom(idate    = date_seq[ix], 
                       ipath    = "/scratch/gent/vo/000/gvo00090/D2D/data/FLEXPART/era_global/particle-o-matic_t0/gglobal/"+str(ryyyy), 
                       ifile_base = ["terabox_NH_AUXTRAJ_", "terabox_SH_AUXTRAJ_"])
        nparticle   = ary.shape[1]
        print("TOTAL: " + str(date_seq[ix]) + " has " + str(nparticle) + " parcels")

        bar = Bar('Processing', suffix='%(percent)d%%', fill="*")
        if mode == "test":
            ntot    = range(1000)
        else:
            ntot    = range(nparticle)
        ## 2) diagnose P, E, H and npart per grid cell
        for i in ntot:
        ## - 2.1) diagnose fluxes and position
            # diagnosis
            diagcodes = diagnoser(parray=ary[:,i,:],
                                  dTH_thresh=dTH_thresh, 
                                  f_dqsdT=f_dqsdT, 
                                  f_dTdqs=f_dTdqs,
                                  hmax_E=hmax_E, 
                                  hmax_H=hmax_H, 
                                  P_dq_min=P_dq_min, 
                                  P_dT_thresh=P_dT_thresh, 
                                  P_RHmin=P_RHmin)
        ## - 2.2) grid to predefined grid
            # write to arrays: npart, H, E, P
            npart_log, ary_heat, ary_evap, ary_prec = gridder(  parray=ary[:,i,:], 
                                                                code=diagcodes,
                                                                fileID=ix, 
                                                                glat=glat, glon=glon, gd_area=gd_area, 
                                                                logger=npart_log, 
                                                                heat=ary_heat, evap=ary_evap, prec=ary_prec)
            ## here comes some multiprocessing to diagnose all particles
            #if __name__ == '__main__':
            #    num_cores = multiprocessing.cpu_count()
            #    diagcodes = Parallel(n_jobs=num_cores)(delayed(diagnoser)(
            #                         ID=i,npart=nparticle,array=ary,
            #                         dTH_thresh=dTH_thresh, f_dqsdT=f_dqsdT, f_dTdqs=f_dTdqs, 
            #                         hmax_E=hmax_E, hmax_H=hmax_H, 
            #                         P_dq_min=P_dq_min, P_dT_thresh=P_dT_thresh, P_RHmin=P_RHmin
            #                         ) for i in range(nparticle))
            ### proceed to call gridding function w/o parallelization
            #for i in range(nparticle):
            #    npart_log, ary_heat, ary_evap, ary_prec = gridder(
            #             ID=i, npart=nparticle, array=ary, code=diagcodes[i],
            #             fileID=ix, glat=glat, glon=glon, gd_area=gd_area, 
            #             logger=npart_log, heat=ary_heat, evap=ary_evap, prec=ary_prec)
            #        
            # update progress bar
            bar.next()
        # finish progress bar
        bar.finish()

    ###########################################################################
    
    #if verbose:
    #    print("\n===========================  CONFIGURATION  =======================")
    #    print("\nayyyy=", ayyyy)
    #    print("\nresolution=", resolution)
    #    print("nlat=", nlat)
    #    print("nlon=", nlon)
                 
    if timethis:
        megatoc = timeit.default_timer()
        print("\n=======    main loop completed, total runtime so far: ",str(round(megatoc-megatic, 2)),"seconds")
    
    ###########################################################################    
    
    if write_netcdf:
            
        ### delete nc file if it is present (avoiding error message)
        try:
            os.remove(opath+sfilename)
        except OSError:
            pass
        
        ### create netCDF4 instance
        nc_f = nc4.Dataset(opath+sfilename,'w', format='NETCDF4')
        
        ### create dimensions ###
        nc_f.createDimension('time', ntime)
        nc_f.createDimension('lat', nlat)
        nc_f.createDimension('lon', nlon)
    
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
    
        nparts[:]           = npart_log[:]
        
        ### close file
        nc_f.close()
        
        print("\n===============================================================")
        print("\n======  Successfully written: "+opath+sfilename+ " !")
        print("\n===============================================================")
        
        
