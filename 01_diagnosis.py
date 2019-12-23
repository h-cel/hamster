#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MAIN FUNCTIONS FOR 01_diagnosis
"""

def str2bol(v):
    """
    ACTION: converts (almost) any string to boolean False/True
    NOTE:   needed for boolean interpretation in parser.add_argument from parsearg in read_cmdargs
    """
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def read_cmdargs():
    """
    ACTION: read dates, thresholds and flags from command line
    RETURN: 'args' contains all
    DEP:    uses argparse
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--ayyyy',      '-ay',  help = "analysis year (YYYY)",                                          type = int,     default = 2002)
    parser.add_argument('--am',         '-am',  help = "analysis month (M)",                                            type = int,     default = 1)
    parser.add_argument('--mode',       '-m',   help = "mode (test,oper)",                                              type = str,     default = "oper")
    parser.add_argument('--expid',      '-id',  help = "experiment ID (string, example versionA)",                      type = str,     default = "FXv")
    parser.add_argument('--diagnosis',  '-dgn', help = "diagnosis method (KAS, SOD, SAJ)",                              type = str,     default = "KAS")
    parser.add_argument('--cprec_dqv',  '-cpq', help = "threshold for detection of P based on delta(qv)",               type = float,   default = 0)
    parser.add_argument('--cprec_rh',   '-cpr', help = "threshold for detection of P based on RH",                      type = float,   default = 80)
    parser.add_argument('--cprec_dtemp','-cpt', help = "threshold for detection of P based on delta(T)",                type = float,   default = 0)
    parser.add_argument('--cevap_cc',   '-cec', help = "threshold for detection of E based on CC criterion",            type = float,   default = 0.7)
    parser.add_argument('--cevap_hgt',  '-ceh', help = "threshold for detection of E using a maximum height",           type = float,   default = 0)
    parser.add_argument('--cheat_cc',   '-chc', help = "threshold for detection of H based on CC criterion",            type = float,   default = 0.7)
    parser.add_argument('--cheat_hgt',  '-chh', help = "threshold for detection of H using a maximum height",           type = float,   default = 0)
    parser.add_argument('--cheat_dtemp','-cht', help = "threshold for detection of H using a minimum delta(T)",         type = float,   default = 0)
    parser.add_argument('--cc_advanced','-cc',  help = "use advanced CC criterion (flag)",                              type = str2bol, default = False,    nargs='?')
    parser.add_argument('--timethis',   '-t',   help = "time the main loop (flag)",                                     type = str2bol, default = False,    nargs='?')
    parser.add_argument('--write_netcdf','-o',  help = "write netcdf output (flag)",                                    type = str2bol, default = True,     nargs='?')
    parser.add_argument('--verbose',    '-v',   help = "verbose output (flag)",                                         type = str2bol, default = True,     nargs='?')
    parser.add_argument('--variable_mass','-vm',help = "use variable mass (flag)",                                      type = str2bol, default = False,    nargs='?')
    parser.add_argument('--gres',       '-r',   help = "output grid resolution (degrees)",                              type = float,   default = 1)
    parser.add_argument('--ryyyy',      '-ry',  help = "run name (here, YYYY, example: 2002, default: ayyyy)",          type = int,     default = parser.parse_args().ayyyy)
    parser.add_argument('--refdate',    '-rd',  help = "reference date (YYYYMMDDHH)",                                   type = str,     default = str(parser.parse_args().ryyyy)+"123118")
    #print(parser.format_help())
    args = parser.parse_args()  # namespace
    return args

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
                print(" Reading " + ifile)
            ary_dim     = pd.read_table(gzip.open(ifile, 'rb'), sep="\s+", header=None, skiprows=1, nrows=1)
            nparticle   = int(ary_dim[0])
            ntrajstep   = int(ary_dim[1])
            nvars       = int(ary_dim[2])
            if verbose:
                print("\t nparticle = ",nparticle, " |  ntrajstep=",ntrajstep,"  | =",nvars)
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

def default_thresholds(cprec_dqv):
    if cprec_dqv == None:
        #if verbose:
        #    print("\n--- INFO: cprec_dqv is calculated based on d(pottemp)-threshold!")
        dummy_dq = 0.2 # this choice doesn't matter too much...
        cprec_dqv = -(1/(calc_pottemp_e(PREF, (5+dummy_dq)/1e3, TREF+15) - 
                       calc_pottemp_e(PREF, 5/1e3, TREF+15)))*dummy_dq/1e3
        #print("cprec_dqv = ", 1e3*cprec_dqv, "g/kg")
    elif cprec_dqv > 0:
        raise SystemExit("------ FATAL ERROR: cprec_dqv should be negative (and in kg/kg)!")
    return cprec_dqv

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

def get_refnpart(refdate, ryyyy, glon, glat):
    """
    INPUT
        - refdate [YYYYMMDDHH] :    reference date (str) used for counting midpoint parcels and scaling
        - glon, glat :              reference grid coordinates      
    ACTION
        - calculates the reference distribution of parcels using the midpoint of parcels at refdate
        - NOTE that this is run specific and needs to be adjusted if FLEXPART runs are setup differently
    DEPEND
        - uses numpy and functions readpom, gridder
    RETURN
        - npart (nlat x nlon) at refdate
    """
    if verbose:
        #print("Reference number of particles: \t" + str(nparticle))
        print(" * Getting reference distribution...")

    ary_npart   = np.zeros(shape=(glat.size,glon.size))
    ary         = readpom( idate    = refdate,
                           ipath    = "/scratch/gent/vo/000/gvo00090/D2D/data/FLEXPART/era_global/particle-o-matic_t0/gglobal/"+str(ryyyy),
                           ifile_base = ["terabox_NH_AUXTRAJ_", "terabox_SH_AUXTRAJ_"])
    nparticle   = ary.shape[1]
    for i in range(nparticle):
        lons, lats, _, _, _, _, _, _, _, _ = readparcel(ary[:,i,:])
        ary_npart[:,:] += gridder(plon=lons, plat=lats, pval=int(1), glon=glon, glat=glat)

    return ary_npart

def scale_mass(ary_val, ary_part, ary_rpart):
    ary_sval    = np.zeros(shape=(ary_val.shape)) 
    for itime in range(ary_val.shape[0]):  
        with np.errstate(divide='ignore', invalid='ignore'):
            ary_sval[itime,:,:] = ary_val[itime,:,:] * np.nan_to_num( ary_rpart[:,:]/ary_part[itime,:,:] ) 
    return(ary_sval)

def writenc(ofile,fdate_seq,glon,glat,ary_prec,ary_evap,ary_heat,ary_npart):
    if verbose:
        print(" * Writing netcdf output...")
                
    # delete nc file if it is present (avoiding error message)
    try:
        os.remove(ofile)
    except OSError:
        pass
        
    # create netCDF4 instance
    nc_f = nc4.Dataset(ofile,'w', format='NETCDF4')
    
    ### create dimensions ###
    nc_f.createDimension('time', len(fdate_seq))
    nc_f.createDimension('lat', glat.size)
    nc_f.createDimension('lon', glon.size)
    
    # create variables
    times               = nc_f.createVariable('time', 'i4', 'time')
    latitudes           = nc_f.createVariable('lat', 'f4', 'lat')
    longitudes          = nc_f.createVariable('lon', 'f4', 'lon')
    heats               = nc_f.createVariable('H', 'f4', ('time','lat','lon'))
    evaps               = nc_f.createVariable('E', 'f4', ('time','lat','lon'))
    precs               = nc_f.createVariable('P', 'f4', ('time','lat','lon'))
    nparts              = nc_f.createVariable('n_part', 'f4', ('time','lat','lon'))
    
    # set attributes
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
    
    # write data
    times[:]            = nc4.date2num(fdate_seq, times.units, times.calendar)
    longitudes[:]       = glon
    latitudes[:]        = glat
    heats[:]            = ary_heat[:]
    evaps[:]            = ary_evap[:]
    precs[:]            = ary_prec[:]     
    nparts[:]           = ary_npart[:]
        
    # close file
    nc_f.close()
        
    print("\n===============================================================")
    print("\n Successfully written: "+ofile+" !")
    print("\n===============================================================")
        
        
############################################################################
#############################    SETTINGS ##################################

def main_diagnosis(
           ryyyy, ayyyy, am,
           ipath, opath,
           mode,
           gres,
           sfnam_base,
           diagnosis,
           cheat_dtemp, # used for E,H,P (if cprec_dqv==None)
           cheat_cc, cevap_cc, # for H, E diagnosis (lower = more strict)
           cevap_hgt, cheat_hgt, # set min ABLh, disabled if 0 
           cprec_dqv, cprec_dtemp, cprec_rh,
           refdate,
           fwrite_netcdf,ftimethis,fcc_advanced,fvariable_mass):

    """
    comments
    
    - more sophisticated ABL criteria as in versionD are highly beneficial to
      avoid E-overestimation over deserts; for H, it does not seem to change
      as much.
    - with the current configuration, there are only 4 parameters:
        
        cheat_dtemp = 1. (Kelvin),
        f_dqdst == cevap_cc,
        cprec_dtemp = 0. (Kelvin), # not a good idea to increase this a lot    
        cprec_rh=80 (%) 
        
        thus, the previosuly introduced dz-Parameter could return,
        with the advantage of being used both for E,H & P 
        (as of now, dz>0 is used for P anyways).
    """

    ## construct precise input and storage paths
    mainpath  = ipath+str(ryyyy)+"/"
    sfilename = str(sfnam_base)+str(ryyyy)[-2:]+"_"+str(ayyyy)+"-"+str(am).zfill(2)+".nc"
    ofile     = opath+sfilename

    ########### LOG W/IN PYTHON SCRIPT by redirecting output #############
    
    if verbose:
        disclaimer()
        print("\n SETTINGS :")
        print("\n PROCESSING: \t", 	ayyyy, "-", str(am).zfill(2))
        print("\n input path: \t", 	ipath)
        print("\n============================================================================================================")
        print(" ! using variable mass: \t" +str(fvariable_mass) )
        if fvariable_mass:
            print(" \t ! reference date for number of particles: \t" +str(refdate) )
        print(" ! writing netcdf output: \t" +str(fwrite_netcdf) )
        if fwrite_netcdf:
            print(" \t ! with grid resolution:: \t", str(gres) )
            print(" \t ! output file: \t", opath+sfilename)
        print(" ! using internal timer: \t" +str(ftimethis) )
        print(" ! mode: \t" +str(mode))
        print(" ! DIAGNOSIS SETTINGS")
        print(" \t ! HEAT: ")
        print(" \t \t  dTH > " + str(cheat_dtemp) )
        print(" \t \t  abs(dqv) < "+str(cheat_cc)+" * (dTH) * ...")
        print(" \t \t  ztra[0] <  max("+str(cheat_hgt)+", hpbl_max) ")
        print(" \t \t  ztra[1] <  max("+str(cheat_hgt)+", hpbl_max) ")
        print(" \t ! + using advanced CC criteria: \t" +str(fcc_advanced) )
        print(" \t ! EVAPORATION: ")
        print(" \t \t  abs(dTH) < "+str(cevap_cc)+" * (dqv) * ...")
        print(" \t \t  ztra[0] <  max("+str(cevap_hgt)+", hpbl_max) ")
        print(" \t \t  ztra[1] <  max("+str(cevap_hgt)+", hpbl_max) ")
        print(" \t ! + using advanced CC criteria: \t" +str(fcc_advanced) )
        print(" \t ! PRECIPITATION: ")
        print(" \t \t  dqv < "+str(cprec_dqv) )
        print(" \t \t  rh[0] > "+str(cprec_rh) )
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

    ## pre-allocate arrays
    ary_heat     = np.zeros(shape=(ntime,glat.size,glon.size))
    ary_evap     = np.zeros(shape=(ntime,glat.size,glon.size))
    ary_prec     = np.zeros(shape=(ntime,glat.size,glon.size))
    ary_npart    = np.zeros(shape=(ntime,glat.size,glon.size))

    # set some default thresholds
    cprec_dqv    = default_thresholds(cprec_dqv) 
    # read in reference distribution of parcels
    if fvariable_mass:
        ary_rnpart   = get_refnpart(refdate=refdate, ryyyy=ryyyy, glon=glon, glat=glat)
    
    ## loop over time to read in files
    if verbose:
        print("\n=== \t Start main program...\n")
    for ix in range(ntime):
        if verbose:
                print("--------------------------------------------------------------------------------------")
        print("Processing "+str(fdate_seq[ix]))
        ## 1) read in all files associated with data --> ary is of dimension (ntrajlen x nparticles x nvars)
        ary = readpom( idate    = date_seq[ix], 
                       ipath    = "/scratch/gent/vo/000/gvo00090/D2D/data/FLEXPART/era_global/particle-o-matic_t0/gglobal/"+str(ryyyy), 
                       ifile_base = ["terabox_NH_AUXTRAJ_", "terabox_SH_AUXTRAJ_"])
        nparticle   = ary.shape[1]
        if verbose:
            print(" TOTAL: " + str(date_seq[ix]) + " has " + str(nparticle) + " parcels")

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
            if diagnosis == 'KAS':
                hpbl_max    = parceldiff(hpbl, 'max')
                dTH         = parceldiff(pottemp, 'diff')
                dTHe        = parceldiff(epottemp, 'diff')
                #dz          = parceldiff(ztra, 'diff')
                if fcc_advanced:
                    dT          = parceldiff(temp, 'diff')
            elif diagnosis == 'SOD':
                hpbl_avg    = parceldiff(hpbl, 'mean')
                dTH         = parceldiff(pottemp, 'diff')

            ## - 2.3) diagnose fluxes

            ## (a) number of parcels
            ary_npart[ix,:,:] += gridder(plon=lons, plat=lats, pval=int(1), glon=glon, glat=glat)

            ##  - 2.3)-KAS: Keune and Schumacher
            if diagnosis == 'KAS':

                ## (b) precipitation
                if ( dq < cprec_dqv and 
                     q2rh(qv[0], pres[0], temp[0]) > cprec_rh  and
                     q2rh(qv[1], pres[1], temp[1]) > cprec_rh ):
                    ary_prec[ix,:,:] += gridder(plon=lons, plat=lats, pval=dq, glon=glon, glat=glat)

                ## (c) evaporation
                if fcc_advanced:
                    if ( ztra[0] <  max(cevap_hgt, hpbl_max)  and
                         ztra[1] <  max(cevap_hgt, hpbl_max)  and
                         (dTHe - dTH) > cheat_dtemp and
                         ( (dT > 0 and dT       < cevap_cc * (dq) * dTdqs(p_hPa=pres[1]/1e2, q_kgkg=qv[1])) or
                           (dT < 0 and abs(dTH) < cevap_cc * (dq) * dTdqs(p_hPa=pres[1]/1e2, q_kgkg=qv[1]))
                         )
                       ):
                        ary_evap[ix,:,:] += gridder(plon=lons, plat=lats, pval=dq, glon=glon, glat=glat)
                else:
                    if ( ztra[0] <  max(cevap_hgt, hpbl_max)  and
                         ztra[1] <  max(cevap_hgt, hpbl_max)  and
                         (dTHe - dTH) > cheat_dtemp and
                         abs(dTH) < cevap_cc * (dq) * dTdqs(p_hPa=pres[1]/1e2, q_kgkg=qv[1]) ):
                        ary_evap[ix,:,:] += gridder(plon=lons, plat=lats, pval=dq, glon=glon, glat=glat)

                ## (d) sensible heat
                if fcc_advanced:
                    if ( ztra[0] <  max(cheat_hgt, hpbl_max) and 
                         ztra[1] <  max(cheat_hgt, hpbl_max) and 
                         (dTH > cheat_dtemp) and 
                         ( (dT > 0 and abs(dq) < cheat_cc * (dT)  * dqsdT(p_hPa=pres[1]/1e2, T_degC=temp[1]-TREF)) or
                           (dT < 0 and abs(dq) < cheat_cc * (dTH) * dqsdT(p_hPa=pres[1]/1e2, T_degC=temp[1]-TREF))
                         )
                       ):
                        ary_heat[ix,:,:] += gridder(plon=lons, plat=lats, pval=dTH, glon=glon, glat=glat) 
                else:
                    if ( ztra[0] <  max(cheat_hgt, hpbl_max) and 
                         ztra[1] <  max(cheat_hgt, hpbl_max) and 
                         (dTH > cheat_dtemp) and 
                         abs(dq) < cheat_cc * (dTH) * dqsdT(p_hPa=pres[1]/1e2, T_degC=temp[1]-TREF) ):
                        ary_heat[ix,:,:] += gridder(plon=lons, plat=lats, pval=dTH, glon=glon, glat=glat) 

            ##  - 2.3)-SOD: Sodemann et al., 2008
            elif diagnosis == 'SOD':
         
                ## (b) precipitation
                if ( dq < 0 and 
                     q2rh((qv[0]+qv[1])/2, (pres[0]+pres[1])/2, (temp[0]+temp[1])/2) > 80 ):
                    ary_prec[ix,:,:] += gridder(plon=lons, plat=lats, pval=dq, glon=glon, glat=glat)

                ## (c) evaporation
                if ( dq > 0.0002 and 
                     (ztra[0]+ztra[1])/2 <  hpbl_avg
                   ):
                    ary_evap[ix,:,:] += gridder(plon=lons, plat=lats, pval=dq, glon=glon, glat=glat)
    
                ## (d) sensible heat (not used originally; analogous to evaporation)
                if ( (dTH > cheat_dtemp) and 
                    (ztra[0]+ztra[1])/2 <  hpbl_avg
                   ):
                    ary_heat[ix,:,:] += gridder(plon=lons, plat=lats, pval=dTH, glon=glon, glat=glat)

            ##  - 2.3)-SAJ: Stohl and James, 2004
            elif diagnosis == 'SAJ':

                ## (b) precipitation
                if ( dq < 0 ):
                    ary_prec[ix,:,:] += gridder(plon=lons, plat=lats, pval=dq, glon=glon, glat=glat)
    
                ## (c) evaporation
                if ( dq > 0 ):
                    ary_evap[ix,:,:] += gridder(plon=lons, plat=lats, pval=dq, glon=glon, glat=glat)

        # Convert units
        if verbose:
            print(" * Converting units...")
        if diagnosis == 'KAS' or diagnosis == 'SOD':
            ary_prec[ix,:,:] = convertunits(ary_prec[ix,:,:], garea, "P")
            ary_evap[ix,:,:] = convertunits(ary_evap[ix,:,:], garea, "E")
            ary_heat[ix,:,:] = convertunits(ary_heat[ix,:,:], garea, "H")
        elif diagnosis =='SAJ':
            # first calculate column sums, assign to E or P according to sign
            colsum = ary_prec[ix,:,:] + ary_evap[ix,:,:]
            colsum_pos = np.zeros_like(colsum)
            colsum_neg = np.zeros_like(colsum)
            colsum_pos[np.where(colsum>0)] = colsum[np.where(colsum>0)]
            colsum_neg[np.where(colsum<0)] = colsum[np.where(colsum<0)]
            ary_evap[ix,:,:] = convertunits(colsum_pos, garea, "E")
            ary_prec[ix,:,:] = convertunits(colsum_neg, garea, "P")

    # Scale with parcel mass
    if fvariable_mass:
        if verbose: 
            print(" * Applying variable mass...")
        ary_prec         = scale_mass(ary_prec, ary_npart, ary_rnpart)
        ary_evap         = scale_mass(ary_evap, ary_npart, ary_rnpart)
        ary_heat         = scale_mass(ary_heat, ary_npart, ary_rnpart)

    if ftimethis:
        megatoc = timeit.default_timer()
        if verbose:
            print("\n=== \t End main program (total runtime so far: ",str(round(megatoc-megatic, 2)),"seconds) \n")
    
    if fwrite_netcdf:
        writenc(ofile,fdate_seq,glon,glat,ary_prec,ary_evap,ary_heat,ary_npart)
