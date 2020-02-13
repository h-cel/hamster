#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
NON-METEOROLOGICAL FUNCTIONS USED BY MODULES
01_attribution, 02_diagnosis, 03_bias-correction

@author: jessica and dominik

To execute interactively: 
> exec(open("./hamsterfunctions.py").read())

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
    parser.add_argument('--tdiagnosis', '-dgn', help = "diagnosis method (KAS, SOD, SAJ)",                              type = str,     default = "KAS")
    parser.add_argument('--ctraj_len',  '-len', help = "threshold for maximum allowed trajectory length in days",       type = int,     default = 10)
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
    parser.add_argument('--fallingdry', '-dry', help = "cut off trajectories falling dry (flag)",                       type = str2bol, default = True,     nargs='?')
    parser.add_argument('--memento',    '-mto', help = "keep track of trajectory history (flag)",                       type = str2bol, default = True,     nargs='?')
    parser.add_argument('--variable_mass','-vm',help = "use variable mass (flag)",                                      type = str2bol, default = False,    nargs='?')
    parser.add_argument('--gres',       '-r',   help = "output grid resolution (degrees)",                              type = float,   default = 1)
    parser.add_argument('--ryyyy',      '-ry',  help = "run name (here, YYYY, example: 2002, default: ayyyy)",          type = int,     default = parser.parse_args().ayyyy)
    parser.add_argument('--refdate',    '-rd',  help = "reference date (YYYYMMDDHH)",                                   type = str,     default = str(parser.parse_args().ryyyy)+"123118")
    #print(parser.format_help())
    args = parser.parse_args()  # namespace
    return args

def printsettings(args):
    if args.tdiagnosis in ['KAS']:
        return(str("Diagnosis following Schumacher & Keune (----) with the following settings: " +
        "[[PRECIPITATION]] cprec_dqv = "+str(args.cprec_dqv)+ ", cprec_rh = " +str(args.cprec_rh)+ ", cprec_dtemp = " +str(args.cprec_dtemp) + ", "
        "[[EVAPORATION]] cevap_cc = "+str(args.cevap_cc)+ ", cevap_hgt = " +str(args.cevap_hgt) + ", "
        "[[SENSIBLE HEAT]] cheat_cc = "+str(args.cheat_cc)+ ", cheat_hgt = " +str(args.cheat_hgt)+ ", cheat_dtemp = " +str(args.cheat_dtemp) + ", "
        "[[OTHERS]]: cc_advanced = "+str(args.cc_advanced)+", variable_mass = "+str(args.variable_mass)+ ", mode = "+str(args.mode)))
    if args.tdiagnosis in ['SOD']:
        return(str({"Diagnosis following Sodemann et al. (2008) with the following settings: " +
        "[[PRECIPITATION]] cprec_dqv = "+str(args.cprec_dqv)+ ", cprec_rh = " +str(args.cprec_rh) + ", " +
        "[[EVAPORATION]] cevap_dqv = 0.2, cevap_hgt < 1.5 * mean ABL, " +
        "[[SENSIBLE HEAT]] cheat_dTH = "+str(args.cheat_dtemp)+ ", cheat_hgt < 1.5 * mean ABL, " +
        "[[OTHERS]]: variable_mass = "+str(args.variable_mass)+ ", mode = "+str(args.mode) + " "
        "Sodemann, H., Schwierz, C., & Wernli, H. (2008). Interannual variability of Greenland winter precipitation sources: Lagrangian moisture diagnostic and North Atlantic Oscillation influence. Journal of Geophysical Research: Atmospheres, 113(D3). http://dx.doi.org/10.1029/2007JD008503"}))
    if args.tdiagnosis in ['SAJ']:
        return(str({"Diagnosis following Stohl and James (2004)." +
        "Stohl, A., & James, P. (2004). A Lagrangian analysis of the atmospheric branch of the global water cycle. Part I: Method description, validation, and demonstration for the August 2002 flooding in central Europe. Journal of Hydrometeorology, 5(4), 656-678. https://doi.org/10.1175/1525-7541(2004)005<0656:ALAOTA>2.0.CO;2"}))


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


def checkpbl(ztra,hpbl,maxhgt):
    #if (ztra < max(np.max(hpbl),maxhgt)).all():
    if (ztra < max(hpbl[1],hpbl[0],maxhgt)).all():
        return True
    else:
        return False

def readparcel(parray):
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

def readmidpoint(parray):
    lats    = parray[:,2]                   # latitude
    lons    = parray[:,1]                   # longitude
    lons[lons>180.0] -= 360                 # transform coordinates from [0 ... 360] to [-180 ... 180]
    lat_mid,lon_mid = midpoint_on_sphere(lats[0],lons[0],lats[1],lons[1]) # use own function to calculate midpoint position
    if (lon_mid>179.5): 
        lon_mid -= 360      # now shift all coords that otherwise would be allocated to +180 deg to - 180
    return lat_mid, lon_mid

def readsparcel(parray):
    # reads standard output that is always needed
    qv      = parray[:,5]                   # specific humidity (kg kg-1)
    temp    = parray[:,8]                   # temperature (K)
    hpbl    = parray[:,7]                   # ABL height (m)
    ztra    = parray[:,3]                   # height (m)
    return qv, temp, ztra, hpbl

def readpres(parray):
    temp    = parray[:,8]                   # temperature (K)
    dens    = parray[:,6]                   # density (kg m-3)
    pres    = calc_pres(dens,temp)          # pressure (Pa)
    return pres

def readepottemp(parray):
    temp    = parray[:,8]                   # temperature (K)
    qv      = parray[:,5]                   # specific humidity (kg kg-1)
    dens    = parray[:,6]                   # density (kg m-3)
    pres    = calc_pres(dens,temp)          # pressure (Pa)
    epottemp= calc_pottemp_e(pres, qv, temp)# equivalent potential temperature (K)
    return epottemp

def readpottemp(parray):
    qv      = parray[:,5]                   # specific humidity (kg kg-1)
    dens    = parray[:,6]                   # density (kg m-3)
    temp    = parray[:,8]                   # temperature (K)
    pres    = calc_pres(dens,temp)          # pressure (Pa)
    pottemp = calc_pottemp(pres, qv, temp)  # potential temperature (K)
    return pottemp

def glanceparcel(parray):

    """
    This is used to take a first glance at parcel properties,
    and as of now, is called only for the first 4 timesteps.
 
    Required for arriving air criterion (which goes back 4 steps).
    However, this could be further optimized by only getting the
    last 4 values of hpbl, as for the other variables here,
    only the last / last 2 values are required at first.

    work-in-progress.

    """

    ## parcel information
    ztra    = parray[:,3]                   # height (m) 
    hpbl    = parray[:,7]                   # ABL height (m)
    temp    = parray[:,8]                   # temperature (K)
    qv      = parray[:,5]                   # specific humidity (kg kg-1)
    dens    = parray[:,6]                   # density (kg m-3)
    pres    = calc_pres(dens,temp)          # pressure (Pa)

    return ztra, hpbl, temp, qv, dens, pres


def parceldiff(pvals, meval):
    # difference 
    if meval in ['diff']:
        dpval   = pvals[0] - pvals[1]
    # mean
    if meval in ['mean']:
        dpval   = (pvals[0]+pvals[1])/2#np.mean(pvals)
    # mean
    if meval in ['max']:
        dpval   = np.max(pvals)
    return(dpval)


def trajparceldiff(pvals, meval):
    # difference 
    if meval in ['diff']:
        dpval   = pvals[:-1] - pvals[1:]
    # mean
    if meval in ['mean']:
        dpval   = (pvals[:-1] + pvals[1:])/2
    # max
    if meval in ['max']:
        dpval   = np.amax(pvals[:-1], pvals[1:])
    return(dpval)


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


def PBL_check(z, h, seth, tdiagnosis):

    if tdiagnosis == 'KAS':
        h[h<seth] = seth
        before_inside = np.logical_or( z[1:] < h[:-1], z[1:]  < h[1:])  
        after_inside  = np.logical_or(z[:-1] < h[:-1], z[:-1] < h[1:])  
        change_inside = np.logical_and(before_inside, after_inside)
    elif tdiagnosis == 'SOD':   
        change_inside = ((z[1:]+z[:-1])/2) < 1.5*((h[1:]+h[:-1])/2)
        # NOTE: factor 1.5 is hardcoded      

    return change_inside


def scale_mass(ary_val, ary_part, ary_rpart):
    ary_sval    = np.zeros(shape=(ary_val.shape)) 
    #for itime in range(ary_val.shape[0]):  
    #    with np.errstate(divide='ignore', invalid='ignore'):
    #        ary_sval[itime,:,:] = ary_val[itime,:,:] * np.nan_to_num( ary_rpart[:,:]/ary_part[itime,:,:] ) 
    with np.errstate(divide='ignore', invalid='ignore'):
         ary_sval[:,:] = ary_val[:,:] * np.nan_to_num( ary_rpart[:,:]/ary_part[:,:] ) 
    return(ary_sval)


def linear_discounter(v, min_gain, min_loss):
    """
    ========= inputs =========
    v:        vector extending into past; v[0] is most recent, v[-1] least recent
    min_gain: minimum gain (dv/dt) to be considered
    min_loss: minimum loss to be considered for discounting, must be <= 0
    ==========================
    """

    ## compute dv/dt, prepare dv_disc
    dv = v[:-1] - v[1:]
    dv_disc = np.zeros(shape=len(dv))
    ## get indices of gains and losses
    idx_gains  = np.where(dv >= min_gain)[0]
    idx_losses = np.where(dv <= min_loss)[0]

    for idx in idx_gains:
        ## determine relevant indices (losses only!)
        idx_consider = idx_losses[np.where(idx_losses < idx)]

        ## skip current iteration if no moisture loss occurs between then and t=0
        if len(idx_consider)==0:
            dv_disc[idx] = dv[idx]
            continue

        ## create vector with fractions,
        frc_vec = v[idx_consider]/v[idx_consider+1]
        ## then multiply current dv with product of vector (order irrelevant)
        dv_disc[idx] = np.prod(frc_vec)*dv[idx]
    return(dv_disc)


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


def mgridder(mlon, mlat, pval,
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
    # get grid index
    ind_lat = np.argmin(np.abs(glat-mlat))    # index on grid # ATTN, works only for 1deg grid
    ind_lon = np.argmin(np.abs(glon-mlon))    # index on grid # ATTN, works only for 1deg grid
    # and assign pval to gridcell (init. with 0's)
    gval    = np.zeros(shape=(glat.size, glon.size))       # shape acc. to pre-allocated result array of dim (ntime, glat.size, glon.size)
    gval[ind_lat,ind_lon]    += pval
    return(gval)

def midpindex(parray,glon,glat):
    lats    = parray[:,2]                   # latitude
    lons    = parray[:,1]                   # longitude
    lons[lons>180.0] -= 360                 # transform coordinates from [0 ... 360] to [-180 ... 180]
    mlat,mlon = midpoint_on_sphere2(lats[0],lons[0],lats[1],lons[1]) # use own function to calculate midpoint position
    if (mlon>179.5):
        mlon -= 360      # now shift all coords that otherwise would be allocated to +180 deg to - 180
    # get grid index
    ind_lat = np.argmin(np.abs(glat-mlat))    # index on grid # ATTN, works only for 1deg grid
    ind_lon = np.argmin(np.abs(glon-mlon))    # index on grid # ATTN, works only for 1deg grid
    return ind_lat, ind_lon

def writeemptync(ofile,fdate_seq,glon,glat,strargs):
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
    nc_f.description    = "01 - " + str(strargs)
    today               = datetime.datetime.now()
    nc_f.history        = "Created " + today.strftime("%d/%m/%Y %H:%M:%S") + " using hamster ((c) Dominik Schumacher and Jessica Keune)"
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
    # close file
    nc_f.close()
    print("\n * Created empty file: "+ofile+" of dimension ("+str(len(fdate_seq))+","+str(glat.size)+","+str(glon.size)+") !")


def writenc(ofile,ix,ary_prec,ary_evap,ary_heat,ary_npart):
    if verbose:
        print(" * Writing to netcdf...")

    nc_f = nc4.Dataset(ofile, 'r+')
    nc_f['P'][ix,:,:]       = ary_prec
    nc_f['E'][ix,:,:]       = ary_evap
    nc_f['H'][ix,:,:]       = ary_heat
    nc_f['n_part'][ix,:,:]  = ary_npart
    nc_f.close()

def writeemptync4D(ofile,fdate_seq,fuptdate_seq,glat,glon,strargs):
                
    # delete nc file if it is present (avoiding error message)
    try:
        os.remove(ofile)
    except OSError:
        pass
        
    # create netCDF4 instance
    nc_f = nc4.Dataset(ofile,'w', format='NETCDF4')
    
    # create dimensions 
    nc_f.createDimension('arrival-time', len(fdate_seq))
    nc_f.createDimension('uptake-time', len(fuptdate_seq))
    nc_f.createDimension('lat', glat.size)
    nc_f.createDimension('lon', glon.size)
    
    # create variables
    atimes              = nc_f.createVariable('arrival-time', 'f4', 'arrival-time')
    utimes              = nc_f.createVariable('uptake-time', 'f4', 'uptake-time')
    latitudes           = nc_f.createVariable('lat', 'f4', 'lat')
    longitudes          = nc_f.createVariable('lon', 'f4', 'lon')
    heats               = nc_f.createVariable('H', 'f4', ('arrival-time','uptake-time','lat','lon'))
    etops               = nc_f.createVariable('E2P', 'f4', ('arrival-time','uptake-time','lat','lon'))
    
    # set attributes
    nc_f.description    = "02 - " + str(strargs)
    today               = datetime.datetime.now()
    nc_f.history        = "Created " + today.strftime("%d/%m/%Y %H:%M:%S") + " using hamster ((c) Jessica Keune and Dominik Schumacher)"
    atimes.units        = 'hours since 1900-01-01 00:00:00'
    atimes.calendar     = 'Standard' 
    utimes.units        = 'hours since 1900-01-01 00:00:00'
    utimes.calendar     = 'Standard' 
    latitudes.units     = 'degrees_north'
    longitudes.units    = 'degrees_east'
    heats.units         = 'W m-2'
    heats.long_name	= 'surface sensible heat flux'
    etops.units         = 'mm'
    etops.long_name	= 'evaporation resulting in precipitation'
  
    # write data
    atimes[:]           = nc4.date2num(fdate_seq, atimes.units, atimes.calendar)[:]
    utimes[:]           = nc4.date2num(fuptdate_seq, utimes.units, utimes.calendar)[:]
    longitudes[:]       = glon[:]
    latitudes[:]        = glat[:]
        
    # close file
    nc_f.close()

    print("\n * Created empty file: "+ofile+" of dimension ("+str(len(fdate_seq))+","+str(len(fuptdate_seq))+","+str(glat.size)+","+str(glon.size)+") !")

        
def writenc4D(ofile,ix,ary_etop,ary_heat):
    if verbose:
        print(" * Writing to netcdf...")

    nc_f = nc4.Dataset(ofile, 'r+')
    nc_f['E2P'][ix,:,:,:]     = ary_etop[:,:,:]
    nc_f['H'][ix,:,:,:]       = ary_heat[:,:,:]
    nc_f.close()
