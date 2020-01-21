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
    for itime in range(ary_val.shape[0]):  
        with np.errstate(divide='ignore', invalid='ignore'):
            ary_sval[itime,:,:] = ary_val[itime,:,:] * np.nan_to_num( ary_rpart[:,:]/ary_part[itime,:,:] ) 
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


def writenc4D(ofile,fdate_seq,fuptdate_seq,glon,glat,ary_etop,ary_heat,ary_npart):
    if verbose:
        print(" * Writing netcdf output...")
                
    # convert date objects to datetime objects if necessary
    if type(fdate_seq[0]) == datetime.date:
        for idt in range(len(fdate_seq)):
            fdate_seq[idt]    = datetime.datetime(fdate_seq[idt].year, fdate_seq[idt].month, fdate_seq[idt].day)
    if type(fuptdate_seq[0]) == datetime.date:
        for idt in range(len(fuptdate_seq)):
            fuptdate_seq[idt] = datetime.datetime(fuptdate_seq[idt].year, fuptdate_seq[idt].month, fuptdate_seq[idt].day)

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
    atimes              = nc_f.createVariable('arrival-time', 'i4', 'arrival-time')
    utimes              = nc_f.createVariable('uptake-time', 'i4', 'uptake-time')
    latitudes           = nc_f.createVariable('lat', 'f4', 'lat')
    longitudes          = nc_f.createVariable('lon', 'f4', 'lon')
    heats               = nc_f.createVariable('H', 'f4', ('arrival-time','uptake-time','lat','lon'))
    etops               = nc_f.createVariable('E2P', 'f4', ('arrival-time','uptake-time','lat','lon'))
    nparts              = nc_f.createVariable('n_part', 'f4', ('arrival-time','uptake-time','lat','lon'))
    
    # set attributes
    nc_f.description    = "FLEXPART: 02_attribution of advected surface sensible heat and evaporative moisture resulting in precipitation"
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
    nparts.units        = 'int'
    nparts.long_name    = 'number of parcels (mid pos.)'
    
    # write data
    atimes[:]           = nc4.date2num(fdate_seq, atimes.units, atimes.calendar)
    utimes[:]           = nc4.date2num(fuptdate_seq, utimes.units, utimes.calendar)
    longitudes[:]       = glon
    latitudes[:]        = glat
    heats[:]            = ary_heat[:]
    etops[:]            = ary_etop[:]     
    nparts[:]           = ary_npart[:]
        
    # close file
    nc_f.close()
        
    print("\n===============================================================")
    print("\n Successfully written: "+ofile+" !")
    print("\n===============================================================")
