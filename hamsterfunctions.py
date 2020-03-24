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
    parser.add_argument('--steps',      '-st',  help = "steps performed (1: 01_diagnosis, ..., 4: all)",                type = int,     default = 4)
    parser.add_argument('--ayyyy',      '-ay',  help = "analysis year (YYYY)",                                          type = int,     default = 2002)
    parser.add_argument('--am',         '-am',  help = "analysis month (M)",                                            type = int,     default = 1)
    parser.add_argument('--ad',         '-ad',  help = "analysis day (D)",                                              type = int,     default = 1)
    parser.add_argument('--mode',       '-m',   help = "mode (test,oper)",                                              type = str,     default = "oper")
    parser.add_argument('--expid',      '-id',  help = "experiment ID (string, example versionA)",                      type = str,     default = "FXv")
    parser.add_argument('--tdiagnosis', '-dgn', help = "diagnosis method (KAS, SOD/SOD2, SAJ)",                         type = str,     default = "KAS")
    parser.add_argument('--maskval',    '-mv',  help = "use <value> from maskfile for masking",                         type = int,     default = 1)
    parser.add_argument('--ctraj_len',  '-len', help = "threshold for maximum allowed trajectory length in days",       type = int,     default = 10)
    parser.add_argument('--cprec_dqv',  '-cpq', help = "threshold for detection of P based on delta(qv)",               type = float,   default = 0)
    parser.add_argument('--cprec_rh',   '-cpr', help = "threshold for detection of P based on RH",                      type = float,   default = 80)
    parser.add_argument('--cprec_dtemp','-cpt', help = "threshold for detection of P based on delta(T)",                type = float,   default = 0)
    parser.add_argument('--cevap_cc',   '-cec', help = "threshold for detection of E based on CC criterion",            type = float,   default = 0.7)
    parser.add_argument('--cevap_hgt',  '-ceh', help = "threshold for detection of E using a maximum height",           type = float,   default = 0)
    parser.add_argument('--cheat_cc',   '-chc', help = "threshold for detection of H based on CC criterion",            type = float,   default = 0.7)
    parser.add_argument('--cheat_hgt',  '-chh', help = "threshold for detection of H using a maximum height",           type = float,   default = 0)
    parser.add_argument('--cheat_dtemp','-cht', help = "threshold for detection of H using a minimum delta(T)",         type = float,   default = 0)
    parser.add_argument('--cpbl_strict','-pbl', help = "1: both within max, 2: one within max, 3: not used",            type = int,     default = 1)
    parser.add_argument('--fjumps',     '-fj',  help = "filter out jumps (flag)",                                       type = str2bol, default = True,    nargs='?')
    parser.add_argument('--fjumpsfull', '-fjf', help = "filter out jumps for full trajectory length (flag)",            type = str2bol, default = False,   nargs='?')
    parser.add_argument('--cjumps',     '-cj',  help = "threshold to filter for jumps [km]",                            type = int,     default = 2000)
    parser.add_argument('--cc_advanced','-cc',  help = "use advanced CC criterion (flag)",                              type = str2bol, default = False,    nargs='?')
    parser.add_argument('--timethis',   '-t',   help = "time the main loop (flag)",                                     type = str2bol, default = False,    nargs='?')
    parser.add_argument('--write_netcdf','-o',  help = "write netcdf output (flag)",                                    type = str2bol, default = True,     nargs='?')
    parser.add_argument('--verbose',    '-v',   help = "verbose output (flag)",                                         type = str2bol, default = True,     nargs='?')
    parser.add_argument('--fallingdry', '-dry', help = "cut off trajectories falling dry (flag)",                       type = str2bol, default = True,     nargs='?')
    parser.add_argument('--memento',    '-mto', help = "keep track of trajectory history (flag)",                       type = str2bol, default = True,     nargs='?')
    parser.add_argument('--variable_mass','-vm',help = "use variable mass (flag)",                                      type = str2bol, default = False,    nargs='?')
    parser.add_argument('--setnegzero', '-sz',  help = "set negative ERA-I E & H fluxes to zero (flag)",                type = str2bol, default = True,     nargs='?')
    parser.add_argument('--debug',      '-d',   help = "debugging option (flag)",                                       type = str2bol, default = False,    nargs='?')
    parser.add_argument('--frankenstein','-fr', help = "PRELIM: pom mask.dat used for P-scaling (flag)",                type = str2bol, default = True,    nargs='?')
    parser.add_argument('--gres',       '-r',   help = "output grid resolution (degrees)",                              type = float,   default = 1)
    parser.add_argument('--ryyyy',      '-ry',  help = "run name (here, YYYY, example: 2002, default: ayyyy)",          type = int,     default = parser.parse_args().ayyyy)
    parser.add_argument('--refdate',    '-rd',  help = "reference date (YYYYMMDDHH)",                                   type = str,     default = str(parser.parse_args().ryyyy)+"123118")
    #print(parser.format_help())
    args = parser.parse_args()  # namespace
    return args

def printsettings(args,step):
    ## 01_DIAGNOSIS
    if step == 1 and args.tdiagnosis in ['KAS']:
        return(str("Diagnosis following Schumacher & Keune (----) with the following settings: " +
        "[[PRECIPITATION]] cprec_dqv = "+str(args.cprec_dqv)+ ", cprec_rh = " +str(args.cprec_rh)+ ", cprec_dtemp = " +str(args.cprec_dtemp) + ", "
        "[[EVAPORATION]] cevap_cc = "+str(args.cevap_cc)+ ", cevap_hgt = " +str(args.cevap_hgt) + ", "
        "[[SENSIBLE HEAT]] cheat_cc = "+str(args.cheat_cc)+ ", cheat_hgt = " +str(args.cheat_hgt)+ ", cheat_dtemp = " +str(args.cheat_dtemp) + ", "
        "[[OTHERS]]: cpbl_strict = "+str(args.cpbl_strict)+", cc_advanced = "+str(args.cc_advanced)+", variable_mass = "+str(args.variable_mass)+ ", mode = "+str(args.mode) 
         +", fjumps = "+str(args.fjumps)+", cjumps = "+str(args.cjumps)))
    if step == 1 and args.tdiagnosis in ['SOD']:
        return(str("Diagnosis following Sodemann et al. (2008) with the following settings: " +
        "[[PRECIPITATION]] cprec_dqv = "+str(args.cprec_dqv)+ ", cprec_rh = " +str(args.cprec_rh) + ", " +
        "[[EVAPORATION]] cevap_dqv = 0.2, cevap_hgt < 1.5 * mean ABL, " +
        "[[SENSIBLE HEAT]] cheat_dTH = "+str(args.cheat_dtemp)+ ", cheat_hgt < 1.5 * mean ABL, " +
        "[[OTHERS]]: variable_mass = "+str(args.variable_mass)+ ", mode = "+str(args.mode)
         + ", fjumps = "+str(args.fjumps)+", cjumps = "+str(args.cjumps)
         + "; REFERENCE: " +
        "Sodemann, H., Schwierz, C., & Wernli, H. (2008). Interannual variability of Greenland winter precipitation sources: Lagrangian moisture diagnostic and North Atlantic Oscillation influence. Journal of Geophysical Research: Atmospheres, 113(D3). http://dx.doi.org/10.1029/2007JD008503"))
    if step == 1 and args.tdiagnosis in ['SOD2']:
        return(str("Diagnosis following Sodemann (2020) with the following settings: " +
        "[[PRECIPITATION]] cprec_dqv = "+str(args.cprec_dqv)+ ", cprec_rh = " +str(args.cprec_rh) + ", " +
        "[[EVAPORATION]] cevap_dqv = 0.1, " +
        "[[SENSIBLE HEAT]] cheat_dTH = "+str(args.cheat_dtemp)+ ", " +
        "[[OTHERS]]: variable_mass = "+str(args.variable_mass)+ ", mode = "+str(args.mode) 
         + ", fjumps = "+str(args.fjumps)+", cjumps = "+str(args.cjumps)
         + "; REFERENCE: " +
        "Sodemann, H. (2020). Beyond Turnover Time: Constraining the Lifetime Distribution of Water Vapor from Simple and Complex Approaches, Journal of the Atmospheric Sciences, 77, 413-433. https://doi.org/10.1175/JAS-D-18-0336.1"))
    if step == 1 and args.tdiagnosis in ['SAJ']:
        return(str("Diagnosis following Stohl and James (2004) with the following settings: " +
         "[[OTHERS]]: variable_mass = "+str(args.variable_mass)+ ", mode = "+str(args.mode)
         + ", fjumps = "+str(args.fjumps)+", cjumps = "+str(args.cjumps)
         + "; REFERENCE: " +
        "Stohl, A., & James, P. (2004). A Lagrangian analysis of the atmospheric branch of the global water cycle. Part I: Method description, validation, and demonstration for the August 2002 flooding in central Europe. Journal of Hydrometeorology, 5(4), 656-678. https://doi.org/10.1175/1525-7541(2004)005<0656:ALAOTA>2.0.CO;2"))
    
    ## 02_ATTRIBUTION
    if step == 2 or step == 3 and args.tdiagnosis in ['KAS']:
        return(str("Diagnosis following Schumacher & Keune (----) with the following settings: " +
        "[[PRECIPITATION]] cprec_dqv = "+str(args.cprec_dqv)+ ", cprec_rh = " +str(args.cprec_rh)+ ", cprec_dtemp = " +str(args.cprec_dtemp) + ", "
        "[[EVAPORATION]] cevap_cc = "+str(args.cevap_cc)+ ", cevap_hgt = " +str(args.cevap_hgt) + ", "
        "[[SENSIBLE HEAT]] cheat_cc = "+str(args.cheat_cc)+ ", cheat_hgt = " +str(args.cheat_hgt)+ ", cheat_dtemp = " +str(args.cheat_dtemp) + ", "
        "[[OTHERS]]: cpbl_strict = "+str(args.cpbl_strict)+", cc_advanced = "+str(args.cc_advanced)+", variable_mass = "+str(args.variable_mass)+ ", mode = "+str(args.mode) 
         +", fjumps = "+str(args.fjumps)+", cjumps = "+str(args.cjumps)+ ", "+
        "[[ATTRIBUTION]]: ctraj_len = "+str(args.ctraj_len)+", fallingdry = "+str(args.fallingdry)+", memento = "+str(args.memento) 
         ))
    if step == 2 or step == 3 and args.tdiagnosis in ['SOD']:
        return(str("Diagnosis following Sodemann et al. (2008) with the following settings: " +
        "[[PRECIPITATION]] cprec_dqv = "+str(args.cprec_dqv)+ ", cprec_rh = " +str(args.cprec_rh) + ", " +
        "[[EVAPORATION]] cevap_dqv = 0.2, cevap_hgt < 1.5 * mean ABL, " +
        "[[SENSIBLE HEAT]] cheat_dTH = "+str(args.cheat_dtemp)+ ", cheat_hgt < 1.5 * mean ABL, " +
        "[[OTHERS]]: variable_mass = "+str(args.variable_mass)+ ", mode = "+str(args.mode)
         + ", fjumps = "+str(args.fjumps)+", cjumps = "+str(args.cjumps) + ", "+
        "[[ATTRIBUTION]]: ctraj_len = "+str(args.ctraj_len)+", fallingdry = "+str(args.fallingdry)+", memento = "+str(args.memento) 
         + "; REFERENCE: " +
        "Sodemann, H., Schwierz, C., & Wernli, H. (2008). Interannual variability of Greenland winter precipitation sources: Lagrangian moisture diagnostic and North Atlantic Oscillation influence. Journal of Geophysical Research: Atmospheres, 113(D3). http://dx.doi.org/10.1029/2007JD008503"
        ))
    if step == 2 or step == 3 and args.tdiagnosis in ['SOD2']:
        return(str("Diagnosis following Sodemann (2020) with the following settings: " +
        "[[PRECIPITATION]] cprec_dqv = "+str(args.cprec_dqv)+ ", cprec_rh = " +str(args.cprec_rh) + ", " +
        "[[EVAPORATION]] cevap_dqv = 0.1, " +
        "[[SENSIBLE HEAT]] cheat_dTH = "+str(args.cheat_dtemp)+ ", " +
        "[[OTHERS]]: variable_mass = "+str(args.variable_mass)+ ", mode = "+str(args.mode) 
         + ", fjumps = "+str(args.fjumps)+", cjumps = "+str(args.cjumps) + ", "+ 
        "[[ATTRIBUTION]]: ctraj_len = "+str(args.ctraj_len)+", fallingdry = "+str(args.fallingdry)+", memento = "+str(args.memento) 
         + "; REFERENCE: " +
        "Sodemann, H. (2020). Beyond Turnover Time: Constraining the Lifetime Distribution of Water Vapor from Simple and Complex Approaches, Journal of the Atmospheric Sciences, 77, 413-433. https://doi.org/10.1175/JAS-D-18-0336.1"))
    if step == 2 or step == 3 and args.tdiagnosis in ['SAJ']:
        return(str("Diagnosis following Stohl and James (2004) with the following settings: " +
         "[[OTHERS]]: variable_mass = "+str(args.variable_mass)+ ", mode = "+str(args.mode)
         + ", fjumps = "+str(args.fjumps)+", cjumps = "+str(args.cjumps) + ", "+ 
        "[[ATTRIBUTION]]: ctraj_len = "+str(args.ctraj_len)+", memento = "+str(args.memento) 
         + "; REFERENCE: " +
        "Stohl, A., & James, P. (2004). A Lagrangian analysis of the atmospheric branch of the global water cycle. Part I: Method description, validation, and demonstration for the August 2002 flooding in central Europe. Journal of Hydrometeorology, 5(4), 656-678. https://doi.org/10.1175/1525-7541(2004)005<0656:ALAOTA>2.0.CO;2"))


def readpom(idate,      # run year
            ipath,      # input data path
            ifile_base, # loop over ifile_base filenames for each date
            verbose=True): # NOTE: temporary solution
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
                print("\t nparticle = ",nparticle, " |  ntrajstep = ",ntrajstep,"  | nvars = ",nvars)
            ary_dat     = pd.read_table(gzip.open(ifile, 'rb'), sep="\s+", header=None, skiprows=2)
            datav       = (np.asarray(ary_dat).flatten('C'))
            if dataar is None:
                dataar      = np.reshape(datav, (ntrajstep,nparticle,nvars), order='F')
            else:
                dataar      = np.append(dataar, np.reshape(datav, (ntrajstep,nparticle,nvars), order='F'), axis=1)
            # flip time axis    (TODO: flip axis depending on forward/backward flag)
    dataar          = dataar[::-1,:,:]

    return(dataar)


def checkpbl(cpbl,ztra,hpbl,maxhgt):
    if (cpbl == 1):
        #if (ztra < max(np.max(hpbl),maxhgt)).all():
        if (ztra < max(hpbl[1],hpbl[0],maxhgt)).all():
            return True
        else:
            return False
    if (cpbl == 2):
        if (ztra < max(hpbl[1],hpbl[0],maxhgt)).any():
            return True
        else:
            return False
    if (cpbl == 3):
        return True

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

def readheights(parray):

    """
    This is basically the same as 'glanceparcel',
    and as of now, is called only for the first 4 timesteps.
 
    NOTE: probably to be updated / removed

    """

    ## parcel information
    ztra    = parray[:,3]                   # height (m) 
    hpbl    = parray[:,7]                   # ABL height (m)

    return ztra, hpbl

def glanceparcel(parray):

    """
    This is used to take a first glance at parcel properties,
    and as of now, is called only for the first 4 timesteps.
 
    Required for arriving air criterion (which goes back 4 steps).
    However, this could be further optimized by only getting the
    last 4 values of hpbl, as for the other variables here,
    only the last / last 2 values are required at first.

    NOTE: beautify or remove!

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
    lat_mid,lon_mid = midpoint_on_sphere2(plat[0],plon[0],plat[1],plon[1]) # calculate midpoint position
    if (lon_mid>179.5):  lon_mid -= 360    # now shift all coords that otherwise would be allocated to +180 deg to - 180
    if (lon_mid<-180.5): lon_mid += 360    # same for the other direction; only correct for 1deg grid (as below)!
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
    elif (mlon<-180.5):
        mlon += 360      # same for the other direction; only correct for 1deg grid (as below)!
    # get grid index
    ind_lat = np.argmin(np.abs(glat-mlat))    # index on grid # ATTN, works only for 1deg grid
    ind_lon = np.argmin(np.abs(glon-mlon))    # index on grid # ATTN, works only for 1deg grid
    return ind_lat, ind_lon

def arrpindex(parray,glon,glat):
    # function to get arrival point index on grid
    lats    = parray[2]                   # latitude
    lons    = parray[1]                   # longitude
    if lons>180.0:
        lons = lons - 360                 # transform coordinates from [0 ... 360] to [-180 ... 180]
    ind_lat = np.argmin(np.abs(glat-lats))
    ind_lon = np.argmin(np.abs(glon-lons))
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
    pnparts             = nc_f.createVariable('P_n_part', 'f4', ('time','lat','lon'))
    enparts             = nc_f.createVariable('E_n_part', 'f4', ('time','lat','lon'))
    hnparts             = nc_f.createVariable('H_n_part', 'f4', ('time','lat','lon'))
    
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
    pnparts.units       = 'int'
    pnparts.long_name   = 'number of P parcels (mid pos.)'
    enparts.units       = 'int'
    enparts.long_name   = 'number of E parcels (mid pos.)'
    hnparts.units       = 'int'
    hnparts.long_name   = 'number of H parcels (mid pos.)'
    
    # write data
    times[:]            = nc4.date2num(fdate_seq, times.units, times.calendar)
    longitudes[:]       = glon
    latitudes[:]        = glat
    # close file
    nc_f.close()
    print("\n * Created empty file: "+ofile+" of dimension ("+str(len(fdate_seq))+","+str(glat.size)+","+str(glon.size)+") !")


def writenc(ofile,ix,ary_prec,ary_evap,ary_heat,ary_npart,ary_pnpart,ary_enpart,ary_hnpart):
    if verbose:
        print(" * Writing to netcdf...")

    nc_f = nc4.Dataset(ofile, 'r+')
    nc_f['P'][ix,:,:]       = ary_prec
    nc_f['E'][ix,:,:]       = ary_evap
    nc_f['H'][ix,:,:]       = ary_heat
    nc_f['n_part'][ix,:,:]  = ary_npart
    nc_f['P_n_part'][ix,:,:]  = ary_pnpart
    nc_f['E_n_part'][ix,:,:]  = ary_enpart
    nc_f['H_n_part'][ix,:,:]  = ary_hnpart
    
    # write data
    times[:]            = nc4.date2num(fdate_seq, times.units, times.calendar)
    longitudes[:]       = glon
    latitudes[:]        = glat
    # close file
    nc_f.close()
    print("\n * Created empty file: "+ofile+" of dimension ("+str(len(fdate_seq))+","+str(glat.size)+","+str(glon.size)+") !")


def writenc(ofile,ix,ary_prec,ary_evap,ary_heat,ary_npart,ary_pnpart,ary_enpart,ary_hnpart):
    if verbose:
        print(" * Writing to netcdf...")

    nc_f = nc4.Dataset(ofile, 'r+')
    nc_f['P'][ix,:,:]       = ary_prec
    nc_f['E'][ix,:,:]       = ary_evap
    nc_f['H'][ix,:,:]       = ary_heat
    nc_f['n_part'][ix,:,:]  = ary_npart
    nc_f['P_n_part'][ix,:,:]  = ary_pnpart
    nc_f['E_n_part'][ix,:,:]  = ary_enpart
    nc_f['H_n_part'][ix,:,:]  = ary_hnpart
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

def basicplot(array, lats, lons, title, coastlines=True, colorbar=True):
    """
    does exactly what the name suggests, 
    produces a very basic plot.
    
    array.shape == lats.size,lons.size
    
    data must be on a regular lat/lon grid!
    """
    
    ## obtain resolution and check if this is a 
    if (abs(lats[1]-lats[0]) == abs(lons[1]-lons[0])):
        res = abs(lats[1]-lats[0])
    else:
        raise SystemExit("------ ERROR: input coordinates not as required, ABORTING!")

    ## create quadrilateral corner coords (dimensions of X & Y one greater than array!)
    cornerlats = np.unique((lats-res/2).tolist()+(lats+res/2).tolist())
    cornerlons = np.unique((lons-res/2).tolist()+(lons+res/2).tolist())
    
    ## prepare figure & plot, add coastlines if desired
    plt.figure()
    ax = plt.axes(projection=ccrs.PlateCarree())
    plot_data = ax.pcolormesh(cornerlons, cornerlats, array, 
                              transform=ccrs.PlateCarree())
    if coastlines:
        coast = NaturalEarthFeature(category='physical', scale='50m',
                                    facecolor='none', name='coastline')
        ax.add_feature(coast, edgecolor='gray')
    if colorbar:
        plt.colorbar(plot_data, ax=ax, orientation = "horizontal")
    plt.title(title)
    plt.show()


def regrid_2Dto1deg(data, lats_in, lons_in, nmin_perpixel=1):
    """
    input: data with axes lats x lons, lons [-180 .. 180]
    action: regrid to regular 1 degree grid
    
    ===> CAUTION: lats, lons must increase monotonically!
    
    DSc, April 2019, alpha version
    """
    
    ## check if input coordinates do increase monotonically...
    if np.all(np.diff(lats_in)>0) and np.all(np.diff(lons_in)>0):
        pass
    else:
        print("---------- WARNING: function cannot handle input coordinates!")
        return(False,False,False)
    
    ## proceed
    lats_rounded     = np.round(lats_in)
    lats_new         = np.unique(lats_rounded)
    lons_rounded     = np.round(lons_in)

    if lons_rounded[-1] == 180 and lons_rounded[0] == -180:
        lons_rounded[lons_rounded==180] = -180  # set 180 to -180 degrees !
    lons_new         = np.unique(lons_rounded)
    
    data_rg = np.zeros(shape=(lats_new.size, lons_new.size))
    
    ## loop through
    for jlat in range(lats_new.size):
        xlat = np.where(lats_rounded==lats_new[jlat])[0]
        for jlon in range(lons_new.size):
            xlon = np.where(lons_rounded==lons_new[jlon])[0]
            
            #######################################
            if xlat.size*xlon.size>nmin_perpixel:
                pass
            else:
                continue
            #######################################
            
            if xlat.size==1 or xlon.size==1: # values at boundaries might be wrong.. (not double-checked)
                data_rg[jlat,jlon] = np.nanmean(data[xlat,xlon])
            else: # feeding in xlat, xlon directly results in weird shapes..
                data_rg[jlat,jlon] = np.nanmean(data[xlat[0]:xlat[-1]+1,xlon[0]:xlon[-1]+1])
                
    return(data_rg, lats_new, lons_new)
    
def freakshow(pommaskpath):
    """
    
    this function is appropriately called 'freakshow',
    because it is 
     - terribly ugly,
     - has only been checked for two ANYAS ecoregions (NGP, AUS)
    
    maskpath: should point to XXXXXXXX_mask.dat produced by particle-o-matic
    
    RETURNS:
        
        1 x 1° gridded mask as obtained from particle-o-matic,
        including coordinates
    """

    ## load mask file
    mask = np.asarray(pd.read_table(pommaskpath, sep="\s", engine='python', header=None))
    
    
    ## data are on a regular 0.2 x 0.2° grid, [90 .. -90], [0 .. 360] 
    dy = 180/mask.shape[0] 
    dx = 360/mask.shape[1]
    mask = np.flip(mask, axis=0) # flip latitudinal axis already
    assume_lats = np.arange(-90, 90+dy, dy) 
    assume_lons = np.arange(  0,   360, dx)
    ## regrid lons from 0 .. 359 to -180 .. 179
    mask_bu        = np.copy(mask)
    assume_lons_bu = np.copy(assume_lons)
    assume_lons[:int(assume_lons.size/2)] = assume_lons_bu[int(assume_lons.size/2):] - 360
    assume_lons[int(assume_lons.size/2):] = assume_lons_bu[:int(assume_lons.size/2)]
    mask[:,:int(assume_lons.size/2)] = mask_bu[:,int(assume_lons.size/2):]
    mask[:,int(assume_lons.size/2):] = mask_bu[:,:int(assume_lons.size/2)]
    
    ## use this crap function for regridding from 0.2 to 1.0°
    mask1deg, lat1deg, lon1deg = regrid_2Dto1deg(data=mask, lats_in=assume_lats, lons_in=assume_lons, nmin_perpixel=1)
    
    ## now ditch some pixels and pray it ends up making sense
    mask1deg[(mask1deg>0) & (mask1deg<0.5)] = 0
    mask1deg[(mask1deg>=0.5)]          = 1

    return(mask1deg, lat1deg, lon1deg)
    

def eraloader_12hourly(var, fullpath, maskpos, uptake_dates, lats, lons):
    """
    only for single years so far!!!
    """
    with nc4.Dataset(fullpath, mode='r') as f: # sshf_12hourly, tp_12hourly
        
        reflats = np.asarray(f['latitude'][:])
        reflons = np.asarray(f['longitude'][:])
        
        reftime = nc4.num2date(f['time'][:], f['time'].units, f['time'].calendar) 
        refdates = np.asarray([datetime.date(rt.year, rt.month, rt.day) for rt in reftime])

        ## first figure out what is needed!
        jbeg = np.min(np.where(refdates == uptake_dates[0])) 
        jend = np.max(np.where(refdates == uptake_dates[-1])) # 12-hourly input here!
        
        array = np.asarray(f[var][jbeg:jend+1,:,:])
        reftime  = reftime[jbeg:jend+1]
        refdates = refdates[jbeg:jend+1]
        
        ## mask positive values (E: condensation, H: downward fluxes, P: do NOT mask)
        if maskpos:
            array[array>0] = np.NaN
            
        ## aggregate to daily
        refudates = np.unique(refdates)
        daily = np.empty(shape=(refudates.size, reflats.size, reflons.size))
        for i in range(refudates.size):
            rud = refudates[i]
            daily[i,:,:] = np.nansum(array[np.where(refdates==rud)[0],:,:], axis=0) 
        if not np.array_equal(refudates, uptake_dates):
            raise SystemExit("---- no good")
    
        ## regrid (flip LATS; LON 0 .. 359 to -180 .. 179)
        reflats = np.flipud(reflats)
        daily  = np.flip(daily, axis=1)
        #-- above: flipping latitudes, below: from 0 .. 359 to -180 .. 179
        daily_bu  = np.copy(daily)
        reflons_bu = np.copy(reflons)
        reflons[:int(reflons.size/2)] = reflons_bu[int(reflons.size/2):] - 360
        reflons[int(reflons.size/2):] = reflons_bu[:int(reflons.size/2)]
        daily[:,:,:int(reflons.size/2)] = daily_bu[:,:,int(reflons.size/2):]
        daily[:,:,int(reflons.size/2):] = daily_bu[:,:,:int(reflons.size/2)]
    
        # check        
        if not (np.array_equal(lats, reflats) or not np.array_equal(lons, reflons)):
            raise SystemExit("---------- FATAL ERROR: regridded ERA-I coordinates don't match!")
            
        ## units... and sign
        if var=='e' and f[var].units == "m of water equivalent":
            daily *= -1e3 # flip sign
        elif var=='tp' and f[var].units == "m":
            daily *= 1e3
        elif var=='sshf' and f[var].units == "J m**-2":
            daily /= -86400 # flip sign
        else:
            raise SystemExit("---- aborted: no can do.")
        
        return daily
    
def alphascreener(alpha, var):
    alphasum = np.nansum(alpha, axis=0)
    plt.figure()
    plt.hist(alphasum.flatten(), bins=np.arange(0.1,2.0,0.01))
    plt.title(var+' :: alphas summed over arrival days -- auto range')
    plt.show()
    plt.figure()
    plt.hist(alphasum.flatten(), bins=np.arange(0.1,2.0,0.01))
    plt.title(var+' :: same as above, manual plotting range')
    plt.ylim(0,2e2)
    plt.show()
    print("--- this value should not really exceed 1.0 ===>", np.nanmax(alphasum))

def nanweight3Dary(array, weights):
    """
    purpose: weight 3-dim array and sum up along spatial dimensions
                if array does not contain data, corresponding weight is
                NOT taken into account as to not distort the average
    input:   3D array (time x lat x lon), weights of same shape
    output:  weighted averages, 1D (timeseries)
    """
    array   = array.reshape(array.shape[0],array.shape[1]*array.shape[2])
    weights = weights.reshape(array.shape)
    weightsum = np.zeros(shape=array.shape[0])
    for ii in range(array.shape[0]):    
        weightsum[ii] = np.nansum(weights[ii,np.where(~np.isnan(weights[ii,]*array[ii,]))])
    return( np.nansum(weights*array, axis=1)/weightsum )
    
def writefinalnc(ofile,fdate_seq,glon,glat,
                 Had, Had_Hs,
                 E2P, E2P_Es, E2P_Ps, E2P_EPs,strargs):
    
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
    heats               = nc_f.createVariable('Had', 'f4', ('time','lat','lon'))
    heats_Hs            = nc_f.createVariable('Had_Hs', 'f4', ('time','lat','lon'))
    evaps               = nc_f.createVariable('E2P', 'f4', ('time','lat','lon'))
    evaps_Es            = nc_f.createVariable('E2P_Es', 'f4', ('time','lat','lon'))
    evaps_Ps            = nc_f.createVariable('E2P_Ps', 'f4', ('time','lat','lon'))
    evaps_EPs           = nc_f.createVariable('E2P_EPs', 'f4', ('time','lat','lon'))
    
 
    # set attributes
    nc_f.description    = "03 - " + str(strargs)
    today               = datetime.datetime.now()
    nc_f.history        = "Created " + today.strftime("%d/%m/%Y %H:%M:%S") + " using hamster ((c) Dominik Schumacher and Jessica Keune)"
    times.units         = 'hours since 1900-01-01 00:00:00'
    times.calendar      = 'Standard' # do NOT use gregorian here!
    latitudes.units     = 'degrees_north'
    longitudes.units    = 'degrees_east'
    heats.units         = 'W m-2'
    heats.long_name     = 'advected surface sensible heat'
    heats_Hs.units      = 'W m-2'
    heats_Hs.long_name  = 'advected surface sensible heat, H-scaled' # this is garbage, I know
    evaps.units         = 'mm'
    evaps.long_name     = 'evaporation resulting in precipitation'
    evaps_Es.units      = 'mm'
    evaps_Es.long_name  = 'evaporation resulting in precipitation, E-scaled'
    evaps_Ps.units      = 'mm'
    evaps_Ps.long_name  = 'evaporation resulting in precipitation, P-scaled'
    evaps_EPs.units     = 'mm'
    evaps_EPs.long_name = 'evaporation resulting in precipitation, E&P-scaled'


    # write data
    times[:]            = nc4.date2num(fdate_seq, times.units, times.calendar)
    latitudes[:]        = glat
    longitudes[:]       = glon
    
    heats[:]      = Had[:]
    heats_Hs[:]   = Had_Hs[:]
    
    evaps[:]      = E2P[:]
    evaps_Es[:]   = E2P_Es[:]
    evaps_Ps[:]   = E2P_Ps[:]
    evaps_EPs[:]  = E2P_EPs[:]
      
    # close file
    nc_f.close()
    
    # print info
    print("\n * Created and wrote to file: "+ofile+" of dimension ("+str(len(fdate_seq))+","+str(glat.size)+","+str(glon.size)+") !")
