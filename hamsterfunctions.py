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
    parser.add_argument('--steps',      '-st',  help = "steps performed (1: diagnosis, 2: attribution, 3: bias correction)",                type = int,     default = 1)
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
    parser.add_argument('--precision',  '-f',   help = "precision for writing netcdf file variables (f4,f8)",           type = str,     default = "f8")
    parser.add_argument('--verbose',    '-v',   help = "verbose output (flag)",                                         type = str2bol, default = True,     nargs='?')
    parser.add_argument('--veryverbose','-vv',  help = "very verbose output (flag)",                                    type = str2bol, default = False,    nargs='?')
    parser.add_argument('--fallingdry', '-dry', help = "cut off trajectories falling dry (flag)",                       type = str2bol, default = True,     nargs='?')
    parser.add_argument('--memento',    '-mto', help = "keep track of trajectory history (flag)",                       type = str2bol, default = True,     nargs='?')
    parser.add_argument('--mattribution','-matt',help= "attribution method (for E2P as of now: random/linear)",         type = str,     default = "linear")
    parser.add_argument('--randomnit',  '-rnit',help = "minimum number of iterations for random attribution",           type = int,     default = 10)
    parser.add_argument('--explainp',   '-exp', help = "trajectory-based upscaling of E2P contributions",               type = str,     default = "none")
    parser.add_argument('--dupscale',   '-dups',help = "daily upscaling of E2P contributions",                          type = str2bol, default = False,    nargs='?')
    parser.add_argument('--mupscale',   '-mups',help = "monthly upscaling of E2P contributions",                        type = str2bol, default = False,    nargs='?')
    parser.add_argument('--variable_mass','-vm',help = "use variable mass (flag)",                                      type = str2bol, default = False,    nargs='?')
    parser.add_argument('--writestats','-ws',   help = "write additional stats to file (02 only; flag)",                type = str2bol, default = False,    nargs='?')
    parser.add_argument('--debug',      '-d',   help = "debugging option (flag)",                                       type = str2bol, default = False,    nargs='?')
    parser.add_argument('--gres',       '-r',   help = "output grid resolution (degrees)",                              type = float,   default = 1)
    parser.add_argument('--ryyyy',      '-ry',  help = "run name (here, YYYY, example: 2002, default: ayyyy)",          type = int,     default = None)
    parser.add_argument('--refdate',    '-rd',  help = "reference date (YYYYMMDDHH)",                                   type = str,     default = None)
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
    if (step == 2) and args.tdiagnosis in ['KAS']:
        return(str("Diagnosis following Schumacher & Keune (----) with the following settings: " +
        "[[PRECIPITATION]] cprec_dqv = "+str(args.cprec_dqv)+ ", cprec_rh = " +str(args.cprec_rh)+ ", cprec_dtemp = " +str(args.cprec_dtemp) + ", "
        "[[EVAPORATION]] cevap_cc = "+str(args.cevap_cc)+ ", cevap_hgt = " +str(args.cevap_hgt) + ", "
        "[[SENSIBLE HEAT]] cheat_cc = "+str(args.cheat_cc)+ ", cheat_hgt = " +str(args.cheat_hgt)+ ", cheat_dtemp = " +str(args.cheat_dtemp) + ", "
        "[[OTHERS]]: cpbl_strict = "+str(args.cpbl_strict)+", cc_advanced = "+str(args.cc_advanced)+", variable_mass = "+str(args.variable_mass)+ ", mode = "+str(args.mode) 
         +", fjumps = "+str(args.fjumps)+", cjumps = "+str(args.cjumps)+ ", "+
        "[[ATTRIBUTION]]: ctraj_len = "+str(args.ctraj_len)+", fallingdry = "+str(args.fallingdry)+", memento = "+str(args.memento) 
         ))
    if (step == 2) and args.tdiagnosis in ['SOD']:
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
    if (step == 2) and args.tdiagnosis in ['SOD2']:
        return(str("Diagnosis following Sodemann (2020) with the following settings: " +
        "[[PRECIPITATION]] cprec_dqv = "+str(args.cprec_dqv)+ ", cprec_rh = " +str(args.cprec_rh) + ", " +
        "[[EVAPORATION]] cevap_dqv = 0.1, " +
        "[[SENSIBLE HEAT]] cheat_dTH = "+str(args.cheat_dtemp)+ ", " +
        "[[OTHERS]]: variable_mass = "+str(args.variable_mass)+ ", mode = "+str(args.mode) 
         + ", fjumps = "+str(args.fjumps)+", cjumps = "+str(args.cjumps) + ", "+ 
        "[[ATTRIBUTION]]: ctraj_len = "+str(args.ctraj_len)+", fallingdry = "+str(args.fallingdry)+", memento = "+str(args.memento) 
         + "; REFERENCE: " +
        "Sodemann, H. (2020). Beyond Turnover Time: Constraining the Lifetime Distribution of Water Vapor from Simple and Complex Approaches, Journal of the Atmospheric Sciences, 77, 413-433. https://doi.org/10.1175/JAS-D-18-0336.1"))
    if (step == 2) and args.tdiagnosis in ['SAJ']:
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


def PBL_check(cpbl_strict, z, hpbl, sethpbl):
    """
    INPUT
        - PBL strictness flag; 1 (moderate), 2 (relaxed), 3 (fully relaxed)
        - parcel altitude
        - PBL heigt at parcel location
        - prescribed PBL height
    ACTION
        - raises any hpbl below sethpbl to this value; no effect if sethpbl=0
        - calculates whether parcel locations 'before' and 'after' each change,
          i.e. during analysis steps, are within PBL
        - depending on cpbl_strict, require both before & after to be
          inside PBL (1), either (2), or none (3), for change locations
    RETURN
        - returns boolean vector for all change locations
          (True if inside PBL), length given by z.size-1
    """
    hpbl[hpbl<sethpbl] = sethpbl
    befor_inside = np.logical_or( z[1:] < hpbl[:-1], z[1:]  < hpbl[1:])
    after_inside = np.logical_or(z[:-1] < hpbl[:-1], z[:-1] < hpbl[1:])
    if cpbl_strict == 1:
        change_inside = np.logical_and(befor_inside, after_inside)
    elif cpbl_strict == 2:
        change_inside = np.logical_or(befor_inside, after_inside)
    elif cpbl_strict == 3:
        change_inside = np.ones(dtype=bool, shape=befor_inside.size)
    return change_inside


def scale_mass(ary_val, ary_part, ary_rpart):
    ary_sval    = np.zeros(shape=(ary_val.shape)) 
    #for itime in range(ary_val.shape[0]):  
    #    with np.errstate(divide='ignore', invalid='ignore'):
    #        ary_sval[itime,:,:] = ary_val[itime,:,:] * np.nan_to_num( ary_rpart[:,:]/ary_part[itime,:,:] ) 
    with np.errstate(divide='ignore', invalid='ignore'):
         ary_sval[:,:] = ary_val[:,:] * np.nan_to_num( ary_rpart[:,:]/ary_part[:,:] ) 
    return(ary_sval)


def linear_discounter(v, min_gain):
    """
    ========= inputs =========
    v:        vector extending into past; v[0] is most recent, v[-1] least recent
    min_gain: minimum gain (dv/dt) to be considered
    ==========================
    """

    ## compute dv/dt, prepare dv_disc
    dv = v[:-1] - v[1:]
    # append initial v as a -fake- uptake 
    # (used to estimate fraction that cannot be attributed)
    dv = np.append(dv,v[-1])
    dv_disc = np.zeros(shape=len(dv))
    ## get indices of gains and losses
    idx_gains  = np.where(dv >= min_gain)[0]
    idx_losses = np.where(dv < 0)[0]

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

def linear_discounter2(qtot):
    dqdt = qtot[:-1] - qtot[1:]
    print(str(dqdt[0]))
    dqdt = dqdt[::-1]
    qtot = qtot[::-1]
    nt   = len(dqdt)
    print(str(dqdt[nt-1]))
    ## scale all humidity changes with subsequent uptakes and losses (rer= rain en route): 
    # using a discounting with the relative gain to the moisture content before the loss
    dqdt_scaled = dqdt
    rer         = abs(dqdt*((dqdt<0)+0))
    for ii in range(nt-2):
        if dqdt[ii]>0:
            for t in range((ii+1),(nt-1)):
                dqdt_scaled[ii] = dqdt_scaled[ii] - dqdt_scaled[ii]/qtot[t] * rer[t]
    dqdt_scaled[np.where(dqdt<=0)] =0 # reset rain en route
    return(dqdt_scaled[::-1])
    ## relative weight of scaled uptake to moisture content before final precip. event
    #fw  = dqdt_scaled / qtot[nt-1]
    #fw[which(fw<0)]=0
    #fw[which(is.na(fw))]=0
    #fw  = fw[::-1]
    #return(fw)

def linear_attribution_p(qv,iupt,explainp):
    dq_disc     = np.zeros(shape=qv.size)
    dq_disc[1:] = linear_discounter(v=qv[1:], min_gain=0)
    ## trajectory-based upscaling
    prec    = abs(qv[0]-qv[1])
    fw_orig = dq_disc/qv[1]
    if explainp=="full":
        # upscaling to 100% of trajectory
        cfac        = qv[1]/np.sum(dq_disc[iupt])
        etop        = prec*fw_orig*cfac
    elif explainp=="max":
        # upscaling to (100-IC)% of trajectory
        cfac        = (qv[1]-dq_disc[-1])/np.sum(dq_disc[iupt])
        etop        = prec*fw_orig*cfac
    elif explainp=="none":
        # no upscaling
        etop        = prec*fw_orig
    return(etop)

def calc_maxatt(qtot, iupt, verbose):
  dqdt = qtot[:-1] - qtot[1:]
  dqdt = np.append(dqdt,qtot[-1])
  nt   = len(dqdt)
  dqdt_max  = np.zeros(shape=nt)
  for ii in iupt[::-1]:
    try:
        imin        = np.argmin(qtot[1:ii])+1
    except:
        imin        = 1
    iatt    = qtot[imin]-round(np.sum(dqdt_max[imin:]),8)
    idqdt   = min(iatt,dqdt[ii]-dqdt_max[ii])
    dqdt_max[ii] = idqdt
  maxatt    = np.sum(dqdt_max)/abs(dqdt[0])
  if maxatt<1 and verbose:
    print(" * Maximum attribution along trajectory: {:.2f}".format(100*maxatt)+"%")
  return(maxatt)


def local_minima(x):
    return np.r_[True, x[1:] < x[:-1]] & np.r_[x[:-1] < x[1:], True]

def random_attribution_p(qtot,iupt,explainp,nmin=1,verbose=True,veryverbose=False):
  qtot = qtot*1000
  # This is only coded for precipitation as of now
  # with:
  # qtot = specific humidity
  # iupt = identified uptake locations
  # explainp = none(default)/full analogue to linear discounting & attribution
  #  - none: 100% is attributed to iupt + initial condition (~not explained)
  #  - full: 100% is attributed to iupt (if possible!)
  # nmin = tuning parameter; ~minimum iterations 
  #  - the higher this value, the more iterations, the uptake locations are covered
  #  - a value of 10 enforces min. 10 iterations
  dqdt  = qtot[:-1] - qtot[1:]
  # append initial condition as artificial uptake
  dqdt  = np.append(dqdt,qtot[-1])
  nt    = len(dqdt)
  if explainp=="none":
    iupt = np.append(iupt,nt-1)
  # indicator for potential uptake locations (1: yes, 0: no)
  pupt  = np.zeros(shape=nt)
  pupt[iupt]=1
  nupt  = len(np.where(pupt==1)[0])
  # adjust minimum number of iterations
  if nmin < nupt:
      nmin  = nupt
  # determine precip. (or max. attr. fraction)
  maxatt    = calc_maxatt(qtot, iupt, verbose)
  if maxatt>=1:
      prec  = dqdt[0]
  else:
      prec  = maxatt*dqdt[0]
  ## starting the random attribution loop
  dqdt_random = np.zeros(shape=nt)
  expl      = 0
  icount    = 0
  while round(expl,8) < round(abs(prec),8):
    i    = random.randint(0,nupt-1) # get a random uptake location number
    ii   = np.where(pupt==1)[0][i]  # uptake location index
    # determine maximum attribution for current uptake location and iteration
    try:
        imin    = np.argmin(qtot[1:ii])+1
    except:
        imin    = 1
    iatt        = qtot[imin]-round(np.sum(dqdt_random[imin:]),8)
    idqdt_max   = min(iatt,dqdt[ii]-dqdt_random[ii])
    # get random value
    rvalue  = random.uniform(0, min(idqdt_max, abs(prec)/nmin))
    if (expl+rvalue) > abs(prec):
      rvalue    = abs(prec)-expl
      if rvalue<0:
          print("OHOH")
    expl  += rvalue
    dqdt_random[ii] += rvalue
    icount += 1
    # safety exit (e.g. sum of dqdt_max cannot fully explain prec)
    if (icount >= 10000*nmin):
        print(" * Stopping at "+str(icount)+" iterations; attributed {:.2f}".format(100*np.sum(dqdt_random)/abs(prec))+"%.")
        break
  if veryverbose:
      print("  *** "+str(icount)+" Iterations for "+str(nupt)+" uptake locations with P={:.4f}".format(dqdt[0])+" g/kg with E2Prandom={:.4f}".format(np.sum(dqdt_random))+ " g/kg (attributed {:.2f}".format(100*np.sum(dqdt_random)/abs(dqdt[0]))+"%).")
  return(dqdt_random/1000)


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

def writeemptync(ofile,fdate_seq,glon,glat,strargs,precision):
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
    times               = nc_f.createVariable('time', 'f8', 'time')
    latitudes           = nc_f.createVariable('lat', 'f8', 'lat')
    longitudes          = nc_f.createVariable('lon', 'f8', 'lon')
    heats               = nc_f.createVariable('H', precision, ('time','lat','lon'))
    evaps               = nc_f.createVariable('E', precision, ('time','lat','lon'))
    precs               = nc_f.createVariable('P', precision, ('time','lat','lon'))
    nparts              = nc_f.createVariable('n_part', 'i4', ('time','lat','lon'))
    pnparts             = nc_f.createVariable('P_n_part', 'i4', ('time','lat','lon'))
    enparts             = nc_f.createVariable('E_n_part', 'i4', ('time','lat','lon'))
    hnparts             = nc_f.createVariable('H_n_part', 'i4', ('time','lat','lon'))
    
    # set attributes
    nc_f.title          = "Diagnosis (01) of FLEXPART fluxes"
    nc_f.description    = "01_diagnosis - " + str(strargs)
    today               = datetime.datetime.now()
    nc_f.history        = "Created " + today.strftime("%d/%m/%Y %H:%M:%S") + " using HAMSTER."
    nc_f.institution    = "Hydro-Climate Extremes Laboratory (H-CEL), Ghent University, Ghent, Belgium"
    nc_f.source         = "HAMSTER v0.1 ((c) Dominik Schumacher and Jessica Keune)" 
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

def writeemptync4D(ofile,fdate_seq,upt_days,glat,glon,strargs,precision):
                
    # delete nc file if it is present (avoiding error message)
    try:
        os.remove(ofile)
    except OSError:
        pass
        
    # create netCDF4 instance
    nc_f = nc4.Dataset(ofile,'w', format='NETCDF4')
    
    # create dimensions 
    nc_f.createDimension('time', len(fdate_seq))
    nc_f.createDimension('level', upt_days.size) # could use len() too
    nc_f.createDimension('lat', glat.size)
    nc_f.createDimension('lon', glon.size)
    
    # create variables
    atimes              = nc_f.createVariable('time', 'f8', 'time')
    utimes              = nc_f.createVariable('level', 'i4', 'level')
    latitudes           = nc_f.createVariable('lat', 'f8', 'lat')
    longitudes          = nc_f.createVariable('lon', 'f8', 'lon')
    heats               = nc_f.createVariable('H', precision, ('time','level','lat','lon'))
    etops               = nc_f.createVariable('E2P', precision, ('time','level','lat','lon'))
    
    # set attributes
    nc_f.title          = "Attribution (02) of sources using FLEXPART output"
    nc_f.description    = "02_attribution - " + str(strargs)
    today               = datetime.datetime.now()
    nc_f.history        = "Created " + today.strftime("%d/%m/%Y %H:%M:%S") + " using HAMSTER."
    nc_f.institution    = "Hydro-Climate Extremes Laboratory (H-CEL), Ghent University, Ghent, Belgium"
    nc_f.source         = "HAMSTER v0.1 ((c) Dominik Schumacher and Jessica Keune)" 
    atimes.units        = 'hours since 1900-01-01 00:00:00'
    atimes.calendar     = 'Standard' 
    utimes.long_name    = 'Difference between uptake and arrival time, in days'
    utimes.units        = 'day'
    latitudes.units     = 'degrees_north'
    longitudes.units    = 'degrees_east'
    heats.units         = 'W m-2'
    heats.long_name	= 'surface sensible heat flux'
    etops.units         = 'mm'
    etops.long_name	= 'evaporation resulting in precipitation'
  
    # write data
    atimes[:]           = nc4.date2num(fdate_seq, atimes.units, atimes.calendar)[:]
    utimes[:]           = upt_days[:]
    longitudes[:]       = glon[:]
    latitudes[:]        = glat[:]
        
    # close file
    nc_f.close()

    print("\n * Created empty file: "+ofile+" of dimension ("+str(len(fdate_seq))+","+str(upt_days.size)+","+str(glat.size)+","+str(glon.size)+") !")
        
def writenc4D(ofile,ix,ary_etop,ary_heat):
    if verbose:
        print(" * Writing to netcdf...")

    nc_f = nc4.Dataset(ofile, 'r+')
    nc_f['E2P'][ix,:,:,:]     = ary_etop[:,:,:]
    nc_f['H'][ix,:,:,:]       = ary_heat[:,:,:]
    nc_f.close()

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
    

def eraloader_12hourly(var, datapath, maskpos, maskneg, uptake_years, uptake_dates, lats, lons):
    """
    quickly adjusted to enable multi-annual support at the cost of reading in
    two entire years, instead of just what is needed.
    """

    uyears = np.unique(uptake_years)

    with nc4.Dataset(datapath+str(uyears[0])+'.nc', mode='r') as f: 
        reflats   = np.asarray(f['latitude'][:])
        reflons   = np.asarray(f['longitude'][:])
        reftime   = nc4.num2date(f['time'][:], f['time'].units, f['time'].calendar) - timedelta(hours=12)
        refdates  = np.asarray([datetime.date(rt.year, rt.month, rt.day) for rt in reftime])
        array     = np.asarray(f[var][:,:,:]) # not ideal to load everything..
        units     = f[var].units

    for ii in range(1,uyears.size):

        with nc4.Dataset(datapath+str(uyears[ii])+'.nc', mode='r') as f:
            reftimeY =  nc4.num2date(f['time'][:], f['time'].units, f['time'].calendar) - timedelta(hours=12)
            reftime  = np.concatenate((reftime, reftimeY))
            refdates = np.concatenate((refdates, np.asarray([datetime.date(rt.year, rt.month, rt.day) for rt in reftimeY])))
            array    = np.concatenate((array, np.asarray(f[var][:,:,:])), axis=0)

    ## first figure out what is needed! NOTE: this is less elegant now due to multi-annual support.
    jbeg = np.min(np.where(refdates == uptake_dates[0])) 
    jend = np.max(np.where(refdates == uptake_dates[-1])) # 12-hourly input here!
        
    array    = array[jbeg:jend+1,:,:]
    reftime  = reftime[jbeg:jend+1]
    refdates = refdates[jbeg:jend+1]
        
    ## mask positive values (E: condensation, H: downward fluxes, P: do NOT mask)
    if maskpos:
        array[array>0] = 0
    ## mask negative values (P only!)
    if maskneg:
        array[array<0] = 0
            
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
    if var=='e' and units == "m of water equivalent":
        daily *= -1e3 # flip sign
    elif var=='tp' and units == "m":
        daily *= 1e3
    elif var=='sshf' and units == "J m**-2":
        daily /= -86400 # flip sign
    else:
        raise SystemExit("---- aborted: no can do.")
    
    return daily

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
                 E2P, E2P_Es, E2P_Ps, E2P_EPs,strargs,precision):
    
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
    times               = nc_f.createVariable('time', 'f8', 'time')
    latitudes           = nc_f.createVariable('lat', 'f8', 'lat')
    longitudes          = nc_f.createVariable('lon', 'f8', 'lon')
    heats               = nc_f.createVariable('Had', precision, ('time','lat','lon'))
    heats_Hs            = nc_f.createVariable('Had_Hs', precision, ('time','lat','lon'))
    evaps               = nc_f.createVariable('E2P', precision, ('time','lat','lon'))
    evaps_Es            = nc_f.createVariable('E2P_Es', precision, ('time','lat','lon'))
    evaps_Ps            = nc_f.createVariable('E2P_Ps', precision, ('time','lat','lon'))
    evaps_EPs           = nc_f.createVariable('E2P_EPs', precision, ('time','lat','lon'))
    
 
    # set attributes
    nc_f.title          = "Bias-corrected source-sink relationships from FLEXPART"
    nc_f.description    = str(strargs)
    today               = datetime.datetime.now()
    nc_f.history        = "Created " + today.strftime("%d/%m/%Y %H:%M:%S") + " using HAMSTER."
    nc_f.institution    = "Hydro-Climate Extremes Laboratory (H-CEL), Ghent University, Ghent, Belgium"
    nc_f.source         = "HAMSTER v0.1 ((c) Dominik Schumacher and Jessica Keune)" 
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

def append2csv(filename, listvals):
    # Open file in append mode
    with open(filename, 'a+', newline='\n') as write_obj:
        csv_writer = csv.writer(write_obj, delimiter='\t', lineterminator='\n')
        csv_writer.writerow(listvals)

def preloop(datetime_bgn, uptdatetime_bgn, timestep,
            ipath, ifile_base, ryyyy,
            mask, mlat, mlon, maskval,
            pidlog, tml,
            verbose):

    ## p1) create required datetime string object
    predatetime_bgn = uptdatetime_bgn + datetime.timedelta(hours=3)
    predatetime_end = datetime_bgn
    predatetime_seq = []
    idatetime       = predatetime_bgn
    while idatetime < predatetime_end:
        predatetime_seq.append(idatetime.strftime('%Y%m%d%H'))
        idatetime += timestep # timestep was defined above
    npretime = len(predatetime_seq)

    if verbose:
        print("\n--------------------------------------------------------------------------------------")
        print("\n ! performing pre-loop to log advected parcels arriving prior to analysis time")
        print("\n ! estimating remaining time for pre-loop ...")

    ## p2) loop through files (.. to log in-ABL hits)
    pretic = timeit.default_timer()
    for pix in range(npretime):
        if verbose and pix==1:
            pretoc = timeit.default_timer()
            print("  ---> "+str(round(npretime*(pretoc-pretic)/60, 2))+" minutes to go, grab a coffee..")

        ## p3) read in all files associated with data --> ary is of dimension (ntrajlen x nparcels x nvars)
        ary = readpom( idate    = predatetime_seq[pix],
                       ipath    = ipath+"/"+str(ryyyy),
                       ifile_base = ifile_base,
                       verbose=False) # NOTE: ugly, but this way, other instances need no change (per default: True)

        nparcel   = ary.shape[1]
        ntot    = range(nparcel)

        ## p4) now loop through parcels
        for i in ntot:
            ## check for arriving parcels
            alat_ind, alon_ind = arrpindex(ary[0,i,:],glon=mlon,glat=mlat)
            if not mask[alat_ind,alon_ind]==maskval:
               continue
            ## read ONLY parcel and ABL heights
            hgt, hpbl = readheights(ary[:4,i,:])

            ## p5) LOG ONLY parcels arriving in PBL (or nocturnal layer)
            if ( hgt[0] < np.max(hpbl[:4]) ):
                ID = int(ary[0,i,0])
                pidlog[ID] = pix - tml # NOTE: tml != npretime (double-check?)
    return( pidlog )

def uptake_locator_KAS(c_hgt, cpbl_strict, hgt, hpbl,
                       dX, dtemp, dA, dB, c_cc, dAdB):
    ## NOTE: for heat, dX==dB, but not for evap!
    is_inpbl    = PBL_check(cpbl_strict, z=hgt, hpbl=hpbl, sethpbl=c_hgt)
    is_uptk     = dX > dtemp
    is_uptkcc   = np.abs(dA) < c_cc * dB * dAdB
    return( np.where(np.logical_and(is_inpbl, np.logical_and(is_uptk, is_uptkcc)))[0] )
    

def convert2daily(xar,ftime,dtime,fagg="mean"):
    dates   = np.asarray([datetime.date(it.year, it.month, it.day) for it in ftime])

    if ftime[0].hour in [0, 6, 12, 18]:
        ## simple fix, subtract 3 hours
        ftime   = np.asarray([t - datetime.timedelta(hours=3) for t in ftime])
        dates   = np.asarray([datetime.date(it.year, it.month, it.day) for it in ftime])
    elif ftime[0].hour in [3, 9, 15, 21]:
        ## NOTE: this is the new norm! retain "old style" for now, though    
        pass
    dtime   = np.unique(dates)

    xtot = np.zeros(shape=(ftime.size, xar.shape[1], xar.shape[2]))

    ## this isn't fast or elegant, but works for literally anything sub-daily
    for i in range(dtime.size):

        iud = dtime[i]
        sel = np.where(dates == iud)[0]

        ## TODO: clean up; there should be a check whether 4 files are present, imo
        if sel.size != 4:
            warnings.warn("\n\n----------------- WARNING: this should NEVER OCCUR; daily aggregation IMPROPER (files missing!)\n\n")

        if fagg=="sum":
            xtot[i,:,:] = np.nansum(xar[sel, :, :], axis=0)
        if fagg=="mean":
            xtot[i,:,:] = np.nanmean(xar[sel, :, :], axis=0)

    return(xtot)

def read_diagdata(opathD,ofile_base,ryyyy,uptake_time,var="E"):
    # get required months
    ayears  = np.asarray([ut.year  for ut in uptake_time])
    amonths = np.asarray([ut.month for ut in uptake_time])
    uset    = np.unique(np.column_stack((ayears, amonths)), axis=0)
    uyears, umonths = uset[:,0], uset[:,1]

    # loop over months
    for jj in range(umonths.size):
        uyr  = str(uyears[jj])
        umon = str(umonths[jj])
        diagfile = str(opathD)+"/"+str(ofile_base)+"_diag_r"+str(ryyyy)[-2:]+"_"+str(uyr)+"-"+umon.zfill(2)+".nc"
        with nc4.Dataset(diagfile, mode="r") as f:
            if var not in ["grid"]:
                ix     = f[var][:]
                timex  = nc4.num2date(f['time'][:], f['time'].units, f['time'].calendar)

        ## concatenate 'em!
        if var=="grid":
            with nc4.Dataset(diagfile, mode="r") as f:
                lats    = f['lat'][:]
                lons    = f['lon'][:]
        else:
            # concatenate 
            if umon == str(umonths[0]):
                x       = np.copy(ix)
                ftime   = np.copy(timex)
            else:
                x       = np.concatenate((x, ix), axis=0)
                ftime   = np.concatenate((ftime, timex))

    # return
    if var=="time":
        return(ftime)
    if var=="grid":
        return(lats, lons)
    if var not in ["grid","time"]:
        return(x)

def gridcheck(lats,totlats,lons,totlons):
    if not np.array_equal(lats, totlats) or not np.array_equal(lons, totlons):
        raise SystemExit("--- ERROR: your grids aren't identical...")

def datecheck(idate,dateseq):
    if idate not in dateseq: 
        raise SystemExit("\n !!! ERROR: INPUT DATA MISSING: date "+str(idate)+" not available as output from 01_diagnosis! Aborting here. !!!\n")

def calc_alpha(top,bot):
    alpha   = np.divide(top,bot)
    ## NOTE: as of now, there is absolutely no check whatsoever concerning
    ## the fractions; if e.g. only 3 6-hourly values are used to generate
    ## daily diagnosis data, this can result in a division by zero above,
    ## so that scaled data blows up to infinity (this actually happened).
    ## hence, check if any alpha clearly exceeds 1, and warn the user
    ## AGAIN that the output cannot be fully trusted (but continue)
    if np.any(alpha>1.0001) or np.any(np.isinf(alpha)):
        print(" \n  \t !!! WARNING: scaling fractions exceed 1 !!!")
        print(" \t Maximum scaling fraction: " + str(np.max(np.nan_to_num(alpha)))+"\n")
    return(alpha)

def udays2udate(atime,utime_srt):
    utime_first = atime[0] - timedelta(days=utime_srt.size-1) # utime_srt.size-1 == trajlen (in days)
    uptake_time = np.asarray([utime_first+timedelta(days=nday) for nday in range(utime_srt.size-1+atime.size)])
    return(uptake_time)

def expand4Darray(myarray,atime,utime_srt,veryverbose):
    utime       = udays2udate(atime,utime_srt)
    myshape     = myarray.shape
    myarray_exp = np.empty(shape=(myshape[0],utime.size,myshape[2],myshape[3]))
    if veryverbose:
        print(" * Expanding array from "+str(myshape)+" to "+str(myarray_exp.shape))
    for iat in range(atime.size):
        myarray_exp[iat,iat:iat+utime_srt.size,:,:] = myarray[iat,:,:,:]
    return(myarray_exp)

def date2year(mydates):
    return( np.asarray([it.year for it in mydates]))
def date2month(mydates):
    return( np.asarray([it.month for it in mydates]))
def cal2date(mydates):
    return( np.asarray([datetime.date(it.year, it.month, it.day) for it in mydates]))

def convert_mm_m3(myarray,areas):
    # we ALWAYS follow the array dimensions order: (anything(s) x lat x lon) here
    if len(areas.shape) > 1:
        # attention: broadcasting with a 2D array in numpy can lead to differences 
        # 1e-9 (f4) or 1e-17 (f8)
        carray = np.multiply(areas,myarray/1e3)
    if len(areas.shape) == 1:
        ## swap axes to enable numpy broadcasting;
        ## (a,b,c,d x b,c,d OK; a,b,c,d x a NOT OK)
        ldim   = len(myarray.shape)-1  
        carray = np.swapaxes(areas*np.moveaxis(myarray/1e3, ldim-1, ldim), ldim-1, ldim)
    return(carray)
    #return( np.multiply(areas,myarray/1e3) ) 
def convert_m3_mm(myarray,areas):
    # we ALWAYS follow the array dimensions order: (anything(s) x lat x lon) here
    if len(areas.shape) > 1:
        # attention: broadcasting with a 2D array in numpy can lead to differences 
        # 1e-9 (f4) or 1e-17 (f8)
        carray = np.nan_to_num(np.divide(myarray*1e3,areas))
    if len(areas.shape) == 1:
        ## swap axes to enable numpy broadcasting;
        ## (a,b,c,d x b,c,d OK; a,b,c,d x a NOT OK)
        ldim   = len(myarray.shape)-1 
        carray = np.swapaxes(np.nan_to_num(np.divide(np.moveaxis(myarray*1e3, ldim-1, ldim), areas)), ldim-1, ldim)
    return(carray) 
