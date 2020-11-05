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
    parser.add_argument('--steps',      '-st',  help = "steps performed (0: flex2traj, 1: diagnosis, 2: attribution, 3: bias correction)", metavar ="", type = int,     default = 1)
    parser.add_argument('--ayyyy',      '-ay',  help = "analysis year (YYYY)",                                          metavar ="", type = int,     default = 2002)
    parser.add_argument('--am',         '-am',  help = "analysis month (M)",                                            metavar ="", type = int,     default = 1)
    parser.add_argument('--ad',         '-ad',  help = "analysis day (D)",                                              metavar ="", type = int,     default = 1)
    parser.add_argument('--mode',       '-m',   help = "mode (test,oper)",                                              metavar ="", type = str,     default = "oper")
    parser.add_argument('--expid',      '-id',  help = "experiment ID (string, example versionA)",                      metavar ="", type = str,     default = "FXv")
    parser.add_argument('--tdiagnosis', '-dgn', help = "diagnosis method (KAS, SOD, SOD2)",                             metavar ="", type = str,     default = "KAS")
    parser.add_argument('--maskval',    '-mv',  help = "use <value> from maskfile for masking",                         metavar ="", type = int,     default = 1)
    parser.add_argument('--ctraj_len',  '-len', help = "threshold for maximum allowed trajectory length in days",       metavar ="", type = int,     default = 10)
    parser.add_argument('--cprec_dqv',  '-cpq', help = "threshold for detection of P based on delta(qv)",               metavar ="", type = float,   default = 0)
    parser.add_argument('--cprec_rh',   '-cpr', help = "threshold for detection of P based on RH",                      metavar ="", type = float,   default = 80)
    parser.add_argument('--cprec_dtemp','-cpt', help = "threshold for detection of P based on delta(T)",                metavar ="", type = float,   default = 0)
    parser.add_argument('--cevap_cc',   '-cec', help = "threshold for detection of E based on CC criterion",            metavar ="", type = float,   default = 0.7)
    parser.add_argument('--cevap_hgt',  '-ceh', help = "threshold for detection of E using a maximum height",           metavar ="", type = float,   default = 0)
    parser.add_argument('--cheat_cc',   '-chc', help = "threshold for detection of H based on CC criterion",            metavar ="", type = float,   default = 0.7)
    parser.add_argument('--cheat_hgt',  '-chh', help = "threshold for detection of H using a maximum height",           metavar ="", type = float,   default = 0)
    parser.add_argument('--cheat_dtemp','-cht', help = "threshold for detection of H using a minimum delta(T)",         metavar ="", type = float,   default = 0)
    parser.add_argument('--cpbl_strict','-pbl', help = "filter for PBL - 1: both within max, 2: one within max, 3: not used", metavar ="", type = int,     default = 1)
    parser.add_argument('--cc_advanced','-cc',  help = "use advanced CC criterion (flag, DEVELOPMENT)",                 metavar ="", type = str2bol, default = False,    nargs='?')
    parser.add_argument('--timethis',   '-t',   help = "time the main loop (flag)",                                     metavar ="", type = str2bol, default = False,    nargs='?')
    parser.add_argument('--write_netcdf','-o',  help = "write netcdf output (flag)",                                    metavar ="", type = str2bol, default = True,     nargs='?')
    parser.add_argument('--write_month','-mo',  help = "write monthly aggreagted netcdf output (03 only; flag)",        metavar ="", type = str2bol, default = False,     nargs='?')
    parser.add_argument('--precision',  '-f',   help = "precision for writing netcdf file variables (f4,f8)",           metavar ="", type = str,     default = "f8")
    parser.add_argument('--verbose',    '-v',   help = "verbose output (flag)",                                         metavar ="", type = str2bol, default = True,     nargs='?')
    parser.add_argument('--veryverbose','-vv',  help = "very verbose output (flag)",                                    metavar ="", type = str2bol, default = False,    nargs='?')
    parser.add_argument('--fallingdry', '-dry', help = "cut off trajectories falling dry (flag)",                       metavar ="", type = str2bol, default = True,     nargs='?')
    parser.add_argument('--memento',    '-mto', help = "keep track of trajectory history (02 only - needed for Had; flag)", metavar ="", type = str2bol, default = True,     nargs='?')
    parser.add_argument('--mattribution','-matt',help= "attribution method (for E2P as of now: random/linear)",         metavar ="", type = str,     default = "linear")
    parser.add_argument('--ratt_nit',   '-rnit',help = "minimum number of iterations for random attribution",           metavar ="", type = int,     default = 10)
    parser.add_argument('--ratt_forcall','-rall',help = "enforcing the attribution to all uptake locations (random att.)", metavar ="", type = str2bol, default = False, nargs='?')
    parser.add_argument('--explainp',   '-exp', help = "trajectory-based upscaling of E2P contributions (02 only: none/max/full)",  metavar ="", type = str,     default = "none")
    parser.add_argument('--dupscale',   '-dups',help = "daily upscaling of E2P contributions (02 only; flag)",          metavar ="", type = str2bol, default = False,    nargs='?')
    parser.add_argument('--mupscale',   '-mups',help = "monthly upscaling of E2P contributions (02 only; flag)",        metavar ="", type = str2bol, default = False,    nargs='?')
    parser.add_argument('--bc_useattp', '-uatt',help = "use precipitation from attribution for bias-correction",        metavar ="", type = str2bol, default = False,    nargs='?')
    parser.add_argument('--bc_time',    '-bct', help = "time scale for bias-correction (daily/monthly)",                metavar ="", type = str,     default = "daily")
    parser.add_argument('--variable_mass','-vm',help = "use variable mass (flag)",                                      metavar ="", type = str2bol, default = False,    nargs='?')
    parser.add_argument('--writestats', '-ws',  help = "write additional stats to file (02 and 03 only; flag)",         metavar ="", type = str2bol, default = False,    nargs='?')
    parser.add_argument('--bc_aggbwtime','-aggbt',help = "aggregate backward time (03 only; flag)",                     metavar ="", type = str2bol, default = True,    nargs='?')
    parser.add_argument('--debug',      '-d',   help = "debugging option (flag)",                                       metavar ="", type = str2bol, default = False,    nargs='?')
    parser.add_argument('--gres',       '-r',   help = "output grid resolution (degrees)",                              metavar ="", type = float,   default = 1)
    parser.add_argument('--ryyyy',      '-ry',  help = "run name (here, YYYY, example: 2002, default: ayyyy)",          metavar ="", type = int,     default = None)
    parser.add_argument('--refdate',    '-rd',  help = "reference date (YYYYMMDDHH)",                                   metavar ="", type = str,     default = None)
    #print(parser.format_help())
    args = parser.parse_args()  # namespace
    # handle None cases already
    if args.ryyyy is None:
        args.ryyyy   = args.ayyyy
    if args.refdate is None:
        args.refdate = str(args.ryyyy)+"123118"
    return args

def printsettings(args,step):
    ## 01_DIAGNOSIS
    if step == 1 and args.tdiagnosis in ['KAS']:
        return(str("Diagnosis following Schumacher & Keune (----) with the following settings: " +
        "[[PRECIPITATION]] cprec_dqv = "+str(args.cprec_dqv)+ ", cprec_rh = " +str(args.cprec_rh)+
        ", cprec_dtemp = " +str(args.cprec_dtemp) + ", "
        "[[EVAPORATION]] cevap_cc = "+str(args.cevap_cc)+ ", cevap_hgt = " +str(args.cevap_hgt) + ", "
        "[[SENSIBLE HEAT]] cheat_cc = "+str(args.cheat_cc)+ ", cheat_hgt = " +str(args.cheat_hgt)+
        ", cheat_dtemp = " +str(args.cheat_dtemp) + ", "
        "[[OTHERS]]: cpbl_strict = "+str(args.cpbl_strict)+", cc_advanced = "+str(args.cc_advanced)+
        ", variable_mass = "+str(args.variable_mass)+ ", mode = "+str(args.mode))) 
    if step == 1 and args.tdiagnosis in ['SOD']:
        return(str("Diagnosis following Sodemann et al. (2008) with the following settings: " +
        "[[PRECIPITATION]] cprec_dqv = "+str(args.cprec_dqv)+ ", cprec_rh = " +str(args.cprec_rh) + ", " +
        "[[EVAPORATION]] cevap_dqv = 0.2, cevap_hgt < 1.5 * mean ABL, " +
        "[[SENSIBLE HEAT]] cheat_dTH = "+str(args.cheat_dtemp)+ ", cheat_hgt < 1.5 * mean ABL, " +
        "[[OTHERS]]: variable_mass = "+str(args.variable_mass)+ ", mode = "+str(args.mode)
         + "; REFERENCE: " +
        "Sodemann, H., Schwierz, C., & Wernli, H. (2008). Interannual variability of Greenland winter precipitation sources: Lagrangian moisture diagnostic and North Atlantic Oscillation influence. Journal of Geophysical Research: Atmospheres, 113(D3). http://dx.doi.org/10.1029/2007JD008503"))
    if step == 1 and args.tdiagnosis in ['SOD2']:
        return(str("Diagnosis following Sodemann (2020) with the following settings: " +
        "[[PRECIPITATION]] cprec_dqv = "+str(args.cprec_dqv)+ ", cprec_rh = " +str(args.cprec_rh) + ", " +
        "[[EVAPORATION]] cevap_dqv = 0.1, " +
        "[[SENSIBLE HEAT]] cheat_dTH = "+str(args.cheat_dtemp)+ ", " +
        "[[OTHERS]]: variable_mass = "+str(args.variable_mass)+ ", mode = "+str(args.mode) 
         + "; REFERENCE: " +
        "Sodemann, H. (2020). Beyond Turnover Time: Constraining the Lifetime Distribution of Water Vapor from Simple and Complex Approaches, Journal of the Atmospheric Sciences, 77, 413-433. https://doi.org/10.1175/JAS-D-18-0336.1"))
    
    ## 02_ATTRIBUTION
    if (step == 2) and args.tdiagnosis in ['KAS']:
        return(str("Diagnosis following Schumacher & Keune (----) with the following settings: " +
        "[[PRECIPITATION]] cprec_dqv = "+str(args.cprec_dqv)+ ", cprec_rh = " +str(args.cprec_rh)+ 
        ", cprec_dtemp = " +str(args.cprec_dtemp) + ", "+
        "[[EVAPORATION]] cevap_cc = "+str(args.cevap_cc)+ ", cevap_hgt = " +str(args.cevap_hgt) + ", "
        "[[SENSIBLE HEAT]] cheat_cc = "+str(args.cheat_cc)+ ", cheat_hgt = " +str(args.cheat_hgt)+ 
        ", cheat_dtemp = " +str(args.cheat_dtemp) + ", "+
        "[[OTHERS]]: cpbl_strict = "+str(args.cpbl_strict)+", cc_advanced = "+str(args.cc_advanced)+
        ", variable_mass = "+str(args.variable_mass)+ ", mode = "+str(args.mode)+", "+
        "[[ATTRIBUTION]]: ctraj_len = "+str(args.ctraj_len)+", fallingdry = "+str(args.fallingdry)+
        ", memento = "+str(args.memento)))
    if (step == 2) and args.tdiagnosis in ['SOD']:
        return(str("Diagnosis following Sodemann et al. (2008) with the following settings: " +
        "[[PRECIPITATION]] cprec_dqv = "+str(args.cprec_dqv)+ ", cprec_rh = " +str(args.cprec_rh) + ", " +
        "[[EVAPORATION]] cevap_dqv = 0.2, cevap_hgt < 1.5 * mean ABL, " +
        "[[SENSIBLE HEAT]] cheat_dTH = "+str(args.cheat_dtemp)+ ", cheat_hgt < 1.5 * mean ABL, " +
        "[[OTHERS]]: variable_mass = "+str(args.variable_mass)+ ", mode = "+str(args.mode)+", "+
        "[[ATTRIBUTION]]: ctraj_len = "+str(args.ctraj_len)+", fallingdry = "+str(args.fallingdry)+
        ", memento = "+str(args.memento)
         + "; REFERENCE: " +
        "Sodemann, H., Schwierz, C., & Wernli, H. (2008). Interannual variability of Greenland winter precipitation sources: Lagrangian moisture diagnostic and North Atlantic Oscillation influence. Journal of Geophysical Research: Atmospheres, 113(D3). http://dx.doi.org/10.1029/2007JD008503"
        ))
    if (step == 2) and args.tdiagnosis in ['SOD2']:
        return(str("Diagnosis following Sodemann (2020) with the following settings: " +
        "[[PRECIPITATION]] cprec_dqv = "+str(args.cprec_dqv)+ ", cprec_rh = " +str(args.cprec_rh) + ", " +
        "[[EVAPORATION]] cevap_dqv = 0.1, " +
        "[[SENSIBLE HEAT]] cheat_dTH = "+str(args.cheat_dtemp)+ ", " +
        "[[OTHERS]]: variable_mass = "+str(args.variable_mass)+ ", mode = "+str(args.mode)+ ", "+ 
        "[[ATTRIBUTION]]: ctraj_len = "+str(args.ctraj_len)+", fallingdry = "+str(args.fallingdry)+
        ", memento = "+str(args.memento)
         + "; REFERENCE: " +
        "Sodemann, H. (2020). Beyond Turnover Time: Constraining the Lifetime Distribution of Water Vapor from Simple and Complex Approaches, Journal of the Atmospheric Sciences, 77, 413-433. https://doi.org/10.1175/JAS-D-18-0336.1"))


def readtraj(idate, # date as string [YYYYMMDDHH]
            ipath,
            ifile_base,
            verbose=True):
    # reads in *h5 data from flex2traj and flips the time axis
    # returns data array of dimension (ntrajlength x nparticles x nvars)
    # Check if file exists /file format
    ifile   = str(ipath+"/"+ifile_base+"_"+idate+".h5")
    if not os.path.isfile(ifile):
        raise SystemExit(ifile + " does not exist!")
    elif os.path.isfile(ifile):
        if verbose:
            print(" Reading " + ifile)
        with h5py.File(ifile, "r") as f:
            dataar      = np.array(f['trajdata'])
    # flip time axis !!!
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
        - uses numpy and functions readtraj, gridder
    RETURN
        - npart (nlat x nlon) at refdate
    """
    if verbose:
        #print("Reference number of particles: \t" + str(nparticle))
        print(" * Getting reference distribution...")

    ary_npart   = np.zeros(shape=(glat.size,glon.size))
    ary         = readtraj(idate    = refdate,
                           ipath    = "/scratch/gent/vo/000/gvo00090/D2D/data/FLEXPART/era_global/flex2traj_t2/"+str(ryyyy),
                           ifile_base = "global")
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

def random_attribution_p(qtot,iupt,explainp,nmin=1,forc_all=False,verbose=True,veryverbose=False):
  qtot = qtot*1000
  # This is only coded for precipitation as of now
  # with:
  # qtot = specific humidity
  # iupt = identified uptake locations
  # explainp = none(default)/full analogue to linear discounting & attribution
  #  - none: maximum is attributed to iupt (can be 100%!)
  #  - max: maximum is attributed to iupt + initial condition (~not explained); enforcing at least one iteration on init. cond.
  #  - full: 100% is attributed to iupt (if possible!)
  # nmin = tuning parameter; ~minimum iterations 
  #  - the higher this value, the more iterations, the uptake locations are covered
  #  - a value of 10 enforces min. 10 iterations
  # forc_all = enforce attribution to all uptake locations (but still random)
  dqdt  = qtot[:-1] - qtot[1:]
  # append initial condition as artificial uptake
  dqdt  = np.append(dqdt,qtot[-1])
  nt    = len(dqdt)
  if explainp=="max":
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
    # enforce attribution to initial cond. if explain==max
    if icount==0 and explainp=="max" and not forc_all:
        ii   = nt-1
    # enfore acttribution to all uptake locations (forc_all==True)    
    elif icount < nupt and forc_all:
        if icount==0 and verbose:
            print("  *** Random attribution with forc_all=True: enforcing at least one attribution to all "+ str(nupt)+ " uptake locations")
        i    = range(nupt)[icount]
        ii   = np.where(pupt==1)[0][i]  # uptake location index
        if veryverbose:
            print("  *** -- enforcing attribution to uptake location " + str(ii))
    else:    
        i    = random.randint(0,nupt-1) # get a random uptake location number
        ii   = np.where(pupt==1)[0][i]  # uptake location index
    # determine maximum attribution for current uptake location and iteration
    try:
        imin    = np.argmin(qtot[1:ii])+1
    except:
        imin    = 1
    iatt        = qtot[imin]-np.sum(dqdt_random[imin:])
    if iatt<0:  # quick fix: this can happen due to precision issues...
        iatt    = 0
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
    # safety check
    if np.any(dqdt_random<0):
        print("\n !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        print(" ABORT: negative values in random attribution... ")
        print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! \n")
        sys.exit()
  # reset for maximum attribution (only needed if veryverbose is set to True)
  if explainp=="max":
      dqdt_random[-1]   = 0
      nupt              -= 1
  if veryverbose:
      print("  *** "+str(icount)+" Iterations for "+str(nupt)+" uptake locations with P={:.4f}".format(dqdt[0])+" g/kg with E2Prandom={:.4f}".format(np.sum(dqdt_random))+ " g/kg (attributed {:.2f}".format(100*np.sum(dqdt_random)/abs(dqdt[0]))+"%).")
  # upscaling to 100% if explain==full  
  if explainp=="full" and maxatt<1:
      explfr    = abs(dqdt[0])/np.sum(dqdt_random)
      dqdt_random *= explfr
      if veryverbose:
        print(" * Upscaling of contributions required...")
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
    ind_lat = np.argmin(np.abs(glat-lat_mid))    # index on grid
    ind_lon = np.argmin(np.abs(glon-lon_mid))    # index on grid
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
    ind_lat = np.argmin(np.abs(glat-mlat))    # index on grid
    ind_lon = np.argmin(np.abs(glon-mlon))    # index on grid
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
    ind_lat = np.argmin(np.abs(glat-mlat))    # index on grid
    ind_lon = np.argmin(np.abs(glon-mlon))    # index on grid
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
    nc_f.source         = "HAMSTER v0.2 ((c) Dominik Schumacher and Jessica Keune)" 
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
    nc_f.source         = "HAMSTER v0.2 ((c) Dominik Schumacher and Jessica Keune)" 
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
   
def checkdim(var):
    # check dimension of variables (has to be consistent) and use 2D, 3D or 4D definitions
    ndims   = len(var.shape)
    if ndims==4:
        mydims  = ('time','level','lat','lon')
    if ndims==3:
        mydims  = ('time','lat','lon')
    if ndims==2:
        mydims  = ('lat','lon')
    return mydims

def writefinalnc(ofile,fdate_seq,udate_seq,glon,glat,
                 Had, Had_Hs,
                 E2P, E2P_Es, E2P_Ps, E2P_EPs,
                 strargs,precision,
                 fwrite_month):
    
    # delete nc file if it is present (avoiding error message)
    try:
        os.remove(ofile)
    except OSError:
        pass
        
    # create netCDF4 instance
    nc_f = nc4.Dataset(ofile,'w', format='NETCDF4')

    ### create dimensions ###
    if fwrite_month:
        nc_f.createDimension('time', 1)
    else:
        nc_f.createDimension('time', len(fdate_seq))
    if not np.any(np.isnan(udate_seq)):    
        nc_f.createDimension('level', len(udate_seq))
    nc_f.createDimension('lat', glat.size)
    nc_f.createDimension('lon', glon.size)

    # create grid + time variables
    times               = nc_f.createVariable('time', 'f8', 'time')
    if not np.any(np.isnan(udate_seq)):    
        utimes              = nc_f.createVariable('level', 'i4', 'level')
    latitudes           = nc_f.createVariable('lat', 'f8', 'lat')
    longitudes          = nc_f.createVariable('lon', 'f8', 'lon')

    # create variables
    heats               = nc_f.createVariable('Had', precision, checkdim(Had))
    heats_Hs            = nc_f.createVariable('Had_Hs', precision, checkdim(Had_Hs))
    evaps               = nc_f.createVariable('E2P', precision, checkdim(E2P))
    evaps_Es            = nc_f.createVariable('E2P_Es', precision, checkdim(E2P_Es))
    evaps_Ps            = nc_f.createVariable('E2P_Ps', precision, checkdim(E2P_Ps))
    evaps_EPs           = nc_f.createVariable('E2P_EPs', precision, checkdim(E2P_EPs))
 
    # set attributes
    nc_f.title          = "Bias-corrected source-sink relationships from FLEXPART"
    nc_f.description    = str(strargs)
    today               = datetime.datetime.now()
    nc_f.history        = "Created " + today.strftime("%d/%m/%Y %H:%M:%S") + " using HAMSTER."
    nc_f.institution    = "Hydro-Climate Extremes Laboratory (H-CEL), Ghent University, Ghent, Belgium"
    nc_f.source         = "HAMSTER v0.2 ((c) Dominik Schumacher and Jessica Keune)" 
    times.units         = 'hours since 1900-01-01 00:00:00'
    times.calendar      = 'Standard' # do NOT use gregorian here!
    if not np.any(np.isnan(udate_seq)):    
        utimes.long_name    = 'Difference between uptake and arrival time, in days'
        utimes.units        = 'day'
    latitudes.units     = 'degrees_north'
    longitudes.units    = 'degrees_east'
    heats.units         = 'W m-2'
    heats.long_name     = 'advected surface sensible heat'
    heats_Hs.units      = 'W m-2'
    heats_Hs.long_name  = 'advected surface sensible heat, H-scaled' # this is garbage, I know
    evaps.units         = 'mm'
    evaps.long_name     = 'evaporation resulting in precipitation'
    evaps_Es.units      = 'mm'
    evaps_Es.long_name  = 'evaporation resulting in precipitation, E-corrected'
    evaps_Ps.units      = 'mm'
    evaps_Ps.long_name  = 'evaporation resulting in precipitation, P-corrected'
    evaps_EPs.units     = 'mm'
    evaps_EPs.long_name = 'evaporation resulting in precipitation, E-and-P-corrected'

    # write data
    if fwrite_month:
        times[:]        = nc4.date2num(fdate_seq[0], times.units, times.calendar)
    else:
        times[:]        = nc4.date2num(fdate_seq, times.units, times.calendar)
    if not np.any(np.isnan(udate_seq)):    
        utimes[:]       = np.arange(-len(udate_seq)+1,1)
    latitudes[:]        = glat
    longitudes[:]       = glon
  
    if fwrite_month:
        heats[:]            = np.nanmean(Had,axis=0,keepdims=True)[:]
        heats_Hs[:]         = np.nanmean(Had_Hs,axis=0,keepdims=True)[:]
        evaps[:]            = np.nansum(E2P,axis=0,keepdims=True)[:]
        evaps_Es[:]         = np.nansum(E2P_Es,axis=0,keepdims=True)[:]
        evaps_Ps[:]         = np.nansum(E2P_Ps,axis=0,keepdims=True)[:]
        evaps_EPs[:]        = np.nansum(E2P_EPs,axis=0,keepdims=True)[:]
    else:
        heats[:]            = Had[:]
        heats_Hs[:]         = Had_Hs[:]
        evaps[:]            = E2P[:]
        evaps_Es[:]         = E2P_Es[:]
        evaps_Ps[:]         = E2P_Ps[:]
        evaps_EPs[:]        = E2P_EPs[:]

    myshape=nc_f['E2P'].shape
    # close file
    nc_f.close()
    
    # print info
    print("\n * Created and wrote to file: "+ofile+" of dimension "+str(myshape)+" !\n")

def append2csv(filename, listvals):
    # Open file in append mode
    with open(filename, 'a+', newline='\n') as write_obj:
        csv_writer = csv.writer(write_obj, delimiter='\t', lineterminator='\n')
        csv_writer.writerow(listvals)

def preloop(datetime_bgn, uptdatetime_bgn, timestep,
            ipath, ifile_base,
            ryyyy,
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
        ary = readtraj( idate    = predatetime_seq[pix],
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
    

def convert2daily(xar,ftime,fagg="mean"):

    if ftime[0].hour in [0, 6, 12, 18]:
        ## simple fix, subtract 3 hours
        ftime   = np.asarray([t - datetime.timedelta(hours=3) for t in ftime])
        dates   = cal2date(ftime)
    elif ftime[0].hour in [3, 9, 15, 21]:
        ## NOTE: this is the new norm! retain "old style" for now, though    
        dates   = cal2date(ftime)
        pass
    dtime   = np.unique(dates)

    xtot = np.zeros(shape=(dtime.size, xar.shape[1], xar.shape[2]))

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
    return(alpha)

def calc_sourcebcf(ref,diag,tscale='daily'):
    # define all positive
    if np.all(ref<=0):
        ref     = -ref
    if np.all(diag<=0):
        diag    = -diag
    # set 0 to nan to avoid 1e300 values
    diag[diag==0]=np.nan
    # calculate bias correction factor
    if tscale=="daily":
        alpha   = np.nan_to_num(np.divide(ref,diag))
    if tscale=="monthly":
        ref_sum = np.nansum(ref,axis=(0))
        diag_sum= np.nansum(diag,axis=(0))
        alpha   = np.nan_to_num(np.divide(ref,diag))
        #alpha   = np.repeat(alpha_H[np.newaxis,:,:,:], 28, axis=0)
    alpha[alpha==np.inf]    = 0
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

def reduce4Darray(myarray,veryverbose):
    myshape     = myarray.shape
    bwtimesteps = myshape[1]-myshape[0]+1
    myarray_red = np.empty(shape=(myshape[0],bwtimesteps,myshape[2],myshape[3]))
    if veryverbose:
        print(" * Reducing array from "+str(myshape)+" to "+str(myarray_red.shape))
    for iat in range(myshape[0]):
        myarray_red[iat,:,:,:]  = myarray[iat,(iat):(iat+bwtimesteps),:,:]
    return(myarray_red)

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

def check_attributedp(pdiag,pattr,veryverbose):
    # define all positive
    if np.all(pdiag<=0):
        pdiag   = -pdiag
    if np.all(pattr<=0):
        pattr   = -pattr
    printwarning= False
    returnval   = False
    pdiag_sum   = np.nansum(pdiag,axis=(1))
    pattr_sum   = np.nansum(pattr[:,:,:,:],axis=(1,2,3))
    if round(np.nansum(pdiag_sum),4) != round(np.nansum(pattr_sum),4):
        print("   --- WARNING: total precipitation from 01_diagnosis and 02_attribution differ")
        print(" \t --- Absolute difference: {:.2f}".format(np.nansum(pdiag_sum)-np.nansum(pattr_sum))+" m3")
        printwarning = True
    if np.any(pdiag_sum-pattr_sum != 0):
        print("   --- WARNING: daily precipitation from 01_diagnosis and 02_attribution differ")
        if veryverbose:
            ndiffs  = len(np.where(pdiag_sum-pattr_sum !=0)[0])
            print(" \t --- "+str(ndiffs)+" days have different precipitation sums. ")
            ndiffs  = len(np.where(pdiag_sum>pattr_sum)[0])
            print(" \t --- "+str(ndiffs)+" days have P(01_diagnosis) > P(02_attribution)")
            ndiffs  = len(np.where(pdiag_sum<pattr_sum)[0])
            print(" \t --- "+str(ndiffs)+" days have P(01_diagnosis) < P(02_attribution)")
        printwarning = True
    if printwarning: 
        print("   --- ATTENTION: Using 02_attribution data for bias correction for consistency.")
        returnval   = True
    #print(pattr_sum)
    #print(np.nansum(pattr_sum))
    #print(pdiag_sum)
    #print(np.nansum(pdiag_sum))
    #print(pdiag_sum-pattr_sum)
    #print(np.nan_to_num(100*(pdiag_sum-pattr_sum)/pattr_sum))
    return(returnval)

def needmonthlyp(pdiag,pref):
    returnval   = False
    # define all positive
    if np.all(pref<=0):
        pref    = -pref
    if np.all(pdiag<=0):
        pdiag   = -pdiag
    tocheck     = np.where(pref>0)[0]
    if np.any(pdiag[tocheck]==0):
        print("   --- WARNING: daily bias correction of precipitation not possible.")
        ndiffs  = len(np.where(pdiag[tocheck]==0)[0])
        print(" \t --- "+str(ndiffs)+" days have no diagnosed but observed precipitation.")
        print("   --- ATTENTION: Using monthly precipitation for bias correction for consistency.")
        returnval   = True
    return(returnval)

def writewarning(wfile):
    with open(wfile,'w') as ifile:
        writer  = csv.writer(ifile, delimiter='\t', lineterminator='\n',                              quoting = csv.QUOTE_NONE, quotechar='',)
        writer.writerow(["WARNING: you're writing out daily data, but (additional) monthly bias correction was performed. Your daily data is thus not representative."])
    print("\n WARNING! \n See: "+wfile+ " ! \n")


def writestats_03(sfile,Pref,P_E2P,P_E2P_Escaled,P_E2P_Pscaled,P_E2P_EPscaled,Had,Had_scaled,xla,xlo,ibgn):
    with open(sfile,'w') as ifile:
        writer  = csv.writer(ifile, delimiter='\t', lineterminator='\n',quoting = csv.QUOTE_NONE, quotechar='',)
        writer.writerow(["* - STATISTICS: "])
        ndays       = Pref[ibgn:,:,:].shape[0]
        writer.writerow(["   --- # DAYS EVALUATED:              {:.0f}".format(ndays)])
        writer.writerow([" "])
        writer.writerow(["* - HEAT ADVECTION STATISTICS: "])
        writer.writerow(["   --- Had_unscaled [W m-2]:          {:.2f}".format(np.nanmean(np.nansum(Had,axis=(1,2))))])
        writer.writerow(["   --- Had_Hscaled [W m-2]:           {:.2f}".format(np.nanmean(np.nansum(Had_scaled,axis=(1,2))))])
        writer.writerow([" "])
        writer.writerow(["* - PRECIPITATION STATISTICS: "])
        writer.writerow(["   --- P_REFERENCE [m3]:              {:.2f}".format(np.nansum(Pref[ibgn:,xla,xlo]))])
        writer.writerow(["   --- P_E2P_unscaled [m3]:           {:.2f}".format(np.nansum(P_E2P))])
        writer.writerow(["   --- P_E2P_Escaled [m3]:            {:.2f}".format(np.nansum(P_E2P_Escaled))])
        writer.writerow(["   --- P_E2P_Pscaled [m3]:            {:.2f}".format(np.nansum(P_E2P_Pscaled))])
        writer.writerow(["   --- P_E2P_EPscaled [m3]:           {:.2f}".format(np.nansum(P_E2P_EPscaled))])
        # some contingency table statistics... 
        writer.writerow([" "])
        writer.writerow(["* - CONTINGENCY TABLE SCORES (PRECIPITATION):"])
        pref_sum    = np.nansum(Pref[ibgn:,xla,xlo],axis=(1))
        pdiag_sum   = np.nansum(P_E2P,axis=(1,2))
        myctab      = contingency_table(pref_sum,pdiag_sum,thresh=0)
        myscores    = calc_ctab_measures(myctab)
        writer.writerow(["   --- DAYS OF FALSE ALARMS:        {:.0f}".format(myctab["b"])])
        writer.writerow(["   --- DAYS OF MISSES:              {:.0f}".format(myctab["c"])])
        writer.writerow(["   --- DAYS OF HITS:                {:.0f}".format(myctab["a"])])
        writer.writerow(["   --- DAYS OF CORRECT NEGATIVES:   {:.0f}".format(myctab["d"])])
        writer.writerow(["   --- SUCCESS RATIO:               {:.2f}".format(myscores["sr"])])
        writer.writerow(["   --- FALSE ALARM RATIO:           {:.2f}".format(myscores["far"])])
        writer.writerow(["   --- FREQUENCY BIAS:              {:.2f}".format(myscores["fbias"])])
        writer.writerow(["   --- PROB. OF DETECTION:          {:.2f}".format(myscores["pod"])])
        writer.writerow(["   --- PROB. OF FALSE DETECTION:    {:.2f}".format(myscores["pofd"])])
        writer.writerow(["   --- PEIRCE'S SKILL SCORE:        {:.2f}".format(myscores["pss"])])


def contingency_table(ref,mod,thresh=0):
    # creates a contingency table based on 1D np.arrays
    ieventobs   = (np.where(ref>thresh)[0])
    ineventobs  = (np.where(ref<=thresh)[0])
    a           = len(np.where(mod[ieventobs]>thresh)[0])     # hits
    b           = len(np.where(mod[ineventobs]>thresh)[0])    # false alarms
    c           = len(np.where(mod[ieventobs]<=thresh)[0])    # misses
    d           = len(np.where(mod[ineventobs]<=thresh)[0])   # correct negatives
    return({"a":a,"b":b,"c":c,"d":d})

def try_div(x,y):
    try: return x/y
    except ZeroDivisionError: return 0

def calc_ctab_measures(cdict):
    # calculates common contingency table scores
    # scores following definitions from https://www.cawcr.gov.au/projects/verification/
    a           = cdict["a"]    # hits
    b           = cdict["b"]    # false alarms
    c           = cdict["c"]    # misses
    d           = cdict["d"]    # correct negatives
    # calculate scores
    acc         = try_div(a+d,a+b+c+d)          # accuracy
    far         = try_div(b,a+b)                # false alarm ratio
    fbias       = try_div(a+b,a+c)              # frequency bias
    pod         = try_div(a,a+c)                # probability of detection (hit rate)
    pofd        = try_div(b,b+d)                # probability of false detection (false alarm rate)
    sr          = try_div(a,a+b)                # success ratio
    ts          = try_div(a,a+c+b)              # threat score (critical success index)
    a_random    = try_div((a+c)*(a+b),a+b+c+d)
    ets         = try_div((a-a_random),(a+b+c+a_random)) # equitable threat score (gilbert skill score)
    pss         = pod-pofd                      # peirce's skill score (true skill statistic)
    odr         = try_div(a*d,c*b)              # odd's ratio      
    return({"acc":acc,"far":far,"fbias":fbias,"pod":pod,"pofd":pofd,"sr":sr,"pss":pss,"odr":odr})

def writestats_02(statfile,tneval,tnnevala,tnevalh,tnnevalh,tnnevalm,tnevalp,tnnevalp,patt,psum,punatt,pmiss):
    with open(statfile,'w') as sfile:
        writer=csv.writer(sfile, delimiter='\t', lineterminator='\n',quoting = csv.QUOTE_NONE, quotechar='',)
        writer.writerow(["* - PARCEL STATISTICS: "])
        writer.writerow(["   --- TOTAL EVALUATED PARCELS:       " , str(tneval)])
        writer.writerow([" "])
        writer.writerow(["   --- # PARCELS ARRIVING INSIDE MASK:" , str(tneval-tnnevala)])
        if tnnevala!=tneval:
            writer.writerow(["   --- # PARCELS EVAL. FOR HEAT-ADV:  " , str(tnevalh)+" ({:.2f}".format(100*tnevalh/(tneval-tnnevala))+"%)"])
        if tnevalh!=0:
            writer.writerow(["   ----- WITHOUT UPTAKES IN THE TRAJ: " , str(tnnevalh)+" ({:.2f}".format(100*tnnevalh/(tnevalh))+"%)"])
        writer.writerow([" "])
        writer.writerow(["   --- # PARCELS MIDPOINT INSIDE MASK:" , str(tneval-tnnevalm)])
        if tnnevalm!=tneval:
            writer.writerow(["   --- # PARCELS EVAL. FOR PRECIP:    " , str(tnevalp)+" ({:.2f}".format(100*tnevalp/(tneval-tnnevalm))+"%)"])
        if tnevalp!=0:
            writer.writerow(["   ----- WITHOUT UPTAKES IN THE TRAJ: " , str(tnnevalp)+" ({:.2f}".format(100*tnnevalp/(tnevalp))+"%)"])
        writer.writerow([" "])
        if psum!=0:
            writer.writerow([" * - PRECIPITATION STATISTICS: "])
            writer.writerow(["   --- ATTRIBUTED FRACTION:             {:.2f}".format(patt/psum)])
            writer.writerow(["   --- UNATTRIBUTED FRACTION (TRAJEC):  {:.2f}".format(punatt/psum)])
            writer.writerow(["   --- UNATTRIBUTED FRACTION (NO-UPT):  {:.2f}".format(pmiss/psum)])

def mask3darray(xarray,xla,xlo):
    marray  = np.zeros(shape=xarray.shape)
    marray[:,xla,xlo] = xarray[:,xla,xlo]
    return(marray)

def writedebugnc(ofile,fdate_seq,udate_seq,glon,glat,mask,
                 Pref,Pdiag,Pattr,Pattr_Es,Pattr_Ps,Pattr_EPs,
                 frac_E2P,
                 frac_Had,
                 alpha_P,alpha_P_Ecorrected,alpha_P_res,
                 alpha_E,alpha_H,
                 strargs,precision):
    
    Prefsum     = np.nansum(Pref,axis=(1,2))
    Pdiagsum    = np.nansum(Pdiag,axis=(1,2))
    Pattrsum    = np.nansum(Pattr,axis=(1,2))
    Pattrsum_Es = np.nansum(Pattr_Es,axis=(1,2))
    Pattrsum_Ps = np.nansum(Pattr_Ps,axis=(1,2))
    Pattrsum_EPs= np.nansum(Pattr_EPs,axis=(1,2))
    malpha_H    = np.max(alpha_H[:,:,:],axis=(1,2))
    malpha_E    = np.max(alpha_E[:,:,:],axis=(1,2))

    # delete nc file if it is present (avoiding error message)
    try:
        os.remove(ofile)
    except OSError:
        pass

    # create netCDF4 instance
    nc_f = nc4.Dataset(ofile,'w', format='NETCDF4')

    ### create dimensions ###
    nc_f.createDimension('time', len(fdate_seq))
    nc_f.createDimension('uptaketime', len(udate_seq))
    nc_f.createDimension('lat', glat.size)
    nc_f.createDimension('lon', glon.size)

    # create variables
    times               = nc_f.createVariable('time', 'f8', 'time')
    utimes              = nc_f.createVariable('uptaketime', 'f8', 'uptaketime')
    latitudes           = nc_f.createVariable('lat', 'f8', 'lat')
    longitudes          = nc_f.createVariable('lon', 'f8', 'lon')
    # Variables
    nc_mask             = nc_f.createVariable('mask', 'i4', ('lat','lon'))
    nc_pref             = nc_f.createVariable('Pref', precision, ('time','lat','lon'))
    nc_pdiag            = nc_f.createVariable('Pdiag', precision, ('time','lat','lon'))
    nc_pattr            = nc_f.createVariable('Pattr', precision, ('time','lat','lon'))
    nc_prefs            = nc_f.createVariable('Pref_sum', precision, ('time'))
    nc_pdiags           = nc_f.createVariable('Pdiag_sum', precision, ('time'))
    nc_pattrs           = nc_f.createVariable('Pattr_sum', precision, ('time'))
    nc_pattrs_es        = nc_f.createVariable('Pattr_Es_sum', precision, ('time'))
    nc_pattrs_ps        = nc_f.createVariable('Pattr_Ps_sum', precision, ('time'))
    nc_pattrs_eps       = nc_f.createVariable('Pattr_EPs_sum', precision, ('time'))
    nc_alphap           = nc_f.createVariable('alpha_P',precision,('time'))
    nc_alphap_ebc       = nc_f.createVariable('alpha_P_Ecorrected',precision,('time'))
    nc_alphap_res       = nc_f.createVariable('alpha_P_res',precision,('time'))
    nc_alphae           = nc_f.createVariable('alpha_E',precision,('uptaketime','lat','lon'))
    nc_alphah           = nc_f.createVariable('alpha_H',precision,('uptaketime','lat','lon'))
    nc_frace2p          = nc_f.createVariable('frac_E2P',precision,('time','uptaketime','lat','lon'))
    nc_frachad          = nc_f.createVariable('frac_Had',precision,('time','uptaketime','lat','lon'))
    nc_malphae          = nc_f.createVariable('max_alpha_E',precision,('uptaketime'))
    nc_malphah          = nc_f.createVariable('max_alpha_H',precision,('uptaketime'))
 
    # set attributes
    nc_f.title          = "Debug-file from 03_biascorrection (HAMSTER)"
    nc_f.description    = str(strargs)
    today               = datetime.datetime.now()
    nc_f.history        = "Created " + today.strftime("%d/%m/%Y %H:%M:%S") + " using HAMSTER."
    nc_f.institution    = "Hydro-Climate Extremes Laboratory (H-CEL), Ghent University, Ghent, Belgium"
    nc_f.source         = "HAMSTER v0.2 ((c) Dominik Schumacher and Jessica Keune)" 
    times.units         = 'hours since 1900-01-01 00:00:00'
    times.calendar      = 'Standard' # do NOT use gregorian here!
    utimes.units        = 'hours since 1900-01-01 00:00:00'
    utimes.calendar     = 'Standard' # do NOT use gregorian here!
    latitudes.units     = 'degrees_north'
    longitudes.units    = 'degrees_east'
    nc_mask.units          = '-'
    nc_pref.units          = 'm3'
    nc_pref.long_name	   = 'reference precipitation'
    nc_pdiag.units         = 'm3'
    nc_pdiag.long_name	   = 'diagnosed precipitation (01_diag)'
    nc_pattr.units         = 'm3'
    nc_pattr.long_name	   = 'attributed precipitation (E2P, 02_attr)'
    nc_prefs.units         = 'm3'
    nc_prefs.long_name	   = 'sum of reference precipitation'
    nc_pdiags.units        = 'm3'
    nc_pdiags.long_name	   = 'sum of diagnosed precipitation (01_diag)'
    nc_pattrs.units        = 'm3'
    nc_pattrs.long_name	   = 'sum of attributed precipitation (E2P, 02_attr)'
    nc_pattrs_es.units     = 'm3'
    nc_pattrs_es.long_name = 'sum of attributed precipitation (E2P_Escaled, 02_attr)'
    nc_pattrs_ps.units     = 'm3'
    nc_pattrs_ps.long_name = 'sum of attributed precipitation (E2P_Pscaled, 02_attr)'
    nc_pattrs_eps.units    = 'm3'
    nc_pattrs_eps.long_name= 'sum of attributed precipitation (E2P_EPscaled, 02_attr)'
    nc_alphap.units        = '-'
    nc_alphap.long_name	   = 'alpha_P'
    nc_alphap_ebc.units    = '-'
    nc_alphap_ebc.long_name= 'alpha_P_Ecorrected'
    nc_alphap_res.units    = '-'
    nc_alphap_res.long_name= 'alpha_P_res'
    nc_alphae.units        = '-'
    nc_alphae.long_name    = 'alpha_E'
    nc_alphah.units        = '-'
    nc_alphah.long_name    = 'alpha_H'
    nc_malphae.units       = '-'
    nc_malphae.long_name   = 'maximum alpha_E'
    nc_malphah.units       = '-'
    nc_malphah.long_name   = 'maximum alpha_H'
    nc_frace2p.units       = '-'
    nc_frace2p.long_name   = 'frac_E2P'
    nc_frachad.units       = '-'
    nc_frachad.long_name   = 'frac_Had'

    # write data
    times[:]            = nc4.date2num(fdate_seq, times.units, times.calendar)
    utimes[:]           = nc4.date2num(udate_seq, utimes.units, utimes.calendar)
    latitudes[:]        = glat
    longitudes[:]       = glon
    nc_pref[:]          = Pref[:]
    nc_mask[:]          = mask[:]
    nc_pdiag[:]         = Pdiag[:]
    nc_pattr[:]         = Pattr[:]
    nc_prefs[:]         = Prefsum[:]
    nc_pdiags[:]        = Pdiagsum[:]
    nc_pattrs[:]        = Pattrsum[:]
    nc_pattrs_es[:]     = Pattrsum_Es[:]
    nc_pattrs_ps[:]     = Pattrsum_Ps[:]
    nc_pattrs_eps[:]    = Pattrsum_EPs[:]
    nc_alphap[:]        = alpha_P[:]
    nc_alphap_ebc[:]    = alpha_P_Ecorrected[:]
    nc_alphap_res[:]    = alpha_P_res[:]
    nc_alphae[:]        = alpha_E[:]
    nc_alphah[:]        = alpha_H[:]
    nc_malphae[:]       = malpha_E[:]
    nc_malphah[:]       = malpha_H[:]
    nc_frace2p[:]       = frac_E2P[:]
    nc_frachad[:]       = frac_Had[:]

    # close file
    nc_f.close()

    # print info
    print("\n * Created and wrote to file: "+ofile+" !")

def maskbymaskval(mask,maskval):
    mymask  = np.copy(mask)
    mymask[np.where(mask!=maskval)]=0
    return(mymask)

def calc_sinkbcf(ref,att,tscale='daily'):
    # define all positive
    if np.all(ref<=0):
        ref     = -ref
    if np.all(att<=0):
        att     = -att
    if tscale=='daily':
        tref    = np.nansum(ref,axis=tuple(range(1,ref.ndim)))
        tatt    = np.nansum(att,axis=tuple(range(1,att.ndim)))
        alpha   = tref/tatt
        alpha[alpha==np.inf]    = 0
        return(alpha)
    if tscale=='monthly':
        tref    = np.nansum(ref)
        tatt    = np.nansum(att)
        alpha   = np.repeat(tref/tatt,ref.shape[0])
        alpha[alpha==np.inf]    = 0
        return(alpha)

def checkpsum(ref,att,verbose):
    if round(np.nansum(ref),4) != round(np.nansum(att),4):
        ident=False
        if verbose:
            print("       !!! Attributed precipitation does not match reference precipitation!")
            print("           * attributed P: "+str(round(np.nansum(ref),4))+" m3")
            print("           * reference P:  "+str(round(np.nansum(att),4))+" m3")
    else:
        ident=True
    return(ident)

def consistencycheck(attr,diag,bcscale,debug):
    if bcscale=="daily":
        aggattr = np.nansum(attr,axis=0)
    if bcscale=="monthly":
        aggattr = np.nansum(attr,axis=(0,1))
        diag    = np.nansum(diag,axis=(0))
    # calculate fractions
    frac    = np.divide(aggattr,diag)
    # print warnings
    if np.any(frac>1.0001) or np.any(np.isinf(frac)):
        print(" \n  \t !!! WARNING: attribution exceeds diagnosis !!!")
        print(" \t !!!          ---> CHECK YOUR DATA !!!")
        print(" \t !!!          ---> Maximum fraction: " + str(np.max(np.nan_to_num(frac))))
        print(" \t !!!          ---> Number of exceedances: " + str(len(frac[np.where(frac>1.0001)])+len(frac[np.where(np.isinf(frac))]))+" \n")
        if debug:
            print(" \t !!!          ---> frac > 1.001 at: "+ str(np.where(frac>1.001)))
            print(" \t !!!          ----------> attr: "+ str(aggattr[np.where(frac>1.001)]))
            print(" \t !!!          ----------> diag: "+ str((aggattr/frac)[np.where(frac>1.001)]))

#################################################################################################

def f2t_read_partposit(ifile, maxn=3e6, verbose=False):
    """
    @action: reads binary outputs from FLEXPART
    @input:  partposit_DATE.gz
    @output: returns a numpy array of dimension nparcels x 13
    @author: Jessica Keune 06/2020
    #modified: Dominik Schumacher, 06/2020 ---> do use pid!
    #modified: Jessica Keune, 10/2020 --> speeeeeedup!!! 
    # ATTN: hardcoded for 60 bytes / parcel from FP-ERA-INT
    """
    nbytes_per_parcel = (8+4+12*4)
    with gzip.open(ifile, 'rb') as strm:
        # skip header
        _       = strm.read(4) # dummy
        _       = struct.unpack('i', strm.read(4))[0] # time
        # grep full binary data set (ATTN: 60 bytes for FP-ERA-Int hardcoded)
        tdata   = strm.read(int(maxn)*nbytes_per_parcel)
        # get number of parcels from length of tdata
        nparc   = math.floor(len(tdata)/(nbytes_per_parcel))
        # decode binary data
        pdata   = struct.unpack((nparc)*'2fi3fi8f', tdata[0:((nparc)*nbytes_per_parcel)])
        flist   = list(pdata)
    strm.close()
    pdata   = np.reshape(flist, newshape=(nparc,15))[:,2:]
    # remove last line if data is bs (pid = -99999)
    if np.any(pdata[:,0]<0):
        pdata = np.delete(pdata, np.where(pdata[:,0]<0), axis=0)
    return(pdata)

def maskgrabber(maskfile, maskvar='mask', latvar='lat', lonvar='lon'):
    # load
    with nc4.Dataset(maskfile, mode='r') as f:
        mask = np.asarray(f[maskvar][:])
        lat = np.asarray(f[latvar][:])
        lon = np.asarray(f[lonvar][:])
    # check if 2dimensional; necessary for ERA-I lsm mask
    if len(mask.shape) == 3:
        mask = mask[0,:,:]
    # lats check (order irrelevant, just must be within [-90,90])
    if not (lat.min()==-90 or lat.max()==90):
        return(None)
    # lons check
    if lon.min()==-180 and lon.max()<180:
        pass
    elif np.array_equal(lon, np.arange(0,360)):
        mask, lon = ncdf_lon360to180(mask, lon, 1)
    else:
        # this case is not implemented
        return(None)
    return(mask, lat, lon)

def ncdf_lon360to180(ary, lons, lonaxis=1):
    # bring lon axis to front to handle any shape of ary
    ary = np.moveaxis(ary, lonaxis, 0)
    ary_bu  = np.copy(ary)
    lons_bu = np.copy(lons)
    # copy lons & data
    lons[:int(lons.size/2)] = lons_bu[int(lons.size/2):] - 360
    lons[int(lons.size/2):] = lons_bu[:int(lons.size/2)]
    ary[:int(lons.size/2)] = ary_bu[int(lons.size/2):]
    ary[int(lons.size/2):] = ary_bu[:int(lons.size/2)]
    # move axis back to where it was
    ary = np.moveaxis(ary, 0, lonaxis)
    return(ary, lons)

def f2t_timelord(ntraj_d, dt_h, tbgn, tend):
    fulltime = []
    fulltime.append(tbgn - datetime.timedelta(days=ntraj_d, hours=dt_h))
    while fulltime[-1] < tend:
        fulltime.append(fulltime[-1]+datetime.timedelta(hours=dt_h))
    # convert to strings in matching format for partposit files
    fulltime_str = [dft.strftime('%Y%m%d%H%M%S') for dft in fulltime]
    return(fulltime_str)

def f2t_loader(ifile, fixlons=True, fixids=True):
    dummy = f2t_read_partposit(ifile, verbose=False)
    ## fix parcel ID's (ATTN: specific to the global FP-ERA-Interim run!)
    if fixids:
        dummy[:,0] = f2t_fixid(IDs=dummy[:,0], verbose=verbose) # fix IDs
    ## fix longitudes
    if fixlons:
        dummy[:,1][dummy[:,1]>=179.5] -= 360
    ## sort array by parcel ID    
    sorted_dummy = dummy[np.argsort(dummy[:,0]),:]
    return(sorted_dummy)

def f2t_fixid(IDs, verbose, thresidx=1997000 ):
    # get duplicate IDs independent of threshidx...
    #u, c    = np.unique(IDs, return_counts=True)
    #dup     = u[c > 1]
    #print("\t --> Number of duplicates (unconditional): "+ str(len(dup)))
    # add 2e6 to the second duplicate...
    #for i in range(len(dup)):
    #    IDs[np.where(IDs==dup[i])[0][1]] += 2e6
    ## simply shift to indices > 2e6
    IDs[thresidx:][IDs[thresidx:]<(2e6-thresidx)] += 2e6
    if verbose:
        ndupl = np.where(IDs>2e6)[0].size
        if ndupl == 0:
            print("        --> NO duplicates present")
        else:
            print("        --> "+str(ndupl)+" duplicate IDs shifted")
    return(IDs)

def f2t_seeker(array2D, mask, val, lat, lon):
    ## check if anything needs to be done at all
    if mask is None:
        return(array2D[:,0][np.where(~np.isnan(array2D[:,0]))])
    ## first, we search potential candidates using rectangular box
    imlat, imlon = np.where(mask==val)
    lat1 = lat[imlat].min() -0.5
    lat2 = lat[imlat].max() +0.5
    lon1 = lon[imlon].min() -0.5
    lon2 = lon[imlon].max() +0.5
    ## go for it (this is pretty fast)
    idx_inbox = np.where( (array2D[:,1] >= lon1) & (array2D[:,1] <= lon2) &
                          (array2D[:,2] >= lat1) & (array2D[:,2] <= lat2) )[0]
    ## now check if *really* in mask (slow!)
    idx = []
    for ii in range(idx_inbox.size):
        jdx = idx_inbox[ii]
        if mask[np.argmin(np.abs(lat-array2D[jdx,2])),
                np.argmin(np.abs(lon-array2D[jdx,1]))] == val:
            idx.append(jdx)
    ## finally, return parcel IDs
    pid = array2D[:,0][np.asarray(idx)]
    return(pid)

def f2t_locator(array2D, pid, tstring):
    ## figure out where parcels are (lines may shift b/w files)
    pidx = np.where(np.isin(array2D[:,0], pid, assume_unique=False))[0] # <----- set True ??
    chosen = np.NaN*np.ones(shape=(len(pid), array2D.shape[1]))
    if not pidx.size == len(pid):
        ## ATTN: this needs to be adjusted for other runs...
        print("---- INFO: not all parcels present in file --> partposit_"+tstring)
        idx_pidok = np.where(np.isin(pid, array2D[pidx,0], assume_unique=False))[0] # <----- set True ??
        chosen[idx_pidok,:] = array2D[pidx,:]
    else:
        chosen[:,:] = array2D[pidx,:]
    return(chosen)

def f2t_constructor(array3D, pid, time_str):
    ## sloppy check
    if not array3D.shape[0] == len(time_str):
        raise IndexError('time_str must match time dimension of array3D!')
    ## prepare large array, loop thru
    trajs = np.empty(shape=(array3D.shape[0], pid.size, array3D.shape[2]))
    for ii in range(array3D.shape[0]):
        ## call locator
        trajs[ii,:,:] = f2t_locator(array2D=array3D[ii,:,:], pid=pid, tstring=time_str[ii])
    return(trajs)

def f2t_saver(odata, outdir, fout, tstring):
    with h5py.File(outdir+'/'+fout+'_'+tstring+'.h5', 'w') as f:
        f.create_dataset("trajdata", data=odata)

def f2t_establisher(partdir, selvars, time_str, ryyyy, mask, maskval, mlat, mlon,
                    outdir, fout, verbose):
    ##-- 1.) load em files
    data = np.empty(shape=(len(time_str),2000001,selvars.size))
    for ii in range(len(time_str)):
         if verbose: print("       "+time_str[ii][:-4], end='')
         ifile = partdir+'/partposit_'+time_str[ii]+'.gz'
         dummy = f2t_loader(ifile, fixlons=True, fixids=True)[:,selvars] # load
         data[ii,:dummy.shape[0]] = dummy[:] # fill only where data available
         data[ii,dummy.shape[0]:] = np.NaN

    ##-- 2.) find IDs within mask
    if verbose: print("       searching IDs", end='')
    pid_inmask = f2t_seeker(array2D=data[-1,:,:],
                            mask=mask, val=maskval, lat=mlat, lon=mlon)

    ##-- 3.) grab data for IDs
    if verbose: print(" | grabbing data for "+str(pid_inmask.size)+" IDs", end='')
    trajs = f2t_constructor(array3D=data, pid=pid_inmask, time_str=time_str)

    ##--4.) save
    if verbose: print(" | writing to file", end='')
    f2t_saver(odata=trajs, outdir=outdir, fout=fout, tstring=time_str[-1][:-4])

    ##--5.) return data & trajs arrays (needed for next files)
    return(data, trajs)

def f2t_ascender(data, partdir, selvars, ryyyy, time_str, mask, maskval,
                 mlat, mlon, outdir, fout, verbose):
    ##--1.) move old data & fill current step with new data  
    if verbose: print("\n      ", time_str[-1][:-4], end='')
    # use loop to avoid RAM spike here
    for ii in range(len(time_str)-1):
        data[ii,:,:] = data[ii+1,:,:]
    # load new data | rely on dummy variable
    ifile = partdir+'/partposit_'+time_str[-1]+'.gz'
    dummy = f2t_loader(ifile, fixlons=True, fixids=True)[:,selvars]
    # insert new data, use NaN for rest
    data[-1,:dummy.shape[0]] = dummy[:]
    data[-1,dummy.shape[0]:] = np.NaN

    ##--2.) find all IDs
    if verbose: print("       searching IDs", end='')
    pid_inmask = f2t_seeker(array2D=data[-1,:,:],
                            mask=mask, val=maskval, lat=mlat, lon=mlon)

    ##--3.) construct new trajectories (trajs recycling option has been removed)
    if verbose: print(" | constructing trajs for "+str(pid_inmask.size)+" IDs from scratch", end='')
    trajs = f2t_constructor(array3D=data, pid=pid_inmask, time_str=time_str[:])

    ##--4.) save
    if verbose: print(" | writing to file", end='')
    f2t_saver(odata=trajs, outdir=outdir, fout=fout, tstring=time_str[-1][:-4]) # omit mins & secs

    ##--5.) return updated data & trajs arrays
    return(data, trajs)

def checknan(x):
    x[x>=9.9e+36]=np.nan
    return(x)

def whereinmask(mask, maskval, masklat, masklon, trajlat, trajlon):
    ## first, we search potential candidates using rectangular box
    imlat, imlon = np.where(mask==maskval)
    lat1 = masklat[imlat].min() -0.5
    lat2 = masklat[imlat].max() +0.5
    lon1 = masklon[imlon].min() -0.5
    lon2 = masklon[imlon].max() +0.5
    ## go for it (this is pretty fast)
    idx_inbox = np.where( (trajlon >= lon1) & (trajlon <= lon2) &
                          (trajlat >= lat1) & (trajlat <= lat2) )[0]
    ## now check if *really* in mask (slow!)
    idx = []
    for ii in range(idx_inbox.size):
        jdx = idx_inbox[ii]
        if mask[np.argmin(np.abs(masklat-trajlat[jdx])),
                np.argmin(np.abs(masklon-trajlon[jdx]))] == maskval:
            idx.append(jdx)
    ## finally, return indices for which traj in mask
    return(np.asarray(idx))

def maxlastn(series, n=4):
    maxy = np.zeros(shape=(n,series.size))
    maxy[0,:]   = series[:]
    for ii in range(1,n):
        maxy[ii,:-ii] = series[ii:]
    return(np.max(maxy, axis=0))

def grabhpbl_partposit(ifile):
    dummy   = f2t_loader(ifile, fixlons=True, fixids=True)[:,[0,9]] # 0: id, 9: hpbl
    return(dummy)

def grabmesomehpbl(ipath, fdatetime_beg, tml, verbose):
    extendarchive = []

    if verbose:
        print("\n--------------------------------------------------------------------------------------")
        print("\n ! performing pre-loop to extend trajectory data & achieve advection month-2-month consistency")
        print("\n ! estimating remaining time for pre-loop ...")
        pretic = timeit.default_timer()

    for qq in range(tml+5):
        if verbose and qq==1:
            pretoc = timeit.default_timer()
            mins   = round((tml+5)*(pretoc-pretic)/60, 2)
            if mins < 1.0:
                print("  ---> "+str(mins*60)+" seconds to go, how about some stretching in the meantime?")
            else:
                print("  ---> "+str(mins)+" minutes to go, might want to grab a coffee...")

        # go back from first timestep..
        timestr = (fdatetime_beg - datetime.timedelta(hours=qq*6+6)).strftime('%Y%m%d%H')+'0000'

        # append data
        ifile = ipath+'/partposit_'+timestr+'.gz'
        extendarchive.append(grabhpbl_partposit(ifile))

    return(extendarchive)
