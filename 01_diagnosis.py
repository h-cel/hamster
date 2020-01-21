#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MAIN FUNCTIONS FOR 01_diagnosis
"""


############################################################################
#############################    SETTINGS ##################################

def main_diagnosis(
           ryyyy, ayyyy, am,
           ipath, ifile_base, 
           opath, ofile_base,
           mode,
           gres,
           tdiagnosis,
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
    ofilename = str(ofile_base)+str(ryyyy)[-2:]+"_"+str(ayyyy)+"-"+str(am).zfill(2)+".nc"
    ofile     = opath+"/"+ofilename

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
            print(" \t ! output file: \t", opath+"/"+ofilename)
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
                       ipath    = ipath+"/"+str(ryyyy), 
                       ifile_base = ifile_base)
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
            if tdiagnosis == 'KAS':
                hpbl_max    = parceldiff(hpbl, 'max')
                dTH         = parceldiff(pottemp, 'diff')
                dTHe        = parceldiff(epottemp, 'diff')
                #dz          = parceldiff(ztra, 'diff')
                if fcc_advanced:
                    dT          = parceldiff(temp, 'diff')
            elif tdiagnosis == 'SOD':
                hpbl_avg    = parceldiff(hpbl, 'mean')
                dTH         = parceldiff(pottemp, 'diff')

            ## - 2.3) diagnose fluxes

            ## (a) number of parcels
            ary_npart[ix,:,:] += gridder(plon=lons, plat=lats, pval=int(1), glon=glon, glat=glat)

            ##  - 2.3)-KAS: Keune and Schumacher
            if tdiagnosis == 'KAS':

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
            elif tdiagnosis == 'SOD':
         
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
            elif tdiagnosis == 'SAJ':

                ## (b) precipitation
                if ( dq < 0 ):
                    ary_prec[ix,:,:] += gridder(plon=lons, plat=lats, pval=dq, glon=glon, glat=glat)
    
                ## (c) evaporation
                if ( dq > 0 ):
                    ary_evap[ix,:,:] += gridder(plon=lons, plat=lats, pval=dq, glon=glon, glat=glat)

        # Convert units
        if verbose:
            print(" * Converting units...")
        if tdiagnosis == 'KAS' or tdiagnosis == 'SOD':
            ary_prec[ix,:,:] = convertunits(ary_prec[ix,:,:], garea, "P")
            ary_evap[ix,:,:] = convertunits(ary_evap[ix,:,:], garea, "E")
            ary_heat[ix,:,:] = convertunits(ary_heat[ix,:,:], garea, "H")
        elif tdiagnosis =='SAJ':
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
