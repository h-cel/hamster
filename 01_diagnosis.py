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
           cpbl_strict,
           refdate,
           fwrite_netcdf,ftimethis,fcc_advanced,fvariable_mass,
           strargs):

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

    ## -- GRID
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

    ## -- WRITE NETCDF OUTPUT (empty, to be filled)
    if fwrite_netcdf:
        writeemptync(ofile,fdate_seq,glon,glat,strargs)

    # set some default thresholds
    cprec_dqv    = default_thresholds(cprec_dqv) 

    # read in reference distribution of parcels
    if fvariable_mass:
        ary_rnpart   = get_refnpart(refdate=refdate, ryyyy=ryyyy, glon=glon, glat=glat)
    
    ## ------- LOOP OVER DATES (FILES) 
    if verbose:
        print("\n=== \t Start main program...\n")

    # pre-allocate arrays
    ary_heat     = np.zeros(shape=(glat.size,glon.size))
    ary_evap     = np.zeros(shape=(glat.size,glon.size))
    ary_prec     = np.zeros(shape=(glat.size,glon.size))
    ary_npart    = np.zeros(shape=(glat.size,glon.size))

    for ix in range(ntime):
        
        if verbose:
            print("--------------------------------------------------------------------------------------")
            print("Processing "+str(fdate_seq[ix]))

        ## READ DATE RELATED TRAJECTORIES -> ary is of dimension (ntrajlen x nparticles x nvars)
        ary         = readpom( idate        = date_seq[ix], 
                               ipath        = ipath+"/"+str(ryyyy), 
                               ifile_base   = ifile_base)
        nparticle   = ary.shape[1]
        if verbose:
            print(" TOTAL: " + str(date_seq[ix]) + " has " + str(nparticle) + " parcels")

        ## TEST mode: less parcels
        if mode == "test":
            ntot    = range(nparticle)#range(1000)
        else:
            ntot    = range(nparticle)

        #smalltic = timeit.default_timer()

        ## ------- LOOP OVER PARCELS TO DIAGNOSE P, E, H (and npart) and assign to grid 
        for i in ntot:

            ## get midpoint at the very beginning
            #lat_mid, lon_mid = readmidpoint(ary[:,i,:])
            lat_ind, lon_ind = midpindex(ary[:,i,:],glon=glon,glat=glat)
            ## log number of parcels
            ary_npart[lat_ind,lon_ind] += int(1)

            ## read only necessary parcel information
            qv, temp, ztra, hpbl    = readsparcel(ary[:,i,:])
            dq = parceldiff(qv, 'diff') 
            dT = parceldiff(temp, 'diff')
            
            ##  KAS: KEUNE AND SCHUMACHER
            if tdiagnosis == 'KAS':

                # precipitation
                if (dq < cprec_dqv):
                    pres            = readpres(ary[:,i,:])
                    if ( q2rh(qv[0], pres[0], temp[0]) > cprec_rh  and
                         q2rh(qv[1], pres[1], temp[1]) > cprec_rh ):
                         ary_prec[lat_ind,lon_ind] += dq

                # evaporation and sensible heat 
                if ( checkpbl(cpbl_strict,ztra,hpbl,cevap_hgt) or checkpbl(cpbl_strict,ztra,hpbl,cheat_hgt) ):
                    if ( dq > 0 or dT > 0):
                        pres                = readpres(ary[:,i,:])
                        pottemp             = readpottemp(ary[:,i,:])
                        dTH                 = parceldiff(pottemp, 'diff')

                        # sensible heat
                        if ( dT > 0 and checkpbl(cpbl_strict,ztra,hpbl,cheat_hgt)):
                            dqmax = cheat_cc * dqsdT(p_hPa=pres[1]/1e2, T_degC=temp[1]-TREF)
                            if fcc_advanced:
                                if  ((dTH > cheat_dtemp) and ( (dT > 0 and abs(dq) < dT*dqmax) or (dT < 0 and abs(dq) < dTHdqmax) )):
                                    ary_heat[lat_ind,lon_ind] += dTH
                            else:
                                if  ((dTH > cheat_dtemp) and abs(dq) < dTH*dqmax ):
                                    ary_heat[lat_ind,lon_ind] += dTH

                        # evaporation
                        if ( dq > 0 and checkpbl(cpbl_strict,ztra,hpbl,cevap_hgt)):
                            epottemp        = readepottemp(ary[:,i,:])
                            dTHe            = parceldiff(epottemp, 'diff')
                            dTmax           = cevap_cc*dTdqs(p_hPa=pres[1]/1e2, q_kgkg=qv[1])
                            if fcc_advanced:
                                if ( (dTHe - dTH) > cheat_dtemp and ( (dT > 0 and dT < dq*dTmax) or (dT < 0 and abs(dTH) < dq*dTmax) )):
                                    ary_evap[lat_ind,lon_ind] += dq
                            else:
                                if ( (dTHe - dTH) > cheat_dtemp and abs(dTH) < dq*dTmax ):
                                    ary_evap[lat_ind,lon_ind] += dq


            ## SOD: SODEMANN ET AL., 2008
            elif tdiagnosis == 'SOD':
        
                # note: not optimized wrt runtime
                pres            = readpres(ary[:,i,:])
                pottemp         = readpottemp(ary[:,i,:])
                dTH             = parceldiff(pottemp, 'diff')

                ## precipitation
                if (dq < 0 and 
                        q2rh((qv[0]+qv[1])/2, (pres[0]+pres[1])/2, (temp[0]+temp[1])/2) > 80):
                    ary_prec[lat_ind,lon_ind] += dq

                ## evaporation
                if (dq > 0.0002 and  
                        (ztra[0]+ztra[1])/2 < 1.5*(hpbl[0]+hpbl[1])/2):
                    ary_evap[lat_ind,lon_ind] += dq
    
                ## sensible heat (not used originally; analogous to evaporation)
                if ((dTH > cheat_dtemp) and 
                        (ztra[0]+ztra[1])/2 < 1.5*(hpbl[0]+hpbl[1])/2):
                    ary_heat[lat_ind,lon_ind] += dTH

            ## SOD2: SODEMANN, 2020; FREMME & SODEMANN, 2019
            elif tdiagnosis == 'SOD2':
        
                # note: not optimized wrt runtime
                pres            = readpres(ary[:,i,:])
                pottemp         = readpottemp(ary[:,i,:])
                dTH             = parceldiff(pottemp, 'diff')

                ## precipitation
                if (dq < 0 and 
                        q2rh((qv[0]+qv[1])/2, (pres[0]+pres[1])/2, (temp[0]+temp[1])/2) > 80):
                    ary_prec[lat_ind,lon_ind] += dq

                ## evaporation
                if (dq > 0.0001): 
                    ary_evap[lat_ind,lon_ind] += dq
    
                ## sensible heat (not used originally; analogous to evaporation)
                if (dTH > cheat_dtemp):
                    ary_heat[lat_ind,lon_ind] += dTH

            ## SAJ: STOHL AND JAMES, 2004
            elif tdiagnosis == 'SAJ':

                ## precipitation AND evaporation
                ary_prec[lat_ind,lon_ind] += dq
                ary_evap[lat_ind,lon_ind] += dq

        #smalltoc = timeit.default_timer()
        #print("=== \t All parcels: ",str(round(smalltoc-smalltic, 2)),"seconds \n")

        # Convert units
        if verbose:
            print(" * Converting units...")
        if tdiagnosis == 'KAS' or tdiagnosis == 'SOD' or tdiagnosis == 'SOD2':
            ary_prec[:,:] = convertunits(ary_prec[:,:], garea, "P")
            ary_evap[:,:] = convertunits(ary_evap[:,:], garea, "E")
            ary_heat[:,:] = convertunits(ary_heat[:,:], garea, "H")
        elif tdiagnosis =='SAJ':
            ary_evap[ary_evap<0]= 0
            ary_prec[ary_prec>0]= 0
            ary_evap[:,:]   = convertunits(ary_evap, garea, "E")
            ary_prec[:,:]   = convertunits(ary_prec, garea, "P")

        # Scale with parcel mass
        if fvariable_mass:
            if verbose: 
                print(" * Applying variable mass...")
            ary_prec[:,:]         = scale_mass(ary_prec[:,:], ary_npart[:,:], ary_rnpart)
            ary_evap[:,:]         = scale_mass(ary_evap[:,:], ary_npart[:,:], ary_rnpart)
            ary_heat[:,:]         = scale_mass(ary_heat[:,:], ary_npart[:,:], ary_rnpart)

        if fwrite_netcdf:
            writenc(ofile,ix,ary_prec[:,:],ary_evap[:,:],ary_heat[:,:],ary_npart[:,:])

        # re-init. arrays
        ary_npart[:,:]  = 0
        ary_prec[:,:]   = 0
        ary_evap[:,:]   = 0
        ary_heat[:,:]   = 0

    if ftimethis:
        megatoc = timeit.default_timer()
        if verbose:
            print("\n=== \t End main program (total runtime so far: ",str(round(megatoc-megatic, 2)),"seconds) \n")

    if verbose:
        if fwrite_netcdf:
            print("\n Successfully written: "+ofile+" !")

        
