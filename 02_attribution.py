#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MAIN FUNCTIONS FOR 02_attribution
"""

############################################################################
#############################    SETTINGS ##################################

def main_attribution(
           ryyyy, ayyyy, am,
           ipath, ifile_base, 
           opath, ofile_base,
           mode,
           gres,
           diagnosis,
           ctraj_len,
           cheat_dtemp, # used for E,H,P (if cprec_dqv==None)
           cheat_cc, cevap_cc, # for H, E diagnosis (lower = more strict)
           cevap_hgt, cheat_hgt, # set min ABLh, disabled if 0 
           cprec_dqv, cprec_dtemp, cprec_rh,
           refdate,
           fwrite_netcdf,ftimethis,fcc_advanced,fvariable_mass):

    if fcc_advanced or fvariable_mass:
        raise SystemExit("---- ABORTED: no can do, master")
    
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
    # NOTE: we begin at 06 UTC...
    datetime_bgn    = datetime.datetime.strptime(str(ayyyy)+"-"+str(am).zfill(2)+"-01-06", "%Y-%m-%d-%H")
    datetime_end    = datetime_bgn + relativedelta(months=1)
    timestep        = datetime.timedelta(hours=6)
    datetime_seq    = []
    fdatetime_seq   = []
    idatetime       = datetime_bgn

    # create arrival datetime string & datetime object
    while idatetime < datetime_end:
        datetime_seq.append(idatetime.strftime('%Y%m%d%H'))
        fdatetime_seq.append(idatetime)
        idatetime += timestep

    # add uptake time dimension
    uptdatetime_bgn = datetime_bgn - datetime.timedelta(days=ctraj_len) - datetime.timedelta(hours=3)
    fuptdatetime_seq   = []
    iuptdatetime       = uptdatetime_bgn
    while iuptdatetime < datetime_end - datetime.timedelta(hours=3):
        fuptdatetime_seq.append(iuptdatetime)
        iuptdatetime += timestep
  
    """
    issue: arrays become HUGE if using 6-hourly res for TWO time axes,
           yet particle mass scaling appears to be implemented
           at 6-hourly resolution as of now.....
    """

    # aggregate to daily, NOTE: arrival at 00 UTC means parcel has arrived on prev day    
    fdate_seq = np.unique([fdt.date() for fdt in fdatetime_seq[:-1]]).tolist()
    fuptdate_seq = np.unique([fdt.date() for fdt in fuptdatetime_seq]).tolist()
    # NOTE: better to keep these as lists to maintain consistency

    ntime           = len(fdatetime_seq)
    nupttime        = len(fuptdatetime_seq)

    ndaytime        = len(fdate_seq)
    ndayupttime     = len(fuptdate_seq)

    # TESTMODE
    if mode == "test":
        ntime            = 1
        datetime_seq     = datetime_seq[:ntime]
        fdatetime_seq    = fdatetime_seq[:ntime]
        ndaytime         = 1
        fdate_seq        = fdate_seq[:ndaytime]
        
        nupttime         = 4*ctraj_len + 1
        fuptdatetime_seq = fuptdatetime_seq[:nupttime]
        ndayupttime      = (ctraj_len) + 1
        fuptdate_seq     = fuptdate_seq[:ndayupttime]

    print("PRELIM: ntime, ndaytime, nupttime, ndayupttime =", ntime, ndaytime, nupttime, ndayupttime)

    # traj max len
    tml = nupttime - ntime

    print("tml=", tml)

    ## pre-allocate arrays
    ary_heat     = np.zeros(shape=(ndaytime,ndayupttime,glat.size,glon.size))
    ary_etop     = np.zeros(shape=(ndaytime,ndayupttime,glat.size,glon.size))
    ary_npart    = np.zeros(shape=(ndaytime,ndayupttime,glat.size,glon.size))

    # set some default thresholds
    cprec_dqv    = default_thresholds(cprec_dqv) 
    # read in reference distribution of parcels
    if fvariable_mass:
        ary_rnpart   = get_refnpart(refdate=refdate, ryyyy=ryyyy, glon=glon, glat=glat)

    ## prepare parcel log to handle trajectories properly   
    pIDlogH = -np.ones(shape=2000001).astype(int) 

    ## prepare uptake indices
    upt_idx = np.asarray([floor(x) for x in np.arange(0,nupttime)/4])
    print(upt_idx)

    ## loop over time to read in files
    if verbose:
        print("\n=== \t Start main program...\n")
    for ix in range(ntime):
        if verbose:
                print("--------------------------------------------------------------------------------------")
        print("Processing "+str(fdatetime_seq[ix]))

        ## 1) read in all files associated with data --> ary is of dimension (ntrajlen x nparticles x nvars)
        ary = readpom( idate    = datetime_seq[ix], 
                       ipath    = ipath+"/"+str(ryyyy), 
                       ifile_base = ifile_base)

        nparticle   = ary.shape[1]
        if verbose:
            print(" TOTAL: " + str(datetime_seq[ix]) + " has " + str(nparticle) + " parcels")

        #bar = Bar('Processing', suffix='%(percent)d%%', fill="*")
        if mode == "test":
            ntot    = range(1000)
        else:
            ntot    = range(nparticle)

        # figure out where to store data (on which arriving day)
        arv_idx = np.where(np.asarray(fdate_seq)==fdatetime_seq[ix].date())[0][0]
 
        ## 2) diagnose P, E, H and npart per grid cell
        for i in ntot:

            ## - 2.0) check how far back trajectory should be evaluated
            # NOTE: this can be moved elsewhere...
            # for Hadv (not needed for E2P)
            ID = int(ary[0,i,0])   
            istepH = pIDlogH[ID]
            if istepH < 0:
                ihf_H = 4*ctraj_len + 2
            else:
                ihf_H = min((ix-istepH+1), 4*ctraj_len + 2)
            if ihf_H < 4*ctraj_len + 2: # might as well use tml+2
                print("*", end="")
            """
                 ======> the above needs to be (double-)checked!
			.. and it would be tight to do this with less lines, IMPROVE
            """

            ## - 2.1) read parcel information
            lons, lats, temp, ztra, qv, hpbl, dens, pres, pottemp, epottemp = readparcel(ary[:tml+2,i,:]) ### tml+2 PRELIM

            ## - 2.2) parcel changes / criteria
            dq          = trajparceldiff(qv[:], 'diff') 
            if diagnosis == 'KAS':
                hpbl_max    = parceldiff(hpbl[:], 'max')
                dTH         = trajparceldiff(pottemp[:], 'diff')
                dTHe        = trajparceldiff(epottemp[:], 'diff')
                #dz          = parceldiff(ztra, 'diff')
                #if fcc_advanced:
                #    dT          = parceldiff(temp[:2], 'diff')
            elif diagnosis == 'SOD':
                hpbl_avg    = parceldiff(hpbl[:], 'mean')
                dTH         = trajparceldiff(pottemp[:], 'diff')

            ## - 2.3) diagnose fluxes

            ## (a) number of parcels
            for itj in range(tml+1): # NOT tml+2 !!!
                ary_npart[arv_idx,upt_idx[ix+tml-itj],:,:] += gridder(plon=lons[itj:itj+2], plat=lats[itj:itj+2], pval=int(1), glon=glon, glat=glat)

            ##  - 2.3)-KAS: Keune and Schumacher
            if diagnosis == 'KAS':

                ## (b) E2P, evaporation resulting in precipitation
                if ( dq[0] < cprec_dqv and 
                     q2rh(qv[0], pres[0], temp[0]) > cprec_rh  and
                     q2rh(qv[1], pres[1], temp[1]) > cprec_rh ):


                    # first things first, check if traj falls dry & adapt ihf_E if so
                    # NOTE: to be discussed if we even use this (imo: ditch it), 
                    #       ===> do include for SOD, though
                    ihf_dry = np.where(qv[:]<= 0.00005)[0]
                    if ihf_dry.size>0:
                       ihf_E = np.min(ihf_dry)
                       print("---- PRELIM INFO: E2P traj fell dry, cut off at", ihf_E)
                    else:
                       ihf_E = tml + 2

                    # SHOULD I MAKE USE OF PARCELDIFF instead of PBL_check???
                    in_PBL     = PBL_check(z=ztra[:ihf_E], h=hpbl[:ihf_E], seth=cevap_hgt, diagnosis=diagnosis)                      
                    evap_uptk  = (dTHe[:ihf_E-1] - dTH[:ihf_E-1]) > cheat_dtemp 
                    evap_plaus = abs(dTH[:ihf_E-1]) < cevap_cc * (dq[:ihf_E-1]) * dTdqs(p_hPa=pres[1:ihf_E]/1e2, q_kgkg=qv[1:ihf_E])
                    evap_idx   = np.where(np.logical_and(in_PBL, np.logical_and(evap_uptk, evap_plaus)))[0]
                    
#                    print("evap_idx, arv_idx = ", evap_idx, arv_idx)
                   
                    if evap_idx.size>0:
                        dq_disc = linear_discounter(v=qv[:ihf_E], min_gain=0, min_loss=0)
                        etop = ((qv[0]-qv[1])/qv[1])*dq_disc ### NOTE: might want to swap sign

#                    print("dq=", dq)
#                    print("dq_disc=", dq_disc)

                    for itj in evap_idx: 
#                        print(itj, upt_idx[ix+tml-itj])
                        ary_etop[arv_idx,upt_idx[ix+tml-itj],:,:] += gridder(plon=lons[itj:itj+2], plat=lats[itj:itj+2], pval=etop[itj], glon=glon, glat=glat)

                ## (c) H, surface sensible heat arriving in PBL (or nocturnal layer)
                if ( ztra[0] < np.max(hpbl[:4]) ):

                    in_PBL     = PBL_check(z=ztra[:ihf_H], h=hpbl[:ihf_H], seth=cheat_hgt, diagnosis=diagnosis)
                    heat_uptk  = dTH[:ihf_H-1] > cheat_dtemp
                    heat_plaus = abs(dq[:ihf_H-1]) < cheat_cc * (dTH[:ihf_H-1]) * dqsdT(p_hPa=pres[1:ihf_H]/1e2, T_degC=temp[1:ihf_H]-TREF)
                    heat_idx   = np.where(np.logical_and(in_PBL, np.logical_and(heat_uptk, heat_plaus)))[0]

                    if heat_idx.size>0:
                        dTH_disc = linear_discounter(v=pottemp[:ihf_H], min_gain=0, min_loss=0)

                    for itj in heat_idx:
                        #print(itj, upt_idx[ix+tml-itj])
                        ary_heat[arv_idx,upt_idx[ix+tml-itj],:,:] += gridder(plon=lons[itj:itj+2], plat=lats[itj:itj+2], pval=dTH_disc[itj], glon=glon, glat=glat)

                    # update parcel log
                    pIDlogH[ID] = ix # NOTE: double-check



            ##  - 2.3)-SOD: Sodemann et al., 2008
#            elif diagnosis == 'SOD':
#         
#                ## (b) precipitation
#                if ( dq < 0 and 
#                     q2rh((qv[0]+qv[1])/2, (pres[0]+pres[1])/2, (temp[0]+temp[1])/2) > 80 ):
#                    ary_prec[ix,:,:] += gridder(plon=lons, plat=lats, pval=dq, glon=glon, glat=glat)
#
#                ## (c) evaporation
#                if ( dq > 0.0002 and 
#                     (ztra[0]+ztra[1])/2 <  hpbl_avg
#                   ):
#                    ary_evap[ix,:,:] += gridder(plon=lons, plat=lats, pval=dq, glon=glon, glat=glat)
#    
#                ## (d) sensible heat (not used originally; analogous to evaporation)
#                if ( (dTH > cheat_dtemp) and 
#                    (ztra[0]+ztra[1])/2 <  hpbl_avg
#                   ):
#                    ary_heat[ix,:,:] += gridder(plon=lons, plat=lats, pval=dTH, glon=glon, glat=glat)
#
#            ##  - 2.3)-SAJ: Stohl and James, 2004
#            elif diagnosis == 'SAJ':
#
#                ## (b) precipitation
#                if ( dq < 0 ):
#                    ary_prec[ix,:,:] += gridder(plon=lons, plat=lats, pval=dq, glon=glon, glat=glat)
#    
#                ## (c) evaporation
#                if ( dq > 0 ):
#                    ary_evap[ix,:,:] += gridder(plon=lons, plat=lats, pval=dq, glon=glon, glat=glat)

        # Convert units
        if verbose:
            print(" * Converting units...")
        if diagnosis == 'KAS' or diagnosis == 'SOD':
            ary_etop[ix,:,:,:] = convertunits(ary_etop[ix,:,:,:], garea, "E")
            ary_heat[ix,:,:,:] = convertunits(ary_heat[ix,:,:,:], garea, "H")
            # NOTE: the above works for 3D arrays too thanks to numpy broadcasting!
#        elif diagnosis =='SAJ':
#            # first calculate column sums, assign to E or P according to sign
#            colsum = ary_prec[ix,:,:] + ary_evap[ix,:,:]
#            colsum_pos = np.zeros_like(colsum)
#            colsum_neg = np.zeros_like(colsum)
#            colsum_pos[np.where(colsum>0)] = colsum[np.where(colsum>0)]
#            colsum_neg[np.where(colsum<0)] = colsum[np.where(colsum<0)]
#            ary_evap[ix,:,:] = convertunits(colsum_pos, garea, "E")
#            ary_prec[ix,:,:] = convertunits(colsum_neg, garea, "P")

    # Scale with parcel mass
#    if fvariable_mass:
#        if verbose: 
#            print(" * Applying variable mass...")
#        ary_prec         = scale_mass(ary_prec, ary_npart, ary_rnpart)
#        ary_evap         = scale_mass(ary_evap, ary_npart, ary_rnpart)
#        ary_heat         = scale_mass(ary_heat, ary_npart, ary_rnpart)

    if ftimethis:
        megatoc = timeit.default_timer()
        if verbose:
            print("\n=== \t End main program (total runtime so far: ",str(round(megatoc-megatic, 2)),"seconds) \n")
    
    if fwrite_netcdf:
        writenc4D(ofile,fdate_seq,fuptdate_seq,glon,glat,ary_etop,ary_heat,ary_npart)
