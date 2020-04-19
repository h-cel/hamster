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
           maskfile,
           maskval,
           tdiagnosis,
           ctraj_len,
           cheat_dtemp, # used for E,H,P (if cprec_dqv==None)
           cheat_cc, cevap_cc, # for H, E diagnosis (lower = more strict)
           cevap_hgt, cheat_hgt, # set min ABLh, disabled if 0 | NOTE: to be unified
           cprec_dqv, cprec_dtemp, cprec_rh,
           fjumps,
           fjumpsfull,
           cjumps,
           refdate,
           fwrite_netcdf,
           precision,
           ftimethis,
           fdry,
           fmemento,
           explainp,fupscale,
           fcc_advanced,fvariable_mass,fwritestats,
           strargs):

    # TODO: add missing features
    if fcc_advanced or fvariable_mass:
        raise SystemExit("---- ABORTED: no can do, not implemented!")
 
    #### OUTPUT FILES
    mainpath  = ipath+str(ryyyy)+"/"
    ## main netcdf output
    ofilename = str(ofile_base)+"_attr_r"+str(ryyyy)[-2:]+"_"+str(ayyyy)+"-"+str(am).zfill(2)+".nc"
    ofile     = opath+"/"+ofilename
    ## additional statistic output files (*.csv)
    # monthly statistics
    sfilename = str(ofile_base)+"_attr_r"+str(ryyyy)[-2:]+"_"+str(ayyyy)+"-"+str(am).zfill(2)+"_stats.csv"
    statfile  = opath+"/"+sfilename
    # trajectory-based precipitation statistics
    if fwritestats:
        pfilename = str(ofile_base)+"_attr_r"+str(ryyyy)[-2:]+"_"+str(ayyyy)+"-"+str(am).zfill(2)+"_pattribution.csv"
        pattfile  = opath+"/"+pfilename
        with open(pattfile,'w') as pfile:
                writer=csv.writer(pfile, delimiter='\t', lineterminator='\n',)
                writer.writerow(["DATE", "F_ATT", "F_POT", "P_DQDT"])
    
    #### INPUT FILES
    ## read netcdf mask
    with nc4.Dataset(maskfile) as f:
        mask = f['mask'][:]
        mlat = f['lat'][:]
        mlon = f['lon'][:]

    #### DISCLAIMER
    if verbose:
        disclaimer()
        print("\n PROCESSING: \t", 	ayyyy, "-", str(am).zfill(2))
        print("\n============================================================================================================")
        print(" ! using input path: \t", 	ipath)
        if fwrite_netcdf:
            print(" ! writing netcdf output: \t" +str(fwrite_netcdf) )
            print(" \t ! with grid resolution: \t", str(gres) )
            print(" \t ! output file: \t", opath+"/"+ofilename)
        print(" ! using internal timer: \t" +str(ftimethis) )
        print(" ! using mode: \t" +str(mode))
        print(" ! additional statistics in: \t"+str(statfile))
        if fwritestats:
            print(" ! precipitation statistics in: \t"+str(pattfile))
        print("\n============================================================================================================")
        print("\n============================================================================================================")

    #### SETTINGS
    ## start timer
    if ftimethis:
        megatic = timeit.default_timer()
    
    ## grids
    glon, glat, garea = makegrid(resolution=gres)
    ## Sanity check: is glon/glat equal to mlon/mlat from maskfile?
    if not np.array_equal(glon,mlon) or not np.array_equal(glat,mlat):
        warnings.warn("\n----------------- WARNING: the grid from the maskfile is not identical to the target grid... please check. Proceeding nevertheless. \n")

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

    # aggregate to daily, NOTE: arrival at 00 UTC means parcel has arrived on prev day    
    fdate_seq = np.unique([fdt.date() for fdt in fdatetime_seq[:-1]]).tolist() # omit last dt object (00:00)
    fuptdate_seq = np.unique([fdt.date() for fdt in fuptdatetime_seq]).tolist()
    # keep a copy of datetime.date formatted list for arv_idx further below
    fdateasdate = np.copy(fdate_seq).tolist() # NOTE: using deepcopy instead of np.copy would be more proper
    # convert these datetime.date objects to datetime.datetime objects for netCDF writing
    for idt in range(len(fdate_seq)):
        fdate_seq[idt]    = datetime.datetime(fdate_seq[idt].year, fdate_seq[idt].month, fdate_seq[idt].day)
    for idt in range(len(fuptdate_seq)):
        fuptdate_seq[idt] = datetime.datetime(fuptdate_seq[idt].year, fuptdate_seq[idt].month, fuptdate_seq[idt].day)
    # NOTE: better to keep these as lists to maintain consistency

    # calculate number of time steps, also aggregated to daily resolution
    ntime           = len(fdatetime_seq)
    nupttime        = len(fuptdatetime_seq)
    ndaytime        = len(fdate_seq)
    ndayupttime     = len(fuptdate_seq)

    ## TESTMODE
    if mode == "test":

        ctraj_len_orig   = ctraj_len
        ctraj_len        = 2 # NOTE: changes here must be accompanied by changes in 01_diagnosis!

        ntime            = 4 #NOTE: use multiples of 4 only, else output is not saved 
        datetime_seq     = datetime_seq[(4*ctraj_len):(4*ctraj_len+ntime)]
        fdatetime_seq    = fdatetime_seq[(4*ctraj_len):(4*ctraj_len+ntime)]
        ndaytime         = int(ntime/4)
        fdate_seq        = fdate_seq[ctraj_len:ctraj_len+ndaytime]
        fdateasdate      = fdateasdate[ctraj_len:ctraj_len+ndaytime]

        nupttime         = ntime + 4*ctraj_len
        # NOTE: this is not coded 'universally' at this point... CAUTION
        fuptdatetime_seq = fuptdatetime_seq[4*ctraj_len_orig:4*ctraj_len_orig+nupttime]
        ndayupttime      = ctraj_len + 1
        fuptdate_seq     = fuptdate_seq[ctraj_len_orig:ctraj_len_orig+ndayupttime]

    ## -- WRITE NETCDF OUTPUT (empty, to be filled)
    if fwrite_netcdf:
        writeemptync4D(ofile,fdate_seq,fuptdate_seq,glat,glon,strargs,precision)

    # traj max len
    tml = nupttime - ntime

    ## prepare parcel log to handle trajectories properly 
    if fmemento: # NOTE: must fill array with negative number whose abs exceeds max traj len  
        pIDlogH = -999*np.ones(shape=2000001).astype(int) 


    ###--- pre-loop to produce independent monthly output
    ## NOTE: this is irrelevant for E2P, but crucial for Had (& Ead)
    ## NOTE: we only need to know if some parcel makes it to the ABL, that's it!
    if fmemento and mode == "oper": # skip if multi-counting somehow desired and/or if testing

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
                if mask[alat_ind,alon_ind]==maskval:
                   continue

                ## read ONLY parcel and ABL heights
                ztra, hpbl = readheights(ary[:4,i,:])

                ## p5) LOG ONLY parcels arriving in PBL (or nocturnal layer)
                if ( ztra[0] < np.max(hpbl[:4]) ):
                    ID = int(ary[0,i,0])
                    pIDlogH[ID] = pix - tml # NOTE: tml != npretime (double-check?)
    
    
    ###--- MAIN LOOP
    ## prepare uptake indices
    upt_idx = np.asarray([floor(x) for x in np.arange(0,nupttime)/4])

    ## prepare STATS
    # number of parcels
    tneval = tnjumps = tnnevala = tnnevalm = tnevalp = tnnevalp = tnevalh = tnnevalh = 0
    # precip. statistics
    psum = patt = punatt = pmiss = 0

    ## loop over time to read in files
    if verbose:
        print("\n=== \t Start main program...\n")
    for ix in range(ntime):
        if verbose:
                print("--------------------------------------------------------------------------------------")
        print("Processing "+str(fdatetime_seq[ix]))

        ## 1) read in all files associated with data --> ary is of dimension (ntrajlen x nparcels x nvars)
        ary = readpom( idate    = datetime_seq[ix], 
                       ipath    = ipath+"/"+str(ryyyy), 
                       ifile_base = ifile_base)

        nparcel   = ary.shape[1]
        ntrajleng   = ary.shape[0]
        if verbose:
            print(" TOTAL: " + str(datetime_seq[ix]) + " has " + str(nparcel) + " parcels")

        if mode == "test":
            ntot    = range(1000)
        else:
            ntot    = range(nparcel)

        # figure out where to store data (on which arriving day)
        arv_idx = np.where(np.asarray(fdateasdate)==(fdatetime_seq[ix]-relativedelta(hours=3)).date())[0][0]

        # pre-allocate arrays (repeat at every 4th step)
        if ix%4==0:
            ary_heat     = np.zeros(shape=(ndayupttime,glat.size,glon.size))
            ary_etop     = np.zeros(shape=(ndayupttime,glat.size,glon.size))
            # upscaling measures (currently has to be per day as well)
            if fupscale:
                ipatt = ipmiss = 0

        # STATS: number of parcels per file
        neval = njumps = nnevala = nnevalm = nevalp = nnevalp = nevalh = nnevalh = 0

        ## 2) diagnose P, E, H and npart per grid cell
        for i in ntot:
           
            ## CHECK FOR JUMPS; disregard entire trajectory if it contains a jump
            if fjumps:
                if fjumpsfull:
                    # checking for the full trajectory length
                    jumps = np.array([])
                    for it in range(ntrajleng-1):
                        jumps = np.append(jumps, dist_on_sphere(ary[it,i,2],ary[it,i,1],ary[it+1,i,2],ary[it+1,i,1]))#lat1,lon1,lat2,lon2
                else:
                    jumps = dist_on_sphere(ary[0,i,2],ary[0,i,1],ary[1,i,2],ary[1,i,1])
                if np.any(jumps > cjumps):
                    njumps += int(1)
                    findjump = np.argwhere(jumps > cjumps)
                    if np.any(findjump > 0):
                        print(" !!! ATTENTION: YOU JUST ENCOUNTERED A JUMP IN THE TAIL OF THE TRAJECTORY !!!")
                    continue

            ## - 2.0) only evaluate if the parcel is in target region
	        ## NOTE: I took only the last two time steps for now; should this be 4?
            ## NOTE2: I am assuming that the mask grid is identical to the target grid for now
            mlat_ind, mlon_ind = midpindex(ary[:2,i,:],glon=mlon,glat=mlat)
            alat_ind, alon_ind = arrpindex(ary[0,i,:],glon=mlon,glat=mlat)
            if not mask[mlat_ind,mlon_ind]==maskval and not mask[alat_ind,alon_ind]==maskval:
                continue
            else:

                ## - 2.1) check how far back trajectory should be evaluated
                # NOTE: this could be moved elsewhere...
                # for Hadv (not needed for E2P):
                if fmemento:
                    ID = int(ary[0,i,0])   
                    istepH = pIDlogH[ID]
                    ihf_H = min((ix-istepH+1), tml + 2) 
                else:
                    ihf_H = tml + 2
                
                ## - 2.2) read only the most basic parcel information
                # NOTE: this could easily be done more efficiently
                ztra, hpbl, temp, qv, dens, pres = glanceparcel(ary[:4,i,:])
                
                # sorry, yet another date for writing the P date to the csv (preliminary).
                # because i wanted to have the hours in there as wel (not assign to day only)
                pdate   = str((fdatetime_seq[ix]-relativedelta(hours=3)).strftime('%Y%m%d%H'))

                ## - 2.3) diagnose fluxes

                ##  - 2.3)-KAS: Keune and Schumacher
                if tdiagnosis == 'KAS':

                    ## (a) E2P, evaporation resulting in precipitation
                    if not mask[mlat_ind,mlon_ind]==maskval:
                        nnevalm += 1
                    else:
                        if ( (qv[0]-qv[1]) < cprec_dqv and 
                             ( (q2rh(qv[0],pres[0],temp[0]) + q2rh(qv[1],pres[1],temp[1]))/2 ) > cprec_rh ):

                            # log some statistics
                            nevalp  += 1
                            psum    += abs(qv[0]-qv[1])

                            # read full parcel information
                            lons, lats, temp, ztra, qv, hpbl, dens, pres, pottemp, epottemp = readparcel(ary[:tml+2,i,:])
                            
                            # calculate all required changes along trajectory
                            dq          = trajparceldiff(qv[:], 'diff')
                            dTH         = trajparceldiff(pottemp[:], 'diff')
                            dTHe        = trajparceldiff(epottemp[:], 'diff')

                            # check if traj falls dry & adjust ihf_E if so
                            ihf_E = tml + 2
                            if fdry:
                                ihf_dry = np.where(qv[1:ihf_E]<= 0.00005)[0] + 1 # omit current time step
                                if ihf_dry.size>0:
                                    ihf_E = np.min(ihf_dry)
                                    
                            # identify evaporative moisture uptakes
                            in_PBL     = PBL_check(z=ztra[:ihf_E], h=hpbl[:ihf_E], seth=cevap_hgt, tdiagnosis=tdiagnosis)                      
                            evap_uptk  = (dTHe[:ihf_E-1] - dTH[:ihf_E-1]) > cheat_dtemp 
                            evap_plaus = np.abs(dTH[:ihf_E-1]) < cevap_cc * (dq[:ihf_E-1]) * dTdqs(p_hPa=pres[1:ihf_E]/1e2, q_kgkg=qv[1:ihf_E])
                            evap_idx   = np.where(np.logical_and(in_PBL, np.logical_and(evap_uptk, evap_plaus)))[0]
                            
                            if evap_idx.size==0:
                                # log some statistics
                                nnevalp     += 1
                                pmiss       += abs(qv[0]-qv[1])
                                # log for upscaling
                                if fupscale:
                                    ipmiss      += abs(qv[0]-qv[1])
                                if fwritestats:
                                    pattdata    = [pdate,str(0),str(0),str(abs(qv[0]-qv[1]))]
                                    append2csv(pattfile,pattdata)
                            
                            if evap_idx.size>0:
                                dq_disc     = np.zeros(shape=qv[:ihf_E].size)
                                dq_disc[1:] = linear_discounter(v=qv[1:ihf_E], min_gain=0)
                                if fwritestats:
                                    pattdata    = [pdate,str(np.sum(dq_disc[evap_idx])/qv[1]),str(1-dq_disc[-1]/qv[1]),str(abs(qv[0]-qv[1]))]
                                    append2csv(pattfile,pattdata)
                                ## trajectory-based upscaling
                                prec    = abs(qv[0]-qv[1])
                                fw_orig = dq_disc/qv[1]
                                if explainp=="full":
                                    # upscaling to 100% of trajectory
                                    cfac        = qv[1]/np.sum(dq_disc[evap_idx])
                                    etop        = prec*fw_orig*cfac
                                elif explainp=="max":
                                    # upscaling to (100-IC)% of trajectory
                                    cfac        = (qv[1]-dq_disc[-1])/np.sum(dq_disc[evap_idx])
                                    etop        = prec*fw_orig*cfac
                                elif explainp=="none":
                                    etop        = prec*fw_orig
                                # log for timestep-based upscaling
                                if fupscale:
                                    ipatt       += np.sum(etop[evap_idx])
                                
                                # log some statistics
                                patt    += np.sum(etop[evap_idx])
                                punatt  += prec-np.sum(etop[evap_idx])

                            for itj in evap_idx: 
                                ary_etop[upt_idx[ix+tml-itj],:,:] += gridder(plon=lons[itj:itj+2], plat=lats[itj:itj+2], pval=etop[itj], glon=glon, glat=glat)

                    ## (b) H, surface sensible heat arriving in PBL (or nocturnal layer)
                    if not mask[alat_ind,alon_ind]==maskval:
                        nnevala += 1
                    else:
                        if ( ihf_H >= 2 and 
                             ztra[0] < np.max(hpbl[:4]) ):
                            
                            # log some statistics
                            nevalh  += 1

                            # read full parcel information #NOTE: redundant when parcel has also (somehow) precipitated
                            lons, lats, temp, ztra, qv, hpbl, dens, pres, pottemp, epottemp = readparcel(ary[:ihf_H,i,:])
                             
                            # calculate all required changes along trajectory
                            dq          = trajparceldiff(qv[:], 'diff')
                            dTH         = trajparceldiff(pottemp[:], 'diff')
                            
                            # identify sensible heat uptakes (NOTE: ihf_H is technically not needed below)
                            in_PBL     = PBL_check(z=ztra[:ihf_H], h=hpbl[:ihf_H], seth=cheat_hgt, tdiagnosis=tdiagnosis)
                            heat_uptk  = dTH[:ihf_H-1] > cheat_dtemp
                            heat_plaus = np.abs(dq[:ihf_H-1]) < cheat_cc * (dTH[:ihf_H-1]) * dqsdT(p_hPa=pres[1:ihf_H]/1e2, T_degC=temp[1:ihf_H]-TREF)
                            heat_idx   = np.where(np.logical_and(in_PBL, np.logical_and(heat_uptk, heat_plaus)))[0]

                            # discount uptakes linearly
                            if heat_idx.size==0:
                                # log some statistics
                                nnevalh += 1
                            if heat_idx.size>0:
                                dTH_disc = linear_discounter(v=pottemp[:ihf_H], min_gain=0)

                            # loop through sensible heat uptakes
                            for itj in heat_idx:
                                #NOTE: hardcoded for writing daily data 
                                ary_heat[upt_idx[ix+tml-itj],:,:] += gridder(plon=lons[itj:itj+2], plat=lats[itj:itj+2], pval=dTH_disc[itj], glon=glon, glat=glat)/4 

                            # update parcel log
                            if fmemento:
                                pIDlogH[ID] = ix # NOTE: double-check

                ##  - 2.3)-SOD: Sodemann et al., 2008
                elif tdiagnosis == 'SOD':
         
                    ## (a) E2P
                    if not mask[mlat_ind,mlon_ind]==maskval:
                        nnevalm += 1
                    else:
                        if ( (qv[0]-qv[1]) < 0 and 
                             ( (q2rh(qv[0],pres[0],temp[0]) + q2rh(qv[1],pres[1],temp[1]))/2 ) > 80 ):

                            # log some statistics
                            nevalp  += 1
                            psum    += abs(qv[0]-qv[1])
                            
                            # read full parcel information
                            lons, lats, temp, ztra, qv, hpbl, dens, pres, pottemp, epottemp = readparcel(ary[:tml+2,i,:])

                            # calculate all required changes along trajectory
                            dq          = trajparceldiff(qv[:], 'diff')
                            # check if traj falls dry & adjust ihf_E if so
                            ihf_E = tml + 2
                            if fdry:
                                ihf_dry = np.where(qv[1:ihf_E]<= 0.00005)[0] + 1 # omit current time step
                                if ihf_dry.size>0:
                                    ihf_E = np.min(ihf_dry)

                            # identify evaporative moisture uptakes
                            in_PBL    = trajparceldiff(ztra[:ihf_E], 'mean') < trajparceldiff(hpbl[:ihf_E], 'mean') 
                            evap_uptk = dq[:ihf_E-1] > 0.0002
                            evap_idx  = np.where(np.logical_and(in_PBL, evap_uptk))[0] 

                            if evap_idx.size==0:
                                # log some statistics
                                nnevalp    += 1
                                pmiss      += abs(qv[0]-qv[1])
                                # log for upscaling
                                if fupscale:
                                    ipmiss      += abs(qv[0]-qv[1])
                                if fwritestats:
                                    pattdata    = [pdate,str(0),str(0),str(abs(qv[0]-qv[1]))]
                                    append2csv(pattfile,pattdata)
                            
                            # discount uptakes linearly, scale with precipitation fraction
                            if evap_idx.size>0:
                                dq_disc     = np.zeros(shape=qv[:ihf_E].size)
                                dq_disc[1:] = linear_discounter(v=qv[1:ihf_E], min_gain=0)
                                if fwritestats:
                                    pattdata    = [pdate,str(np.sum(dq_disc[evap_idx])/qv[1]),str(1-dq_disc[-1]/qv[1]),str(abs(qv[0]-qv[1]))]
                                    append2csv(pattfile,pattdata)
                                ## trajectory-based upscaling
                                prec    = abs(qv[0]-qv[1])
                                fw_orig = dq_disc/qv[1]
                                if explainp=="full":
                                    # upscaling to 100% of trajectory
                                    cfac        = qv[1]/np.sum(dq_disc[evap_idx])
                                    etop        = prec*fw_orig*cfac
                                elif explainp=="max":
                                    # upscaling to (100-IC)% of trajectory
                                    cfac        = (qv[1]-dq_disc[-1])/np.sum(dq_disc[evap_idx])
                                    etop        = prec*fw_orig*cfac
                                elif explainp=="none":
                                    etop        = prec*fw_orig
                                # log for timestep-based upscaling
                                if fupscale:
                                    ipatt       += np.sum(etop[evap_idx])
                                
                                # log some statistics
                                patt    += np.sum(etop[evap_idx])
                                punatt  += prec-np.sum(etop[evap_idx])

                            # loop through evaporative uptakes
                            for itj in evap_idx:
                                ary_etop[upt_idx[ix+tml-itj],:,:] += gridder(plon=lons[itj:itj+2], plat=lats[itj:itj+2], pval=etop[itj], glon=glon, glat=glat)
  
    
                    ## (b) H, surface sensible heat (not used originally; analogous to evaporation)
                    if not mask[alat_ind,alon_ind]==maskval:
                        nnevala += 1
                    else:
                        if ( ihf_H >= 2 and
                             ztra[0] < np.max(hpbl[:4]) ):
                            
                            # log some statistics
                            nevalh  += 1

                            # read full parcel information #NOTE: redundant when parcel has also (somehow) precipitated
                            lons, lats, temp, ztra, qv, hpbl, dens, pres, pottemp, epottemp = readparcel(ary[:ihf_H,i,:])

                            # calculate all required changes along trajectory
                            dTH         = trajparceldiff(pottemp[:], 'diff')

                            # identify sensible heat uptakes #NOTE: same as for KAS, ihf_H not needed here (again)
                            in_PBL    = trajparceldiff(ztra[:ihf_H], 'mean') < trajparceldiff(hpbl[:ihf_H], 'mean') 
                            heat_uptk = dTH[:ihf_H-1] > cheat_dtemp
                            heat_idx  = np.where(np.logical_and(in_PBL, heat_uptk))[0]     

                            # discount uptakes linearly
                            if heat_idx.size==0:
                                # log some statistics
                                nnevalh += 1
                            if heat_idx.size>0:
                                dTH_disc = linear_discounter(v=pottemp[:ihf_H], min_gain=0)

                            # loop through sensible heat uptakes
                            for itj in heat_idx:
                                #NOTE: hardcoded for writing daily data 
                                ary_heat[upt_idx[ix+tml-itj],:,:] += gridder(plon=lons[itj:itj+2], plat=lats[itj:itj+2], pval=dTH_disc[itj], glon=glon, glat=glat)/4 

                            # update parcel log
                            if fmemento:
                                pIDlogH[ID] = ix # NOTE: double-check


                ##  SOD2: SODEMANN, 2020; FREMME & SODEMANN, 2019
                elif tdiagnosis == 'SOD2':
         
                    ## (a) E2P
                    if not mask[mlat_ind,mlon_ind]==maskval:
                        nnevalm += 1
                    else:
                        if ( (qv[0]-qv[1]) < 0 and 
                             ( (q2rh(qv[0],pres[0],temp[0]) + q2rh(qv[1],pres[1],temp[1]))/2 ) > 80 ):
                            
                            # log some statistics
                            nevalp  += 1
                            psum    += abs(qv[0]-qv[1])

                            # read full parcel information
                            lons, lats, temp, ztra, qv, hpbl, dens, pres, pottemp, epottemp = readparcel(ary[:tml+2,i,:])

                            # calculate all required changes along trajectory
                            dq          = trajparceldiff(qv[:], 'diff')
                            
                            # check if traj falls dry & adjust ihf_E if so
                            ihf_E = tml + 2
                            if fdry:
                                ihf_dry = np.where(qv[1:ihf_E]<= 0.00005)[0] + 1 # omit current time step
                                if ihf_dry.size>0:
                                    ihf_E = np.min(ihf_dry)

                            # identify evaporative moisture uptakes
                            evap_uptk = dq[:ihf_E-1] > 0.0001
                            evap_idx  = np.where(evap_uptk)[0] 

                            if evap_idx.size==0:
                                # log some statistics
                                nnevalp    += 1
                                pmiss      += abs(qv[0]-qv[1])
                                # log for upscaling
                                if fupscale:
                                    ipmiss      += abs(qv[0]-qv[1])
                                if fwritestats:
                                    pattdata    = [pdate,str(0),str(0),str(abs(qv[0]-qv[1]))]
                                    append2csv(pattfile,pattdata)
                            
                            # discount uptakes linearly, scale with precipitation fraction
                            if evap_idx.size>0:
                                dq_disc     = np.zeros(shape=qv[:ihf_E].size)
                                dq_disc[1:] = linear_discounter(v=qv[1:ihf_E], min_gain=0)
                                if fwritestats:
                                    pattdata    = [pdate,str(np.sum(dq_disc[evap_idx])/qv[1]),str(1-dq_disc[-1]/qv[1]),str(abs(qv[0]-qv[1]))]
                                    append2csv(pattfile,pattdata)
                                ## trajectory-based upscaling
                                prec    = abs(qv[0]-qv[1])
                                fw_orig = dq_disc/qv[1]
                                if explainp=="full":
                                    # upscaling to 100% of trajectory
                                    cfac        = qv[1]/np.sum(dq_disc[evap_idx])
                                    etop        = prec*fw_orig*cfac
                                elif explainp=="max":
                                    # upscaling to (100-IC)% of trajectory
                                    cfac        = (qv[1]-dq_disc[-1])/np.sum(dq_disc[evap_idx])
                                    etop        = prec*fw_orig*cfac
                                elif explainp=="none":
                                    etop        = prec*fw_orig
                                # log for timestep-based upscaling
                                if fupscale:
                                    ipatt       += np.sum(etop[evap_idx])
                                
                                # log some statistics
                                patt    += np.sum(etop[evap_idx])
                                punatt  += prec-np.sum(etop[evap_idx])

                            # loop through evaporative uptakes
                            for itj in evap_idx:
                                ary_etop[upt_idx[ix+tml-itj],:,:] += gridder(plon=lons[itj:itj+2], plat=lats[itj:itj+2], pval=etop[itj], glon=glon, glat=glat)
  
    
                    ## (b) H, surface sensible heat (not used originally; analogous to evaporation)
                    if not mask[alat_ind,alon_ind]==maskval:
                        nnevala += 1
                    else:
                        if ( ihf_H >= 2 and
                             ztra[0] < np.max(hpbl[:4]) ):
                            
                            # log some statistics
                            nevalh  += 1
                            
                            # read full parcel information #NOTE: redundant when parcel has also (somehow) precipitated
                            lons, lats, temp, ztra, qv, hpbl, dens, pres, pottemp, epottemp = readparcel(ary[:ihf_H,i,:])

                            # calculate all required changes along trajectory
                            dTH         = trajparceldiff(pottemp[:], 'diff')

                            # identify sensible heat uptakes #NOTE: same as for KAS, ihf_H not needed here (again)
                            heat_uptk = dTH[:ihf_H-1] > cheat_dtemp
                            heat_idx  = np.where(heat_uptk)[0]     

                            # discount uptakes linearly
                            if heat_idx.size==0:
                                # log some statistics
                                nnevalh += 1
                            if heat_idx.size>0:
                                dTH_disc = linear_discounter(v=pottemp[:ihf_H], min_gain=0)

                            # loop through sensible heat uptakes
                            for itj in heat_idx:
                                #NOTE: hardcoded for writing daily data 
                                ary_heat[upt_idx[ix+tml-itj],:,:] += gridder(plon=lons[itj:itj+2], plat=lats[itj:itj+2], pval=dTH_disc[itj], glon=glon, glat=glat)/4 

                            # update parcel log
                            if fmemento:
                                pIDlogH[ID] = ix # NOTE: double-check


                ##  - 2.3)-SAJ: Stohl and James, 2004
                elif tdiagnosis == 'SAJ':

                    if not mask[mlat_ind,mlon_ind]==maskval:
                        nnevalm += 1
                        continue
                    
                    ## (a) E-P based on ALL parcels residing over target region, no precipitation-criterion used
                    # read full parcel information (but only what is needed; ihf_H)
                    lons, lats, temp, ztra, qv, hpbl, dens, pres, pottemp, epottemp = readparcel(ary[:ihf_H,i,:])

                    # calculate all required changes along trajectory
                    dq          = trajparceldiff(qv[:], 'diff')

                    # add any moisture change along trajectory to respective column sum
                    for itj in range(ihf_H-1):
                        ary_etop[upt_idx[ix+tml-itj],:,:] += gridder(plon=lons[itj:itj+2], plat=lats[itj:itj+2], pval=dq[itj], glon=glon, glat=glat)

                    # update parcel log
                    if fmemento:
                        pIDlogH[ID] = ix # NOTE: making use of heat parcel log for E-P


        neval   = len(ntot)
        if verbose:
            if fjumps:
                print(" STATS: Encountered " + str(njumps) + " ({:.2f}".format(100*njumps/neval) +"%) jumps.")
            print(" STATS: Evaluated "+str(neval-nnevala)+" ({:.2f}".format(100*(neval-nnevala)/(neval)) +"%) arriving parcels inside mask (advection).")
            if nnevalh!=0:
                print(" --- ATTENTION: "+str(nnevalh)+"/"+str(neval-nnevala)+" arriving parcels are not associated with any heat uptakes...")
            print(" STATS: Evaluated "+str(neval-nnevalm)+" ({:.2f}".format(100*(neval-nnevalm)/(neval)) +"%) midpoint parcels inside mask (precipitation).")
            print(" STATS: Evaluated "+str(nevalp)+" ({:.2f}".format(100*(nevalp)/(neval-nnevalm)) +"%) precipitating parcels.")
            if nnevalp!=0:
                print(" --- ATTENTION: "+str(nnevalp)+"/"+str(nevalp)+" precipitating parcels are not associated with any evap uptakes...")
        
        ## SOME DAILY CALCULATIONS
        if ( (ix+1)%4==0 ):
            # DAILY UPSCALING of E2P, taking into account the missing trajectories (i.e. the ones without any uptakes)
            if fupscale and nnevalp!=0:
                if ipatt==0:
                    warnings.warn(" --- WARNING: there were no trajectories with uptakes, so upscaling is impossible...")
                else:
                    upsfac              = 1+(ipmiss/ipatt)
                    ary_etop[:,:,:]     = upsfac*ary_etop[:,:,:]
                    # corrections for final statistics
                    patt                += -np.sum(ipatt) + np.sum(ary_etop)
                    pmiss               += -np.sum(ipmiss)
                    if verbose:
                        print(" * Upscaling... (factor: {:.4f}".format(upsfac)+")")
            # Convert units
            if verbose:
                print(" * Converting units...")
            ary_etop[:,:,:] = convertunits(ary_etop[:,:,:], garea, "E")
            ary_heat[:,:,:] = convertunits(ary_heat[:,:,:], garea, "H")
    
            if fwrite_netcdf:
                writenc4D(ofile,arv_idx,ary_etop,ary_heat)
        
        ## STATS summary
        tneval  += neval
        tnjumps += njumps
        tnnevala+= nnevala
        tnnevalm+= nnevalm
        tnevalp += nevalp
        tnnevalp+= nnevalp
        tnevalh += nevalh
        tnnevalh+= nnevalh

    if ftimethis:
        megatoc = timeit.default_timer()
        if verbose:
            print("\n=== \t End main program (total runtime so far: ",str(round(megatoc-megatic, 2)),"seconds) \n")
    
    if verbose:
        if fwrite_netcdf:
            print("\n Successfully written: "+ofile+" !\n")
    
    with open(statfile,'w') as sfile:
        writer=csv.writer(sfile, delimiter='\t', lineterminator='\n',quoting = csv.QUOTE_NONE, quotechar='',)
        writer.writerow(["* - PARCEL STATISTICS: "])
        writer.writerow(["   --- TOTAL EVALUATED PARCELS:       " , str(tneval)])
        writer.writerow(["   --- # PARCELS FILTERED OUT (JUMPS):" , str(tnjumps)])
        writer.writerow([" "])
        writer.writerow(["   --- # PARCELS ARRIVING INSIDE MASK:" , str(tneval-tnnevala)])
        writer.writerow(["   --- # PARCELS EVAL. FOR HEAT-ADV:  " , str(tnevalh)+" ({:.2f}".format(100*tnevalh/(tneval-tnnevala))+"%)"])
        if tnevalh!=0:
            writer.writerow(["   ----- WITHOUT UPTAKES IN THE TRAJ: " , str(tnnevalh)+" ({:.2f}".format(100*tnnevalh/(tnevalh))+"%)"])
        writer.writerow([" "])
        writer.writerow(["   --- # PARCELS MIDPOINT INSIDE MASK:" , str(tneval-tnnevalm)])
        writer.writerow(["   --- # PARCELS EVAL. FOR PRECIP:    " , str(tnevalp)+" ({:.2f}".format(100*tnevalp/(tneval-tnnevalm))+"%)"])
        if tnevalp!=0:
            writer.writerow(["   ----- WITHOUT UPTAKES IN THE TRAJ: " , str(tnnevalp)+" ({:.2f}".format(100*tnnevalp/(tnevalp))+"%)"])
        writer.writerow([" "])
        if psum!=0:
            writer.writerow([" * - PRECIPITATION STATISTICS: "])
            writer.writerow(["   --- ATTRIBUTED FRACTION:             {:.2f}".format(patt/psum)])
            writer.writerow(["   --- UNATTRIBUTED FRACTION (TRAJEC):  {:.2f}".format(punatt/psum)])
            writer.writerow(["   --- UNATTRIBUTED FRACTION (NO-UPT):  {:.2f}".format(pmiss/psum)])
    if verbose: 
        with open(statfile, 'r') as sfile:
            print(sfile.read())
