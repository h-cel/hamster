#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MAIN FUNCTIONS FOR 02_attribution
"""

############################################################################
#############################    SETTINGS ##################################

def main_attribution(
           ryyyy, ayyyy, am, ad,
           ipath, ifile_base, ifile_format, ipath_f2t,
           opath, ofile_base,
           mode,
           gres,
           maskfile,
           maskval,
           verbose,
           veryverbose,
           tdiagnosis,
           ctraj_len,
           cheat_dtemp, # used for E,H,P (if cprec_dqv==None)
           cheat_cc, cevap_cc, # for H, E diagnosis (lower = more strict)
           cevap_hgt, cheat_hgt, # set min ABLh, disabled if 0 | NOTE: to be unified
           cprec_dqv, cprec_dtemp, cprec_rh,
           cpbl_strict,
           fjumps,
           fjumpsfull,
           cjumps,
           refdate,
           fwrite_netcdf,
           precision,
           ftimethis,
           fdry,
           fmemento,
           mattribution,
           crandomnit,
           randatt_forcall,
           explainp,fdupscale,fmupscale,
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
    if fwritestats and mattribution=="linear":
        pfilename = str(ofile_base)+"_attr_r"+str(ryyyy)[-2:]+"_"+str(ayyyy)+"-"+str(am).zfill(2)+"_p-linear-attribution.csv"
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
        print("\n============================================================================================================\n")
    ## Resets & consistency checks
    if fwritestats and mattribution=="random":
        print(" ! Option <writestats> not yet available for attribution method random. Continuing anyhow...")
    if mode=="oper" and precision=="f4":
        precision = "f8"
        print(" ! Single precision should only be used for testing. Reset to double-precision.")
    if verbose:
        print(" ! using input path: \t", 	ipath)
        print(" ! using internal timer: \t" +str(ftimethis) )
        print(" ! using mode: \t" +str(mode))
        print(" ! using attribution method (P): \t" +str(mattribution))
        print(" ! using trajectory-based upscaling: \t" +str(explainp))
        if mattribution == "random":
            print(" ! using minimum iterations: \t" +str(crandomnit))
        print(" ! using daily upscaling: \t" +str(fdupscale))
        print(" ! using monthly upscaling: \t" +str(fmupscale))
        if fwrite_netcdf:
            print(" ! writing netcdf output: \t" +str(fwrite_netcdf) )
            print(" \t ! with grid resolution: \t", str(gres) )
            print(" \t ! output file: \t", opath+"/"+ofilename)
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
        print("\n--- WARNING: the grid from the maskfile is not identical to the target grid... please check. Proceeding nevertheless. \n")

    ## -- DATES
    # NOTE: we begin at 06 UTC...
    datetime_bgn    = datetime.datetime.strptime(str(ayyyy)+"-"+str(am).zfill(2)+"-"+str(ad).zfill(2)+"-06", "%Y-%m-%d-%H") 
    # get end date (always 00 UTC of the 1st of the next month)
    nayyyy          = (datetime_bgn + relativedelta(months=1)).strftime('%Y')
    nam             = (datetime_bgn + relativedelta(months=1)).strftime('%m')
    datetime_end    = datetime.datetime.strptime(str(nayyyy)+"-"+str(nam).zfill(2)+"-01-00",  "%Y-%m-%d-%H")
    timestep        = datetime.timedelta(hours=6)
    datetime_seq    = []
    fdatetime_seq   = []
    idatetime       = datetime_bgn

    # create arrival datetime string & datetime object
    while idatetime <= datetime_end:
        datetime_seq.append(idatetime.strftime('%Y%m%d%H'))
        fdatetime_seq.append(idatetime)
        idatetime += timestep

    # figure out when first uptakes occur (for pre-loop)
    uptdatetime_bgn = datetime_bgn - datetime.timedelta(days=ctraj_len) - datetime.timedelta(hours=3)

    # aggregate to daily, NOTE: arrival at 00 UTC means parcel has arrived on prev day    
    fdate_seq = np.unique([fdt.date() for fdt in fdatetime_seq[:-1]]).tolist() # omit last dt object (00:00)
    # keep a copy of datetime.date formatted list for arv_idx further below
    fdateasdate = np.copy(fdate_seq).tolist() # NOTE: using deepcopy instead of np.copy would be more proper
    # convert these datetime.date objects to datetime.datetime objects for netCDF writing
    for idt in range(len(fdate_seq)):
        fdate_seq[idt]    = datetime.datetime(fdate_seq[idt].year, fdate_seq[idt].month, fdate_seq[idt].day)
    # NOTE: better to keep these as lists to maintain consistency

    # calculate number of time steps, also aggregated to daily resolution
    ntime           = len(fdatetime_seq)
    ndaytime        = len(fdate_seq)

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

    ## -- WRITE NETCDF OUTPUT (empty, to be filled)
    if fwrite_netcdf:
        writeemptync4D(ofile,fdate_seq,np.arange(-ctraj_len,1),glat,glon,strargs,precision)

    # traj max len, expressed in input data (6-hourly) steps
    tml = int(4*ctraj_len) # hardcoded for 6-hourly input
    # compact form of max traj len in days (used for array filling w/ shortened uptake dim)
    ctl = ctraj_len

    ### MEMENTO --- pre-loop to produce independent monthly output
    ## NOTE: this is irrelevant for E2P, but crucial for Had (& Ead)
    ## NOTE: we only need to know if some parcel makes it to the ABL, that's it!
    ## NOTE: must fill array with negative number whose abs exceeds max traj len  
    if fmemento: 
        pidlog = -999*np.ones(shape=2100000).astype(int) 
        
        if mode == "oper": # skip if multi-counting somehow desired and/or if testing
#            pidlog = preloop(datetime_bgn, uptdatetime_bgn, timestep,
#                         ipath, ifile_base, ifile_format,
#                         ryyyy,
#                         mask, mlat, mlon, maskval,
#                         pidlog, tml,
#                         verbose)
            ###--- PRELOOP v2
            # NOTE: code further below this function call here could also be moved to
            #       first main loop iteration, so that no dim checking necessary
            ntrajstep, _, _ = readtrajdim(idate    = datetime_seq[0],
                                          ipath    = ipath+"/"+str(ryyyy),
                                          ifile_base = ifile_base,
                                          ifile_format = ifile_format,
                                          verbose=False)

            if ntrajstep < tml+2+4:
                # only do this if data really isn't already 'there'
                extendarchive = grabmesomehpbl(pposbasepath=ipath_f2t, ryyyy=ryyyy,
                                               fdatetime_beg=fdatetime_seq[0], tml=tml,
                                               verbose=verbose)
            else:
                if verbose: print("\n=== \t INFO: no pre-loop needed, trajectories are long enough")

    ###--- MAIN LOOP

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

        # safety exit hardcoded for FLEXPART–ERA-Interim runs
        if ryyyy==ayyyy and am==12 and ix==range(ntime)[-1]:
            continue

        ## 1) read in all files associated with data --> ary is of dimension (ntrajlen x nparcels x nvars)
        ary = readtraj(idate    = datetime_seq[ix], 
                       ipath    = ipath+"/"+str(ryyyy), 
                       ifile_base = ifile_base, 
                       ifile_format = ifile_format,
                       verbose=verbose)

        nparcel     = ary.shape[1]
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
            ary_heat     = np.zeros(shape=(ctl+1,glat.size,glon.size))
            ary_etop     = np.zeros(shape=(ctl+1,glat.size,glon.size))
            # upscaling measures (currently has to be per day as well)
            if fdupscale:
                ipatt = ipmiss = 0

        # STATS: number of parcels per file
        neval = njumps = nnevala = nnevalm = nevalp = nnevalp = nevalh = nnevalh = 0

        # grab extended trajectory data
        if fmemento:
            # NOTE: pom data can come with duplicate IDs; remove to avoid (some) trouble
            if not np.unique(ary[0,:,0]).size == ary[0,:,0].size:
                print("\t INFO: duplicates detected, original pom array shape=",ary.shape)
                _, ikeep = np.unique(ary[0,:,0], return_index=True)
                ary = ary[:,ikeep,:]
                print("\t INFO: duplicates eliminated, new pom array shape=",ary.shape)
                nparcel   = ary.shape[1]   # update
                ntot      = range(nparcel)
            # NOTE: yet another ******* pom fix ... to be removed!
            #       (flex2traj output is not affected by this)
            IDs   = ary[0,:,0]
            thresidx=int((9/10)*ary.shape[1]) # this should do the trick
            # simply shift to indices > 2e6
            IDs[thresidx:][IDs[thresidx:]<3000] += 2e6

            # extract what is needed from extendarchive if trajs 'too short'
            if ary.shape[0] < tml+2+4 and ix < ctraj_len*4:
                extendtrajs = np.empty(shape=(4, nparcel, 2))
                for pp in range(4):
                    allIDs     = extendarchive[-(4-pp+ix)][:,0]
                    extendtrajs[pp,:,0] = extendarchive[-(4-pp+ix)][:,0][np.where(np.isin(allIDs, ary[0,:,0]))] # ID
                    extendtrajs[pp,:,1] = extendarchive[-(4-pp+ix)][:,1][np.where(np.isin(allIDs, ary[0,:,0]))] # hpbl

        ## 2) diagnose P, E, H and npart per grid cell
        for i in ntot:
           
            ## CHECK FOR JUMPS; disregard entire trajectory if it contains a jump
            if fjumps and ifile_format=='dat.gz' and (ary[0,i,0] < 3000 or ary[0,i,0] > 2e6):
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
            if not mask[alat_ind,alon_ind]==maskval:
                nnevalm     += 1
                nnevala     += 1
                continue
            else:

                ## - 2.1) check how far back trajectory should be evaluated
                ID = int(ary[0,i,0])
                if fmemento and ary[0,i,3] < np.max(ary[:4,i,7]):
                    if ix < ctraj_len*4: # rely on (extended) traj data

                        # check if parcel has been inside before
                        is_inmask = whereinmask(mask=mask, maskval=maskval, masklat=mlat, masklon=mlon,
                                                trajlat=ary[:(tml+2),i,2], trajlon=ary[:(tml+2),i,1])

                        # check if parcel was 'in PBL'
                        hgt = ary[:(tml+2),i,3] # consistent with max traj len
                        if ary.shape[0] < tml+2+4:
                            longhpbl = np.concatenate((ary[:(tml+2),i,7], extendtrajs[:,i,1]))
                        else:
                            longhpbl = ary[:(tml+2+4),i,7]
                        is_inpbl = np.where(hgt < maxlastn(longhpbl, n=4)[:-4])[0] # omit last 4

                        # check where parcel inside and 'in PBL'
                        is_arrv  = np.intersect1d(is_inmask, is_inpbl)

                        # now determine ihf_H
                        if is_arrv.size > 1:
                            ihf_H = is_arrv[is_arrv>0].min() + 1
                        elif is_arrv == 0:
                            ihf_H = tml + 2 # use max traj len
                        else:
                            raise RuntimeError('--- FATAL ERROR: Schrödingers cat situation; is parcel inside and in PBL, or not?')

#                        ## checking mode
#                        ihf_H_orig = min((ix-pidlog[ID]+1), tml + 2)
#                        if not ihf_H == ihf_H_orig:
#                            print("DISCREPANCIES DETECTED!!!!!!! if ID < 3'000 and using pom data, this is to be expected..")
#                            print("ihf_H, ihf_H_orig=", ihf_H, ihf_H_orig)
#                            print("ID=",ID)
#                            _ = input("proceed?")

                    else: # fully rely on log from now
                        istep = pidlog[ID]
                        ihf_H = min((ix-istep+1), tml + 2)
                else:
                    ihf_H = tml + 2
                
                ## - 2.2) read only the most basic parcel information
                # NOTE: this could easily be done more efficiently
                hgt, hpbl, temp, qv, dens, pres = glanceparcel(ary[:4,i,:])
                
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

                            # prec
                            prec    = abs(qv[0]-qv[1])
                            # log some statistics
                            nevalp  += 1
                            psum    += prec

                            # read full parcel information
                            lons, lats, temp, hgt, qv, hpbl, dens, pres, pottemp, epottemp = readparcel(ary[:tml+2,i,:])
                            
                            # calculate all required changes along trajectory
                            dq          = trajparceldiff(qv[:], 'diff')
                            dTH         = trajparceldiff(pottemp[:], 'diff')
                            dTHe        = trajparceldiff(epottemp[:], 'diff')

                            # evaluate only until trajectory falls dry
                            ihf_E = tml + 2
                            if fdry and np.any(qv[1:ihf_E]<= 0.00005):
                                ihf_E = np.min(np.where(qv[1:ihf_E]<= 0.00005)[0] + 1)
                                    
                            # identify evaporative moisture uptakes
                            evap_idx = uptake_locator_KAS(cevap_hgt, cpbl_strict, hgt[:ihf_E], hpbl[:ihf_E],
                                                          (dTHe-dTH)[:ihf_E-1], cheat_dtemp, dTH[:ihf_E-1], dq[:ihf_E-1], cevap_cc,
                                                          dTdqs(p_hPa=pres[1:ihf_E]/1e2, q_kgkg=qv[1:ihf_E]))
                           
                            # WRITE STATS
                            if fwritestats and mattribution=="linear":
                                if evap_idx.size==0:
                                    pattdata    = [pdate,str(0),str(0),str(prec)]
                                elif evap_idx.size>0:
                                    etop    = linear_attribution_p(qv[:ihf_E],iupt=evap_idx,explainp="none")
                                    pattdata= [pdate,str(np.sum(etop[evap_idx]/prec)),str(1-etop[-1]/prec),str(prec)]
                                append2csv(pattfile,pattdata)
                           
                            # TRAJECTORY WITHOUT UPTAKES
                            if evap_idx.size==0:
                                # log some statistics (for upscaling)
                                nnevalp     += 1
                                pmiss       += prec
                                if fdupscale:
                                    ipmiss      += prec
                            
                            # ATTRIBUTION
                            if evap_idx.size>0:
                                if mattribution=="linear":
                                    etop    = linear_attribution_p(qv[:ihf_E],iupt=evap_idx,explainp=explainp)
                                elif mattribution=="random":
                                    etop    = random_attribution_p(qtot=qv[:ihf_E],iupt=evap_idx,explainp=explainp,
                                            nmin=crandomnit,forc_all=randatt_forcall,
                                            verbose=verbose,veryverbose=veryverbose)
                                for itj in evap_idx:
                                    ary_etop[ctl-(itj+3-ix%4)//4,:,:] += gridder(plon=lons[itj:itj+2], plat=lats[itj:itj+2], pval=etop[itj], glon=glon, glat=glat)
                                # log some statistics (for upscaling)
                                patt    += np.sum(etop[evap_idx])
                                punatt  += prec-np.sum(etop[evap_idx])
                                if fdupscale:
                                    ipatt       += np.sum(etop[evap_idx])
                                    ipmiss      += prec-np.sum(etop[evap_idx])

                    ## (b) H, surface sensible heat arriving in PBL (or nocturnal layer)
                    if not mask[alat_ind,alon_ind]==maskval:
                        nnevala += 1
                    else:
                        if ( ihf_H >= 2 and 
                             hgt[0] < np.max(hpbl[:4]) ):
                            
                            # log some statistics
                            nevalh  += 1

                            # read full parcel information #NOTE: redundant when parcel has also (somehow) precipitated
                            lons, lats, temp, hgt, qv, hpbl, dens, pres, pottemp, epottemp = readparcel(ary[:ihf_H,i,:])
                             
                            # calculate all required changes along trajectory
                            dq          = trajparceldiff(qv[:], 'diff')
                            dTH         = trajparceldiff(pottemp[:], 'diff')
                            
                            # identify sensible heat uptakes (NOTE: ihf_H is technically not needed below)
                            heat_idx = uptake_locator_KAS(cheat_hgt, cpbl_strict, hgt[:ihf_H], hpbl[:ihf_H],
                                                          dTH[:ihf_H-1], cheat_dtemp, dq[:ihf_H-1], dTH[:ihf_H-1], cheat_cc,
                                                          dqsdT(p_hPa=pres[1:ihf_H]/1e2, T_degC=temp[1:ihf_H]-TREF))

                            # discount uptakes linearly
                            if heat_idx.size==0:
                                # log some statistics
                                nnevalh += 1
                            if heat_idx.size>0:
                                dTH_disc = linear_discounter(v=pottemp[:ihf_H], min_gain=0)

                            # loop through sensible heat uptakes
                            for itj in heat_idx:
                                #NOTE: hardcoded for writing daily data 
                                ary_heat[ctl-(itj+3-ix%4)//4,:,:] += gridder(plon=lons[itj:itj+2], plat=lats[itj:itj+2], pval=dTH_disc[itj], glon=glon, glat=glat)/4 

                            # update parcel log
                            if fmemento:
                                pidlog[ID] = ix # NOTE: double-check

                ##  - 2.3)-SOD: Sodemann et al., 2008
                elif tdiagnosis == 'SOD':
         
                    ## (a) E2P
                    if not mask[mlat_ind,mlon_ind]==maskval:
                        nnevalm += 1
                    else:
                        if ( (qv[0]-qv[1]) < 0 and 
                             ( (q2rh(qv[0],pres[0],temp[0]) + q2rh(qv[1],pres[1],temp[1]))/2 ) > 80 ):

                            # prec
                            prec    = abs(qv[0]-qv[1])
                            # log some statistics
                            nevalp  += 1
                            psum    += prec
                            
                            # read full parcel information
                            lons, lats, temp, hgt, qv, hpbl, dens, pres, pottemp, epottemp = readparcel(ary[:tml+2,i,:])

                            # calculate all required changes along trajectory
                            dq          = trajparceldiff(qv[:], 'diff')
                            # evaluate only until trajectory falls dry
                            ihf_E = tml + 2
                            if fdry and np.any(qv[1:ihf_E]<= 0.00005):
                                ihf_E = np.min(np.where(qv[1:ihf_E]<= 0.00005)[0] + 1)

                            # identify uptake locations
                            is_inpbl    = trajparceldiff(hgt[:ihf_E], 'mean') < 1.5*trajparceldiff(hpbl[:ihf_E], 'mean')
                            is_uptk     = dq[:ihf_E-1] > 0.0002
                            evap_idx    = np.where(np.logical_and(is_inpbl, is_uptk))[0] 
                               
                            # WRITE STATS
                            if fwritestats and mattribution=="linear":
                                if evap_idx.size==0:
                                    pattdata    = [pdate,str(0),str(0),str(prec)]
                                elif evap_idx.size>0:
                                    etop    = linear_attribution_p(qv[:ihf_E],iupt=evap_idx,explainp="none")
                                    pattdata= [pdate,str(np.sum(etop[evap_idx]/prec)),str(1-etop[-1]/prec),str(prec)]
                                append2csv(pattfile,pattdata)

                            ## TRAJECTORY WITHOUT UPTAKE LOCATIONS
                            if evap_idx.size==0:
                                # log some statistics (for upscaling)
                                nnevalp    += 1
                                pmiss      += prec
                                if fdupscale:
                                    ipmiss      += prec
                            
                            # ATTRIBUTION
                            if evap_idx.size>0:
                                if mattribution=="linear":
                                    etop    = linear_attribution_p(qv[:ihf_E],iupt=evap_idx,explainp=explainp)
                                elif mattribution=="random":
                                    etop    = random_attribution_p(qtot=qv[:ihf_E],iupt=evap_idx,explainp=explainp,
                                            nmin=crandomnit, forc_all=randatt_forcall,
                                            verbose=verbose,veryverbose=veryverbose)
                                    print("Attributed fraction: "+str(np.sum(etop[evap_idx])/prec))
                                for itj in evap_idx:
                                    ary_etop[ctl-(itj+3-ix%4)//4,:,:] += gridder(plon=lons[itj:itj+2], plat=lats[itj:itj+2], pval=etop[itj], glon=glon, glat=glat)
                                # log some statistics (for upscaling)
                                patt    += np.sum(etop[evap_idx])
                                punatt  += prec-np.sum(etop[evap_idx])
                                if fdupscale:
                                    ipatt       += np.sum(etop[evap_idx])
                                    ipmiss      += prec-np.sum(etop[evap_idx]) 
    
                    ## (b) H, surface sensible heat (not used originally; analogous to evaporation)
                    if not mask[alat_ind,alon_ind]==maskval:
                        nnevala += 1
                    else:
                        if ( ihf_H >= 2 and
                             hgt[0] < np.max(hpbl[:4]) ):
                            
                            # log some statistics
                            nevalh  += 1

                            # read full parcel information #NOTE: redundant when parcel has also (somehow) precipitated
                            lons, lats, temp, hgt, qv, hpbl, dens, pres, pottemp, epottemp = readparcel(ary[:ihf_H,i,:])

                            # calculate all required changes along trajectory
                            dTH         = trajparceldiff(pottemp[:], 'diff')

                            # identify sensible heat uptakes #NOTE: same as for KAS, ihf_H not needed here (again)
                            is_inpbl    = trajparceldiff(hgt[:ihf_H], 'mean') < 1.5*trajparceldiff(hpbl[:ihf_H], 'mean')
                            is_uptk     = dTH[:ihf_H-1] > cheat_dtemp
                            heat_idx    = np.where(np.logical_and(is_inpbl, is_uptk))[0]     

                            # discount uptakes linearly
                            if heat_idx.size==0:
                                # log some statistics
                                nnevalh += 1
                            if heat_idx.size>0:
                                dTH_disc = linear_discounter(v=pottemp[:ihf_H], min_gain=0)

                            # loop through sensible heat uptakes
                            for itj in heat_idx:
                                #NOTE: hardcoded for writing daily data 
                                ary_heat[ctl-(itj+3-ix%4)//4,:,:] += gridder(plon=lons[itj:itj+2], plat=lats[itj:itj+2], pval=dTH_disc[itj], glon=glon, glat=glat)/4 

                            # update parcel log
                            if fmemento:
                                pidlog[ID] = ix # NOTE: double-check


                ##  SOD2: SODEMANN, 2020; FREMME & SODEMANN, 2019
                elif tdiagnosis == 'SOD2':
         
                    ## (a) E2P
                    if not mask[mlat_ind,mlon_ind]==maskval:
                        nnevalm += 1
                    else:
                        if ( (qv[0]-qv[1]) < 0 and 
                             ( (q2rh(qv[0],pres[0],temp[0]) + q2rh(qv[1],pres[1],temp[1]))/2 ) > 80 ):

                            # prec
                            prec    = abs(qv[0]-qv[1])
                            # log some statistics
                            nevalp  += 1
                            psum    += prec

                            # read full parcel information
                            lons, lats, temp, hgt, qv, hpbl, dens, pres, pottemp, epottemp = readparcel(ary[:tml+2,i,:])

                            # calculate all required changes along trajectory
                            dq          = trajparceldiff(qv[:], 'diff')
                            
                            # evaluate only until trajectory falls dry
                            ihf_E = tml + 2
                            if fdry and np.any(qv[1:ihf_E]<= 0.00005):
                                ihf_E   = np.min(np.where(qv[1:ihf_E]<= 0.00005)[0] + 1)

                            # identify evaporative moisture uptakes
                            is_uptk     = dq[:ihf_E-1] > 0.0001
                            evap_idx    = np.where(is_uptk)[0] 

                            # WRITE STATS
                            if fwritestats and mattribution=="linear":
                                if evap_idx.size==0:
                                    pattdata    = [pdate,str(0),str(0),str(prec)]
                                elif evap_idx.size>0:
                                    etop    = linear_attribution_p(qv[:ihf_E],iupt=evap_idx,explainp="none")
                                    pattdata= [pdate,str(np.sum(etop[evap_idx]/prec)),str(1-etop[-1]/prec),str(prec)]
                                append2csv(pattfile,pattdata)
 
                            # TRAJECTORY WITHOUT UPTAKE LOCATIONS
                            if evap_idx.size==0:
                                # log some statistics (for upscaling)
                                nnevalp    += 1
                                pmiss      += prec
                                if fdupscale:
                                    ipmiss      += prec
                            
                            # ATTRIBUTION
                            if evap_idx.size>0:
                                if mattribution=="linear":
                                    etop    = linear_attribution_p(qv[:ihf_E],iupt=evap_idx,explainp=explainp)
                                elif mattribution=="random":
                                    etop    = random_attribution_p(qtot=qv[:ihf_E],iupt=evap_idx,explainp=explainp,
                                            nmin=crandomnit,forc_all=randatt_forcall,
                                            verbose=verbose,veryverbose=veryverbose)
                                for itj in evap_idx:
                                    ary_etop[ctl-(itj+3-ix%4)//4,:,:] += gridder(plon=lons[itj:itj+2], plat=lats[itj:itj+2], pval=etop[itj], glon=glon, glat=glat)
                                # log some statistics (for upscaling)
                                patt    += np.sum(etop[evap_idx])
                                punatt  += prec-np.sum(etop[evap_idx])
                                if fdupscale:
                                    ipatt       += np.sum(etop[evap_idx])
                                    ipmiss      += prec-np.sum(etop[evap_idx]) 

                    ## (b) H, surface sensible heat (not used originally; analogous to evaporation)
                    if not mask[alat_ind,alon_ind]==maskval:
                        nnevala += 1
                    else:
                        if ( ihf_H >= 2 and
                             hgt[0] < np.max(hpbl[:4]) ):
                            
                            # log some statistics
                            nevalh  += 1
                            
                            # read full parcel information #NOTE: redundant when parcel has also (somehow) precipitated
                            lons, lats, temp, hgt, qv, hpbl, dens, pres, pottemp, epottemp = readparcel(ary[:ihf_H,i,:])

                            # calculate all required changes along trajectory
                            dTH         = trajparceldiff(pottemp[:], 'diff')

                            # identify sensible heat uptakes #NOTE: same as for KAS, ihf_H not needed here (again)
                            is_uptk     = dTH[:ihf_H-1] > cheat_dtemp
                            heat_idx    = np.where(is_uptk)[0]     

                            # discount uptakes linearly
                            if heat_idx.size==0:
                                # log some statistics
                                nnevalh += 1
                            if heat_idx.size>0:
                                dTH_disc = linear_discounter(v=pottemp[:ihf_H], min_gain=0)

                            # loop through sensible heat uptakes
                            for itj in heat_idx:
                                #NOTE: hardcoded for writing daily data 
                                ary_heat[ctl-(itj+3-ix%4)//4,:,:] += gridder(plon=lons[itj:itj+2], plat=lats[itj:itj+2], pval=dTH_disc[itj], glon=glon, glat=glat)/4 

                            # update parcel log
                            if fmemento:
                                pidlog[ID] = ix # NOTE: double-check

        if fjumps:
            neval   = len(ntot)-njumps
        else: 
            neval   = len(ntot)
        if verbose:
            if fjumps:
                print(" STATS: Encountered " + str(njumps) + " ({:.2f}".format(100*njumps/neval) +"%) jumps.")
            print(" STATS: Evaluated "+str(neval-nnevala)+" ({:.2f}".format(100*(neval-nnevala)/(neval)) +"%) arriving parcels inside mask (advection).")
            if nnevalh!=0:
                print(" --- ATTENTION: "+str(nnevalh)+"/"+str(neval-nnevala)+" arriving parcels are not associated with any heat uptakes...")
            print(" STATS: Evaluated "+str(neval-nnevalm)+" ({:.2f}".format(100*(neval-nnevalm)/(neval)) +"%) midpoint parcels inside mask (precipitation).")
            if nnevalm==neval:
                print(" STATS: Evaluated "+str(nevalp)+" precipitating parcels.")
            else:
                print(" STATS: Evaluated "+str(nevalp)+" ({:.2f}".format(100*(nevalp)/(neval-nnevalm)) +"%) precipitating parcels.")
            if nnevalp!=0:
                print(" --- ATTENTION: "+str(nnevalp)+"/"+str(nevalp)+" precipitating parcels are not associated with any evap uptakes...")
        
        ## SOME DAILY CALCULATIONS
        if ( (ix+1)%4==0 ):
            # DAILY UPSCALING of E2P, taking into account the missing trajectories (i.e. the ones without any uptakes)
            if fdupscale and (nnevalp!=0 or ipmiss!=0):
                if ipatt==0:
                    print(" \n--- WARNING: there were no trajectories with uptakes, so upscaling is impossible...\n ")
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
    
    # MONTHLY UPSCALING of E2P, taking into account the missing trajectories (i.e. the ones without any uptakes)
    if fmupscale and (pmiss!=0 or punatt!=0):
        if patt==0:
            print(" \n--- WARNING: there were no trajectories with uptakes, so upscaling is impossible...\n")
        else:
            upsfac              = 1+((pmiss+punatt)/patt)
            # load full etop array and upscale
            fdata       = nc4.Dataset(ofile,'r+')
            uns_etop    = fdata.variables['E2P'][:]
            ups_etop    = upsfac*uns_etop
            fdata.variables['E2P'][:]= ups_etop
            fdata.close()
            # corrections for final statistics
            patt                += np.sum(pmiss) +np.sum(punatt)
            pmiss               += -np.sum(pmiss)
            punatt              += -np.sum(punatt)
            if verbose:
                print(" * Monthly upscaling for unattributed precipitation... (factor: {:.4f}".format(upsfac)+")")

    if ftimethis:
        megatoc = timeit.default_timer()
        if verbose:
            print("\n=== \t End main program (total runtime so far: ",str(round(megatoc-megatic, 2)),"seconds) \n")
    
    if verbose:
        if fwrite_netcdf:
            print("\n Successfully written: "+ofile+" !\n")

    writestats_02(statfile,tneval,tnjumps,tnnevala,tnevalh,tnnevalh,tnnevalm,tnevalp,tnnevalp,patt,psum,punatt,pmiss) 
    if verbose: 
        with open(statfile, 'r') as sfile:
            print(sfile.read())
