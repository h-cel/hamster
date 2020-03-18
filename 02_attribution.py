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
           refdate,
           fwrite_netcdf,ftimethis,
           fdry,fmemento,fcc_advanced,fvariable_mass,
           strargs):

    # TODO: add missing features
    if fcc_advanced or fvariable_mass:
        raise SystemExit("---- ABORTED: no can do, not implemented!")
 
    ## construct precise input and storage paths
    mainpath  = ipath+str(ryyyy)+"/"
    ofilename = str(ofile_base)+"_attr_r"+str(ryyyy)[-2:]+"_"+str(ayyyy)+"-"+str(am).zfill(2)+".nc"
    ofile     = opath+"/"+ofilename
    
    ## read netcdf mask
    with nc4.Dataset(maskfile) as f:
        mask = f['mask'][:]
        mlat = f['lat'][:]
        mlon = f['lon'][:]

    ########### LOG W/IN PYTHON SCRIPT by redirecting output #############
    
    if verbose:
        disclaimer()
        print("\n PROCESSING: \t", 	ayyyy, "-", str(am).zfill(2))
        print("\n============================================================================================================")
        print(" ! using input path: \t", 	ipath)
        print(" ! using variable mass: \t" +str(fvariable_mass) )
        if fvariable_mass:
            print(" \t ! reference date for number of particles: \t" +str(refdate) )
        if fwrite_netcdf:
            print(" ! writing netcdf output: \t" +str(fwrite_netcdf) )
            print(" \t ! with grid resolution:: \t", str(gres) )
            print(" \t ! output file: \t", opath+"/"+ofilename)
        print(" ! using internal timer: \t" +str(ftimethis) )
        print(" ! using mode: \t" +str(mode))
        print("\n============================================================================================================")
        print("\n============================================================================================================")

    ## start timer
    if ftimethis:
        megatic = timeit.default_timer()
    
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

    # TESTMODE
    if mode == "test":
        ntime            = 4 #NOTE: use multiples of 4 only, else output is not saved 
        datetime_seq     = datetime_seq[:ntime]
        fdatetime_seq    = fdatetime_seq[:ntime]
        ndaytime         = 1
        fdate_seq        = fdate_seq[:ndaytime]
        
        nupttime         = 4*ctraj_len + 1
        fuptdatetime_seq = fuptdatetime_seq[:nupttime]
        ndayupttime      = (ctraj_len) + 1
        fuptdate_seq     = fuptdate_seq[:ndayupttime]

    ## -- WRITE NETCDF OUTPUT (empty, to be filled)
    if fwrite_netcdf:
        writeemptync4D(ofile,fdate_seq,fuptdate_seq,glat,glon,strargs)

    # traj max len
    tml = nupttime - ntime

    # set some default thresholds
    cprec_dqv    = default_thresholds(cprec_dqv) 
    # read in reference distribution of parcels
    if fvariable_mass:
        ary_rnpart   = get_refnpart(refdate=refdate, ryyyy=ryyyy, glon=glon, glat=glat)

    ## prepare parcel log to handle trajectories properly 
    if fmemento: # NOTE: must fill array with negative number whose abs exceeds max traj len  
        pIDlogH = -999*np.ones(shape=2000001).astype(int) 

    ## prepare uptake indices
    upt_idx = np.asarray([floor(x) for x in np.arange(0,nupttime)/4])

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
        arv_idx = np.where(np.asarray(fdateasdate)==(fdatetime_seq[ix]-relativedelta(hours=3)).date())[0][0]

        # pre-allocate arrays (repeat at every 4th step)
        if ix%4==0:
            ary_heat     = np.zeros(shape=(ndayupttime,glat.size,glon.size))
            ary_etop     = np.zeros(shape=(ndayupttime,glat.size,glon.size))

        ## 2) diagnose P, E, H and npart per grid cell
        for i in ntot:
            
            ## - 2.0) only evaluate if the parcel is in target region
            ## NOTE: I took only the last two time steps for now; should this be 4?
            ## NOTE2: I am assuming that the mask grid is identical to the target grid for now
            lat_ind, lon_ind = midpindex(ary[:2,i,:],glon=mlon,glat=mlat)
            if mask[lat_ind,lon_ind]!=maskval:
                pass
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

                ## - 2.3) diagnose fluxes

                ##  - 2.3)-KAS: Keune and Schumacher
                if tdiagnosis == 'KAS':

                    ## (a) E2P, evaporation resulting in precipitation
                    if ( (qv[0]-qv[1]) < cprec_dqv and 
                         q2rh(qv[0], pres[0], temp[0]) > cprec_rh  and
                         q2rh(qv[1], pres[1], temp[1]) > cprec_rh ):

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
                        
                        if evap_idx.size>0:
                            dq_disc     = np.zeros(shape=qv[:ihf_E].size-1)
                            dq_disc[1:] = linear_discounter(v=qv[1:ihf_E], min_gain=0, min_loss=0)
                            etop        = ((qv[0]-qv[1])/qv[1])*dq_disc

                        for itj in evap_idx: 
                            ary_etop[upt_idx[ix+tml-itj],:,:] += gridder(plon=lons[itj:itj+2], plat=lats[itj:itj+2], pval=etop[itj], glon=glon, glat=glat)

                    ## (b) H, surface sensible heat arriving in PBL (or nocturnal layer)
                    if ( ztra[0] < np.max(hpbl[:4]) ):

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
                        if heat_idx.size>0:
                            dTH_disc = linear_discounter(v=pottemp[:ihf_H], min_gain=0, min_loss=0)

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
                    if ( (qv[0]-qv[1]) < 0 and 
                         q2rh((qv[0]+qv[1])/2, (pres[0]+pres[1])/2, (temp[0]+temp[1])/2) > 80 ):

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

                        # discount uptakes linearly, scale with precipitation fraction
                        if evap_idx.size>0:
                            dq_disc     = np.zeros(shape=qv[:ihf_E].size-1)
                            dq_disc[1:] = linear_discounter(v=qv[1:ihf_E], min_gain=0, min_loss=0)
                            etop        = ((qv[0]-qv[1])/qv[1])*dq_disc 

                        # loop through evaporative uptakes
                        for itj in evap_idx:
                            ary_etop[upt_idx[ix+tml-itj],:,:] += gridder(plon=lons[itj:itj+2], plat=lats[itj:itj+2], pval=etop[itj], glon=glon, glat=glat)
  
    
                    ## (b) H, surface sensible heat (not used originally; analogous to evaporation)
                    if ( ztra[0] < np.max(hpbl[:4]) ):

                        # read full parcel information #NOTE: redundant when parcel has also (somehow) precipitated
                        lons, lats, temp, ztra, qv, hpbl, dens, pres, pottemp, epottemp = readparcel(ary[:ihf_H,i,:])

                        # calculate all required changes along trajectory
                        dTH         = trajparceldiff(pottemp[:], 'diff')

                        # identify sensible heat uptakes #NOTE: same as for KAS, ihf_H not needed here (again)
                        in_PBL    = trajparceldiff(ztra[:ihf_H], 'mean') < trajparceldiff(hpbl[:ihf_H], 'mean') 
                        heat_uptk = dTH[:ihf_H-1] > cheat_dtemp
                        heat_idx  = np.where(np.logical_and(in_PBL, heat_uptk))[0]     

                        # discount uptakes linearly
                        if heat_idx.size>0:
                            dTH_disc = linear_discounter(v=pottemp[:ihf_H], min_gain=0, min_loss=0)

                        # loop through sensible heat uptakes
                        for itj in heat_idx:
                            #NOTE: hardcoded for writing daily data 
                            ary_heat[upt_idx[ix+tml-itj],:,:] += gridder(plon=lons[itj:itj+2], plat=lats[itj:itj+2], pval=dTH_disc[itj], glon=glon, glat=glat)/4 

                        # update parcel log
                        if fmemento:
                            pIDlogH[ID] = ix # NOTE: double-check


                ##  - 2.3)-SAJ: Stohl and James, 2004
                elif tdiagnosis == 'SAJ':

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


        # Convert units, but only after the last time step of each day
        if ( (ix+1)%4==0 ):
            if verbose:
                print(" * Converting units...")
            ary_etop[:,:,:] = convertunits(ary_etop[:,:,:], garea, "E")
            ary_heat[:,:,:] = convertunits(ary_heat[:,:,:], garea, "H")
    
            if fwrite_netcdf:
                writenc4D(ofile,arv_idx,ary_etop,ary_heat)
        
#TODO: decide what to do with parcel mass scaling here
#    # Scale with parcel mass
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
    
    if verbose:
        if fwrite_netcdf:
            print("\n Successfully written: "+ofile+" !")
