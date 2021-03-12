#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MAIN FUNCTIONS FOR 01_diagnosis
"""


############################################################################
#############################    SETTINGS ##################################

def main_diagnosis(
           ryyyy, ayyyy, am, ad,
           ipath, ifile_base,
           opath, ofile_base,
           mode,
           gres,
           verbose,
           veryverbose,
           fproc_npart,
           # E criteria
           cevap_dqv, fevap_drh, cevap_drh, cevap_hgt,
           # P criteria
           cprec_dqv, cprec_rh,
           # H criteria
           cheat_dtemp, fheat_drh, cheat_drh, cheat_hgt, fheat_rdq, cheat_rdq,
           # pbl and height criteria
           cpbl_method, cpbl_strict, cpbl_factor,
           refdate,
           fwrite_netcdf,
           precision,
           ftimethis,fvariable_mass,
           strargs):

    ## Perform consistency checks
    if mode=="oper" and precision=="f4":
        precision = "f8"
        print("Single precision should only be used for testing. Reset to double-precision.")
    if fvariable_mass and not fproc_npart:
        fproc_npart = True
        print("Have to process all parcels for variable mass...")

    ## Construct precise input and storage paths
    mainpath  = ipath+str(ryyyy)+"/"
    ofilename = str(ofile_base)+"_diag_r"+str(ryyyy)[-2:]+"_"+str(ayyyy)+"-"+str(am).zfill(2)+".nc"
    ofile     = opath+"/"+ofilename
    
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
            print(" \t ! with grid resolution: \t", str(gres) )
            print(" \t ! output file: \t", opath+"/"+ofilename)
        print(" ! using internal timer: \t" +str(ftimethis) )
        print(" ! using mode: \t" +str(mode))
        print("\n============================================================================================================")
        print("\n============================================================================================================")

    ## Start timer
    if ftimethis:
        megatic = timeit.default_timer()

    ## Prepare grid
    glon, glat, garea = makegrid(resolution=gres)

    ## Handle dates
    date_bgn        = datetime.datetime.strptime(str(ayyyy)+"-"+str(am).zfill(2)+"-"+str(ad).zfill(2), "%Y-%m-%d")
    # get end date (always 00 UTC of the 1st of the next month)
    nayyyy          = (date_bgn + relativedelta(months=1)).strftime('%Y')
    nam             = (date_bgn + relativedelta(months=1)).strftime('%m')
    date_end        = datetime.datetime.strptime(str(nayyyy)+"-"+str(nam).zfill(2)+"-01-00",  "%Y-%m-%d-%H")
    timestep        = datetime.timedelta(hours=6)
    date_seq        = []
    fdate_seq       = []
    mfdate_seq      = []
    idate           = date_bgn + timestep
    while idate <= date_end:
        date_seq.append(idate.strftime('%Y%m%d%H'))
        fdate_seq.append(idate)
        mfdate_seq.append(idate-timestep/2) # -dt/2 for backward run
        idate   += timestep
    ntime           = len(date_seq)

    ##-- TESTMODE
    if mode == "test":
        ntime       = 12
        date_seq    = date_seq[0:ntime]
        fdate_seq   = fdate_seq[0:ntime]
        mfdate_seq  = mfdate_seq[0:ntime]

    ## Create empty netcdf file (to be filled)
    if fwrite_netcdf:
        writeemptync(ofile,mfdate_seq,glon,glat,strargs,precision)

    # Read in reference distribution of parcels
    if fvariable_mass:
        print(" \n !!! WARNING !!! With this version, variable mass can only be applied to 01_diagnosis -- it cannot be used consistently for all steps yet! \n")
        ary_rnpart   = get_refnpart(refdate=refdate, ryyyy=ryyyy, glon=glon, glat=glat)
    
    ## Pre-allocate arrays
    ary_heat     = np.zeros(shape=(glat.size,glon.size))
    ary_evap     = np.zeros(shape=(glat.size,glon.size))
    ary_prec     = np.zeros(shape=(glat.size,glon.size))
    ary_npart    = np.zeros(shape=(glat.size,glon.size))
    ary_pnpart   = np.zeros(shape=(glat.size,glon.size))
    ary_enpart   = np.zeros(shape=(glat.size,glon.size))
    ary_hnpart   = np.zeros(shape=(glat.size,glon.size))

    ##-- LOOP THROUGH FILES
    if verbose:
        print("\n=== \t Start main program: 01_diagnosis...\n")

    for ix in range(ntime):
        
        if verbose:
            print("--------------------------------------------------------------------------------------")
            print("Processing "+str(fdate_seq[ix]))

        ## Read date related trajectories -> ary is of dimension (ntrajlen x nparticles x nvars)
        ary         = readtraj(idate        = date_seq[ix], 
                               ipath        = ipath+"/"+str(ryyyy), 
                               ifile_base   = ifile_base,
                               verbose      = verbose)
        ary         = calc_allvars(ary)
        dq          = trajparceldiff(ary[:,:,5], "diff")
        mrh         = np.apply_over_axes(np.mean, ary[:,:,10], 0)
        dTH         = trajparceldiff(ary[:,:,11], "diff")

        nparticle   = ary.shape[1]
        if verbose:
            print(" TOTAL: " + str(date_seq[ix]) + " has " + str(nparticle) + " parcels")

        ## TESTMODE: less parcels
        if mode == "test":
            ntot    = range(10000,10100)
        else:
            ntot    = range(nparticle)

        #smalltic = timeit.default_timer()

        ##-- LOOP OVER PARCELS TO DIAGNOSE P, E, H (and npart) and assign to grid
        if fproc_npart:
            # get midpoint indices on grid from ary
            lary        = [y for y in (np.moveaxis(ary, 1, 0))] # convert to list for first dimension (parcels) to be able to use map
            res         = np.asarray(list(map(lambda p: midpindex(p, glon=glon, glat=glat), lary)))
            lat_ind     = res[:,0]
            lon_ind     = res[:,1]
            for i in ntot:
                ## log number of parcels
                ary_npart[lat_ind[i],lon_ind[i]] += int(1)
        
        ## Precipitation
        # filter for P
        fdqv        = np.where(dq[0,:]<cprec_dqv)
        frh         = np.where(mrh[0,:]>cprec_rh)
        isprec      = np.intersect1d(fdqv,frh)
        p_ary       = ary[:,isprec,:]
        # get midpoint indices on grid from ary
        lary        = [y for y in (np.moveaxis(p_ary, 1, 0))] # convert to list for first dimension (parcels) to be able to use map
        pres        = np.asarray(list(map(lambda p: midpindex(p, glon=glon, glat=glat), lary)))
        plat_ind    = pres[:,0]
        plon_ind    = pres[:,1]
        for i in range(p_ary.shape[1]):
            ary_prec[plat_ind[i],plon_ind[i]]   += dq[:,isprec[i]][0]
            ary_pnpart[plat_ind[i],plon_ind[i]] += int(1)
            
        ## Evaporation
        # filter for E
        fdqv        = np.where(dq[0,:]>cevap_dqv)
        lary        = [y for y in (np.moveaxis(ary[:,:,[3,7]], 1, 0))] # convert to list for first dimension (parcels) to be able to use map
        fpbl        = np.where(np.asarray(list(map(lambda p: pblcheck2(p, cpbl_strict, cevap_hgt, cpbl_factor, cpbl_method), lary)))[:,0])
        lary        = [y for y in (np.moveaxis(ary[:,:,[10]], 1, 0))] # convert to list for first dimension (parcels) to be able to use map
        fdrh        = np.where(np.asarray(list(map(lambda p: drhcheck(p, checkit=fevap_drh, maxdrh=cevap_drh), lary)))[:,0])
        isevap      = reduce(np.intersect1d, (fdqv, fpbl, fdrh))
        e_ary       = ary[:,isevap,:]
        # get midpoint indices on grid from ary
        lary        = [y for y in (np.moveaxis(e_ary, 1, 0))] # convert to list for first dimension (parcels) to be able to use map
        eres        = np.asarray(list(map(lambda p: midpindex(p, glon=glon, glat=glat), lary)))
        elat_ind    = eres[:,0]
        elon_ind    = eres[:,1]
        for i in range(e_ary.shape[1]):
                ary_evap[elat_ind[i],elon_ind[i]]   += dq[:,isevap[i]][0]
                ary_enpart[elat_ind[i],elon_ind[i]] += int(1)

        ## Sensible heat
        # filter for H
        fdTH        = np.where(dTH[0,:]>cheat_dtemp)
        lary        = [y for y in (np.moveaxis(ary[:,:,[3,7]], 1, 0))] # convert to list for first dimension (parcels) to be able to use map
        fpbl        = np.where(np.asarray(list(map(lambda p: pblcheck2(p, cpbl_strict, cheat_hgt, cpbl_factor, cpbl_method), lary)))[:,0])
        lary        = [y for y in (np.moveaxis(ary[:,:,[10]], 1, 0))] # convert to list for first dimension (parcels) to be able to use map
        fdrh        = np.where(np.asarray(list(map(lambda p: drhcheck(p, checkit=fheat_drh, maxdrh=cheat_drh), lary)))[:,0])
        lary        = [y for y in (np.moveaxis(ary[:,:,[5]], 1, 0))] # convert to list for first dimension (parcels) to be able to use map
        frdq        = np.where(np.asarray(list(map(lambda p: rdqvcheck(p, checkit=fheat_rdq, maxrdqv=cheat_rdq), lary)))[:,0])
        isheat      = reduce(np.intersect1d, (fdTH, fpbl, fdrh, frdq))
        h_ary       = ary[:,isheat,:]
        # get midpoint indices on grid from ary
        lary        = [y for y in (np.moveaxis(h_ary, 1, 0))] # convert to list for first dimension (parcels) to be able to use map
        hres        = np.asarray(list(map(lambda p: midpindex(p, glon=glon, glat=glat), lary)))
        hlat_ind    = hres[:,0]
        hlon_ind    = hres[:,1]
        for i in range(h_ary.shape[1]):
                ary_heat[hlat_ind[i],hlon_ind[i]]   += dTH[:,isheat[i]][0]
                ary_hnpart[hlat_ind[i],hlon_ind[i]] += int(1)

        #smalltoc = timeit.default_timer()
        #print("=== \t All parcels: ",str(round(smalltoc-smalltic, 2)),"seconds \n")

        ## Convert units
        if verbose:
            print(" * Converting units...")
        ary_prec[:,:] = convertunits(ary_prec[:,:], garea, "P")
        ary_evap[:,:] = convertunits(ary_evap[:,:], garea, "E")
        ary_heat[:,:] = convertunits(ary_heat[:,:], garea, "H")

        ## Scale with parcel mass
        if fvariable_mass:
            print(" !!! WARNING !!! With this version, variable mass can only be applied to 01_diagnosis -- it cannot be used consistently for all steps yet!")
            if verbose: 
                print(" * Applying variable mass...")
            ary_prec[:,:]         = scale_mass(ary_prec[:,:], ary_npart[:,:], ary_rnpart)
            ary_evap[:,:]         = scale_mass(ary_evap[:,:], ary_npart[:,:], ary_rnpart)
            ary_heat[:,:]         = scale_mass(ary_heat[:,:], ary_npart[:,:], ary_rnpart)

        if fwrite_netcdf:
            writenc(ofile,ix,ary_prec[:,:],ary_evap[:,:],ary_heat[:,:],ary_npart[:,:],ary_pnpart[:,:],ary_enpart[:,:],ary_hnpart[:,:])

        ## re-init. arrays
        ary_npart[:,:]  = 0
        ary_pnpart[:,:] = 0
        ary_enpart[:,:] = 0
        ary_hnpart[:,:] = 0
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
