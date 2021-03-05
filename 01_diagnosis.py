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
           tdiagnosis,
           cheat_dtemp, # used for E,H,P (if cprec_dqv==None)
           cevap_dqv,
           cevap_hgt, cheat_hgt, # set min ABLh, disabled if 0 
           cprec_dqv, cprec_rh,
           cpbl_strict,
           # pbl and height criteria
           cpbl_method,
           cpbl_strict,
           cpbl_factor,
           refdate,
           fwrite_netcdf,
           precision,
           ftimethis,fvariable_mass,
           strargs):

    
    # Consistency checks
    if mode=="oper" and precision=="f4":
        precision = "f8"
        print("Single precision should only be used for testing. Reset to double-precision.")

    ## construct precise input and storage paths
    mainpath  = ipath+str(ryyyy)+"/"
    ofilename = str(ofile_base)+"_diag_r"+str(ryyyy)[-2:]+"_"+str(ayyyy)+"-"+str(am).zfill(2)+".nc"
    ofile     = opath+"/"+ofilename

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
            print(" \t ! with grid resolution: \t", str(gres) )
            print(" \t ! output file: \t", opath+"/"+ofilename)
        print(" ! using internal timer: \t" +str(ftimethis) )
        print(" ! using mode: \t" +str(mode))
        print("\n============================================================================================================")
        print("\n============================================================================================================")

    ## start timer
    if ftimethis:
        megatic = timeit.default_timer()

    ## -- GRID
    glon, glat, garea = makegrid(resolution=gres)

    ## -- DATES
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

    # TESTMODE
    if mode == "test":
        ntime       = 12
        date_seq    = date_seq[0:ntime]
        fdate_seq   = fdate_seq[0:ntime]
        mfdate_seq  = mfdate_seq[0:ntime]

    ## -- WRITE NETCDF OUTPUT (empty, to be filled)
    if fwrite_netcdf:
        writeemptync(ofile,mfdate_seq,glon,glat,strargs,precision)

    # read in reference distribution of parcels
    if fvariable_mass:
        print(" \n !!! WARNING !!! With this version, variable mass can only be applied to 01_diagnosis -- it cannot be used consistently for all steps yet! \n")
        ary_rnpart   = get_refnpart(refdate=refdate, ryyyy=ryyyy, glon=glon, glat=glat)
    
    ## ------- LOOP OVER DATES (FILES) 
    if verbose:
        print("\n=== \t Start main program: 01_diagnosis...\n")

    # pre-allocate arrays
    ary_heat     = np.zeros(shape=(glat.size,glon.size))
    ary_evap     = np.zeros(shape=(glat.size,glon.size))
    ary_prec     = np.zeros(shape=(glat.size,glon.size))
    ary_npart    = np.zeros(shape=(glat.size,glon.size))
    ary_pnpart   = np.zeros(shape=(glat.size,glon.size))
    ary_enpart   = np.zeros(shape=(glat.size,glon.size))
    ary_hnpart   = np.zeros(shape=(glat.size,glon.size))

    for ix in range(ntime):
        
        if verbose:
            print("--------------------------------------------------------------------------------------")
            print("Processing "+str(fdate_seq[ix]))

        ## READ DATE RELATED TRAJECTORIES -> ary is of dimension (ntrajlen x nparticles x nvars)
        ary         = readtraj(idate        = date_seq[ix], 
                               ipath        = ipath+"/"+str(ryyyy), 
                               ifile_base   = ifile_base,
                               verbose      = verbose)
        nparticle   = ary.shape[1]
        if verbose:
            print(" TOTAL: " + str(date_seq[ix]) + " has " + str(nparticle) + " parcels")

        ## TEST mode: less parcels
        if mode == "test":
            ntot    = range(10000,10100)
        else:
            ntot    = range(nparticle)

        #smalltic = timeit.default_timer()
       
        # log variables
       
        ## ------- LOOP OVER PARCELS TO DIAGNOSE P, E, H (and npart) and assign to grid 
        for i in ntot:

            ## get midpoint at the very beginning
            #lat_mid, lon_mid = readmidpoint(ary[:,i,:])
            lat_ind, lon_ind = midpindex(ary[:2,i,:],glon=glon,glat=glat) # :2 for clarity, not needed
            ## log number of parcels
            ary_npart[lat_ind,lon_ind] += int(1)

            ## read only necessary parcel information
            qv, temp, hgt, hpbl     = readsparcel(ary[:2,i,:]) # load only ('last') 2 steps
            dq                      = trajparceldiff(qv, 'diff') 
            pottemp                 = readpottemp(ary[:2,i,:])
            dTH                     = trajparceldiff(pottemp, 'diff')
            pres                    = readpres(ary[:2,i,:])

            ## PRECIPITATION
            if (dq < cprec_dqv and 
                    ( (q2rh(qv[0],pres[0],temp[0]) + q2rh(qv[1],pres[1],temp[1]))/2 ) > cprec_rh ):
                ary_prec[lat_ind,lon_ind] += dq
                ary_pnpart[lat_ind,lon_ind] += int(1)
            
            ## ALLPBL
            if tdiagnosis == 'ALLPBL':

                ## evaporation
                if ( pblcheck(cpbl_strict,hgt,hpbl,cevap_hgt,cpbl_factor,cpbl_method) and dq>cevap_dqv ):
                    ary_evap[lat_ind,lon_ind]  += dq
                    ary_enpart[lat_ind,lon_ind] += int(1)

                ## sensible heat
                if ( pblcheck(cpbl_strict,hgt,hpbl,cheat_hgt,cpbl_factor,cpbl_method) and dTH>cheat_dtemp ):
                    ary_heat[lat_ind,lon_ind]  += dTH
                    ary_hnpart[lat_ind,lon_ind] += int(1)

            ## SOD: SODEMANN ET AL., 2008
            elif tdiagnosis == 'SOD':
        
                ## evaporation
                if (dq > 0.0002 and  
                        (hgt[0]+hgt[1])/2 < 1.5*(hpbl[0]+hpbl[1])/2):
                    ary_evap[lat_ind,lon_ind] += dq
                    ary_enpart[lat_ind,lon_ind] += int(1)
    
                ## sensible heat (not used originally; analogous to evaporation)
                if ((dTH > cheat_dtemp) and 
                        (hgt[0]+hgt[1])/2 < 1.5*(hpbl[0]+hpbl[1])/2):
                    ary_heat[lat_ind,lon_ind] += dTH
                    ary_hnpart[lat_ind,lon_ind] += int(1)

            ## SOD2: SODEMANN, 2020; FREMME & SODEMANN, 2019
            elif tdiagnosis == 'SOD2':
        
                ## evaporation
                if (dq > 0.0001): 
                    ary_evap[lat_ind,lon_ind] += dq
                    ary_enpart[lat_ind,lon_ind] += int(1)
    
                ## sensible heat (not used originally; analogous to evaporation)
                if (dTH > cheat_dtemp):
                    ary_heat[lat_ind,lon_ind] += dTH
                    ary_hnpart[lat_ind,lon_ind] += int(1)
        
        #smalltoc = timeit.default_timer()
        #print("=== \t All parcels: ",str(round(smalltoc-smalltic, 2)),"seconds \n")

        # Convert units
        if verbose:
            print(" * Converting units...")
        ary_prec[:,:] = convertunits(ary_prec[:,:], garea, "P")
        ary_evap[:,:] = convertunits(ary_evap[:,:], garea, "E")
        ary_heat[:,:] = convertunits(ary_heat[:,:], garea, "H")

        # Scale with parcel mass
        if fvariable_mass:
            print(" !!! WARNING !!! With this version, variable mass can only be applied to 01_diagnosis -- it cannot be used consistently for all steps yet!")
            if verbose: 
                print(" * Applying variable mass...")
            ary_prec[:,:]         = scale_mass(ary_prec[:,:], ary_npart[:,:], ary_rnpart)
            ary_evap[:,:]         = scale_mass(ary_evap[:,:], ary_npart[:,:], ary_rnpart)
            ary_heat[:,:]         = scale_mass(ary_heat[:,:], ary_npart[:,:], ary_rnpart)

        if fwrite_netcdf:
            writenc(ofile,ix,ary_prec[:,:],ary_evap[:,:],ary_heat[:,:],ary_npart[:,:],ary_pnpart[:,:],ary_enpart[:,:],ary_hnpart[:,:])

        # re-init. arrays
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

        
