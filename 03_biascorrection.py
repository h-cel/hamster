#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MAIN FUNCTION FOR 03_biascorrection
"""

############################################################################
#############################    SETTINGS ##################################

def main_biascorrection(
           ryyyy, ayyyy, am,
           opathA, # attribution (output)
           opathD, # diagnosis (output)
           ipathR, # reference data (input)
           opath, ofile_base,           # output
           mode,
           maskfile,
           maskval,
           verbose,
           veryverbose,
           fdebug,
           fwrite_netcdf,
           precision,
           strargs):
    
    ## SOME PRELIMINARY SETTINGS TO REDUCE OUTPUT
    ## suppressing warnings, such as
    #  invalid value encountered in true_divide
    #  invalid value encountered in multiply 
    if not fdebug:
        np.seterr(divide='ignore', invalid='ignore')

    ## construct precise input and storage paths
    attrfile  = opathA+"/"+str(ofile_base)+"_attr_r"+str(ryyyy)[-2:]+"_"+str(ayyyy)+"-"+str(am).zfill(2)+".nc"
    ofilename = str(ofile_base)+"_biascor-attr_r"+str(ryyyy)[-2:]+"_"+str(ayyyy)+"-"+str(am).zfill(2)+".nc"
    ofile     = opath+"/"+ofilename

    #### DISCLAIMER
    if verbose:
        disclaimer()
        print("\n PROCESSING: \t", 	ayyyy, "-", str(am).zfill(2)+"\n")
    ## Resets & consistency checks
    if mode=="oper" and precision=="f4":
        precision = "f8"
        print(" ! Single precision should only be used for testing. Reset to double-precision.")
    if verbose:
        print(" ! using input paths: \t")
        print("\t"+str(opathD))
        print("\t"+str(opathA))
        print(" ! using reference data from: \t")
        print("\t"+str(ipathR))
        print(" ! using mode: \t" +str(mode))
        print(" ! writing netcdf output: \t")
        print("\t"+str(ofile))
        print("\n============================================================================================================")
        print("\n")

    ##--1. load attribution data; grab all uptake days ############################    
    if verbose: 
        print(" * Reading attribution data...")

    with nc4.Dataset(attrfile, mode="r") as f:
        E2Psrt       = np.asarray(f['E2P'][:])
        Hadsrt       = np.asarray(f['H'][:])
        arrival_time = nc4.num2date(f['time'][:], f['time'].units, f['time'].calendar)
        utime_srt    = np.asarray(f['level'][:])
        lats         = np.asarray(f['lat'][:])
        lons         = np.asarray(f['lon'][:])
        areas        = 1e6*np.nan_to_num(gridded_area_exact(lats, res=abs(lats[1]-lats[0]), nlon=lons.size))

    ## expand uptake dimension to conform to original approach
    utime_first = arrival_time[0] - timedelta(days=utime_srt.size-1) # utime_srt.size-1 == trajlen (in days)
    uptake_time = np.asarray([utime_first+timedelta(days=nday) for nday in range(utime_srt.size-1+arrival_time.size)])

    date_bgn = datetime.date(uptake_time[0].year, uptake_time[0].month, uptake_time[0].day)
    date_end = datetime.date(uptake_time[-1].year, uptake_time[-1].month, uptake_time[-1].day)

    ## 'reshape' arrays accordingly
    E2P = np.empty(shape=(arrival_time.size,uptake_time.size,lats.size,lons.size))
    Had = np.empty(shape=(arrival_time.size,uptake_time.size,lats.size,lons.size))
    for iat in range(arrival_time.size):
        E2P[iat,iat:iat+utime_srt.size,:,:] = E2Psrt[iat,:,:,:]
        Had[iat,iat:iat+utime_srt.size,:,:] = Hadsrt[iat,:,:,:]
    del(E2Psrt, Hadsrt) # not needed anymore

    
    ##--2. load diagnosis data ####################################################
    if verbose: 
        print(" * Reading diagnosis data...")
    
    # get required months
    ayears  = np.asarray([ut.year  for ut in uptake_time])
    amonths = np.asarray([ut.month for ut in uptake_time])
    uyears = (np.unique(np.column_stack((ayears, amonths)), axis=0))[:,0]

    # read concatenated data
    totlats, totlons    = read01data(opathD,ofile_base,ryyyy,uptake_time,var="grid")
    gridcheck(lats,totlats,lons,totlons)
    ftime               = read01data(opathD,ofile_base,ryyyy,uptake_time,var="time")
    E                   = read01data(opathD,ofile_base,ryyyy,uptake_time,var="E")
    P                   = read01data(opathD,ofile_base,ryyyy,uptake_time,var="P")
    H                   = read01data(opathD,ofile_base,ryyyy,uptake_time,var="H")
    
    ## must check if data comes in daily resolution; fix if not
    dates   = np.asarray([datetime.date(it.year, it.month, it.day) for it in ftime])
    udates = np.unique(dates)

    if udates.size != dates.size:
        Etot    = convert2daily(E,ftime,udates,fagg="sum")
        Ptot    = convert2daily(P,ftime,udates,fagg="sum")
        Htot    = convert2daily(H,ftime,udates,fagg="mean")
    else: ## NOTE: preparing for DAILY OUTPUT from 01_diagnosis
        ## neat, we can just proceed
        Etot = E
        Ptot = P
        Htot = H
        
    ## only keep what is really needed
    datecheck(date_bgn,udates)
    ibgn = np.where(udates==date_bgn)[0][0]
    iend = np.where(udates==date_end)[0][0]
    Etot = Etot[ibgn:iend+1]
    Ptot = Ptot[ibgn:iend+1]
    Htot = Htot[ibgn:iend+1]
    datestot = udates[ibgn:iend+1]  
    
    ## make sure we grabbed the right data
    uptake_dates = np.asarray([datetime.date(iut.year, iut.month, iut.day) for iut in uptake_time])
    if not np.array_equal(uptake_dates, datestot):
        raise SystemExit("---- hold your horses; datetime matching failed!")
    
    ## clean up a bit
    #del(E, P, H, Ex, Px, Hx)
    del(E, P, H)

    ##--3. load reference data ####################################################
    """
    this part is STRICTLY CODED FOR (12-hourly) ERA-INTERIM only (so far),
    and HARDCODED too
    """
    if verbose: 
        print(" * Reading reference data...")
    
    Eref = eraloader_12hourly(var='e',
                     datapath=ipathR+"/evap_12hourly/E_1deg_",
                     maskpos=True,
                     maskneg=False,
                     uptake_years=uyears,
                     uptake_dates=uptake_dates, lats=lats, lons=lons)
        
    Href = eraloader_12hourly(var='sshf',
                     datapath=ipathR+"/sshf_12hourly/H_1deg_",
                     maskpos=True,
                     maskneg=False,
                     uptake_years=uyears,
                     uptake_dates=uptake_dates, lats=lats, lons=lons)
    
    Pref = -eraloader_12hourly(var='tp',
                     datapath=ipathR+"/tp_12hourly/P_1deg_",
                     maskpos=False, # do NOT set this to True!
                     maskneg=True,
                     uptake_years=uyears,
                     uptake_dates=uptake_dates, lats=lats, lons=lons)
    
    ##--4. scale ##################################################################
    if verbose: 
        print(" * Starting bias correction...")
    
    ## P-scaling requires arrival region mask
    with nc4.Dataset(maskfile) as f:
        mask = f['mask'][:]
        mlat = f['lat'][:]
        mlon = f['lon'][:]   
    
    ## area-weight arrival region precipitation (FLEXPART & REF)
    if verbose: print("---- INFO: area-weighting precipitation data...")
    xla, xlo    = np.where(mask==maskval) # P[:,xla,xlo] is merely a 2D array... ;)
    ibgn        = np.where(uptake_time==arrival_time[0])[0][0] # only arrival days!
    # convert from mm to m3
    PrefTS      = np.nansum(areas[xla,xlo]*Pref[ibgn:,xla,xlo], axis=1)/1e3
    PtotTS      = np.nansum(areas[xla,xlo]*Ptot[ibgn:,xla,xlo], axis=1)/1e3

    ## here is where the magic (part I) happens.
    #******************************************************************************
    alpha_Had = np.divide(Had, Htot)
    alpha_E2P = np.divide(E2P, Etot)
    #******************************************************************************
   
    ## NOTE: as of now, there is absolutely no check whatsoever concerning
    ## the fractions; if e.g. only 3 6-hourly values are used to generate
    ## daily diagnosis data, this can result in a division by zero above,
    ## so that scaled data blows up to infinity (this actually happened).
    ## hence, check if any alpha clearly exceeds 1, and warn the user
    ## AGAIN that the output cannot be fully trusted (but continue)
    if ( (np.any(alpha_Had > 1.0001) or np.any(np.isinf(alpha_Had))) or 
         (np.any(alpha_E2P > 1.0001) or np.any(np.isinf(alpha_E2P))) ):
        warnings.warn("\n\n----------------- WARNING: scaling fractions exceed 1, might encounter infinity!\n\n")
 
    
    ## here comes the magic (part II); plugging in reference dat; DONE
    #******************************************************************************
    if verbose: 
        print("   --- Bias correction using source data...")
    Had_Hscaled  = np.multiply(alpha_Had, Href)
    E2P_Escaled  = np.multiply(alpha_E2P, Eref)
    #******************************************************************************
    
    ## for P-scaling, a bit more effort is required: 
    # 1. figure out relative P bias per arrival day (NOT per uptake day!)
    Pratio =  PrefTS / PtotTS # make sure this stays positive
    Pratio[Pratio==np.inf] = 0 # replace inf by 0 (happens if FLEX-P is zero)
    
    # 2.) now check how much E2P changed due to E-scaling already
    ## convert from mm to m3 as for P before, take areas into account
    E2P_Escaled_ts  = np.nansum(areas[:,0]*np.moveaxis(np.nansum(E2P_Escaled,axis=1), 1, 2), axis=(1,2))/1e3
    E2P_ts          = np.nansum(areas[:,0]*np.moveaxis(np.nansum(E2P,axis=1), 1, 2), axis=(1,2))/1e3
    f_Escaled       = np.divide(E2P_Escaled_ts, E2P_ts)
    ### Jessica: I believe the part below should be correct. It gives 1e-8 differences in E2P_EPscaled though, as a result of different f_Escaled
    ### ... to test, just uncomment the next lines
    #print(f_Escaled)
    #E2P_Escaled_ts  = np.nansum(np.multiply(areas,E2P_Escaled/1e3),axis=(1,2,3))
    #E2P_ts          = np.nansum(np.multiply(areas,E2P/1e3),axis=(1,2,3))
    #f_Escaled       = np.divide(E2P_Escaled_ts, E2P_ts)
    #print(f_Escaled)
    
    # 3.) alright, now calculate how much more scaling is needed to match P too
    f_remain = np.divide(Pratio, f_Escaled)
    
    #******************************************************************************
    ## swap axes to enable numpy broadcasting; 
    ## (a,b,c,d x b,c,d OK; a,b,c,d x a NOT OK)
    ## swap back and then store    
    if verbose: 
        print("   --- Bias correction using sink data...")
    E2P_Pscaled  = np.swapaxes(Pratio * np.swapaxes(E2P, 0, 3), 0, 3) 
    if verbose: 
        print("   --- Bias correction using source and sink data...")
    E2P_EPscaled = np.swapaxes(f_remain * np.swapaxes(E2P_Escaled, 0, 3), 0, 3) 
    #******************************************************************************
    
    
    ##--5. aggregate ##############################################################
    if verbose: 
        print(" * Writing final output... ")
    
    ## uptake time dimension is no longer needed!
    Had          = np.nansum(Had, axis=1)
    Had_scaled   = np.nansum(Had_Hscaled, axis=1)
    
    E2P          = np.nansum(E2P, axis=1)
    E2P_Escaled  = np.nansum(E2P_Escaled, axis=1)
    E2P_Pscaled  = np.nansum(E2P_Pscaled, axis=1)
    E2P_EPscaled = np.nansum(E2P_EPscaled, axis=1)
    
    ### plots that help to debug
    if fdebug: 
        debugmask(mask,maskval,mlat,mlon)
        alphascreener(alpha_Had, var='Had')
        alphascreener(alpha_E2P, var='E2P')
        plotpratio(f_remain,Pratio)
        basicplot(np.nanmean(Had, axis=0), lats, lons, 
                  title="raw Had, daily mean")
        basicplot(np.nanmean(Had_scaled, axis=0), lats, lons, 
                  title="H-scaled Had, daily mean")
        basicplot(np.nanmean(E2P, axis=0), lats, lons, 
                  title="raw E2P, daily mean")
        basicplot(np.nanmean(E2P_Escaled, axis=0), lats, lons, 
                  title="E-scaled E2P, daily mean")    
        basicplot(np.nanmean(E2P_Pscaled, axis=0), lats, lons, 
                  title="P-scaled E2P, daily mean")
        basicplot(np.nanmean(E2P_EPscaled, axis=0), lats, lons, 
                  title="E-P-scaled E2P, daily mean")
    
    ##--6. save output ############################################################
    if fwrite_netcdf:
        if verbose: 
            print(" * Writing final output... ")
        
        writefinalnc(ofile=ofile, fdate_seq=arrival_time, glon=lons, glat=lats, Had=Had, Had_Hs=Had_scaled, 
                 E2P=E2P, E2P_Es=E2P_Escaled, E2P_Ps=E2P_Pscaled, E2P_EPs=E2P_EPscaled, strargs=strargs, precision=precision)
