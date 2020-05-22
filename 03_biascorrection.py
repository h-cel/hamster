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
           fdebug,
           fwrite_netcdf,
           precision,
           strargs):
    
    # Consistency checks
    if mode=="oper" and precision=="f4":
        precision = "f8"
        print("Single precision should only be used for testing. Reset to double-precision.")

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

    ##--1. load attribution data; grab all uptake days ############################    
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
    ## NOTE: this was first coded with the intention of being able to deal with
    ## basically any sort of attribution input, i.e. multi-monthly data..
    ## should be cleaned up since we will, as discussed on 04-03-2020, stick
    ## to the monthly format

    ayears  =  np.asarray([ut.year  for ut in uptake_time])
    amonths =  np.asarray([ut.month for ut in uptake_time])
    uset = np.unique(np.column_stack((ayears, amonths)), axis=0)
    uyears, umonths = uset[:,0], uset[:,1]
    
    for jj in range(umonths.size):
        uyr  = str(uyears[jj])
        umon = str(umonths[jj])
        diagfile = str(opathD)+"/"+str(ofile_base)+"_diag_r"+str(ryyyy)[-2:]+"_"+str(uyr)+"-"+umon.zfill(2)+".nc"
        with nc4.Dataset(diagfile, mode="r") as f:
            Ex     = f['E'][:]
            Px     = f['P'][:]
            Hx     = f['H'][:]
            timex  = nc4.num2date(f['time'][:], f['time'].units, f['time'].calendar)
        
        ## concatenate 'em!
        if umon == str(umonths[0]):
            E = np.copy(Ex)
            P = np.copy(Px)
            H = np.copy(Hx)
            time = np.copy(timex)
        else:
            E = np.concatenate((E, Ex), axis=0)
            P = np.concatenate((P, Px), axis=0)
            H = np.concatenate((H, Hx), axis=0)
            time = np.concatenate((time, timex))
    
        ## check if coordinates are the same
        with nc4.Dataset(diagfile, mode="r") as f:
            totlats = f['lat'][:]
            totlons = f['lon'][:]
            
        if not np.array_equal(lats, totlats) or not np.array_equal(lons, totlons):
            raise SystemExit("--- ERROR: your grids aren't identical...")
    
    ## must check if data comes in daily resolution; fix if not
    dates   = np.asarray([datetime.date(it.year, it.month, it.day) for it in time])
    udates = np.unique(dates)
    if udates.size != dates.size:
        ## might need to redo dates in case uptakes are stored on arrival times...
        if time[0].hour in [0, 6, 12, 18]:
            ## simple fix, subtract 3 hours
            time = np.asarray([t - datetime.timedelta(hours=3) for t in time])
            dates   = np.asarray([datetime.date(it.year, it.month, it.day) for it in time])
            udates = np.unique(dates)
        elif time[0].hour in [3, 9, 15, 21]:
            ## NOTE: this is the new norm! retain "old style" for now, though    
            pass 
        
        Etot = np.zeros(shape=(time.size, lats.size, lons.size))
        Ptot = np.zeros(shape=(time.size, lats.size, lons.size))
        Htot = np.zeros(shape=(time.size, lats.size, lons.size))
        
        ## this isn't fast or elegant, but works for literally anything sub-daily
        for i in range(udates.size):
            
            iud = udates[i]       
            sel = np.where(dates == iud)[0]
     
            ## TODO: clean up; there should be a check whether 4 files are present, imo
            if sel.size != 4:
                warnings.warn("\n\n----------------- WARNING: this should NEVER OCCUR; daily aggregation IMPROPER (files missing!)\n\n")
           
            Etot[i,:,:] = np.nansum(E[sel, :, :], axis=0) # nan likely unnecessary
            Ptot[i,:,:] = np.nansum(P[sel, :, :], axis=0)
            Htot[i,:,:] = np.nanmean(H[sel, :, :], axis=0) # NOTE: mean here!
    
    else: ## NOTE: preparing for DAILY OUTPUT from 01_diagnosis
        ## neat, we can just proceed
        Etot = E
        Ptot = P
        Htot = H
        
    ## only keep what is really needed
    if date_bgn not in udates: 
        raise SystemExit("\n !!! ERROR: INPUT DATA MISSING: date "+str(date_bgn)+" not available as output from 01_diagnosis! Aborting here. !!!\n")
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
    del(E, P, H, Ex, Px, Hx)

    ##--3. load reference data ####################################################
    """
    this part is STRICTLY CODED FOR (12-hourly) ERA-INTERIM only (so far),
    and HARDCODED too
    """
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
    
    ## P-scaling requires arrival region mask
    with nc4.Dataset(maskfile) as f:
        mask = f['mask'][:]
        mlat = f['lat'][:]
        mlon = f['lon'][:]   
    if fdebug: 
        basicplot(mask, mlat, mlon, title="content of mask file (all values plotted)")
        mask[(mask>0) & (mask!=maskval)] = 0
        basicplot(mask, mlat, mlon, title="mask used for bias-correction")
    
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
 
    ## have a look if you're curious
    if fdebug:
        alphascreener(alpha_Had, var='Had')
        alphascreener(alpha_E2P, var='E2P')
    
    ## here comes the magic (part II); plugging in reference dat; DONE
    #******************************************************************************
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
    
    if fdebug:
        plt.figure
        plt.plot(f_remain, label="remaining scaling factor")
        plt.plot(Pratio, label="P ratio (REF/FLEX)")
        plt.legend()
        plt.show()    
    
    #******************************************************************************
    ## swap axes to enable numpy broadcasting; 
    ## (a,b,c,d x b,c,d OK; a,b,c,d x a NOT OK)
    ## swap back and then store    
    E2P_Pscaled  = np.swapaxes(Pratio * np.swapaxes(E2P, 0, 3), 0, 3) 
    E2P_EPscaled = np.swapaxes(f_remain * np.swapaxes(E2P_Escaled, 0, 3), 0, 3) 
    #******************************************************************************
    
    
    ##--5. aggregate ##############################################################
    
    ## uptake time dimension is no longer needed!
    Had          = np.nansum(Had, axis=1)
    Had_scaled   = np.nansum(Had_Hscaled, axis=1)
    
    E2P          = np.nansum(E2P, axis=1)
    E2P_Escaled  = np.nansum(E2P_Escaled, axis=1)
    E2P_Pscaled  = np.nansum(E2P_Pscaled, axis=1)
    E2P_EPscaled = np.nansum(E2P_EPscaled, axis=1)
    
    if fdebug:
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
        writefinalnc(ofile=ofile, fdate_seq=arrival_time, glon=lons, glat=lats, Had=Had, Had_Hs=Had_scaled, 
                 E2P=E2P, E2P_Es=E2P_Escaled, E2P_Ps=E2P_Pscaled, E2P_EPs=E2P_EPscaled, strargs=strargs, precision=precision)
