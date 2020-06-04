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
           fuseattp,
           fdebug,
           fwrite_netcdf,
           fwritestats,
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
    ## additional statistic output files includes P validation data (*.csv)
    if fwritestats:
        sfilename = str(ofile_base)+"_biascor-attr_r"+str(ryyyy)[-2:]+"_"+str(ayyyy)+"-"+str(am).zfill(2)+"_stats.csv"
        sfile     = opath+"/"+sfilename

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
        print(" ! using attribution data to bias-correct P: \t" +str(fuseattp))
        print(" ! writing netcdf output: \t")
        print("\t"+str(ofile))
        if fwritestats:
            print(" ! precipitation statistics in: \t")
            print("\t"+str(sfile))
        print("\n============================================================================================================")
        print("\n")

    ##--1. load attribution data; grab all uptake days ############################    
    if verbose: 
        print(" * Reading attribution data...")
        if veryverbose:
            print("   --- file: "+str(attrfile))

    with nc4.Dataset(attrfile, mode="r") as f:
        E2Psrt       = np.asarray(f['E2P'][:])
        Hadsrt       = np.asarray(f['H'][:])
        arrival_time = nc4.num2date(f['time'][:], f['time'].units, f['time'].calendar)
        utime_srt    = np.asarray(f['level'][:])
        uptake_time  = udays2udate(arrival_time,utime_srt)
        uptake_dates = cal2date(uptake_time)
        uyears       = np.unique(date2year(uptake_time))
        lats         = np.asarray(f['lat'][:])
        lons         = np.asarray(f['lon'][:])
        areas        = 1e6*np.nan_to_num(gridded_area_exact(lats, res=abs(lats[1]-lats[0]), nlon=lons.size))[:,0]
    # expand uptake dimension to dates (instead of backward days)
    E2P = expand4Darray(E2Psrt,arrival_time,utime_srt,veryverbose)
    Had = expand4Darray(Hadsrt,arrival_time,utime_srt,veryverbose)
    # convert water fluxes from mm-->m3
    E2P = convert_mm_m3(E2P, areas)

    # clean up
    del(E2Psrt, Hadsrt)
    
    ##--2. load diagnosis data ####################################################
    if verbose: 
        print(" * Reading diagnosis data...")
    
    # read concatenated data
    totlats, totlons    = read_diagdata(opathD,ofile_base,ryyyy,uptake_time,var="grid")
    gridcheck(lats,totlats,lons,totlons)
    ftime               = read_diagdata(opathD,ofile_base,ryyyy,uptake_time,var="time")
    fdays               = np.unique(cal2date(ftime))
    E                   = read_diagdata(opathD,ofile_base,ryyyy,uptake_time,var="E")
    P                   = read_diagdata(opathD,ofile_base,ryyyy,uptake_time,var="P")
    H                   = read_diagdata(opathD,ofile_base,ryyyy,uptake_time,var="H")
    # convert water fluxes from mm-->m3 to avoid area weighting in between
    E                   = convert_mm_m3(E, areas)
    P                   = convert_mm_m3(P, areas)

    # make sure we use daily aggregates
    if fdays.size != ftime.size:
        Etot    = convert2daily(E,ftime,fagg="sum")
        Ptot    = convert2daily(P,ftime,fagg="sum")
        Htot    = convert2daily(H,ftime,fagg="mean")
    else: 
        Etot = E
        Ptot = P
        Htot = H
        
    ## only keep what is really needed (P is stored analogous to E and H for consistency)
    datecheck(uptake_dates[0],fdays)
    ibgn = np.where(fdays==uptake_dates[0])[0][0]
    iend = np.where(fdays==uptake_dates[-1])[0][0]
    Etot = Etot[ibgn:iend+1,:,:]
    Ptot = Ptot[ibgn:iend+1,:,:]
    Htot = Htot[ibgn:iend+1,:,:]
    fdates = fdays[ibgn:iend+1]  
    ## make sure we grabbed the right data
    if not np.array_equal(uptake_dates, fdates):
        raise SystemExit("---- hold your horses; datetime matching failed!")

    ## clean up
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
    # convert water fluxes from mm-->m3 to avoid area weighting in between
    Eref = convert_mm_m3(Eref, areas)
        
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
    # convert water fluxes from mm-->m3 to avoid area weighting in between
    Pref = convert_mm_m3(Pref, areas)
     
    ##--4. biascorrection #########################################################
    if verbose: 
        print(" * Starting bias correction...")
    
    ## P-scaling requires arrival region mask
    with nc4.Dataset(maskfile) as f:
        mask = f['mask'][:]
        mlat = f['lat'][:]
        mlon = f['lon'][:]   
    xla, xlo    = np.where(mask==maskval) # P[:,xla,xlo] is merely a 2D array... ;)
    ibgn        = np.where(uptake_time==arrival_time[0])[0][0] # only arrival days!
    
    ## preliminary checks
    if not fuseattp:
        # re-evaluate precip. data to check if it can be used
        fuseattp = checkprec(pdiag=Ptot[ibgn:,xla,xlo],pattr=E2P,veryverbose=veryverbose)
    
    #******************************************************************************
    ## (i) BIAS CORRECTING THE SOURCE
    #******************************************************************************
    if verbose: 
        print("   --- Bias correction using source data...")
    # calculate source fraction contribution (FLEXPART data only)
    alpha_Had = calc_alpha(Had, Htot)
    alpha_E2P = calc_alpha(E2P, Etot)
    # and use reference data for bias-correcting the totals 
    Had_Hscaled  = np.multiply(alpha_Had, Href)
    E2P_Escaled  = np.multiply(alpha_E2P, Eref)
    
    #******************************************************************************
    ## (ii) BIAS CORRETING THE SINK (P only)
    #******************************************************************************
    if verbose: 
        print("   --- Bias correction using sink data...")
    # sum up precpitation (arrival days) over mask only
    PrefTS      = np.nansum(Pref[ibgn:,xla,xlo], axis=1)
    if fuseattp:
        PtotTS  = -np.nansum(E2P,axis=(1,2,3))
    else:    
        PtotTS  = np.nansum(Ptot[ibgn:,xla,xlo], axis=1)
    # switch to monthly bias correction of P if necessary
    # attention: still writing out daily data though (days won't match!)
    fusemonthly = needmonthlyp(pdiag=PtotTS,pref=PrefTS)    
    if fusemonthly:
        ndays   = PtotTS.shape[0] 
        PtotTS  = np.repeat(np.nansum(PtotTS),ndays)
        PrefTS  = np.repeat(np.nansum(PrefTS),ndays)
    # calculate bias correction fractor
    Pratio      = PrefTS / PtotTS # make sure this stays positive
    Pratio[Pratio==np.inf] = 0 # replace inf by 0 (happens if FLEX-P is zero)
    E2P_Pscaled = np.swapaxes(Pratio * np.swapaxes(E2P, 0, 3), 0, 3) 
    if round(np.nansum(E2P_Pscaled),4) != round(np.nansum(-Pref[ibgn:,xla,xlo]),4):
        print("  --- OOOPS... something must be wrong in the biascorrection of P.")
    
    #******************************************************************************
    ## (iii) BIAS CORRETING THE SOURCE AND THE SINK (P only)
    #******************************************************************************
    if verbose: 
        print("   --- Bias correction using source and sink data...")
    ## step 1: check how much E2P changed due to E-scaling already
    E2P_Escaled_ts  = np.nansum(E2P_Escaled,axis=(1,2,3))
    E2P_ts          = np.nansum(E2P,axis=(1,2,3))
    f_Escaled       = np.divide(E2P_Escaled_ts, E2P_ts)
    # step 2: calculate how much more scaling is needed to match P too 
    f_remain = np.divide(Pratio, f_Escaled)
    E2P_EPscaled = np.swapaxes(f_remain * np.swapaxes(E2P_Escaled, 0, 3), 0, 3) 
    if round(np.nansum(E2P_EPscaled),4) != round(np.nansum(-Pref[ibgn:,xla,xlo]),4):
        print("  --- OOOPS... something must be wrong in the biascorrection of E or P.")
    
    ##--5. aggregate ##############################################################
    ## aggregate over uptake time (uptake time dimension is no longer needed!)
    Had          = np.nansum(Had, axis=1)
    Had_scaled   = np.nansum(Had_Hscaled, axis=1)
    E2P          = np.nansum(E2P, axis=1)
    E2P_Escaled  = np.nansum(E2P_Escaled, axis=1)
    E2P_Pscaled  = np.nansum(E2P_Pscaled, axis=1)
    E2P_EPscaled = np.nansum(E2P_EPscaled, axis=1)
    if fwritestats:
        # write some additional statistics about P-biascorrection before converting back to mm
        writestats_03(sfile,Pref,E2P,E2P_Escaled,E2P_Pscaled,E2P_EPscaled,xla,xlo,ibgn)

    # and convert water fluxes back from m3 --> mm
    E2P          = convert_m3_mm(E2P,areas)
    E2P_Escaled  = convert_m3_mm(E2P_Escaled,areas)
    E2P_Pscaled  = convert_m3_mm(E2P_Pscaled,areas)
    E2P_EPscaled = convert_m3_mm(E2P_EPscaled,areas)
    
    ##--6. debugging needed? ######################################################
    ### plots that help to debug
    if fdebug:
        print("   --- Creating debugging file")
        # this is where debugging function will be called
    
    ##--7. save output ############################################################
    if verbose: 
        print(" * Writing final output... ")
        
    if fwrite_netcdf:
        # get attributes from attribution file and modify
        attrdesc    = getattr(nc4.Dataset(attrfile),"description")
        biasdesc    = attrdesc.replace("02_attribution","03_biascorrection") 

        # write to netcdf
        writefinalnc(ofile=ofile, fdate_seq=arrival_time, glon=lons, glat=lats, Had=Had, Had_Hs=Had_scaled, 
                 E2P=E2P, E2P_Es=E2P_Escaled, E2P_Ps=E2P_Pscaled, E2P_EPs=E2P_EPscaled, strargs=biasdesc, precision=precision)
