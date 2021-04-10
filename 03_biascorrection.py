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
           bcscale,
           pref_data,
           eref_data,
           href_data,
           faggbwtime,
           fdebug,
           fwrite_netcdf,
           fwrite_month,
           fwritestats,
           precision,
           strargs):
   
    ## SOME PRELIMINARY SETTINGS TO REDUCE OUTPUT
    ## suppressing warnings, such as
    #  invalid value encountered in true_divide
    #  invalid value encountered in multiply 
    if not fdebug:
        np.seterr(divide='ignore', invalid='ignore')
    # default values
    fwritewarning = False

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
        E2Psrt       = np.asarray(checknan(f['E2P'][:]))
        Hadsrt       = np.asarray(checknan(f['Had'][:]))
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
    P                   = -read_diagdata(opathD,ofile_base,ryyyy,uptake_time,var="P")
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
    if verbose: 
        print(" * Reading reference data...")
    
    print("Reading E")
    if eref_data == "eraint":
        Eref, reflats, reflons = eraloader_12hourly(var='e',
                     datapath=ipathR+"/evap_12hourly/E_1deg_",
                     maskpos=True,
                     maskneg=False,
                     uptake_years=uyears,
                     uptake_dates=uptake_dates, lats=lats, lons=lons)
    elif eref_data == "others":    
        # attention: data has to be on the correct grid and daily (or subdaily that can be summed up) and with the correct sign (all positive)
        Eref, reflats, reflons = get_reference_data(ipathR+"/gleam+oaflux+eraint_merge_1deg_e/", "evaporation", uptake_dates)
    gridcheck(totlats,reflats,totlons,reflons)

    # convert water fluxes from mm-->m3 to avoid area weighting in between
    Eref = convert_mm_m3(Eref, areas)
        
    print("Reading H")
    if href_data == "eraint":
        Href, reflats, reflons = eraloader_12hourly(var='sshf',
                     datapath=ipathR+"/sshf_12hourly/H_1deg_",
                     maskpos=True,
                     maskneg=False,
                     uptake_years=uyears,
                     uptake_dates=uptake_dates, lats=lats, lons=lons)
    elif href_data == "others":    
        # attention: data has to be on the correct grid and daily (or subdaily that can be summed up) and with the correct sign (all positive)
        Href, reflats, reflons = get_reference_data(ipathR+"/oaflux+eraint_merge_1deg_h/", "sensible heat flux", uptake_dates)
    gridcheck(totlats,reflats,totlons,reflons)
    
    print("Reading P")
    if pref_data == "eraint":
        Pref, reflats, reflons = eraloader_12hourly(var='tp',
                     datapath=ipathR+"/tp_12hourly/P_1deg_",
                     maskpos=False, # do NOT set this to True!
                     maskneg=True,
                     uptake_years=uyears,
                     uptake_dates=uptake_dates, lats=lats, lons=lons)
    elif pref_data == "others":
        # attention: data has to be on the correct grid and daily (or subdaily that can be summed up) and with the correct sign (all positive)
        Pref, reflats, reflons = get_reference_data(ipathR+"/mswep_1deg_p/", "precipitation", uptake_dates)
    gridcheck(totlats,reflats,totlons,reflons)

    # convert water fluxes from mm-->m3 to avoid area weighting in between
    Pref = convert_mm_m3(Pref, areas)
     
    ##--4. biascorrection #########################################################
    if verbose: 
        print(" * Starting bias correction...")
    
    ## P-scaling requires arrival region mask
    mask, mlat, mlon = maskgrabber(maskfile)
    # currently, only identical grids (mask = attribution = reference data) are supported...
    gridcheck(mlat,totlats,mlon,totlons)

    xla, xlo    = np.where(mask==maskval) # P[:,xla,xlo] is merely a 2D array... ;)
    ibgn        = np.where(uptake_time==arrival_time[0])[0][0] # only arrival days!
    
    ## preliminary checks
    if not fuseattp:
        # re-evaluate precip. data to check if it can be used (need daily data here because of upscaling in 02)
        fuseattp = check_attributedp(pdiag=Ptot[ibgn:,xla,xlo],pattr=E2P,veryverbose=veryverbose)
    
    #******************************************************************************
    ## (i) BIAS CORRECTING THE SOURCE
    #******************************************************************************
    if verbose: 
        print("   --- Bias correction using source data...")
    # quick consistency check
    consistencycheck(Had, Htot, bcscale, fdebug)
    consistencycheck(E2P, Etot, bcscale, fdebug)
    # calculate bias correction factor
    alpha_H     = calc_sourcebcf(ref=Href, diag=Htot, tscale=bcscale)
    alpha_E     = calc_sourcebcf(ref=Eref, diag=Etot, tscale=bcscale)
    # apply bias correction factor
    Had_Hcorrtd = np.multiply(alpha_H, Had)
    E2P_Ecorrtd = np.multiply(alpha_E, E2P)
    
    #******************************************************************************
    ## (ii) BIAS CORRECTING THE SINK (P only)
    #******************************************************************************
    if verbose: 
        print("   --- Bias correction using sink data...")
    # calculate (daily) bias correction factor
    if fuseattp:
        alpha_P     = calc_sinkbcf(ref=Pref[ibgn:,xla,xlo], att=E2P, tscale=bcscale)
        # perform monthly bias correction if necessary
        if np.all( np.nan_to_num(alpha_P) == 0):
            print("        * Monthly bias correction needed to match reference precipitation...")
            alpha_P     = calc_sinkbcf(ref=Pref[ibgn:,xla,xlo], att=E2P, tscale='monthly')
            fwritewarning= True
    else:    
        alpha_P     = calc_sinkbcf(ref=Pref[ibgn:,xla,xlo], att=Ptot[ibgn:,xla,xlo], tscale=bcscale)
        # perform monthly bias correction if necessary
        if np.all( np.nan_to_num(alpha_P) == 0):
            print("        * Monthly bias correction needed to match reference precipitation...")
            alpha_P     = calc_sinkbcf(ref=Pref[ibgn:,xla,xlo], att=Ptot[ibgn:,xla,xlo], tscale='monthly')
            fwritewarning= True
    # apply bias correction factor
    E2P_Pcorrtd     = np.swapaxes(alpha_P * np.swapaxes(E2P, 0, 3), 0, 3) 
    
    # additionally perform monthly bias correction of P if necessary
    if not checkpsum(Pref[ibgn:,xla,xlo], E2P_Pcorrtd, verbose=False):
        print("        * Additional monthly bias correction needed to match reference precipitation...")
        alpha_P     = calc_sinkbcf(ref=Pref[ibgn:,xla,xlo], att=E2P_Pcorrtd, tscale='monthly')
        E2P_Pcorrtd = np.swapaxes(alpha_P * np.swapaxes(E2P_Pcorrtd, 0, 3), 0, 3) 
        fwritewarning= True
    checkpsum(Pref[ibgn:,xla,xlo], E2P_Pcorrtd, verbose=verbose)
    
    #******************************************************************************
    ## (iii) BIAS CORRECTING THE SOURCE AND THE SINK (P only)
    #******************************************************************************
    if verbose: 
        print("   --- Bias correction using source and sink data...")
    # step 1: check how much E2P changed due to source-correction already
    alpha_P_Ecor    = calc_sinkbcf(ref=E2P_Ecorrtd, att=E2P, tscale=bcscale)
    # step 2: calculate how much more correction is needed to match sink 
    alpha_P         = calc_sinkbcf(ref=Pref[ibgn:,xla,xlo], att=E2P_Pcorrtd, tscale=bcscale)
    # perform monthly bias correction if necessary
    if np.all( np.nan_to_num(alpha_P) == 0):
        print("        * Monthly bias correction needed to match reference precipitation...")
        alpha_P     = calc_sinkbcf(ref=Pref[ibgn:,xla,xlo], att=E2P_Pcorrtd, tscale='monthly')
        fwritewarning= True
    # step 3: calculate adjusted bias correction factor
    alpha_P_res     = np.divide(alpha_P, alpha_P_Ecor)
    E2P_EPcorrtd    = np.swapaxes(alpha_P_res * np.swapaxes(E2P_Ecorrtd, 0, 3), 0, 3) 
    
    # additionally perform monthly bias correction of P if necessary
    if not checkpsum(Pref[ibgn:,xla,xlo], E2P_EPcorrtd, verbose=False):
        print("        * Additional monthly bias correction needed to match reference precipitation...")
        alpha_P_res = calc_sinkbcf(ref=Pref[ibgn:,xla,xlo], att=E2P_EPcorrtd, tscale='monthly')
        E2P_EPcorrtd= np.swapaxes(alpha_P_res * np.swapaxes(E2P_EPcorrtd, 0, 3), 0, 3)
        fwritewarning= True
    checkpsum(Pref[ibgn:,xla,xlo], E2P_EPcorrtd, verbose=verbose)

    # save some data in case debugging is needed
    if fdebug:
        frac_E2P = calc_alpha(E2P,Etot)
        frac_Had = calc_alpha(Had,Htot)
    
    ##--5. aggregate ##############################################################
    ## aggregate over uptake time (uptake time dimension is no longer needed!)
    aHad          = np.nansum(Had, axis=1)
    aHad_Hcorrtd  = np.nansum(Had_Hcorrtd, axis=1)
    aE2P          = np.nansum(E2P, axis=1)
    aE2P_Ecorrtd  = np.nansum(E2P_Ecorrtd, axis=1)
    aE2P_Pcorrtd  = np.nansum(E2P_Pcorrtd, axis=1)
    aE2P_EPcorrtd = np.nansum(E2P_EPcorrtd, axis=1)
    # free up memory if backward time not needed anymore... 
    if faggbwtime:
        del(Had,Had_Hcorrtd,E2P,E2P_Ecorrtd,E2P_Pcorrtd,E2P_EPcorrtd)

    if fwritestats:
        # write some additional statistics about P-biascorrection before converting back to mm
        writestats_03(sfile,Pref,aE2P,aE2P_Ecorrtd,aE2P_Pcorrtd,aE2P_EPcorrtd,aHad,aHad_Hcorrtd,xla,xlo,ibgn)

    ##--6. unit conversion ##############################################################
    # and convert water fluxes back from m3 --> mm
    if not faggbwtime:
        E2P           = convert_m3_mm(E2P,areas)
        E2P_Ecorrtd   = convert_m3_mm(E2P_Ecorrtd,areas)
        E2P_Pcorrtd   = convert_m3_mm(E2P_Pcorrtd,areas)
        E2P_EPcorrtd  = convert_m3_mm(E2P_EPcorrtd,areas)
    if fdebug or faggbwtime:    
        aE2P          = convert_m3_mm(aE2P,areas)
        aE2P_Ecorrtd  = convert_m3_mm(aE2P_Ecorrtd,areas)
        aE2P_Pcorrtd  = convert_m3_mm(aE2P_Pcorrtd,areas)
        aE2P_EPcorrtd = convert_m3_mm(aE2P_EPcorrtd,areas)
    
    ##--7. debugging needed? ######################################################
    if fdebug:
        print(" * Creating debugging file")
        writedebugnc(opath+"/debug.nc",arrival_time,uptake_time,lons,lats,maskbymaskval(mask,maskval),
                mask3darray(Pref[ibgn:,:,:],xla,xlo),mask3darray(Ptot[ibgn:,:,:],xla,xlo),
                convert_mm_m3(aE2P,areas),convert_mm_m3(aE2P_Ecorrtd,areas),
                convert_mm_m3(aE2P_Pcorrtd,areas),convert_mm_m3(aE2P_EPcorrtd,areas),
                np.nan_to_num(frac_E2P),
                np.nan_to_num(frac_Had),
                alpha_P,np.nan_to_num(alpha_P_Ecor),np.nan_to_num(alpha_P_res),
                np.nan_to_num(alpha_E),np.nan_to_num(alpha_H),
                strargs,precision)
    
    ##--8. write final output ############################################################
    if verbose: 
        print(" * Writing final output... ")
        
    if fwrite_netcdf:
        # get attributes from attribution file and modify
        attrdesc    = getattr(nc4.Dataset(attrfile),"description") + "; " + strargs
        biasdesc    = attrdesc.replace("02_attribution","03_biascorrection") 

        # write to netcdf
        if faggbwtime:
            writefinalnc(ofile=ofile, 
                        fdate_seq=arrival_time, udate_seq=np.nan, 
                        glon=lons, glat=lats, 
                        Had=aHad, 
                        Had_Hs=aHad_Hcorrtd, 
                        E2P=aE2P, 
                        E2P_Es=aE2P_Ecorrtd, 
                        E2P_Ps=aE2P_Pcorrtd, 
                        E2P_EPs=aE2P_EPcorrtd, 
                        strargs=biasdesc, 
                        precision=precision,
                        fwrite_month=fwrite_month)
        if not faggbwtime:
            writefinalnc(ofile=ofile, 
                        fdate_seq=arrival_time, udate_seq=utime_srt, 
                        glon=lons, glat=lats, 
                        Had=reduce4Darray(Had,veryverbose), 
                        Had_Hs=reduce4Darray(Had_Hcorrtd,veryverbose), 
                        E2P=reduce4Darray(E2P,veryverbose), 
                        E2P_Es=reduce4Darray(E2P_Ecorrtd,veryverbose), 
                        E2P_Ps=reduce4Darray(E2P_Pcorrtd,veryverbose), 
                        E2P_EPs=reduce4Darray(E2P_EPcorrtd,veryverbose), 
                        strargs=biasdesc, 
                        precision=precision,
                        fwrite_month=fwrite_month)
    if fwritewarning:
        wfile = opath+"/"+str(ofile_base)+"_biascor-attr_r"+str(ryyyy)[-2:]+"_"+str(ayyyy)+"-"+str(am).zfill(2)+"_WARNING.csv"
        writewarning(wfile)

    if os.path.exists(ofile):
        print("Removing "+str(attrfile)+" ...")
        os.remove(attrfile)
