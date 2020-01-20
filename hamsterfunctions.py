#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
NON-METEOROLOGICAL FUNCTIONS USED BY MODULES
01_attribution, 02_diagnosis, 03_bias-correction

@author: jessica and dominik

To execute interactively: 
> exec(open("./hamsterfunctions.py").read())

"""

def PBL_check(z, h, seth, diagnosis):

    if diagnosis == 'KAS':
        h[h<seth] = seth
        before_inside = np.logical_or( z[1:] < h[:-1], z[1:]  < h[1:])  
        after_inside  = np.logical_or(z[:-1] < h[:-1], z[:-1] < h[1:])  
        change_inside = np.logical_and(before_inside, after_inside)
    elif diagnosis == 'SOD':   
        change_inside = ((z[1:]+z[:-1])/2) < 1.5*((h[1:]+h[:-1])/2)
        # NOTE: factor 1.5 is hardcoded      

    return change_inside

def trajparceldiff(pvals, meval):
    # difference 
    if meval in ['diff']:
        dpval   = pvals[:-1] - pvals[1:]
    # mean
    if meval in ['mean']:
        dpval   = (pvals[:-1] + pvals[1:])/2
    # max
    if meval in ['max']:
        dpval   = np.amax(pvals[:-1], pvals[1:])
    return(dpval)


def linear_discounter(v, min_gain, min_loss):
    """
    ========= inputs =========
    v:        vector extending into past; v[0] is most recent, v[-1] least recent
    min_gain: minimum gain (dv/dt) to be considered
    min_loss: minimum loss to be considered for discounting, must be <= 0
    ==========================
    """

    ## compute dv/dt, prepare dv_disc
    dv = v[:-1] - v[1:]
    dv_disc = np.zeros(shape=len(dv))
    ## get indices of gains and losses
    idx_gains  = np.where(dv >= min_gain)[0]
    idx_losses = np.where(dv <= min_loss)[0]

    for idx in idx_gains:
        ## determine relevant indices (losses only!)
        idx_consider = idx_losses[np.where(idx_losses < idx)]

        ## skip current iteration if no moisture loss occurs between then and t=0
        if len(idx_consider)==0:
            dv_disc[idx] = dv[idx]
            continue

        ## create vector with fractions,
        frc_vec = v[idx_consider]/v[idx_consider+1]
        ## then multiply current dv with product of vector (order irrelevant)
        dv_disc[idx] = np.prod(frc_vec)*dv[idx]
    return(dv_disc)


def writenc4D(ofile,fdate_seq,fuptdate_seq,glon,glat,ary_etop,ary_heat,ary_npart):
    if verbose:
        print(" * Writing netcdf output...")
                
    # convert date objects to datetime objects if necessary
    if type(fdate_seq[0]) == datetime.date:
        for idt in range(len(fdate_seq)):
            fdate_seq[idt]    = datetime.datetime(fdate_seq[idt].year, fdate_seq[idt].month, fdate_seq[idt].day)
    if type(fuptdate_seq[0]) == datetime.date:
        for idt in range(len(fuptdate_seq)):
            fuptdate_seq[idt] = datetime.datetime(fuptdate_seq[idt].year, fuptdate_seq[idt].month, fuptdate_seq[idt].day)

    # delete nc file if it is present (avoiding error message)
    try:
        os.remove(ofile)
    except OSError:
        pass
        
    # create netCDF4 instance
    nc_f = nc4.Dataset(ofile,'w', format='NETCDF4')
    
    # create dimensions 
    nc_f.createDimension('arrival-time', len(fdate_seq))
    nc_f.createDimension('uptake-time', len(fuptdate_seq))
    nc_f.createDimension('lat', glat.size)
    nc_f.createDimension('lon', glon.size)
    
    # create variables
    atimes              = nc_f.createVariable('arrival-time', 'i4', 'arrival-time')
    utimes              = nc_f.createVariable('uptake-time', 'i4', 'uptake-time')
    latitudes           = nc_f.createVariable('lat', 'f4', 'lat')
    longitudes          = nc_f.createVariable('lon', 'f4', 'lon')
    heats               = nc_f.createVariable('H', 'f4', ('arrival-time','uptake-time','lat','lon'))
    etops               = nc_f.createVariable('E2P', 'f4', ('arrival-time','uptake-time','lat','lon'))
    nparts              = nc_f.createVariable('n_part', 'f4', ('arrival-time','uptake-time','lat','lon'))
    
    # set attributes
    nc_f.description    = "FLEXPART: 02_attribution of advected surface sensible heat and evaporative moisture resulting in precipitation"
    atimes.units        = 'hours since 1900-01-01 00:00:00'
    atimes.calendar     = 'Standard' 
    utimes.units        = 'hours since 1900-01-01 00:00:00'
    utimes.calendar     = 'Standard' 
    latitudes.units     = 'degrees_north'
    longitudes.units    = 'degrees_east'
    heats.units         = 'W m-2'
    heats.long_name	= 'surface sensible heat flux'
    etops.units         = 'mm'
    etops.long_name	= 'evaporation resulting in precipitation'
    nparts.units        = 'int'
    nparts.long_name    = 'number of parcels (mid pos.)'
    
    # write data
    atimes[:]           = nc4.date2num(fdate_seq, atimes.units, atimes.calendar)
    utimes[:]           = nc4.date2num(fuptdate_seq, utimes.units, utimes.calendar)
    longitudes[:]       = glon
    latitudes[:]        = glat
    heats[:]            = ary_heat[:]
    etops[:]            = ary_etop[:]     
    nparts[:]           = ary_npart[:]
        
    # close file
    nc_f.close()
        
    print("\n===============================================================")
    print("\n Successfully written: "+ofile+" !")
    print("\n===============================================================")
