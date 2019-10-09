#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MAIN SCRIPT MOISTURE AND HEAT DIAGNOSIS
@author: dominik and jessica

To execute interactively: 
> exec(open("./main.py").read())

"""

## METEOROLOGY FUNCTIONS

def calc_pres(dens,temp):
    return(dens*RSPECIFIC*temp)

def dist_on_sphere(lat1,lon1,lat2,lon2):
    """
    INPUT 
        - lon1 :    longitude of point 1
        - lon2 :    longitude of point 2
        - lat1 :    latitude of point 1
        - lat2 :    latitude of point 2

    ACTION
    This function calculates the linear distance between two points
    defined by lat1,lon1 and lat2,lon2 (in °, lat [-90,90], lon[-180,180])
    on Earth's surface, relying on the following assumptions:
        - Earth is a perfect sphere,
        - the radius is precisely 6371 km
    The calculation is based on a transformation to cartesian coordinates,
    followed by computing the dot product between two vectors spanned
    by the two points on Earth's surface and the center, respectively.

    DEPENDENCIES
        all from module math: sin, cos, acos
        
    RETURNS:
        - dist :    distance between locations in kilometers.

    """

    ## first off, must adapt lat/lon for spherical coordinates
    lat1 = 90 - lat1
    lat2 = 90 - lat2
    lon1 = 180 + lon1
    lon2 = 180 + lon2

    ## second, must obtain angles in radian from (arc) degree
    la1 = (PI/180)*lat1
    la2 = (PI/180)*lat2
    lo1 = (PI/180)*lon1
    lo2 = (PI/180)*lon2

    ## third, convert to cartesian coordinates
    ## use unit sphere for simplicity (r=1)
    x1 = sin(la1) * cos(lo1)
    y1 = sin(la1) * sin(lo1)
    z1 = cos(la1)

    x2 = sin(la2) * cos(lo2)
    y2 = sin(la2) * sin(lo2)
    z2 = cos(la2)

    ## fourth, calculate dot product and angle
    dotp = (x1*x2 + y1*y2 + z1*z2)
    angle = acos(dotp)

    ## fifth, calculate distance
    dist = angle * EARTHRADIUS   # pi cancels out

    return(dist)

def gridded_area_exact(lats_centr, res, nlon):
    """
    INPUT
        - lats_centr :  latitute of the center [deg N], between -90 to +90 (float or 1D np ary)
        - res :         regular (!) grid resolution [deg]
        - nlon :        to return array of dim (nlat x nlon)

    RETURNS
        - area :        EXACT gridded area, shape as defined by input [km^2] as array of dimension (nlon x nlat)

    ACTION
        based on the grid spacing inferred by lats & lons,
        areas of all grid cells specified by indices and coordinates (lats)
        are computed and summed up in the end.
        Since the spacing between grid cell border longitudes ALWAYS remains
        the same, that is, the resolution in degrees, this function does NOT
        rely on any longitudinal input.

    CAUTION
        1.) this will NOT work properly close to the poles!
        2.) this is based on input coordinates referring to CENTROIDS
        3.) only for regular grids (independent from longitude)
        4.) will produce CRAP if / results in integer division (Python 2.7) !

    NOTES
        - original formula, using degrees (caution -- np.sin requires radian)
            A = (pi/180)R^2 |sin(lat1)-sin(lat2)| |lon1-lon2|
          obtained from:
            pmel.noaa.gov/maillists/tmap/ferret_users/fu_2004/msg00023.html

    """

    ## make use of numpy vectorization
    lats1 = (lats_centr+(res/2))*np.pi/180 # np.sin requires radians
    lats2 = (lats_centr-(res/2))*np.pi/180
    areas = (np.pi/180)*(EARTHRADIUS**2)*np.abs(np.sin(lats1)-np.sin(lats2))*res

    ## overwrite any areas of 0 (at the poles) with np.NaN to prevent problems
    try:
        areas[np.where(areas==0.)] = np.NaN # only works for arrays
    except TypeError:
        pass # simply ignore if it's a float
    # return array of dimension nlat x nlon
    ary_area = np.swapaxes(np.tile(areas, (nlon,1)), 0,1)
    return(ary_area)

def midpoint_on_sphere(lat1,lon1,lat2,lon2):

    """
    INPUT
        - lon1 :    longitude of point 1
        - lon2 :    longitude of point 2
        - lat1 :    latitude of point 1
        - lat2 :    latitude of point 2
        * must come as lons [-180 .. 180], lats [-90 .. 90]
    
    WARNING 
        the coordinate transformation is probably correct,
        but the last transformations were done empirically
        (and not based on logic). 
        A number of checks has not revealed any issues.

    DEPENDENCIES:
        numpy, various trigonometric functions from module math
        
    OUTPUT: 
        - lon_mid :     geographical center longitude (lat/lon)
        - lat_mid :     geographical center latitude (lat/lon)

    """

    ## first off, must adapt lat/lon for spherical coordinates
    lat1 = 90 - lat1
    lat2 = 90 - lat2
    lon1 = 180 + lon1
    lon2 = 180 + lon2

    ## second, must obtain angles in radian from (arc) degree
    lat1 = (PI/180)*lat1
    lat2 = (PI/180)*lat2
    lon1 = (PI/180)*lon1
    lon2 = (PI/180)*lon2

    ## third, convert to cartesian coordinates
    ## use unit sphere for simplicity (r=1)
    x1 = sin(lat1) * cos(lon1)
    y1 = sin(lat1) * sin(lon1)
    z1 = cos(lat1)

    x2 = sin(lat2) * cos(lon2)
    y2 = sin(lat2) * sin(lon2)
    z2 = cos(lat2)

    ## fourth, add vector to obtain "middle" vector, scale back to length 1
    middle_vec = [x1+x2, y1+y2, z1+z2]
    mv = np.asarray(middle_vec)/sqrt(sum(np.asarray(middle_vec)**2))

    ## TRANSFORM BACK TO RADIAN
    theta = atan(sqrt((mv[0]**2)+(mv[1]**2))/(mv[2]))
    phi   = atan2(mv[1],mv[0]) # is "same" as atan(y/x), but NEED atan2 (sign!)

    ## convert these mid coordinates back to (arc) degree & shift
    if theta > 0:
        lat_mid = 90 - theta*(180/PI)
    else:
        lat_mid = -(90 + theta*(180/PI))
    if phi > 0:
        lon_mid = (phi*(180/PI)-180)
    else:
        lon_mid = 180 + phi*(180/PI)

    if (lon_mid>179.5): lon_mid -= 360    # now shift all coords that otherwise would be allocated to +180 deg to - 180

    return(lat_mid, lon_mid)

def q2rh(q_kgkg,p_Pa,T_K):
    """
    INPUT
        - q_kgkg: kg kg-1,  specific humidity, float or vector
        - p_Pa: Pa,         pressure, float or vector
        - T_K: K,           temperature, float or vector
    
    ACTION
        Converting specific to relative humidity following Bolton (1980), see:
        https://www.eol.ucar.edu/projects/ceop/dm/documents/refdata_report/eqns.html

    OUTPUTS
        - RH: %        relative humidity in percent.
    """
    # convert q into e, vapor pressure
    e = q_kgkg*p_Pa/(0.622+0.378*q_kgkg)

    # compute saturation vapor pressure, must convert T to °C
    es = 611.2*np.exp(17.67*(T_K-TREF)/(T_K-TREF+243.5))

    # return relative humidity
    return(1e2*e/es)

def dqsdT(p_hPa, T_degC):
    """
    INPUTS
        - p_hPa : hPa,              pressure
        - T_degC : degree Celsius,  temperature
    
    ACTION
        1. CC eq approximation by Bolton (1980) differentiated by T, obtained from
        https://www.wolframalpha.com/input/?i=differentiate+a*exp(b*x%2F(x%2Bc))
        2. convert saturation vapor pressure to saturation specific humidity
        NOTE: pressures (e=water vapor, p=total) must be in hPa
        returns specific humidity in kg/kg
        https://archive.eol.ucar.edu/projects/ceop/dm/documents/refdata_report/eqns.html
    
    OUTPUTS
        - dqdT: kg kg-1 K-1     the change of specific humidity per temperature change ((kg/kg)/K)
    """
    ## 1. calculate des/dT 
    dedT = (6.112*17.67*243.5*np.exp(17.67*T_degC/(T_degC+243.5))/(243.5+T_degC)**2)
    ## 2. convert to q
    dqdT = (0.622 * dedT)/(p_hPa - (0.378 * dedT))
    return(dqdT)
    
def dTdqs(p_hPa, q_kgkg):
    """
    INPUTS
        - p_hPa : hPa,          pressure
        - q_kgkg : kg kg-1,     specific humidity
    
    ACTION
        1. CC eq approximation by Bolton (1980) differentiated by T, obtained from
        https://www.wolframalpha.com/input/?i=differentiate+a*exp(b*x%2F(x%2Bc))
        2. convert saturation vapor pressure to saturation specific humidity
        NOTE: pressures (e=water vapor, p=total) must be in hPa
        returns specific humidity in kg/kg
        https://archive.eol.ucar.edu/projects/ceop/dm/documents/refdata_report/eqns.html
        (copied the above from dqsdT function; same thing, just the inverse.)
    
    OUTPUTS
        - dTdq: K kg kg-1,      the change of temperature per specific humidity change (K/(kg/kg)
    """
    ## 1. convert q to e; for this, Bolton approx. q=(0.622 * e)/(p - (0.378 * e)
    ##    was solved for e using WolframAlpha
    ## https://www.wolframalpha.com/input/?i=q+%3D+(0.622+*+e)%2F(p+-+(0.378+*+e))+solve+for+e
    e_hPa = 2.6455*p_hPa*q_kgkg/(q_kgkg+1.6455)
    ## 2. now plug into this differentiated CC eq (solved for T) to find slope
    ### MY VERSION ###
    #dTde = 17.67*243.5/(e_hPa*(np.log(e_hPa/6.112)-17.67)**2)
    ### DIFFERENTIATED BOLTON (1980) EQ. 11 ###
    dTde = (243.5*19.48-440.8)/(e_hPa*(19.48-np.log(e_hPa))**2)
    ## 3. finally, convert back to Kelvin per specific humidity...
    dTdq = dTde*(e_hPa/((0.622 * e_hPa)/(p_hPa - (0.378 * e_hPa))))
    return(dTdq)
    
def calc_theta_e(p_Pa, q_kgkg, T_K):
    """
    INPUTS
        - p_Pa : hPa,           pressure
        - q_kgkg : kg kg-1,     specific humidity
        - T_K : K,              temperature
    
    ACTION
        Dewpoint temperature is calculated to approximate the temperature at
        the lifting condensation level, which is then used to find theta there,
        and at last calculate theta-e.
        The procedure is essentially identical to what is used by MetPy,
        with the exception of making use of eq. 6.5 instead of eq. B39
        (see Davies-Jones, 2009; doi=10.1175/2009MWR2774.1)
    
    OUTPUTS
        - theta-e: K,   equivalent potential temperature
    """
    ## convert specific humidity to mixing ratio (exact)
    r_kgkg = -q_kgkg/(q_kgkg-1) # needs q in kg/kg
    r_gkg  = r_kgkg*1e3
    
    ## calculate theta using Bolton eq. 7
    p_hPa = p_Pa/1e2
    theta = T_K*(1000/p_hPa)**(0.2854*(1-0.00028*r_gkg))
    
    ## convert q into e, vapor pressure (Bolton 1980)
    e_Pa = q_kgkg*p_Pa/(0.622+0.378*q_kgkg)
    
    ## calculate dewpoint temperature according to inverted Bolton 1980 eq.,
    ## https://archive.eol.ucar.edu/projects/ceop/dm/documents/refdata_report/eqns.html
    e_hPa = e_Pa / 1e2
    T_D = np.log(e_hPa/6.112)*243.5/(17.67-np.log(e_hPa/6.112)) + 273.15
    
    ## temperature at lifting condensation level (all temperatures in KElVIN; eq. 15)
    T_L = (1. / (1. / (T_D - 56.) + np.log(T_K / T_D) / 800.)) + 56.
    
    ## theta at lifting condensation level (eq. 3.19, Davies-Jones 2009)
    theta_DL = theta*(((theta/T_L)**(0.28*r_kgkg))*(1+(r_kgkg/EPSILON))**KAPPAD)
    
    ## theta at lifting condensation level (eq. 6.5, Davies-Jones 2009)
    theta_E = theta_DL*np.exp((L0s-L1s*(T_L-C)+K2*r_kgkg)*r_kgkg/(CPD*T_L))
    
    return(theta_E)
    
def calc_theta(p_Pa, q_kgkg, T_K):
    """
    INPUTS
        - p_Pa : hPa,           pressure
        - q_kgkg : kg kg-1,     specific humidity
        - T_K : K,              temperature
    
    ACTION
        The potential temperature is calculated using Bolton's eq. 7;
        to this end, specific humidity is converted into mixing ratio
    
    OUTPUTS
        - theta: K,     potential temperature
    """
    ## convert specific humidity to mixing ratio (exact)
    r_kgkg = -q_kgkg/(q_kgkg-1) # needs q in kg/kg
    r_gkg  = r_kgkg*1e3
    
    ## calculate theta using Bolton eq. 7
    p_hPa = p_Pa/1e2
    return( T_K*(1000/p_hPa)**(0.2854*(1-0.00028*r_gkg)) )
    
def moist_ascender(p_Pa, q_kgkg, T_K):
    """    
    INPUTS
        - p_Pa : hPa,           pressure
        - q_kgkg : kg kg-1,     specific humidity
        - T_K : K,              temperature
    
    ACTION
        Calculate lifting condensation level (LCL) temperature analogously to
        theta_e,
        then approximate lapse rate at LCL using the approximation from:
            http://glossary.ametsoc.org/wiki/Saturation-adiabatic_lapse_rate
    
    OUTPUTS
        - T_LCL: K,     temperature at LCL
        - MALR: K m-1,  moist adiabatic lapse rate
    
    """
    ## convert q into e, vapor pressure (Bolton 1980)
    e_Pa = q_kgkg*p_Pa/(0.622+0.378*q_kgkg)
    
    ## calculate dewpoint temperature according to inverted Bolton 1980 eq.,
    ## https://archive.eol.ucar.edu/projects/ceop/dm/documents/refdata_report/eqns.html
    e_hPa = e_Pa / 1e2
    T_D = np.log(e_hPa/6.112)*243.5/(17.67-np.log(e_hPa/6.112)) + 273.15
    
    ## temperature at lifting condensation level (all temperatures in KElVIN; eq. 15)
    T_LCL = (1. / (1. / (T_D - 56.) + np.log(T_K / T_D) / 800.)) + 56.
    
    ## now that we have T_LCL, calculate moist adiabatic lapse rate...
    
    ## convert specific humidity to mixing ratio (exact)
    r_kgkg = -q_kgkg/(q_kgkg-1) # needs q in kg/kg
    ## now, finally, approximate lapse rate
    MALR = GRAV*( (RSPECIFIC*(T_LCL**2)) + (HV*r_kgkg*T_LCL) )/( (CPD*RSPECIFIC*T_LCL**2) + (HV**2*r_kgkg*EPS))
    
    return(T_LCL, -MALR) # T_LCL in Kelvin, MALR in K/m

