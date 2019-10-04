#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  2 21:52:42 2019

@author: dominik

===============================================================================

copied from cluster, 19-09-2019

===============================================================================

"""

###########################################################################
#############################    MODULES ##################################

import gzip
import pandas as pd
import numpy as np
import os, fnmatch
import timeit
import netCDF4 as nc4
import sys
import random
import multiprocessing    
from joblib import Parallel, delayed
from copy import deepcopy
from datetime import datetime, timedelta
from math import sin,cos,acos,atan,atan2,sqrt

###############################################################################
###############################################################################

def absoluteFilePaths(directory, flist):
    ## returns absolute file paths for given file list
    paths = list()
    for dirpath,_,filenames in os.walk(directory):
        for f in filenames:
            if f in flist: 
                paths.append(os.path.abspath(os.path.join(dirpath, f)))
    return(paths)

def dist_on_sphere(lat1,lon1,lat2,lon2):
    """
    This function calculates the linear distance between two points
    defined by lat1,lon1 and lat2,lon2 (in °, lat [-90,90], lon[-180,180])
    on Earth's surface, relying on the following assumptions:

        - Earth is a perfect sphere,
        - the radius is precisely 6371 km

    The calculation is based on a transformation to cartesian coordinates,
    followed by computing the dot product between two vectors spanned
    by the two points on Earth's surface and the center, respectively.

    ============ required functions ================

    all from module math:
        - sin
        - cos
        - acos
        
    RETURNS:
        
        distance between locations in kilometers.

    """
    

    ## set Earth radius, define pi
    R  = 6371
    pi = 3.141592654

    ## first off, must adapt lat/lon for spherical coordinates
    lat1 = 90 - lat1
    lat2 = 90 - lat2
    lon1 = 180 + lon1
    lon2 = 180 + lon2

    ## second, must obtain angles in radian from (arc) degree
    la1 = (pi/180)*lat1
    la2 = (pi/180)*lat2
    lo1 = (pi/180)*lon1
    lo2 = (pi/180)*lon2

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

    #len1  = sqrt(x1**2 + y1**2 + z1**2)
    #len2  = sqrt(x2**2 + y2**2 + z2**2)

    angle = acos(dotp)

    ## fifth, calculate distance
    dist = angle * R   # pi cancels out

    return(dist)

def gridded_area_exact(lats_centr, res, R):
    """
    ============================== INFO =======================================
    INPUT
        lats_centr:    lats [deg North], between -90 to +90 (float or 1D np ary)
         res:        (regular!) grid resolution [deg]
         R:          Earth radius [km]
    OUTPUT
        area:        EXACT gridded area, shape as defined by input [km^2]
        ACTION
        based on the grid spacing inferred by lats & lons,
        areas of all grid cells specified by indices and coordinates (lats)
        are computed and summed up in the end.
        Since the spacing between grid cell border longitudes ALWAYS remains
        the same, that is, the resolution in degrees, this function does NOT
        rely on any longitudinal input.
    ***************************************************************************
        CAUTION
        1.) this will NOT work properly close to the poles!
        2.) this is based on input coordinates referring to CENTROIDS
        3.) only for regular grids (independent from longitude)
        4.) will produce CRAP if / results in integer division (Python 2.7) !
    ***************************************************************************
    NOTES
        - replaced ugly loop by vectorized expressions,
            speed-up based on 10000 iterations: x50
        - original formula, using degrees (caution -- np.sin requires radian)
            A = (pi/180)R^2 |sin(lat1)-sin(lat2)| |lon1-lon2|
          obtained from:
            pmel.noaa.gov/maillists/tmap/ferret_users/fu_2004/msg00023.html

    ========================== DSc, January 2018 ==============================
    """

    ## make use of numpy vectorization
    lats1 = (lats_centr+(res/2))*np.pi/180 # np.sin requires radians
    lats2 = (lats_centr-(res/2))*np.pi/180
    areas = (np.pi/180)*(R**2)*np.abs(np.sin(lats1)-np.sin(lats2))*res

    ## overwrite any areas of 0 (at the poles) with np.NaN to prevent problems
    try:
        areas[np.where(areas==0.)] = np.NaN # only works for arrays
    except TypeError:
        pass # simply ignore if it's a float
    return(areas)

def midpoint_on_sphere(lat1,lon1,lat2,lon2):

    """
    
    INPUT: must come as lons [-180 .. 180], lats [-90 .. 90]
    
    WARNING: the coordinate transformation is probably correct,
             but the last transformations were done empirically
             (and not based on logic). 
             A number of checks has not revealed any issues.

    DEPENDENCIES:

        numpy, various trigonometric functions from module math
        
    OUTPUT: geographical center coordinates (lat/lon)

    """

    ## define pi
    pi = 3.141592654

    ## first off, must adapt lat/lon for spherical coordinates
    lat1 = 90 - lat1
    lat2 = 90 - lat2
    lon1 = 180 + lon1
    lon2 = 180 + lon2

    ## second, must obtain angles in radian from (arc) degree
    lat1 = (pi/180)*lat1
    lat2 = (pi/180)*lat2
    lon1 = (pi/180)*lon1
    lon2 = (pi/180)*lon2

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
        lat_mid = 90 - theta*(180/pi)
    else:
        lat_mid = -(90 + theta*(180/pi))
    if phi > 0:
        lon_mid = (phi*(180/pi)-180)
    else:
        lon_mid = 180 + phi*(180/pi)

    return(lat_mid, lon_mid)

def q2rh(q,p,T):
    """
    INPUTS
    -------
        q: kg/kg,    float or vector
        p: Pa,       float or vector
        T: K,        float or vector
    
    ACTION
    ------
        Converting specific to relative humidity following Bolton (1980), see:
        https://www.eol.ucar.edu/projects/ceop/dm/documents/refdata_report/eqns.html

    OUTPUTS
    -------
        returns RH in percent.
    """
    # convert q into e, vapor pressure
    e = q*p/(0.622+0.378*q)

    # compute saturation vapor pressure, must convert T to °C
    Tref = 273.15
    es = 611.2*np.exp(17.67*(T-Tref)/(T-Tref+243.5))

    # return relative humidity
    return(1e2*e/es)

def dqsdT(p_hPa, T_degC):
    """
    INPUTS
    -------
        - pressure [>>>hPa<<<], 
        - temperature [degree Celsius]
    
    ACTION
    -------
        1. CC eq approximation by Bolton (1980) differentiated by T, obtained from
        https://www.wolframalpha.com/input/?i=differentiate+a*exp(b*x%2F(x%2Bc))
        
        2. convert saturation vapor pressure to saturation specific humidity
        NOTE: pressures (e=water vapor, p=total) must be in hPa
        returns specific humidity in kg/kg
        https://archive.eol.ucar.edu/projects/ceop/dm/documents/refdata_report/eqns.html
    
    OUTPUTS
    -------
        - dqdT: the change of specific humidity per temperature change ((kg/kg)/K)
    """
    ## 1. calculate des/dT 
    dedT = (6.112*17.67*243.5*np.exp(17.67*T_degC/(T_degC+243.5))/(243.5+T_degC)**2)
    ## 2. convert to q
    dqdT = (0.622 * dedT)/(p_hPa - (0.378 * dedT))
    return(dqdT)
    
def dTdqs(p_hPa, qv_kgkg):
    """
    INPUTS
    -------
        - pressure [>>>hPa<<<], 
        - specific humidity [kg/kg]
    
    ACTION
    -------
        1. CC eq approximation by Bolton (1980) differentiated by T, obtained from
        https://www.wolframalpha.com/input/?i=differentiate+a*exp(b*x%2F(x%2Bc))
        
        2. convert saturation vapor pressure to saturation specific humidity
        NOTE: pressures (e=water vapor, p=total) must be in hPa
        returns specific humidity in kg/kg
        https://archive.eol.ucar.edu/projects/ceop/dm/documents/refdata_report/eqns.html
        
        (copied the above from dqsdT function; same thing, just the inverse.)
    
    OUTPUTS
    -------
        - dTdq: the change of temperature per specific humidity change (K/(kg/kg)
    """
    ## 1. convert q to e; for this, Bolton approx. q=(0.622 * e)/(p - (0.378 * e)
    ##    was solved for e using WolframAlpha
    ## https://www.wolframalpha.com/input/?i=q+%3D+(0.622+*+e)%2F(p+-+(0.378+*+e))+solve+for+e
    e_hPa = 2.6455*p_hPa*qv_kgkg/(qv_kgkg+1.6455)
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
    -------
        - pressure [Pa], 
        - specific humidity [kg/kg],
        - temperature [K]
    
    ACTION
    -------
        Dewpoint temperature is calculated to approximate the temperature at
        the lifting condensation level, which is then used to find theta there,
        and at last calculate theta-e.
        The procedure is essentially identical to what is used by MetPy,
        with the exception of making use of eq. 6.5 instead of eq. B39
        (see Davies-Jones, 2009; doi=10.1175/2009MWR2774.1)
    
    OUTPUTS
    -------
        - theta-e: equivalent potential temperature [K]
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
    epsilon = 0.6620
    kappa_d = 0.2854
    theta_DL = theta*(((theta/T_L)**(0.28*r_kgkg))*(1+(r_kgkg/epsilon))**kappa_d)
    
    ## theta at lifting condensation level (eq. 6.5, Davies-Jones 2009)
    C  = 273.15
    L0s = 2.56313e6 
    L1s = 1754
    K2 = 1.137e6
    cpd = 1005.7
    theta_E = theta_DL*np.exp((L0s-L1s*(T_L-C)+K2*r_kgkg)*r_kgkg/(cpd*T_L))
    
    return(theta_E)
    
def calc_theta(p_Pa, q_kgkg, T_K):
    """
    INPUTS
    -------
        - pressure [Pa], 
        - specific humidity [kg/kg],
        - temperature [K]
    
    ACTION
    -------
        The potential temperature is calculated using Bolton's eq. 7;
        to this end, specific humidity is converted into mixing ratio
    
    OUTPUTS
    -------
        - theta: potential temperature [K]
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
    -------
        - pressure [Pa], 
        - specific humidity [kg/kg],
        - temperature [K]
    
    ACTION
    -------
        Calculate lifting condensation level (LCL) temperature analogously to
        theta_e,
        then approximate lapse rate at LCL using the approximation from:
            http://glossary.ametsoc.org/wiki/Saturation-adiabatic_lapse_rate
    
    OUTPUTS
    -------
        - T_LCL: temperature at LCL (K)
        - MALR: moist adiabatic lapse rate (K/m)
    
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
    g   = 9.8076
    Hv  = 2.5e6
    Rsd = 287.057 # round(8.3144598 / (28.9645/1e3), 3)
    Rsw = 461.5   # specific gas constant of water vapor 
    eps = Rsd/Rsw
    cpd = 1005.7
    
    ## convert specific humidity to mixing ratio (exact)
    r_kgkg = -q_kgkg/(q_kgkg-1) # needs q in kg/kg
    ## now, finally, approximate lapse rate
    MALR = g*( (Rsd*(T_LCL**2)) + (Hv*r_kgkg*T_LCL) )/( (cpd*Rsd*T_LCL**2) + (Hv**2*r_kgkg*eps))
    
    return(T_LCL, -MALR) # T_LCL in Kelvin, MALR in K/m

def diagnoser(ID, npart, array, 
              dTH_thresh=1.0, f_dqsdT=1.0, f_dTdqs=1.0, sample_E_upto=0, sample_H_upto=0, 
              P_dq_min=None, P_dT_thresh=0, P_RHmin=80):
    
    """
    NO description yet.. :/
    """
        
    ###### CONSTANTS #######
    R_specific = 287.057    # round(8.3144598 / (28.9645/1e3), 3)
    DALR       = (-9.8/1e3) # dry adiabatic lapse rate

    ## determine indices, copy data to small nparticle loop arrays
    ii1 = ID +   npart
    ii2 = ID + 2*npart
    ii3 = ID + 3*npart
    #ii4 = ID + 4*npart
    ii5 = ID + 5*npart
    ii6 = ID + 6*npart
    ii7 = ID + 7*npart
    ii8 = ID + 8*npart

    #######################################################################
    #pID   = array[ ID,:2]  # ID
    lats  = array[ii2,:2]  # LAT (>> NOT << lon )
    lons  = np.copy(array[ii1,:2])  # LON; making copy prevents error msg
    ######################## for ERA-global ###############################
    ## as described in TECHNICAL NOTE
    if ID < 0.002*npart: # only check first few files
        if dist_on_sphere(lats[0],lons[0],lats[1],lons[1])>1620:
            return(0)
    #######################################################################
    temp  = array[ii8,:2]  # temperature (K)
    ztra  = array[ii3,:2]  # height (m)
    #topo  = array[ii4,:2]  # topography (m) 
    qv    = array[ii5,:2]  # specific humidity (kg/kg)
    hmix  = array[ii7,:2]  # ABL height (m)
    dens  = array[ii6,:2]  # density (kg/m^3)
    #######################################################################

    ## calculate everything else that is needed at least once from here
    dq = qv[0] - qv[1]
    BLh_max = np.max(hmix[:2]) 

    #######################################################################
    #######################################################################
#        import metpy.calc as mpcalc
#        from metpy.units import units
    pres = dens*R_specific*temp #* units('Pa')
#        td = mpcalc.dewpoint_from_specific_humidity(qv * units('kg/kg'), temp * units('K'), pressure)
#        potT = mpcalc.potential_temperature(pres, temp * units('K'))
#        equivpotT = mpcalc.equivalent_potential_temperature(pres, temp * units('K'), td)
    potT = calc_theta(pres, qv, temp)  
    equivpotT = calc_theta_e(pres, qv, temp)
    dTH  = potT[0] - potT[1]
    dTHe = (equivpotT[0]-equivpotT[1])
    #######################################################################
    ####################################################################### 
    
    ####################
    counter = 0 # not elegant, but it works
    ####################
    
    """
    dz = (ztra[0]+topo[0]) - (ztra[1]+topo[1]) # WHY WOULD THIS BE CORRECT?
    """
    dz = ztra[0] - ztra[1]
    
    ########### PRECIPITATION ############
    if ( 
          dq < P_dq_min and
          q2rh(qv[0], dens[0]*R_specific*temp[0], temp[0])>P_RHmin and 
          q2rh(qv[1], dens[1]*R_specific*temp[1], temp[1])>P_RHmin 
          ):
            counter += 1
    elif ( 
        dq < P_dq_min and
        dz > 0 # not really a parameter, but rather: is the parcel ascending?
        ): 
        ## moist adiabatic lapse rate, T @ Lifting Condensation Level
        T_LCL, MALR = moist_ascender(p_Pa=pres[1], q_kgkg=qv[1], T_K=temp[1])
        ## calculate whether LCL was reached during ascent
        dz_reachLCL = (T_LCL-temp[1])/DALR
        if dz > dz_reachLCL:
            dz_rem = dz - dz_reachLCL    # remaining ascent
            T_moi = T_LCL  + dz_rem*MALR # hypoth. T if partially moist adiabatic
            T_dry = temp[1]+ dz*DALR     # hypoth. T if purely dry adibatic
            ## check if temperatures are somewhat consistent
            if ( # yes, yet another if
                # check if temperature change closer to what we expect in case
                # of condensation occurring during ascent
                abs(temp[0]-T_moi) < abs(temp[0]-T_dry) and 
                (temp[0]-T_moi) < P_dT_thresh 
                # above: if air is 'too warm', ABL 'detrainment' in large-scale
                # subsidence areas more likely (results in pseudo-precip)
                ): 
                counter += 1

    ########### EVAPORATION & SENSIBLE HEAT ############
    if (
        (ztra[0] < hmix[0]) and
        (ztra[1] < hmix[1]) and
        (dTH > dTH_thresh) and 
        ((dTHe - dTH) > dTH_thresh) and         
        abs(dq) < f_dqsdT*(dTH)*dqsdT(p_hPa=dens[1]*R_specific*temp[1]/1e2, T_degC=temp[1]-273.15) and
        abs(dTH) < f_dTdqs*(dq)*dTdqs(p_hPa=dens[1]*R_specific*temp[1]/1e2, qv_kgkg=qv[1])
        ):
        ### EVAP-HEAT ###
        print(".")
        counter += 6
        return(counter)
        
    ########### SENSIBLE HEAT ############
    if (   
        (ztra[0] <  max(sample_H_upto, BLh_max)) and 
        (ztra[1] <  max(sample_H_upto, BLh_max)) and 
        (dTH > dTH_thresh) and 
        abs(dq) < f_dqsdT*(dTH)*dqsdT(p_hPa=dens[1]*R_specific*temp[1]/1e2, T_degC=temp[1]-273.15)
        ):
        ### HEAT ###
        counter += 4
    
    ########### EVAPORATION ############
    if ( 
        (ztra[0] <  max(sample_E_upto, BLh_max)) and 
        (ztra[1] <  max(sample_E_upto, BLh_max)) and
        ((dTHe - dTH) > dTH_thresh) and
         abs(dTH) < f_dTdqs*(dq)*dTdqs(p_hPa=dens[1]*R_specific*temp[1]/1e2, qv_kgkg=qv[1])
       ):
        ### EVAP ###
        counter += 2

    return(counter)
    
def gridder(ID, npart, array, code,
            fileID, glat, glon, gd_area, 
            logger, heat, evap, prec):
    
    """
    NO description yet.. :/
    """

    #if ID%1e5==0: print("..")  
    
    ###### CONSTANTS #######
    R_specific = 287.057 # round(8.3144598 / (28.9645/1e3), 3)
    cpd        = 1005.7
    mass_part  = 2548740090557.712 # 5.09198e18/1997842
    
    ## determine indices, copy data to small npart loop arrays
    ii1 = ID +   npart
    ii2 = ID + 2*npart

    #######################################################################
    lats  = array[ii2,:2]  # LAT (>> NOT << lon )
    lons  = array[ii1,:2]  # LON

    ## enforce proper coordinates for midpoint function [-180.0 .. 180.0]
    lons[lons>180.0] -= 360
    
    ## use own function to obtain midpoint position  
    lat_interm,lon_interm = midpoint_on_sphere(lats[0],lons[0],lats[1],lons[1])
            
    ## now shift all coords that otherwise would be allocated to +180 deg to - 180
    if (lon_interm>179.5): lon_interm -= 360
    
    ## finally, obtain intermediate position indices
    lat_ix = np.argmin(np.abs(glat-lat_interm))
    lon_ix = np.argmin(np.abs(glon-lon_interm))
    
    ###################################################################
    logger[fileID,lat_ix,lon_ix] += 1 ## store npart per pixel
    ###################################################################
    
    if code > 0:
        ii5 = ID + 5*npart
        qv    = array[ii5,:2]  # specific humidity (kg/kg)
        dq = qv[0] - qv[1]
    
        if code>=4:
            #ii3 = ID + 3*npart
            #ii4 = ID + 4*npart
            ii6 = ID + 6*npart
            ii8 = ID + 8*npart
            temp  = array[ii8,:2]  # temperature (K)
            #ztra  = array[ii3,:2]  # height (m)
            #topo  = array[ii4,:2]  # topography (m) 
            dens  = array[ii6,:2]  # density (kg/m^3)
            pres = dens*R_specific*temp # Pa
            theta = calc_theta(pres, qv, temp)  
            dT = (theta[0] - theta[1])*cpd # <<-------------------------------- assuming dry air; check this!
        #######################################################################
    
        if code >= 4:
            heat[fileID,lat_ix,lon_ix] += mass_part*dT/(1e6*gd_area[lat_ix]*6*3600) # Wm-2
        if code in [2,3,6,7]: # I know, so amateur-like. sorry.    
            evap[fileID,lat_ix,lon_ix] += mass_part*dq/(1e6*gd_area[lat_ix]) # mm
        if code%2-1==0:
            prec[fileID,lat_ix,lon_ix] += mass_part*dq/(1e6*gd_area[lat_ix]) # mm
            
    #-- DONE!
    return(logger, heat, evap, prec)


############################################################################
#############################    SETTINGS ##################################

def readNmore(
           runyr, ayear, amonth,
           inpath="/media/gvo00090_scratch/vsc42561/tools/particle-o-matic/era_global/terabox/",
           outpath="/tmp/", sfnam_base="FXvG_r",           
           dTH_thresh=1.0, # used for E,H,P (if P_dq_min==None)
           f_dqsdT=0.7, f_dTdqs=0.7, # for H, E diagnosis (lower = more strict)
           sample_E_upto=0, sample_H_upto=0, # set min ABLh, disabled if 0 
           P_dq_min=None, P_dT_thresh=0, P_RHmin=80): # P settings

    """
    comments
    
    - more sophisticated ABL criteria as in versionD are highly beneficial to
      avoid E-overestimation over deserts; for H, it does not seem to change
      as much.
    - with the current configuration, there are only 4 parameters:
        
        dTH_thresh = 1. (Kelvin),
        f_dqdst == f_dTdqs,
        P_dT_thresh = 0. (Kelvin), # not a good idea to increase this a lot    
        P_RHmin=80 (%) 
        
        thus, the previosuly introduced dz-Parameter could return,
        with the advantage of being used both for E,H & P 
        (as of now, dz>0 is used for P anyways).
    """

    ###########################################################################    
    save_output = True   #<=================
    log_this    = False
    verbose     = True
    timethis    = True
    ###########################################################################
    
    ## construct precise input and storage paths
    basepath  = inpath + "/NH/year"
    basepath2 = inpath + "/SH/year"    
    runyr     = str(runyr)
    mainpath  = basepath+runyr+"/output_temp/" # remove output_temp for data on DATA (not SCRATCH)
    mainpath2 = basepath2+runyr+"/output_temp/"
    sfilename = sfnam_base+runyr[-2:]+"_"+str(ayear)+"-"+str(amonth)+".nc"
    
    ########### LOG W/IN PYTHON SCRIPT by redirecting output #############
    if log_this:    
        new_target = open(outpath+sfilename[:-3]+'.txt', 'w')
        old_target, sys.stdout = sys.stdout, new_target
    
    if verbose:
        print("\n====================================================================")
        print("=========================--- WELCOME ---============================")
        print("\n\n----- this is how the file will be stored:\n", sfilename)
        print("\n====================================================================\n")
        
    ##########################    EXPERIMENTAL    #############################
    if P_dq_min == None:
        if verbose:
            print("\n--- INFO: P_dq_min is calculated based on d(theta)-threshold!")
        dummy_dq = 0.2 # this choice doesn't matter too much...
        P_dq_min = -(1/(calc_theta_e(1013.25e2, (5+dummy_dq)/1e3, 273.15+15) - 
                       calc_theta_e(1013.25e2, 5/1e3, 273.15+15)))*dummy_dq/1e3
        print("P_dq_min = ", 1e3*P_dq_min, "g/kg")
    elif P_dq_min > 0:
        raise SystemExit("------ FATAL ERROR: P_dq_min should be negative (and in kg/kg)!")
    ###########################################################################
        
    ###########################################################################
    #############################     SETUP      ##############################
    
    ## start timer
    if timethis:
        megatic = timeit.default_timer()
    
    ## prepare grid
    resolution = 1. # in degrees
    glat = np.arange(-90,90+resolution,resolution) # lats from -90 to + 90 deg
    glon = np.arange(-180,180,resolution)
    nlat = glat.size
    nlon = glon.size
    
    ###########################################################################
    
    ## prepare filelist
    filelist = fnmatch.filter(os.listdir(mainpath), '*.dat.gz')
    filelist.sort() # probably unnecessary
       
    ## translate into datetime 
    alldates = []
    for file in filelist:
        filedate = file[-17:-7]
        alldates.append(datetime(int(filedate[:4]),  int(filedate[4:6]),
                      int(filedate[6:8]), int(filedate[8:10]))) # check if files change!
    #fdate_end = filelist[-1][-17:-7] # not needed
    
    ## create start date to find right files
    if amonth==11 and int(ayear)!=int(runyr):
        ## exception for incomplete November (could go a bit further; OCT incomplete)
        date_bgn = datetime(year=ayear, month=amonth, day=20, hour=6)
    else:
        date_bgn = datetime(year=ayear, month=amonth, day=1, hour=6)
    
    ## same thing for end date; for DECEMBER, 'trick' doesn't work
    if amonth!=12: # these are easy
        date_end = datetime(year=ayear, month=amonth+1, day=1, hour=0)
    elif amonth==12 and int(ayear)==int(runyr):
        date_end = datetime(year=ayear, month=amonth, day=31, hour=18)
    else:
        date_end = datetime(year=ayear+1, month=1, day=1, hour=0)
    
    fID_bgn = np.where(np.asarray(alldates)==date_bgn)[0][0]
    fID_end = np.where(np.asarray(alldates)==date_end)[0][0] + 1
    
    ## crop filelist accordingly
    filelist_orig = deepcopy(filelist)
    filelist = filelist[fID_bgn:fID_end]
    
    dt_dates = []
    for somedate in alldates[fID_bgn:fID_end]: 
        dt_dates.append(somedate - timedelta(hours=3)) 
    dt_dates = np.asarray(dt_dates)
    
    if verbose:
        print("\nfirst file --- ", filelist[0])
        print("last file  --- ", filelist[-1])
        
        print("\nfirst uptake time --- ", dt_dates[0])
        print("last uptake time  --- ", dt_dates[-1])
    
    ## complete file paths
    filepaths = absoluteFilePaths(mainpath, filelist)
    filepaths.sort()
    
    ###########################################################################
    ######################   REPEAT FOR SOUTHERN HEMISPHERE    ################
       
    ## prepare filelist
    filelist2 = fnmatch.filter(os.listdir(mainpath2), '*.dat.gz')
    filelist2.sort() # probably unnecessary
    
    if len(filelist2) != len(filelist_orig):
        raise SystemExit("---- ERROR: no can do")
    
    filelist2 = filelist2[fID_bgn:fID_end]
        
    if verbose:
        print("\n---- now for the southern hemisphere:")
        print("\nfirst file --- ", filelist2[0])
        print("last file  --- ", filelist2[-1])
    
    ## complete file paths
    filepaths2 = absoluteFilePaths(mainpath2, filelist2)
    filepaths2.sort()
    
    
    ###########################################################################
    ##########    unzip first file, read first two lines to obtain dims  ######
    
    firstfile = filelist_orig[0]
    os.chdir(mainpath)
    rint = random.randint(10000,100000)
    os.system(str("cp "+firstfile+" "+firstfile[:-7]+"_TMP-"+str(rint)+".DAT.gz"))
    os.system(str("gunzip "+firstfile[:-7]+"_TMP-"+str(rint)+".DAT.gz"))
    
    with open(mainpath+firstfile[:-7]+"_TMP-"+str(rint)+".DAT", 'r') as f:
        frst_line = f.readline()
        scnd_line = f.readline()
    
    if verbose:
        print("\nworking in the following dir: ", mainpath)
        print("name of first file (used for grabbing dimensions):", firstfile)
        print("\n1st line of first file:", frst_line[:-1])
        print("2nd line of first file:", scnd_line[:-1])
    
    ## clean up
    os.system(str("rm "+firstfile[:-7]+"_TMP-"+str(rint)+".DAT"))
    
    ## now extract dimensions.......
    if len(scnd_line)==12: # added for SH output: more than 1 mio particles!
        nparticle = int(scnd_line[:7])
        ntrajstep = int(scnd_line[8])
        ncolumn   = int(scnd_line[10])
    else:
        nparticle = int(scnd_line[:6])
        ntrajstep = int(scnd_line[7])
        ncolumn   = int(scnd_line[9])
    
    if verbose:
        print("\nDimensions determined as follows:")
        print("----------------------------------------------------------------------")
        print("nparticle = ",nparticle, " |  ntrajstep=",ntrajstep,"  | ncolumn=",ncolumn)
        print("----------------------------------------------------------------------")
        print("\n--- NOTE: nparticle varies between timesteps and hence not needed")
        print("\n--- INFO: trajectory is evaluated only for steps now and -6h")
    
    ###########################################################################
    ###########################################################################
    
    if verbose:
        print("\n===========================  CONFIGURATION  =======================")
        print("\nayear=", ayear)
        print("\nresolution=", resolution)
        print("nlat=", nlat)
        print("nlon=", nlon)
    
        print("\nncolumn=", ncolumn)
        print("nday=", len(filelist)/4)
        print("ntrajstep=", ntrajstep)
        print("\n===================================================================")
    
    ###########################################################################
    ###########################################################################
    
    ## grab gridded areas prior to loop to save some CPU time
    gd_area = gridded_area_exact(glat, res=1.0, R=6371)
    
    ## pre-allocate arrays
    ary_heat     = np.zeros(shape=(len(filelist),nlat,nlon))
    ary_evap     = np.zeros(shape=(len(filelist),nlat,nlon))
    ary_prec     = np.zeros(shape=(len(filelist),nlat,nlon))
    npart_log    = np.zeros(shape=(len(filelist),nlat,nlon), dtype=int)
    
    if verbose:
        print("\n======================================================================")
        print("=======    main loop initiated! ")
    
    ###########################################################################
    ##########################      MAIN -- LOOP ##############################
    for ix in range(len(filelist)):
    
        if verbose:
            print(ix, ".. ", end="")
    
        ## grab file and read it
        file = filelist[ix]
        
        ary_pd = pd.read_table(gzip.open(filepaths[ix], 'rb'), sep="\s+", header=None, skiprows=2) #, header=1)
        ## now can proceed with  >> FLIPPING << data (reverse time axis)
        ary = np.fliplr(np.asarray(ary_pd)) # flips axis 1
    
        ## must ALWAYS determine nparticle, as it changes from file to file
        if ary.shape[0]%ncolumn == 0:
            nparticle = int(ary.shape[0] / ncolumn)
        else:
            raise ValueError("Array shape incompatible with code")
            
        ## here comes some multiprocessing to diagnose all particles
        if __name__ == '__main__':
            num_cores = multiprocessing.cpu_count()
            diagcodes = Parallel(n_jobs=num_cores)(delayed(diagnoser)(
                                 ID=i,npart=nparticle,array=ary,
                                 dTH_thresh=dTH_thresh, f_dqsdT=f_dqsdT, f_dTdqs=f_dTdqs, 
                                 sample_E_upto=sample_E_upto, sample_H_upto=sample_H_upto, 
                                 P_dq_min=P_dq_min, P_dT_thresh=P_dT_thresh, P_RHmin=P_RHmin
                                 ) for i in range(nparticle))

        ## proceed to call gridding function w/o parallelization
        for i in range(nparticle):
            npart_log, ary_heat, ary_evap, ary_prec = gridder(
                     ID=i, npart=nparticle, array=ary, code=diagcodes[i],
                     fileID=ix, glat=glat, glon=glon, gd_area=gd_area, 
                     logger=npart_log, heat=ary_heat, evap=ary_evap, prec=ary_prec)
                
                 
        #######################################################################
        #######################################################################
        #######################################################################
        
        #print("---- repeating procedure for SH..")
        
        ## nomenclature here for SH could easily win an award for worst coding practice
        ary2_pd = pd.read_table(gzip.open(filepaths2[ix], 'rb'), sep="\s+", header=None, skiprows=2) #, header=1)
        ary2 = np.fliplr(np.asarray(ary2_pd)) # flips axis 1
    
        if ary2.shape[0]%ncolumn == 0:
            nparticle2 = int(ary2.shape[0] / ncolumn)
        else:
            raise ValueError("Array shape incompatible with code")
            
        if __name__ == '__main__':
            #num_cores = multiprocessing.cpu_count()
            diagcodes2 = Parallel(n_jobs=num_cores)(delayed(diagnoser)(
                                  ID=i,npart=nparticle2,array=ary2,
                                  dTH_thresh=dTH_thresh, f_dqsdT=f_dqsdT, f_dTdqs=f_dTdqs, 
                                  sample_E_upto=sample_E_upto, sample_H_upto=sample_H_upto, 
                                  P_dq_min=P_dq_min, P_dT_thresh=P_dT_thresh, P_RHmin=P_RHmin
                                  ) for i in range(nparticle2))
            
        ## proceed to call gridding function w/o parallelization
        for i in range(nparticle2):
            npart_log, ary_heat, ary_evap, ary_prec = gridder(
                     ID=i, npart=nparticle2, array=ary2, code=diagcodes2[i],
                     fileID=ix, glat=glat, glon=glon, gd_area=gd_area, 
                     logger=npart_log, heat=ary_heat, evap=ary_evap, prec=ary_prec)
            
        #--- end nparticle loop ----#
    
    #--- end trajstep loop ----#
    if timethis:
        megatoc = timeit.default_timer()
        print("\n=======    main loop completed, total runtime so far: ",str(round(megatoc-megatic, 2)),"seconds")
        print("======================================================================")
        
    
    ###########################################################################    
    
    if save_output:
            
        ### delete nc file if it is present (avoiding error message)
        try:
            os.remove(outpath+sfilename)
        except OSError:
            pass
        
        ### create netCDF4 instance
        nc_f = nc4.Dataset(outpath+sfilename,'w', format='NETCDF4') #'w' stands for write
        
        ### create dimensions ###
        nc_f.createDimension('time', dt_dates.size)
        nc_f.createDimension('lat', nlat)
        nc_f.createDimension('lon', nlon)
    
        ### create variables ###
        times      = nc_f.createVariable('time', 'i4', 'time')
        latitudes  = nc_f.createVariable('lat', 'f4', 'lat')
        longitudes = nc_f.createVariable('lon', 'f4', 'lon')
        
        heats     = nc_f.createVariable('H', 'f4', ('time','lat','lon'))
        evaps     = nc_f.createVariable('E', 'f4', ('time','lat','lon'))
        precs     = nc_f.createVariable('P', 'f4', ('time','lat','lon'))
        
        nparts     = nc_f.createVariable('n_part', 'f4', ('time','lat','lon'))
    
        ### set attributes ###
        nc_f.description = "FLEXPART-based upward land surface fluxes & precipitation"
        times.units       = 'hours since 1900-01-01 00:00:00'
        times.calendar    = 'Standard' # do NOT use gregorian here!
        latitudes.units   = 'degrees_north'
        longitudes.units  = 'degrees_east'
    
        heats.units     = 'W/m2'
        heats.long_name	= 'surface sensible heat flux'
        evaps.units     = 'mm'
        evaps.long_name	= 'evaporation'
        precs.units     = 'mm'
        precs.long_name	= 'precipitation'
    	    
        nparts.units     = 'int'
        nparts.long_name = 'number of parcels (intermed. pos.)'
    
        ### insert data ###
        times[:]       = nc4.date2num(dt_dates, times.units, times.calendar)
        longitudes[:]  = glon
        latitudes[:]   = glat
        
        heats[:] = ary_heat[:]
        evaps[:] = ary_evap[:]
        precs[:] = ary_prec[:]     
    
        nparts[:] = npart_log[:]
        
        ### close file, done ###
        nc_f.close()
        
        print("\n===============================================================")
        print("======  file",sfilename," written to:")
        print("======  ",outpath)
        print("===============================================================")
        
        
    if log_this:
        ## stop writing to text file
        new_target.close()
        sys.stdout = old_target
        print("\n\n\n==================------------------ ALL DONE !!!")
