#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Collection of functions for HAMSTER
# 
# This file is part of HAMSTER, 
# originally created by Dominik Schumacher, Jessica Keune, Diego G. Miralles
# at the Hydro-Climate Extremes Lab, Department of Environment, Ghent University
# 
# https://github.com/h-cel/hamster
# 
# HAMSTER is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation v3.
#
# HAMSTER is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with HAMSTER. If not, see <http://www.gnu.org/licenses/>.
#

import argparse
import calendar
import csv
import fnmatch
import gzip
import imp
import math
import os
import random
import re
import struct
import sys
import time
import timeit
import warnings
from datetime import date, datetime, timedelta
import datetime as datetime
from functools import reduce
from math import acos, atan, atan2, cos, floor, sin, sqrt

import h5py
import netCDF4 as nc4
import numpy as np
import pandas as pd
from dateutil.relativedelta import relativedelta

from hamsterfunctions import *

def get_currentversion():
    version_file = os.path.join("../.","VERSION")
    with open(version_file) as vfile:
        version = vfile.readlines()[0].strip()
    return(version)

def disclaimer():
    print(
        "\n============================================================================================================"
    )
    print(
        "\n============================================================================================================"
    )
    print(os.system("figlet hamster"))
    print(" * Copyright 2021                                                      *")
    print(" * Dominik Schumacher, Jessica Keune                                   *")
    print(" *                                                                     *")
    print(" * This program is part of HAMSTER v"+str(get_currentversion())+"                              *")
    print(" *                                                                     *")
    print(" * HAMSTER is free under the terms of the GNU General Public license   *")
    print(" * version 3 as published by the Free Software Foundation              *")
    print(" * HAMSTER is distributed WITHOUT ANY WARRANTY! (see LICENSE)          *")
    print(
        "\n============================================================================================================"
    )

def check_paths(pfile, path):
    try: 
        fpath = getattr(pfile, path)
    except:
        fpath = ""
    return fpath

###------------------------------------------------------------------------------
###------------------------------------------------------------------------------
## CONSTANTS
###------------------------------------------------------------------------------
###------------------------------------------------------------------------------

RSPECIFIC = 287.057  # specific gas constant
RSPECIFICW = 461.5  # specific gas constant for water vapor
EPS = RSPECIFIC / RSPECIFICW
LAPSERATE_DRY = -9.8 / 1e3  # dry adiabatic lapse rate
EARTHRADIUS = 6371  # earth radius
PI = 3.141592654  # pi
GRAV = 9.8076  #
CPD = 1005.7  # specific heat of dry air [J kg-1 K-1]
HV = 2.5e6
TREF = 273.15  # reference temperature [K]
PREF = 1013.25e2  # reference pressure of standard atmosphere [Pa]
EPSILON = 0.6620
KAPPAD = 0.2854
C = 273.15
L0s = 2.56313e6
L1s = 1754
K2 = 1.137e6

# RUN-SPECIFIC CONSTANTS
PMASS = 2548740090557.712  # average parcel mass


###------------------------------------------------------------------------------
###------------------------------------------------------------------------------
## METEOROLOGY FUNCTIONS
###------------------------------------------------------------------------------
###------------------------------------------------------------------------------


def calc_pres(rho_kgm3, q_kgkg, T_K):
    """
    INPUTS
        - rho_kgm3 : kg m-3,    density
        - q_kgkg : kg kg-1,     specific humidity
        - T_K : K,              temperature
    ACTION
        Pressure is obtained using the ideal gas law, to which end the
        virtual temperature is calculated, allowing to use the specific
        gas constant of dry air
    OUTPUTS
        - pres_Pa: Pa,             pressure
    """
    ## convert specific humidity to mixing ratio (exact)
    r_kgkg = -q_kgkg / (q_kgkg - 1)
    ## now calculate virtual temperature
    Tv_K = T_K * (1 + r_kgkg / EPSILON) / (1 + r_kgkg)
    ## use ideal gas to obtain pressure
    return rho_kgm3 * RSPECIFIC * Tv_K


def dist_on_sphere(lat1, lon1, lat2, lon2):
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
    la1 = (PI / 180) * lat1
    la2 = (PI / 180) * lat2
    lo1 = (PI / 180) * lon1
    lo2 = (PI / 180) * lon2

    ## third, convert to cartesian coordinates
    ## use unit sphere for simplicity (r=1)
    x1 = sin(la1) * cos(lo1)
    y1 = sin(la1) * sin(lo1)
    z1 = cos(la1)

    x2 = sin(la2) * cos(lo2)
    y2 = sin(la2) * sin(lo2)
    z2 = cos(la2)

    ## fourth, calculate dot product and angle
    dotp = x1 * x2 + y1 * y2 + z1 * z2
    angle = acos(dotp)

    ## fifth, calculate distance
    dist = angle * EARTHRADIUS  # pi cancels out

    return dist


def dist_on_sphere2(lat1, lon1, lat2, lon2):

    lon1 = math.radians(lon1)
    lon2 = math.radians(lon2)
    dlon = lon2 - lon1

    lat1 = math.radians(lat1)
    lat2 = math.radians(lat2)
    dlat = lat2 - lat1

    a = (
        math.sin(dlat / 2) ** 2
        + math.cos(lat1) * math.cos(lat2) * math.sin(dlon / 2) ** 2
    )
    c = 2 * math.atan2(sqrt(a), sqrt(1 - a))

    dist = c * EARTHRADIUS

    return dist


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
    lats1 = (lats_centr + (res / 2)) * np.pi / 180  # np.sin requires radians
    lats2 = (lats_centr - (res / 2)) * np.pi / 180
    areas = (
        (np.pi / 180) * (EARTHRADIUS ** 2) * np.abs(np.sin(lats1) - np.sin(lats2)) * res
    )

    ## overwrite any areas of 0 (at the poles) with np.NaN to prevent problems
    try:
        areas[np.where(areas == 0.0)] = np.NaN  # only works for arrays
    except TypeError:
        pass  # simply ignore if it's a float
    # return array of dimension nlat x nlon
    ary_area = np.swapaxes(np.tile(areas, (nlon, 1)), 0, 1)
    return ary_area


def midpoint_on_sphere(lat1, lon1, lat2, lon2):

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
    lat1 = (PI / 180) * lat1
    lat2 = (PI / 180) * lat2
    lon1 = (PI / 180) * lon1
    lon2 = (PI / 180) * lon2

    ## third, convert to cartesian coordinates
    ## use unit sphere for simplicity (r=1)
    x1 = sin(lat1) * cos(lon1)
    y1 = sin(lat1) * sin(lon1)
    z1 = cos(lat1)

    x2 = sin(lat2) * cos(lon2)
    y2 = sin(lat2) * sin(lon2)
    z2 = cos(lat2)

    ## fourth, add vector to obtain "middle" vector, scale back to length 1
    middle_vec = [x1 + x2, y1 + y2, z1 + z2]
    mv = np.asarray(middle_vec) / sqrt(sum(np.asarray(middle_vec) ** 2))

    ## TRANSFORM BACK TO RADIAN
    potT = atan(sqrt((mv[0] ** 2) + (mv[1] ** 2)) / (mv[2]))
    phi = atan2(mv[1], mv[0])  # is "same" as atan(y/x), but NEED atan2 (sign!)

    ## convert these mid coordinates back to (arc) degree & shift
    if potT > 0:
        lat_mid = 90 - potT * (180 / PI)
    else:
        lat_mid = -(90 + potT * (180 / PI))
    if phi > 0:
        lon_mid = phi * (180 / PI) - 180
    else:
        lon_mid = 180 + phi * (180 / PI)

    if lon_mid > 179.5:
        lon_mid -= 360  # now shift all coords that otherwise would be allocated to +180 deg to - 180
    return (lat_mid, lon_mid)


def midpoint_on_sphere2(lat1, lon1, lat2, lon2):
    # Input values as degrees
    # following http://www.movable-type.co.uk/scripts/latlong.html
    lat1 = math.radians(lat1)
    lon1 = math.radians(lon1)
    lat2 = math.radians(lat2)
    lon2 = math.radians(lon2)
    bx = math.cos(lat2) * math.cos(lon2 - lon1)
    by = math.cos(lat2) * math.sin(lon2 - lon1)
    latm = math.atan2(
        math.sin(lat1) + math.sin(lat2),
        math.sqrt((math.cos(lat1) + bx) * (math.cos(lat1) + bx) + by ** 2),
    )
    lonm = lon1 + math.atan2(by, math.cos(lat1) + bx)
    londeg = math.degrees(lonm)
    if londeg > 180:
        londeg -= 360  # now shift all coords that otherwise would be allocated to +180 deg to - 180
    latdeg = math.degrees(latm)
    return latdeg, londeg


def q2rh(q_kgkg, p_Pa, T_K):
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
    e = q_kgkg * p_Pa / (0.622 + 0.378 * q_kgkg)

    # compute saturation vapor pressure, must convert T to °C
    es = 611.2 * np.exp(17.67 * (T_K - TREF) / (T_K - TREF + 243.5))

    # return relative humidity
    return 1e2 * e / es


def calc_pottemp_e(p_Pa, q_kgkg, T_K):
    """
    INPUTS
        - p_Pa : hPa,           pressure
        - q_kgkg : kg kg-1,     specific humidity
        - T_K : K,              temperature
    
    ACTION
        Dewpoint temperature is calculated to approximate the temperature at
        the lifting condensation level, which is then used to find potT there,
        and at last calculate potT-e.
        The procedure is essentially identical to what is used by MetPy,
        with the exception of making use of eq. 6.5 instead of eq. B39
        (see Davies-Jones, 2009; doi=10.1175/2009MWR2774.1)
    
    OUTPUTS
        - potT-e: K,   equivalent potential temperature
    """
    ## convert specific humidity to mixing ratio (exact)
    r_kgkg = -q_kgkg / (q_kgkg - 1)  # needs q in kg/kg
    r_gkg = r_kgkg * 1e3

    ## calculate potT using Bolton eq. 7
    p_hPa = p_Pa / 1e2
    potT = T_K * (1000 / p_hPa) ** (0.2854 * (1 - 0.00028 * r_gkg))

    ## convert q into e, vapor pressure (Bolton 1980)
    e_Pa = q_kgkg * p_Pa / (0.622 + 0.378 * q_kgkg)

    ## calculate dewpoint temperature according to inverted Bolton 1980 eq.,
    ## https://archive.eol.ucar.edu/projects/ceop/dm/documents/refdata_report/eqns.html
    e_hPa = e_Pa / 1e2
    T_D = np.log(e_hPa / 6.112) * 243.5 / (17.67 - np.log(e_hPa / 6.112)) + 273.15

    ## temperature at lifting condensation level (all temperatures in KElVIN; eq. 15)
    T_L = (1.0 / (1.0 / (T_D - 56.0) + np.log(T_K / T_D) / 800.0)) + 56.0

    ## potT at lifting condensation level (eq. 3.19, Davies-Jones 2009)
    potT_DL = potT * (
        ((potT / T_L) ** (0.28 * r_kgkg)) * (1 + (r_kgkg / EPSILON)) ** KAPPAD
    )

    ## potT at lifting condensation level (eq. 6.5, Davies-Jones 2009)
    potT_E = potT_DL * np.exp(
        (L0s - L1s * (T_L - C) + K2 * r_kgkg) * r_kgkg / (CPD * T_L)
    )

    return potT_E


def calc_pottemp(p_Pa, q_kgkg, T_K):
    """
    INPUTS
        - p_Pa : hPa,           pressure
        - q_kgkg : kg kg-1,     specific humidity
        - T_K : K,              temperature
    
    ACTION
        The potential temperature is calculated using Bolton's eq. 7;
        to this end, specific humidity is converted into mixing ratio
    
    OUTPUTS
        - potT: K,     potential temperature
    """
    ## convert specific humidity to mixing ratio (exact)
    r_kgkg = -q_kgkg / (q_kgkg - 1)  # needs q in kg/kg
    r_gkg = r_kgkg * 1e3

    ## calculate potT using Bolton eq. 7
    p_hPa = p_Pa / 1e2
    return T_K * (1000 / p_hPa) ** (0.2854 * (1 - 0.00028 * r_gkg))


def moist_ascender(p_Pa, q_kgkg, T_K):
    """    
    INPUTS
        - p_Pa : hPa,           pressure
        - q_kgkg : kg kg-1,     specific humidity
        - T_K : K,              temperature
    
    ACTION
        Calculate lifting condensation level (LCL) temperature analogously to
        potT_e,
        then approximate lapse rate at LCL using the approximation from:
            http://glossary.ametsoc.org/wiki/Saturation-adiabatic_lapse_rate
    
    OUTPUTS
        - T_LCL: K,     temperature at LCL
        - MALR: K m-1,  moist adiabatic lapse rate
    
    """
    ## convert q into e, vapor pressure (Bolton 1980)
    e_Pa = q_kgkg * p_Pa / (0.622 + 0.378 * q_kgkg)

    ## calculate dewpoint temperature according to inverted Bolton 1980 eq.,
    ## https://archive.eol.ucar.edu/projects/ceop/dm/documents/refdata_report/eqns.html
    e_hPa = e_Pa / 1e2
    T_D = np.log(e_hPa / 6.112) * 243.5 / (17.67 - np.log(e_hPa / 6.112)) + 273.15

    ## temperature at lifting condensation level (all temperatures in KElVIN; eq. 15)
    T_LCL = (1.0 / (1.0 / (T_D - 56.0) + np.log(T_K / T_D) / 800.0)) + 56.0

    ## now that we have T_LCL, calculate moist adiabatic lapse rate...

    ## convert specific humidity to mixing ratio (exact)
    r_kgkg = -q_kgkg / (q_kgkg - 1)  # needs q in kg/kg
    ## now, finally, approximate lapse rate
    MALR = (
        GRAV
        * ((RSPECIFIC * (T_LCL ** 2)) + (HV * r_kgkg * T_LCL))
        / ((CPD * RSPECIFIC * T_LCL ** 2) + (HV ** 2 * r_kgkg * EPS))
    )

    return (T_LCL, -MALR)  # T_LCL in Kelvin, MALR in K/m


def makegrid(resolution):
    glat = np.arange(-90, 90 + resolution, resolution)
    glon = np.arange(-180, 180, resolution)
    gd_area = gridded_area_exact(glat, res=resolution, nlon=glon.size)
    return glon, glat, gd_area


###------------------------------------------------------------------------------
###------------------------------------------------------------------------------
## HAMSTER-SPECIFIC FUNCTIONS
###------------------------------------------------------------------------------
###------------------------------------------------------------------------------
def str2bol(v):
    """
    ACTION: converts (almost) any string to boolean False/True
    NOTE:   needed for boolean interpretation in parser.add_argument from parsearg in read_cmdargs
    """
    if isinstance(v, bool):
        return v
    if v.lower() in ("yes", "true", "t", "y", "1"):
        return True
    elif v.lower() in ("no", "false", "f", "n", "0"):
        return False
    else:
        raise argparse.ArgumentTypeError("Boolean value expected.")


def read_cmdargs():
    """
    ACTION: read dates, thresholds and flags from command line
    RETURN: 'args' contains all
    DEP:    uses argparse
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--pathfile",
        "-pf",
        help="name of pathfile (default: paths.txt)",
        metavar="",
        type=str,
        default="paths.txt",
    )
    parser.add_argument(
        "--steps",
        "-st",
        help="steps performed (0: flex2traj, 1: diagnosis, 2: attribution, 3: bias correction)",
        metavar="",
        type=int,
        default=1,
    )
    parser.add_argument(
        "--ayyyy",
        "-ay",
        help="analysis year (YYYY)",
        metavar="",
        type=int,
        default=2002,
    )
    parser.add_argument(
        "--am", "-am", help="analysis month (M)", metavar="", type=int, default=1
    )
    parser.add_argument(
        "--ad", "-ad", help="analysis day (D)", metavar="", type=int, default=1
    )
    parser.add_argument(
        "--mode", "-m", help="mode (test,oper)", metavar="", type=str, default="oper"
    )
    parser.add_argument(
        "--expid",
        "-id",
        help="experiment ID (string, example versionA)",
        metavar="",
        type=str,
        default="FXv",
    )
    parser.add_argument(
        "--maskval",
        "-mv",
        help="use <value> from maskfile for masking",
        metavar="",
        type=int,
        default=1,
    )
    parser.add_argument(
        "--ctraj_len",
        "-len",
        help="threshold for maximum allowed trajectory length in days",
        metavar="",
        type=int,
        default=10,
    )
    parser.add_argument(
        "--cprec_dqv",
        "-cpq",
        help="threshold for detection of P based on delta(qv)",
        metavar="",
        type=float,
        default=0,
    )
    parser.add_argument(
        "--fprec",
        "-fp",
        help="flag: filter for P (01 only)",
        metavar="",
        type=str2bol,
        default=True,
        nargs="?",
    )
    parser.add_argument(
        "--cprec_rh",
        "-cpr",
        help="threshold for detection of P based on RH",
        metavar="",
        type=float,
        default=80,
    )
    parser.add_argument(
        "--fevap",
        "-fe",
        help="flag: filter for E (01 only)",
        metavar="",
        type=str2bol,
        default=True,
        nargs="?",
    )
    parser.add_argument(
        "--cevap_dqv",
        "-ceq",
        help="threshold for detection of E based on delta(qv)",
        metavar="",
        type=float,
        default=0,
    )
    parser.add_argument(
        "--cevap_hgt",
        "-ceh",
        help="threshold for detection of E using a maximum height",
        metavar="",
        type=float,
        default=0,
    )
    parser.add_argument(
        "--fevap_drh",
        "-fer",
        help="flag: check for maximum delta(RH) for detection of E",
        metavar="",
        type=str2bol,
        default=False,
        nargs="?",
    )
    parser.add_argument(
        "--cevap_drh",
        "-cer",
        help="threshold for detection of E using a maximum delta(RH)",
        metavar="",
        type=float,
        default=15,
    )
    parser.add_argument(
        "--fheat",
        "-fh",
        help="flag: filter for H (01 only)",
        metavar="",
        type=str2bol,
        default=True,
        nargs="?",
    )
    parser.add_argument(
        "--cheat_hgt",
        "-chh",
        help="threshold for detection of H using a maximum height",
        metavar="",
        type=float,
        default=0,
    )
    parser.add_argument(
        "--cheat_dtemp",
        "-cht",
        help="threshold for detection of H using a minimum delta(T)",
        metavar="",
        type=float,
        default=0,
    )
    parser.add_argument(
        "--fheat_drh",
        "-fhr",
        help="flag: check for maximum delta(RH) for detection of H",
        metavar="",
        type=str2bol,
        default=False,
        nargs="?",
    )
    parser.add_argument(
        "--cheat_drh",
        "-chr",
        help="threshold for detection of H using a maximum delta(RH)",
        metavar="",
        type=float,
        default=15,
    )
    parser.add_argument(
        "--fheat_rdq",
        "-fhq",
        help="flag: check for maximum relative delta(Q) for detection of H",
        metavar="",
        type=str2bol,
        default=False,
        nargs="?",
    )
    parser.add_argument(
        "--cheat_rdq",
        "-chq",
        help="threshold for detection of H using a maximum relative d(Q) [%]",
        metavar="",
        type=float,
        default=10,
    )
    parser.add_argument(
        "--cpbl_factor",
        "-pblf",
        help="factor for PBL relaxation",
        metavar="",
        type=float,
        default=1,
    )
    parser.add_argument(
        "--cpbl_method",
        "-pblm",
        help="filter for PBL: mean, max, actual heights between 2 points",
        metavar="",
        type=str,
        default="max",
    )
    parser.add_argument(
        "--cpbl_strict",
        "-pbls",
        help="filter for PBL: 0/1/2 locations within max PBL (0: no filter)",
        metavar="",
        type=int,
        default=1,
    )
    parser.add_argument(
        "--timethis",
        "-t",
        help="time the main loop (flag)",
        metavar="",
        type=str2bol,
        default=False,
        nargs="?",
    )
    parser.add_argument(
        "--write_netcdf",
        "-o",
        help="write netcdf output (flag)",
        metavar="",
        type=str2bol,
        default=True,
        nargs="?",
    )
    parser.add_argument(
        "--pref_data",
        "-pref",
        help="preciptation bias correction data set (eraint/others)",
        metavar="",
        type=str,
        default="eraint",
    )
    parser.add_argument(
        "--eref_data",
        "-eref",
        help="evaporation bias correction data set (eraint/others)",
        metavar="",
        type=str,
        default="eraint",
    )
    parser.add_argument(
        "--href_data",
        "-href",
        help="sensible heat bias correction data set (eraint/others)",
        metavar="",
        type=str,
        default="eraint",
    )
    parser.add_argument(
        "--write_month",
        "-mo",
        help="write monthly aggreagted netcdf output (03 only; flag)",
        metavar="",
        type=str2bol,
        default=False,
        nargs="?",
    )
    parser.add_argument(
        "--precision",
        "-f",
        help="precision for writing netcdf file variables (f4,f8)",
        metavar="",
        type=str,
        default="f8",
    )
    parser.add_argument(
        "--verbose",
        "-v",
        help="verbose output (flag)",
        metavar="",
        type=str2bol,
        default=True,
        nargs="?",
    )
    parser.add_argument(
        "--veryverbose",
        "-vv",
        help="very verbose output (flag)",
        metavar="",
        type=str2bol,
        default=False,
        nargs="?",
    )
    parser.add_argument(
        "--fallingdry",
        "-dry",
        help="cut off trajectories falling dry (flag)",
        metavar="",
        type=str2bol,
        default=False,
        nargs="?",
    )
    parser.add_argument(
        "--memento",
        "-mto",
        help="keep track of trajectory history (02 only - needed for Had; flag)",
        metavar="",
        type=str2bol,
        default=True,
        nargs="?",
    )
    parser.add_argument(
        "--mattribution",
        "-matt",
        help="attribution method (for E2P as of now: random/linear)",
        metavar="",
        type=str,
        default="linear",
    )
    parser.add_argument(
        "--ratt_nit",
        "-rnit",
        help="minimum number of iterations for random attribution",
        metavar="",
        type=int,
        default=10,
    )
    parser.add_argument(
        "--ratt_forcall",
        "-rall",
        help="enforcing the attribution to all uptake locations (random att.)",
        metavar="",
        type=str2bol,
        default=False,
        nargs="?",
    )
    parser.add_argument(
        "--ratt_wloc",
        "-rwloc",
        help="weight probability of locations according to their maximum potential contribution (random att.)",
        metavar="",
        type=str2bol,
        default=True,
        nargs="?",
    )
    parser.add_argument(
        "--explainp",
        "-exp",
        help="trajectory-based upscaling of E2P contributions (02 only: none/max/full)",
        metavar="",
        type=str,
        default="none",
    )
    parser.add_argument(
        "--dupscale",
        "-dups",
        help="daily upscaling of E2P contributions (02 only; flag)",
        metavar="",
        type=str2bol,
        default=False,
        nargs="?",
    )
    parser.add_argument(
        "--mupscale",
        "-mups",
        help="monthly upscaling of E2P contributions (02 only; flag)",
        metavar="",
        type=str2bol,
        default=False,
        nargs="?",
    )
    parser.add_argument(
        "--bc_useattp",
        "-uatt",
        help="use precipitation from attribution for bias-correction",
        metavar="",
        type=str2bol,
        default=False,
        nargs="?",
    )
    parser.add_argument(
        "--bc_time",
        "-bct",
        help="time scale for bias-correction (daily/monthly)",
        metavar="",
        type=str,
        default="daily",
    )
    parser.add_argument(
        "--variable_mass",
        "-vm",
        help="use variable mass (flag)",
        metavar="",
        type=str2bol,
        default=False,
        nargs="?",
    )
    parser.add_argument(
        "--writestats",
        "-ws",
        help="write additional stats to file (02 and 03 only; flag)",
        metavar="",
        type=str2bol,
        default=False,
        nargs="?",
    )
    parser.add_argument(
        "--bc_aggbwtime",
        "-aggbt",
        help="aggregate backward time (03 only; flag)",
        metavar="",
        type=str2bol,
        default=True,
        nargs="?",
    )
    parser.add_argument(
        "--bc_e2p",
        "-e2p",
        help="bias correction: write E2P (raw) to file (flag)",
        metavar="",
        type=str2bol,
        default=True,
        nargs="?",
    )
    parser.add_argument(
        "--bc_e2p_p",
        "-e2pp",
        help="bias correction: write E2P_Ps to file (flag)",
        metavar="",
        type=str2bol,
        default=True,
        nargs="?",
    )
    parser.add_argument(
        "--bc_e2p_e",
        "-e2pe",
        help="bias correction: write E2P_Es to file (flag)",
        metavar="",
        type=str2bol,
        default=False,
        nargs="?",
    )
    parser.add_argument(
        "--bc_e2p_ep",
        "-e2pep",
        help="bias correction: write E2P_EPs to file (flag)",
        metavar="",
        type=str2bol,
        default=True,
        nargs="?",
    )
    parser.add_argument(
        "--bc_t2p_ep",
        "-t2pep",
        help="bias correction: calculate transpiration to P (T2P_EPs, flag)",
        metavar="",
        type=str2bol,
        default=False,
        nargs="?",
    )
    parser.add_argument(
        "--bc_had",
        "-had",
        help="bias correction: write Had (raw) to file (flag)",
        metavar="",
        type=str2bol,
        default=True,
        nargs="?",
    )
    parser.add_argument(
        "--bc_had_h",
        "-hadh",
        help="bias correction: write Had_Hs to file (flag)",
        metavar="",
        type=str2bol,
        default=True,
        nargs="?",
    )
    parser.add_argument(
        "--debug",
        "-d",
        help="debugging option (flag)",
        metavar="",
        type=str2bol,
        default=False,
        nargs="?",
    )
    parser.add_argument(
        "--gres",
        "-r",
        help="output grid resolution (degrees)",
        metavar="",
        type=float,
        default=1,
    )
    parser.add_argument(
        "--ryyyy",
        "-ry",
        help="run name (here, YYYY, example: 2002, default: ayyyy)",
        metavar="",
        type=int,
        default=None,
    )
    parser.add_argument(
        "--refdate",
        "-rd",
        help="reference date (YYYYMMDDHH)",
        metavar="",
        type=str,
        default=None,
    )
    parser.add_argument(
        "--waiter",
        "-wt",
        help="random waiter to avoid simulatenous access to files",
        metavar="",
        type=str2bol,
        default=True,
        nargs="?",
    )
    parser.add_argument(
        "--fproc_npart",
        "-fpn",
        help="01: process all parcels (n parcel evaluation)",
        metavar="",
        type=str2bol,
        default=True,
        nargs="?",
    )
    # print(parser.format_help())
    args = parser.parse_args()  # namespace
    # handle None cases already
    if args.ryyyy is None:
        args.ryyyy = args.ayyyy
    if args.refdate is None:
        args.refdate = str(args.ryyyy) + "123118"
    return args


def printsettings(args):

    str0 = str(
        "[[DATES]] ayyyy = "
        + str(args.ayyyy)
        + ", am = "
        + str(args.am)
        + ", ad = "
        + str(args.ad)
        + ", ryyyy = "
        + str(args.ryyyy)
        + " "
        + "[[GRID]] gres = "
        + str(args.gres)
        + " "
        + "[[PATHS & MASK]] pathfile = "
        + str(args.pathfile)
        + ", maskval = "
        + str(args.maskval)
        + " "
        + "[[EXPERIMENT ID]] expid = "
        + str(args.expid)
    )

    str1 = str(
        "Diagnosis with the following settings: "
        + "[[PRECIPITATION]] cprec_dqv = "
        + str(args.cprec_dqv)
        + ", cprec_rh = "
        + str(args.cprec_rh)
        + " "
        + "[[EVAPORATION]] cevap_hgt = "
        + str(args.cevap_hgt)
        + ", cevap_dqv = "
        + str(args.cevap_dqv)
        + ", fevap_drh = "
        + str(args.fevap_drh)
        + ", cevap_drh = "
        + str(args.cevap_drh)
        + " "
        + "[[SENSIBLE HEAT]] cheat_hgt = "
        + str(args.cheat_hgt)
        + ", cheat_dtemp = "
        + str(args.cheat_dtemp)
        + ", fheat_drh = "
        + str(args.fheat_drh)
        + ", cheat_drh = "
        + str(args.cheat_drh)
        + " "
        + "[[OTHERS]]: cpbl_strict = "
        + str(args.cpbl_strict)
        + ", cpbl_method = "
        + str(args.cpbl_method)
        + ", cpbl_factor = "
        + str(args.cpbl_factor)
        + ", variable_mass = "
        + str(args.variable_mass)
        + ", mode = "
        + str(args.mode)
    )

    str2 = str(
        "Attribution with the following settings: "
        + "[[ATTRIBUTION]]: ctraj_len = "
        + str(args.ctraj_len)
        + ", fallingdry = "
        + str(args.fallingdry)
        + ", memento (H) = "
        + str(args.memento)
        + ", attribution = "
        + str(args.mattribution)
        + ", explainp (P) = "
        + str(args.explainp)
        + ", dupscale (P) = "
        + str(args.dupscale)
        + ", mupscale (P) = "
        + str(args.mupscale)
    )
    if args.mattribution == "random":
        str2 = str2 + str(
            " [random attribution settings] ratt_nit = "
            + str(args.ratt_nit)
            + ", ratt_forcall = "
            + str(args.ratt_forcall)
        )

    str3 = str(
        "Bias correction with the following settings: "
        + "[[BIAS CORRECTION]]: bc_time = "
        + str(args.bc_time)
        + ", bc_useattp = "
        + str(args.bc_useattp)
        + ", bc_aggbwtime = "
        + str(args.bc_aggbwtime)
        + ", write_month = "
        + str(args.write_month)
    )

    if args.steps == 1:
        return str0 + "; " + str1
    if args.steps == 2:
        return str0 + "; " + str1 + "; " + str2
    if args.steps == 3:
        return (
            str3  # copying ncdf description from 2 in step 3; thus only passing BC info
        )


def readtraj(idate, ipath, ifile_base, verbose=True):  # date as string [YYYYMMDDHH]
    # reads in *h5 data from flex2traj and flips the time axis
    # returns data array of dimension (ntrajlength x nparticles x nvars)
    # Check if file exists /file format
    ifile = str(ipath + "/" + ifile_base + "_" + idate + ".h5")
    if not os.path.isfile(ifile):
        raise SystemExit(ifile + " does not exist!")
    elif os.path.isfile(ifile):
        if verbose:
            print(" Reading " + ifile)
        with h5py.File(ifile, "r") as f:
            dataar = np.array(f["trajdata"])
    # flip time axis !!!
    dataar = dataar[::-1, :, :]
    return dataar


def checkpbl(cpbl, ztra, hpbl, maxhgt):
    if cpbl == 0:
        # do not check for PBL; use everything
        return True
    if cpbl == 1:
        # at least one location is within the max PBL
        if (ztra < max(hpbl[1], hpbl[0], maxhgt)).any():
            return True
        else:
            return False
    if cpbl == 2:
        # both locations are within the max PBL
        if (ztra < max(hpbl[1], hpbl[0], maxhgt)).all():
            return True
        else:
            return False


def readparcel(parray):
    lats = parray[:, 2]  # latitude
    lons = parray[:, 1]  # longitude
    lons[
        lons > 180.0
    ] -= 360  # transform coordinates from [0 ... 360] to [-180 ... 180]
    temp = parray[:, 8]  # temperature (K)
    ztra = parray[:, 3]  # height (m)
    # topo   = parray[:,4]                   # topography (m)
    qv = parray[:, 5]  # specific humidity (kg kg-1)
    hpbl = parray[:, 7]  # ABL height (m)
    dens = parray[:, 6]  # density (kg m-3)
    pres = calc_pres(dens, qv, temp)  # pressure (Pa)
    pottemp = calc_pottemp(pres, qv, temp)  # potential temperature (K)
    epottemp = calc_pottemp_e(pres, qv, temp)  # equivalent potential temperature (K)

    return lons, lats, temp, ztra, qv, hpbl, dens, pres, pottemp, epottemp


def calc_allvars(parray):
    # calculate all additional variables required for heat + moisture analysis; for all parcels!
    # hardcoded on 3D array and variable structure...
    # pressure (10th variable, index 9)
    parray = np.dstack(
        (parray, calc_pres(parray[:, :, 6], parray[:, :, 5], parray[:, :, 8]))
    )
    # relative humidity (11th variable, index 10)
    parray = np.dstack(
        (parray, q2rh(parray[:, :, 5], parray[:, :, 9], parray[:, :, 8]))
    )
    # potential temperature (12th variable, index 11)
    parray = np.dstack(
        (parray, calc_pottemp(parray[:, :, 9], parray[:, :, 5], parray[:, :, 8]))
    )
    return parray


def readmidpoint(parray):
    lats = parray[:, 2]  # latitude
    lons = parray[:, 1]  # longitude
    lons[
        lons > 180.0
    ] -= 360  # transform coordinates from [0 ... 360] to [-180 ... 180]
    lat_mid, lon_mid = midpoint_on_sphere(
        lats[0], lons[0], lats[1], lons[1]
    )  # use own function to calculate midpoint position
    if lon_mid > 179.5:
        lon_mid -= 360  # now shift all coords that otherwise would be allocated to +180 deg to - 180
    return lat_mid, lon_mid


def readsparcel(parray):
    # reads standard output that is always needed
    qv = parray[:, 5]  # specific humidity (kg kg-1)
    temp = parray[:, 8]  # temperature (K)
    hpbl = parray[:, 7]  # ABL height (m)
    ztra = parray[:, 3]  # height (m)
    return qv, temp, ztra, hpbl


def readpres(parray):
    temp = parray[:, 8]  # temperature (K)
    dens = parray[:, 6]  # density (kg m-3)
    qv = parray[:, 5]  # specific humidity (kg kg-1)
    pres = calc_pres(dens, qv, temp)  # pressure (Pa)
    return pres


def readepottemp(parray):
    temp = parray[:, 8]  # temperature (K)
    qv = parray[:, 5]  # specific humidity (kg kg-1)
    dens = parray[:, 6]  # density (kg m-3)
    pres = calc_pres(dens, qv, temp)  # pressure (Pa)
    epottemp = calc_pottemp_e(pres, qv, temp)  # equivalent potential temperature (K)
    return epottemp


def readpottemp(parray):
    qv = parray[:, 5]  # specific humidity (kg kg-1)
    dens = parray[:, 6]  # density (kg m-3)
    temp = parray[:, 8]  # temperature (K)
    pres = calc_pres(dens, qv, temp)  # pressure (Pa)
    pottemp = calc_pottemp(pres, qv, temp)  # potential temperature (K)
    return pottemp


def readheights(parray):

    """
    This is basically the same as 'glanceparcel',
    and as of now, is called only for the first 4 timesteps.
 
    NOTE: probably to be updated / removed

    """

    ## parcel information
    ztra = parray[:, 3]  # height (m)
    hpbl = parray[:, 7]  # ABL height (m)

    return ztra, hpbl


def glanceparcel(parray):

    """
    This is used to take a first glance at parcel properties,
    and as of now, is called only for the first 4 timesteps.
 
    Required for arriving air criterion (which goes back 4 steps).
    However, this could be further optimized by only getting the
    last 4 values of hpbl, as for the other variables here,
    only the last / last 2 values are required at first.

    NOTE: beautify or remove!

    """

    ## parcel information
    ztra = parray[:, 3]  # height (m)
    hpbl = parray[:, 7]  # ABL height (m)
    temp = parray[:, 8]  # temperature (K)
    qv = parray[:, 5]  # specific humidity (kg kg-1)
    dens = parray[:, 6]  # density (kg m-3)
    pres = calc_pres(dens, qv, temp)  # pressure (Pa)

    return ztra, hpbl, temp, qv, dens, pres


def trajparceldiff(pvals, meval):
    # difference
    if meval in ["diff"]:
        dpval = pvals[:-1] - pvals[1:]
    # mean
    if meval in ["mean"]:
        dpval = (pvals[:-1] + pvals[1:]) / 2
    # max
    if meval in ["max"]:
        dpval = np.amax(pvals[:-1], pvals[1:])
    return dpval


def get_refnpart(refdate, ryyyy, glon, glat):
    """
    INPUT
        - refdate [YYYYMMDDHH] :    reference date (str) used for counting midpoint parcels and scaling
        - glon, glat :              reference grid coordinates      
    ACTION
        - calculates the reference distribution of parcels using the midpoint of parcels at refdate
        - NOTE that this is run specific and needs to be adjusted if FLEXPART runs are setup differently
    DEPEND
        - uses numpy and functions readtraj, gridder
    RETURN
        - npart (nlat x nlon) at refdate
    """
    if verbose:
        # print("Reference number of particles: \t" + str(nparticle))
        print(" * Getting reference distribution...")

    ary_npart = np.zeros(shape=(glat.size, glon.size))
    ary = readtraj(
        idate=refdate,
        ipath="/scratch/gent/vo/000/gvo00090/D2D/data/FLEXPART/era_global/flex2traj_t2/"
        + str(ryyyy),
        ifile_base="global",
    )
    nparticle = ary.shape[1]
    for i in range(nparticle):
        lons, lats, _, _, _, _, _, _, _, _ = readparcel(ary[:, i, :])
        ary_npart[:, :] += gridder(
            plon=lons, plat=lats, pval=int(1), glon=glon, glat=glat
        )

    return ary_npart


def convertunits(ary_val, garea, var):
    """
    INPUT
        - aryval
    ACTION
        - calculates grid cell values 
    RETURN
        - returns P and E as mm
        - returns H as W m-2
    """
    if var in ["P", "E"]:
        return PMASS * ary_val / (1e6 * garea)
    if var in ["H"]:
        return PMASS * ary_val * CPD / (1e6 * garea * 6 * 3600)


def movingmax(x, n=2):
    if len(x) == 2:
        return max(x[0], x[1])
    else:
        return np.array([np.max(x[i : i + n]) for i in range(len(x) - (n - 1))])


def movingmean(x, n=2):
    if len(x) == 2:
        return mean(x[0], x[1])
    else:
        return np.array([np.mean(x[i : i + n]) for i in range(len(x) - (n - 1))])


def pblcheck(ary, cpbl_strict, minh, fpbl, method):
    # returns boolean vector for all change locations (z.size-1)
    # manually tweak PBL heights to account for minimum heights (attn; if fpbl != 1; the heights are adjusted)
    z = ary[:, 0]
    # (i) skip everything if no pbl check required
    if cpbl_strict == 0:
        return np.ones(dtype=bool, shape=z.size - 1)
    else:
        hpbl = ary[:, 1]
        hpbl[hpbl < minh] = minh
        if method == "mean":
            before_inside = z[1:] <= fpbl * movingmean(hpbl, n=2)
            after_inside = z[:-1] <= fpbl * movingmean(hpbl, n=2)
        elif method == "max":
            before_inside = z[1:] <= fpbl * movingmax(hpbl, n=2)
            after_inside = z[:-1] <= fpbl * movingmax(hpbl, n=2)
        elif method == "actual":
            before_inside = z[1:] <= fpbl * hpbl[1:]
            after_inside = z[:-1] <= fpbl * hpbl[:-1]
        if cpbl_strict == 2:
            # both inside (and)
            return np.logical_and(before_inside, after_inside)
        elif cpbl_strict == 1:
            # one inside (or)
            return np.logical_or(before_inside, after_inside)


def pblcheck_diag(z, hpbl, cpbl_strict, minh, fpbl, method):
    hpbl[hpbl < minh] = minh
    if cpbl_strict == 0:
        return np.asarray([range(z.shape[1])])
    if method == "actual":
        fpbl1 = z[0, :] <= fpbl * hpbl[0, :]
        fpbl2 = z[1, :] <= fpbl * hpbl[1, :]
    elif method == "mean":
        mpbl = np.apply_over_axes(np.mean, hpbl[:, :], 0)[0, :]
        fpbl1 = z[0, :] <= fpbl * mpbl
        fpbl2 = z[1, :] <= fpbl * mpbl
    elif method == "max":
        mpbl = np.apply_over_axes(np.max, hpbl[:, :], 0)[0, :]
        fpbl1 = z[0, :] <= fpbl * mpbl
        fpbl2 = z[1, :] <= fpbl * mpbl
    if cpbl_strict == 1:
        fpbl = np.where(np.logical_or(fpbl1, fpbl2))
    if cpbl_strict == 2:
        fpbl = np.where(np.logical_and(fpbl1, fpbl2))
    return fpbl


def drhcheck(rh, checkit, maxdrh):
    if not checkit:
        retvals = np.ones(dtype=bool, shape=rh.size - 1)
    elif checkit:
        drh = trajparceldiff(rh, "diff")
        retvals = np.abs(drh) <= maxdrh
    return retvals


def drhcheck_diag(ary2d, checkit, maxdrh):
    if not checkit:
        fdrh = np.asarray([range(len(ary2d[0, :]))])
    else:
        fdrh = np.where(abs(ary2d[0, :] - ary2d[1, :]) <= maxdrh)
    return fdrh


def rdqvcheck(qv, checkit, maxrdqv):
    # checks if the absolute humidity changed by less than maxrdqv %
    if not checkit:
        retvals = np.ones(dtype=bool, shape=qv.size - 1)
    elif checkit:
        dqv = trajparceldiff(qv, "diff")
        rdqv = abs(dqv) / qv[:-1]  # using most recent Q as a reference here
        retvals = 100 * np.abs(rdqv) <= maxrdqv
    return retvals


def rdqcheck_diag(ary2d, checkit, maxrdqv):
    # checks if the absolute humidity changed by less than maxrdqv %
    if not checkit:
        frdq = np.asarray([range(len(ary2d[0, :]))])
    elif checkit:
        pchange = 100 * abs(ary2d[0, :] - ary2d[1, :]) / ary2d[0, :]  # %
        frdq = np.where(pchange <= maxrdqv)
    return frdq


def filter_for_evap_parcels(
    eary,
    dq,
    cpbl_method,
    cpbl_strict,
    cpbl_factor,
    cevap_hgt,
    fevap_drh,
    cevap_drh,
    cevap_dqv,
    veryverbose,
):
    # check for dq > cevap_dqv
    fdqv = np.where(dq[0, :] > cevap_dqv)
    eary = eary[:, fdqv[0], :]
    # check for drh
    fdrh = drhcheck_diag(eary[:, :, 10], checkit=fevap_drh, maxdrh=cevap_drh)
    eary = eary[:, fdrh[0], :]
    # check for pbl
    fpbl = pblcheck_diag(
        eary[:, :, 3], eary[:, :, 7], cpbl_strict, cevap_hgt, cpbl_factor, cpbl_method
    )
    if veryverbose:
        print(" * Filter for E: " + str(fdqv[0].size) + " parcels after dqv-filter")
        print("                 " + str(fdrh[0].size) + " parcels after drh-filter")
        print("                 " + str(fpbl[0].size) + " parcels after pbl-filter")
    return fdqv[0][fdrh[0][fpbl[0]]]


def filter_for_heat_parcels(
    hary,
    dTH,
    cpbl_method,
    cpbl_strict,
    cpbl_factor,
    cheat_hgt,
    fheat_drh,
    cheat_drh,
    cheat_dtemp,
    fheat_rdq,
    cheat_rdq,
    veryverbose,
):
    # check for dTH > cheat_dtemp
    fdTH = np.where(dTH[0, :] > cheat_dtemp)
    hary = hary[:, fdTH[0], :]
    # check for drh
    fdrh = drhcheck_diag(hary[:, :, 10], checkit=fheat_drh, maxdrh=cheat_drh)
    hary = hary[:, fdrh[0], :]
    # check for rdq
    frdq = rdqcheck_diag(hary[:, :, 5], checkit=fheat_rdq, maxrdqv=cheat_rdq)
    hary = hary[:, frdq[0], :]
    # check for pbl
    fpbl = pblcheck_diag(
        hary[:, :, 3], hary[:, :, 7], cpbl_strict, cheat_hgt, cpbl_factor, cpbl_method
    )
    if veryverbose:
        print(" * Filter for H: " + str(fdTH[0].size) + " parcels after dTH-filter")
        print("                 " + str(fdrh[0].size) + " parcels after drh-filter")
        print("                 " + str(frdq[0].size) + " parcels after rdq-filter")
        print("                 " + str(fpbl[0].size) + " parcels after pbl-filter")
    return fdTH[0][fdrh[0][frdq[0][fpbl[0]]]]


def scale_mass(ary_val, ary_part, ary_rpart):
    ary_sval = np.zeros(shape=(ary_val.shape))
    # for itime in range(ary_val.shape[0]):
    #    with np.errstate(divide='ignore', invalid='ignore'):
    #        ary_sval[itime,:,:] = ary_val[itime,:,:] * np.nan_to_num( ary_rpart[:,:]/ary_part[itime,:,:] )
    with np.errstate(divide="ignore", invalid="ignore"):
        ary_sval[:, :] = ary_val[:, :] * np.nan_to_num(ary_rpart[:, :] / ary_part[:, :])
    return ary_sval


def linear_discounter(v, min_gain):
    """
    ========= inputs =========
    v:        vector extending into past; v[0] is most recent, v[-1] least recent
    min_gain: minimum gain (dv/dt) to be considered
    ==========================
    """

    ## compute dv/dt, prepare dv_disc
    dv = v[:-1] - v[1:]
    # append initial v as a -fake- uptake
    # (used to estimate fraction that cannot be attributed)
    dv = np.append(dv, v[-1])
    dv_disc = np.zeros(shape=len(dv))
    ## get indices of gains and losses
    idx_gains = np.where(dv >= min_gain)[0]
    idx_losses = np.where(dv < 0)[0]

    for idx in idx_gains:
        ## determine relevant indices (losses only!)
        idx_consider = idx_losses[np.where(idx_losses < idx)]

        ## skip current iteration if no moisture loss occurs between then and t=0
        if len(idx_consider) == 0:
            dv_disc[idx] = dv[idx]
            continue

        ## create vector with fractions,
        frc_vec = v[idx_consider] / v[idx_consider + 1]
        ## then multiply current dv with product of vector (order irrelevant)
        dv_disc[idx] = np.prod(frc_vec) * dv[idx]
    return dv_disc


def linear_attribution_p(qv, iupt, explainp):
    dq_disc = np.zeros(shape=qv.size)
    dq_disc[1:] = linear_discounter(v=qv[1:], min_gain=0)
    ## trajectory-based upscaling
    prec = abs(qv[0] - qv[1])
    fw_orig = dq_disc / qv[1]
    if explainp == "full":
        # upscaling to 100% of trajectory
        cfac = qv[1] / np.sum(dq_disc[iupt])
        etop = prec * fw_orig * cfac
    elif explainp == "max":
        # upscaling to (100-IC)% of trajectory
        cfac = (qv[1] - dq_disc[-1]) / np.sum(dq_disc[iupt])
        etop = prec * fw_orig * cfac
    elif explainp == "none":
        # no upscaling
        etop = prec * fw_orig
    return etop


def calc_maxatt(qtot, iupt, verbose):
    dqdt = qtot[:-1] - qtot[1:]
    dqdt = np.append(dqdt, qtot[-1])
    nt = len(dqdt)
    dqdt_max = np.zeros(shape=nt)
    for ii in iupt[::-1]:
        try:
            imin = np.argmin(qtot[1:ii]) + 1
        except:
            imin = 1
        iatt = qtot[imin] - round(np.sum(dqdt_max[imin:]), 8)
        idqdt = min(iatt, dqdt[ii] - dqdt_max[ii])
        dqdt_max[ii] = idqdt
    maxatt = np.sum(dqdt_max) / abs(dqdt[0])
    if maxatt < 1 and verbose:
        print(
            " * Maximum attribution along trajectory: {:.2f}".format(100 * maxatt) + "%"
        )
    return maxatt

def calc_maxcon(qtot, iupt, verbose):
    dqdt = qtot[:-1] - qtot[1:]
    dqdt = np.append(dqdt, qtot[-1])
    nt = len(dqdt)
    dqdt_max = np.zeros(shape=nt)
    for ii in iupt[::-1]:
        try:
            imin = np.argmin(qtot[1:ii]) + 1
        except:
            imin = 1
        dqdt_max[ii] = min(qtot[imin], dqdt[ii])
    return dqdt_max

def local_minima(x):
    return np.r_[True, x[1:] < x[:-1]] & np.r_[x[:-1] < x[1:], True]


def random_attribution_p(
    qtot, iupt, explainp, nmin=1, forc_all=False, weight_locations=True, verbose=True, veryverbose=False
):
    qtot = qtot * 1000
    # This is only coded for precipitation as of now
    # with:
    # qtot = specific humidity
    # iupt = identified uptake locations
    # explainp = none(default)/full analogue to linear discounting & attribution
    #  - none: maximum is attributed to iupt (can be 100%!)
    #  - max: maximum is attributed to iupt + initial condition (~not explained); enforcing at least one iteration on init. cond.
    #  - full: 100% is attributed to iupt (if possible!)
    # nmin = tuning parameter; ~minimum iterations
    #  - the higher this value, the more iterations, the uptake locations are covered
    #  - a value of 10 enforces min. 10 iterations
    # forc_all = enforce attribution to all uptake locations (but still random)
    # weight_locations = weighting of location picks according to their maximum potential contribution (True/False; False = Default)
    dqdt = qtot[:-1] - qtot[1:]
    # append initial condition as artificial uptake
    dqdt = np.append(dqdt, qtot[-1])
    nt = len(dqdt)
    if explainp == "max":
        iupt = np.append(iupt, nt - 1)
    # indicator for potential uptake locations (1: yes, 0: no)
    pupt = np.zeros(shape=nt)
    pupt[iupt] = 1
    nupt = len(np.where(pupt == 1)[0])
    # adjust minimum number of iterations
    if nmin < nupt:
        nmin = nupt
    # determine precip. (or max. attr. fraction)
    maxatt = calc_maxatt(qtot, iupt, verbose)
    if maxatt >= 1:
        prec = dqdt[0]
    else:
        prec = maxatt * dqdt[0]
    # location weights?
    maxcon = calc_maxcon(qtot, iupt, verbose)
    if weight_locations:
        # using dqdt with maximum contribution constraint
        lmax = calc_maxcon(qtot, iupt, verbose=True)
        lweights = lmax[iupt]/np.sum(lmax[iupt])
    else:
        lweights = np.repeat(1/nupt, nupt)
    ## starting the random attribution loop
    dqdt_random = np.zeros(shape=nt)
    expl = 0
    icount = 0
    while np.round(expl, 8) < np.round(abs(prec), 8):
        # enforce attribution to initial cond. if explain==max
        if icount == 0 and explainp == "max" and not forc_all:
            ii = nt - 1
        # enfore acttribution to all uptake locations (forc_all==True)
        elif icount < nupt and forc_all:
            if icount == 0 and verbose:
                print(
                    "  *** Random attribution with forc_all=True: enforcing at least one attribution to all "
                    + str(nupt)
                    + " uptake locations"
                )
            i = range(nupt)[icount]
            ii = np.where(pupt == 1)[0][i]  # uptake location index
            if veryverbose:
                print("  *** -- enforcing attribution to uptake location " + str(ii))
        else:
            ii = random.choices(iupt, weights=lweights, k=1)
        # determine maximum attribution for current uptake location and iteration
        try:
            imin = np.argmin(qtot[1:ii]) + 1
        except:
            imin = 1
        iatt = qtot[imin] - np.sum(dqdt_random[imin:])
        if iatt < 0:  # quick fix: this can happen due to precision issues...
            iatt = 0
        idqdt_max = min(iatt, dqdt[ii] - dqdt_random[ii])
        # get random value
        rvalue = random.uniform(0, min(idqdt_max, abs(prec) / nmin))
        if (expl + rvalue) > abs(prec):
            rvalue = abs(prec) - expl
            if rvalue < 0:
                print("OHOH")
        expl += rvalue
        dqdt_random[ii] += rvalue
        icount += 1
        # safety exit (e.g. sum of dqdt_max cannot fully explain prec)
        if icount >= 10000 * nmin:
            print(
                " * Stopping at "
                + str(icount)
                + " iterations; attributed {:.2f}".format(
                    100 * np.sum(dqdt_random) / abs(prec)
                )
                + "%."
            )
            break
        # safety check
        if np.any(dqdt_random < 0):
            print("\n !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
            print(" ABORT: negative values in random attribution... ")
            print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! \n")
            sys.exit()
    # reset for maximum attribution (only needed if veryverbose is set to True)
    if explainp == "max":
        dqdt_random[-1] = 0
        nupt -= 1
    if veryverbose:
        print(
            "  *** "
            + str(icount)
            + " Iterations for "
            + str(nupt)
            + " uptake locations with P={:.4f}".format(dqdt[0])
            + " g/kg with E2Prandom={:.4f}".format(np.sum(dqdt_random))
            + " g/kg (attributed {:.2f}".format(
                100 * np.sum(dqdt_random) / abs(dqdt[0])
            )
            + "%)."
        )
    # upscaling to 100% if explain==full
    if explainp == "full" and maxatt < 1:
        explfr = abs(dqdt[0]) / np.sum(dqdt_random)
        dqdt_random *= explfr
        if veryverbose:
            print(" * Upscaling of contributions required...")
            print(
                "  *** "
                + str(icount)
                + " Iterations for "
                + str(nupt)
                + " uptake locations with P={:.4f}".format(dqdt[0])
                + " g/kg with E2Prandom={:.4f}".format(np.sum(dqdt_random))
                + " g/kg (attributed {:.2f}".format(
                    100 * np.sum(dqdt_random) / abs(dqdt[0])
                )
                + "%)."
            )
    return dqdt_random / 1000


def gridder(plon, plat, pval, glat, glon):
    """
    INPUT
        - plon, plat: parcel longitutde and latitude
        - glon, glat: grid longitude and latitutde
        - pval      : parcel value to be assigned to grid
    ACTION
        1. calculated midpoint of two coordinates
        2. assigns val to gridcell corresponding to midpoint
    RETURN
        - array of dimension (glat.size x glon.size) with 0's and one value assigned
    """
    # 1. calculate midpoint
    lat_mid, lon_mid = midpoint_on_sphere2(
        plat[0], plon[0], plat[1], plon[1]
    )  # calculate midpoint position
    if lon_mid > 179.5:
        lon_mid -= 360  # now shift all coords that otherwise would be allocated to +180 deg to - 180
    if lon_mid < -180.5:
        lon_mid += (
            360  # same for the other direction; only correct for 1deg grid (as below)!
        )
    # 2. get grid index
    ind_lat = np.argmin(np.abs(glat - lat_mid))  # index on grid
    ind_lon = np.argmin(np.abs(glon - lon_mid))  # index on grid
    # and assign pval to gridcell (init. with 0's)
    gval = np.zeros(
        shape=(glat.size, glon.size)
    )  # shape acc. to pre-allocated result array of dim (ntime, glat.size, glon.size)
    gval[ind_lat, ind_lon] += pval
    return gval


def gridall(lon, lat, val, glon, glat):
    # create pandas data frame to sum up over each grid cell
    mydf = pd.DataFrame({"lon": lon, "lat": lat, "v": val})
    uniq = mydf.groupby(["lon", "lat"]).sum()
    # extract values from pandas df (all 1D vectors; unique x,y combin.)
    v = np.asarray(uniq["v"])
    x = np.asarray(uniq.index.get_level_values(0))
    y = np.asarray(uniq.index.get_level_values(1))
    # write into field
    gvalues = np.zeros(shape=(glat.size, glon.size))
    gvalues[y, x] = v
    return gvalues


def midpindex(parray, glon, glat):
    lats = parray[:, 2]  # latitude
    lons = parray[:, 1]  # longitude
    lons[
        lons > 180.0
    ] -= 360  # transform coordinates from [0 ... 360] to [-180 ... 180]
    mlat, mlon = midpoint_on_sphere2(
        lats[0], lons[0], lats[1], lons[1]
    )  # use own function to calculate midpoint position
    if mlon > 179.5:
        mlon -= 360  # now shift all coords that otherwise would be allocated to +180 deg to - 180
    elif mlon < -180.5:
        mlon += (
            360  # same for the other direction; only correct for 1deg grid (as below)!
        )
    # get grid index
    ind_lat = np.argmin(np.abs(glat - mlat))  # index on grid
    ind_lon = np.argmin(np.abs(glon - mlon))  # index on grid
    return ind_lat, ind_lon


def arrpindex(parray, glon, glat):
    # function to get arrival point index on grid
    lats = parray[2]  # latitude
    lons = parray[1]  # longitude
    if lons > 180.0:
        lons = lons - 360  # transform coordinates from [0 ... 360] to [-180 ... 180]
    ind_lat = np.argmin(np.abs(glat - lats))
    ind_lon = np.argmin(np.abs(glon - lons))
    return ind_lat, ind_lon


def get_all_midpindices(ary, glon, glat):
    # gets midpoint indices for all parcels from 3D array of dimension nsteps x nparcels x nvars
    lary = [
        y for y in (np.moveaxis(ary, 1, 0))
    ]  # convert to list for first dimension (parcels) to be able to use map
    imidi = np.asarray(list(map(lambda p: midpindex(p, glon=glon, glat=glat), lary)))
    return imidi


def writeemptync(
    ofile,
    fdate_seq,
    glon,
    glat,
    strargs,
    precision,
    fproc_npart,
    fprec,
    fevap,
    fheat,
    currentversion=get_currentversion(),
):
    # delete nc file if it is present (avoiding error message)
    try:
        os.remove(ofile)
    except OSError:
        pass

    # create netCDF4 instance
    nc_f = nc4.Dataset(ofile, "w", format="NETCDF4")

    ### create dimensions ###
    nc_f.createDimension("time", len(fdate_seq))
    nc_f.createDimension("lat", glat.size)
    nc_f.createDimension("lon", glon.size)

    # create variables
    times = nc_f.createVariable("time", "f8", "time")
    latitudes = nc_f.createVariable("lat", "f8", "lat")
    longitudes = nc_f.createVariable("lon", "f8", "lon")
    if fheat:
        heats = nc_f.createVariable(
            "H",
            precision,
            ("time", "lat", "lon"),
            fill_value=nc4.default_fillvals[precision],
        )
        hnparts = nc_f.createVariable(
            "H_n_part",
            "i4",
            ("time", "lat", "lon"),
            fill_value=nc4.default_fillvals["i4"],
        )
    if fevap:
        evaps = nc_f.createVariable(
            "E",
            precision,
            ("time", "lat", "lon"),
            fill_value=nc4.default_fillvals[precision],
        )
        enparts = nc_f.createVariable(
            "E_n_part",
            "i4",
            ("time", "lat", "lon"),
            fill_value=nc4.default_fillvals["i4"],
        )
    if fprec:
        precs = nc_f.createVariable(
            "P",
            precision,
            ("time", "lat", "lon"),
            fill_value=nc4.default_fillvals[precision],
        )
        pnparts = nc_f.createVariable(
            "P_n_part",
            "i4",
            ("time", "lat", "lon"),
            fill_value=nc4.default_fillvals["i4"],
        )
    if fproc_npart:
        nparts = nc_f.createVariable(
            "n_part",
            "i4",
            ("time", "lat", "lon"),
            fill_value=nc4.default_fillvals["i4"],
        )

    # set attributes
    nc_f.title = "Diagnosis (01) of FLEXPART fluxes"
    nc_f.description = "01_diagnosis - " + str(strargs)
    today = datetime.datetime.now()
    nc_f.history = "Created " + today.strftime("%d/%m/%Y %H:%M:%S") + " using HAMSTER."
    nc_f.institution = (
        "Hydro-Climate Extremes Laboratory (H-CEL), Ghent University, Ghent, Belgium"
    )
    nc_f.source = (
        "HAMSTER "
        + str(currentversion)
        + " ((c) Dominik Schumacher and Jessica Keune): https://github.ugent.be/jkeune/hamster"
    )
    times.units = "hours since 1900-01-01 00:00:00"
    times.calendar = "Standard"  # do NOT use gregorian here!
    latitudes.units = "degrees_north"
    longitudes.units = "degrees_east"
    if fheat:
        heats.units = "W m-2"
        heats.long_name = "surface sensible heat flux"
        hnparts.units = "int"
        hnparts.long_name = "number of H parcels (mid pos.)"
    if fevap:
        evaps.units = "mm"
        evaps.long_name = "evaporation"
        enparts.units = "int"
        enparts.long_name = "number of E parcels (mid pos.)"
    if fprec:
        precs.units = "mm"
        precs.long_name = "precipitation"
        pnparts.units = "int"
        pnparts.long_name = "number of P parcels (mid pos.)"
    if fproc_npart:
        nparts.units = "int"
        nparts.long_name = "number of parcels (mid pos.)"

    # write data
    times[:] = nc4.date2num(fdate_seq, times.units, times.calendar)
    longitudes[:] = glon
    latitudes[:] = glat
    # close file
    nc_f.close()
    print(
        "\n * Created empty file: "
        + ofile
        + " of dimension ("
        + str(len(fdate_seq))
        + ","
        + str(glat.size)
        + ","
        + str(glon.size)
        + ") !"
    )


def writenc(ofile, ix, iarray, ivar, verbose=True):
    if verbose:
        print(" * Writing " + ivar + " to netcdf...")
    nc_f = nc4.Dataset(ofile, "r+")
    nc_f[ivar][ix, :, :] = iarray
    nc_f.close()


def writeemptync4D(
    ofile, fdate_seq, upt_days, glat, glon, strargs, precision, currentversion=get_currentversion()
):

    # delete nc file if it is present (avoiding error message)
    try:
        os.remove(ofile)
    except OSError:
        pass

    # create netCDF4 instance
    nc_f = nc4.Dataset(ofile, "w", format="NETCDF4")

    # create dimensions
    nc_f.createDimension("time", len(fdate_seq))
    nc_f.createDimension("level", upt_days.size)  # could use len() too
    nc_f.createDimension("lat", glat.size)
    nc_f.createDimension("lon", glon.size)

    # create variables
    atimes = nc_f.createVariable("time", "f8", "time")
    utimes = nc_f.createVariable("level", "i4", "level")
    latitudes = nc_f.createVariable("lat", "f8", "lat")
    longitudes = nc_f.createVariable("lon", "f8", "lon")
    heats = nc_f.createVariable(
        "Had",
        precision,
        ("time", "level", "lat", "lon"),
        fill_value=nc4.default_fillvals[precision],
    )
    etops = nc_f.createVariable(
        "E2P",
        precision,
        ("time", "level", "lat", "lon"),
        fill_value=nc4.default_fillvals[precision],
    )

    # set attributes
    nc_f.title = "Attribution (02) of sources using FLEXPART output"
    nc_f.description = "02_attribution - " + str(strargs)
    today = datetime.datetime.now()
    nc_f.history = "Created " + today.strftime("%d/%m/%Y %H:%M:%S") + " using HAMSTER."
    nc_f.institution = (
        "Hydro-Climate Extremes Laboratory (H-CEL), Ghent University, Ghent, Belgium"
    )
    nc_f.source = (
        "HAMSTER "
        + str(currentversion)
        + " ((c) Dominik Schumacher and Jessica Keune): https://github.ugent.be/jkeune/hamster"
    )
    atimes.units = "hours since 1900-01-01 00:00:00"
    atimes.calendar = "Standard"
    utimes.long_name = "Difference between uptake and arrival time, in days"
    utimes.units = "day"
    latitudes.units = "degrees_north"
    longitudes.units = "degrees_east"
    heats.units = "W m-2"
    heats.long_name = "surface sensible heat flux"
    etops.units = "mm"
    etops.long_name = "evaporation resulting in precipitation"

    # write data
    atimes[:] = nc4.date2num(fdate_seq, atimes.units, atimes.calendar)[:]
    utimes[:] = upt_days[:]
    longitudes[:] = glon[:]
    latitudes[:] = glat[:]

    # close file
    nc_f.close()

    print(
        "\n * Created empty file: "
        + ofile
        + " of dimension ("
        + str(len(fdate_seq))
        + ","
        + str(upt_days.size)
        + ","
        + str(glat.size)
        + ","
        + str(glon.size)
        + ") !"
    )


def writenc4D(ofile, ix, ary_etop, ary_heat, verbose):
    if verbose:
        print(" * Writing to netcdf...")

    nc_f = nc4.Dataset(ofile, "r+")
    nc_f["E2P"][ix, :, :, :] = ary_etop[:, :, :]
    nc_f["Had"][ix, :, :, :] = ary_heat[:, :, :]
    nc_f.close()


def get_reference_data(refpath, ivar, idates):
    # ATTENTION: this function assumes that
    # (i)   all data is already on the correct grid (resolution, orientation!)
    # (ii)  all data comes at daily (or subdaily) time steps; in case it's subdaily, these are *summed up* to daily (likely works for E and P, but not for H!)
    # (iii) there are no duplicate files and dates in the refpath folder
    # (iv)  the units HAVE to be as follows: E [mm]; P [mm]; H [W m-2] !

    seldates = np.asarray([datetime.date(rt.year, rt.month, rt.day) for rt in idates])
    selyears = np.unique([dft.strftime("%Y") for dft in idates])

    # get files from refpath that have selyears in the filename
    ifiles = [
        fnmatch.filter(os.listdir(refpath), "*" + str(iyear) + "*.nc")
        for iyear in selyears
    ]

    # read grid
    with nc4.Dataset(os.path.join(refpath, str(ifiles[0][0])), mode="r") as f:
        reflats = np.asarray(f["lat"][:])
        reflons = np.asarray(f["lon"][:])

    # initialize empty array of correct size (seldates x latitude x longitude)
    refdata = np.zeros(shape=(len(seldates), len(reflats), len(reflons)))
    dcounter = np.zeros(shape=(len(seldates)))

    for ifile in ifiles:
        reffile = os.path.join(refpath, ifile[0])
        print("     * Reading " + str(reffile))
        with nc4.Dataset(reffile, mode="r") as f:
            reftime = nc4.num2date(f["time"][:], f["time"].units, f["time"].calendar)
            refdates = np.asarray(
                [datetime.date(rt.year, rt.month, rt.day) for rt in reftime]
            )
            idata = np.asarray(f[ivar][:, :, :])
            ifiledates = np.intersect1d(refdates, seldates)
        for i in range(len(ifiledates)):
            isel = np.where(seldates == ifiledates[i])[0][0]
            iref = np.where(refdates == ifiledates[i])[0][0]
            dcounter[isel] += 1
            refdata[isel, :, :] += idata[iref, :, :]

    if np.any(dcounter == 0):
        print(
            " * The following dates are missing in the reference data set (for variable "
            + str(ivar)
            + "): "
            + str(seldates[np.where(dcounter == 0)[0]])
        )
        raise SystemExit(
            "---------- FATAL ERROR: DATES MISSING IN THE REFERENCES DATA SET!!!"
        )

    return refdata, reflats, reflons


def eraloader_12hourly(
    var, datapath, maskpos, maskneg, uptake_years, uptake_dates, lats, lons
):
    """
    quickly adjusted to enable multi-annual support at the cost of reading in
    two entire years, instead of just what is needed.
    """

    uyears = np.unique(uptake_years)

    with nc4.Dataset(datapath + str(uyears[0]) + ".nc", mode="r") as f:
        print("     * Reading " + str(datapath + str(uyears[0]) + ".nc"))
        reflats = np.asarray(f["latitude"][:])
        reflons = np.asarray(f["longitude"][:])
        reftime = nc4.num2date(
            f["time"][:], f["time"].units, f["time"].calendar
        ) - timedelta(hours=12)
        refdates = np.asarray(
            [datetime.date(rt.year, rt.month, rt.day) for rt in reftime]
        )
        array = np.asarray(f[var][:, :, :])  # not ideal to load everything..
        units = f[var].units

    for ii in range(1, uyears.size):

        print("     * Reading " + str(datapath + str(uyears[ii]) + ".nc"))
        with nc4.Dataset(datapath + str(uyears[ii]) + ".nc", mode="r") as f:
            reftimeY = nc4.num2date(
                f["time"][:], f["time"].units, f["time"].calendar
            ) - timedelta(hours=12)
            reftime = np.concatenate((reftime, reftimeY))
            refdates = np.concatenate(
                (
                    refdates,
                    np.asarray(
                        [datetime.date(rt.year, rt.month, rt.day) for rt in reftimeY]
                    ),
                )
            )
            array = np.concatenate((array, np.asarray(f[var][:, :, :])), axis=0)

    ## first figure out what is needed! NOTE: this is less elegant now due to multi-annual support.
    jbeg = np.min(np.where(refdates == uptake_dates[0]))
    jend = np.max(np.where(refdates == uptake_dates[-1]))  # 12-hourly input here!

    array = array[jbeg : jend + 1, :, :]
    reftime = reftime[jbeg : jend + 1]
    refdates = refdates[jbeg : jend + 1]

    ## mask positive values (E: condensation, H: downward fluxes, P: do NOT mask)
    if maskpos:
        array[array > 0] = 0
    ## mask negative values (P only!)
    if maskneg:
        array[array < 0] = 0

    ## aggregate to daily
    refudates = np.unique(refdates)

    daily = np.empty(shape=(refudates.size, reflats.size, reflons.size))
    for i in range(refudates.size):
        rud = refudates[i]
        daily[i, :, :] = np.nansum(array[np.where(refdates == rud)[0], :, :], axis=0)
    if not np.array_equal(refudates, uptake_dates):
        raise SystemExit("---- no good")

    ## regrid (flip LATS; LON 0 .. 359 to -180 .. 179)
    reflats = np.flipud(reflats)
    daily = np.flip(daily, axis=1)
    # -- above: flipping latitudes, below: from 0 .. 359 to -180 .. 179
    daily_bu = np.copy(daily)
    reflons_bu = np.copy(reflons)
    reflons[: int(reflons.size / 2)] = reflons_bu[int(reflons.size / 2) :] - 360
    reflons[int(reflons.size / 2) :] = reflons_bu[: int(reflons.size / 2)]
    daily[:, :, : int(reflons.size / 2)] = daily_bu[:, :, int(reflons.size / 2) :]
    daily[:, :, int(reflons.size / 2) :] = daily_bu[:, :, : int(reflons.size / 2)]

    # check
    if not (np.array_equal(lats, reflats) or not np.array_equal(lons, reflons)):
        raise SystemExit(
            "---------- FATAL ERROR: regridded ERA-I coordinates don't match!"
        )

    ## units... and sign
    if var == "e" and units == "m of water equivalent":
        daily *= -1e3  # flip sign
    elif var == "tp" and units == "m":
        daily *= 1e3
    elif var == "sshf" and units == "J m**-2":
        daily /= -86400  # flip sign
    else:
        raise SystemExit("---- aborted: no can do.")

    return daily, reflats, reflons


def checkdim(var):
    # check dimension of variables (has to be consistent) and use 2D, 3D or 4D definitions
    ndims = len(var.shape)
    if ndims == 4:
        mydims = ("time", "level", "lat", "lon")
    if ndims == 3:
        mydims = ("time", "lat", "lon")
    if ndims == 2:
        mydims = ("lat", "lon")
    return mydims


def writefinalnc(
    ofile,
    fdate_seq,
    udate_seq,
    glon,
    glat,
    Had,
    Had_Hs,
    E2P,
    E2P_Es,
    E2P_Ps,
    E2P_EPs,
    T2P_EPs,
    strargs,
    precision,
    fwrite_month,
    fbc_had,
    fbc_had_h,
    fbc_e2p,
    fbc_e2p_p,
    fbc_e2p_e,
    fbc_e2p_ep,
    fbc_t2p_ep,
    currentversion=get_currentversion(),
):

    # delete nc file if it is present (avoiding error message)
    try:
        os.remove(ofile)
    except OSError:
        pass

    # create netCDF4 instance
    nc_f = nc4.Dataset(ofile, "w", format="NETCDF4")

    ### create dimensions ###
    if fwrite_month:
        nc_f.createDimension("time", 1)
    else:
        nc_f.createDimension("time", len(fdate_seq))
    if not np.any(np.isnan(udate_seq)):
        nc_f.createDimension("level", len(udate_seq))
    nc_f.createDimension("lat", glat.size)
    nc_f.createDimension("lon", glon.size)

    # create grid + time variables
    times = nc_f.createVariable("time", "f8", "time")
    if not np.any(np.isnan(udate_seq)):
        utimes = nc_f.createVariable("level", "i4", "level")
    latitudes = nc_f.createVariable("lat", "f8", "lat")
    longitudes = nc_f.createVariable("lon", "f8", "lon")

    # create variables
    if fbc_had:
        heats = nc_f.createVariable(
            "Had", precision, checkdim(Had), fill_value=nc4.default_fillvals[precision]
        )
    if fbc_had_h:
        heats_Hs = nc_f.createVariable(
            "Had_Hs",
            precision,
            checkdim(Had_Hs),
            fill_value=nc4.default_fillvals[precision],
        )
    if fbc_e2p:
        evaps = nc_f.createVariable(
            "E2P", precision, checkdim(E2P), fill_value=nc4.default_fillvals[precision]
        )
    if fbc_e2p_e:
        evaps_Es = nc_f.createVariable(
            "E2P_Es",
            precision,
            checkdim(E2P_Es),
            fill_value=nc4.default_fillvals[precision],
        )
    if fbc_e2p_p:
        evaps_Ps = nc_f.createVariable(
            "E2P_Ps",
            precision,
            checkdim(E2P_Ps),
            fill_value=nc4.default_fillvals[precision],
        )
    if fbc_e2p_ep:
        evaps_EPs = nc_f.createVariable(
            "E2P_EPs",
            precision,
            checkdim(E2P_EPs),
            fill_value=nc4.default_fillvals[precision],
        )
    if fbc_t2p_ep:
        trans_EPs = nc_f.createVariable(
            "T2P_EPs",
            precision,
            checkdim(T2P_EPs),
            fill_value=nc4.default_fillvals[precision],
        )

    # set attributes
    nc_f.title = "Bias-corrected source-sink relationships from FLEXPART"
    nc_f.description = str(strargs)
    today = datetime.datetime.now()
    nc_f.history = "Created " + today.strftime("%d/%m/%Y %H:%M:%S") + " using HAMSTER."
    nc_f.institution = (
        "Hydro-Climate Extremes Laboratory (H-CEL), Ghent University, Ghent, Belgium"
    )
    nc_f.source = (
        "HAMSTER "
        + str(currentversion)
        + " ((c) Dominik Schumacher and Jessica Keune): https://github.ugent.be/jkeune/hamster"
    )
    times.units = "hours since 1900-01-01 00:00:00"
    times.calendar = "Standard"  # do NOT use gregorian here!
    if not np.any(np.isnan(udate_seq)):
        utimes.long_name = "Difference between uptake and arrival time, in days"
        utimes.units = "day"
    latitudes.units = "degrees_north"
    longitudes.units = "degrees_east"
    if fbc_had:
        heats.units = "W m-2"
        heats.long_name = "advected surface sensible heat"
    if fbc_had_h:
        heats_Hs.units = "W m-2"
        heats_Hs.long_name = (
            "advected surface sensible heat, H-scaled"  # this is garbage, I know
        )
    if fbc_e2p:
        evaps.units = "mm"
        evaps.long_name = "evaporation resulting in precipitation"
    if fbc_e2p_e:
        evaps_Es.units = "mm"
        evaps_Es.long_name = "evaporation resulting in precipitation, E-corrected"
    if fbc_e2p_p:
        evaps_Ps.units = "mm"
        evaps_Ps.long_name = "evaporation resulting in precipitation, P-corrected"
    if fbc_e2p_ep:
        evaps_EPs.units = "mm"
        evaps_EPs.long_name = (
            "evaporation resulting in precipitation, E-and-P-corrected"
        )
    if fbc_t2p_ep:
        trans_EPs.units = "mm"
        trans_EPs.long_name = (
            "transpiration resulting in precipitation, E-and-P-corrected"
        )

    # write data
    if fwrite_month:
        times[:] = nc4.date2num(fdate_seq[0], times.units, times.calendar)
    else:
        times[:] = nc4.date2num(fdate_seq, times.units, times.calendar)
    if not np.any(np.isnan(udate_seq)):
        utimes[:] = np.arange(-len(udate_seq) + 1, 1)
    latitudes[:] = glat
    longitudes[:] = glon

    if fwrite_month:
        if fbc_had:
            heats[:] = np.nanmean(Had, axis=0, keepdims=True)[:]
        if fbc_had_h:
            heats_Hs[:] = np.nanmean(Had_Hs, axis=0, keepdims=True)[:]
        if fbc_e2p:
            evaps[:] = np.nansum(E2P, axis=0, keepdims=True)[:]
        if fbc_e2p_e:
            evaps_Es[:] = np.nansum(E2P_Es, axis=0, keepdims=True)[:]
        if fbc_e2p_p:
            evaps_Ps[:] = np.nansum(E2P_Ps, axis=0, keepdims=True)[:]
        if fbc_e2p_ep:
            evaps_EPs[:] = np.nansum(E2P_EPs, axis=0, keepdims=True)[:]
        if fbc_t2p_ep:
            trans_EPs[:] = np.nansum(T2P_EPs, axis=0, keepdims=True)[:]
    else:
        if fbc_had:
            heats[:] = Had[:]
        if fbc_had_h:
            heats_Hs[:] = Had_Hs[:]
        if fbc_e2p:
            evaps[:] = E2P[:]
        if fbc_e2p_e:
            evaps_Es[:] = E2P_Es[:]
        if fbc_e2p_p:
            evaps_Ps[:] = E2P_Ps[:]
        if fbc_e2p_ep:
            evaps_EPs[:] = E2P_EPs[:]
        if fbc_t2p_ep:
            trans_EPs[:] = T2P_EPs[:]

    # close file
    nc_f.close()

    # print info
    print("\n * Created and wrote to file: " + ofile + " !\n")


def append2csv(filename, listvals):
    # Open file in append mode
    with open(filename, "a+", newline="\n") as write_obj:
        csv_writer = csv.writer(write_obj, delimiter="\t", lineterminator="\n")
        csv_writer.writerow(listvals)


def convert2daily(xar, ftime, fagg="mean"):

    if ftime[0].hour in [0, 6, 12, 18]:
        ## simple fix, subtract 3 hours
        ftime = np.asarray([t - datetime.timedelta(hours=3) for t in ftime])
        dates = cal2date(ftime)
    elif ftime[0].hour in [3, 9, 15, 21]:
        ## NOTE: this is the new norm! retain "old style" for now, though
        dates = cal2date(ftime)
        pass
    dtime = np.unique(dates)

    xtot = np.zeros(shape=(dtime.size, xar.shape[1], xar.shape[2]))

    ## this isn't fast or elegant, but works for literally anything sub-daily
    for i in range(dtime.size):

        iud = dtime[i]
        sel = np.where(dates == iud)[0]

        ## TODO: clean up; there should be a check whether 4 files are present, imo
        if sel.size != 4:
            warnings.warn(
                "\n\n----------------- WARNING: this should NEVER OCCUR; daily aggregation IMPROPER (files missing!)\n\n"
            )

        if fagg == "sum":
            xtot[i, :, :] = np.nansum(xar[sel, :, :], axis=0)
        if fagg == "mean":
            xtot[i, :, :] = np.nanmean(xar[sel, :, :], axis=0)

    return xtot


def read_diagdata(opathD, ofile_base, ryyyy, uptake_time, var="E"):
    # get required months
    ayears = np.asarray([ut.year for ut in uptake_time])
    amonths = np.asarray([ut.month for ut in uptake_time])
    uset = np.unique(np.column_stack((ayears, amonths)), axis=0)
    uyears, umonths = uset[:, 0], uset[:, 1]

    # loop over months
    for jj in range(umonths.size):
        uyr = str(uyears[jj])
        umon = str(umonths[jj])
        diagfile = (
            str(opathD)
            + "/"
            + str(ofile_base)
            + "_diag_r"
            + str(ryyyy)[-2:]
            + "_"
            + str(uyr)
            + "-"
            + umon.zfill(2)
            + ".nc"
        )
        with nc4.Dataset(diagfile, mode="r") as f:
            if var not in ["grid"]:
                ix = f[var][:]
                timex = nc4.num2date(f["time"][:], f["time"].units, f["time"].calendar)

        ## concatenate 'em!
        if var == "grid":
            with nc4.Dataset(diagfile, mode="r") as f:
                lats = f["lat"][:]
                lons = f["lon"][:]
        else:
            # concatenate
            if umon == str(umonths[0]):
                x = np.copy(ix)
                ftime = np.copy(timex)
            else:
                x = np.concatenate((x, ix), axis=0)
                ftime = np.concatenate((ftime, timex))

    # return
    if var == "time":
        return ftime
    if var == "grid":
        return (lats, lons)
    if var not in ["grid", "time"]:
        return x


def gridcheck(lats, totlats, lons, totlons):
    if not np.array_equal(lats, totlats) or not np.array_equal(lons, totlons):
        raise SystemExit("--- ERROR: your grids aren't identical...")


def datecheck(idate, dateseq):
    if idate not in dateseq:
        raise SystemExit(
            "\n !!! ERROR: INPUT DATA MISSING: date "
            + str(idate)
            + " not available as output from 01_diagnosis! Aborting here. !!!\n"
        )


def calc_alpha(top, bot):
    alpha = np.divide(top, bot)
    return alpha


def calc_sourcebcf(ref, diag, tscale="daily"):
    # define all positive
    if np.all(ref <= 0):
        ref = -ref
    if np.all(diag <= 0):
        diag = -diag
    # set 0 to nan to avoid 1e300 values
    diag[diag == 0] = np.nan
    # calculate bias correction factor
    if tscale == "daily":
        alpha = np.nan_to_num(np.divide(ref, diag))
    if tscale == "monthly":
        ref_sum = np.nansum(ref, axis=(0))
        diag_sum = np.nansum(diag, axis=(0))
        alpha = np.nan_to_num(np.divide(ref, diag))
        # alpha   = np.repeat(alpha_H[np.newaxis,:,:,:], 28, axis=0)
    alpha[alpha == np.inf] = 0
    return alpha


def udays2udate(atime, utime_srt):
    utime_first = atime[0] - timedelta(
        days=utime_srt.size - 1
    )  # utime_srt.size-1 == trajlen (in days)
    uptake_time = np.asarray(
        [
            utime_first + timedelta(days=nday)
            for nday in range(utime_srt.size - 1 + atime.size)
        ]
    )
    return uptake_time


def expand4Darray(myarray, atime, utime_srt, veryverbose):
    utime = udays2udate(atime, utime_srt)
    myshape = myarray.shape
    myarray_exp = np.empty(shape=(myshape[0], utime.size, myshape[2], myshape[3]))
    if veryverbose:
        print(
            " * Expanding array from " + str(myshape) + " to " + str(myarray_exp.shape)
        )
    for iat in range(atime.size):
        myarray_exp[iat, iat : iat + utime_srt.size, :, :] = myarray[iat, :, :, :]
    return myarray_exp


def reduce4Darray(myarray, veryverbose):
    myshape = myarray.shape
    bwtimesteps = myshape[1] - myshape[0] + 1
    myarray_red = np.empty(shape=(myshape[0], bwtimesteps, myshape[2], myshape[3]))
    if veryverbose:
        print(
            " * Reducing array from " + str(myshape) + " to " + str(myarray_red.shape)
        )
    for iat in range(myshape[0]):
        myarray_red[iat, :, :, :] = myarray[iat, (iat) : (iat + bwtimesteps), :, :]
    return myarray_red


def date2year(mydates):
    return np.asarray([it.year for it in mydates])


def date2month(mydates):
    return np.asarray([it.month for it in mydates])


def cal2date(mydates):
    return np.asarray([datetime.date(it.year, it.month, it.day) for it in mydates])


def convert_mm_m3(myarray, areas):
    # we ALWAYS follow the array dimensions order: (anything(s) x lat x lon) here
    if len(areas.shape) > 1:
        # attention: broadcasting with a 2D array in numpy can lead to differences
        # 1e-9 (f4) or 1e-17 (f8)
        carray = np.multiply(areas, myarray / 1e3)
    if len(areas.shape) == 1:
        ## swap axes to enable numpy broadcasting;
        ## (a,b,c,d x b,c,d OK; a,b,c,d x a NOT OK)
        ldim = len(myarray.shape) - 1
        carray = np.swapaxes(
            areas * np.moveaxis(myarray / 1e3, ldim - 1, ldim), ldim - 1, ldim
        )
    return carray


def convert_m3_mm(myarray, areas):
    # we ALWAYS follow the array dimensions order: (anything(s) x lat x lon) here
    if len(areas.shape) > 1:
        # attention: broadcasting with a 2D array in numpy can lead to differences
        # 1e-9 (f4) or 1e-17 (f8)
        carray = np.nan_to_num(np.divide(myarray * 1e3, areas))
    if len(areas.shape) == 1:
        ## swap axes to enable numpy broadcasting;
        ## (a,b,c,d x b,c,d OK; a,b,c,d x a NOT OK)
        ldim = len(myarray.shape) - 1
        carray = np.swapaxes(
            np.nan_to_num(np.divide(np.moveaxis(myarray * 1e3, ldim - 1, ldim), areas)),
            ldim - 1,
            ldim,
        )
    return carray


def check_attributedp(pdiag, pattr, veryverbose):
    # define all positive
    if np.all(pdiag <= 0):
        pdiag = -pdiag
    if np.all(pattr <= 0):
        pattr = -pattr
    printwarning = False
    returnval = False
    pdiag_sum = np.nansum(pdiag, axis=(1))
    pattr_sum = np.nansum(pattr[:, :, :, :], axis=(1, 2, 3))
    if round(np.nansum(pdiag_sum), 4) != round(np.nansum(pattr_sum), 4):
        print(
            "   --- WARNING: total precipitation from 01_diagnosis and 02_attribution differ"
        )
        print(
            " \t --- Absolute difference: {:.2f}".format(
                np.nansum(pdiag_sum) - np.nansum(pattr_sum)
            )
            + " m3"
        )
        printwarning = True
    if np.any(pdiag_sum - pattr_sum != 0):
        print(
            "   --- WARNING: daily precipitation from 01_diagnosis and 02_attribution differ"
        )
        if veryverbose:
            ndiffs = len(np.where(pdiag_sum - pattr_sum != 0)[0])
            print(
                " \t --- " + str(ndiffs) + " days have different precipitation sums. "
            )
            ndiffs = len(np.where(pdiag_sum > pattr_sum)[0])
            print(
                " \t --- "
                + str(ndiffs)
                + " days have P(01_diagnosis) > P(02_attribution)"
            )
            ndiffs = len(np.where(pdiag_sum < pattr_sum)[0])
            print(
                " \t --- "
                + str(ndiffs)
                + " days have P(01_diagnosis) < P(02_attribution)"
            )
        printwarning = True
    if printwarning:
        print(
            "   --- ATTENTION: Using 02_attribution data for bias correction for consistency."
        )
        returnval = True
    # print(pattr_sum)
    # print(np.nansum(pattr_sum))
    # print(pdiag_sum)
    # print(np.nansum(pdiag_sum))
    # print(pdiag_sum-pattr_sum)
    # print(np.nan_to_num(100*(pdiag_sum-pattr_sum)/pattr_sum))
    return returnval


def needmonthlyp(pdiag, pref):
    returnval = False
    # define all positive
    if np.all(pref <= 0):
        pref = -pref
    if np.all(pdiag <= 0):
        pdiag = -pdiag
    tocheck = np.where(pref > 0)[0]
    if np.any(pdiag[tocheck] == 0):
        print("   --- WARNING: daily bias correction of precipitation not possible.")
        ndiffs = len(np.where(pdiag[tocheck] == 0)[0])
        print(
            " \t --- "
            + str(ndiffs)
            + " days have no diagnosed but observed precipitation."
        )
        print(
            "   --- ATTENTION: Using monthly precipitation for bias correction for consistency."
        )
        returnval = True
    return returnval


def writewarning(wfile):
    with open(wfile, "w") as ifile:
        writer = csv.writer(
            ifile,
            delimiter="\t",
            lineterminator="\n",
            quoting=csv.QUOTE_NONE,
            quotechar="",
        )
        writer.writerow(
            [
                "WARNING: you're writing out daily data, but (additional) monthly bias correction was performed. Your daily data is thus not representative."
            ]
        )
    print("\n WARNING! \n See: " + wfile + " ! \n")


def writestats_03(
    sfile,
    Pref,
    P_E2P,
    P_E2P_Escaled,
    P_E2P_Pscaled,
    P_E2P_EPscaled,
    Had,
    Had_scaled,
    xla,
    xlo,
    ibgn,
):
    with open(sfile, "w") as ifile:
        writer = csv.writer(
            ifile,
            delimiter="\t",
            lineterminator="\n",
            quoting=csv.QUOTE_NONE,
            quotechar="",
        )
        writer.writerow(["* - STATISTICS: "])
        ndays = Pref[ibgn:, :, :].shape[0]
        writer.writerow(["   --- # DAYS EVALUATED:              {:.0f}".format(ndays)])
        writer.writerow([" "])
        writer.writerow(["* - HEAT ADVECTION STATISTICS: "])
        writer.writerow(
            [
                "   --- Had_unscaled [W m-2]:          {:.2f}".format(
                    np.nanmean(np.nansum(Had, axis=(1, 2)))
                )
            ]
        )
        writer.writerow(
            [
                "   --- Had_Hscaled [W m-2]:           {:.2f}".format(
                    np.nanmean(np.nansum(Had_scaled, axis=(1, 2)))
                )
            ]
        )
        writer.writerow([" "])
        writer.writerow(["* - PRECIPITATION STATISTICS: "])
        writer.writerow(
            [
                "   --- P_REFERENCE [m3]:              {:.2f}".format(
                    np.nansum(Pref[ibgn:, xla, xlo])
                )
            ]
        )
        writer.writerow(
            ["   --- P_E2P_unscaled [m3]:           {:.2f}".format(np.nansum(P_E2P))]
        )
        writer.writerow(
            [
                "   --- P_E2P_Escaled [m3]:            {:.2f}".format(
                    np.nansum(P_E2P_Escaled)
                )
            ]
        )
        writer.writerow(
            [
                "   --- P_E2P_Pscaled [m3]:            {:.2f}".format(
                    np.nansum(P_E2P_Pscaled)
                )
            ]
        )
        writer.writerow(
            [
                "   --- P_E2P_EPscaled [m3]:           {:.2f}".format(
                    np.nansum(P_E2P_EPscaled)
                )
            ]
        )
        # some contingency table statistics...
        writer.writerow([" "])
        writer.writerow(["* - CONTINGENCY TABLE SCORES (PRECIPITATION):"])
        pref_sum = np.nansum(Pref[ibgn:, xla, xlo], axis=(1))
        pdiag_sum = np.nansum(P_E2P, axis=(1, 2))
        myctab = contingency_table(pref_sum, pdiag_sum, thresh=0)
        myscores = calc_ctab_measures(myctab)
        writer.writerow(
            ["   --- DAYS OF FALSE ALARMS:        {:.0f}".format(myctab["b"])]
        )
        writer.writerow(
            ["   --- DAYS OF MISSES:              {:.0f}".format(myctab["c"])]
        )
        writer.writerow(
            ["   --- DAYS OF HITS:                {:.0f}".format(myctab["a"])]
        )
        writer.writerow(
            ["   --- DAYS OF CORRECT NEGATIVES:   {:.0f}".format(myctab["d"])]
        )
        writer.writerow(
            ["   --- SUCCESS RATIO:               {:.2f}".format(myscores["sr"])]
        )
        writer.writerow(
            ["   --- FALSE ALARM RATIO:           {:.2f}".format(myscores["far"])]
        )
        writer.writerow(
            ["   --- FREQUENCY BIAS:              {:.2f}".format(myscores["fbias"])]
        )
        writer.writerow(
            ["   --- PROB. OF DETECTION:          {:.2f}".format(myscores["pod"])]
        )
        writer.writerow(
            ["   --- PROB. OF FALSE DETECTION:    {:.2f}".format(myscores["pofd"])]
        )
        writer.writerow(
            ["   --- PEIRCE'S SKILL SCORE:        {:.2f}".format(myscores["pss"])]
        )


def contingency_table(ref, mod, thresh=0):
    # creates a contingency table based on 1D np.arrays
    ieventobs = np.where(ref > thresh)[0]
    ineventobs = np.where(ref <= thresh)[0]
    a = len(np.where(mod[ieventobs] > thresh)[0])  # hits
    b = len(np.where(mod[ineventobs] > thresh)[0])  # false alarms
    c = len(np.where(mod[ieventobs] <= thresh)[0])  # misses
    d = len(np.where(mod[ineventobs] <= thresh)[0])  # correct negatives
    return {"a": a, "b": b, "c": c, "d": d}


def try_div(x, y):
    try:
        return x / y
    except ZeroDivisionError:
        return 0


def calc_ctab_measures(cdict):
    # calculates common contingency table scores
    # scores following definitions from https://www.cawcr.gov.au/projects/verification/
    a = cdict["a"]  # hits
    b = cdict["b"]  # false alarms
    c = cdict["c"]  # misses
    d = cdict["d"]  # correct negatives
    # calculate scores
    acc = try_div(a + d, a + b + c + d)  # accuracy
    far = try_div(b, a + b)  # false alarm ratio
    fbias = try_div(a + b, a + c)  # frequency bias
    pod = try_div(a, a + c)  # probability of detection (hit rate)
    pofd = try_div(b, b + d)  # probability of false detection (false alarm rate)
    sr = try_div(a, a + b)  # success ratio
    ts = try_div(a, a + c + b)  # threat score (critical success index)
    a_random = try_div((a + c) * (a + b), a + b + c + d)
    ets = try_div(
        (a - a_random), (a + b + c + a_random)
    )  # equitable threat score (gilbert skill score)
    pss = pod - pofd  # peirce's skill score (true skill statistic)
    odr = try_div(a * d, c * b)  # odd's ratio
    return {
        "acc": acc,
        "far": far,
        "fbias": fbias,
        "pod": pod,
        "pofd": pofd,
        "sr": sr,
        "pss": pss,
        "odr": odr,
    }


def writestats_02(
    statfile,
    tneval,
    tnnevala,
    tnevalh,
    tnnevalh,
    tnnevalm,
    tnevalp,
    tnnevalp,
    patt,
    psum,
    punatt,
    pmiss,
):
    with open(statfile, "w") as sfile:
        writer = csv.writer(
            sfile,
            delimiter="\t",
            lineterminator="\n",
            quoting=csv.QUOTE_NONE,
            quotechar="",
        )
        writer.writerow(["* - PARCEL STATISTICS: "])
        writer.writerow(["   --- TOTAL EVALUATED PARCELS:       ", str(tneval)])
        writer.writerow([" "])
        writer.writerow(
            ["   --- # PARCELS ARRIVING INSIDE MASK:", str(tneval - tnnevala)]
        )
        if tnnevala != tneval:
            writer.writerow(
                [
                    "   --- # PARCELS EVAL. FOR HEAT-ADV:  ",
                    str(tnevalh)
                    + " ({:.2f}".format(100 * tnevalh / (tneval - tnnevala))
                    + "%)",
                ]
            )
        if tnevalh != 0:
            writer.writerow(
                [
                    "   ----- WITHOUT UPTAKES IN THE TRAJ: ",
                    str(tnnevalh)
                    + " ({:.2f}".format(100 * tnnevalh / (tnevalh))
                    + "%)",
                ]
            )
        writer.writerow([" "])
        writer.writerow(
            ["   --- # PARCELS MIDPOINT INSIDE MASK:", str(tneval - tnnevalm)]
        )
        if tnnevalm != tneval:
            writer.writerow(
                [
                    "   --- # PARCELS EVAL. FOR PRECIP:    ",
                    str(tnevalp)
                    + " ({:.2f}".format(100 * tnevalp / (tneval - tnnevalm))
                    + "%)",
                ]
            )
        if tnevalp != 0:
            writer.writerow(
                [
                    "   ----- WITHOUT UPTAKES IN THE TRAJ: ",
                    str(tnnevalp)
                    + " ({:.2f}".format(100 * tnnevalp / (tnevalp))
                    + "%)",
                ]
            )
        writer.writerow([" "])
        if psum != 0:
            writer.writerow([" * - PRECIPITATION STATISTICS: "])
            writer.writerow(
                ["   --- ATTRIBUTED FRACTION:             {:.2f}".format(patt / psum)]
            )
            writer.writerow(
                ["   --- UNATTRIBUTED FRACTION (TRAJEC):  {:.2f}".format(punatt / psum)]
            )
            writer.writerow(
                ["   --- UNATTRIBUTED FRACTION (NO-UPT):  {:.2f}".format(pmiss / psum)]
            )


def mask3darray(xarray, xla, xlo):
    marray = np.zeros(shape=xarray.shape)
    marray[:, xla, xlo] = xarray[:, xla, xlo]
    return marray


def writedebugnc(
    ofile,
    fdate_seq,
    udate_seq,
    glon,
    glat,
    mask,
    Pref,
    Pdiag,
    Pattr,
    Pattr_Es,
    Pattr_Ps,
    Pattr_EPs,
    frac_E2P,
    frac_Had,
    alpha_P,
    alpha_P_Ecorrected,
    alpha_P_res,
    alpha_E,
    alpha_H,
    strargs,
    precision,
    currentversion=get_currentversion(),
):

    Prefsum = np.nansum(Pref, axis=(1, 2))
    Pdiagsum = np.nansum(Pdiag, axis=(1, 2))
    Pattrsum = np.nansum(Pattr, axis=(1, 2))
    Pattrsum_Es = np.nansum(Pattr_Es, axis=(1, 2))
    Pattrsum_Ps = np.nansum(Pattr_Ps, axis=(1, 2))
    Pattrsum_EPs = np.nansum(Pattr_EPs, axis=(1, 2))
    malpha_H = np.max(alpha_H[:, :, :], axis=(1, 2))
    malpha_E = np.max(alpha_E[:, :, :], axis=(1, 2))

    # delete nc file if it is present (avoiding error message)
    try:
        os.remove(ofile)
    except OSError:
        pass

    # create netCDF4 instance
    nc_f = nc4.Dataset(ofile, "w", format="NETCDF4")

    ### create dimensions ###
    nc_f.createDimension("time", len(fdate_seq))
    nc_f.createDimension("uptaketime", len(udate_seq))
    nc_f.createDimension("lat", glat.size)
    nc_f.createDimension("lon", glon.size)

    # create variables
    times = nc_f.createVariable("time", "f8", "time")
    utimes = nc_f.createVariable("uptaketime", "f8", "uptaketime")
    latitudes = nc_f.createVariable("lat", "f8", "lat")
    longitudes = nc_f.createVariable("lon", "f8", "lon")
    # Variables
    nc_mask = nc_f.createVariable(
        "mask", "i4", ("lat", "lon"), fill_value=nc4.default_fillvals["i4"]
    )
    nc_pref = nc_f.createVariable(
        "Pref",
        precision,
        ("time", "lat", "lon"),
        fill_value=nc4.default_fillvals[precision],
    )
    nc_pdiag = nc_f.createVariable(
        "Pdiag",
        precision,
        ("time", "lat", "lon"),
        fill_value=nc4.default_fillvals[precision],
    )
    nc_pattr = nc_f.createVariable(
        "Pattr",
        precision,
        ("time", "lat", "lon"),
        fill_value=nc4.default_fillvals[precision],
    )
    nc_prefs = nc_f.createVariable(
        "Pref_sum", precision, ("time"), fill_value=nc4.default_fillvals[precision]
    )
    nc_pdiags = nc_f.createVariable(
        "Pdiag_sum", precision, ("time"), fill_value=nc4.default_fillvals[precision]
    )
    nc_pattrs = nc_f.createVariable(
        "Pattr_sum", precision, ("time"), fill_value=nc4.default_fillvals[precision]
    )
    nc_pattrs_es = nc_f.createVariable(
        "Pattr_Es_sum", precision, ("time"), fill_value=nc4.default_fillvals[precision]
    )
    nc_pattrs_ps = nc_f.createVariable(
        "Pattr_Ps_sum", precision, ("time"), fill_value=nc4.default_fillvals[precision]
    )
    nc_pattrs_eps = nc_f.createVariable(
        "Pattr_EPs_sum", precision, ("time"), fill_value=nc4.default_fillvals[precision]
    )
    nc_alphap = nc_f.createVariable(
        "alpha_P", precision, ("time"), fill_value=nc4.default_fillvals[precision]
    )
    nc_alphap_ebc = nc_f.createVariable(
        "alpha_P_Ecorrected",
        precision,
        ("time"),
        fill_value=nc4.default_fillvals[precision],
    )
    nc_alphap_res = nc_f.createVariable(
        "alpha_P_res", precision, ("time"), fill_value=nc4.default_fillvals[precision]
    )
    nc_alphae = nc_f.createVariable(
        "alpha_E",
        precision,
        ("uptaketime", "lat", "lon"),
        fill_value=nc4.default_fillvals[precision],
    )
    nc_alphah = nc_f.createVariable(
        "alpha_H",
        precision,
        ("uptaketime", "lat", "lon"),
        fill_value=nc4.default_fillvals[precision],
    )
    nc_frace2p = nc_f.createVariable(
        "frac_E2P",
        precision,
        ("time", "uptaketime", "lat", "lon"),
        fill_value=nc4.default_fillvals[precision],
    )
    nc_frachad = nc_f.createVariable(
        "frac_Had",
        precision,
        ("time", "uptaketime", "lat", "lon"),
        fill_value=nc4.default_fillvals[precision],
    )
    nc_malphae = nc_f.createVariable(
        "max_alpha_E",
        precision,
        ("uptaketime"),
        fill_value=nc4.default_fillvals[precision],
    )
    nc_malphah = nc_f.createVariable(
        "max_alpha_H",
        precision,
        ("uptaketime"),
        fill_value=nc4.default_fillvals[precision],
    )

    # set attributes
    nc_f.title = "Debug-file from 03_biascorrection (HAMSTER)"
    nc_f.description = str(strargs)
    today = datetime.datetime.now()
    nc_f.history = "Created " + today.strftime("%d/%m/%Y %H:%M:%S") + " using HAMSTER."
    nc_f.institution = (
        "Hydro-Climate Extremes Laboratory (H-CEL), Ghent University, Ghent, Belgium"
    )
    nc_f.source = (
        "HAMSTER "
        + str(currentversion)
        + " ((c) Dominik Schumacher and Jessica Keune): https://github.ugent.be/jkeune/hamster"
    )
    times.units = "hours since 1900-01-01 00:00:00"
    times.calendar = "Standard"  # do NOT use gregorian here!
    utimes.units = "hours since 1900-01-01 00:00:00"
    utimes.calendar = "Standard"  # do NOT use gregorian here!
    latitudes.units = "degrees_north"
    longitudes.units = "degrees_east"
    nc_mask.units = "-"
    nc_pref.units = "m3"
    nc_pref.long_name = "reference precipitation"
    nc_pdiag.units = "m3"
    nc_pdiag.long_name = "diagnosed precipitation (01_diag)"
    nc_pattr.units = "m3"
    nc_pattr.long_name = "attributed precipitation (E2P, 02_attr)"
    nc_prefs.units = "m3"
    nc_prefs.long_name = "sum of reference precipitation"
    nc_pdiags.units = "m3"
    nc_pdiags.long_name = "sum of diagnosed precipitation (01_diag)"
    nc_pattrs.units = "m3"
    nc_pattrs.long_name = "sum of attributed precipitation (E2P, 02_attr)"
    nc_pattrs_es.units = "m3"
    nc_pattrs_es.long_name = "sum of attributed precipitation (E2P_Escaled, 02_attr)"
    nc_pattrs_ps.units = "m3"
    nc_pattrs_ps.long_name = "sum of attributed precipitation (E2P_Pscaled, 02_attr)"
    nc_pattrs_eps.units = "m3"
    nc_pattrs_eps.long_name = "sum of attributed precipitation (E2P_EPscaled, 02_attr)"
    nc_alphap.units = "-"
    nc_alphap.long_name = "alpha_P"
    nc_alphap_ebc.units = "-"
    nc_alphap_ebc.long_name = "alpha_P_Ecorrected"
    nc_alphap_res.units = "-"
    nc_alphap_res.long_name = "alpha_P_res"
    nc_alphae.units = "-"
    nc_alphae.long_name = "alpha_E"
    nc_alphah.units = "-"
    nc_alphah.long_name = "alpha_H"
    nc_malphae.units = "-"
    nc_malphae.long_name = "maximum alpha_E"
    nc_malphah.units = "-"
    nc_malphah.long_name = "maximum alpha_H"
    nc_frace2p.units = "-"
    nc_frace2p.long_name = "frac_E2P"
    nc_frachad.units = "-"
    nc_frachad.long_name = "frac_Had"

    # write data
    times[:] = nc4.date2num(fdate_seq, times.units, times.calendar)
    utimes[:] = nc4.date2num(udate_seq, utimes.units, utimes.calendar)
    latitudes[:] = glat
    longitudes[:] = glon
    nc_pref[:] = Pref[:]
    nc_mask[:] = mask[:]
    nc_pdiag[:] = Pdiag[:]
    nc_pattr[:] = Pattr[:]
    nc_prefs[:] = Prefsum[:]
    nc_pdiags[:] = Pdiagsum[:]
    nc_pattrs[:] = Pattrsum[:]
    nc_pattrs_es[:] = Pattrsum_Es[:]
    nc_pattrs_ps[:] = Pattrsum_Ps[:]
    nc_pattrs_eps[:] = Pattrsum_EPs[:]
    nc_alphap[:] = alpha_P[:]
    nc_alphap_ebc[:] = alpha_P_Ecorrected[:]
    nc_alphap_res[:] = alpha_P_res[:]
    nc_alphae[:] = alpha_E[:]
    nc_alphah[:] = alpha_H[:]
    nc_malphae[:] = malpha_E[:]
    nc_malphah[:] = malpha_H[:]
    nc_frace2p[:] = frac_E2P[:]
    nc_frachad[:] = frac_Had[:]

    # close file
    nc_f.close()

    # print info
    print("\n * Created and wrote to file: " + ofile + " !")


def append_attrfrac_netcdf(ofile, attrfrac):
    nc_f = nc4.Dataset(ofile, "r+")
    attrdesc = (
        getattr(nc_f, "description")
        + "; [[STATS]] attributed fraction (P) = "
        + attrfrac
    )
    nc_f.description = attrdesc
    nc_f.close()


def maskbymaskval(mask, maskval):
    mymask = np.copy(mask)
    mymask[np.where(mask != maskval)] = 0
    return mymask


def calc_sinkbcf(ref, att, tscale="daily"):
    # define all positive
    if np.all(ref <= 0):
        ref = -ref
    if np.all(att <= 0):
        att = -att
    if tscale == "daily":
        tref = np.nansum(ref, axis=tuple(range(1, ref.ndim)))
        tatt = np.nansum(att, axis=tuple(range(1, att.ndim)))
        alpha = tref / tatt
        alpha[alpha == np.inf] = 0
        return alpha
    if tscale == "monthly":
        tref = np.nansum(ref)
        tatt = np.nansum(att)
        alpha = np.repeat(tref / tatt, ref.shape[0])
        alpha[alpha == np.inf] = 0
        return alpha


def checkpsum(ref, att, verbose):
    if round(np.nansum(ref), 4) != round(np.nansum(att), 4):
        ident = False
        if verbose:
            print(
                "       !!! Attributed precipitation does not match reference precipitation!"
            )
            print("           * attributed P: " + str(round(np.nansum(ref), 4)) + " m3")
            print("           * reference P:  " + str(round(np.nansum(att), 4)) + " m3")
    else:
        ident = True
    return ident


def consistencycheck(attr, diag, bcscale, debug):
    if bcscale == "daily":
        aggattr = np.nansum(attr, axis=0)
    if bcscale == "monthly":
        aggattr = np.nansum(attr, axis=(0, 1))
        diag = np.nansum(diag, axis=(0))
    # calculate fractions
    frac = np.divide(aggattr, diag)
    # print warnings
    if np.any(frac > 1.0001) or np.any(np.isinf(frac)):
        print(" \n  \t !!! WARNING: attribution exceeds diagnosis !!!")
        print(" \t !!!          ---> CHECK YOUR DATA !!!")
        print(
            " \t !!!          ---> Maximum fraction: "
            + str(np.max(np.nan_to_num(frac)))
        )
        print(
            " \t !!!          ---> Number of exceedances: "
            + str(
                len(frac[np.where(frac > 1.0001)]) + len(frac[np.where(np.isinf(frac))])
            )
            + " \n"
        )
        if debug:
            print(
                " \t !!!          ---> frac > 1.001 at: " + str(np.where(frac > 1.001))
            )
            print(
                " \t !!!          ----------> attr: "
                + str(aggattr[np.where(frac > 1.001)])
            )
            print(
                " \t !!!          ----------> diag: "
                + str((aggattr / frac)[np.where(frac > 1.001)])
            )


#################################################################################################


def f2t_read_partposit(ifile, maxn=3e6, verbose=False):
    """
    @action: reads binary outputs from FLEXPART
    @input:  partposit_DATE.gz
    @output: returns a numpy array of dimension nparcels x 13
    @author: Jessica Keune 06/2020
    #modified: Dominik Schumacher, 06/2020 ---> do use pid!
    #modified: Jessica Keune, 10/2020 --> speeeeeedup!!! 
    # ATTN: hardcoded for 60 bytes / parcel from FP-ERA-INT
    """
    nbytes_per_parcel = 8 + 4 + 12 * 4
    with gzip.open(ifile, "rb") as strm:
        # skip header
        _ = strm.read(4)  # dummy
        _ = struct.unpack("i", strm.read(4))[0]  # time
        # grep full binary data set (ATTN: 60 bytes for FP-ERA-Int hardcoded)
        tdata = strm.read(int(maxn) * nbytes_per_parcel)
        # get number of parcels from length of tdata
        nparc = math.floor(len(tdata) / (nbytes_per_parcel))
        # decode binary data
        pdata = struct.unpack(
            (nparc) * "2fi3fi8f", tdata[0 : ((nparc) * nbytes_per_parcel)]
        )
        flist = list(pdata)
    strm.close()
    pdata = np.reshape(flist, newshape=(nparc, 15))[:, 2:]
    # remove last line if data is bs (pid = -99999)
    if np.any(pdata[:, 0] < 0):
        pdata = np.delete(pdata, np.where(pdata[:, 0] < 0), axis=0)
    return pdata


def maskgrabber(maskfile, maskvar="mask", latvar="lat", lonvar="lon"):
    # load
    with nc4.Dataset(maskfile, mode="r") as f:
        mask = np.asarray(f[maskvar][:])
        lat = np.asarray(f[latvar][:])
        lon = np.asarray(f[lonvar][:])
    # check if 2dimensional; necessary for ERA-I lsm mask
    if len(mask.shape) == 3:
        mask = mask[0, :, :]
    # lats check (order irrelevant, just must be within [-90,90])
    if not (lat.min() == -90 or lat.max() == 90):
        return None
    # lats check - now check order (currently required: -90 --> 90) and adjust if needed
    if not latsok(lat):
        mask, lat = ncdf_fliplats(mask, lat, lataxis=0)
    # lons check
    if lon.min() == -180 and lon.max() < 180:
        pass
    elif np.array_equal(lon, np.arange(0, 360)):
        mask, lon = ncdf_lon360to180(mask, lon, 1)
    else:
        # this case is not implemented
        return None
    return (mask, lat, lon)


def ncdf_lon360to180(ary, lons, lonaxis=1):
    # bring lon axis to front to handle any shape of ary
    ary = np.moveaxis(ary, lonaxis, 0)
    ary_bu = np.copy(ary)
    lons_bu = np.copy(lons)
    # copy lons & data
    lons[: int(lons.size / 2)] = lons_bu[int(lons.size / 2) :] - 360
    lons[int(lons.size / 2) :] = lons_bu[: int(lons.size / 2)]
    ary[: int(lons.size / 2)] = ary_bu[int(lons.size / 2) :]
    ary[int(lons.size / 2) :] = ary_bu[: int(lons.size / 2)]
    # move axis back to where it was
    ary = np.moveaxis(ary, 0, lonaxis)
    return (ary, lons)


def latsok(lats):
    # returns True if no flip needed; returns False if flip needed
    # standard for now: -90 --> 90 (contrary to ERA-Int)
    if lats[0] < lats[-1]:
        return True
    elif lats[0] > lats[-1]:
        return False


def ncdf_fliplats(ary, lats, lataxis=0):
    flip_ary = np.flip(ary, axis=lataxis)
    flip_lats = np.flipud(lats)
    return flip_ary, flip_lats


def writemasknc(mask, mlat, mlon, ofile="mask.nc", currentversion=get_currentversion()):
    # create netCDF4 instance
    nc_f = nc4.Dataset(ofile, "w", format="NETCDF4")
    ### create dimensions ###
    nc_f.createDimension("lat", mlat.size)
    nc_f.createDimension("lon", mlon.size)
    # create variables
    latitudes = nc_f.createVariable("lat", "f8", "lat")
    longitudes = nc_f.createVariable("lon", "f8", "lon")
    ncmask = nc_f.createVariable(
        "mask", "i4", ("lat", "lon"), fill_value=nc4.default_fillvals["i4"]
    )
    # set attributes
    nc_f.title = "HAMSTER: mask"
    today = datetime.datetime.now()
    nc_f.history = "Created " + today.strftime("%d/%m/%Y %H:%M:%S") + " using HAMSTER."
    nc_f.institution = (
        "Hydro-Climate Extremes Laboratory (H-CEL), Ghent University, Ghent, Belgium"
    )
    nc_f.source = (
        "HAMSTER " + str(currentversion) + " ((c) Dominik Schumacher and Jessica Keune)"
    )
    latitudes.units = "degrees_north"
    longitudes.units = "degrees_east"
    ncmask.units = "-"
    # write data
    longitudes[:] = mlon
    latitudes[:] = mlat
    ncmask[:] = mask
    # close file
    nc_f.close()
    print("\n * Created " + str(ofile) + " !")


def extendmask(mask, mlat, mlon, maskval, nx=5, ny=5, debug=False):
    imlat, imlon = np.where(mask == maskval)
    lat1 = mlat[imlat].min() - ny
    lat2 = mlat[imlat].max() + ny
    lon1 = mlon[imlon].min() - nx
    lon2 = mlon[imlon].max() + nx
    if lon1 < -180:
        lon1 += 360
    if lon2 > 180:
        lon2 -= 360
    # super ugly, but does what it should..
    ilats = np.where((mlat >= lat1) & (mlat <= lat2))
    if lon1 < lon2:
        ilons = np.where((mlon >= lon1) & (mlon <= lon2))
    else:
        ilons = np.where((mlon >= lon1) | (mlon <= lon2))
    extmask = np.zeros(shape=mask.shape)
    dmask = np.zeros(shape=mask.shape)
    dmask[ilats, :] += 1
    dmask[:, ilons] += 1
    extmask[np.where(dmask == 2)] = maskval
    if debug:
        writemasknc(extmask, mlat, mlon, "extmask.nc")
    return extmask


def nextmonth(ddate):
    nyyyy = (ddate + relativedelta(months=1)).strftime("%Y")
    nmm = (ddate + relativedelta(months=1)).strftime("%m")
    return (nyyyy, nmm)


def datetime2date(datetimeseq):
    dateseq = np.unique([idate.date() for idate in datetimeseq]).tolist()
    return dateseq


def timelord(startdate, enddate, timestep, ret="", fformat="%Y%m%d%H"):
    datetime_seq = []
    fdatetime_seq = []
    ffdatetime_seq = []
    idatetime = startdate
    # create datetime string & datetime object
    while idatetime <= enddate:
        datetime_seq.append(idatetime.strftime(fformat))
        fdatetime_seq.append(idatetime)
        ffdatetime_seq.append(idatetime.strftime("%Y%m%d%H") + "0000")
        idatetime += timestep
    if ret == "string":
        return datetime_seq
    elif ret == "datetime":
        return fdatetime_seq
    elif ret == "fileformat":
        return ffdatetime_seq
    else:
        return (datetime_seq, fdatetime_seq, ffdatetime_seq)


def f2t_timelord(ntraj_d, dt_h, tbgn, tend):
    fulltime = []
    fulltime.append(tbgn - datetime.timedelta(days=ntraj_d, hours=dt_h))
    while fulltime[-1] < tend:
        fulltime.append(fulltime[-1] + datetime.timedelta(hours=dt_h))
    # convert to strings in matching format for partposit files
    fulltime_str = [dft.strftime("%Y%m%d%H%M%S") for dft in fulltime]
    return fulltime_str


def f2t_loader(ifile, fixlons=True, fixids=True, verbose=True):
    dummy = f2t_read_partposit(ifile, verbose=False)
    ## fix parcel ID's (ATTN: specific to the global FP-ERA-Interim run!)
    if fixids:
        dummy[:, 0] = f2t_fixid(IDs=dummy[:, 0], verbose=verbose)  # fix IDs
    ## fix longitudes
    if fixlons:
        dummy[:, 1][dummy[:, 1] >= 179.5] -= 360
    ## sort array by parcel ID
    sorted_dummy = dummy[np.argsort(dummy[:, 0]), :]
    return sorted_dummy


def f2t_fixid(IDs, verbose, thresidx=1997000):
    # get duplicate IDs independent of threshidx...
    # u, c    = np.unique(IDs, return_counts=True)
    # dup     = u[c > 1]
    # print("\t --> Number of duplicates (unconditional): "+ str(len(dup)))
    # add 2e6 to the second duplicate...
    # for i in range(len(dup)):
    #    IDs[np.where(IDs==dup[i])[0][1]] += 2e6
    ## simply shift to indices > 2e6
    IDs[thresidx:][IDs[thresidx:] < (2e6 - thresidx)] += 2e6
    if verbose:
        ndupl = np.where(IDs > 2e6)[0].size
        if ndupl == 0:
            print("        --> NO duplicates present")
        else:
            print("        --> " + str(ndupl) + " duplicate IDs shifted")
    return IDs


def is_parcel_in_mask(plat, plon, mlat, mlon, mask, maskval):
    imlat = np.argmin(np.abs(mlat - plat))
    imlon = np.argmin(np.abs(mlon - plon))
    if mask[imlat, imlon] == maskval:
        return True
    else:
        return False


def find_potential_parcels(plon, plat, mlon, mlat, mask, maskval):
    ## search for potential candidates using rectangular box using min/max of mask extent
    imlat, imlon = np.where(mask == maskval)
    dlon = abs(mlon[0] - mlon[1]) / 2
    dlat = abs(mlat[0] - mlat[1]) / 2
    lat1 = mlat[imlat].min() - dlat
    lat2 = mlat[imlat].max() + dlat
    lon1 = mlon[imlon].min() - dlon
    lon2 = mlon[imlon].max() + dlon
    ## go for it (this is pretty fast)
    idx_inbox = np.where(
        (plon >= lon1) & (plon <= lon2) & (plat >= lat1) & (plat <= lat2)
    )[0]
    return idx_inbox  # ATTN: returns index only (not pid)


def f2t_seeker(array3D, mask, val, lat, lon):
    ## check if anything needs to be done at all
    if mask is None:
        return array3D[-1, :, 0][np.where(~np.isnan(array3D[-1, :, 0]))]
    ## use an extended mask to make sure we get everything (arriving + midpoint parcels!)
    extmask = extendmask(
        mask=mask, mlat=lat, mlon=lon, maskval=val, nx=5, ny=5, debug=False
    )
    ## first, we search potential candidates arriving in extended mask, using a rectangular box
    idx_inbox = find_potential_parcels(
        plon=array3D[-1, :, 1],
        plat=array3D[-1, :, 2],
        mlon=lon,
        mlat=lat,
        mask=extmask,
        maskval=val,
    )
    # work with subset of array3D for all potential candidates
    carray3D = f2t_constructor(
        array3D=array3D[-2:, :, :], pid=array3D[-1, idx_inbox, 0], time_str=["-2", "-1"]
    )
    ## now check if *really* in mask (slow!)
    pid = []
    for ii in range(idx_inbox.size):
        # check if arrival point in mask
        if is_parcel_in_mask(
            plon=carray3D[-1, ii, 1],
            plat=carray3D[-1, ii, 2],
            mlon=lon,
            mlat=lat,
            mask=mask,
            maskval=val,
        ):
            pid.append(carray3D[-1, ii, 0])
        else:
            if np.any(np.isnan(carray3D[:, ii, 0:3])):
                # parcel disappeared --> skip
                continue
            # check if midpoint is mask
            # (only calc. if arrival point is not in mask, to speed up process)
            midlat, midlon = midpoint_on_sphere2(
                carray3D[-1, ii, 2],
                carray3D[-1, ii, 1],
                carray3D[-2, ii, 2],
                carray3D[-2, ii, 1],
            )
            if is_parcel_in_mask(
                plon=midlon, plat=midlat, mlon=lon, mlat=lat, mask=mask, maskval=val
            ):
                pid.append(carray3D[-1, ii, 0])
    ## finally, return parcel IDs
    return np.asarray(pid)


def f2t_locator(array2D, pid, tstring):
    ## figure out where parcels are (lines may shift b/w files)
    pidx = np.where(np.isin(array2D[:, 0], pid, assume_unique=False))[
        0
    ]  # <----- set True ??
    chosen = np.NaN * np.ones(shape=(len(pid), array2D.shape[1]))
    if not pidx.size == len(pid):
        ## ATTN: this needs to be adjusted for other runs...
        print("---- INFO: not all parcels present in file --> partposit_" + tstring)
        idx_pidok = np.where(np.isin(pid, array2D[pidx, 0], assume_unique=False))[
            0
        ]  # <----- set True ??
        chosen[idx_pidok, :] = array2D[pidx, :]
    else:
        chosen[:, :] = array2D[pidx, :]
    return chosen


def f2t_constructor(array3D, pid, time_str):
    ## sloppy check
    if not array3D.shape[0] == len(time_str):
        raise IndexError("time_str must match time dimension of array3D!")
    ## prepare large array, loop thru
    trajs = np.empty(shape=(array3D.shape[0], pid.size, array3D.shape[2]))
    for ii in range(array3D.shape[0]):
        ## call locator
        trajs[ii, :, :] = f2t_locator(
            array2D=array3D[ii, :, :], pid=pid, tstring=time_str[ii]
        )
    return trajs


def f2t_saver(odata, outdir, fout, tstring):
    with h5py.File(outdir + "/" + fout + "_" + tstring + ".h5", "w") as f:
        f.create_dataset("trajdata", data=odata)


def f2t_establisher(
    partdir, selvars, time_str, ryyyy, mask, maskval, mlat, mlon, outdir, fout, verbose
):
    ##-- 1.) load em files
    data = np.empty(shape=(len(time_str), 2000001, selvars.size))
    for ii in range(len(time_str)):
        if verbose:
            print("       " + time_str[ii][:-4], end="")
        ifile = partdir + "/partposit_" + time_str[ii] + ".gz"
        dummy = f2t_loader(ifile, fixlons=True, fixids=True, verbose=verbose)[
            :, selvars
        ]  # load
        data[ii, : dummy.shape[0]] = dummy[:]  # fill only where data available
        data[ii, dummy.shape[0] :] = np.NaN

    ##-- 2.) find IDs within mask
    if verbose:
        print("       searching IDs", end="")
    pid_inmask = f2t_seeker(
        array3D=data[:, :, :], mask=mask, val=maskval, lat=mlat, lon=mlon
    )

    ##-- 3.) grab data for IDs
    if verbose:
        print(" | grabbing data for " + str(pid_inmask.size) + " IDs", end="")
    trajs = f2t_constructor(array3D=data, pid=pid_inmask, time_str=time_str)

    ##--4.) save
    if verbose:
        print(" | writing to file", end="")
    f2t_saver(odata=trajs, outdir=outdir, fout=fout, tstring=time_str[-1][:-4])

    ##--5.) return data & trajs arrays (needed for next files)
    return (data, trajs)


def f2t_ascender(
    data,
    partdir,
    selvars,
    ryyyy,
    time_str,
    mask,
    maskval,
    mlat,
    mlon,
    outdir,
    fout,
    verbose,
):
    ##--1.) move old data & fill current step with new data
    if verbose:
        print("\n      ", time_str[-1][:-4], end="")
    # use loop to avoid RAM spike here
    for ii in range(len(time_str) - 1):
        data[ii, :, :] = data[ii + 1, :, :]
    # load new data | rely on dummy variable
    ifile = partdir + "/partposit_" + time_str[-1] + ".gz"
    dummy = f2t_loader(ifile, fixlons=True, fixids=True, verbose=verbose)[:, selvars]
    # insert new data, use NaN for rest
    data[-1, : dummy.shape[0]] = dummy[:]
    data[-1, dummy.shape[0] :] = np.NaN

    ##--2.) find all IDs
    if verbose:
        print("       searching IDs", end="")
    pid_inmask = f2t_seeker(
        array3D=data[:, :, :], mask=mask, val=maskval, lat=mlat, lon=mlon
    )

    ##--3.) construct new trajectories (trajs recycling option has been removed)
    if verbose:
        print(
            " | constructing trajs for " + str(pid_inmask.size) + " IDs from scratch",
            end="",
        )
    trajs = f2t_constructor(array3D=data, pid=pid_inmask, time_str=time_str[:])

    ##--4.) save
    if verbose:
        print(" | writing to file", end="")
    f2t_saver(
        odata=trajs, outdir=outdir, fout=fout, tstring=time_str[-1][:-4]
    )  # omit mins & secs

    ##--5.) return updated data & trajs arrays
    return (data, trajs)


def checknan(x):
    x[x >= 9.9e36] = np.nan
    return x


def whereinmask(mask, maskval, masklat, masklon, trajlat, trajlon):
    ## first, we search potential candidates using rectangular box
    imlat, imlon = np.where(mask == maskval)
    lat1 = masklat[imlat].min() - 0.5
    lat2 = masklat[imlat].max() + 0.5
    lon1 = masklon[imlon].min() - 0.5
    lon2 = masklon[imlon].max() + 0.5
    ## go for it (this is pretty fast)
    idx_inbox = np.where(
        (trajlon >= lon1) & (trajlon <= lon2) & (trajlat >= lat1) & (trajlat <= lat2)
    )[0]
    ## now check if *really* in mask (slow!)
    idx = []
    for ii in range(idx_inbox.size):
        jdx = idx_inbox[ii]
        if (
            mask[
                np.argmin(np.abs(masklat - trajlat[jdx])),
                np.argmin(np.abs(masklon - trajlon[jdx])),
            ]
            == maskval
        ):
            idx.append(jdx)
    ## finally, return indices for which traj in mask
    return np.asarray(idx)


def maxlastn(series, n=4):
    maxy = np.zeros(shape=(n, series.size))
    maxy[0, :] = series[:]
    for ii in range(1, n):
        maxy[ii, :-ii] = series[ii:]
    return np.max(maxy, axis=0)


def grabhpbl_partposit(ifile):
    dummy = f2t_loader(ifile, fixlons=True, fixids=True, verbose=False)[
        :, [0, 9]
    ]  # 0: id, 9: hpbl
    return dummy


def grabmesomehpbl(filelist, verbose):
    extendarchive = []

    if verbose:
        print(
            "\n--------------------------------------------------------------------------------------"
        )
        print(
            "\n ! performing pre-loop to extend trajectory data & achieve advection month-2-month consistency"
        )
        print("\n ! estimating remaining time for pre-loop ...")
        pretic = timeit.default_timer()

    for it in range(len(filelist)):
        if verbose and it == 1:
            pretoc = timeit.default_timer()
            mins = round((len(filelist)) * (pretoc - pretic) / 60, 2)
            if mins < 1.0:
                print(
                    "  ---> "
                    + str(mins * 60)
                    + " seconds to go, how about some stretching in the meantime?"
                )
            else:
                print(
                    "  ---> "
                    + str(mins)
                    + " minutes to go, might want to grab a coffee..."
                )

        # append data
        extendarchive.append(grabhpbl_partposit(filelist[it]))

    return extendarchive
