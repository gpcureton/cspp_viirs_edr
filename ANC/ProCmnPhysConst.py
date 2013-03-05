#!/usr/bin/env python
# encoding: utf-8
"""
ProCmnPhysConst.py

DESCRIPTION:  Fundamental Physical Constants in SI units.
              National Institute of Standards and Technology (NIST)
              CODATA 2002 constants:
              http://physics.nist.gov/cuu/Constants/index.html
              Plus SI Prefixes
              http://physics.nist.gov/cuu/Units/prefixes.html

This file is based on ProCmnPhysConst.h from ADL 4.1.

Created by Geoff Cureton on 2009-04-04.
Copyright (c) 2011 University of Wisconsin SSEC. All rights reserved.
"""

file_Date = '$Date$'
file_Revision = '$Revision$'
file_Author = '$Author$'
file_HeadURL = '$HeadURL$'
file_Id = '$Id$'

__author__ = 'G.P. Cureton <geoff.cureton@ssec.wisc.edu>'
__version__ = '$Id$'
__docformat__ = 'Epytext'


# Value of PI/2 in Float32 format
#FLOAT32_PIO2 = 1.570795;

# Value of PI in Float32 format
#FLOAT32_PI = 3.141591;


#PI =       M_PI;                # math constant pi
#PIO2 =     M_PI_2;              # pi/2
#PIO4 =     M_PI_4;              # pi/4
#TREPIO2 =  3.0L*M_PI/2.0L;      # 3*pi/2
#TWOPI =    2.0L*M_PI;           # 2*pi

# used to convert degrees to radians
#DEG2RAD =  M_PI/180.0L;         # (pi/180)

# used to convert radians to degrees
#RAD2DEG =  180.0L/M_PI;         # (180/pi)

# used to convert degrees to arcsec
DEG2ARCSEC =  3600.0;         # seconds in a degree

# Universal Gas Constant - used to compute virtual temp.  J/(kg*K)
DRYGAS = 287.05;

# Average Atmospheric Pressure at Sea Level
AVG_PRESS_SEALVL = 1013.25e0;

# Gravity (m/s^2)
GRAVITY = 9.80665;

# Earth radius (m)
EARTH_RADIUS_METERS = 6371007.181;
EARTH_RADIUS_METERS_INVERSE = 1 / EARTH_RADIUS_METERS;

# The following constant definitions are the WGS84 Earth Ellipsoid
# constants.  The main reference for the WGS84 ellipsoid is
# NIMA Physical Geodesy web page: 164.214.2.59/GandG/wgs-84/egm96.htm

EQUAT_RAD = 6.37813700000000e+6; # Equatoral rad., meters WGS84
POLAR_RAD = 6.35675231424518e+6; # Polar radius, meters WGS84
ECCEN_SQ = 6.69437999014132e-3;  # Eccentricity Squared WGS84
FLATFAC = 3.35281066474748e-3;   # Flattening Factor WGS84

# USAF Orbit Analyst Manuals, circa 1978.  1 - eccen_sqr
DETIC2CENTRIC = 9.93305620009859e-1;
CENTRIC2DETIC = 1.00673949674228e+0;

# The following constant definitions are for time conversions

TAI2IET      = 1.0e+06;         # Conversion factor
MIN_IN_HOUR  = 60.0e+0;         # Number of minutes in an hour
SEC_IN_HOUR  = 3600.0e+0;       # Number of seconds in an hour
MJD_CONV_FAC = 2.4000005e+6;    # Factor to convert AJD to MJD
SEC_IN_DAY   = 8.64e+04;        # Number of seconds in a day
UJD58        = 2.43620450e+06;  # Jan 1 1958  UJD format
JAN012030    = 2.272147232e+09; # Jan 1 2030  TAI format
TJD_CONV_FAC = 32.184e+0;       # Factor to convert TAI to TJD
DEG_IN_HOUR  = 15.0e+0;         # Number of degrees in an hour

# The following constant definitions are for polarstereographic dataset
MINUS30      = -0.523598775598299e0;  # -30 degrees in radians
PLUS30       =  0.523598775598299e0;  # 30 degrees in radians

# The following constant definitions are for nwp ancillary granulation 
# declare constant for calculation of water vapor mixing ratio (r)
GAS   = 621.97;                #-- ratio of the molecular weight
                               #-- of water vapor to dry air
                               #-- in units of grams/kilogram
    
# declare just a few of the more popular of the twenty SI prefixes
MICRO =   0.000001;  #-- scale by 1/1000000th
MILLI =   0.001;     #-- scale by 1/1000th
CENTI =   0.01;      #-- scale by 1/100th
DECI  =   0.1;       #-- scale by 1/10th
DEKA  =   10.0;      #-- scale by 10x
HECTO =   100.0;     #-- scale by 100x 

# Kelvin/Celsius conversion factor
TCOEFF=273.15;

# Constant used to generate surface reflectance
# multiplier to convert pascal to atmospheres (1 atm = 101325 pascal)
PRESS_CONV = 1.0 / 101325.0; 

# Standard Atmosphere Surface Pressure
STDPSL=1013.0;

# Moist air adiabatic lapse rate is 6.5 K/Km (equivalent to 6.5 C/Km)
# Converted value would be .0065 C/m
MOIST_AIR_LAPSE_RATE = 6.5/1000;

# Constant used to convert atm-cm to Dobson units
ATM_CM2DOBSON = 1000.0;

