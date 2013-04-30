#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
awips2_gtm_edr
~~~~~~~~~~~~~~

Convert HDF5 format Ground-Track Mercator (GTM) Imagery EDR files to NetCDF4,
compatible with AWIPS2.

From code written by DJHoese, Apr2013

:copyright: 2013 by University of Wisconsin Regents, see AUTHORS for more details
:license: GPLv3, see LICENSE for more details
"""
__author__ = 'rayg'
__docformat__ = 'reStructuredText'

import os
import sys
import re
import logging
import unittest
import datetime

import h5py
from netCDF4 import Dataset
import numpy as np

LOG = logging.getLogger(__name__)


def h5path(elf, path):
    if not path:
        return elf
    if isinstance(path, str):
        path = h5path(elf, path.split('/'))
    rx = re.compile(path[0])
    for k,v in elf.iteritems():
        if rx.match(k):
            return h5path(v, path[1:])
    return None

EDR_PATH = r'All_Data/VIIRS-.*-EDR_All/BrightnessTemperature'
# where to find Latitude and Longitude
GEO_PATH = r'All_Data/VIIRS-.*GTM-EDR_GEO_All'
TIME_PATH = r'Data_Products/VIIRS-.*GTM-EDR-GEO/VIIRS-.*-EDR-GEO_Aggr'
GRING_PATH = r'Data_Products/VIIRS-.*GTM-EDR-GEO/VIIRS-.*-EDR-GEO_Gran_0'
g_ring_lat = geo_file["Data_Products"]["VIIRS-IMG-GTM-EDR-GEO"]["VIIRS-IMG-GTM-EDR-GEO_Gran_0"].attrs["G-Ring_Latitude"][:].astype(numpy.float64)

# Time helpers
class UTC(datetime.tzinfo):
    """Time zone class for UTC
    """
    ZERO = datetime.timedelta(0)
    def utcoffset(self, dt):
        return self.ZERO
    def tzname(self, dt):
        return "UTC"
    def dst(self, dt):
        return self.ZERO

def utc_now():
    return datetime.datetime.utcnow().replace(tzinfo=UTC())


class Granule(object):
    """

    """
    def __init__(self, edr_path, geo_path):
        self.edr = h5py.File(edr_path, 'r')
        self.geo = h5py.File(geo_path, 'r')

    def data(self):
        return h5path(self.edr, EDR_PATH)[:]

    def lat_lon(self):
        base = h5path(self.geo, GEO_PATH)
        return base['Latitude'][:], base['Longitude'][:]

    def start_end(self):
        base = h5path(self.geo, TIME_PATH)
        start_date = base.attrs["AggregateBeginningDate"][0,0]
        end_date = base.attrs["AggregateEndingDate"][0,0]
        start_time = base.attrs["AggregateBeginningTime"][0,0]
        end_time = base.attrs["AggregateEndingTime"][0,0]
        return as_datetime(start_date, start_time), as_datetime(end_date, end_time)

    def gring_lat_lon(self):
        # FIXME this doesn't deal with aggregates having more than one granule
        base = h5path(self.geo, GRING_PATH)
        return base.attrs["G-Ring_Latitude"][:].astype(numpy.float64), base.attrs["G-Ring_Longitude"][:].astype(numpy.float64)

    # FIXME deal with scale factors

    def geo_group(self):
        for k,v in self.geo['All_Data'].iteritems():
            if RE_GEO.match(k):
                return v
        return None

    def time_range(self):
        for k,v in self.geo['Data_Products'].iteritems():
            if RE_GEO.match(k):
                return v






# Read the input files
img_file = h5py.File("./VI4BO_npp_d20130116_t0944041_e0945402_b06328_c20130305060723097136_noaa_ops.h5")
geo_file = h5py.File("./GIGTO_npp_d20130116_t0944041_e0945402_b06328_c20130305060441825409_noaa_ops.h5")

img_data = img_file["All_Data"]["VIIRS-I4-IMG-EDR_All"]["BrightnessTemperature"][:,:]
img_factors = img_file["All_Data"]["VIIRS-I4-IMG-EDR_All"]["BrightnessFactors"][:]
lat_data = geo_file["All_Data"]["VIIRS-IMG-GTM-EDR-GEO_All"]["Latitude"][:,:]
lon_data = geo_file["All_Data"]["VIIRS-IMG-GTM-EDR-GEO_All"]["Longitude"][:,:]
start_date = geo_file["Data_Products"]["VIIRS-IMG-GTM-EDR-GEO"]["VIIRS-IMG-GTM-EDR-GEO_Aggr"].attrs["AggregateBeginningDate"][0,0]
start_time = geo_file["Data_Products"]["VIIRS-IMG-GTM-EDR-GEO"]["VIIRS-IMG-GTM-EDR-GEO_Aggr"].attrs["AggregateBeginningTime"][0,0]
end_date = geo_file["Data_Products"]["VIIRS-IMG-GTM-EDR-GEO"]["VIIRS-IMG-GTM-EDR-GEO_Aggr"].attrs["AggregateEndingDate"][0,0]
end_time = geo_file["Data_Products"]["VIIRS-IMG-GTM-EDR-GEO"]["VIIRS-IMG-GTM-EDR-GEO_Aggr"].attrs["AggregateEndingTime"][0,0]
g_ring_lat = geo_file["Data_Products"]["VIIRS-IMG-GTM-EDR-GEO"]["VIIRS-IMG-GTM-EDR-GEO_Gran_0"].attrs["G-Ring_Latitude"][:].astype(numpy.float64)
g_ring_lon = geo_file["Data_Products"]["VIIRS-IMG-GTM-EDR-GEO"]["VIIRS-IMG-GTM-EDR-GEO_Gran_0"].attrs["G-Ring_Longitude"][:].astype(numpy.float64)

# Make data manipulations
mid_idx = lat_data.shape[-1]/2
lat_envelope = lat_data[:,(0,mid_idx,-1)]
lon_envelope = lon_data[:,(0,mid_idx,-1)]

# Write the output file
nc_file = Dataset("VIIRS_I4_IMG_EDR_TIPB99_KNES_npp_s201301160944041_e201301160945402_c20130305060441825409.nc", mode="w")
#nc_file = Dataset("viirs_img_edr_20121001.nc", mode="w")

# Create Dimensions
row_dim_name = "AlongTrack-%d" % img_data.shape[0]
fac_dim_name = "Granule-%d" % img_factors.shape[0]
col_dim_name = "CrossTrack-%d" % img_data.shape[1]
env_dim_name = "CrossTrack-%d" % 3
nc_file.createDimension(row_dim_name, img_data.shape[0])
nc_file.createDimension(fac_dim_name, img_factors.shape[0])
nc_file.createDimension(col_dim_name, img_data.shape[1])
nc_file.createDimension(env_dim_name, 3)

# Create Global Attributes
nc_file.time_coverage_start = datetime.datetime.strptime(start_date + start_time.split(".")[0], "%Y%m%d%H%M%S").strftime("%Y-%m-%dT%H:%M:%SZ")
nc_file.time_coverage_end   = datetime.datetime.strptime(end_date + end_time.split(".")[0], "%Y%m%d%H%M%S").strftime("%Y-%m-%dT%H:%M:%SZ")
nc_file.date_created = utc_now().strftime("%Y-%m-%dT%H:%M:%SZ")
print g_ring_lat.dtype
print g_ring_lat
print g_ring_lon
g_ring_lat = numpy.append(g_ring_lat, g_ring_lat[0])
g_ring_lon = numpy.append(g_ring_lon, g_ring_lon[0])
nc_file.setncattr("g_ring_latitude", g_ring_lat)
nc_file.setncattr("g_ring_longitude", g_ring_lon)

# Create and write Variables
# Image data
bt_var = nc_file.createVariable("BrightnessTemperature@VIIRS-I4-IMG-EDR", 'u2', dimensions=(row_dim_name,col_dim_name))
bt_var[:,:] = img_data
bt_var.setncattr("missing_value", "65535 65534 65533 65532 65531 65530 65529 65528")
# XXX: Need missing_valuename?

# Scaling Factors
bt_factors_var = nc_file.createVariable("BrightnessFactors@VIIRS-I4-IMG-EDR", 'f4', dimensions=(fac_dim_name,))
# Factors should be 2 but this aggregate has more
bt_factors_var[:] = img_factors

# Navigation
lat_var = nc_file.createVariable("Latitude@VIIRS-IMG-GTM-EDR-GEO", 'f4', dimensions=(row_dim_name,env_dim_name))
lat_var[:,:] = lat_envelope
lat_var.setncattr("missing_value", "-999.9 -999.8 -999.5 -999.4 -999.3")
lon_var = nc_file.createVariable("Longitude@VIIRS-IMG-GTM-EDR-GEO", 'f4', dimensions=(row_dim_name,env_dim_name))
lon_var[:,:] = lon_envelope
lon_var.setncattr("missing_value", "-999.9 -999.8 -999.5 -999.4 -999.3")

nc_file.sync()
nc_file.close()

