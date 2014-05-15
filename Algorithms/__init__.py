#!/usr/bin/env python
# encoding: utf-8
"""
This module contains data required for any VIIRS EDR algorithms.

Created by Geoff Cureton on 2013-02-27.
Copyright (c) 2013 University of Wisconsin SSEC. All rights reserved.
"""

file_Date = '$Date$'
file_Revision = '$Revision$'
file_Author = '$Author$'
file_HeadURL = '$HeadURL$'
file_Id = '$Id$'

__author__ = 'G.P. Cureton <geoff.cureton@ssec.wisc.edu>'
__version__ = '$Id$'
__docformat__ = 'Epytext'

import CloudMaskIP
import AerosolOpticalThicknessIP
import SeaSurfaceTemperatureEDR
import SurfaceReflectanceIP
import VegetationIndexEDR
import SurfaceTypeEDR
#import LandSurfaceTemperatureEDR
import dummy_multiprocessing

modules = {}
modules['VCM'] = 'CloudMaskIP'
modules['AOT'] = 'AerosolOpticalThicknessIP'
modules['SST'] = 'SeaSurfaceTemperatureEDR'
modules['SRFREF'] = 'SurfaceReflectanceIP'
modules['VI'] = 'VegetationIndexEDR'
modules['ST'] = 'SurfaceTypeEDR'
#modules['LST'] = 'LandSurfaceTemperatureEDR'
modules['ALL'] = 'VegetationIndexEDR'
modules['MPC'] = 'dummy_multiprocessing'

meta_algorithms = ['ATMOS','LAND','OCEAN','ALL']

# What are the cross granule requirements for each algorithm?
crossGranules = {}
crossGranules['VCM'] = 1
crossGranules['AOT'] = 1
crossGranules['SST'] = 0
crossGranules['SRFREF'] = 0
crossGranules['VI'] = 0
crossGranules['ST'] = 0
#crossGranules['LST'] = 1
crossGranules['MPC'] = 1
crossGranules['ALL'] = 0
crossGranules['ATMOS'] = 0
crossGranules['LAND'] = 0
crossGranules['OCEAN'] = 0

# What previous alg products are *explicitly* required for each algorithm.
prerequisites = {}
prerequisites['VCM'] = []
prerequisites['AOT'] = ['VCM']
prerequisites['SST'] = ['VCM','AOT']
prerequisites['SRFREF'] = ['VCM','AOT']
prerequisites['VI'] = ['VCM','AOT','SRFREF']
prerequisites['ST'] = []
#prerequisites['LST'] = []
prerequisites['ATMOS'] = ['VCM','AOT']
prerequisites['LAND'] = ['VCM','AOT','SRFREF','VI']
prerequisites['OCEAN'] = ['VCM','AOT','SST']
prerequisites['ALL'] = ['VCM','AOT','SST','SRFREF','VI']
prerequisites['MPC'] = []

geo_hdf5_prefix={}
geo_hdf5_prefix['VIIRS-MOD-GEO'] = 'GMODO'
geo_hdf5_prefix['VIIRS-MOD-GEO-TC'] = 'GMTCO'
geo_hdf5_prefix['VIIRS-IMG-GEO'] = 'GIMGO'
geo_hdf5_prefix['VIIRS-IMG-GEO-TC'] = 'GITCO'

sdr_hdf5_prefix={}
sdr_hdf5_prefix['VIIRS-I1-SDR'] = 'SVI01'
sdr_hdf5_prefix['VIIRS-I2-SDR'] = 'SVI02'
sdr_hdf5_prefix['VIIRS-I3-SDR'] = 'SVI03'
sdr_hdf5_prefix['VIIRS-I4-SDR'] = 'SVI04'
sdr_hdf5_prefix['VIIRS-I5-SDR'] = 'SVI05'
sdr_hdf5_prefix['VIIRS-M1-SDR'] = 'SVM01'
sdr_hdf5_prefix['VIIRS-M2-SDR'] = 'SVM02'
sdr_hdf5_prefix['VIIRS-M3-SDR'] = 'SVM03'
sdr_hdf5_prefix['VIIRS-M4-SDR'] = 'SVM04'
sdr_hdf5_prefix['VIIRS-M5-SDR'] = 'SVM05'
sdr_hdf5_prefix['VIIRS-M6-SDR'] = 'SVM06'
sdr_hdf5_prefix['VIIRS-M7-SDR'] = 'SVM07'
sdr_hdf5_prefix['VIIRS-M8-SDR'] = 'SVM08'
sdr_hdf5_prefix['VIIRS-M9-SDR'] = 'SVM09'
sdr_hdf5_prefix['VIIRS-M10-SDR'] = 'SVM10'
sdr_hdf5_prefix['VIIRS-M11-SDR'] = 'SVM11'
sdr_hdf5_prefix['VIIRS-M12-SDR'] = 'SVM12'
sdr_hdf5_prefix['VIIRS-M13-SDR'] = 'SVM13'
sdr_hdf5_prefix['VIIRS-M14-SDR'] = 'SVM14'
sdr_hdf5_prefix['VIIRS-M15-SDR'] = 'SVM15'
sdr_hdf5_prefix['VIIRS-M16-SDR'] = 'SVM16'

edr_hdf5_prefix={}
edr_hdf5_prefix['VCM'] = ['IICMO','AVAFO']
edr_hdf5_prefix['AOT'] = ['IVAOT','GAERO','VAOOO','IVAMI','VSUMO']
edr_hdf5_prefix['SST'] = ['VSSTO']
edr_hdf5_prefix['SRFREF'] = ['IVISR']
edr_hdf5_prefix['VI'] = ['VIVIO']
edr_hdf5_prefix['ST'] = ['VSTYO']
#edr_hdf5_prefix['LST'] = ['VLSTO']
edr_hdf5_prefix['MPC'] = ['MPC']

edr_hdf5_Gran_0={}
edr_hdf5_Gran_0['VCM']    = ['/Data_Products/VIIRS-CM-IP/VIIRS-CM-IP_Gran_0',\
                             '/Data_Products/VIIRS-AF-EDR/VIIRS-AF-EDR_Gran_0']
edr_hdf5_Gran_0['AOT']    = ['/Data_Products/VIIRS-Aeros-Opt-Thick-IP/VIIRS-Aeros-Opt-Thick-IP_Gran_0',\
                             '/Data_Products/VIIRS-Aeros-EDR-GEO/VIIRS-Aeros-EDR-GEO_Gran_0',\
                             '/Data_Products/VIIRS-Aeros-EDR/VIIRS-Aeros-EDR_Gran_0',\
                             '/Data_Products/VIIRS-Aeros-Modl-Info-IP/VIIRS-Aeros-Modl-Info-IP_Gran_0',\
                             '/Data_Products/VIIRS-SusMat-EDR/VIIRS-SusMat-EDR_Gran_0']
edr_hdf5_Gran_0['SST']    = ['/Data_Products/VIIRS-SST-EDR/VIIRS-SST-EDR_Gran_0']
edr_hdf5_Gran_0['SRFREF'] = ['/Data_Products/VIIRS-Surf-Refl-IP/VIIRS-Surf-Refl-IP_Gran_0']
edr_hdf5_Gran_0['VI']     = ['/Data_Products/VIIRS-VI-EDR/VIIRS-VI-EDR_Gran_0']
edr_hdf5_Gran_0['ST']     = ['/Data_Products/VIIRS-ST-EDR/VIIRS-ST-EDR_Gran_0']
#edr_hdf5_Gran_0['LST']    = ['/Data_Products/VIIRS-LST-EDR/VIIRS-LST-EDR_Gran_0']
edr_hdf5_Gran_0['MPC']    = ['/Data_Products/VIIRS-DUMMY/VIIRS-DUMMY_Gran_0']
