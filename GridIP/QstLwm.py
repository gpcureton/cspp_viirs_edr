#!/usr/bin/env python
# encoding: utf-8
"""
QstLwm.py

 * DESCRIPTION:  Class to granulate the QstLwm data product 

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



import os, sys, logging, traceback
from os import path,uname,environ
import string
import re
import uuid
from glob import glob
from time import time
from datetime import datetime,timedelta

from scipy import round_

import numpy as np
from numpy import ma
import copy
from bisect import bisect_left,bisect_right

import ctypes
from numpy.ctypeslib import ndpointer

import ViirsData

from NCEPtoBlob import NCEPclass

# skim and convert routines for reading .asc metadata fields of interest
import adl_blob
import adl_asc
from adl_asc import skim_dir, contiguous_granule_groups, granule_groups_contain, effective_anc_contains,_eliminate_duplicates,_is_contiguous, RDR_REQUIRED_KEYS, POLARWANDER_REQUIRED_KEYS
from adl_common import ADL_HOME, CSPP_RT_ANC_PATH, CSPP_RT_ANC_CACHE_DIR, COMMON_LOG_CHECK_TABLE

# every module should have a LOG object
try :
    sourcename= file_Id.split(" ")
    LOG = logging.getLogger(sourcename[1])
except :
    LOG = logging.getLogger('QstLwm')

from Utils import getURID, getAscLine, getAscStructs, shipOutToFile
from Utils import index, find_lt, find_le, find_gt, find_ge

class QstLwm() :

    def __init__(self,inDir=None, sdrEndian=None, ancEndian=None):
        self.collectionShortName = 'VIIRS-GridIP-VIIRS-Qst-Lwm-Mod-Gran'
        self.xmlName = 'VIIRS_GRIDIP_VIIRS_QST_LWM_MOD_GRAN.xml'
        self.blobDatasetName = 'qstlwm'
        self.dataType = 'uint8'
        self.sourceType = 'DEM-IGBP'
        self.sourceList = ['']
        self.trimObj = ViirsData.ViirsTrimTable()

        if inDir is None :
            self.inDir = path.abspath(path.curdir)
        else :
            self.inDir = inDir

        if sdrEndian is None :
            self.sdrEndian = adl_blob.LITTLE_ENDIAN
        else :
            self.sdrEndian = sdrEndian

        if ancEndian is None :
            self.ancEndian = adl_blob.LITTLE_ENDIAN
        else :
            self.ancEndian = ancEndian

        # QST Land Water Mask type value enumerations to be used with the qstlwm field
        self.QSTLWM_list = [
            'EVERGREEN_NEEDLELEAF_FOREST', 'EVERGREEN_BROADLEAF_FORESTS', 'DECIDUOUS_NEEDLELEAF_FORESTS', \
            'DECIDUOUS_BROADLEAF_FORESTS', 'MIXED_FORESTS', 'CLOSED_SHRUBLANDS', 'OPEN_SHRUBLANDS', \
            'WOODY_SAVANNAS', 'SAVANNAS', 'GRASSLANDS', 'PERMANENT_WETLANDS', 'CROPLANDS', \
            'URBAN_AND_BUILDUP_LANDS', 'CROPLAND_NATURAL_VEGETATION', 'SNOW_AND_ICE', \
            'BARREN_OR_SPARSELY_VEGETATED', \
            'OCEAN_SEA', 'INLAND_WATER', 'COASTAL_WATER', \
            'UNCLASSIFIED_LAND', 'FILL_VALUE'
        ]

        self.QSTLWM_dict = {
            'EVERGREEN_NEEDLELEAF_FOREST'  : 1,
            'EVERGREEN_BROADLEAF_FORESTS'  : 2,
            'DECIDUOUS_NEEDLELEAF_FORESTS' : 3,
            'DECIDUOUS_BROADLEAF_FORESTS'  : 4,
            'MIXED_FORESTS'                : 5,
            'CLOSED_SHRUBLANDS'            : 6,
            'OPEN_SHRUBLANDS'              : 7,
            'WOODY_SAVANNAS'               : 8,
            'SAVANNAS'                     : 9,
            'GRASSLANDS'                   : 10,
            'PERMANENT_WETLANDS'           : 11,
            'CROPLANDS'                    : 12,
            'URBAN_AND_BUILDUP_LANDS'      : 13, 
            'CROPLAND_NATURAL_VEGETATION'  : 14,
            'SNOW_AND_ICE'                 : 15,
            'BARREN_OR_SPARSELY_VEGETATED' : 16,
            'OCEAN_SEA'                    : 17,
            'INLAND_WATER'                 : 18,
            'COASTAL_WATER'                : 19,
            'UNCLASSIFIED_LAND'            : 20,
            'FILL_VALUE'                   : 255
        }

        self.DEM_to_QSTLWM = np.array([
                self.QSTLWM_dict['OCEAN_SEA'],                    # DEM 0 : QSTLWM 17
                self.QSTLWM_dict['EVERGREEN_NEEDLELEAF_FOREST'],  # DEM 1 : QSTLWM  1
                self.QSTLWM_dict['COASTAL_WATER'],                # DEM 2 : QSTLWM 19
                self.QSTLWM_dict['INLAND_WATER'],                 # DEM 3 : QSTLWM 18
                self.QSTLWM_dict['SAVANNAS'],                     # DEM 4 : QSTLWM  9
                self.QSTLWM_dict['OCEAN_SEA'],                    # DEM 5 : QSTLWM 17
                self.QSTLWM_dict['OCEAN_SEA'],                    # DEM 6 : QSTLWM 17
                self.QSTLWM_dict['OCEAN_SEA']                     # DEM 7 : QSTLWM 17
            ],dtype=('uint8'))


    def subset(self,latMinList,latMaxList,lonMinList,lonMaxList,latCrnList,lonCrnList):
        '''Subsets the global elevation dataset to cover the required geolocation range.'''
        pass


    def granulate(self,geoDicts):
        '''Granulates the input DEM files.'''
        pass


    def shipOutToFile(self):
        ''' Pass the current class instance to this Utils method to generate 
            a blob/asc file pair from the input ancillary data object.'''

        shipOutToFile(self)
