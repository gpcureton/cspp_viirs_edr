#!/usr/bin/env python
# encoding: utf-8
"""
QuarterlySurfaceType.py

 * DESCRIPTION:  Class to granulate the IGBP Quarterly Surface Type data product 

Created by Geoff Cureton on 2013-03-05.
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
    LOG = logging.getLogger('QuarterlySurfaceType')

from Utils import getURID, getAscLine, getAscStructs, shipOutToFile
from Utils import index, find_lt, find_le, find_gt, find_ge

class QuarterlySurfaceType() :

    def __init__(self,inDir=None, sdrEndian=None, ancEndian=None):
        self.collectionShortName = 'VIIRS-GridIP-VIIRS-Qst-Mod-Gran'
        self.xmlName = 'VIIRS_GRIDIP_VIIRS_QST_MOD_GRAN.xml'
        self.blobDatasetName = 'igbp'
        self.dataType = 'uint8'
        self.sourceType = 'IGBP'
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

        # QST type value enumerations to be used with the igbp field
        self.IGBP_dict = {
            'IGBP_EVERNEEDLE'    : 1,
            'IGBP_EVERBROAD'     : 2,
            'IGBP_DECIDNEEDLE'   : 3,
            'IGBP_DECIDBROAD'    : 4,
            'IGBP_MIXEDFOREST'   : 5,
            'IGBP_CLOSEDSHRUBS'  : 6,
            'IGBP_OPENSHRUBS'    : 7,
            'IGBP_WOODYSAVANNA'  : 8,
            'IGBP_SAVANNA'       : 9,
            'IGBP_GRASSLAND'     : 10,
            'IGBP_WETLANDS'      : 11,
            'IGBP_CROPLANDS'     : 12,
            'IGBP_URBAN'         : 13,
            'IGBP_CROPNATMOSAIC' : 14,
            'IGBP_SNOWICE'       : 15,
            'IGBP_BARREN'        : 16,
            'IGBP_WATERBODIES'   : 17,
            'IGBP_UNCLASSIFIED'  : 30,
            'IGBP_FILL'          : 31,
            'IGBP_MIN'           : 1,
            'IGBP_MAX'           : 17
        }


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
