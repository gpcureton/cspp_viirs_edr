#!/usr/bin/env python
# encoding: utf-8
"""
TerrainGeopotentialHeight.py

 * DESCRIPTION:  Class to granulate the ViirsAncTerrainHeightType data product
                 from a Digital Elevation Model (DEM).

Created by Geoff Cureton on 2013-02-25.
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

# every module should have a LOG object
try :
    sourcename= file_Id.split(" ")
    LOG = logging.getLogger(sourcename[1])
except :
    LOG = logging.getLogger('TerrainGeopotentialHeight')

# import the superclass
from ANC import ANC

class TerrainGeopotentialHeight() :

    # Terrain height
    # FIXME : VIIRS-ANC-Surf-Ht-Mod-Gran is supposed to come from the Terrain-Eco-ANC-Tile 
    # FIXME :   gridded data. The granulated terrain height from this source is provided in
    # FIXME :   the VIIRS-MOD-RGEO-TC from the VIIRS SDR controller. If we get the non 
    # FIXME :   terrain correction geolocation (VIIRS-MOD-RGEO) we can use the surface
    # FIXME :   geopotential height, which is approximately the same.

    def __init__(self):
        self.collectionShortName = 'VIIRS-ANC-Surf-Ht-Mod-Gran'
        self.xmlName = 'VIIRS_ANC_SURF_HT_MOD_GRAN.xml'
        self.blobDatasetName = 'terrainGeopotentialHeight'
        self.dataType = 'int16'

    def ingest(self,ancBlob):
        '''
        Ingest the ancillary dataset.
        '''
        LOG.info("Ingesting of %s is to be handled differently." % (self.collectionShortName))        
        self.gridData = None
        #return getattr(ancBlob,self.blobDatasetName).astype(self.dataType)


    def granulate():
        '''
        Granulate the ancillary dataset.
        '''
        pass
