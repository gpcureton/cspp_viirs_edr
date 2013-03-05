#!/usr/bin/env python
# encoding: utf-8
"""
AerosolOpticalThicknessIP.py

 * DESCRIPTION:  Class containing data relevent to the VIIRS Aerosol Optical Thickness IP. 

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

# every module should have a LOG object
LOG = logging.getLogger('AerosolOpticalThicknessIP')

ANC_collectionShortNames = [
                           'VIIRS-ANC-Preci-Wtr-Mod-Gran',
                           'VIIRS-ANC-Temp-Surf2M-Mod-Gran',
                           'VIIRS-ANC-Wind-Speed-Mod-Gran',
                           'VIIRS-ANC-Wind-Direction-Mod-Gran',
                           'VIIRS-ANC-Press-Surf-Mod-Gran',
                           'VIIRS-ANC-Tot-Col-Mod-Gran'
                          ]

GridIP_collectionShortNames = [
                            'VIIRS-GridIP-VIIRS-Qst-Lwm-Mod-Gran',
                            'VIIRS-GridIP-VIIRS-Snow-Ice-Cover-Mod-Gran',
                            'VIIRS-GridIP-VIIRS-Nbar-Ndvi-Mod-Gran'
                          ]

controllerBinary = 'ProEdrViirsAerosolController.exe'

# Attribute paths for Aerosol EDR and IP
attributePaths = {}
attributePaths['MOD_GEO_TC'] = 'Data_Products/VIIRS-MOD-GEO-TC/VIIRS-MOD-GEO-TC_Gran_0/N_Granule_ID'
attributePaths['AOT_IP'] = 'Data_Products/VIIRS-Aeros-Opt-Thick-IP/VIIRS-Aeros-Opt-Thick-IP_Gran_0/N_Granule_ID'
attributePaths['AOT_EDR'] = 'Data_Products/VIIRS-Aeros-EDR/VIIRS-Aeros-EDR_Gran_0/N_Granule_ID'
attributePaths['SUSMAT_EDR'] = 'Data_Products/VIIRS-SusMat-EDR/VIIRS-SusMat-EDR_Gran_0/N_Granule_ID'

MOD_GEO_TC_GRANULE_ID_ATTR_PATH = 'Data_Products/VIIRS-MOD-GEO-TC/VIIRS-MOD-GEO-TC_Gran_0/N_Granule_ID'
AOT_IP_GRANULE_ID_ATTR_PATH = 'Data_Products/VIIRS-Aeros-Opt-Thick-IP/VIIRS-Aeros-Opt-Thick-IP_Gran_0/N_Granule_ID'
AOT_EDR_GRANULE_ID_ATTR_PATH = 'Data_Products/VIIRS-Aeros-EDR/VIIRS-Aeros-EDR_Gran_0/N_Granule_ID'
SUSMAT_EDR_GRANULE_ID_ATTR_PATH = 'Data_Products/VIIRS-SusMat-EDR/VIIRS-SusMat-EDR_Gran_0/N_Granule_ID'

# XML template for ProEdrViirsAerosolController.exe
xmlTemplate = """<InfTkConfig>
  <idpProcessName>ProEdrViirsAerosolController.exe</idpProcessName>
  <siSoftwareId></siSoftwareId>
  <isRestart>FALSE</isRestart>
  <useExtSiWriter>FALSE</useExtSiWriter>
  <debugLogLevel>LOW</debugLogLevel>
  <debugLevel>DBG_LOW</debugLevel>
  <dbgDest>D_FILE</dbgDest>
  <enablePerf>FALSE</enablePerf>
  <perfPath>${WORK_DIR}/perf</perfPath>
  <dbgPath>${WORK_DIR}/log</dbgPath>
  <initData>
     <domain>OPS</domain>
     <subDomain>SUBDOMAIN</subDomain>
     <startMode>INF_STARTMODE_COLD</startMode>
     <executionMode>INF_EXEMODE_PRIMARY</executionMode>
     <healthTimeoutPeriod>30</healthTimeoutPeriod>
  </initData>
  <lockinMem>FALSE</lockinMem>
  <rootDir>${WORK_DIR}</rootDir>
  <inputPath>${WORK_DIR}</inputPath>
  <outputPath>${WORK_DIR}</outputPath>
  <dataStartIET>0000000000000000</dataStartIET>
  <dataEndIET>1111111111111111</dataEndIET>
  <actualScans>0</actualScans>
  <previousActualScans>0</previousActualScans>
  <nextActualScans>0</nextActualScans> 
  <usingMetadata>TRUE</usingMetadata>
  <configGuideName>ProEdrViirsAerosolController_GuideList.cfg</configGuideName>

  <task>
    <taskType>EDR</taskType>
    <taskDetails1>%(N_Granule_ID)s</taskDetails1>
    <taskDetails2>%(N_Granule_Version)s</taskDetails2>
    <taskDetails3>NPP</taskDetails3>
    <taskDetails4>VIIRS</taskDetails4>
  </task>

</InfTkConfig>
"""

def aggregate(fileGlob):
    '''
    Aggregate the files given by the input file glob.
    '''
    pass

