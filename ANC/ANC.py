#!/usr/bin/env python
# encoding: utf-8
"""
ANC.py

 * DESCRIPTION:  Top level class for the ANC module. 

Created by Geoff Cureton on 2013-02-28.
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
LOG = logging.getLogger('ANC')

class ANC():

    def ingest(self,ancBlob):
        '''
        Ingest the ancillary dataset.
        '''
        return getattr(ancBlob,blobDatasetName).astype('float')


    def granulate(self):
        '''
        Granulate the ancillary dataset.
        '''
        pass


    def toFile(self):
        '''
        Copy the granulated ANC dataset to a blob file, and create the appropriate 
        asc file.
        '''
        pass
