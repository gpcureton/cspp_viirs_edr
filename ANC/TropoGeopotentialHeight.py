#!/usr/bin/env python
# encoding: utf-8
"""
TropoGeopotentialHeight.py

 * DESCRIPTION:  Class to granulate the ViirsAncTropoGeopotentialHeight
 *               data product 

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
sourcename= __name__.split(" ")
LOG = logging.getLogger(sourcename[1])

def ingest():
    '''
    Ingest the ancillary dataset.
    '''
    print "Ingesting the ancillary %s" % (sourcename[1])

def granulate():
    '''
    Granulate the ancillary dataset.
    '''
    print "Granulating the ancillary %s" % (sourcename[1])
