#!/usr/bin/env python
# encoding: utf-8
"""
adl_viirs_edr_masks.py

Purpose: Run the VIIRS EDR Masks Controller using ADL 3.1. 

Input:
    * One or more HDF5 VIIRS SDR input files, aggregated or single-granule.
      You need one granule before and one granule after a given granule to process.
    * Static IGBP,LWM,NDVI files.
    * A work directory, typically, empty, in which to unpack the granules and generate the output.
      If the work directory specified does not exist, it will be created.

Output:
    * ADL VIIRS Cloud Mask Intermediate Product (IP) blob files, and associated HDF5 output files
    * ADL VIIRS Active Fires Intermediate Product (IP) blob files, and associated HDF5 output files 
    * Associated .asc metadata files.
    * Unless --debug is specified, all blob and asc files, and the log, perf and ancillary directories 
      will be deleted.

Details:
    * If you input a series of granules, the software will scan the work directory.
      Thus, for N input SDR granules, you may get up to N-2 output EDR granules if 
      they are all contiguous.
    * It is ambiguous to provide several copies of the same granule in the work 
      directory; this will result in an error abort.
    * The unpacker gives each unpacking of a given granule its own.

Preconditions:
    * Requires ADL_HOME, CSPP_RT_ANC_PATH, CSPP_RT_ANC_CACHE_DIR environment variables are set.
    * Requires that any needed LD_LIBRARY_PATH is set.
    * Requires that DSTATICDATA is set.

Optional:
    * Environment variables CSPP_RT_ANC_PATH, CSPP_RT_ANC_CACHE_DIR, if static ancillary data is not 
      placed within ADL_HOME.

Minimum commandline:

    python adl_viirs_edr_masks.py  --input_files=INPUTFILES

where...

    INPUTFILES: The fully qualified path to the input files. May be a directory or a file glob.


Created by Geoff Cureton on 2011-09-30.
Copyright (c) 2011 University of Wisconsin SSEC. All rights reserved.
"""

file_Date = '$Date$'
file_Revision = '$Revision$'
file_Author = '$Author$'
file_HeadURL = '$HeadURL$'
file_Id = '$Id$'

__author__ = 'Geoff Cureton <geoff.cureton@ssec.wisc.edu>'
__version__ = '$Id$'
__docformat__ = 'Epytext'


import os, sys, logging, traceback
from os import path,uname,environ
import string
import re
import uuid
import shlex, subprocess
from subprocess import CalledProcessError, call
from shutil import rmtree
from glob import glob
import optparse as optparse
from time import time
from datetime import datetime,timedelta

from scipy import round_

import numpy as np
from numpy import ma
import copy
from bisect import bisect_left,bisect_right

import ctypes
from numpy.ctypeslib import ndpointer

import tables as pytables
from tables import exceptions as pyEx

import ViirsData

from thermo import rh_to_mr
rh_to_mr_vec = np.vectorize(rh_to_mr)

from NCEPtoBlob import NCEPclass
#from NAAPStoBlob import NAAPSclass
#import pyhdf.SD as SD
from HDF4File import HDF4File

# skim and convert routines for reading .asc metadata fields of interest
import adl_blob
import adl_asc
from adl_asc import skim_dir, contiguous_granule_groups, granule_groups_contain, effective_anc_contains,_eliminate_duplicates,_is_contiguous, RDR_REQUIRED_KEYS, POLARWANDER_REQUIRED_KEYS

# ancillary search and unpacker common routines
# We need [ 'CSPP_RT_HOME', 'ADL_HOME', 'CSPP_RT_ANC_TILE_PATH', 'CSPP_RT_ANC_CACHE_DIR', 'CSPP_RT_ANC_PATH' ] environment 
# variables to be set...
from adl_common import sh, anc_files_needed, link_ancillary_to_work_dir, unpack, env, h5_xdr_inventory, get_return_code, check_env
from adl_common import ADL_HOME, CSPP_RT_ANC_PATH, CSPP_RT_ANC_CACHE_DIR, COMMON_LOG_CHECK_TABLE

# log file scanning
import adl_log

# cache loading and searching
#import adl_anc_retrieval

# N_GEO_Ref fix for ADL
#from adl_geo_ref import write_geo_ref

# post-processing on blob+asc products
#from adl_post_process import repack_products,aggregate_products,add_geo_attribute_to_aggregates
#from adl_post_process import SHORTNAME_2_PRODUCTID

# every module should have a LOG object
sourcename= file_Id.split(" ")
LOG = logging.getLogger(sourcename[1])
from adl_common import configure_logging
from adl_common import _test_logging as test_logging

# locations of executables in ADL
ADL_VIIRS_MASKS_EDR=os.path.abspath(os.path.join(ADL_HOME, 'bin', 'ProEdrViirsMasksController.exe'))
ADL_VIIRS_AEROSOL_EDR=os.path.abspath(os.path.join(ADL_HOME, 'bin', 'ProEdrViirsAerosolController.exe'))

# we're not likely to succeed in processing using geolocation smaller than this many bytes
MINIMUM_SDR_BLOB_SIZE = 81000000

# maximum delay between granule end time and next granule start time to consider them contiguous
MAX_CONTIGUOUS_DELTA=timedelta(seconds = 2)

# Attribute paths for Masks EDR and IP
MOD_GEO_TC_GRANULE_ID_ATTR_PATH = 'Data_Products/VIIRS-MOD-GEO-TC/VIIRS-MOD-GEO-TC_Gran_0/N_Granule_ID'
CM_GRANULE_ID_ATTR_PATH = 'Data_Products/VIIRS-CM-IP/VIIRS-CM-IP_Gran_0/N_Granule_ID'
AF_GRANULE_ID_ATTR_PATH = 'Data_Products/VIIRS-AF-EDR/VIIRS-AF-EDR_Gran_0/N_Granule_ID'

# Attribute paths for Aerosol EDR and IP
AOT_IP_GRANULE_ID_ATTR_PATH = 'Data_Products/VIIRS-Aeros-Opt-Thick-IP/VIIRS-Aeros-Opt-Thick-IP_Gran_0/N_Granule_ID'
AOT_EDR_GRANULE_ID_ATTR_PATH = 'Data_Products/VIIRS-Aeros-EDR/VIIRS-Aeros-EDR_Gran_0/N_Granule_ID'
SUSMAT_EDR_GRANULE_ID_ATTR_PATH = 'Data_Products/VIIRS-SusMat-EDR/VIIRS-SusMat-EDR_Gran_0/N_Granule_ID'

CSPP_RT_ANC_HOME = os.path.abspath(os.getenv('CSPP_RT_ANC_HOME'))



###################################################
#                  Global Data                    #
###################################################

environ['TZ'] = 'UTC'
hexPat = '[\\dA-Fa-f]'

def set_sdr_endian(inputEndianness) :
    global sdrEndian 
    if inputEndianness=='big' :
        sdrEndian = adl_blob.BIG_ENDIAN
    elif inputEndianness=='little' :
        sdrEndian = adl_blob.LITTLE_ENDIAN
    else :
        LOG.error('Invalid value for the VIIRS SDR endianness : %s ' % (inputEndianness))

def set_anc_endian(inputEndianness) :
    global ancEndian 
    if inputEndianness=='big' :
        ancEndian = adl_blob.BIG_ENDIAN
    elif inputEndianness=='little' :
        ancEndian = adl_blob.LITTLE_ENDIAN
    else :
        LOG.error('Invalid value for the VIIRS ancillary endianness : %s ' % (inputEndianness))

# QST type value enumerations to be used with the igbp field
IGBP_dict = {
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

# QST Land Water Mask type value enumerations to be used with the qstlwm field
QSTLWM_list = [
    'EVERGREEN_NEEDLELEAF_FOREST', 'EVERGREEN_BROADLEAF_FORESTS', 'DECIDUOUS_NEEDLELEAF_FORESTS', 'DECIDUOUS_BROADLEAF_FORESTS', 'MIXED_FORESTS', 'CLOSED_SHRUBLANDS', 'OPEN_SHRUBLANDS', 'WOODY_SAVANNAS', 'SAVANNAS', 'GRASSLANDS', 'PERMANENT_WETLANDS', 'CROPLANDS', 'URBAN_AND_BUILDUP_LANDS', 'CROPLAND_NATURAL_VEGETATION', 'SNOW_AND_ICE', 'BARREN_OR_SPARSELY_VEGETATED', 'OCEAN_SEA', 'INLAND_WATER', 'COASTAL_WATER', 'UNCLASSIFIED_LAND', 'FILL_VALUE'
]
QSTLWM_dict = {
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

# Digital Elevation Model (DEM) land sea mask types
DEM_list = ['DEM_SHALLOW_OCEAN','DEM_LAND','DEM_COASTLINE','DEM_SHALLOW_INLAND_WATER','DEM_EPHEMERAL_WATER','DEM_DEEP_INLAND_WATER','DEM_MOD_CONT_OCEAN','DEM_DEEP_OCEAN']
DEM_dict = {
    'DEM_SHALLOW_OCEAN'        : 0,
    'DEM_LAND'                 : 1,
    'DEM_COASTLINE'            : 2,
    'DEM_SHALLOW_INLAND_WATER' : 3,
    'DEM_EPHEMERAL_WATER'      : 4,
    'DEM_DEEP_INLAND_WATER'    : 5,
    'DEM_MOD_CONT_OCEAN'       : 6,
    'DEM_DEEP_OCEAN'           : 7
}

DEM_to_QSTLWM = np.array([
    QSTLWM_dict['OCEAN_SEA'],                    # DEM 0 : QSTLWM 17
    QSTLWM_dict['EVERGREEN_NEEDLELEAF_FOREST'],  # DEM 1 : QSTLWM  1
    QSTLWM_dict['COASTAL_WATER'],                # DEM 2 : QSTLWM 19
    QSTLWM_dict['INLAND_WATER'],                 # DEM 3 : QSTLWM 18
    QSTLWM_dict['SAVANNAS'],                     # DEM 4 : QSTLWM  9
    QSTLWM_dict['OCEAN_SEA'],                    # DEM 5 : QSTLWM 17
    QSTLWM_dict['OCEAN_SEA'],                    # DEM 6 : QSTLWM 17
    QSTLWM_dict['OCEAN_SEA']                     # DEM 7 : QSTLWM 17
],dtype=('uint8'))


def index(a, x):
    '''Locate the leftmost value exactly equal to x'''
    i = bisect_left(a, x)
    if i != len(a) and a[i] == x:
        return i
    raise ValueError

def find_lt(a, x):
    '''Find rightmost value less than x'''
    i = bisect_left(a, x)
    if i:
        return a[i-1]
    raise ValueError

def find_le(a, x):
    '''Find rightmost value less than or equal to x'''
    i = bisect_right(a, x)
    if i:
        return a[i-1]
    raise ValueError

def find_gt(a, x):
    '''Find leftmost value greater than x'''
    i = bisect_right(a, x)
    if i != len(a):
        return a[i]
    raise ValueError

def find_ge(a, x):
    '''Find leftmost item greater than or equal to x'''
    i = bisect_left(a, x)
    if i != len(a):
        return a[i]
    raise ValueError

def toJulianDate(inTime):
    '''
    Takes time in regular yyyymmdd format and returns time string in Julian yyyyddd format.
    '''
    inTime = str(inTime)
    try :
        return time.strftime("%Y%j",time.strptime(inTime, "%Y%m%d"))
    except :
        LOG.error("Incorrect data format (%s). Should conform to yyyymmdd." % (inTime))
        return 1

def fromJulianDate(inTime):
    '''
    Takes time in Julian yyyyddd format and returns time string in regular yyyymmdd format .
    '''
    inTime = str(inTime)
    try :
        return time.strftime("%Y%m%d",time.strptime(inTime,"%Y%j"))
    except :
        LOG.error("Incorrect data format (%s). Should conform to yyyyddd." % (inTime))
        return 1

#----------------------------------------------------------------------------
# Finds the places where the boundary points that will make up a polygon
# cross the dateline.
#
# This method is heavily based on the AltNN NNfind_crossings() method
#----------------------------------------------------------------------------
def findDatelineCrossings(latCrnList,lonCrnList):

    #-------------------------------------------------------------------------
    # NOTE:  This loop will find the place(s) where the boundary crosses 180
    # degrees longitude.  It will also record the index after the crossing
    # for the first two crossings.
    # 
    # NOTE:  Since the last point in the boundary is equal to the first point
    # in the boundary, there is no chance of a crossing between the last
    # and first points.
    #
    # initialize the number of crossings to zero
    # for loop over the boundary
    #    if the longitudes cross the 180 degree line, then
    #       increment the number of crossings
    #       if this is first crossing, then
    #          save the index after the crossing
    #       else if this is the second crossing
    #          save the index after the second crossing
    #       endif
    #    endif
    # end for loop
    #-------------------------------------------------------------------------

    status = 0
    numCrosses = 0

    # For an ascending granule, the corner points are numbered [0,1,3,2], from the southeast
    # corner moving anti-clockwise.

    LOG.debug("latCrnList = %r " % (latCrnList))
    LOG.debug("lonCrnList = %r " % (lonCrnList))

    for idx1,idx2 in zip([1,3,2],[0,1,3]):
        
        # Convert the longitudes to radians, and calculate the 
        # absolute difference
        lon1 = np.radians(lonCrnList[idx1])
        lon2 = np.radians(lonCrnList[idx2])
        lonDiff = np.fabs( lon1 - lon2 )
        
        if ( np.fabs(lonDiff) > np.pi ):

            # We have a crossing, incrememnt the number of crossings
            numCrosses += 1
            
            if(numCrosses == 1):

                # This was the first crossing
                cross1Idx_ = idx1

            elif(numCrosses == 2):

                # This was the second crossing
                cross2Idx_ = idx1

            else :

                # we should never get here
                return -1

    num180Crossings_ = numCrosses

    '''
    # now determine the minimum and maximum latitude
    maxLat_ = latCrnList[0]
    minLat_ = maxLat_

    for idx in [1,3,2]:
        if(latCrnList[idx] > maxLat_):
            # if current lat is bigger than maxLat_, make the current point the
            # maximum
            maxLat_ = latCrnList[idx]

        if(latCrnList[idx] < minLat_):
            # if current lat is smaller than minLat_, make the current point the
            # minimum
            minLat_ = latCrnList[idx]

    return num180Crossings_,minLat_,maxLat_
    '''

    return num180Crossings_

def isDatelineCrossed(latCrnList,lonCrnList):

    dateLineCrossed = []
    ascendingNode = []
    descendingNode = []

    LOG.debug("Determining granule dateline crossings...")

    for granLats,granLons in zip(latCrnList,lonCrnList):
        LOG.debug("Processing granule...")
        isDatelineCrosser = False
        isAscendingNode = False
        isDecendingNode = False

        # TODO : Use the ADL method of determining dateline crossings
        numCrossings = findDatelineCrossings(granLats,granLons)
        LOG.debug("Number of dataline crossings %d" % (numCrossings))

        # Ascending node ? ...
        if (granLats[2] > granLats[0]):
            LOG.debug("Ascending node...")
            isAscendingNode = True
            # Dateline crosser ? ...
            if (granLons[0] < granLons[3]):
                LOG.debug("Dateline crosser...\n")
                isDatelineCrosser = True

        # Descending node ? ...
        if (granLats[0] > granLats[2]):
            LOG.debug("Descending node...")
            isDecendingNode = True
            # Dateline crosser ? ...
            if (granLons[1] < granLons[2]):
                LOG.debug("Dateline crosser...\n")
                isDatelineCrosser = True

        dateLineCrossed.append(isDatelineCrosser)
        ascendingNode.append(isAscendingNode)
        descendingNode.append(isDecendingNode)

    return dateLineCrossed,ascendingNode,descendingNode


def _create_input_file_globs(inputFiles):
    '''
    Determine the correct input file path and globs
    '''
    input_path = os.path.abspath(inputFiles)
    if os.path.isdir(input_path) :
        input_dir = input_path
        input_files = None
    else :
        input_dir = os.path.dirname(input_path)
        input_files = os.path.basename(input_path)

    LOG.debug("input_path = %s" %(input_path))
    LOG.debug("input_dir = %s" %(input_dir))
    LOG.debug("input_files = %s" %(input_files))

    inputGlobs = {"GEO":None,\
                  "MOD":None,\
                  "IMG":None}

    charsToKill = string.ascii_letters + string.digits + "."

    if (input_files is None):
        inputGlobs['GEO'] = 'GMTCO_npp*.h5'
        inputGlobs['MOD'] = 'SVM*_npp*.h5'
        inputGlobs['IMG'] = 'SVI*_npp*.h5'
    elif ((('GMTCO' in input_files) or ('SVM' in input_files) or ('SVI' in input_files)) and ('*' in input_files)) :
        fileGlob = string.rstrip(string.lstrip(input_files,charsToKill),charsToKill)
        LOG.debug("fileGlob = %s" %(fileGlob))
        inputGlobs['GEO'] = "GMTCO%s.h5" %(fileGlob)
        inputGlobs['MOD'] = "SVM*%s.h5" %(fileGlob)
        inputGlobs['IMG'] = "SVI*%s.h5" %(fileGlob)
        for fileType in ['GEO','MOD','IMG']:
            inputGlobs[fileType] = string.replace(inputGlobs[fileType],"**","*")
    elif os.path.isfile(input_path) :
        fileGlob = string.rstrip(string.lstrip(string.split(input_files,"b")[0],charsToKill),charsToKill)
        LOG.debug("fileGlob = %s" %(fileGlob))
        inputGlobs['GEO'] = "GMTCO%s*.h5" %(fileGlob)
        inputGlobs['MOD'] = "SVM*%s*.h5" %(fileGlob)
        inputGlobs['IMG'] = "SVI*%s*.h5" %(fileGlob)
        for fileType in ['GEO','MOD','IMG']:
            inputGlobs[fileType] = string.replace(inputGlobs[fileType],"**","*")

    return input_dir,inputGlobs


def _test_sdr_granules(collectionShortName, work_dir='.'):
    "list granules we'd generate XML for"
    granules_to_process =  list(sift_metadata_for_viirs_sdr(collectionShortName ,work_dir))
    from pprint import pprint
    pprint([x['N_Granule_ID'] for x in granules_to_process])
    return granules_to_process

def _skim_viirs_sdr(collectionShortName,work_dir):
    "skim for VIIRS SDR data meeting minimum requirements"
    for info in skim_dir(work_dir, N_Collection_Short_Name=collectionShortName):
        print info
        path = info['BlobPath']
        if not os.path.isfile(path) or os.stat(path).st_size < MINIMUM_SDR_BLOB_SIZE:
            LOG.info('Ignoring %s, invalid blob. For direct broadcast this can happen at start or end of pass.' % (path))
        else:
            yield info

def _contiguous_granule_groups(granules, tolerance=MAX_CONTIGUOUS_DELTA, larger_granules_preferred=False):
    """
    given a sequence of granule dictionaries, yield a sequence of contiguous granule groups as tuples
    tolerance, if provided, is a datetime.timedelta object representing max differance between endtime and starttime

    This is a custom version of adl_asc.contiguous_granule_groups(), which keys off of 'StartTime', rather than
    'ObservedStartTime' as is done here.
    """
    
    # sort granules into start time order and eliminate exact duplicates
    # FUTURE: is lex-compare sufficient for A2/A1/etc
    #start_time_key = lambda x: (x['StartTime'], x.get('N_Granule_Version', None))
    start_time_key = lambda x: (x['ObservedStartTime'], x.get('N_Granule_Version', None))
    granlist = _eliminate_duplicates(sorted(granules, key = start_time_key),larger_granules_preferred=larger_granules_preferred)
    granset = set(x['N_Granule_ID'] for x in granlist)

    # it's ambiguous if we have a work directory with multiple different blobs for any given granule
    if len(granlist) != len(granset):
        LOG.error('Aborting. Multiple granule copies in work directory: %d files in directory representing only %d granules' % (len(granlist), len(granset)))
        LOG.error('To solve this, only unpack granules once in given work directory')
        return

    # pair off in start-time order
    # build a set containing contiguous granules
    # when we find a break in sequence, yield the set as an ordered tuple and start over
    seq = {}
    for a,b in zip(granlist[:-1], granlist[1:]):
        if a['N_Granule_ID']==b['N_Granule_ID']:
            LOG.error('Granule %r has been unpacked to this directory multiple times!' % (a['N_Granule_ID']))
            return
        if _is_contiguous(a, b, tolerance):
            seq[a['URID']] = a
            seq[b['URID']] = b
        else:
            LOG.info('contiguous sequence has %d granules' % (len(seq)))
            yield tuple(sorted(seq.values(), key=start_time_key))
            seq.clear()
    # leftovers! yum!
    if seq:
        LOG.info('contiguous sequence has %d granules' % (len(seq)))
        yield tuple(sorted(seq.values(), key=start_time_key))


def sift_metadata_for_viirs_sdr(collectionShortName, crossGran=None, work_dir='.'):
    """
    Search through the ASC metadata in a directory, grouping in StartTime order.
    Look for back-to-back granules and break into contiguous sequences.
    Yield sequence of granules we can process.
    """
    LOG.debug('Collecting information for VIIRS SDRs')

    work_dir = os.path.abspath(work_dir)

    geoGroupList = list(_contiguous_granule_groups(skim_dir(work_dir, N_Collection_Short_Name=collectionShortName)))

    if len(geoGroupList)==0:
        LOG.warn('No VIIRS geolocation files were found for collection shortname %s!'%(collectionShortName))

    LOG.debug('Sifting VIIRS SDR data for processing opportunities')
    for group in _contiguous_granule_groups(skim_dir(work_dir, N_Collection_Short_Name=collectionShortName)):
        ##- for VIIRS, we can process everything but the first and last granule
        ##- for CrIS, use [4:-4]
        LOG.debug('contiguous granule group: %r' % (group,))

        if not crossGran :
            startGran,endGran = None,None
        else :
            startGran,endGran = crossGran,-1*crossGran

        for gran in group[startGran:endGran]:
            if not granule_groups_contain(geoGroupList, gran):
                LOG.info("Insufficient VIIRS SDR coverage to process %s @ %s (%s) - skipping" % (gran['N_Granule_ID'], gran['StartTime'], gran['URID']))
                continue
                #pass
            LOG.info('Processing opportunity: %r at %s with uuid %s' % (gran['N_Granule_ID'], gran['StartTime'], gran['URID']))
            yield gran
 
# XML template for ProEdrViirsMasksController.exe
XML_TMPL_VIIRS_MASKS_EDR = """<InfTkConfig>
  <idpProcessName>ProEdrViirsMasksController.exe</idpProcessName>
  <siSoftwareId />
  <isRestart>FALSE</isRestart>
  <useExtSiWriter>FALSE</useExtSiWriter>
  <debugLogLevel>LOW</debugLogLevel>
  <debugLevel>DBG_HIGH</debugLevel>
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
  <actualScans>47</actualScans>
  <previousActualScans>48</previousActualScans>
  <nextActualScans>48</nextActualScans> 
  <usingMetadata>TRUE</usingMetadata>
  <configGuideName>ProEdrViirsMasksController_GuideList.cfg</configGuideName>

  <task>
    <taskType>EDR</taskType>
    <taskDetails1>%(N_Granule_ID)s</taskDetails1>
    <taskDetails2>%(N_Granule_Version)s</taskDetails2>
    <taskDetails3>NPP</taskDetails3>
    <taskDetails4>VIIRS</taskDetails4>
  </task>
</InfTkConfig>
"""

# XML template for ProEdrViirsAerosolController.exe
XML_TMPL_VIIRS_AEROSOL_EDR = """<InfTkConfig>
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

XML_TMPL_VIIRS_MASKS_EDR_ADL41 = """<InfTkConfig>
  <idpProcessName>ProEdrViirsMasksController.exe</idpProcessName>
  <siSoftwareId></siSoftwareId>
  <isRestart>FALSE</isRestart>
  <useExtSiWriter>FALSE</useExtSiWriter>
  <debugLogLevel>NORMAL</debugLogLevel>
  <debugLevel>DBG_LOW</debugLevel>
  <dbgDest>D_FILE</dbgDest>
  <enablePerf>FALSE</enablePerf>
  <perfPath>${ADL_HOME}/log</perfPath>
  <dbgPath>${ADL_HOME}/log</dbgPath>
  <initData>
     <domain>OPS</domain>
     <subDomain>SUBDOMAIN</subDomain>
     <startMode>INF_STARTMODE_COLD</startMode>
     <executionMode>INF_EXEMODE_PRIMARY</executionMode>
     <healthTimeoutPeriod>30</healthTimeoutPeriod>
  </initData>
  <lockinMem>FALSE</lockinMem>
  <rootDir>${WORK_DIR}/log</rootDir>
  <inputPath>${ADL_HOME}/data/input/withMetadata/ProEdrViirsMasksControllerInputs:${WORK_DIR}</inputPath>
  <outputPath>${ADL_HOME}/data/output/withMetadata/ProEdrViirsMasksControllerOutputs</outputPath>
  <dataStartIET>0</dataStartIET>
  <dataEndIET>0</dataEndIET>
  <actualScans>0</actualScans>
  <previousActualScans>0</previousActualScans>
  <nextActualScans>0</nextActualScans> 
  <usingMetadata>TRUE</usingMetadata>
  <configGuideName>ProEdrViirsMasksController_GuideListsrm.cfg</configGuideName>

  <task>
    <taskType>EDR</taskType>
    <taskDetails1>NPP001212025477</taskDetails1>
    <taskDetails2>A1</taskDetails2>
    <taskDetails3>NPP</taskDetails3>
    <taskDetails4>VIIRS</taskDetails4>
  </task>

</InfTkConfig>

"""

def generate_viirs_masks_edr_xml(work_dir, granule_seq):
    "generate XML files for VIIRS Masks EDR granule generation"
    to_process = []
    for gran in granule_seq:
        name = gran['N_Granule_ID']
        fnxml = 'edr_viirs_masks_%s.xml' % (name)
        LOG.debug('writing XML file %r' % (fnxml))
        fpxml = file(os.path.join(work_dir, fnxml), 'wt')
        fpxml.write(XML_TMPL_VIIRS_MASKS_EDR % gran)
        to_process.append([name,fnxml])
    return to_process


def generate_viirs_aerosol_edr_xml(work_dir, granule_seq):
    "generate XML files for VIIRS Masks EDR granule generation"
    to_process = []
    for gran in granule_seq:
        name = gran['N_Granule_ID']
        fnxml = 'edr_viirs_aerosol_%s.xml' % (name)
        LOG.debug('writing XML file %r' % (fnxml))
        fpxml = file(os.path.join(work_dir, fnxml), 'wt')
        fpxml.write(XML_TMPL_VIIRS_AEROSOL_EDR % gran)
        to_process.append([name,fnxml])
    return to_process

def _get_geo_info(inDir,collectionShortName):
    ''' Return a list of dictionaries summarising the metadata for the geolocation files '''

    LOG.info('Geolocation dir: %s' % (inDir))
    LOG.debug('Collection short name is : %s' % (collectionShortName))

    geoDicts = []
    geoGroups = adl_asc.contiguous_granule_groups(adl_asc.skim_dir(inDir,N_Collection_Short_Name=collectionShortName))

    for group in geoGroups :
        for granDict in group :
            geoDicts.append(granDict)

    if (geoDicts == []) :
        group = adl_asc.skim_dir(inDir,N_Collection_Short_Name=collectionShortName)
        while True:
            try :
                thisDict = group.next()
                geoDicts.append(thisDict)
            except StopIteration :
                break

    for dicts in geoDicts :
        for key,item in dicts.items() :
            LOG.info("%s : %s" % (key,repr(item)))

    return geoDicts

def _get_geo_Arrays(geoDicts):

    # Determine the latitude and longitude ranges of the geolocation granules

    latitudeList = []
    longitudeList = []
    latMinList = []
    latMaxList = []
    lonMinList = []
    lonMaxList = []
    latCrnList = []
    lonCrnList = []
    scanModeList = []

    adlHome = os.getenv('ADL_HOME')

    # Moderate resolution trim table arrays. These are 
    # bool arrays, and the trim pixels are set to True.

    trimObj = ViirsData.ViirsTrimTable()
    modTrimMask = trimObj.createModTrimArray(nscans=48,trimType=bool)

    for geoDict in geoDicts :

        URID = geoDict['URID']
        geo_Collection_ShortName = geoDict['N_Collection_Short_Name']
        N_Granule_ID = geoDict['N_Granule_ID']
        ObservedStartTimeObj = geoDict['ObservedStartTime']

        LOG.debug("\n###########################")
        LOG.debug("  Geolocation Information  ")
        LOG.debug("###########################")
        LOG.debug("N_Granule_ID : %r" % (N_Granule_ID))
        LOG.debug("ObservedStartTime : %s" % (ObservedStartTimeObj.__str__()))
        LOG.debug("N_Collection_Short_Name : %s" %(geo_Collection_ShortName))
        LOG.debug("URID : %r" % (URID))
        LOG.debug("###########################\n")



        geoAscFileName = geoDict['_filename']
        geoBlobFileName = string.replace(geoAscFileName,'asc',geo_Collection_ShortName)

        LOG.debug("geolocation asc file :  %s" % (geoAscFileName))
        LOG.debug("geolocation blob file : %s" % (geoBlobFileName))

        # Do we have terrain corrected geolocation?

        terrainCorrectedGeo = True if 'GEO-TC' in geo_Collection_ShortName else False

        # Do we have long or short style geolocation field names?

        if (geo_Collection_ShortName=='VIIRS-MOD-GEO-TC' or geo_Collection_ShortName=='VIIRS-MOD-RGEO') :
            longFormGeoNames = True
            LOG.debug("We have long form geolocation names")
        elif (geo_Collection_ShortName=='VIIRS-MOD-GEO' or geo_Collection_ShortName=='VIIRS-MOD-RGEO-TC') :
            LOG.debug("We have short form geolocation names")
            longFormGeoNames = False
        else :
            LOG.error("Invalid geolocation shortname")
            return -1

        # Get the geolocation xml file

        geoXmlFile = "%s.xml" % (string.replace(geo_Collection_ShortName,'-','_'))
        geoXmlFile = path.join(adlHome,'xml/VIIRS',geoXmlFile)
        if os.path.exists(geoXmlFile):
            LOG.debug("We are using for %s: %s,%s" %(geo_Collection_ShortName,geoXmlFile,geoBlobFileName))

        # Open the geolocation blob and get the latitude and longitude

        endian=sdrEndian

        geoBlobObj = adl_blob.map(geoXmlFile,geoBlobFileName,endian=endian)
        geoBlobArrObj = geoBlobObj.as_arrays()

        if longFormGeoNames :
            latitude = getattr(geoBlobArrObj,'latitude').astype('float')
            longitude = getattr(geoBlobArrObj,'longitude').astype('float')
        else :
            latitude = getattr(geoBlobArrObj,'lat').astype('float')
            longitude = getattr(geoBlobArrObj,'lon').astype('float')
        
        scanMode = getattr(geoBlobArrObj,'scan_mode').astype('uint8')
        LOG.debug("Scan Mode = %r" % (scanMode))

        latitude = ma.masked_less(latitude,-800.)
        latMin,latMax = np.min(latitude),np.max(latitude)
        latRange = latMax-latMin

        longitude = ma.masked_less(longitude,-800.)
        lonMin,lonMax = np.min(longitude),np.max(longitude)
        lonRange = lonMax-lonMin

        LOG.debug("min,max,range of latitide: %f %f %f" % (latMin,latMax,latRange))
        LOG.debug("min,max,range of longitude: %f %f %f" % (lonMin,lonMax,lonRange))

        # Determine the latitude and longitude fill masks, so we can restore the 
        # fill values after we have scaled...

        latMask = latitude.mask
        lonMask = longitude.mask

        # Check if the geolocation is in radians, convert to degrees

        if (lonRange < 2.*np.pi) :
            LOG.info("Geolocation is in radians, convert to degrees...")
            latitude = np.degrees(latitude)
            longitude = np.degrees(longitude)
        
            latMin,latMax = np.min(latitude),np.max(latitude)
            latRange = latMax-latMin
            lonMin,lonMax = np.min(longitude),np.max(longitude)
            lonRange = lonMax-lonMin

            LOG.debug("New min,max,range of latitide: %f %f %f" % (latMin,latMax,latRange))
            LOG.debug("New min,max,range of longitude: %f %f %f" % (lonMin,lonMax,lonRange))

        # Restore fill values to masked pixels in geolocation

        geoFillValue = trimObj.sdrTypeFill['VDNE_FLOAT64_FILL'][latitude.dtype.name]
        latitude = ma.array(latitude,mask=latMask,fill_value=geoFillValue)
        latitude = latitude.filled()

        geoFillValue = trimObj.sdrTypeFill['VDNE_FLOAT64_FILL'][longitude.dtype.name]
        longitude = ma.array(longitude,mask=lonMask,fill_value=geoFillValue)
        longitude = longitude.filled()

        # Shift the longitudes to be between -180 and 180 degrees
        ### FIXME ###
        if lonMax > 180. :
            LOG.debug("\nFinal min,max,range of longitude: %f %f %f" % (lonMin,lonMax,lonRange))
            # Scale to restore -ve longitues, not necessarily # FIXME
            dateLineIdx = np.where(longitude>180.)
            LOG.debug("dateLineIdx = %r" % (dateLineIdx))
            longitude[dateLineIdx] -= 360.
            lonMax = np.max(ma.array(longitude,mask=lonMask))
            lonMin = np.min(ma.array(longitude,mask=lonMask))
            lonRange = lonMax-lonMin
            LOG.debug("\nFinal min,max,range of longitude: %f %f %f" % (lonMin,lonMax,lonRange))
        ### FIXME ###

        latitudeList.append(latitude)
        longitudeList.append(longitude)
        latMinList.append(latMin)
        latMaxList.append(latMax)
        lonMinList.append(lonMin)
        lonMaxList.append(lonMax)
        scanModeList.append(scanMode)

        # Record the corners, taking care to exclude any bad scans...
        firstGoodScan = np.where(scanMode<=2)[0][0]
        lastGoodScan = np.where(scanMode<=2)[0][-1]
        firstGoodRow = firstGoodScan * 16
        lastGoodRow = lastGoodScan * 16 + 15

        latCrnList.append([latitude[firstGoodRow,0],latitude[firstGoodRow,-1],latitude[lastGoodRow,0],latitude[lastGoodRow,-1]])
        lonCrnList.append([longitude[firstGoodRow,0],longitude[firstGoodRow,-1],longitude[lastGoodRow,0],longitude[lastGoodRow,-1]])

    return latitudeList,longitudeList,latMinList,latMaxList,lonMinList,lonMaxList,latCrnList,lonCrnList


def _subset_IGBP(latMinList,latMaxList,lonMinList,lonMaxList,latCrnList,lonCrnList):
    '''Subsets the IGBP global ecosystem dataset to cover the required geolocation range.'''

    #CSPP_RT_ANC_HOME = os.getenv('CSPP_RT_ANC_HOME')

    IGBP_dLat = 60.*(1./3600.)
    IGBP_dLon = 60.*(1./3600.)

    # Get the subset of NDVI global dataset.
    IGBP_fileName = path.join(CSPP_RT_ANC_HOME,'IGBP/IGBP.EcoMap.v1.0.2004.129.v004.h5')

    try :
        IGBPobj = pytables.openFile(IGBP_fileName)
        IGBP_node = IGBPobj.getNode('/IGBP_Land_Cover_Type')
    except :
        LOG.error("Problem opening IGBP file %s"%(IGBP_fileName))
        return -1

    dateLineCrossed,ascendingNode,descendingNode = isDatelineCrossed(latCrnList,lonCrnList)
    LOG.debug("dateLineCross is %r" % (dateLineCrossed))
    LOG.debug("ascendingNode is %r" % (ascendingNode))
    LOG.debug("descendingNode is %r\n" % (descendingNode))

    try :
        IGBP_gridLats = IGBPobj.getNode('/Latitude')[:]
        IGBP_gridLons = IGBPobj.getNode('/Longitude')[:]

        latMin = min(latMinList)
        latMax = max(latMaxList)
        lonMin = min(lonMinList)
        lonMax = max(lonMaxList)

        IGBP_latMask = np.equal((IGBP_gridLats<(latMax+IGBP_dLat)),(IGBP_gridLats>(latMin-IGBP_dLat)))
        IGBP_lonMask = np.equal((IGBP_gridLons<(lonMax+IGBP_dLon)),(IGBP_gridLons>(lonMin-IGBP_dLon)))

        IGBP_latIdx = np.where(IGBP_latMask==True)[0]
        IGBP_lonIdx = np.where(IGBP_lonMask==True)[0]

        IGBP_latMinIdx = IGBP_latIdx[0]
        IGBP_latMaxIdx = IGBP_latIdx[-1]
        IGBP_lonMinIdx = IGBP_lonIdx[0]
        IGBP_lonMaxIdx = IGBP_lonIdx[-1]

        lat_subset = IGBP_gridLats[IGBP_latMinIdx:IGBP_latMaxIdx+1]

        if True in dateLineCrossed :
            posLonCrn = np.min(ma.masked_less_equal(np.array(lonCrnList),0.))
            negLonCrn = np.max(ma.masked_outside(np.array(lonCrnList),-800.,0.))
            posIdx = index(IGBP_gridLons,find_lt(IGBP_gridLons,posLonCrn))
            negIdx = index(IGBP_gridLons,find_gt(IGBP_gridLons,negLonCrn))

            posBlock = IGBP_node[IGBP_latMinIdx:IGBP_latMaxIdx+1,posIdx:]
            negBlock = IGBP_node[IGBP_latMinIdx:IGBP_latMaxIdx+1,:negIdx]

            IGBP_subset = np.concatenate((posBlock,negBlock),axis=1)

            posLons_subset = IGBP_gridLons[posIdx:]
            negLons_subset = IGBP_gridLons[:negIdx]
            lon_subset = np.concatenate((posLons_subset,negLons_subset))

        else :

            IGBP_subset = IGBP_node[IGBP_latMinIdx:IGBP_latMaxIdx+1,IGBP_lonMinIdx:IGBP_lonMaxIdx+1]
            lon_subset = IGBP_gridLons[IGBP_lonMinIdx:IGBP_lonMaxIdx+1]

        IGBP_node.close()
        IGBPobj.close()

    except Exception, err :

        LOG.debug("EXCEPTION: %s" % (err))
        IGBP_node.close()
        IGBPobj.close()

    # Fix erroneous number for water bodies from 0 to 17
    waterIdx = np.where(IGBP_subset==0)
    IGBP_subset[waterIdx] = 17

    return IGBP_subset.astype('uint8'),lat_subset,lon_subset,IGBP_fileName

def _granulate_IGBP(geoDicts,inDir):
    '''Granulates the input IGBP files.'''

    # Moderate resolution trim table arrays. These are 
    # bool arrays, and the trim pixels are set to True.

    trimObj = ViirsData.ViirsTrimTable()
    modTrimMask = trimObj.createModTrimArray(nscans=48,trimType=bool)

    # Get a bunch of information about the geolocation
    latitudeList,longitudeList,latMinList,latMaxList,lonMinList,lonMaxList,latCrnList,lonCrnList = _get_geo_Arrays(geoDicts)

    LOG.debug("\nGranules -->")
    LOG.debug("latMin : %r" % (latMinList))
    LOG.debug("latMax : %r" % (latMaxList))
    LOG.debug("lonMin : %r" % (lonMinList))
    LOG.debug("lonMax : %r" % (lonMaxList))
    LOG.debug("latCrnList : %r" % (latCrnList))
    LOG.debug("lonCrnList : %r\n" % (lonCrnList))

    IGBP_subset,lat_subset,lon_subset,IGBP_fileName = \
            _subset_IGBP(latMinList,latMaxList,lonMinList,lonMaxList,latCrnList,lonCrnList)

    gridLon,gridLat = np.meshgrid(lon_subset,lat_subset[::-1])
    IGBP_subset = IGBP_subset[::-1,:]

    IGBP_type = IGBP_subset.dtype

    IGBP_list = []

    for latitude,longitude,geoDict in zip(latitudeList,longitudeList,geoDicts):

        N_Granule_ID = geoDict['N_Granule_ID']

        LOG.debug("\nGranulating %s ..." % ('IGBP'))
        LOG.debug("latitide,longitude shapes: %s, %s"%(str(latitude.shape) , str(longitude.shape)))
        LOG.debug("IGBP_subset.shape = %s" % (str(IGBP_subset.shape)))
        LOG.debug("gridLat.shape = %s" % (str(gridLat.shape)))
        LOG.debug("gridLon.shape = %s" % (str(gridLon.shape)))

        LOG.debug("min of IGBP_subset  = %s"%(np.min(IGBP_subset)))
        LOG.debug("max of IGBP_subset  = %s"%(np.max(IGBP_subset)))

        data,dataIdx = _grid2Gran(np.ravel(latitude),
                                  np.ravel(longitude),
                                  IGBP_subset.astype(np.float64),
                                  gridLat.astype(np.float64),
                                  gridLon.astype(np.float64))

        data = data.reshape(latitude.shape)
        dataIdx = dataIdx.reshape(latitude.shape)
        LOG.debug("Shape of first granulated %s data is %s" % ('IGBP',np.shape(data)))
        LOG.debug("Shape of first granulated %s dataIdx is %s" % ('IGBP',np.shape(dataIdx)))

        # Convert granulated data back to original type...

        data = data.astype(IGBP_type)

        # Fill the required pixel trim rows in the granulated NCEP data with 
        # the ONBOARD_PT_FILL value for the correct data type

        fillValue = trimObj.sdrTypeFill['ONBOARD_PT_FILL'][data.dtype.name]        
        data = ma.array(data,mask=modTrimMask,fill_value=fillValue)
        LOG.debug("min of IGBP granule = %d"%(np.min(data)))
        LOG.debug("max of IGBP granule = %d"%(np.max(data)))

        data = data.filled()

        IGBP_list.append(data)

    return IGBP_list


def _subset_DEM(latMinList,latMaxList,lonMinList,lonMaxList,latCrnList,lonCrnList):
    '''Subsets the global elevation dataset to cover the required geolocation range.'''

    #CSPP_RT_ANC_HOME = os.getenv('CSPP_RT_ANC_HOME')

    DEM_dLat = 30.*(1./3600.)
    DEM_dLon = 30.*(1./3600.)

    DEM_gridLats = -1. * (np.arange(21600.) * DEM_dLat - 90.)
    DEM_gridLons = np.arange(43200.) * DEM_dLon - 180.

    # Get the subset of DEM global dataset.
    DEM_fileName = path.join(CSPP_RT_ANC_HOME,'LSM/dem30ARC_Global_LandWater_uncompressed.h5')

    dateLineCrossed,ascendingNode,descendingNode = isDatelineCrossed(latCrnList,lonCrnList)
    LOG.debug("dateLineCross is %r" % (dateLineCrossed))
    LOG.debug("ascendingNode is %r" % (ascendingNode))
    LOG.debug("descendingNode is %r\n" % (descendingNode))

    latMin = min(latMinList)
    latMax = max(latMaxList)
    lonMin = min(lonMinList)
    lonMax = max(lonMaxList)

    DEM_latMask = np.equal((DEM_gridLats<(latMax+DEM_dLat)),(DEM_gridLats>(latMin-DEM_dLat)))
    DEM_lonMask = np.equal((DEM_gridLons<(lonMax+DEM_dLon)),(DEM_gridLons>(lonMin-DEM_dLon)))

    DEM_latIdx = np.where(DEM_latMask==True)[0]
    DEM_lonIdx = np.where(DEM_lonMask==True)[0]

    DEM_latMinIdx = DEM_latIdx[0]
    DEM_latMaxIdx = DEM_latIdx[-1]
    DEM_lonMinIdx = DEM_lonIdx[0]
    DEM_lonMaxIdx = DEM_lonIdx[-1]

    LOG.debug("Opening DEM file %s" % (DEM_fileName))
    # TODO : Use original HDF4 file which contains elevation and LWM.
    DEMobj = pytables.openFile(DEM_fileName,'r')
    DEM_node = DEMobj.getNode('/demGRID/Data Fields/LandWater')

    try :

        lat_subset = DEM_gridLats[DEM_latMinIdx:DEM_latMaxIdx+1]

        if True in dateLineCrossed :
            posLonCrn = np.min(ma.masked_less_equal(np.array(lonCrnList),0.))
            negLonCrn = np.max(ma.masked_outside(np.array(lonCrnList),-800.,0.))
            posIdx = index(DEM_gridLons,find_lt(DEM_gridLons,posLonCrn))
            negIdx = index(DEM_gridLons,find_gt(DEM_gridLons,negLonCrn))

            posBlock = DEM_node[DEM_latMinIdx:DEM_latMaxIdx+1,posIdx:]
            negBlock = DEM_node[DEM_latMinIdx:DEM_latMaxIdx+1,:negIdx]

            DEM_subset = np.concatenate((posBlock,negBlock),axis=1)

            posLons_subset = DEM_gridLons[posIdx:]
            negLons_subset = DEM_gridLons[:negIdx]
            lon_subset = np.concatenate((posLons_subset,negLons_subset))

        else :

            DEM_subset = DEM_node[DEM_latMinIdx:DEM_latMaxIdx+1,DEM_lonMinIdx:DEM_lonMaxIdx+1]
            lon_subset = DEM_gridLons[DEM_lonMinIdx:DEM_lonMaxIdx+1]

        DEM_node.close()
        DEMobj.close()

    except Exception, err :

        LOG.debug("EXCEPTION: %s" % (err))

        DEM_node.close()
        DEMobj.close()

    return DEM_subset.astype('uint8'),lat_subset,lon_subset,DEM_fileName


def _granulate_DEM(geoDicts,inDir):
    '''Granulates the input DEM files.'''

    # Moderate resolution trim table arrays. These are 
    # bool arrays, and the trim pixels are set to True.

    trimObj = ViirsData.ViirsTrimTable()
    modTrimMask = trimObj.createModTrimArray(nscans=48,trimType=bool)

    # Get a bunch of information about the geolocation
    latitudeList,longitudeList,latMinList,latMaxList,lonMinList,lonMaxList,latCrnList,lonCrnList = _get_geo_Arrays(geoDicts)

    LOG.debug("\nGranules -->")
    LOG.debug("latMin : %r" % (latMinList))
    LOG.debug("latMax : %r" % (latMaxList))
    LOG.debug("lonMin : %r" % (lonMinList))
    LOG.debug("lonMax : %r" % (lonMaxList))
    LOG.debug("latCrnList : %r" % (latCrnList))
    LOG.debug("lonCrnList : %r\n" % (lonCrnList))

    DEM_subset,lat_subset,lon_subset,DEM_fileName = \
            _subset_DEM(latMinList,latMaxList,lonMinList,lonMaxList,latCrnList,lonCrnList)

    gridLon,gridLat = np.meshgrid(lon_subset,lat_subset[::-1])
    DEM_subset = DEM_subset[::-1,:]

    DEM_type = DEM_subset.dtype

    DEM_list = []

    for latitude,longitude,geoDict in zip(latitudeList,longitudeList,geoDicts):

        N_Granule_ID = geoDict['N_Granule_ID']

        LOG.debug("\nGranulating %s ..." % ('DEM'))
        LOG.debug("latitide,longitude shapes: %s, %s"%(str(latitude.shape) , str(longitude.shape)))
        LOG.debug("DEM_subset.shape = %s" % (str(DEM_subset.shape)))
        LOG.debug("gridLat.shape = %s" % (str(gridLat.shape)))
        LOG.debug("gridLon.shape = %s" % (str(gridLon.shape)))

        LOG.debug("min of DEM_subset  = %s"%(np.min(DEM_subset)))
        LOG.debug("max of DEM_subset  = %s"%(np.max(DEM_subset)))

        data,dataIdx = _grid2Gran(np.ravel(latitude),
                                  np.ravel(longitude),
                                  DEM_subset.astype(np.float64),
                                  gridLat.astype(np.float64),
                                  gridLon.astype(np.float64))

        data = data.reshape(latitude.shape)
        dataIdx = dataIdx.reshape(latitude.shape)
        LOG.debug("Shape of first granulated %s data is %s" % ('DEM',np.shape(data)))
        LOG.debug("Shape of first granulated %s dataIdx is %s" % ('DEM',np.shape(dataIdx)))

        # Convert granulated data back to original type...

        data = data.astype(DEM_type)

        # Fill the required pixel trim rows in the granulated NCEP data with 
        # the ONBOARD_PT_FILL value for the correct data type

        fillValue = trimObj.sdrTypeFill['ONBOARD_PT_FILL'][data.dtype.name]        
        data = ma.array(data,mask=modTrimMask,fill_value=fillValue)
        LOG.debug("min of DEM granule = %d"%(np.min(data)))
        LOG.debug("max of DEM granule = %d"%(np.max(data)))

        data = data.filled()

        DEM_list.append(data)

    return DEM_list


def _subset_NDVI(latMinList,latMaxList,lonMinList,lonMaxList,latCrnList,lonCrnList,inDir,geoDict):
    '''Subsets the global NDVI dataset to cover the required geolocation range.'''

    #CSPP_RT_ANC_HOME = os.getenv('CSPP_RT_ANC_HOME')

    NDVI_dLat = 60.*(1./3600.)
    NDVI_dLon = 60.*(1./3600.)

    # Get the subset of NDVI global dataset.

    NDVIdays = np.array([1, 17, 33, 49, 65, 81, 97, 113, 129, 145, \
                         161, 177, 193, 209, 225, 241, 257, 273, 289, 305, 321, 337, 353])

    startTimeObj = geoDict['ObservedStartTime']
    julianDay = int(startTimeObj.strftime('%j'))
    LOG.debug("Julian day = %d" % (julianDay))

    lowerDay = find_le(NDVIdays,julianDay)
    lowerIdx = index(NDVIdays,lowerDay)
    LOG.debug("lowerdDay, lowerIdx = %d, %d" % (lowerDay,lowerIdx))

    if lowerIdx == 22 :
        NDVIday = NDVIdays[22]
        NDVIidx = 22
        remain = -1
    else :
        #upperIdx = lowerIdx+1
        #remain = (NDVIdays[upperIdx]-1) - julianDay
        #print "remain = ",remain
        #if remain < 8 :
            #NDVIday = NDVIdays[upperIdx]
            #NDVIidx = upperIdx
        #else :
            #NDVIday = NDVIdays[lowerIdx]
            #NDVIidx = lowerIdx
        NDVIday = NDVIdays[lowerIdx]
        NDVIidx = lowerIdx

    NDVI_fileName = path.join(CSPP_RT_ANC_HOME,'NDVI/NDVI.FM.c004.v2.0.WS.00-04.%03d.h5'%(NDVIday))

    LOG.debug("NDVI file : %s" % (NDVI_fileName))

    NDVIobj = pytables.openFile(NDVI_fileName)
    NDVI_node = NDVIobj.getNode('/NDVI')

    dateLineCrossed,ascendingNode,descendingNode = isDatelineCrossed(latCrnList,lonCrnList)
    LOG.debug("dateLineCross is %r" % (dateLineCrossed))
    LOG.debug("ascendingNode is %r" % (ascendingNode))
    LOG.debug("descendingNode is %r\n" % (descendingNode))

    try :
        NDVI_gridLats = NDVIobj.getNode('/Latitude')[:]
        NDVI_gridLons = NDVIobj.getNode('/Longitude')[:]

        NDVI_scaleFactor = NDVI_node.attrs['scale_factor']
        NDVI_offset = NDVI_node.attrs['add_offset']
        NDVI_fillValue = NDVI_node.attrs['_FillValue']

        latMin = min(latMinList)
        latMax = max(latMaxList)
        lonMin = min(lonMinList)
        lonMax = max(lonMaxList)

        NDVI_latMask = np.equal((NDVI_gridLats<(latMax+NDVI_dLat)),(NDVI_gridLats>(latMin-NDVI_dLat)))
        NDVI_lonMask = np.equal((NDVI_gridLons<(lonMax+NDVI_dLon)),(NDVI_gridLons>(lonMin-NDVI_dLon)))

        NDVI_latIdx = np.where(NDVI_latMask==True)[0]
        NDVI_lonIdx = np.where(NDVI_lonMask==True)[0]

        NDVI_latMinIdx = NDVI_latIdx[0]
        NDVI_latMaxIdx = NDVI_latIdx[-1]
        NDVI_lonMinIdx = NDVI_lonIdx[0]
        NDVI_lonMaxIdx = NDVI_lonIdx[-1]

        lat_subset = NDVI_gridLats[NDVI_latMinIdx:NDVI_latMaxIdx+1]

        if True in dateLineCrossed :
            posLonCrn = np.min(ma.masked_less_equal(np.array(lonCrnList),0.))
            negLonCrn = np.max(ma.masked_outside(np.array(lonCrnList),-800.,0.))
            posIdx = index(NDVI_gridLons,find_lt(NDVI_gridLons,posLonCrn))
            negIdx = index(NDVI_gridLons,find_gt(NDVI_gridLons,negLonCrn))

            posBlock = NDVI_node[NDVI_latMinIdx:NDVI_latMaxIdx+1,posIdx:]
            negBlock = NDVI_node[NDVI_latMinIdx:NDVI_latMaxIdx+1,:negIdx]

            NDVI_subset = np.concatenate((posBlock,negBlock),axis=1)

            posLons_subset = NDVI_gridLons[posIdx:]
            negLons_subset = NDVI_gridLons[:negIdx]
            lon_subset = np.concatenate((posLons_subset,negLons_subset))

        else :

            NDVI_subset = NDVI_node[NDVI_latMinIdx:NDVI_latMaxIdx+1,NDVI_lonMinIdx:NDVI_lonMaxIdx+1]
            lon_subset = NDVI_gridLons[NDVI_lonMinIdx:NDVI_lonMaxIdx+1]

        NDVI_mask = ma.masked_equal(NDVI_subset,NDVI_fillValue).mask
        NDVI_subset = NDVI_subset * NDVI_scaleFactor + NDVI_offset
        NDVI_subset = ma.array(NDVI_subset,mask=NDVI_mask,fill_value=-999.9)
        NDVI_subset = NDVI_subset.filled()

        NDVI_node.close()
        NDVIobj.close()

    except Exception, err :

        LOG.debug("EXCEPTION: %s" % (err))

        NDVI_node.close()
        NDVIobj.close()

    return NDVI_subset,lat_subset,lon_subset,NDVI_fileName


def _granulate_NDVI(inDir,geoDicts):
    '''Granulates the input NDVI files.'''

    global ancEndian 

    
    #CSPP_RT_ANC_HOME = os.getenv('CSPP_RT_ANC_HOME')
    ADL_ASC_TEMPLATES = path.join(CSPP_RT_ANC_HOME,'asc_templates')

    masksCollShortNames = 'VIIRS-GridIP-VIIRS-Nbar-Ndvi-Mod-Gran'

    GridIP_shortNameToBlobName = {}
    GridIP_dataType = {}
    GridIP_shortNameToXmlName = {}
    GridIP_shortNameToBlobName['VIIRS-GridIP-VIIRS-Nbar-Ndvi-Mod-Gran'] = 'nbarNdvi'
    GridIP_dataType['VIIRS-GridIP-VIIRS-Nbar-Ndvi-Mod-Gran'] = 'float64'
    GridIP_shortNameToXmlName['VIIRS-GridIP-VIIRS-Nbar-Ndvi-Mod-Gran'] = 'VIIRS_GRIDIP_VIIRS_NBAR_NDVI_MOD_GRAN.xml'

    # Moderate resolution trim table arrays. These are 
    # bool arrays, and the trim pixels are set to True.

    trimObj = ViirsData.ViirsTrimTable()
    modTrimMask = trimObj.createModTrimArray(nscans=48,trimType=bool)

    # Get a bunch of information about the geolocation
    latitudeList,longitudeList,latMinList,latMaxList,lonMinList,lonMaxList,latCrnList,lonCrnList = _get_geo_Arrays(geoDicts)

    LOG.debug("\nGranules -->")
    LOG.debug("latMin : %r" % (latMinList))
    LOG.debug("latMax : %r" % (latMaxList))
    LOG.debug("lonMin : %r" % (lonMinList))
    LOG.debug("lonMax : %r" % (lonMaxList))
    LOG.debug("latCrnList : %r" % (latCrnList))
    LOG.debug("lonCrnList : %r\n" % (lonCrnList))

    NDVI_subset,lat_subset,lon_subset,NDVI_fileName = \
            _subset_NDVI(latMinList,latMaxList,lonMinList,lonMaxList,latCrnList,lonCrnList,inDir,geoDicts[0])

    gridLon,gridLat = np.meshgrid(lon_subset,lat_subset[::-1])
    NDVI_subset = NDVI_subset[::-1,:]

    for latitude,longitude,geoDict in zip(latitudeList,longitudeList,geoDicts):

        N_Granule_ID = geoDict['N_Granule_ID']

        LOG.debug("\nGranulating %s ..." % (masksCollShortNames))
        LOG.debug("latitide,longitude shapes: %s, %s"%(str(latitude.shape) , str(longitude.shape)))
        LOG.debug("NDVI_subset.shape = %s" % (str(NDVI_subset.shape)))
        LOG.debug("gridLat.shape = %s" % (str(gridLat.shape)))
        LOG.debug("gridLon.shape = %s" % (str(gridLon.shape)))

        LOG.debug("min of NDVI_subset = %r"%(np.min(NDVI_subset)))
        LOG.debug("max of NDVI_subset = %r"%(np.max(NDVI_subset)))

        data,dataIdx = _grid2Gran(np.ravel(latitude),
                                  np.ravel(longitude),
                                  NDVI_subset.astype(np.float64),
                                  gridLat.astype(np.float64),
                                  gridLon.astype(np.float64))

        data = data.reshape(latitude.shape)
        dataIdx = dataIdx.reshape(latitude.shape)
        LOG.debug("Shape of first granulated %s data is %s" % (masksCollShortNames,np.shape(data)))
        LOG.debug("Shape of first granulated %s dataIdx is %s" % (masksCollShortNames,np.shape(dataIdx)))

        # Fill the required pixel trim rows in the granulated NCEP data with 
        # the ONBOARD_PT_FILL value for the correct data type

        fillValue = trimObj.sdrTypeFill['ONBOARD_PT_FILL'][data.dtype.name]        
        data = ma.array(data,mask=modTrimMask,fill_value=fillValue)
        LOG.debug("min of NDVI granule = %d"%(np.min(data)))
        LOG.debug("max of NDVI granule = %d"%(np.max(data)))

        data = data.filled()

        # Create new NDVI ancillary blob, and copy granulated data to it

        endian = ancEndian
        xmlName = path.join(ADL_HOME,'xml/VIIRS',GridIP_shortNameToXmlName[masksCollShortNames])

        # Create a new URID to be used in making the asc filenames

        URID_dict = _getURID()

        URID = URID_dict['URID']
        creationDate_nousecStr = URID_dict['creationDate_nousecStr']
        creationDateStr = URID_dict['creationDateStr']

        # Create a new directory in the input directory for the new ancillary
        # asc and blob files

        blobDir = inDir

        ascFileName = path.join(blobDir,URID+'.asc')
        blobName = path.join(blobDir,URID+'.'+masksCollShortNames)

        LOG.debug("ascFileName : %s" % (ascFileName))
        LOG.debug("blobName : %s" % (blobName))

        # Create a new ancillary blob, and copy the data to it.

        newNDVIblobObj = adl_blob.create(xmlName, blobName, endian=endian, overwrite=True)
        newNDVIblobArrObj = newNDVIblobObj.as_arrays()

        blobData = getattr(newNDVIblobArrObj,GridIP_shortNameToBlobName[masksCollShortNames])
        blobData[:,:] = data[:,:]

        # Make a new NDVI asc file from the template, and substitute for the various tags

        ascTemplateFileName = path.join(ADL_ASC_TEMPLATES,"VIIRS-ANC_Template.asc")

        LOG.debug("Creating new asc file\n%s\nfrom template\n%s" % (ascFileName,ascTemplateFileName))
        
        ANC_fileList = [NDVI_fileName]
        for idx in range(len(ANC_fileList)) :
            ANC_fileList[idx] = path.basename(ANC_fileList[idx])
        ANC_fileList.sort()
        ancGroupRecipe = '    ("N_Anc_Filename" STRING EQ "%s")'
        ancFileStr = "%s" % ("\n").join([ancGroupRecipe % (str(files)) for files in ANC_fileList])

        # Parse the geolocation asc file to get struct information which will be 
        # written to the ancillary asc files

        geoAscFileName = path.join(inDir,geoDict['URID']+".asc")
        LOG.debug("\nOpening %s..." % (geoAscFileName))

        geoAscFile = open(geoAscFileName,'rt')

        #RangeDateTimeStr =  _getAscLine(geoAscFile,"RangeDateTime")
        RangeDateTimeStr =  _getAscLine(geoAscFile,"ObservedDateTime")
        RangeDateTimeStr =  string.replace(RangeDateTimeStr,"ObservedDateTime","RangeDateTime")
        GRingLatitudeStr =  _getAscStructs(geoAscFile,"GRingLatitude",12)
        GRingLongitudeStr =  _getAscStructs(geoAscFile,"GRingLongitude",12)

        geoAscFile.close()

        ascTemplateFile = open(ascTemplateFileName,"rt") # Open template file for reading
        ascFile = open(ascFileName,"wt") # create a new text file

        LOG.debug("Template file %s is %r with mode %s" %(ascTemplateFileName,'not open' if ascTemplateFile.closed else 'open',ascTemplateFile.mode))

        LOG.debug("New file %s is %r with mode %s" %(ascFileName,'not open' if ascFile.closed else 'open',ascFile.mode))

        for line in ascTemplateFile.readlines():
           line = line.replace("CSPP_URID",URID)
           line = line.replace("CSPP_CREATIONDATETIME_NOUSEC",creationDate_nousecStr)
           line = line.replace("CSPP_ANC_BLOB_FULLPATH",blobName)
           line = line.replace("CSPP_ANC_COLLECTION_SHORT_NAME",masksCollShortNames)
           line = line.replace("CSPP_GRANULE_ID",N_Granule_ID)
           line = line.replace("CSPP_CREATIONDATETIME",creationDateStr)
           line = line.replace("  CSPP_RANGE_DATE_TIME",RangeDateTimeStr)
           line = line.replace("  CSPP_GRINGLATITUDE",GRingLatitudeStr)
           line = line.replace("  CSPP_GRINGLONGITUDE",GRingLongitudeStr)
           line = line.replace("    CSPP_ANC_SOURCE_FILES",ancFileStr)
           ascFile.write(line) 

        ascFile.close()
        ascTemplateFile.close()


def _getURID() :
    '''
    Create a new URID to be used in making the asc filenames
    '''
    
    URID_dict = {}

    URID_timeObj = datetime.utcnow()
    
    creationDateStr = URID_timeObj.strftime("%Y-%m-%d %H:%M:%S.%f")
    creationDate_nousecStr = URID_timeObj.strftime("%Y-%m-%d %H:%M:%S.000000")
    
    tv_sec = int(URID_timeObj.strftime("%s"))
    tv_usec = int(URID_timeObj.strftime("%f"))
    hostId_ = uuid.getnode()
    thisAddress = id(URID_timeObj)
    
    l = tv_sec + tv_usec + hostId_ + thisAddress
    
    URID = '-'.join( ('{0:08x}'.format(tv_sec)[:8],
                      '{0:05x}'.format(tv_usec)[:5],
                      '{0:08x}'.format(hostId_)[:8],
                      '{0:08x}'.format(l)[:8]) )
    
    URID_dict['creationDateStr'] = creationDateStr
    URID_dict['creationDate_nousecStr'] = creationDate_nousecStr
    URID_dict['tv_sec'] = tv_sec
    URID_dict['tv_usec'] = tv_usec
    URID_dict['hostId_'] = hostId_
    URID_dict['thisAddress'] = thisAddress
    URID_dict['URID'] = URID
    
    return URID_dict


def _setupAuxillaryFiles(inDir):
    '''
    Create asc files for the various auxillary files (tunable parameters etc...) 
    in the inDir, and create links to the binary auxillary files in inDir.
    '''

    CSPP_RT_HOME = os.getenv('CSPP_RT_HOME')
    #CSPP_RT_ANC_HOME = os.getenv('CSPP_RT_ANC_HOME')
    CSPP_RT_ANC_CACHE_DIR = os.getenv('CSPP_RT_ANC_CACHE_DIR')

    ANC_SCRIPTS_PATH = path.join(CSPP_RT_HOME,'viirs')
    ADL_ASC_TEMPLATES = path.join(CSPP_RT_ANC_HOME,'asc_templates')

    

    #auxillaryCollShortNames = ['VIIRS-CM-IP-AC-Int','VIIRS-AF-EDR-AC-Int','VIIRS-Aeros-EDR-AC-Int','NAAPS-ANC-Int','AOT-ANC','VIIRS-AOT-LUT','VIIRS-AOT-Sunglint-LUT','VIIRS-AF-EDR-DQTT','VIIRS-Aeros-EDR-DQTT','VIIRS-SusMat-EDR-DQTT']
    #auxillaryAscTemplateFile = ['VIIRS-CM-IP-AC-Template.asc','VIIRS-AF-EDR-AC-Template.asc','VIIRS-Aeros-EDR-AC-Template.asc','NAAPS-ANC-Inc-Template.asc','AOT-ANC-Template.asc',']
    #auxillaryBlobTemplateFile = ['template.VIIRS-CM-IP-AC','template.VIIRS-AF-EDR-AC','template.VIIRS-Aeros-EDR-AC','template.NAAPS-ANC-Int','template.AOT-ANC']
    #auxillaryPaths = ['ViirsEdrMasks_Aux','ViirsEdrMasks_Aux','ViirsEdrMasks_Aux','NAAPS-ANC-Int','ViirsEdrMasks_Aux']

    #auxillaryCollShortNames = ['VIIRS-CM-IP-AC',
                               #'VIIRS-AF-EDR-AC',
                               #'VIIRS-AF-EDR-DQTT',
                               #'VIIRS-Aeros-EDR-AC',
                               #'VIIRS-Aeros-EDR-DQTT',
                               #'NAAPS-ANC-Int',
                               #'AOT-ANC',
                               #'VIIRS-AOT-LUT',
                               #'VIIRS-AOT-Sunglint-LUT',
                               #'VIIRS-SusMat-EDR-DQTT']

    auxillaryCollShortNames = ['VIIRS-CM-IP-AC-Int',
                               'VIIRS-AF-EDR-AC-Int',
                               'VIIRS-AF-EDR-DQTT-Int',
                               'VIIRS-Aeros-EDR-AC-Int',
                               'VIIRS-Aeros-EDR-DQTT-Int',
                               'NAAPS-ANC-Int',
                               'AOT-ANC',
                               'VIIRS-AOT-LUT',
                               'VIIRS-AOT-Sunglint-LUT',
                               'VIIRS-SusMat-EDR-DQTT-Int']

    auxillaryAscTemplateFile = ['VIIRS-CM-IP-AC-Int_Template.asc',
                                'VIIRS-AF-EDR-AC-Int_Template.asc',
                                'VIIRS-AF-EDR-DQTT-Int_Template.asc',
                                'VIIRS-Aeros-EDR-AC-Int_Template.asc',
                                'VIIRS-Aeros-EDR-DQTT-Int_Template.asc',
                                'NAAPS-ANC-Int_Template.asc',
                                'AOT-ANC_Template.asc',
                                'VIIRS-AOT-LUT_Template.asc',
                                'VIIRS-AOT-Sunglint-LUT_Template.asc',
                                'VIIRS-SusMat-EDR-DQTT-Int_Template.asc']

    auxillaryBlobTemplateFile = ['template.VIIRS-CM-IP-AC-Int',
                                 'template.VIIRS-AF-EDR-AC-Int',
                                 'template.VIIRS-AF-EDR-DQTT-Int',
                                 'template.VIIRS-Aeros-EDR-AC-Int',
                                 'template.VIIRS-Aeros-EDR-DQTT-Int',
                                 'template.NAAPS-ANC-Int',
                                 'template.AOT-ANC',
                                 'template.VIIRS-AOT-LUT',
                                 'template.VIIRS-AOT-Sunglint-LUT',
                                 'template.VIIRS-SusMat-EDR-DQTT-Int']

    auxillaryPaths = ['ViirsEdrMasks_Aux',
                      'ViirsEdrMasks_Aux',
                      'ViirsEdrMasks_Aux',
                      'ViirsEdrMasks_Aux',
                      'ViirsEdrMasks_Aux',
                      'NAAPS-ANC-Int',
                      'ViirsEdrMasks_Aux',
                      'ViirsEdrMasks_Aux',
                      'ViirsEdrMasks_Aux',
                      'ViirsEdrMasks_Aux']


    auxillarySourceFiles = []
    # Number of characters in a URID, plus the trailing "."...
    charsInUrid = 32+1

    for templatePath,blobTempFileName in zip(auxillaryPaths,auxillaryBlobTemplateFile) :
        blobTempFileName = path.join(CSPP_RT_ANC_HOME,templatePath,blobTempFileName)
        if os.path.islink(blobTempFileName) :
            #auxillarySourceFile = string.split(os.path.basename(os.readlink(blobTempFileName)),'.')[1]
            auxillarySourceFile = os.path.basename(os.readlink(blobTempFileName))[charsInUrid:]
            auxillarySourceFiles.append(auxillarySourceFile)
        else :
            auxillarySourceFile = os.path.basename(blobTempFileName)
            auxillarySourceFiles.append(auxillarySourceFile)

    for shortName,auxillarySourceFile in zip(auxillaryCollShortNames,auxillarySourceFiles) :
        LOG.debug("%s --> %s" %(shortName,auxillarySourceFile))



    for shortName,ascTempFileName,blobTempFileName,templatePath,auxillarySourceFile in zip(auxillaryCollShortNames,auxillaryAscTemplateFile,auxillaryBlobTemplateFile,auxillaryPaths,auxillarySourceFiles):

        LOG.info("Creating new %s asc file from template %s" % (shortName,ascTempFileName))

        # Create a new URID to be used in making the asc filenames

        URID_dict = _getURID()

        URID = URID_dict['URID']
        creationDate_nousecStr = URID_dict['creationDate_nousecStr']
        creationDateStr = URID_dict['creationDateStr']

        # The names for the new asc and blob files
        ascFileName = path.join(inDir,URID+'.asc')
        blobFileName = path.join(inDir,string.replace(blobTempFileName,'template',URID))

        # Make a new asc file from the template, and substitute for the various tags

        ascTempFileName = path.join(ADL_ASC_TEMPLATES,ascTempFileName)

        LOG.info("Creating new asc file\n%s\nfrom template\n%s" % (ascFileName,ascTempFileName))

        ascTemplateFile = open(ascTempFileName,"rt") # Open template file for reading
        ascFile = open(ascFileName,"wt") # create a new text file

        LOG.info("Template file %s is %r with mode %s" %(ascTempFileName,'not open' if ascTemplateFile.closed else 'open',ascTemplateFile.mode))

        LOG.info("New file %s is %r with mode %s" %(ascFileName,'not open' if ascFile.closed else 'open',ascFile.mode))

        for line in ascTemplateFile.readlines():
           line = line.replace("CSPP_URID",URID)
           line = line.replace("CSPP_CREATIONDATETIME_NOUSEC",creationDate_nousecStr)
           line = line.replace("CSPP_AUX_BLOB_FULLPATH",blobFileName)
           line = line.replace("CSPP_CREATIONDATETIME",creationDateStr)
           line = line.replace("CSPP_AUX_SOURCE_FILE",auxillarySourceFile)
           ascFile.write(line) 

        ascFile.close()
        ascTemplateFile.close()

        # Create a link between the binary template file and working directory

        blobTempFileName = path.join(CSPP_RT_ANC_HOME,templatePath,blobTempFileName)
        LOG.info("Creating the link %s -> %s" %(blobFileName,blobTempFileName))

        if not os.path.exists(blobFileName):
            LOG.debug('%r -> %r' % (blobFileName, blobTempFileName))
            os.symlink(blobTempFileName, blobFileName)
        else:
            LOG.info('%r already exists; continuing' % blobFileName)
        try:
            LOG.debug('testing %r' % blobFileName)
            s = os.stat(blobFileName)
        except OSError as oops:
            LOG.error("link at %r is broken" % blobFileName)
            raise


def _getGRC(inDir,geoDicts):
    '''
    Setup three dummy GRC files with the appropriate asc files
    '''

    CSPP_RT_HOME = os.getenv('CSPP_RT_HOME')
    #CSPP_RT_ANC_HOME = os.getenv('CSPP_RT_ANC_HOME')
    CSPP_RT_ANC_CACHE_DIR = os.getenv('CSPP_RT_ANC_CACHE_DIR')

    ANC_SCRIPTS_PATH = path.join(CSPP_RT_HOME,'viirs')
    ADL_ASC_TEMPLATES = path.join(CSPP_RT_ANC_HOME,'asc_templates')

    

    # Get a bunch of information about the geolocation
    latitudeList,longitudeList,latMinList,latMaxList,lonMinList,lonMaxList,latCrnList,lonCrnList = _get_geo_Arrays(geoDicts)

    LOG.debug("\nGranules -->")
    LOG.debug("latMin : %r" % (latMinList))
    LOG.debug("latMax : %r" % (latMaxList))
    LOG.debug("lonMin : %r" % (lonMinList))
    LOG.debug("lonMax : %r" % (lonMaxList))
    LOG.debug("latCrnList : %r" % (latCrnList))
    LOG.debug("lonCrnList : %r\n" % (lonCrnList))

    masksCollShortNames = ['VIIRS-MOD-GRC-TC', 'VIIRS-MOD-GRC']
    xmlNames = ['VIIRS_MOD_GRC_TC.xml', 'VIIRS_MOD_GRC.xml']

    for shortName,xmlName in zip(masksCollShortNames,xmlNames):

        for latitude,longitude,geoDict in zip(latitudeList,longitudeList,geoDicts):

            N_Granule_ID = geoDict['N_Granule_ID']
            LOG.info("Creating new %s asc file for granule ID %s" % (shortName,N_Granule_ID))

            LOG.debug("geoDict keys are... %r" % (geoDict.keys()))

            # Create a new URID to be used in making the asc filenames

            URID_dict = _getURID()

            URID = URID_dict['URID']
            creationDate_nousecStr = URID_dict['creationDate_nousecStr']
            creationDateStr = URID_dict['creationDateStr']

            blobDir = inDir

            ascFileName = path.join(blobDir,URID+'.asc')
            blobName = path.join(blobDir,URID+'.'+shortName)

            LOG.info("ascFileName : %s" % (ascFileName))
            LOG.info("blobName : %s" % (blobName))

            # Create new VIIRS-MOD-GRC blob

            endian = ancEndian
            xmlPath = path.join(ADL_HOME,'xml/VIIRS',xmlName)
            #xmlName = path.join(ADL_HOME,'xml/VIIRS','VIIRS_MOD_GRC.xml')
            #xmlName = path.join(ADL_HOME,'xml/VIIRS','VIIRS_MOD_GRC_TC.xml')

            newGRCblobObj = adl_blob.create(xmlPath, blobName, endian=endian, overwrite=True)

            # Make a new VIIRS-MOD-GRC asc file from the template, and substitute for the various tags

            ascTemplateFileName = path.join(ADL_ASC_TEMPLATES,"VIIRS-MOD-GRC_Template.asc")

            LOG.debug("Creating new asc file\n%s\nfrom template\n%s" % (ascFileName,ascTemplateFileName))

            # Parse the geolocation asc file to get struct information which will be 
            # written to the ancillary asc files

            geoAscFileName = path.join(inDir,geoDict['URID']+".asc")
            LOG.debug("\nOpening %s..." % (geoAscFileName))

            geoAscFile = open(geoAscFileName,'rt')

            ObservedDateTime =  _getAscLine(geoAscFile,"ObservedDateTime")
            #RangeDateTimeStr =  _getAscLine(geoAscFile,"RangeDateTime")
            RangeDateTimeStr =  ObservedDateTime
            RangeDateTimeStr =  string.replace(RangeDateTimeStr,"ObservedDateTime","RangeDateTime")

            GRingLatitudeStr =  _getAscStructs(geoAscFile,"GRingLatitude",12)
            GRingLongitudeStr =  _getAscStructs(geoAscFile,"GRingLongitude",12)

            BeginningOrbitNumber =  _getAscLine(geoAscFile,"BeginningOrbitNumber")

            N_Nadir_Latitude_Max  = _getAscLine(geoAscFile,"N_Nadir_Latitude_Max")
            N_Nadir_Latitude_Min  = _getAscLine(geoAscFile,"N_Nadir_Latitude_Min")
            N_Nadir_Longitude_Max = _getAscLine(geoAscFile,"N_Nadir_Longitude_Max")
            N_Nadir_Longitude_Min = _getAscLine(geoAscFile,"N_Nadir_Longitude_Min")

            North_Bounding_Coordinate = _getAscLine(geoAscFile,"North_Bounding_Coordinate")
            South_Bounding_Coordinate = _getAscLine(geoAscFile,"South_Bounding_Coordinate")
            East_Bounding_Coordinate  = _getAscLine(geoAscFile,"East_Bounding_Coordinate")
            West_Bounding_Coordinate  = _getAscLine(geoAscFile,"West_Bounding_Coordinate")

            N_Day_Night_Flag  = _getAscLine(geoAscFile,"N_Day_Night_Flag")
            Ascending_Descending_Indicator  = _getAscLine(geoAscFile,"Ascending/Descending_Indicator")

            geoAscFile.close()

            ascTemplateFile = open(ascTemplateFileName,"rt") # Open template file for reading
            ascFile = open(ascFileName,"wt") # create a new text file

            LOG.debug("Template file %s is %r with mode %s" %(ascTemplateFileName,'not open' if ascTemplateFile.closed else 'open',ascTemplateFile.mode))

            LOG.debug("New file %s is %r with mode %s" %(ascFileName,'not open' if ascFile.closed else 'open',ascFile.mode))

            for line in ascTemplateFile.readlines():
               line = line.replace("CSPP_URID",URID)
               line = line.replace("CSPP_CREATIONDATETIME_NOUSEC",creationDate_nousecStr)
               line = line.replace("CSPP_ANC_BLOB_FULLPATH",blobName)
               line = line.replace("CSPP_ANC_COLLECTION_SHORT_NAME",shortName)
               line = line.replace("CSPP_GRANULE_ID",N_Granule_ID)
               line = line.replace("CSPP_BEGINNING_ORBIT_NUMBER",BeginningOrbitNumber)
               line = line.replace("CSPP_CREATIONDATETIME",creationDateStr)
               line = line.replace("  CSPP_OBSERVED_DATE_TIME",ObservedDateTime)
               line = line.replace("  CSPP_RANGE_DATE_TIME",RangeDateTimeStr)
               line = line.replace("  CSPP_GRINGLATITUDE",GRingLatitudeStr)
               line = line.replace("  CSPP_GRINGLONGITUDE",GRingLongitudeStr)
               #line = line.replace("    CSPP_RT_ANC_SOURCE_FILES",ancFileStr)

               line = line.replace("  CSPP_NADIR_LATITUDE_MAX",N_Nadir_Latitude_Max)
               line = line.replace("  CSPP_NADIR_LATITUDE_MIN",N_Nadir_Latitude_Min)
               line = line.replace("  CSPP_NADIR_LONGITUDE_MAX",N_Nadir_Longitude_Max)
               line = line.replace("  CSPP_NADIR_LONGITUDE_MIN",N_Nadir_Longitude_Min)

               line = line.replace("  CSPP_NORTH_BOUNDING_COORD",North_Bounding_Coordinate)
               line = line.replace("  CSPP_SOUTH_BOUNDING_COORD",South_Bounding_Coordinate)
               line = line.replace("  CSPP_EAST_BOUNDING_COORD",East_Bounding_Coordinate)
               line = line.replace("  CSPP_WEST_BOUNDING_COORD",West_Bounding_Coordinate)

               line = line.replace("  CSPP_DAY_NIGHT_FLAG",N_Day_Night_Flag)
               line = line.replace("  CSPP_ASCENDING_DESCENDING_INDICATOR",Ascending_Descending_Indicator)

               ascFile.write(line) 

            ascFile.close()
            ascTemplateFile.close()


def _QSTLWM(LWM_list,IGBP_list,geoDicts,inDir):
    '''
    Combines the LWM and IGBP to make the QSTLWM
    '''
    
    '''
    LWM  = np.array([ 0, 1, 2, 3, 4, 5, 6, 7, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],dtype='uint8')
    IGBP = np.array([17,17,17,17,17,17,17,17, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16],dtype='uint8')
    QSTLWM_exp = np.array([17,19,19,18, 9,17,17,17, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16],dtype='uint8')

    '''

    global ancEndian 

    
    #CSPP_RT_ANC_HOME = os.getenv('CSPP_RT_ANC_HOME')
    ADL_ASC_TEMPLATES = path.join(CSPP_RT_ANC_HOME,'asc_templates')

    masksCollShortNames = 'VIIRS-GridIP-VIIRS-Qst-Lwm-Mod-Gran'

    GridIP_shortNameToBlobName = {}
    GridIP_dataType = {}
    GridIP_shortNameToXmlName = {}
    GridIP_shortNameToBlobName[masksCollShortNames] = 'qstlwm'
    GridIP_dataType[masksCollShortNames] = 'uint8'
    GridIP_shortNameToXmlName[masksCollShortNames] = 'VIIRS_GRIDIP_VIIRS_QST_LWM_MOD_GRAN.xml'

    # Moderate resolution trim table arrays. These are 
    # bool arrays, and the trim pixels are set to True.

    trimObj = ViirsData.ViirsTrimTable()
    modTrimMask = trimObj.createModTrimArray(nscans=48,trimType=bool)


    for LWM,IGBP,geoDict in zip(LWM_list,IGBP_list,geoDicts):

        N_Granule_ID = geoDict['N_Granule_ID']

        LOG.debug("\nMaking %s ..." % (masksCollShortNames))

        LOG.debug("min of LWM  = %d" % (np.min(LWM)))
        LOG.debug("max of LWM  = %d" % (np.max(LWM)))
        LOG.debug("min of IGBP = %d" % (np.min(IGBP)))
        LOG.debug("max of IGBP = %d" % (np.max(IGBP)))

        # Construct the IGBP and LWM masks

        IGBP_mask = ma.masked_greater(IGBP,17).mask
        LWM_mask = ma.masked_greater(LWM,7).mask

        # Determine the combined mask

        totalMask = np.ravel(np.bitwise_or(IGBP_mask,LWM_mask))

        # Mask and compress the IGBP and LWM datasets

        IGBP_reduced = ma.array(np.ravel(IGBP),mask=totalMask).compressed()
        LWM_reduced = ma.array(np.ravel(LWM),mask=totalMask).compressed()

        DEM_LAND_mask = (LWM_reduced == DEM_dict['DEM_LAND'])
        DEM_water_mask = np.bitwise_not(DEM_LAND_mask)
        IGBP_waterBodies_mask = (IGBP_reduced == IGBP_dict['IGBP_WATERBODIES'])
        
        QSTLWM_coastalWater_mask = np.logical_and(DEM_LAND_mask,IGBP_waterBodies_mask)
        QSTLWM_IGBP_land_mask = np.invert(IGBP_waterBodies_mask)
        
        QSTLWM_IGBP_idx = np.where(QSTLWM_IGBP_land_mask)
        QSTLWM_coastalWater_idx = np.where(QSTLWM_coastalWater_mask)
        QSTLWM_LSM_idx = np.where(DEM_water_mask)
        
        QSTLWM_reduced = np.ones(IGBP_reduced.shape,dtype='uint8') * QSTLWM_dict['FILL_VALUE']
        
        QSTLWM_reduced[QSTLWM_IGBP_idx] = IGBP_reduced[QSTLWM_IGBP_idx]
        QSTLWM_reduced[QSTLWM_coastalWater_idx] = np.array([QSTLWM_dict['COASTAL_WATER']],dtype='uint8')[0]
        QSTLWM_reduced[QSTLWM_LSM_idx] = DEM_to_QSTLWM[LWM_reduced[QSTLWM_LSM_idx]]

        # Restore to original shape...
        QSTLWM = np.ones(IGBP.shape,dtype='uint8') * QSTLWM_dict['FILL_VALUE']
        QSTLWM = np.ravel(QSTLWM)
        validIdx = np.where(totalMask==False)
        QSTLWM[validIdx] = QSTLWM_reduced
        QSTLWM = QSTLWM.reshape(IGBP.shape)

        # Fill the required pixel trim rows in the granulated NCEP data with 
        # the ONBOARD_PT_FILL value for the correct data type

        fillValue = trimObj.sdrTypeFill['ONBOARD_PT_FILL'][QSTLWM.dtype.name]        
        QSTLWM = ma.array(QSTLWM,mask=modTrimMask,fill_value=fillValue)
        LOG.debug("min of QSTLWM granule = %d" % (np.min(QSTLWM)))
        LOG.debug("max of QSTLWM granule = %d" % (np.max(QSTLWM)))
        QSTLWM = QSTLWM.filled()

        # Create new NCEP ancillary blob, and copy granulated data to it

        endian = ancEndian
        xmlName = path.join(ADL_HOME,'xml/VIIRS',GridIP_shortNameToXmlName[masksCollShortNames])

        # Create a new URID to be used in making the asc filenames

        URID_dict = _getURID()

        URID = URID_dict['URID']
        creationDate_nousecStr = URID_dict['creationDate_nousecStr']
        creationDateStr = URID_dict['creationDateStr']

        # Create a new directory in the input directory for the new ancillary
        # asc and blob files

        blobDir = inDir

        ascFileName = path.join(blobDir,URID+'.asc')
        blobName = path.join(blobDir,URID+'.'+masksCollShortNames)

        LOG.debug("ascFileName : %s" % (ascFileName))
        LOG.debug("blobName : %s" % (blobName))

        # Create a new ancillary blob, and copy the data to it.
        newQSTLWMblobObj = adl_blob.create(xmlName, blobName, endian=endian, overwrite=True)
        newQSTLWMblobArrObj = newQSTLWMblobObj.as_arrays()

        blobData = getattr(newQSTLWMblobArrObj,GridIP_shortNameToBlobName[masksCollShortNames])
        blobData[:,:] = QSTLWM[:,:]

        # Make a new QSTLWM asc file from the template, and substitute for the various tags

        ascTemplateFileName = path.join(ADL_ASC_TEMPLATES,"VIIRS-ANC_Template.asc")

        LOG.debug("Creating new asc file\n%s\nfrom template\n%s" % (ascFileName,ascTemplateFileName))
        
        ANC_fileList = ['IGBP.EcoMap.v1.0.2004.129.v004.h5','dem30ARC_Global_LandWater.h5']
        for idx in range(len(ANC_fileList)) :
            ANC_fileList[idx] = path.basename(ANC_fileList[idx])
        ANC_fileList.sort()
        ancGroupRecipe = '    ("N_Anc_Filename" STRING EQ "%s")'
        ancFileStr = "%s" % ("\n").join([ancGroupRecipe % (str(files)) for files in ANC_fileList])

        # Parse the geolocation asc file to get struct information which will be 
        # written to the ancillary asc files

        geoAscFileName = path.join(inDir,geoDict['URID']+".asc")
        LOG.debug("\nOpening %s..." % (geoAscFileName))

        geoAscFile = open(geoAscFileName,'rt')

        #RangeDateTimeStr =  _getAscLine(geoAscFile,"RangeDateTime")
        RangeDateTimeStr =  _getAscLine(geoAscFile,"ObservedDateTime")
        RangeDateTimeStr =  string.replace(RangeDateTimeStr,"ObservedDateTime","RangeDateTime")

        GRingLatitudeStr =  _getAscStructs(geoAscFile,"GRingLatitude",12)
        GRingLongitudeStr =  _getAscStructs(geoAscFile,"GRingLongitude",12)

        geoAscFile.close()

        ascTemplateFile = open(ascTemplateFileName,"rt") # Open template file for reading
        ascFile = open(ascFileName,"wt") # create a new text file

        LOG.debug("Template file %s is %r with mode %s" %(ascTemplateFileName,'not open' if ascTemplateFile.closed else 'open',ascTemplateFile.mode))

        LOG.debug("New file %s is %r with mode %s" %(ascFileName,'not open' if ascFile.closed else 'open',ascFile.mode))

        for line in ascTemplateFile.readlines():
           line = line.replace("CSPP_URID",URID)
           line = line.replace("CSPP_CREATIONDATETIME_NOUSEC",creationDate_nousecStr)
           line = line.replace("CSPP_ANC_BLOB_FULLPATH",blobName)
           line = line.replace("CSPP_ANC_COLLECTION_SHORT_NAME",masksCollShortNames)
           line = line.replace("CSPP_GRANULE_ID",N_Granule_ID)
           line = line.replace("CSPP_CREATIONDATETIME",creationDateStr)
           line = line.replace("  CSPP_RANGE_DATE_TIME",RangeDateTimeStr)
           line = line.replace("  CSPP_GRINGLATITUDE",GRingLatitudeStr)
           line = line.replace("  CSPP_GRINGLONGITUDE",GRingLongitudeStr)
           line = line.replace("    CSPP_ANC_SOURCE_FILES",ancFileStr)
           ascFile.write(line) 

        ascFile.close()
        ascTemplateFile.close()

    return QSTLWM


def _granulate_NISE(latitude,longitude,LSM,NISE_fileName):
    '''Granulates the NISE Snow and Ice files.'''

    #latitude = ma.masked_less(latitude,-800.)
    latMin,latMax = np.min(latitude),np.max(latitude)
    latRange = latMax-latMin

    #longitude = ma.masked_less(longitude,-800.)
    lonMin,lonMax = np.min(longitude),np.max(longitude)
    lonRange = lonMax-lonMin

    LOG.debug("\nmin,max,range of latitide: %f %f %f" % (latMin,latMax,latRange))
    LOG.debug("min,max,range of longitude: %f %f %f" % (lonMin,lonMax,lonRange))

    LOG.debug("\nGranulating NISE file %s" % (NISE_fileName))

    #################
    # Using HDF5
    #################

    #NISEobj = pytables.openFile(NISE_fileName,'r')
    #nHemi = NISEobj.getNode('/Northern Hemisphere/Data Fields/Extent')[:,:]
    #sHemi = NISEobj.getNode('/Southern Hemisphere/Data Fields/Extent')[:,:]
    #NISEobj.close()

    #################
    # Using HDF4
    #################

    LOG.debug("Opening the NISE file %s" % (NISE_fileName))
    fileObj = HDF4File(NISE_fileName)

    northDsetName = "Northern Hemisphere/Data Fields/Extent"
    southDsetName = "Southern Hemisphere/Data Fields/Extent"
    #northDsetName = "Northern Hemisphere/Data Fields/Age"
    #southDsetName = "Southern Hemisphere/Data Fields/Age"

    LOG.debug("Retrieving NISE HDF4 path '%s'" % (northDsetName))
    nHemi = fileObj.getDataset(northDsetName)

    LOG.debug("Retrieving NISE HDF4 path '%s'" % (southDsetName))
    sHemi = fileObj.getDataset(southDsetName)

    #################

    xsize, ysize = nHemi.shape    

    nHemi = nHemi.flatten()
    sHemi = sHemi.flatten()
    
    pi = np.pi

    rg = 6371.228 / 25.067525

    r0 = (ysize-1) / 2.0
    s0 = (xsize-1) / 2.0

    gridrows, gridcols = latitude.shape

    posLatMask = ma.masked_less(latitude,0.).mask
    negLatMask = ma.masked_greater_equal(latitude,0.).mask

    # For positive latitudes...
    
    phi = np.radians(ma.array(latitude,mask=posLatMask,fill_value=0))
    lam = np.radians(ma.array(longitude,mask=posLatMask,fill_value=0))

    rho =  2. * rg * np.sin((pi / 4.0) - (phi / 2.0))
    r = r0 + rho * np.sin(lam)
    s = s0 + rho * np.cos(lam)          

    i_pos = round_(r, 0).astype('int') - 1
    j_pos = round_(s, 0).astype('int') - 1


    # For negative latitudes...

    phi = np.radians(ma.array(latitude,mask=negLatMask,fill_value=0))
    lam = np.radians(ma.array(longitude,mask=negLatMask,fill_value=0))

    rho = 2. * rg * np.cos((pi / 4.0) - (phi / 2.0))
    r = r0 + rho * np.sin(lam)
    s = s0 - rho * np.cos(lam)

    i_neg = round_(r, 0).astype('int') - 1
    j_neg = round_(s, 0).astype('int') - 1

    ###
    # Combine the +ve and -ve latitudes
    ###

    posIdx = np.where(posLatMask==False)
    negIdx = np.where(negLatMask==False)

    # convert i_pos and j_pos to an index into the raveled nHemi..
    k_pos = (j_pos * xsize) + i_pos

    # convert i_neg and j_neg to an index into the raveled sHemi..
    k_neg = (j_neg * xsize) + i_neg

    nise_val = 999 * np.ones(latitude.shape,dtype='int')

    nise_val[posIdx] = nHemi[k_pos[posIdx]]
    nise_val[negIdx] = sHemi[k_neg[negIdx]]

    # From ADL/include/ProEdrViirsCMIPGbl.h ...

    CM_LAND = 1
    CM_COASTAL = 5
    CM_SEA_WATER = 3
    CM_IN_WATER = 2
    CM_SNOW = 1
    CM_NO_SNOW = 0

    # Permanent Ice
    permIceMask = ma.masked_equal(nise_val,101).mask

    # Dry Snow
    drySnowMask = ma.masked_equal(nise_val,103).mask

    # Wet Snow
    wetSnowMask = ma.masked_equal(nise_val,104).mask

    # Land and Coastal
    LSM_Land_Coastal = (ma.masked_equal(LSM,CM_LAND) * ma.masked_equal(LSM,CM_COASTAL)).mask

    # Sea Water and Inland Water
    LSM_SeaWater_InlandWater = (ma.masked_equal(LSM,CM_SEA_WATER) * ma.masked_equal(LSM,CM_IN_WATER)).mask

    # Transient sea ice > 25% concentration, AND (Land OR Coastal)
    transSeaIceCoastalMask = ma.masked_inside(nise_val,(25+1),(101-1)).mask
    
    # Transient sea ice > 25% concentration, AND Permanent Ice AND (Sea OR Inland Water)
    transSeaIceOceanMask = ma.masked_inside(nise_val,(25+1),(102-1)).mask

    snowMask = permIceMask + drySnowMask + wetSnowMask \
            + (transSeaIceCoastalMask * LSM_Land_Coastal) \
            + (transSeaIceOceanMask * LSM_SeaWater_InlandWater)

    snowIceMask  = np.ones(nise_val.shape,dtype='int') * CM_NO_SNOW
    snowIceMask = ma.array(snowIceMask,mask=snowMask,fill_value=CM_SNOW)
    snowIceMask = snowIceMask.filled()

    # Moderate resolution trim table arrays. These are 
    # bool arrays, and the trim pixels are set to True.

    trimObj = ViirsData.ViirsTrimTable()
    modTrimMask = trimObj.createModTrimArray(nscans=48,trimType=bool)

    # Fill the required pixel trim rows in the granulated NISE data with 
    # the ONBOARD_PT_FILL value for the correct data type

    fillValue = trimObj.sdrTypeFill['ONBOARD_PT_FILL']['uint8']
    snowIceMask = ma.array(snowIceMask,mask=modTrimMask,fill_value=fillValue)
    snowIceMask = snowIceMask.filled()

    return snowIceMask


def _granulate_NISE_list(inDir,geoDicts,DEM_granules,NISEfiles):
    '''Granulates the input NISE snow and ice files.'''

    global ancEndian 

    #CSPP_RT_ANC_HOME = os.getenv('CSPP_RT_ANC_HOME')
    
    ADL_ASC_TEMPLATES = path.join(CSPP_RT_ANC_HOME,'asc_templates')

    masksCollShortNames = 'VIIRS-GridIP-VIIRS-Snow-Ice-Cover-Mod-Gran'

    GridIP_shortNameToBlobName = {}
    GridIP_dataType = {}
    GridIP_shortNameToXmlName = {}
    GridIP_shortNameToBlobName['VIIRS-GridIP-VIIRS-Snow-Ice-Cover-Mod-Gran'] = 'snowIceCover'
    GridIP_dataType['VIIRS-GridIP-VIIRS-Snow-Ice-Cover-Mod-Gran'] = 'uint8'
    GridIP_shortNameToXmlName['VIIRS-GridIP-VIIRS-Snow-Ice-Cover-Mod-Gran'] = 'VIIRS_GRIDIP_VIIRS_SNOW_ICE_COVER_MOD_GRAN.xml'

    # Moderate resolution trim table arrays. These are 
    # bool arrays, and the trim pixels are set to True.

    trimObj = ViirsData.ViirsTrimTable()
    modTrimMask = trimObj.createModTrimArray(nscans=48,trimType=bool)

    # Get a bunch of information about the geolocation
    latitudeList,longitudeList,latMinList,latMaxList,lonMinList,lonMaxList,latCrnList,lonCrnList = _get_geo_Arrays(geoDicts)

    LOG.debug("\nGranules -->")
    LOG.debug("latMin : %r" % (latMinList))
    LOG.debug("latMax : %r" % (latMaxList))
    LOG.debug("lonMin : %r" % (lonMinList))
    LOG.debug("lonMax : %r" % (lonMaxList))
    LOG.debug("latCrnList : %r" % (latCrnList))
    LOG.debug("lonCrnList : %r\n" % (lonCrnList))

    LOG.debug("\nNISE files --> %r" % (NISEfiles))

    # Loop through the geolocation files and granulate...

    for latitude,longitude,LSM,geoDict in zip(latitudeList,longitudeList,DEM_granules,geoDicts):

        N_Granule_ID = geoDict['N_Granule_ID']

        LOG.info("Granulating %s ..." % (masksCollShortNames))

        ####################################################
        # Granulate the NISE data using this geolocation...
        ####################################################

        data = _granulate_NISE(latitude,longitude,LSM,NISEfiles[0])

        # Fill the required pixel trim rows in the granulated NCEP data with 
        # the ONBOARD_PT_FILL value for the correct data type

        #fillValue = trimObj.sdrTypeFill['ONBOARD_PT_FILL'][data.dtype.name]        
        #data = ma.array(data,mask=modTrimMask,fill_value=fillValue)
        #print "min of NISE granule = ",np.min(data)
        #print "max of NISE granule = ",np.max(data)
        #data = data.filled()

        # Create new NISE ancillary blob, and copy granulated data to it

        endian = ancEndian
        xmlName = path.join(ADL_HOME,'xml/VIIRS',GridIP_shortNameToXmlName[masksCollShortNames])

        # Create a new URID to be used in making the asc filenames

        URID_dict = _getURID()

        URID = URID_dict['URID']
        creationDate_nousecStr = URID_dict['creationDate_nousecStr']
        creationDateStr = URID_dict['creationDateStr']

        # Create a new directory in the input directory for the new ancillary
        # asc and blob files

        blobDir = inDir


        ascFileName = path.join(blobDir,URID+'.asc')
        blobName = path.join(blobDir,URID+'.'+masksCollShortNames)

        LOG.debug("ascFileName : %s" % (ascFileName))
        LOG.debug("blobName : %s" % (blobName))

        # Create a new ancillary blob, and copy the data to it.

        newNISEblobObj = adl_blob.create(xmlName, blobName, endian=endian, overwrite=True)
        newNISEblobArrObj = newNISEblobObj.as_arrays()

        blobData = getattr(newNISEblobArrObj,GridIP_shortNameToBlobName[masksCollShortNames])
        blobData[:,:] = data[:,:]

        # Make a new NCEP asc file from the template, and substitute for the various tags

        ascTemplateFileName = path.join(ADL_ASC_TEMPLATES,"VIIRS-ANC_Template.asc")

        LOG.debug("Creating new asc file\n%s\nfrom template\n%s" % (ascFileName,ascTemplateFileName))
        
        ANC_fileList = copy.copy(NISEfiles)
        for idx in range(len(ANC_fileList)) :
            ANC_fileList[idx] = path.basename(ANC_fileList[idx])
        ANC_fileList.sort()
        ancGroupRecipe = '    ("N_Anc_Filename" STRING EQ "%s")'
        ancFileStr = "%s" % ("\n").join([ancGroupRecipe % (str(files)) for files in ANC_fileList])

        # Parse the geolocation asc file to get struct information which will be 
        # written to the ancillary asc files

        geoAscFileName = path.join(inDir,geoDict['URID']+".asc")
        LOG.debug("\nOpening %s..." % (geoAscFileName))

        geoAscFile = open(geoAscFileName,'rt')

        #RangeDateTimeStr =  _getAscLine(geoAscFile,"RangeDateTime")
        RangeDateTimeStr =  _getAscLine(geoAscFile,"ObservedDateTime")
        RangeDateTimeStr =  string.replace(RangeDateTimeStr,"ObservedDateTime","RangeDateTime")

        GRingLatitudeStr =  _getAscStructs(geoAscFile,"GRingLatitude",12)
        GRingLongitudeStr =  _getAscStructs(geoAscFile,"GRingLongitude",12)

        geoAscFile.close()

        ascTemplateFile = open(ascTemplateFileName,"rt") # Open template file for reading
        ascFile = open(ascFileName,"wt") # create a new text file

        LOG.debug("Template file %s is %r with mode %s" %(ascTemplateFileName,'not open' if ascTemplateFile.closed else 'open',ascTemplateFile.mode))

        LOG.debug("New file %s is %r with mode %s" %(ascFileName,'not open' if ascFile.closed else 'open',ascFile.mode))

        for line in ascTemplateFile.readlines():
           line = line.replace("CSPP_URID",URID)
           line = line.replace("CSPP_CREATIONDATETIME_NOUSEC",creationDate_nousecStr)
           line = line.replace("CSPP_ANC_BLOB_FULLPATH",blobName)
           line = line.replace("CSPP_ANC_COLLECTION_SHORT_NAME",masksCollShortNames)
           line = line.replace("CSPP_GRANULE_ID",N_Granule_ID)
           line = line.replace("CSPP_CREATIONDATETIME",creationDateStr)
           line = line.replace("  CSPP_RANGE_DATE_TIME",RangeDateTimeStr)
           line = line.replace("  CSPP_GRINGLATITUDE",GRingLatitudeStr)
           line = line.replace("  CSPP_GRINGLONGITUDE",GRingLongitudeStr)
           line = line.replace("    CSPP_ANC_SOURCE_FILES",ancFileStr)
           ascFile.write(line) 

        ascFile.close()
        ascTemplateFile.close()


def _retrieve_grib_files(geoDicts):
    ''' Download the GRIB files which cover the dates of the geolocation files.'''

    CSPP_RT_HOME = os.getenv('CSPP_RT_HOME')
    ANC_SCRIPTS_PATH = path.join(CSPP_RT_HOME,'viirs')
    CSPP_RT_ANC_CACHE_DIR = os.getenv('CSPP_RT_ANC_CACHE_DIR')

    # FIXME : Fix rounding up of the seconds if the decisecond>=9.5 
    gribFiles = []
    for geoDict in geoDicts:
        #timeObj = geoDict['StartTime']
        timeObj = geoDict['ObservedStartTime']
        dateStamp = timeObj.strftime("%Y%m%d")
        seconds = repr(int(round(timeObj.second + float(timeObj.microsecond)/1000000.)))
        deciSeconds = int(round(float(timeObj.microsecond)/100000.))
        deciSeconds = repr(0 if deciSeconds > 9 else deciSeconds)
        startTimeStamp = "%s%s" % (timeObj.strftime("%H%M%S"),deciSeconds)

        #timeObj = geoDict['EndTime']
        timeObj = geoDict['ObservedEndTime']
        seconds = repr(int(round(timeObj.second + float(timeObj.microsecond)/1000000.)))
        deciSeconds = int(round(float(timeObj.microsecond)/100000.))
        deciSeconds = repr(0 if deciSeconds > 9 else deciSeconds)
        endTimeStamp = "%s%s" % (timeObj.strftime("%H%M%S"),deciSeconds)

        timeObj = geoDict['UnpackTime']
        unpackTimeStamp = timeObj.strftime("%Y%m%d%H%M%S%f")

        granuleName = "GMODO_npp_d%s_t%s_e%s_b00014_c%s.h5" % (dateStamp,startTimeStamp,endTimeStamp,unpackTimeStamp)

        try :
            LOG.info('Retrieving NCEP files for %s ...' % (granuleName))
            cmdStr = '%s/cspp_retrieve_gdas_gfs.csh %s' % (ANC_SCRIPTS_PATH,granuleName)
            LOG.info('\t%s' % (cmdStr))
            args = shlex.split(cmdStr)
            LOG.debug('\t%s' % (repr(args)))

            procRetVal = 0
            procObj = subprocess.Popen(args,bufsize=0, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            procObj.wait()
            procRetVal = procObj.returncode

            procOutput = procObj.stdout.readlines()
            #procOutput = procObj.stdout.read()

            # FIXME : How to get this output to have linebreaks when using readlines()
            #LOG.debug(procOutput)
            
            for lines in procOutput:
                if "GDAS/GFS file" in lines :
                    lines = string.replace(lines,'GDAS/GFS file 1: ','')
                    lines = string.replace(lines,'GDAS/GFS file 2: ','')
                    lines = string.replace(lines,'\n','')
                    gribFiles.append(lines)

            # TODO : On error, jump to a cleanup routine
            if not (procRetVal == 0) :
                LOG.error('Retrieval of ancillary files failed for %s.' % (granuleName))
                #sys.exit(procRetVal)

        except Exception, err:
            LOG.warn( "%s" % (str(err)))

    # Get a unique list of grib files that were fetched
    gribFiles.sort()
    gribFiles = dict(map(lambda i: (i,1),gribFiles)).keys()
    gribFiles.sort()

    for gribFile in gribFiles :
        LOG.info('Retrieved gribfile: %r' % (gribFile))

    return gribFiles


def _create_NCEP_gridBlobs(gribFiles):
    '''Converts NCEP GRIB files into NCEP blobs'''

    blobFiles = []

    CSPP_RT_HOME = os.getenv('CSPP_RT_HOME')
    ANC_SCRIPTS_PATH = path.join(CSPP_RT_HOME,'viirs')
    CSPP_RT_ANC_CACHE_DIR = os.getenv('CSPP_RT_ANC_CACHE_DIR')
    csppPython = os.getenv('PY')
    

    for files in gribFiles :
        gribPath = path.dirname(files)
        gribFile = path.basename(files)
        gribBlob = "%s_blob.le" % (gribFile)
        gribBlob = path.join(gribPath,gribBlob)
        LOG.debug("Creating NCEP GRIB blob %s" % (gribBlob))

        if not os.path.exists(gribBlob):
            try :
                LOG.info('Transcoding %s to %s ...' % (files,gribBlob))
                NCEPscript = '%s/NCEPtoBlob.py' % (ANC_SCRIPTS_PATH)
                NCEPxml = path.join(ADL_HOME,'xml/ANC/NCEP_ANC_Int.xml')
                cmdStr = '%s %s -x %s -g %s -e little -o %s' % (csppPython,NCEPscript,NCEPxml,files,gribBlob)
                LOG.info('\t%s' % (cmdStr))
                args = shlex.split(cmdStr)
                LOG.debug('\t%s' % (repr(args)))

                procRetVal = 0
                procObj = subprocess.Popen(args,bufsize=0, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                procObj.wait()
                procRetVal = procObj.returncode

                procOutput = procObj.stdout.read()

                LOG.debug("%s" % (procOutput))
                
                if not (procRetVal == 0) :
                    LOG.error('Transcoding of ancillary files failed for %s.' % (files))
                    sys.exit(procRetVal)
                else :
                    blobFiles.append(gribBlob)

            except Exception, err:
                LOG.warn( "%s" % (str(err)))
        else :
            LOG.info('Gridded NCEP blob files %s exists, skipping.' % (gribBlob))
            blobFiles.append(gribBlob)

    return blobFiles


def _create_NCEP_gridBlobs_alt(gribFiles):
    '''Converts NCEP GRIB files into NCEP blobs'''

    from copy import deepcopy

    blobFiles = []

    CSPP_RT_HOME = os.getenv('CSPP_RT_HOME')
    ANC_SCRIPTS_PATH = path.join(CSPP_RT_HOME,'viirs')
    CSPP_RT_ANC_CACHE_DIR = os.getenv('CSPP_RT_ANC_CACHE_DIR')
    csppPython = os.getenv('PY')
    

    for files in gribFiles :
        gribPath = path.dirname(files)
        gribFile = path.basename(files)
        gribBlob = "%s_blob.le" % (gribFile)
        gribBlob = path.join(gribPath,gribBlob)
        LOG.debug("Creating NCEP GRIB blob %s" % (gribBlob))

        if not os.path.exists(gribBlob):
            try :
                LOG.info('Transcoding %s to %s ...' % (files,gribBlob))
                NCEPxml = path.join(ADL_HOME,'xml/ANC/NCEP_ANC_Int.xml')

                # Create the grib object and populate with the grib file data
                NCEPobj = NCEPclass(gribFile=files)

                # Convert surface pressure from Pa to mb or hPa ...
                NCEPobj.NCEPmessages['surfacePressure'].data /= 100.

                # Convert total column ozone from DU or kg m**-2 to Atm.cm ...
                NCEPobj.NCEPmessages['totalColumnOzone'].data /= 1000.

                # Convert total precipitable water kg m^{-2} to cm ...
                NCEPobj.NCEPmessages['totalPrecipitableWater'].data /= 10.

                # Convert specific humidity in kg.kg^{-1} to water vapor mixing ratio in g.kg^{-1}
                moistureObj = NCEPobj.NCEPmessages['waterVaporMixingRatioLayers'].messageLevelData
                temperatureObj = NCEPobj.NCEPmessages['temperatureLayers'].messageLevelData

                # Compute the 100mb mixing ratio in g/kg
                if moistureObj['100'].name == 'Specific humidity':
                    specHumidity_100mb = moistureObj['100'].data
                    mixingRatio_100mb = 1000. * specHumidity_100mb/(1. - specHumidity_100mb)
                elif  moistureObj['100'].name == 'Relative humidity':
                    relativeHumidity_100mb = moistureObj['100'].data
                    temperature_100mb = temperatureObj['100'].data
                    mixingRatio_100mb = rh_to_mr_vec(relativeHumidity_100mb,100.,temperature_100mb)
                else :
                    pass

                for level,levelIdx in NCEPobj.NCEP_LAYER_LEVELS.items() : 
                    levelStrIdx = level[:-2]
                    pressure = NCEPobj.NCEP_LAYER_VALUES[levelIdx]

                    if pressure < 100. :
                        # Copy the 100mb message object to this pressure, and assign the correct
                        # mixing ratio
                        moistureObj[levelStrIdx] = deepcopy(moistureObj['100'])
                        moistureObj[levelStrIdx].level = levelStrIdx

                        # Compute the mixing ratio in g/kg
                        mixingRatio = np.maximum(mixingRatio_100mb,0.003) * ((pressure/100.)**3.)
                        mixingRatio = np.maximum(mixingRatio_100mb,0.003)
                    else :
                        # Compute the mixing ratio in g/kg
                        if moistureObj[levelStrIdx].name == 'Specific humidity':
                            specHumidity = moistureObj[levelStrIdx].data 
                            mixingRatio = 1000. * specHumidity/(1. - specHumidity)
                        elif  moistureObj[levelStrIdx].name == 'Relative humidity':
                            relativeHumidity = moistureObj[levelStrIdx].data
                            temperature = temperatureObj[levelStrIdx].data
                            mixingRatio = rh_to_mr_vec(relativeHumidity,pressure,temperature)
                        else :
                            pass

                    moistureObj[levelStrIdx].data = mixingRatio

                # Write the contents of the NCEPobj object to an ADL blob file
                endian = adl_blob.LITTLE_ENDIAN
                procRetVal = NCEPclass.NCEPgribToBlob_interpNew(NCEPobj,NCEPxml,gribBlob,endian=endian)

                #blobFiles.append(gribBlob)
                if not (procRetVal == 0) :
                    LOG.error('Transcoding of ancillary files failed for %s.' % (files))
                    sys.exit(procRetVal)
                else :
                    LOG.info('Finished creating NCEP GRIB blob %s' % (gribBlob))
                    blobFiles.append(gribBlob)

            except Exception, err:
                LOG.warn( "%s" % (str(err)))
        else :
            LOG.info('Gridded NCEP blob files %s exists, skipping.' % (gribBlob))
            blobFiles.append(gribBlob)

    LOG.info('Returning NCEP GRIB blob file names %r' % (blobFiles))
    return blobFiles


def _retrieve_NAAPS_grib_files(geoDicts):
    ''' Download the NAAPS GRIB files which cover the dates of the geolocation files.'''

    CSPP_RT_HOME = os.getenv('CSPP_RT_HOME')
    ANC_SCRIPTS_PATH = path.join(CSPP_RT_HOME,'viirs')
    CSPP_RT_ANC_CACHE_DIR = os.getenv('CSPP_RT_ANC_CACHE_DIR')

    # FIXME : Fix rounding up of the seconds if the decisecond>=9.5 
    gribFiles = []
    for geoDict in geoDicts:
        #timeObj = geoDict['StartTime']
        timeObj = geoDict['ObservedStartTime']
        dateStamp = timeObj.strftime("%Y%m%d")
        seconds = repr(int(round(timeObj.second + float(timeObj.microsecond)/1000000.)))
        deciSeconds = int(round(float(timeObj.microsecond)/100000.))
        deciSeconds = repr(0 if deciSeconds > 9 else deciSeconds)
        startTimeStamp = "%s%s" % (timeObj.strftime("%H%M%S"),deciSeconds)

        #timeObj = geoDict['EndTime']
        timeObj = geoDict['ObservedEndTime']
        seconds = repr(int(round(timeObj.second + float(timeObj.microsecond)/1000000.)))
        deciSeconds = int(round(float(timeObj.microsecond)/100000.))
        deciSeconds = repr(0 if deciSeconds > 9 else deciSeconds)
        endTimeStamp = "%s%s" % (timeObj.strftime("%H%M%S"),deciSeconds)

        timeObj = geoDict['UnpackTime']
        unpackTimeStamp = timeObj.strftime("%Y%m%d%H%M%S%f")

        granuleName = "GMODO_npp_d%s_t%s_e%s_b00014_c%s.h5" % (dateStamp,startTimeStamp,endTimeStamp,unpackTimeStamp)

        try :
            LOG.info('Retrieving NAAPS files for %s ...' % (granuleName))
            cmdStr = '%s/cspp_retrieve_gdas_gfs.csh %s' % (ANC_SCRIPTS_PATH,granuleName)
            LOG.info('\t%s' % (cmdStr))
            args = shlex.split(cmdStr)
            LOG.debug('\t%s' % (repr(args)))

            procRetVal = 0
            procObj = subprocess.Popen(args,bufsize=0, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            procObj.wait()
            procRetVal = procObj.returncode

            procOutput = procObj.stdout.readlines()
            #procOutput = procObj.stdout.read()

            # FIXME : How to get this output to have linebreaks when using readlines()
            #LOG.debug(procOutput)
            
            for lines in procOutput:
                if "GDAS/GFS file" in lines :
                    lines = string.replace(lines,'GDAS/GFS file 1: ','')
                    lines = string.replace(lines,'GDAS/GFS file 2: ','')
                    lines = string.replace(lines,'\n','')
                    gribFiles.append(lines)

            # TODO : On error, jump to a cleanup routine
            if not (procRetVal == 0) :
                LOG.error('Retrieval of ancillary files failed for %s.' % (granuleName))
                #sys.exit(procRetVal)

        except Exception, err:
            LOG.warn( "%s" % (str(err)))

    # Get a unique list of grib files that were fetched
    gribFiles.sort()
    gribFiles = dict(map(lambda i: (i,1),gribFiles)).keys()
    gribFiles.sort()

    for gribFile in gribFiles :
        LOG.info('Retrieved gribfile: %r' % (gribFile))

    return gribFiles


def _create_NAAPS_gridBlobs(gribFiles):
    '''Converts NAAPS GRIB files into NAAPS blobs'''

    from copy import deepcopy

    blobFiles = []

    CSPP_RT_HOME = os.getenv('CSPP_RT_HOME')
    ANC_SCRIPTS_PATH = path.join(CSPP_RT_HOME,'viirs')
    CSPP_RT_ANC_CACHE_DIR = os.getenv('CSPP_RT_ANC_CACHE_DIR')
    #CSPP_RT_ANC_HOME = os.getenv('CSPP_RT_ANC_HOME')
    csppPython = os.getenv('PY')
    

    '''
    for files in gribFiles :
        gribPath = path.dirname(files)
        gribFile = path.basename(files)
        gribBlob = "%s_blob.le" % (gribFile)
        gribBlob = path.join(gribPath,gribBlob)
        LOG.debug("Creating NAAPS GRIB blob %s" % (gribBlob))

        if not os.path.exists(gribBlob):
            try :
                LOG.info('Transcoding %s to %s ...' % (files,gribBlob))
                NAAPSxml = path.join(ADL_HOME,'xml/ANC/NAAPS_ANC_Int.xml')

                # Create the grib object and populate with the grib file data
                NAAPSobj = NAAPSclass(gribFile=files)

                # Write the contents of the NAAPSobj object to an ADL blob file
                endian = adl_blob.LITTLE_ENDIAN
                procRetVal = NAAPSclass.NAAPSgribToBlob_interpNew(NAAPSobj,NAAPSxml,gribBlob,endian=endian)

                #blobFiles.append(gribBlob)
                if not (procRetVal == 0) :
                    LOG.error('Transcoding of ancillary files failed for %s.' % (files))
                    sys.exit(procRetVal)
                else :
                    LOG.info('Finished creating NAAPS GRIB blob %s' % (gribBlob))
                    blobFiles.append(gribBlob)

            except Exception, err:
                LOG.warn( "%s" % (str(err)))
        else :
            LOG.info('Gridded NAAPS blob files %s exists, skipping.' % (gribBlob))
            blobFiles.append(gribBlob)
    '''

    blobFiles = []

    LOG.info('Returning NAAPS GRIB blob file names %r' % (blobFiles))
    return blobFiles


def _grid2Gran(dataLat, dataLon, gridData, gridLat, gridLon):
    '''Granulates a gridded dataset using an input geolocation'''

    nData = np.int64(dataLat.size)
    gridRows = np.int32(gridLat.shape[0])
    gridCols = np.int32(gridLat.shape[1])

    data = np.ones(np.shape(dataLat),dtype=np.float64)* -999.9
    dataIdx  = np.ones(np.shape(dataLat),dtype=np.int64) * -99999

    CSPP_RT_HOME = os.getenv('CSPP_RT_HOME')
    ANC_SCRIPTS_PATH = path.join(CSPP_RT_HOME,'viirs')

    libFile = path.join(ANC_SCRIPTS_PATH,'libgriddingAndGranulation.so.1.0.1')
    LOG.debug("Gridding and granulation library file: %s" % (libFile))
    lib = ctypes.cdll.LoadLibrary(libFile)
    grid2gran = lib.grid2gran
    grid2gran.restype = None
    grid2gran.argtypes = [
            ndpointer(ctypes.c_double,ndim=1,shape=(nData),flags='C_CONTIGUOUS'),
            ndpointer(ctypes.c_double,ndim=1,shape=(nData),flags='C_CONTIGUOUS'),
            ndpointer(ctypes.c_double,ndim=1,shape=(nData),flags='C_CONTIGUOUS'),
            ctypes.c_int64,
            ndpointer(ctypes.c_double,ndim=2,shape=(gridRows,gridCols),flags='C_CONTIGUOUS'),
            ndpointer(ctypes.c_double,ndim=2,shape=(gridRows,gridCols),flags='C_CONTIGUOUS'),
            ndpointer(ctypes.c_double,ndim=2,shape=(gridRows,gridCols),flags='C_CONTIGUOUS'),
            ndpointer(ctypes.c_int64,ndim=1,shape=(nData),flags='C_CONTIGUOUS'),
            ctypes.c_int32,
            ctypes.c_int32
            ]

    t1 = time()

    '''
    int snapGrid_ctypes(double *lat, 
                    double *lon, 
                    double *data, 
                    long nData, 
                    double *gridLat,
                    double *gridLon,
                    double *gridData,
                    long *gridDataIdx,
                    int nGridRows,
                    int nGridCols
                    )
    '''

    LOG.info("Calling C routine grid2gran()...")

    retVal = grid2gran(dataLat,
                       dataLon,
                       data,
                       nData,
                       gridLat,
                       gridLon,
                       gridData,
                       dataIdx,
                       gridRows,
                       gridCols)

    LOG.info("Returning from C routine grid2gran()")

    t2 = time()
    elapsedTime = t2-t1
    LOG.info("Granulation took %f seconds for %d points" % (elapsedTime,nData))

    return data,dataIdx

def _getAscLine(fileObj,searchString):
    dataStr = ''
    try :
        while True :
            line = fileObj.readline()

            if searchString in line : 
                dataStr = "%s" % (string.replace(line,'\n',''));
                break

        fileObj.seek(0)

    except Exception, err:
        geoAscFile.close()

    return dataStr

def _getAscStructs(fileObj,searchString,linesOfContext):

    dataList = []
    data_count = 0
    dataFound = False

    try :
        while True :
            line = fileObj.readline()

            if searchString in line : 
                dataFound = True

            if dataFound :
                dataStr = "%s" % (string.replace(line,'\n',''));
                dataList.append(dataStr)
                data_count += 1
            else :
                pass

            if (data_count == linesOfContext) :
                break

        fileObj.seek(0)

    except Exception, err:
        fileObj.close()
        return -1

    dataStr=''
    dataStr = "%s" % ("\n").join(['%s' % (str(lines)) for lines in dataList])

    return dataStr

def _granulate_NCEP_gridBlobs(inDir,geoDicts, gridBlobFiles):
    '''Granulates the input gridded blob files into the required NCEP granulated datasets.'''

    global ancEndian 

    CSPP_RT_HOME = os.getenv('CSPP_RT_HOME')
    ANC_SCRIPTS_PATH = path.join(CSPP_RT_HOME,'viirs')
    CSPP_RT_ANC_CACHE_DIR = os.getenv('CSPP_RT_ANC_CACHE_DIR')
    #CSPP_RT_ANC_HOME = os.getenv('CSPP_RT_ANC_HOME')
    csppPython = os.getenv('PY')
    

    ADL_ASC_TEMPLATES = path.join(CSPP_RT_ANC_HOME,'asc_templates')
    
    # Collection shortnames of the required NCEP ancillary datasets
    # FIXME : Poll ADL/cfg/ProEdrViirsCM_CFG.xml for this information

    masksCollShortNamesDict = {'VCM' : [
                           'VIIRS-ANC-Preci-Wtr-Mod-Gran',
                           'VIIRS-ANC-Temp-Surf2M-Mod-Gran',
                           'VIIRS-ANC-Wind-Speed-Mod-Gran',
                           'VIIRS-ANC-Surf-Ht-Mod-Gran'
                          ],
                          'AOT' : [
                            'VIIRS-ANC-Preci-Wtr-Mod-Gran',
                            'VIIRS-ANC-Temp-Surf2M-Mod-Gran',
                            'VIIRS-ANC-Wind-Speed-Mod-Gran',
                            'VIIRS-ANC-Wind-Direction-Mod-Gran',
                            'VIIRS-ANC-Press-Surf-Mod-Gran',
                            'VIIRS-ANC-Tot-Col-Mod-Gran'
                              ]
                          }

    masksCollShortNames = [
                           'VIIRS-ANC-Preci-Wtr-Mod-Gran',
                           'VIIRS-ANC-Temp-Surf2M-Mod-Gran',
                           'VIIRS-ANC-Wind-Speed-Mod-Gran',
                           'VIIRS-ANC-Surf-Ht-Mod-Gran'
                          ]

    # Dictionary relating the required NCEP collection short names and the 
    # NCEP gridded blob dataset names

    NCEP_shortNameToBlobName = {}
    NCEP_shortNameToBlobName['VIIRS-ANC-Preci-Wtr-Mod-Gran'] = 'totalPrecipitableWater'
    NCEP_shortNameToBlobName['VIIRS-ANC-Temp-Surf2M-Mod-Gran'] = 'surfaceTemperature'
    NCEP_shortNameToBlobName['VIIRS-ANC-Wind-Speed-Mod-Gran'] = ['uComponentOfWind', 'vComponentOfWind']
    NCEP_shortNameToBlobName['VIIRS-ANC-Surf-Ht-Mod-Gran'] = 'surfaceGeopotentialHeight'

    NCEP_shortNameToBlobName['VIIRS-ANC-Press-Surf-Mod-Gran'] = 'surfacePressure'
    NCEP_shortNameToBlobName['VIIRS-ANC-Tot-Col-Mod-Gran'] = 'totalColumnOzone'

    NCEP_shortNameToXmlName = {}
    NCEP_shortNameToXmlName['VIIRS-ANC-Preci-Wtr-Mod-Gran'] = 'VIIRS_ANC_PRECI_WTR_MOD_GRAN.xml'
    NCEP_shortNameToXmlName['VIIRS-ANC-Temp-Surf2M-Mod-Gran'] = 'VIIRS_ANC_TEMP_SURF2M_MOD_GRAN.xml'
    NCEP_shortNameToXmlName['VIIRS-ANC-Wind-Speed-Mod-Gran'] = 'VIIRS_ANC_WIND_SPEED_MOD_GRAN.xml'
    NCEP_shortNameToXmlName['VIIRS-ANC-Surf-Ht-Mod-Gran'] = 'VIIRS_ANC_SURF_HT_MOD_GRAN.xml'

    NCEP_shortNameToXmlName['VIIRS-ANC-Press-Surf-Mod-Gran'] = 'VIIRS_ANC_PRESS_SURF_MOD_GRAN.xml'
    NCEP_shortNameToXmlName['VIIRS-ANC-Tot-Col-Mod-Gran'] = 'VIIRS_ANC_TOT_COL_MOD_GRAN.xml'

    # Moderate resolution trim table arrays. These are 
    # bool arrays, and the trim pixels are set to True.

    trimObj = ViirsData.ViirsTrimTable()
    modTrimMask = trimObj.createModTrimArray(nscans=48,trimType=bool)

    # Open the NCEP gridded blob file

    ncepXmlFile = path.join(ADL_HOME,'xml/ANC/NCEP_ANC_Int.xml')
    # FIXME : Should be using two NCEP blob files, and averaging
    gridBlobFile = gridBlobFiles[0]

    if os.path.exists(ncepXmlFile):
        LOG.info("We are using for %s: %s,%s" %('NCEP-ANC-Int',ncepXmlFile,gridBlobFile))
    
    endian = ancEndian

    ncepBlobObj = adl_blob.map(ncepXmlFile,gridBlobFile, endian=endian)
    ncepBlobArrObj = ncepBlobObj.as_arrays()

    LOG.debug("%s...\n%r" % (gridBlobFile,ncepBlobArrObj._fields))

    # Save the global grids of the required datasets into a dictionary...

    NCEP_globalGridData = {}

    # Precipitable water
    blobDsetName = NCEP_shortNameToBlobName['VIIRS-ANC-Preci-Wtr-Mod-Gran']
    NCEP_globalGridData['VIIRS-ANC-Preci-Wtr-Mod-Gran'] = \
            getattr(ncepBlobArrObj,blobDsetName).astype('float')
    LOG.info("Shape of dataset %s is %s" % (blobDsetName,np.shape(NCEP_globalGridData['VIIRS-ANC-Preci-Wtr-Mod-Gran'])))

    # 2m Surface temperature
    blobDsetName = NCEP_shortNameToBlobName['VIIRS-ANC-Temp-Surf2M-Mod-Gran']
    NCEP_globalGridData['VIIRS-ANC-Temp-Surf2M-Mod-Gran'] = \
            getattr(ncepBlobArrObj,blobDsetName).astype('float')
    LOG.info("Shape of dataset %s is %s" % (blobDsetName,np.shape(NCEP_globalGridData['VIIRS-ANC-Temp-Surf2M-Mod-Gran'])))

    # Omnidirectional wind speed
    blobDsetName = NCEP_shortNameToBlobName['VIIRS-ANC-Wind-Speed-Mod-Gran'][0]
    uWind = getattr(ncepBlobArrObj,blobDsetName).astype('float')
    blobDsetName = NCEP_shortNameToBlobName['VIIRS-ANC-Wind-Speed-Mod-Gran'][1]
    vWind = getattr(ncepBlobArrObj,blobDsetName).astype('float')
    windSpeed = np.sqrt(uWind*uWind + vWind*vWind)
    NCEP_globalGridData['VIIRS-ANC-Wind-Speed-Mod-Gran'] = windSpeed
    LOG.info("Shape of dataset %s is %s" % ('windSpeed',np.shape(NCEP_globalGridData['VIIRS-ANC-Wind-Speed-Mod-Gran'])))

    # Terrain height
    # FIXME : VIIRS-ANC-Surf-Ht-Mod-Gran is supposed to come from the Terrain-Eco-ANC-Tile 
    # FIXME :   gridded data. The granulated terrain height from this source is provided in
    # FIXME :   the VIIRS-MOD-RGEO-TC from the VIIRS SDR controller. If we get the non 
    # FIXME :   terrain correction geolocation (VIIRS-MOD-RGEO) we can use the surface
    # FIXME :   geopotential height, which is approximately the same.
    blobDsetName = NCEP_shortNameToBlobName['VIIRS-ANC-Surf-Ht-Mod-Gran']
    NCEP_globalGridData['VIIRS-ANC-Surf-Ht-Mod-Gran'] = \
            getattr(ncepBlobArrObj,blobDsetName).astype('float')
    LOG.info("Shape of dataset %s is %s" % (blobDsetName,np.shape(NCEP_globalGridData['VIIRS-ANC-Surf-Ht-Mod-Gran'])))

    # Contruct a default 0.5 degree grid...

    degInc = 0.5
    grids = np.mgrid[-90.:90.+degInc:degInc,-180.:180.:degInc]
    gridLat,gridLon = grids[0],grids[1]

    # Loop through the geolocation files and granulate...

    for dicts in geoDicts :

        URID = dicts['URID']
        geo_Collection_ShortName = dicts['N_Collection_Short_Name']
        N_Granule_ID = dicts['N_Granule_ID']
        ObservedStartTimeObj = dicts['ObservedStartTime']
        geoFiles = glob('%s/%s*' % (inDir,URID))
        geoFiles.sort()

        LOG.debug("\n###########################")
        LOG.debug("  Geolocation Information  ")
        LOG.debug("###########################")
        LOG.debug("N_Granule_ID : %r" % (N_Granule_ID))
        LOG.debug("ObservedStartTime : %s" % (ObservedStartTimeObj.__str__()))
        LOG.debug("N_Collection_Short_Name : %s" %(geo_Collection_ShortName))
        LOG.debug("URID : %r" % (URID))
        LOG.debug("geoFiles : %r" % (geoFiles))
        LOG.debug("###########################\n")

        # Do we have terrain corrected geolocation?

        terrainCorrectedGeo = True if 'GEO-TC' in geo_Collection_ShortName else False

        # Do we have long or short style geolocation field names?

        if (geo_Collection_ShortName=='VIIRS-MOD-GEO-TC' or geo_Collection_ShortName=='VIIRS-MOD-RGEO') :
            longFormGeoNames = True
            LOG.debug("We have long form geolocation names")
        elif (geo_Collection_ShortName=='VIIRS-MOD-GEO' or geo_Collection_ShortName=='VIIRS-MOD-RGEO-TC') :
            LOG.debug("We have short form geolocation names")
            longFormGeoNames = False
        else :
            LOG.error("Invalid geolocation shortname: %s",geo_Collection_ShortName)
            return -1

        # Get the geolocation xml file

        geoXmlFile = "%s.xml" % (string.replace(geo_Collection_ShortName,'-','_'))
        geoXmlFile = path.join(ADL_HOME,'xml/VIIRS',geoXmlFile)
        if os.path.exists(geoXmlFile):
            LOG.debug("We are using for %s: %s,%s" %(geo_Collection_ShortName,geoXmlFile,geoFiles[0]))

        # Open the geolocation blob and get the latitude and longitude

        endian=sdrEndian

        geoBlobObj = adl_blob.map(geoXmlFile,geoFiles[0], endian=endian)
        geoBlobArrObj = geoBlobObj.as_arrays()

        # If we have the terrain corrected geolocation, get the terrain height

        if terrainCorrectedGeo :
            terrainHeight = geoBlobArrObj.height[:,:]

        # Get scan_mode to find any bad scans

        scanMode = geoBlobArrObj.scan_mode[:]
        badScanIdx = np.where(scanMode==254)[0]
        LOG.debug("Bad Scans: %r" % (badScanIdx))

        # Detemine the min, max and range of the latitude and longitude, 
        # taking care to exclude any fill values.

        if longFormGeoNames :
            latitude = getattr(geoBlobArrObj,'latitude').astype('float')
            longitude = getattr(geoBlobArrObj,'longitude').astype('float')
        else :
            latitude = getattr(geoBlobArrObj,'lat').astype('float')
            longitude = getattr(geoBlobArrObj,'lon').astype('float')
        
        if terrainCorrectedGeo :
            terrain = getattr(geoBlobArrObj,'height').astype('float')

        latitude = ma.masked_less(latitude,-800.)
        latMin,latMax = np.min(latitude),np.max(latitude)
        latRange = latMax-latMin

        longitude = ma.masked_less(longitude,-800.)
        lonMin,lonMax = np.min(longitude),np.max(longitude)
        lonRange = lonMax-lonMin

        LOG.debug("min,max,range of latitide: %f %f %f" % (latMin,latMax,latRange))
        LOG.debug("min,max,range of longitude: %f %f %f" % (lonMin,lonMax,lonRange))

        # Determine the latitude and longitude fill masks, so we can restore the 
        # fill values after we have scaled...

        latMask = latitude.mask
        lonMask = longitude.mask

        # Check if the geolocation is in radians, convert to degrees
        # FIXME : This information is likely conveyed by whether the 
        # FIXME :     geolocation short-name is *-GEO-TC (degrees) or
        # FIXME :     *-RGEO_TC (radians).
        if (lonRange < 2.*np.pi) :
            LOG.debug("Geolocation is in radians, convert to degrees...")
            latitude = np.degrees(latitude)
            longitude = np.degrees(longitude)
        
            latMin,latMax = np.min(latitude),np.max(latitude)
            latRange = latMax-latMin
            lonMin,lonMax = np.min(longitude),np.max(longitude)
            lonRange = lonMax-lonMin

            LOG.debug("New min,max,range of latitide: %f %f %f" % (latMin,latMax,latRange))
            LOG.debug("New min,max,range of longitude: %f %f %f" % (lonMin,lonMax,lonRange))

        # Restore fill values to masked pixels in geolocation

        geoFillValue = trimObj.sdrTypeFill['VDNE_FLOAT64_FILL'][latitude.dtype.name]
        latitude = ma.array(latitude,mask=latMask,fill_value=geoFillValue)
        latitude = latitude.filled()

        geoFillValue = trimObj.sdrTypeFill['VDNE_FLOAT64_FILL'][longitude.dtype.name]
        longitude = ma.array(longitude,mask=lonMask,fill_value=geoFillValue)
        longitude = longitude.filled()

        # Parse the geolocation asc file to get struct information which will be 
        # written to the ancillary asc files

        geoAscFileName = path.join(inDir,URID+".asc")
        LOG.debug("\nOpening %s..." % (geoAscFileName))

        geoAscFile = open(geoAscFileName,'rt')

        #RangeDateTimeStr =  _getAscLine(geoAscFile,"RangeDateTime")
        RangeDateTimeStr =  _getAscLine(geoAscFile,"ObservedDateTime")
        RangeDateTimeStr =  string.replace(RangeDateTimeStr,"ObservedDateTime","RangeDateTime")
        GRingLatitudeStr =  _getAscStructs(geoAscFile,"GRingLatitude",12)
        GRingLongitudeStr =  _getAscStructs(geoAscFile,"GRingLongitude",12)

        geoAscFile.close()

        # Loop through the required NCEP datasets and create the blobs.
        # FIXME : Handle pathological geolocation cases

        firstGranule = True

        for dSet in masksCollShortNames :
        
            LOG.debug("\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n")
            LOG.debug("Processing dataset %s for %s" % (NCEP_shortNameToBlobName[dSet],dSet))

            if (dSet=='VIIRS-ANC-Surf-Ht-Mod-Gran') and terrainCorrectedGeo :

                data[:,:] = terrain[:,:]

            else :

                # FIXME : Account for dateline and pole crossings...

                # Massage the NCEP data array a bit...
                NCEP_anc = np.array(NCEP_globalGridData[dSet])[::-1,:]
                NCEP_anc = np.roll(NCEP_anc,360)

                if (firstGranule) :

                    LOG.debug("\nGranulating %s ..." % (dSet))
                    LOG.debug("latitide,longitude shapes: %s, %s"%(str(latitude.shape) , str(longitude.shape)))
                    LOG.debug("NCEP_anc.shape = %s" % (str(NCEP_anc.shape)))
                    LOG.debug("gridLat.shape = %s" % (str(gridLat.shape)))
                    LOG.debug("gridLon.shape = %s" % (str(gridLon.shape)))

                    LOG.debug("min of NCEP_anc  = %r"%(np.min(NCEP_anc)))
                    LOG.debug("max of NCEP_anc  = %r"%(np.max(NCEP_anc)))

                    data,dataIdx = _grid2Gran(np.ravel(latitude),
                                              np.ravel(longitude),
                                              NCEP_anc.astype(np.float64),
                                              gridLat,
                                              gridLon)

                    data = data.reshape(latitude.shape)
                    dataIdx = dataIdx.reshape(latitude.shape)
                    firstGranule = False
                    LOG.debug("Shape of first granulated %s data is %s" % (dSet,np.shape(data)))
                    LOG.debug("Shape of first granulated %s dataIdx is %s" % (dSet,np.shape(dataIdx)))

                else :

                    LOG.debug("Granulating %s using existing data indices." % (dSet))
                    NCEP_anc = np.ravel(NCEP_anc)
                    data = np.ravel(NCEP_anc)[np.ravel(dataIdx)]
                    data = data.reshape(latitude.shape)
                    LOG.debug("Shape of subsequent granulated %s is %s" % (dSet,np.shape(data)))

            # Fill the required pixel trim rows in the granulated NCEP data with 
            # the ONBOARD_PT_FILL value for the correct data type

            fillValue = trimObj.sdrTypeFill['ONBOARD_PT_FILL'][data.dtype.name]        
            data = ma.array(data,mask=modTrimMask,fill_value=fillValue)
            data = data.filled()

            # Create new NCEP ancillary blob, and copy granulated data to it

            endian = ancEndian
            xmlName = path.join(ADL_HOME,'xml/VIIRS',NCEP_shortNameToXmlName[dSet])

            # Create a new URID to be used in making the asc filenames

            URID_dict = _getURID()

            URID = URID_dict['URID']
            creationDate_nousecStr = URID_dict['creationDate_nousecStr']
            creationDateStr = URID_dict['creationDateStr']

            # Create a new directory in the input directory for the new ancillary
            # asc and blob files

            blobDir = inDir

            ascFileName = path.join(blobDir,URID+'.asc')
            blobName = path.join(blobDir,URID+'.'+dSet)

            LOG.debug("ascFileName : %s" % (ascFileName))
            LOG.debug("blobName : %s" % (blobName))

            # Create a new ancillary blob, and copy the data to it.
            newNCEPblobObj = adl_blob.create(xmlName, blobName, endian=endian, overwrite=True)
            newNCEPblobArrObj = newNCEPblobObj.as_arrays()

            blobData = getattr(newNCEPblobArrObj,'data')
            blobData[:,:] = data[:,:]

            # Make a new NCEP asc file from the template, and substitute for the various tags

            ascTemplateFileName = path.join(ADL_ASC_TEMPLATES,"VIIRS-ANC_Template.asc")

            LOG.debug("Creating new asc file\n%s\nfrom template\n%s" % (ascFileName,ascTemplateFileName))
            
            ANC_fileList = gridBlobFiles
            for idx in range(len(ANC_fileList)) :
                ANC_fileList[idx] = path.basename(ANC_fileList[idx])
            ANC_fileList.sort()
            ancGroupRecipe = '    ("N_Anc_Filename" STRING EQ "%s")'
            ancFileStr = "%s" % ("\n").join([ancGroupRecipe % (str(files)) for files in ANC_fileList])

            LOG.debug("RangeDateTimeStr = %s\n" % (RangeDateTimeStr))
            LOG.debug("GRingLatitudeStr = \n%s\n" % (GRingLatitudeStr))
            LOG.debug("GRingLongitudeStr = \n%s\n" % (GRingLongitudeStr))


            ascTemplateFile = open(ascTemplateFileName,"rt") # Open template file for reading
            ascFile = open(ascFileName,"wt") # create a new text file

            LOG.debug("Template file %s is %r with mode %s" %(ascTemplateFileName,'not open' if ascTemplateFile.closed else 'open',ascTemplateFile.mode))

            LOG.debug("New file %s is %r with mode %s" %(ascFileName,'not open' if ascFile.closed else 'open',ascFile.mode))

            for line in ascTemplateFile.readlines():
               line = line.replace("CSPP_URID",URID)
               line = line.replace("CSPP_CREATIONDATETIME_NOUSEC",creationDate_nousecStr)
               line = line.replace("CSPP_ANC_BLOB_FULLPATH",blobName)
               line = line.replace("CSPP_ANC_COLLECTION_SHORT_NAME",dSet)
               line = line.replace("CSPP_GRANULE_ID",N_Granule_ID)
               line = line.replace("CSPP_CREATIONDATETIME",creationDateStr)
               line = line.replace("  CSPP_RANGE_DATE_TIME",RangeDateTimeStr)
               line = line.replace("  CSPP_GRINGLATITUDE",GRingLatitudeStr)
               line = line.replace("  CSPP_GRINGLONGITUDE",GRingLongitudeStr)
               line = line.replace("    CSPP_ANC_SOURCE_FILES",ancFileStr)
               ascFile.write(line) 

            ascFile.close()
            ascTemplateFile.close()


def _granulate_NCEP_gridBlobs_alt(inDir,geoDicts, gridBlobFiles):
    '''Granulates the input gridded blob files into the required NCEP granulated datasets.'''

    global ancEndian 

    CSPP_RT_HOME = os.getenv('CSPP_RT_HOME')
    ANC_SCRIPTS_PATH = path.join(CSPP_RT_HOME,'viirs')
    CSPP_RT_ANC_CACHE_DIR = os.getenv('CSPP_RT_ANC_CACHE_DIR')
    #CSPP_RT_ANC_HOME = os.getenv('CSPP_RT_ANC_HOME')
    csppPython = os.getenv('PY')
    

    ADL_ASC_TEMPLATES = path.join(CSPP_RT_ANC_HOME,'asc_templates')
    
    # Collection shortnames of the required NCEP ancillary datasets
    # FIXME : Poll ADL/cfg/ProEdrViirsCM_CFG.xml for this information

    masksCollShortNamesDict = {'VCM' : [
                           'VIIRS-ANC-Preci-Wtr-Mod-Gran',
                           'VIIRS-ANC-Temp-Surf2M-Mod-Gran',
                           'VIIRS-ANC-Wind-Speed-Mod-Gran',
                           'VIIRS-ANC-Surf-Ht-Mod-Gran'
                          ],
                          'AOT' : [
                            'VIIRS-ANC-Preci-Wtr-Mod-Gran',
                            'VIIRS-ANC-Temp-Surf2M-Mod-Gran',
                            'VIIRS-ANC-Wind-Speed-Mod-Gran',
                            'VIIRS-ANC-Wind-Direction-Mod-Gran',
                            'VIIRS-ANC-Press-Surf-Mod-Gran',
                            'VIIRS-ANC-Tot-Col-Mod-Gran'
                              ]
                          }

    masksCollShortNames = [
                           'VIIRS-ANC-Preci-Wtr-Mod-Gran',
                           'VIIRS-ANC-Temp-Surf2M-Mod-Gran',
                           'VIIRS-ANC-Wind-Speed-Mod-Gran',
                           'VIIRS-ANC-Surf-Ht-Mod-Gran',
                           'VIIRS-ANC-Wind-Direction-Mod-Gran',
                           'VIIRS-ANC-Press-Surf-Mod-Gran',
                           'VIIRS-ANC-Tot-Col-Mod-Gran'
                          ]

    # Dictionary relating the required NCEP collection short names and the 
    # NCEP gridded blob dataset names

    NCEP_shortNameToBlobName = {}
    NCEP_shortNameToBlobName['VIIRS-ANC-Preci-Wtr-Mod-Gran'] = 'totalPrecipitableWater'
    NCEP_shortNameToBlobName['VIIRS-ANC-Temp-Surf2M-Mod-Gran'] = 'surfaceTemperature'
    NCEP_shortNameToBlobName['VIIRS-ANC-Wind-Speed-Mod-Gran'] = ['uComponentOfWind', 'vComponentOfWind']
    NCEP_shortNameToBlobName['VIIRS-ANC-Surf-Ht-Mod-Gran'] = 'surfaceGeopotentialHeight'
    NCEP_shortNameToBlobName['VIIRS-ANC-Wind-Direction-Mod-Gran'] = 'windDirection'  # Dummy Entry
    NCEP_shortNameToBlobName['VIIRS-ANC-Press-Surf-Mod-Gran'] = 'surfacePressure'
    NCEP_shortNameToBlobName['VIIRS-ANC-Tot-Col-Mod-Gran'] = 'totalColumnOzone'

    NCEP_shortNameToXmlName = {}
    NCEP_shortNameToXmlName['VIIRS-ANC-Preci-Wtr-Mod-Gran'] = 'VIIRS_ANC_PRECI_WTR_MOD_GRAN.xml'
    NCEP_shortNameToXmlName['VIIRS-ANC-Temp-Surf2M-Mod-Gran'] = 'VIIRS_ANC_TEMP_SURF2M_MOD_GRAN.xml'
    NCEP_shortNameToXmlName['VIIRS-ANC-Wind-Speed-Mod-Gran'] = 'VIIRS_ANC_WIND_SPEED_MOD_GRAN.xml'
    NCEP_shortNameToXmlName['VIIRS-ANC-Surf-Ht-Mod-Gran'] = 'VIIRS_ANC_SURF_HT_MOD_GRAN.xml'
    NCEP_shortNameToXmlName['VIIRS-ANC-Wind-Direction-Mod-Gran'] = 'VIIRS_ANC_WIND_DIRECTION_MOD_GRAN.xml'
    NCEP_shortNameToXmlName['VIIRS-ANC-Press-Surf-Mod-Gran'] = 'VIIRS_ANC_PRESS_SURF_MOD_GRAN.xml'
    NCEP_shortNameToXmlName['VIIRS-ANC-Tot-Col-Mod-Gran'] = 'VIIRS_ANC_TOT_COL_MOD_GRAN.xml'

    # Moderate resolution trim table arrays. These are 
    # bool arrays, and the trim pixels are set to True.

    trimObj = ViirsData.ViirsTrimTable()
    modTrimMask = trimObj.createModTrimArray(nscans=48,trimType=bool)

    # Open the NCEP gridded blob file

    ncepXmlFile = path.join(ADL_HOME,'xml/ANC/NCEP_ANC_Int.xml')
    # FIXME : Should be using two NCEP blob files, and averaging
    gridBlobFile = gridBlobFiles[0]

    if os.path.exists(ncepXmlFile):
        LOG.info("We are using for %s: %s,%s" %('NCEP-ANC-Int',ncepXmlFile,gridBlobFile))
    
    endian = ancEndian

    ncepBlobObj = adl_blob.map(ncepXmlFile,gridBlobFile, endian=endian)
    ncepBlobArrObj = ncepBlobObj.as_arrays()

    LOG.debug("%s...\n%r" % (gridBlobFile,ncepBlobArrObj._fields))

    # Save the global grids of the required datasets into a dictionary...

    NCEP_globalGridData = {}

    # Precipitable water
    blobDsetName = NCEP_shortNameToBlobName['VIIRS-ANC-Preci-Wtr-Mod-Gran']
    NCEP_globalGridData['VIIRS-ANC-Preci-Wtr-Mod-Gran'] = \
            getattr(ncepBlobArrObj,blobDsetName).astype('float')
    LOG.info("Shape of dataset %s is %s" % (blobDsetName,np.shape(NCEP_globalGridData['VIIRS-ANC-Preci-Wtr-Mod-Gran'])))

    # 2m Surface temperature
    blobDsetName = NCEP_shortNameToBlobName['VIIRS-ANC-Temp-Surf2M-Mod-Gran']
    NCEP_globalGridData['VIIRS-ANC-Temp-Surf2M-Mod-Gran'] = \
            getattr(ncepBlobArrObj,blobDsetName).astype('float')
    LOG.info("Shape of dataset %s is %s" % (blobDsetName,np.shape(NCEP_globalGridData['VIIRS-ANC-Temp-Surf2M-Mod-Gran'])))

    # Omnidirectional wind speed
    blobDsetName = NCEP_shortNameToBlobName['VIIRS-ANC-Wind-Speed-Mod-Gran'][0]
    uWind = getattr(ncepBlobArrObj,blobDsetName).astype('float')
    blobDsetName = NCEP_shortNameToBlobName['VIIRS-ANC-Wind-Speed-Mod-Gran'][1]
    vWind = getattr(ncepBlobArrObj,blobDsetName).astype('float')
    windSpeed = np.sqrt(uWind*uWind + vWind*vWind)
    NCEP_globalGridData['VIIRS-ANC-Wind-Speed-Mod-Gran'] = windSpeed
    LOG.info("Shape of dataset %s is %s" % ('windSpeed',np.shape(NCEP_globalGridData['VIIRS-ANC-Wind-Speed-Mod-Gran'])))

    # Wind direction
    windDirection_rad  = np.pi/2. - np.arctan2(vWind, uWind)
    negWindDirIdx = np.where(windDirection_rad < 0.)
    windDirection_rad[negWindDirIdx] = windDirection_rad[negWindDirIdx] + 2.*np.pi
    windDirection = np.degrees(windDirection_rad)
    NCEP_globalGridData['VIIRS-ANC-Wind-Direction-Mod-Gran'] = windDirection
    LOG.info("Shape of dataset %s is %s" % ('windSpeed',np.shape(NCEP_globalGridData['VIIRS-ANC-Wind-Direction-Mod-Gran'])))

    # Terrain height
    # FIXME : VIIRS-ANC-Surf-Ht-Mod-Gran is supposed to come from the Terrain-Eco-ANC-Tile 
    # FIXME :   gridded data. The granulated terrain height from this source is provided in
    # FIXME :   the VIIRS-MOD-RGEO-TC from the VIIRS SDR controller. If we get the non 
    # FIXME :   terrain correction geolocation (VIIRS-MOD-RGEO) we can use the surface
    # FIXME :   geopotential height, which is approximately the same.
    blobDsetName = NCEP_shortNameToBlobName['VIIRS-ANC-Surf-Ht-Mod-Gran']
    NCEP_globalGridData['VIIRS-ANC-Surf-Ht-Mod-Gran'] = \
            getattr(ncepBlobArrObj,blobDsetName).astype('float')
    LOG.info("Shape of dataset %s is %s" % (blobDsetName,np.shape(NCEP_globalGridData['VIIRS-ANC-Surf-Ht-Mod-Gran'])))

    # Surface Pressure
    blobDsetName = NCEP_shortNameToBlobName['VIIRS-ANC-Press-Surf-Mod-Gran']
    NCEP_globalGridData['VIIRS-ANC-Press-Surf-Mod-Gran'] = \
            getattr(ncepBlobArrObj,blobDsetName).astype('float')
    LOG.info("Shape of dataset %s is %s" % (blobDsetName,np.shape(NCEP_globalGridData['VIIRS-ANC-Press-Surf-Mod-Gran'])))
    # TODO : Terrain correction of surface pressure in ADL/ANC/VIIRS/SurfPres/src/ProAncViirsGranulateSurfPres.cpp

    # Total Column Ozone
    blobDsetName = NCEP_shortNameToBlobName['VIIRS-ANC-Tot-Col-Mod-Gran']
    NCEP_globalGridData['VIIRS-ANC-Tot-Col-Mod-Gran'] = \
            getattr(ncepBlobArrObj,blobDsetName).astype('float')
    LOG.info("Shape of dataset %s is %s" % (blobDsetName,np.shape(NCEP_globalGridData['VIIRS-ANC-Tot-Col-Mod-Gran'])))

    # Contruct a default 0.5 degree grid...

    degInc = 0.5
    grids = np.mgrid[-90.:90.+degInc:degInc,-180.:180.:degInc]
    gridLat,gridLon = grids[0],grids[1]

    # Loop through the geolocation files and granulate...

    for dicts in geoDicts :

        URID = dicts['URID']
        geo_Collection_ShortName = dicts['N_Collection_Short_Name']
        N_Granule_ID = dicts['N_Granule_ID']
        ObservedStartTimeObj = dicts['ObservedStartTime']
        geoFiles = glob('%s/%s*' % (inDir,URID))
        geoFiles.sort()

        LOG.debug("\n###########################")
        LOG.debug("  Geolocation Information  ")
        LOG.debug("###########################")
        LOG.debug("N_Granule_ID : %r" % (N_Granule_ID))
        LOG.debug("ObservedStartTime : %s" % (ObservedStartTimeObj.__str__()))
        LOG.debug("N_Collection_Short_Name : %s" %(geo_Collection_ShortName))
        LOG.debug("URID : %r" % (URID))
        LOG.debug("geoFiles : %r" % (geoFiles))
        LOG.debug("###########################\n")

        # Do we have terrain corrected geolocation?

        terrainCorrectedGeo = True if 'GEO-TC' in geo_Collection_ShortName else False

        # Do we have long or short style geolocation field names?

        if (geo_Collection_ShortName=='VIIRS-MOD-GEO-TC' or geo_Collection_ShortName=='VIIRS-MOD-RGEO') :
            longFormGeoNames = True
            LOG.debug("We have long form geolocation names")
        elif (geo_Collection_ShortName=='VIIRS-MOD-GEO' or geo_Collection_ShortName=='VIIRS-MOD-RGEO-TC') :
            LOG.debug("We have short form geolocation names")
            longFormGeoNames = False
        else :
            LOG.error("Invalid geolocation shortname: %s",geo_Collection_ShortName)
            return -1

        # Get the geolocation xml file

        geoXmlFile = "%s.xml" % (string.replace(geo_Collection_ShortName,'-','_'))
        geoXmlFile = path.join(ADL_HOME,'xml/VIIRS',geoXmlFile)
        if os.path.exists(geoXmlFile):
            LOG.debug("We are using for %s: %s,%s" %(geo_Collection_ShortName,geoXmlFile,geoFiles[0]))

        # Open the geolocation blob and get the latitude and longitude

        endian=sdrEndian

        geoBlobObj = adl_blob.map(geoXmlFile,geoFiles[0], endian=endian)
        geoBlobArrObj = geoBlobObj.as_arrays()

        # If we have the terrain corrected geolocation, get the terrain height

        if terrainCorrectedGeo :
            terrainHeight = geoBlobArrObj.height[:,:]

        # Get scan_mode to find any bad scans

        scanMode = geoBlobArrObj.scan_mode[:]
        badScanIdx = np.where(scanMode==254)[0]
        LOG.debug("Bad Scans: %r" % (badScanIdx))

        # Detemine the min, max and range of the latitude and longitude, 
        # taking care to exclude any fill values.

        if longFormGeoNames :
            latitude = getattr(geoBlobArrObj,'latitude').astype('float')
            longitude = getattr(geoBlobArrObj,'longitude').astype('float')
        else :
            latitude = getattr(geoBlobArrObj,'lat').astype('float')
            longitude = getattr(geoBlobArrObj,'lon').astype('float')
        
        if terrainCorrectedGeo :
            terrain = getattr(geoBlobArrObj,'height').astype('float')

        latitude = ma.masked_less(latitude,-800.)
        latMin,latMax = np.min(latitude),np.max(latitude)
        latRange = latMax-latMin

        longitude = ma.masked_less(longitude,-800.)
        lonMin,lonMax = np.min(longitude),np.max(longitude)
        lonRange = lonMax-lonMin

        LOG.debug("min,max,range of latitide: %f %f %f" % (latMin,latMax,latRange))
        LOG.debug("min,max,range of longitude: %f %f %f" % (lonMin,lonMax,lonRange))

        # Determine the latitude and longitude fill masks, so we can restore the 
        # fill values after we have scaled...

        latMask = latitude.mask
        lonMask = longitude.mask

        # Check if the geolocation is in radians, convert to degrees
        # FIXME : This information is likely conveyed by whether the 
        # FIXME :     geolocation short-name is *-GEO-TC (degrees) or
        # FIXME :     *-RGEO_TC (radians).
        if (lonRange < 2.*np.pi) :
            LOG.debug("Geolocation is in radians, convert to degrees...")
            latitude = np.degrees(latitude)
            longitude = np.degrees(longitude)
        
            latMin,latMax = np.min(latitude),np.max(latitude)
            latRange = latMax-latMin
            lonMin,lonMax = np.min(longitude),np.max(longitude)
            lonRange = lonMax-lonMin

            LOG.debug("New min,max,range of latitide: %f %f %f" % (latMin,latMax,latRange))
            LOG.debug("New min,max,range of longitude: %f %f %f" % (lonMin,lonMax,lonRange))

        # Restore fill values to masked pixels in geolocation

        geoFillValue = trimObj.sdrTypeFill['VDNE_FLOAT64_FILL'][latitude.dtype.name]
        latitude = ma.array(latitude,mask=latMask,fill_value=geoFillValue)
        latitude = latitude.filled()

        geoFillValue = trimObj.sdrTypeFill['VDNE_FLOAT64_FILL'][longitude.dtype.name]
        longitude = ma.array(longitude,mask=lonMask,fill_value=geoFillValue)
        longitude = longitude.filled()

        # Parse the geolocation asc file to get struct information which will be 
        # written to the ancillary asc files

        geoAscFileName = path.join(inDir,URID+".asc")
        LOG.debug("\nOpening %s..." % (geoAscFileName))

        geoAscFile = open(geoAscFileName,'rt')

        #RangeDateTimeStr =  _getAscLine(geoAscFile,"RangeDateTime")
        RangeDateTimeStr =  _getAscLine(geoAscFile,"ObservedDateTime")
        RangeDateTimeStr =  string.replace(RangeDateTimeStr,"ObservedDateTime","RangeDateTime")
        GRingLatitudeStr =  _getAscStructs(geoAscFile,"GRingLatitude",12)
        GRingLongitudeStr =  _getAscStructs(geoAscFile,"GRingLongitude",12)

        geoAscFile.close()

        # Loop through the required NCEP datasets and create the blobs.
        # FIXME : Handle pathological geolocation cases

        firstGranule = True

        for dSet in masksCollShortNames :
        
            LOG.debug("\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n")
            LOG.debug("Processing dataset %s for %s" % (NCEP_shortNameToBlobName[dSet],dSet))

            if (dSet=='VIIRS-ANC-Surf-Ht-Mod-Gran') and terrainCorrectedGeo :

                data[:,:] = terrain[:,:]

            else :

                # FIXME : Account for dateline and pole crossings...

                # Massage the NCEP data array a bit...
                NCEP_anc = np.array(NCEP_globalGridData[dSet])[::-1,:]
                NCEP_anc = np.roll(NCEP_anc,360)

                if (firstGranule) :

                    LOG.debug("\nGranulating %s ..." % (dSet))
                    LOG.debug("latitide,longitude shapes: %s, %s"%(str(latitude.shape) , str(longitude.shape)))
                    LOG.debug("NCEP_anc.shape = %s" % (str(NCEP_anc.shape)))
                    LOG.debug("gridLat.shape = %s" % (str(gridLat.shape)))
                    LOG.debug("gridLon.shape = %s" % (str(gridLon.shape)))

                    LOG.debug("min of NCEP_anc  = "%(np.min(NCEP_anc)))
                    LOG.debug("max of NCEP_anc  = "%(np.max(NCEP_anc)))

                    data,dataIdx = _grid2Gran(np.ravel(latitude),
                                              np.ravel(longitude),
                                              NCEP_anc.astype(np.float64),
                                              gridLat,
                                              gridLon)

                    data = data.reshape(latitude.shape)
                    dataIdx = dataIdx.reshape(latitude.shape)
                    firstGranule = False
                    LOG.debug("Shape of first granulated %s data is %s" % (dSet,np.shape(data)))
                    LOG.debug("Shape of first granulated %s dataIdx is %s" % (dSet,np.shape(dataIdx)))

                else :

                    LOG.debug("Granulating %s using existing data indices." % (dSet))
                    NCEP_anc = np.ravel(NCEP_anc)
                    data = np.ravel(NCEP_anc)[np.ravel(dataIdx)]
                    data = data.reshape(latitude.shape)
                    LOG.debug("Shape of subsequent granulated %s is %s" % (dSet,np.shape(data)))

            # Fill the required pixel trim rows in the granulated NCEP data with 
            # the ONBOARD_PT_FILL value for the correct data type

            fillValue = trimObj.sdrTypeFill['ONBOARD_PT_FILL'][data.dtype.name]        
            data = ma.array(data,mask=modTrimMask,fill_value=fillValue)
            data = data.filled()

            # Create new NCEP ancillary blob, and copy granulated data to it

            endian = ancEndian
            xmlName = path.join(ADL_HOME,'xml/VIIRS',NCEP_shortNameToXmlName[dSet])

            # Create a new URID to be used in making the asc filenames

            URID_dict = _getURID()

            URID = URID_dict['URID']
            creationDate_nousecStr = URID_dict['creationDate_nousecStr']
            creationDateStr = URID_dict['creationDateStr']

            # Create a new directory in the input directory for the new ancillary
            # asc and blob files

            blobDir = inDir

            ascFileName = path.join(blobDir,URID+'.asc')
            blobName = path.join(blobDir,URID+'.'+dSet)

            LOG.debug("ascFileName : %s" % (ascFileName))
            LOG.debug("blobName : %s" % (blobName))

            # Create a new ancillary blob, and copy the data to it.
            newNCEPblobObj = adl_blob.create(xmlName, blobName, endian=endian, overwrite=True)
            newNCEPblobArrObj = newNCEPblobObj.as_arrays()

            blobData = getattr(newNCEPblobArrObj,'data')
            blobData[:,:] = data[:,:]

            # Make a new NCEP asc file from the template, and substitute for the various tags

            ascTemplateFileName = path.join(ADL_ASC_TEMPLATES,"VIIRS-ANC_Template.asc")

            LOG.debug("Creating new asc file\n%s\nfrom template\n%s" % (ascFileName,ascTemplateFileName))
            
            ANC_fileList = gridBlobFiles
            for idx in range(len(ANC_fileList)) :
                ANC_fileList[idx] = path.basename(ANC_fileList[idx])
            ANC_fileList.sort()
            ancGroupRecipe = '    ("N_Anc_Filename" STRING EQ "%s")'
            ancFileStr = "%s" % ("\n").join([ancGroupRecipe % (str(files)) for files in ANC_fileList])

            LOG.debug("RangeDateTimeStr = %s\n" % (RangeDateTimeStr))
            LOG.debug("GRingLatitudeStr = \n%s\n" % (GRingLatitudeStr))
            LOG.debug("GRingLongitudeStr = \n%s\n" % (GRingLongitudeStr))


            ascTemplateFile = open(ascTemplateFileName,"rt") # Open template file for reading
            ascFile = open(ascFileName,"wt") # create a new text file

            LOG.debug("Template file %s is %r with mode %s" %(ascTemplateFileName,'not open' if ascTemplateFile.closed else 'open',ascTemplateFile.mode))

            LOG.debug("New file %s is %r with mode %s" %(ascFileName,'not open' if ascFile.closed else 'open',ascFile.mode))

            for line in ascTemplateFile.readlines():
               line = line.replace("CSPP_URID",URID)
               line = line.replace("CSPP_CREATIONDATETIME_NOUSEC",creationDate_nousecStr)
               line = line.replace("CSPP_ANC_BLOB_FULLPATH",blobName)
               line = line.replace("CSPP_ANC_COLLECTION_SHORT_NAME",dSet)
               line = line.replace("CSPP_GRANULE_ID",N_Granule_ID)
               line = line.replace("CSPP_CREATIONDATETIME",creationDateStr)
               line = line.replace("  CSPP_RANGE_DATE_TIME",RangeDateTimeStr)
               line = line.replace("  CSPP_GRINGLATITUDE",GRingLatitudeStr)
               line = line.replace("  CSPP_GRINGLONGITUDE",GRingLongitudeStr)
               line = line.replace("    CSPP_ANC_SOURCE_FILES",ancFileStr)
               ascFile.write(line) 

            ascFile.close()
            ascTemplateFile.close()


def _granulate_NAAPS_gridBlobs(inDir,geoDicts, gridBlobFiles):
    '''Granulates the input gridded blob files into the required NAAPS granulated datasets.'''

    global ancEndian 

    CSPP_RT_HOME = os.getenv('CSPP_RT_HOME')
    ANC_SCRIPTS_PATH = path.join(CSPP_RT_HOME,'masks')
    CSPP_RT_ANC_CACHE_DIR = os.getenv('CSPP_RT_ANC_CACHE_DIR')
    #CSPP_RT_ANC_HOME = os.getenv('CSPP_RT_ANC_HOME')
    csppPython = os.getenv('PY')
    

    ADL_ASC_TEMPLATES = path.join(CSPP_RT_ANC_HOME,'asc_templates')
    
    # Collection shortnames of the required NAAPS ancillary datasets
    # FIXME : Poll ADL/cfg/ProEdrViirsCM_CFG.xml for this information

    masksCollShortNames = [
                           'VIIRS-ANC-Optical-Depth-Mod-Gran'
                          ]

    # Dictionary relating the required NAAPS collection short names and the 
    # NAAPS gridded blob dataset names

    NAAPS_shortNameToBlobName = {}
    NAAPS_shortNameToBlobName['VIIRS-ANC-Optical-Depth-Mod-Gran'] = 'aotGrid'

    NAAPS_shortNameToXmlName = {}
    NAAPS_shortNameToXmlName['VIIRS-ANC-Optical-Depth-Mod-Gran'] = 'VIIRS_ANC_OPTICAL_DEPTH_MOD_GRAN.xml'

    # Moderate resolution trim table arrays. These are 
    # bool arrays, and the trim pixels are set to True.

    trimObj = ViirsData.ViirsTrimTable()
    modTrimMask = trimObj.createModTrimArray(nscans=48,trimType=bool)

    # Open the NAAPS gridded blob file

    naapsXmlFile = path.join(ADL_HOME,'xml/ANC/NAAPS_ANC_Int.xml')
    # FIXME : Should be using two NAAPS blob files, and averaging
    gridBlobFile = gridBlobFiles[0]

    if os.path.exists(naapsXmlFile):
        LOG.info("We are using for %s: %s,%s" %('NAAPS-ANC-Int',naapsXmlFile,gridBlobFile))
    
    endian = ancEndian

    naapsBlobObj = adl_blob.map(naapsXmlFile,gridBlobFile, endian=adl_blob.BIG_ENDIAN)
    naapsBlobArrObj = naapsBlobObj.as_arrays()

    LOG.debug("%s...\n%r" % (gridBlobFile,naapsBlobArrObj._fields))

    # Save the global grids of the required datasets into a dictionary...

    NAAPS_globalGridData = {}

    # Precipitable water
    blobDsetName = NAAPS_shortNameToBlobName['VIIRS-ANC-Optical-Depth-Mod-Gran']
    NAAPS_globalGridData['VIIRS-ANC-Optical-Depth-Mod-Gran'] = \
            getattr(naapsBlobArrObj,blobDsetName).astype('float')
    LOG.info("Shape of dataset %s is %s" % (blobDsetName,np.shape(NAAPS_globalGridData['VIIRS-ANC-Optical-Depth-Mod-Gran'])))


    # Contruct a default 0.5 degree grid...

    degInc = 0.5
    grids = np.mgrid[-90.:90.+degInc:degInc,-180.:180.:degInc]
    gridLat,gridLon = grids[0],grids[1]

    # Loop through the geolocation files and granulate...

    for dicts in geoDicts :

        URID = dicts['URID']
        geo_Collection_ShortName = dicts['N_Collection_Short_Name']
        N_Granule_ID = dicts['N_Granule_ID']
        ObservedStartTimeObj = dicts['ObservedStartTime']
        geoFiles = glob('%s/%s*' % (inDir,URID))
        geoFiles.sort()

        LOG.debug("\n###########################")
        LOG.debug("  Geolocation Information  ")
        LOG.debug("###########################")
        LOG.debug("N_Granule_ID : %r" % (N_Granule_ID))
        LOG.debug("ObservedStartTime : %s" % (ObservedStartTimeObj.__str__()))
        LOG.debug("N_Collection_Short_Name : %s" %(geo_Collection_ShortName))
        LOG.debug("URID : %r" % (URID))
        LOG.debug("geoFiles : %r" % (geoFiles))
        LOG.debug("###########################\n")

        # Do we have terrain corrected geolocation?

        terrainCorrectedGeo = True if 'GEO-TC' in geo_Collection_ShortName else False

        # Do we have long or short style geolocation field names?

        if (geo_Collection_ShortName=='VIIRS-MOD-GEO-TC' or geo_Collection_ShortName=='VIIRS-MOD-RGEO') :
            longFormGeoNames = True
            LOG.debug("We have long form geolocation names")
        elif (geo_Collection_ShortName=='VIIRS-MOD-GEO' or geo_Collection_ShortName=='VIIRS-MOD-RGEO-TC') :
            LOG.debug("We have short form geolocation names")
            longFormGeoNames = False
        else :
            LOG.error("Invalid geolocation shortname: %s"% (geo_Collection_ShortName))
            return -1

        # Get the geolocation xml file

        geoXmlFile = "%s.xml" % (string.replace(geo_Collection_ShortName,'-','_'))
        geoXmlFile = path.join(ADL_HOME,'xml/VIIRS',geoXmlFile)
        if os.path.exists(geoXmlFile):
            LOG.debug("We are using for %s: %s,%s" %(geo_Collection_ShortName,geoXmlFile,geoFiles[0]))

        # Open the geolocation blob and get the latitude and longitude

        endian=sdrEndian

        geoBlobObj = adl_blob.map(geoXmlFile,geoFiles[0], endian=endian)
        geoBlobArrObj = geoBlobObj.as_arrays()

        # If we have the terrain corrected geolocation, get the terrain height

        if terrainCorrectedGeo :
            terrainHeight = geoBlobArrObj.height[:,:]

        # Get scan_mode to find any bad scans

        scanMode = geoBlobArrObj.scan_mode[:]
        badScanIdx = np.where(scanMode==254)[0]
        LOG.debug("Bad Scans: %r" % (badScanIdx))

        # Detemine the min, max and range of the latitude and longitude, 
        # taking care to exclude any fill values.

        if longFormGeoNames :
            latitude = getattr(geoBlobArrObj,'latitude').astype('float')
            longitude = getattr(geoBlobArrObj,'longitude').astype('float')
        else :
            latitude = getattr(geoBlobArrObj,'lat').astype('float')
            longitude = getattr(geoBlobArrObj,'lon').astype('float')
        
        latitude = ma.masked_less(latitude,-800.)
        latMin,latMax = np.min(latitude),np.max(latitude)
        latRange = latMax-latMin

        longitude = ma.masked_less(longitude,-800.)
        lonMin,lonMax = np.min(longitude),np.max(longitude)
        lonRange = lonMax-lonMin

        LOG.debug("min,max,range of latitide: %f %f %f" % (latMin,latMax,latRange))
        LOG.debug("min,max,range of longitude: %f %f %f" % (lonMin,lonMax,lonRange))

        # Determine the latitude and longitude fill masks, so we can restore the 
        # fill values after we have scaled...

        latMask = latitude.mask
        lonMask = longitude.mask

        # Check if the geolocation is in radians, convert to degrees
        # FIXME : This information is likely conveyed by whether the 
        # FIXME :     geolocation short-name is *-GEO-TC (degrees) or
        # FIXME :     *-RGEO_TC (radians).
        if (lonRange < 2.*np.pi) :
            LOG.debug("Geolocation is in radians, convert to degrees...")
            latitude = np.degrees(latitude)
            longitude = np.degrees(longitude)
        
            latMin,latMax = np.min(latitude),np.max(latitude)
            latRange = latMax-latMin
            lonMin,lonMax = np.min(longitude),np.max(longitude)
            lonRange = lonMax-lonMin

            LOG.debug("New min,max,range of latitide: %f %f %f" % (latMin,latMax,latRange))
            LOG.debug("New min,max,range of longitude: %f %f %f" % (lonMin,lonMax,lonRange))

        # Restore fill values to masked pixels in geolocation

        geoFillValue = trimObj.sdrTypeFill['VDNE_FLOAT64_FILL'][latitude.dtype.name]
        latitude = ma.array(latitude,mask=latMask,fill_value=geoFillValue)
        latitude = latitude.filled()

        geoFillValue = trimObj.sdrTypeFill['VDNE_FLOAT64_FILL'][longitude.dtype.name]
        longitude = ma.array(longitude,mask=lonMask,fill_value=geoFillValue)
        longitude = longitude.filled()

        # Parse the geolocation asc file to get struct information which will be 
        # written to the ancillary asc files

        geoAscFileName = path.join(inDir,URID+".asc")
        LOG.debug("\nOpening %s..." % (geoAscFileName))

        geoAscFile = open(geoAscFileName,'rt')

        #RangeDateTimeStr =  _getAscLine(geoAscFile,"RangeDateTime")
        RangeDateTimeStr =  _getAscLine(geoAscFile,"ObservedDateTime")
        RangeDateTimeStr =  string.replace(RangeDateTimeStr,"ObservedDateTime","RangeDateTime")
        GRingLatitudeStr =  _getAscStructs(geoAscFile,"GRingLatitude",12)
        GRingLongitudeStr =  _getAscStructs(geoAscFile,"GRingLongitude",12)

        geoAscFile.close()

        # Loop through the required NAAPS datasets and create the blobs.
        # FIXME : Handle pathological geolocation cases

        firstGranule = True

        for dSet in masksCollShortNames :
        
            LOG.debug("\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n")
            LOG.debug("Processing dataset %s for %s" % (NAAPS_shortNameToBlobName[dSet],dSet))

            # FIXME : Account for dateline and pole crossings...

            # Massage the NAAPS data array a bit...
            NAAPS_anc = np.array(NAAPS_globalGridData[dSet])[::-1,:]
            NAAPS_anc = np.roll(NAAPS_anc,360)

            if (firstGranule) :

                LOG.debug("\nGranulating %s ..." % (dSet))
                LOG.debug("latitide,longitude shapes: %s, %s" %(str(latitude.shape) , str(longitude.shape)))
                LOG.debug("NAAPS_anc.shape = %s" % (str(NAAPS_anc.shape)))
                LOG.debug("gridLat.shape = %s" % (str(gridLat.shape)))
                LOG.debug("gridLon.shape = %s" % (str(gridLon.shape)))

                LOG.debug("min of NAAPS_anc  = "%(np.min(NAAPS_anc)))
                LOG.debug("max of NAAPS_anc  = "%(np.max(NAAPS_anc)))

                data,dataIdx = _grid2Gran(np.ravel(latitude),
                                          np.ravel(longitude),
                                          NAAPS_anc.astype(np.float64),
                                          gridLat,
                                          gridLon)

                data = data.reshape(latitude.shape)
                dataIdx = dataIdx.reshape(latitude.shape)
                firstGranule = False
                LOG.debug("Shape of first granulated %s data is %s" % (dSet,np.shape(data)))
                LOG.debug("Shape of first granulated %s dataIdx is %s" % (dSet,np.shape(dataIdx)))

            else :

                LOG.debug("Granulating %s using existing data indices." % (dSet))
                NAAPS_anc = np.ravel(NAAPS_anc)
                data = np.ravel(NAAPS_anc)[np.ravel(dataIdx)]
                data = data.reshape(latitude.shape)
                LOG.debug("Shape of subsequent granulated %s is %s" % (dSet,np.shape(data)))

            # Fill the required pixel trim rows in the granulated NAAPS data with 
            # the ONBOARD_PT_FILL value for the correct data type

            fillValue = trimObj.sdrTypeFill['ONBOARD_PT_FILL'][data.dtype.name]        
            data = ma.array(data,mask=modTrimMask,fill_value=fillValue)
            data = data.filled()

            # Create new NAAPS ancillary blob, and copy granulated data to it

            endian = ancEndian
            xmlName = path.join(ADL_HOME,'xml/VIIRS',NAAPS_shortNameToXmlName[dSet])

            # Create a new URID to be used in making the asc filenames

            URID_dict = _getURID()

            URID = URID_dict['URID']
            creationDate_nousecStr = URID_dict['creationDate_nousecStr']
            creationDateStr = URID_dict['creationDateStr']

            # Create a new directory in the input directory for the new ancillary
            # asc and blob files

            blobDir = inDir

            ascFileName = path.join(blobDir,URID+'.asc')
            blobName = path.join(blobDir,URID+'.'+dSet)

            LOG.debug("ascFileName : %s" % (ascFileName))
            LOG.debug("blobName : %s" % (blobName))

            # Create a new ancillary blob, and copy the data to it.
            newNAAPSblobObj = adl_blob.create(xmlName, blobName, endian=endian, overwrite=True)
            newNAAPSblobArrObj = newNAAPSblobObj.as_arrays()

            blobData = getattr(newNAAPSblobArrObj,'faot550')
            blobData[:,:] = data[:,:]

            # Make a new NAAPS asc file from the template, and substitute for the various tags

            ascTemplateFileName = path.join(ADL_ASC_TEMPLATES,"VIIRS-ANC_Template.asc")

            LOG.debug("Creating new asc file\n%s\nfrom template\n%s" % (ascFileName,ascTemplateFileName))
            
            ANC_fileList = gridBlobFiles
            for idx in range(len(ANC_fileList)) :
                ANC_fileList[idx] = path.basename(ANC_fileList[idx])
            ANC_fileList.sort()
            ancGroupRecipe = '    ("N_Anc_Filename" STRING EQ "%s")'
            ancFileStr = "%s" % ("\n").join([ancGroupRecipe % (str(files)) for files in ANC_fileList])

            LOG.debug("RangeDateTimeStr = %s\n" % (RangeDateTimeStr))
            LOG.debug("GRingLatitudeStr = \n%s\n" % (GRingLatitudeStr))
            LOG.debug("GRingLongitudeStr = \n%s\n" % (GRingLongitudeStr))


            ascTemplateFile = open(ascTemplateFileName,"rt") # Open template file for reading
            ascFile = open(ascFileName,"wt") # create a new text file

            LOG.debug("Template file %s is %r with mode %s" %(ascTemplateFileName,'not open' if ascTemplateFile.closed else 'open',ascTemplateFile.mode))

            LOG.debug("New file %s is %r with mode %s" %(ascFileName,'not open' if ascFile.closed else 'open',ascFile.mode))

            for line in ascTemplateFile.readlines():
               line = line.replace("CSPP_URID",URID)
               line = line.replace("CSPP_CREATIONDATETIME_NOUSEC",creationDate_nousecStr)
               line = line.replace("CSPP_ANC_BLOB_FULLPATH",blobName)
               line = line.replace("CSPP_ANC_COLLECTION_SHORT_NAME",dSet)
               line = line.replace("CSPP_GRANULE_ID",N_Granule_ID)
               line = line.replace("CSPP_CREATIONDATETIME",creationDateStr)
               line = line.replace("  CSPP_RANGE_DATE_TIME",RangeDateTimeStr)
               line = line.replace("  CSPP_GRINGLATITUDE",GRingLatitudeStr)
               line = line.replace("  CSPP_GRINGLONGITUDE",GRingLongitudeStr)
               line = line.replace("    CSPP_ANC_SOURCE_FILES",ancFileStr)
               ascFile.write(line) 

            ascFile.close()
            ascTemplateFile.close()


def _retrieve_NISE_files_and_convert(geoDicts):
    ''' Download the NISE Snow/Ice files which cover the dates of the geolocation files.'''

    CSPP_RT_HOME = os.getenv('CSPP_RT_HOME')
    ANC_SCRIPTS_PATH = path.join(CSPP_RT_HOME,'viirs')
    CSPP_RT_ANC_CACHE_DIR = os.getenv('CSPP_RT_ANC_CACHE_DIR')

    niseFiles = []

    for geoDict in geoDicts:
        timeObj = geoDict['ObservedStartTime']
        dateStamp = timeObj.strftime("%Y%j")

        try :
            LOG.info('Retrieving NISE files for %s ...' % (dateStamp))
            cmdStr = '%s/get_anc_cspp_nise.csh %s' % (ANC_SCRIPTS_PATH,dateStamp)
            LOG.info('\t%s' % (cmdStr))
            args = shlex.split(cmdStr)
            LOG.debug('\t%s' % (repr(args)))

            procRetVal = 0
            procObj = subprocess.Popen(args,bufsize=0, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            procObj.wait()
            procRetVal = procObj.returncode

            procOutput = procObj.stdout.readlines()
            #procOutput = procObj.stdout.read()

            # FIXME : How to get this output to have linebreaks when using readlines()
            #LOG.debug(procOutput)
            
            for lines in procOutput:
                if CSPP_RT_ANC_CACHE_DIR in lines :
                    lines = string.replace(lines,'\n','')
                    niseFiles.append(lines)

            # TODO : On error, jump to a cleanup routine
            if not (procRetVal == 0) :
                LOG.error('Retrieval of NISE files failed for %s.' % (dateStamp))
                #sys.exit(procRetVal)

        except Exception, err:
            LOG.warn( "%s" % (str(err)))

    # Uniqify the list of NISE files
    niseFiles = list(set(niseFiles))
    niseFiles.sort()

    # Convert the NISE HDFEOS file to HDF5
    newNiseFiles = []
    for files in niseFiles:
        niseDir = path.dirname(files)
        niseName = path.basename(files)
        new_niseName = string.replace(niseName,'HDFEOS','h5')
        new_niseName = path.join(niseDir,new_niseName)

        if not os.path.exists(new_niseName):

            try :
                LOG.info('Converting HDFEOS NISE file\n\t%s\nto HDF5...\n\t%s' % (files,new_niseName))
                cmdStr = '%s/common/COTS/hdf5/hdf5-1.8.4/bin/h4toh5 %s %s' % (CSPP_RT_HOME,files,new_niseName)
                LOG.info('\t%s' % (cmdStr))
                args = shlex.split(cmdStr)
                LOG.debug('\t%s' % (repr(args)))

                procRetVal = 0
                procObj = subprocess.Popen(args,bufsize=0, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                procObj.wait()
                procRetVal = procObj.returncode
                procOutput = procObj.stdout.read()
                LOG.debug(procOutput)

                newNiseFiles.append(new_niseName)

            except Exception, err:
                LOG.warn( "%s" % (str(err)))
        else :
            LOG.info('NISE file %s exists, skipping.' % (new_niseName))
            newNiseFiles.append(new_niseName)

    return newNiseFiles


def _retrieve_NISE_files(geoDicts):
    ''' Download the NISE Snow/Ice files which cover the dates of the geolocation files.'''

    CSPP_RT_HOME = os.getenv('CSPP_RT_HOME')
    ANC_SCRIPTS_PATH = path.join(CSPP_RT_HOME,'viirs')
    CSPP_RT_ANC_CACHE_DIR = os.getenv('CSPP_RT_ANC_CACHE_DIR')

    niseFiles = []

    for geoDict in geoDicts:
        timeObj = geoDict['ObservedStartTime']
        dateStamp = timeObj.strftime("%Y%j")

        try :
            LOG.info('Retrieving NISE files for %s ...' % (dateStamp))
            cmdStr = '%s/get_anc_cspp_nise.csh %s' % (ANC_SCRIPTS_PATH,dateStamp)
            LOG.info('\t%s' % (cmdStr))
            args = shlex.split(cmdStr)
            LOG.debug('\t%s' % (repr(args)))

            procRetVal = 0
            procObj = subprocess.Popen(args,bufsize=0, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            procObj.wait()
            procRetVal = procObj.returncode

            procOutput = procObj.stdout.readlines()
            #procOutput = procObj.stdout.read()

            # FIXME : How to get this output to have linebreaks when using readlines()
            #LOG.debug(procOutput)
            
            for lines in procOutput:
                if CSPP_RT_ANC_CACHE_DIR in lines :
                    lines = string.replace(lines,'\n','')
                    niseFiles.append(lines)

            # TODO : On error, jump to a cleanup routine
            if not (procRetVal == 0) :
                LOG.error('Retrieval of NISE files failed for %s.' % (dateStamp))
                #sys.exit(procRetVal)

        except Exception, err:
            LOG.warn( "%s" % (str(err)))

    # Uniqify the list of NISE files
    niseFiles = list(set(niseFiles))
    niseFiles.sort()

    return niseFiles


# Look through new log files for completed messages
# based on CrIS SDR original
def check_log_files(work_dir, pid, xml):
    """
    Find the log file
    Look for success
    Return True if found
    Display log message and hint if problem occurred
    """
    # retrieve exe name and log path from lw file
    logpat = os.path.join(work_dir, "log", "*%d*.log" % pid)

    n_err=0
    for logFile in glob(logpat):
        LOG.debug( "Checking Log file " +logFile +" for errors.")
        n_err += adl_log.scan_log_file(COMMON_LOG_CHECK_TABLE, logFile)

    if n_err == 0 :
        LOG.debug(" Log: "+ logFile + " "+ xml + " Completed successfully" )
        return True
    else :
        LOG.error( " Log: "+ logFile +" " + xml + " Completed unsuccessfully" )
        return False


def run_xml_files(work_dir, xml_files_to_process, setup_only=False, **additional_env):
    """Run each VIIRS EDR MASKS XML input in sequence
    return the list of granule IDs which crashed, and list of granule ids which did not create output
    """
    crashed_runs = set()
    no_output_runs = set()
    geo_problem_runs = set()
    bad_log_runs = set()
    first = True

    # obtain pre-existing granule list
    modGeoTCPattern = os.path.join(work_dir, 'GMTCO*.h5')
    cmaskPattern = os.path.join(work_dir, 'IICMO*.h5')
    activeFiresPattern = os.path.join(work_dir, 'AVAFO*.h5')
    # TODO : Enable for VIIRS AOT
    #aotIpPattern = os.path.join(work_dir, 'IVAOT*.h5')
    #aotEdrPattern = os.path.join(work_dir, 'VAOOO*.h5')
    #suspMatEdrPattern = os.path.join(work_dir, 'VSUMO*.h5')

    # prior_granules dicts contain (N_GranuleID,HDF5File) key,value pairs.
    # *ID dicts contain (HDF5File,N_GranuleID) key,value pairs.
    
    # Get the N_GranuleID and filename of existing cmask files in the work_dir ? ...
    cmask_prior_granules, cmask_ID = h5_xdr_inventory(cmaskPattern, CM_GRANULE_ID_ATTR_PATH)
    cmask_prior_granules = set(cmask_prior_granules.keys())

    # Get the (N_GranuleID,hdfFileName) pairs for the existing Active fires files
    activeFires_prior_granules, afires_ID = h5_xdr_inventory(activeFiresPattern, AF_GRANULE_ID_ATTR_PATH)
    activeFires_prior_granules = set(activeFires_prior_granules.keys())

    # Get the (N_GranuleID,hdfFileName) pairs for the existing Aerosol IP files
    # TODO : Enable for VIIRS AOT
    #aerosolIP_prior_granules, aotIp_ID = h5_xdr_inventory(aotIpPattern, AOT_IP_GRANULE_ID_ATTR_PATH)
    #aerosolIP_prior_granules = set(aerosolIP_prior_granules.keys())

    # Get the (N_GranuleID,hdfFileName) pairs for the existing Aerosol EDR files
    # TODO : Enable for VIIRS AOT
    #aerosolEDR_prior_granules, aotEdr_ID = h5_xdr_inventory(aotEdrPattern, AOT_EDR_GRANULE_ID_ATTR_PATH)
    #aerosolEDR_prior_granules = set(aerosolEDR_prior_granules.keys())

    # Get the (N_GranuleID,hdfFileName) pairs for the existing Suspended Matter EDR files
    # TODO : Enable for VIIRS AOT
    #suspMatEDR_prior_granules, suspMatEdr_ID = h5_xdr_inventory(suspMatEdrPattern, SUSMAT_EDR_GRANULE_ID_ATTR_PATH)
    #suspMatEDR_prior_granules = set(suspMatEDR_prior_granules.keys())

    for granule_id, xml in xml_files_to_process:

        t1 = time()
        
        cmd = [ADL_VIIRS_MASKS_EDR, xml]
        # TODO : Enable for VIIRS AOT
        #cmd = [ADL_VIIRS_AEROSOL_EDR, xml]
        
        if setup_only:
            print ' '.join(cmd)
        else:
            LOG.debug('executing "%s"' % ' '.join(cmd))
            LOG.debug('additional environment variables: %s' % repr(additional_env))
            try:
                pid = sh(cmd, env=env(**additional_env), cwd=work_dir)
                LOG.debug("%r ran as pid %d" % (cmd, pid))
                if not check_log_files(work_dir, pid, xml):
                    bad_log_runs.add(granule_id)

            except CalledProcessError as oops:
                LOG.debug(traceback.format_exc())
                LOG.error('ProEdrViirsMasksController.exe failed on %r: %r. Continuing...' % (xml, oops))
                crashed_runs.add(granule_id)
            first = False

            # check new IICMO output granules
            cmask_new_granules, cmask_ID = h5_xdr_inventory(cmaskPattern, CM_GRANULE_ID_ATTR_PATH, state=cmask_ID)
            LOG.debug('new IICMO granules after this run: %s' % (repr(cmask_new_granules)))
            if granule_id not in cmask_new_granules:
                LOG.warning('no IICMO HDF5 output for %s' % (granule_id))
                no_output_runs.add(granule_id)
            else:
                filename = cmask_new_granules[granule_id]

            # check new AVAFO output granules
            afires_new_granules, afires_ID = h5_xdr_inventory(activeFiresPattern, AF_GRANULE_ID_ATTR_PATH, state=afires_ID)
            LOG.debug('new AVAFO granules after this run: %s' % (repr(afires_new_granules)))
            if granule_id not in afires_new_granules:
                LOG.warning('no AVAFO HDF5 output for %s' % (granule_id))
                no_output_runs.add(granule_id)
            else:
                filename = afires_new_granules[granule_id]

            # TODO : Enable for VIIRS AOT
            # check new IVAOT output granules
            #aotIp_new_granules, aotIp_ID = h5_xdr_inventory(aotIpPattern, AOT_IP_GRANULE_ID_ATTR_PATH, state=aotIp_ID)
            #LOG.debug('new IVAOT granules after this run: %s' % repr(aotIp_new_granules))
            #if granule_id not in aotIp_new_granules:
                #LOG.warning('no IVAOT HDF5 output for %s' % granule_id)
                #no_output_runs.add(granule_id)
            #else:
                #filename = aotIp_new_granules[granule_id]

            # TODO : Enable for VIIRS AOT
            # check new VAOOO output granules
            #aotEdr_new_granules, aotEdr_ID = h5_xdr_inventory(aotEdrPattern, AOT_EDR_GRANULE_ID_ATTR_PATH, state=aotEdr_ID)
            #LOG.debug('new VAOOO granules after this run: %s' % repr(aotEdr_new_granules))
            #if granule_id not in aotEdr_new_granules:
                #LOG.warning('no VAOOO HDF5 output for %s' % granule_id)
                #no_output_runs.add(granule_id)
            #else:
                #filename = aotEdr_new_granules[granule_id]

            # TODO : Enable for VIIRS AOT
            # check new VSUMO output granules
            #suspMatEdr_new_granules, suspMatEdr_ID = h5_xdr_inventory(suspMatEdrPattern, SUSMAT_EDR_GRANULE_ID_ATTR_PATH, state=suspMatEdr_ID)
            #LOG.debug('new VSUMO granules after this run: %s' % repr(suspMatEdr_new_granules))
            #if granule_id not in suspMatEdr_new_granules:
                #LOG.warning('no VSUMO HDF5 output for %s' % granule_id)
                #no_output_runs.add(granule_id)
            #else:
                #filename = suspMatEdr_new_granules[granule_id]

        t2 = time()
        LOG.info ( "Controller ran in %f seconds." % (t2-t1))

    LOG.debug("cmask_ID.values() = \n%r" % (cmask_ID.values()))
    LOG.debug("set(cmask_ID.values()) = \n%r" % (set(cmask_ID.values())))
    LOG.debug("cmask_prior_granules = \n%r" % (cmask_prior_granules))
    cmask_granules_made = set(cmask_ID.values()) - cmask_prior_granules

    LOG.debug("afires_ID.values() = \n%r" % (afires_ID.values()))
    LOG.debug("set(afires_ID.values()) = \n%r" % (set(afires_ID.values())))
    LOG.debug("activeFires_prior_granules = \n%r" % (activeFires_prior_granules))
    activeFires_granules_made = set(afires_ID.values()) - activeFires_prior_granules

    # TODO : Enable for VIIRS AOT
    #LOG.debug("aotIp_ID.values() = \n%r" % (aotIp_ID.values()))
    #LOG.debug("set(aotIp_ID.values()) = \n%r" % (set(aotIp_ID.values())))
    #LOG.debug("aerosolIP_prior_granules = \n%r" % (aerosolIP_prior_granules))
    #aerosolIP_granules_made = set(aotIp_ID.values()) - aerosolIP_prior_granules

    # TODO : Enable for VIIRS AOT
    #LOG.debug("aotEdr_ID.values() = \n%r" % (aotEdr_ID.values()))
    #LOG.debug("set(aotEdr_ID.values()) = \n%r" % (set(aotEdr_ID.values())))
    #LOG.debug("aerosolEDR_prior_granules = \n%r" % (aerosolEDR_prior_granules))
    #aerosolEDR_granules_made = set(aotEdr_ID.values()) - aerosolEDR_prior_granules

    # TODO : Enable for VIIRS AOT
    #LOG.debug("suspMatEdr_ID.values() = \n%r" % (suspMatEdr_ID.values()))
    #LOG.debug("set(suspMatEdr_ID.values()) = \n%r" % (set(suspMatEdr_ID.values())))
    #LOG.debug("suspMatEDR_prior_granules = \n%r" % (suspMatEDR_prior_granules))
    #suspMatEDR_granules_made = set(suspMatEdr_ID.values()) - suspMatEDR_prior_granules


    LOG.info('cmask granules created: %s' %( ', '.join(list(cmask_granules_made))))
    LOG.info('activeFires granules created: %s' % (', '.join(list(activeFires_granules_made))))
    #LOG.info('aerosolIP granules created: %s' % ', '.join(list(aerosolIP_granules_made)))
    #LOG.info('aerosolEDR granules created: %s' % ', '.join(list(aerosolEDR_granules_made)))
    #LOG.info('suspMatEDR granules created: %s' % ', '.join(list(suspMatEDR_granules_made)))

    if no_output_runs:
        LOG.info('granules that failed to generate output: %s' % (', '.join(no_output_runs)))
    if geo_problem_runs:
        LOG.warning('granules which had no N_Geo_Ref: %s' % ', '.join(geo_problem_runs))
    if crashed_runs:
        LOG.warning('granules that crashed ADL: %s' % (', '.join(crashed_runs)))
    if bad_log_runs:
        LOG.warning('granules that produced logs indicating problems: %s' % (', '.join(bad_log_runs)))
    if not cmask_granules_made:
        LOG.warning('no HDF5 SDRs were created')

    return crashed_runs, no_output_runs, geo_problem_runs, bad_log_runs


def _cleanup(work_dir, xml_glob, log_dir_glob, *more_dirs):
    """upon successful run, clean out work directory"""

    LOG.info("cleaning up work directory...")
    for fn in glob(os.path.join(work_dir, '????????-?????-????????-????????.*')):
        LOG.debug('removing %s' % (fn))
        os.unlink(fn)

    LOG.info("Removing task xml files...")
    for fn in glob(os.path.join(work_dir, xml_glob)):
        LOG.debug('removing task file %s' % (fn))
        os.unlink(fn)

    LOG.info("Removing log directories %s ..."%(log_dir_glob))
    for dirname in glob(os.path.join(work_dir,log_dir_glob)):
        LOG.debug('removing logs in %s' % (dirname))
        try :
            rmtree(dirname, ignore_errors=False)
        except Exception, err:
            LOG.warn( "%s" % (str(err)))

    # FIXME : Cannot seem to remove the log directory, complains that the dir is not empty.
    LOG.info("Removing other directories ...")
    for dirname in more_dirs:
        fullDirName = os.path.join(work_dir,dirname)
        LOG.debug('removing %s' % (fullDirName))
        try :
            rmtree(fullDirName, ignore_errors=False)
        except Exception, err:
            LOG.warn( "%s" % (str(err)))


def main():

    endianChoices = ['little','big']

    description = '''Run the ADL VIIRS EDR Masks Controller.'''
    usage = "usage: %prog [mandatory args] [options]"
    version = "%prog "+__version__

    parser = optparse.OptionParser(description=description,usage=usage,version=version)

    # Mandatory arguments

    mandatoryGroup = optparse.OptionGroup(parser, "Mandatory Arguments",
                        "At a minimum these arguments must be specified")

    mandatoryGroup.add_option('-i','--input_files',
                      action="store",
                      dest="inputFiles",
                      type="string",
                      help="The fully qualified path to the input files. May be a directory or a file glob.")

    parser.add_option_group(mandatoryGroup)

    # Optional arguments 

    optionalGroup = optparse.OptionGroup(parser, "Extra Options",
                         "These options may be used to customize behaviour of this program.")

    optionalGroup.add_option('-w','--work_directory',
                      action="store",
                      dest="work_dir",
                      type="string",
                      default='.',
                      help="The directory which all activity will occur in, defaults to the current directory.")
    
    optionalGroup.add_option('--skip_sdr_unpack',
                      action="store_true",
                      dest="skipSdrUnpack",
                      help="Skip the unpacking of the VIIRS SDR HDF5 files.")

    optionalGroup.add_option('--skip_aux_linking',
                      action="store_true",
                      dest="skipAuxLinking",
                      help="Skip the the linking to auxillary files.")

    optionalGroup.add_option('--skip_ancillary',
                      action="store_true",
                      dest="skipAncillary",
                      help="Skip the retrieval and granulation of ancillary data.")

    optionalGroup.add_option('--skip_algorithm',
                      action="store_true",
                      dest="skipAlgorithm",
                      help="Skip running the VIIRS Masks algorithm.")

    optionalGroup.add_option('--debug',
                      action="store_true",
                      dest="cspp_debug",
                      default=False,
                      help="Enable debug mode on ADL and avoid cleaning workspace")

    optionalGroup.add_option('--sdr_endianness',
                      action="store",
                      dest="sdr_Endianness",
                      type="choice",
                      default='little',
                      choices=endianChoices,
                      help='''The input VIIRS SDR endianness.\n\n
                              Possible values are...
                              %s. [default: 'little']
                           ''' % (endianChoices.__str__()[1:-1]))
    
    optionalGroup.add_option('--anc_endianness',
                      action="store",
                      dest="anc_Endianness",
                      type="choice",
                      default='little',
                      choices=endianChoices,
                      help='''The input VIIRS ancillary endianness.\n\n
                              Possible values are...
                              %s. [default: 'little']
                           ''' % (endianChoices.__str__()[1:-1]))
    
    optionalGroup.add_option('-v', '--verbose',
                      dest='verbosity',
                      action="count",
                      default=0,
                      help='each occurrence increases verbosity 1 level from ERROR: -v=WARNING -vv=INFO -vvv=DEBUG')

    parser.add_option_group(optionalGroup)

    # Parse the arguments from the command line

    (options, args) = parser.parse_args()

    # Check that all of the mandatory options are given. If one or more 
    # are missing, print error message and exit...

    mandatories = ['inputFiles']
    mand_errors = ["Missing mandatory argument [-i input_files --input_files=input_files]"]
    isMissingMand = False
    for m,m_err in zip(mandatories,mand_errors):
        if not options.__dict__[m]:
            isMissingMand = True
            print m_err
    if isMissingMand :
        parser.error("Incomplete mandatory arguments, aborting...")

    # Set up the logging

    levels = [logging.ERROR, logging.WARN, logging.INFO, logging.DEBUG]
    configure_logging(level = levels[min(options.verbosity,3)])

    # Test the logging

    #test_logging()
    #sys.exit(0)

    # Determine the correct input file path and glob

    input_dir,inputGlobs = _create_input_file_globs(options.inputFiles)
    #########################################################################################
    # Determine the correct input file path and glob

  

    if inputGlobs['GEO'] is None or inputGlobs['MOD'] is None or inputGlobs['IMG'] is None :
        LOG.error("No input files found, aborting...")
        sys.exit(1)
    #########################################################################################
    # Set the work directory

    work_dir = os.path.abspath(options.work_dir)
    LOG.debug('Setting the work directory to %r' % (work_dir))

    # Check the environment variables, and whether we can write to the working directory
    check_env(work_dir)
    
    # create work directory
    if not os.path.isdir(work_dir):
        LOG.info('creating directory %s' % (work_dir))
        os.makedirs(work_dir)
    log_dir = os.path.join(work_dir, 'log')
    if not os.path.isdir(log_dir):
        LOG.debug('creating directory %s' % (log_dir))
        os.makedirs(log_dir)
    perf_dir = os.path.join(work_dir, 'perf')
    if not os.path.isdir(perf_dir):
        LOG.debug('creating directory %s' % (perf_dir))
        os.makedirs(perf_dir)
    anc_dir = os.path.join(work_dir, 'linked_data')
    if not os.path.isdir(anc_dir):
        LOG.debug('creating directory %s' % (anc_dir))
        os.makedirs(anc_dir)

    # Unpack HDF5 VIIRS SDRs in the input directory to the work directory
    unpacking_problems = 0
    if not options.skipSdrUnpack :
        h5_names = glob(path.join(input_dir,inputGlobs['GEO'])) \
                 + glob(path.join(input_dir,inputGlobs['MOD'])) \
                 + glob(path.join(input_dir,inputGlobs['IMG']))

        t1 = time()
        for fn in h5_names:
            try:
                LOG.info('Unpacking %r ...' % (fn))
                unpack(work_dir, fn)
            except CalledProcessError as oops:
                LOG.debug(traceback.format_exc())
                LOG.error('ADL_Unpacker failed on %r: %r . Continuing' % (fn, oops))
                unpacking_problems += 1
        t2 = time()
        LOG.info ("Unpacking of VIIRS SDR files took %f seconds." % (t2-t1))
    else :
        LOG.info('Skipping SDR unpacking, assuming all VIIRS SDR blob and asc files are present.')

    LOG.debug("Unpacking problems = %d" % (unpacking_problems))

    # Read through ascii metadata and build up information table
    LOG.info('Sifting through metadata to find VIIRS SDR processing candidates')
    geolocationShortNames = ['VIIRS-MOD-RGEO-TC','VIIRS-MOD-RGEO','VIIRS-MOD-GEO-TC','VIIRS-MOD-GEO']
    for geoType in geolocationShortNames :
        LOG.info("Searching for VIIRS geolocation %s..." % (geoType))
        anc_granules_to_process = sorted(list(sift_metadata_for_viirs_sdr(geoType,crossGran=None,work_dir=work_dir)))
        granules_to_process     = sorted(list(sift_metadata_for_viirs_sdr(geoType,crossGran=1,   work_dir=work_dir)))
        
        if granules_to_process :
            LOG.debug("\tgranules_to_process() has %d objects..."%(len(granules_to_process)))
            LOG.info(', '.join(x['N_Granule_ID'] for x in granules_to_process))
            LOG.debug("\tanc_granules_to_process() has %d objects..."%(len(anc_granules_to_process)))
            LOG.info(', '.join(x['N_Granule_ID'] for x in anc_granules_to_process))
            break
        else :
            LOG.warn("\tNo granules for VIIRS geolocation %s..." % (geoType))

    if not granules_to_process:
        LOG.error("Error: Found no granules to process!")
        num_xml_files_to_process = 0
        num_no_output_runs = 0
        noncritical_problem = False
        environment_error = False
        return get_return_code(unpacking_problems, num_xml_files_to_process, num_no_output_runs, noncritical_problem, environment_error)

    # Set the VIIRS SDR endianness from the input option...

    global sdrEndian
    set_sdr_endian(options.sdr_Endianness)

    # Set the VIIRS ancillary endianness from the input option...

    global ancEndian
    set_anc_endian(options.anc_Endianness)

    # Get some information about the geolocation files

    LOG.debug("\nGetting geolocation information...")

    LOG.debug("\n%13s%28s%29s" % ('N_Granule_ID','ObservedStartTime','ObservedEndTime'))
    for dicts in granules_to_process :
        LOG.debug("%15s%30s%30s"%(dicts['N_Granule_ID'],dicts['ObservedStartTime'],dicts['ObservedEndTime']))

    # Expand any user specifiers in the various paths

    #CSPP_RT_ANC_HOME = os.getenv('CSPP_RT_ANC_HOME')
    CSPP_RT_ANC_PATH = os.getenv('CSPP_RT_ANC_PATH')
    CSPP_RT_ANC_CACHE_DIR = os.getenv('CSPP_RT_ANC_CACHE_DIR')
    CSPP_RT_ANC_TILE_PATH = os.getenv('CSPP_RT_ANC_TILE_PATH')
    LD_LIBRARY_PATH = os.getenv('LD_LIBRARY_PATH')
    DSTATICDATA = os.getenv('DSTATICDATA')
    LOG.debug("\nCSPP_RT_ANC_HOME:     %s" % (CSPP_RT_ANC_HOME))
    LOG.debug("CSPP_RT_ANC_PATH:       %s" % (CSPP_RT_ANC_PATH))
    LOG.debug("CSPP_RT_ANC_CACHE_DIR:  %s" % (CSPP_RT_ANC_CACHE_DIR))
    LOG.debug("CSPP_RT_ANC_TILE_PATH:  %s" % (CSPP_RT_ANC_TILE_PATH))
    LOG.debug("LD_LIBRARY_PATH:     %s" % (LD_LIBRARY_PATH))
    LOG.debug("DSTATICDATA:         %s" % (DSTATICDATA))

    # Link in auxillary files

    if not options.skipAuxLinking :
        LOG.info('Linking in the VIIRS EDR auxillary files...')
        _setupAuxillaryFiles(work_dir)
    else :
        LOG.info('Skipping linking in the VIIRS EDR auxillary files.')

    # Retrieve and granulate the required ancillary data...

    if not options.skipAncillary :

        t1 = time()

        LOG.info('Retrieving and granulating ancillary data...')

        # TODO : Clarify what this "allow_cache_update" is for?
        allow_cache_update = True
        # Search for dynamic ancillary data (GRIB and NISE) for the granules we intend to process.
        # - If it's not in the cache already, download it.
        # - Softlink those into the workspace
        # FUTURE: This currently doesn't operate correctly without this, since we rely on the cache scripts
        #   to tell us where to find the ancillary data whether it was downloaded or not.
        if allow_cache_update:
            LOG.info("Downloading GRIB and NISE ancillary into cache")
            gribFiles = _retrieve_grib_files(anc_granules_to_process)
            LOG.debug('dynamic ancillary GRIB files: %s' % repr(gribFiles))
            if (gribFiles == []) :
                LOG.error('Failed to find or retrieve any GRIB files, aborting...')
                return -1

            niseFiles = _retrieve_NISE_files(anc_granules_to_process)
            LOG.debug('dynamic ancillary NISE files: %s' % repr(niseFiles))
            if (niseFiles == []) :
                LOG.error('Failed to find or retrieve any NISE files, aborting...')
                return -1

            # TODO : Obtain a working NAAPS retrieval script from Kathy Strabala
            # TODO : Enable for VIIRS AOT
            #naapsFiles = _retrieve_NAAPS_files(anc_granules_to_process)
            #LOG.debug('dynamic ancillary NAAPS files: %s' % repr(naapsFiles))
            #if (naapsFiles == []) :
                #LOG.error('Failed to find or retrieve any NAAPS files, aborting...')
                #return -1

            all_dyn_anc = list(gribFiles) + list(niseFiles) # + list(naapsFiles)
            LOG.debug('dynamic ancillary files: %s' % repr(all_dyn_anc))


        # Setup GRC files

        _getGRC(work_dir,anc_granules_to_process)

        # Transcode the NCEP GRIB files into NCEP global grid blob files

        #gridBlobFiles = _create_NCEP_gridBlobs(gribFiles)
        gridBlobFiles = _create_NCEP_gridBlobs_alt(gribFiles)

        LOG.debug("gridBlobFiles: %r" % (gridBlobFiles))

        # Granulate the global grid NCEP blob files

        granBlobFiles = _granulate_NCEP_gridBlobs(work_dir,anc_granules_to_process,gridBlobFiles)
        #granBlobFiles = _granulate_NCEP_gridBlobs_alt(work_dir,anc_granules_to_process,gridBlobFiles)

        LOG.debug("granBlobFiles: %r" % (granBlobFiles))

        # Transcode the NAAPS GRIB files into NAAPS global grid blob files
        #gridBlobFiles = _create_NAAPS_gridBlobs(naapsFiles)
        #gridBlobFiles = glob(path.join(work_dir,'*.NAAPS-ANC-Int'))
        #LOG.debug("gridBlobFiles: %r" % (gridBlobFiles)

        # Granulate the global grid NAAPS blob files
        #granBlobFiles = _granulate_NAAPS_gridBlobs(work_dir,anc_granules_to_process,gridBlobFiles)
        #LOG.debug("granBlobFiles: %r" % (granBlobFiles)

        # Granulate the global DEM file

        DEM_granules = _granulate_DEM(anc_granules_to_process,work_dir)

        # Granulate the global IGBP file

        IGBP_granules = _granulate_IGBP(anc_granules_to_process,work_dir)

        # Combine the IGBP and the DEM to make the QSTLWM

        _QSTLWM(DEM_granules,IGBP_granules,anc_granules_to_process,work_dir)

        # Granulate the global NISE HDF5 files

        granNiseFiles = _granulate_NISE_list(work_dir,anc_granules_to_process,DEM_granules,niseFiles)

        # Granulate the global NDVI files

        _granulate_NDVI(work_dir,anc_granules_to_process)

        t2 = time()
        LOG.info ( "Generation of VIIRS Masks ancillary data took %f seconds." % (t2-t1))

    else :

        LOG.info('Skipping retrieval and granulation of ancillary data.')


    # Specify the algorithm we want to run via the Lw XML file.

    if not options.skipAlgorithm :

        # build XML configuration files for jobs that can be run
        LOG.debug("Building XML files for %d granules" % len(granules_to_process))
        
        # Generate the VIIRS Cloud Mask/Active Fires LW xml files for each N_Granule_ID...
        xml_files_to_process = generate_viirs_masks_edr_xml(work_dir, granules_to_process)

        # TODO : Enable for VIIRS AOT
        # Generate the VIIRS Aerosol Optical Thickness LW xml files for each N_Granule_ID...
        #xml_files_to_process = generate_viirs_aerosol_edr_xml(work_dir, granules_to_process)

        LOG.info('%d granules to process: \n%s' % (len(xml_files_to_process), ''.join(name+' -> '+xmlfile+'\n' for (name,xmlfile) in xml_files_to_process)))

        LOG.info("Running VIIRS EDR Masks on XML files")
        crashed_runs, no_output_runs, geo_problem_runs, bad_log_runs = run_xml_files(work_dir, xml_files_to_process, setup_only = False, WORK_DIR = work_dir, LINKED_ANCILLARY = work_dir, ADL_HOME=ADL_HOME)

        ## considered a noncritical problem if there were any crashed runs, runs that produced no output,
        ## runs where Geo failed, or runs where ADL logs indicated a problem
        noncritical_problem = crashed_runs or no_output_runs or geo_problem_runs or bad_log_runs

        ## print final disposition message and get return code
        environment_error=False
        rc = get_return_code(unpacking_problems, len(xml_files_to_process), len(no_output_runs), noncritical_problem, environment_error)

        ## if no errors or only non-critical errors: clean up
        if rc == 0 and not options.cspp_debug:
            LOG.debug("Cleaning up workspace...")
            _cleanup(work_dir, 'edr_viirs_masks*.xml', 'ProEdrViirsMasksController.exe_*', log_dir, anc_dir, perf_dir)

            # TODO : Enable for VIIRS AOT
            #_cleanup(work_dir, 'edr_viirs_aerosol*.xml', 'ProEdrViirsAerosolController.exe_*', log_dir, anc_dir, perf_dir)

        return rc
    
    else :

        LOG.info('Skipping execution of VIIRS Masks Controller.')


if __name__=='__main__':
    sys.exit(main())  
