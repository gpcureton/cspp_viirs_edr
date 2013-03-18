#!/usr/bin/env python
# encoding: utf-8
"""
adl_viirs_edr.py

Purpose: Run one or more of the VIIRS EDR Controllers using ADL.

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
ADL_VIIRS_MASKS_EDR=path.abspath(path.join(ADL_HOME, 'bin', 'ProEdrViirsMasksController.exe'))
ADL_VIIRS_AEROSOL_EDR=path.abspath(path.join(ADL_HOME, 'bin', 'ProEdrViirsAerosolController.exe'))
ADL_VIIRS_SST_EDR=path.abspath(path.join(ADL_HOME, 'bin', 'ProEdrViirsSstController.exe'))

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

CSPP_RT_ANC_HOME = path.abspath(os.getenv('CSPP_RT_ANC_HOME'))



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
#IGBP_dict = {
    #'IGBP_EVERNEEDLE'    : 1,
    #'IGBP_EVERBROAD'     : 2,
    #'IGBP_DECIDNEEDLE'   : 3,
    #'IGBP_DECIDBROAD'    : 4,
    #'IGBP_MIXEDFOREST'   : 5,
    #'IGBP_CLOSEDSHRUBS'  : 6,
    #'IGBP_OPENSHRUBS'    : 7,
    #'IGBP_WOODYSAVANNA'  : 8,
    #'IGBP_SAVANNA'       : 9,
    #'IGBP_GRASSLAND'     : 10,
    #'IGBP_WETLANDS'      : 11,
    #'IGBP_CROPLANDS'     : 12,
    #'IGBP_URBAN'         : 13,
    #'IGBP_CROPNATMOSAIC' : 14,
    #'IGBP_SNOWICE'       : 15,
    #'IGBP_BARREN'        : 16,
    #'IGBP_WATERBODIES'   : 17,
    #'IGBP_UNCLASSIFIED'  : 30,
    #'IGBP_FILL'          : 31,
    #'IGBP_MIN'           : 1,
    #'IGBP_MAX'           : 17
#}

# QST Land Water Mask type value enumerations to be used with the qstlwm field
#QSTLWM_list = [
    #'EVERGREEN_NEEDLELEAF_FOREST', 'EVERGREEN_BROADLEAF_FORESTS', 'DECIDUOUS_NEEDLELEAF_FORESTS', 'DECIDUOUS_BROADLEAF_FORESTS', 'MIXED_FORESTS', 'CLOSED_SHRUBLANDS', 'OPEN_SHRUBLANDS', 'WOODY_SAVANNAS', 'SAVANNAS', 'GRASSLANDS', 'PERMANENT_WETLANDS', 'CROPLANDS', 'URBAN_AND_BUILDUP_LANDS', 'CROPLAND_NATURAL_VEGETATION', 'SNOW_AND_ICE', 'BARREN_OR_SPARSELY_VEGETATED', 'OCEAN_SEA', 'INLAND_WATER', 'COASTAL_WATER', 'UNCLASSIFIED_LAND', 'FILL_VALUE'
#]
#QSTLWM_dict = {
    #'EVERGREEN_NEEDLELEAF_FOREST'  : 1,
    #'EVERGREEN_BROADLEAF_FORESTS'  : 2,
    #'DECIDUOUS_NEEDLELEAF_FORESTS' : 3,
    #'DECIDUOUS_BROADLEAF_FORESTS'  : 4,
    #'MIXED_FORESTS'                : 5,
    #'CLOSED_SHRUBLANDS'            : 6,
    #'OPEN_SHRUBLANDS'              : 7,
    #'WOODY_SAVANNAS'               : 8,
    #'SAVANNAS'                     : 9,
    #'GRASSLANDS'                   : 10,
    #'PERMANENT_WETLANDS'           : 11,
    #'CROPLANDS'                    : 12,
    #'URBAN_AND_BUILDUP_LANDS'      : 13, 
    #'CROPLAND_NATURAL_VEGETATION'  : 14,
    #'SNOW_AND_ICE'                 : 15,
    #'BARREN_OR_SPARSELY_VEGETATED' : 16,
    #'OCEAN_SEA'                    : 17,
    #'INLAND_WATER'                 : 18,
    #'COASTAL_WATER'                : 19,
    #'UNCLASSIFIED_LAND'            : 20,
    #'FILL_VALUE'                   : 255
#}

# Digital Elevation Model (DEM) land sea mask types
#DEM_list = ['DEM_SHALLOW_OCEAN','DEM_LAND','DEM_COASTLINE','DEM_SHALLOW_INLAND_WATER','DEM_EPHEMERAL_WATER','DEM_DEEP_INLAND_WATER','DEM_MOD_CONT_OCEAN','DEM_DEEP_OCEAN']
#DEM_dict = {
    #'DEM_SHALLOW_OCEAN'        : 0,
    #'DEM_LAND'                 : 1,
    #'DEM_COASTLINE'            : 2,
    #'DEM_SHALLOW_INLAND_WATER' : 3,
    #'DEM_EPHEMERAL_WATER'      : 4,
    #'DEM_DEEP_INLAND_WATER'    : 5,
    #'DEM_MOD_CONT_OCEAN'       : 6,
    #'DEM_DEEP_OCEAN'           : 7
#}

#DEM_to_QSTLWM = np.array([
    #QSTLWM_dict['OCEAN_SEA'],                    # DEM 0 : QSTLWM 17
    #QSTLWM_dict['EVERGREEN_NEEDLELEAF_FOREST'],  # DEM 1 : QSTLWM  1
    #QSTLWM_dict['COASTAL_WATER'],                # DEM 2 : QSTLWM 19
    #QSTLWM_dict['INLAND_WATER'],                 # DEM 3 : QSTLWM 18
    #QSTLWM_dict['SAVANNAS'],                     # DEM 4 : QSTLWM  9
    #QSTLWM_dict['OCEAN_SEA'],                    # DEM 5 : QSTLWM 17
    #QSTLWM_dict['OCEAN_SEA'],                    # DEM 6 : QSTLWM 17
    #QSTLWM_dict['OCEAN_SEA']                     # DEM 7 : QSTLWM 17
#],dtype=('uint8'))


#def index(a, x):
    #'''Locate the leftmost value exactly equal to x'''
    #i = bisect_left(a, x)
    #if i != len(a) and a[i] == x:
        #return i
    #raise ValueError

#def find_lt(a, x):
    #'''Find rightmost value less than x'''
    #i = bisect_left(a, x)
    #if i:
        #return a[i-1]
    #raise ValueError

#def find_le(a, x):
    #'''Find rightmost value less than or equal to x'''
    #i = bisect_right(a, x)
    #if i:
        #return a[i-1]
    #raise ValueError

#def find_gt(a, x):
    #'''Find leftmost value greater than x'''
    #i = bisect_right(a, x)
    #if i != len(a):
        #return a[i]
    #raise ValueError

#def find_ge(a, x):
    #'''Find leftmost item greater than or equal to x'''
    #i = bisect_left(a, x)
    #if i != len(a):
        #return a[i]
    #raise ValueError

#def toJulianDate(inTime):
    #'''
    #Takes time in regular yyyymmdd format and returns time string in Julian yyyyddd format.
    #'''
    #inTime = str(inTime)
    #try :
        #return time.strftime("%Y%j",time.strptime(inTime, "%Y%m%d"))
    #except :
        #LOG.error("Incorrect data format (%s). Should conform to yyyymmdd." % (inTime))
        #return 1

#def fromJulianDate(inTime):
    #'''
    #Takes time in Julian yyyyddd format and returns time string in regular yyyymmdd format .
    #'''
    #inTime = str(inTime)
    #try :
        #return time.strftime("%Y%m%d",time.strptime(inTime,"%Y%j"))
    #except :
        #LOG.error("Incorrect data format (%s). Should conform to yyyyddd." % (inTime))
        #return 1

#----------------------------------------------------------------------------
# Finds the places where the boundary points that will make up a polygon
# cross the dateline.
#
# This method is heavily based on the AltNN NNfind_crossings() method
#----------------------------------------------------------------------------
'''
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

    return num180Crossings_
'''

#def isDatelineCrossed(latCrnList,lonCrnList):

    #dateLineCrossed = []
    #ascendingNode = []
    #descendingNode = []

    #LOG.debug("Determining granule dateline crossings...")

    #for granLats,granLons in zip(latCrnList,lonCrnList):
        #LOG.debug("Processing granule...")
        #isDatelineCrosser = False
        #isAscendingNode = False
        #isDecendingNode = False

        ## TODO : Use the ADL method of determining dateline crossings
        #numCrossings = findDatelineCrossings(granLats,granLons)
        #LOG.debug("Number of dataline crossings %d" % (numCrossings))

        ## Ascending node ? ...
        #if (granLats[2] > granLats[0]):
            #LOG.debug("Ascending node...")
            #isAscendingNode = True
            ## Dateline crosser ? ...
            #if (granLons[0] < granLons[3]):
                #LOG.debug("Dateline crosser...\n")
                #isDatelineCrosser = True

        ## Descending node ? ...
        #if (granLats[0] > granLats[2]):
            #LOG.debug("Descending node...")
            #isDecendingNode = True
            ## Dateline crosser ? ...
            #if (granLons[1] < granLons[2]):
                #LOG.debug("Dateline crosser...\n")
                #isDatelineCrosser = True

        #dateLineCrossed.append(isDatelineCrosser)
        #ascendingNode.append(isAscendingNode)
        #descendingNode.append(isDecendingNode)

    #return dateLineCrossed,ascendingNode,descendingNode


def _create_input_file_globs(inputFiles):
    '''
    Determine the correct input file path and globs
    '''
    input_path = path.abspath(inputFiles)
    if path.isdir(input_path) :
        input_dir = input_path
        input_files = None
    else :
        input_dir = path.dirname(input_path)
        input_files = path.basename(input_path)

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
    elif path.isfile(input_path) :
        fileGlob = string.rstrip(string.lstrip(string.split(input_files,"b")[0],charsToKill),charsToKill)
        LOG.debug("fileGlob = %s" %(fileGlob))
        inputGlobs['GEO'] = "GMTCO%s*.h5" %(fileGlob)
        inputGlobs['MOD'] = "SVM*%s*.h5" %(fileGlob)
        inputGlobs['IMG'] = "SVI*%s*.h5" %(fileGlob)
        for fileType in ['GEO','MOD','IMG']:
            inputGlobs[fileType] = string.replace(inputGlobs[fileType],"**","*")

    return input_dir,inputGlobs


def _skim_viirs_sdr(collectionShortName,work_dir):
    "skim for VIIRS SDR data meeting minimum requirements"
    for info in skim_dir(work_dir, N_Collection_Short_Name=collectionShortName):
        print info
        path = info['BlobPath']
        if not path.isfile(path) or os.stat(path).st_size < MINIMUM_SDR_BLOB_SIZE:
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

    work_dir = path.abspath(work_dir)

    geoGroupList = list(_contiguous_granule_groups(skim_dir(work_dir, N_Collection_Short_Name=collectionShortName)))

    if len(geoGroupList)==0:
        return

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
#XML_TMPL_VIIRS_MASKS_EDR = """<InfTkConfig>
  #<idpProcessName>ProEdrViirsMasksController.exe</idpProcessName>
  #<siSoftwareId />
  #<isRestart>FALSE</isRestart>
  #<useExtSiWriter>FALSE</useExtSiWriter>
  #<debugLogLevel>LOW</debugLogLevel>
  #<debugLevel>DBG_HIGH</debugLevel>
  #<dbgDest>D_FILE</dbgDest>
  #<enablePerf>FALSE</enablePerf>
  #<perfPath>${WORK_DIR}/perf</perfPath>
  #<dbgPath>${WORK_DIR}/log</dbgPath>
  #<initData>
     #<domain>OPS</domain>
     #<subDomain>SUBDOMAIN</subDomain>
     #<startMode>INF_STARTMODE_COLD</startMode>
     #<executionMode>INF_EXEMODE_PRIMARY</executionMode>
     #<healthTimeoutPeriod>30</healthTimeoutPeriod>
  #</initData>
  #<lockinMem>FALSE</lockinMem>
  #<rootDir>${WORK_DIR}</rootDir>
  #<inputPath>${WORK_DIR}</inputPath>
  #<outputPath>${WORK_DIR}</outputPath>
  #<dataStartIET>0000000000000000</dataStartIET>
  #<dataEndIET>1111111111111111</dataEndIET>
  #<actualScans>47</actualScans>
  #<previousActualScans>48</previousActualScans>
  #<nextActualScans>48</nextActualScans> 
  #<usingMetadata>TRUE</usingMetadata>
  #<configGuideName>ProEdrViirsMasksController_GuideList.cfg</configGuideName>

  #<task>
    #<taskType>EDR</taskType>
    #<taskDetails1>%(N_Granule_ID)s</taskDetails1>
    #<taskDetails2>%(N_Granule_Version)s</taskDetails2>
    #<taskDetails3>NPP</taskDetails3>
    #<taskDetails4>VIIRS</taskDetails4>
  #</task>
#</InfTkConfig>
#"""

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

#def generate_viirs_masks_edr_xml(work_dir, granule_seq):
    #"generate XML files for VIIRS Masks EDR granule generation"
    #to_process = []
    #for gran in granule_seq:
        #name = gran['N_Granule_ID']
        #fnxml = 'edr_viirs_masks_%s.xml' % (name)
        #LOG.debug('writing XML file %r' % (fnxml))
        #fpxml = file(path.join(work_dir, fnxml), 'wt')
        #fpxml.write(XML_TMPL_VIIRS_MASKS_EDR % gran)
        #to_process.append([name,fnxml])
    #return to_process


def generate_viirs_aerosol_edr_xml(work_dir, granule_seq):
    "generate XML files for VIIRS Masks EDR granule generation"
    to_process = []
    for gran in granule_seq:
        name = gran['N_Granule_ID']
        fnxml = 'edr_viirs_aerosol_%s.xml' % (name)
        LOG.debug('writing XML file %r' % (fnxml))
        fpxml = file(path.join(work_dir, fnxml), 'wt')
        fpxml.write(XML_TMPL_VIIRS_AEROSOL_EDR % gran)
        to_process.append([name,fnxml])
    return to_process


def _granulate_GridIP(inDir,geoDicts):
    '''Granulates the input gridded static data into the required GridIP granulated datasets.'''

    import GridIP
    import Algorithms
    global ancEndian 

    CSPP_RT_HOME = os.getenv('CSPP_RT_HOME')
    ANC_SCRIPTS_PATH = path.join(CSPP_RT_HOME,'viirs')
    CSPP_RT_ANC_CACHE_DIR = os.getenv('CSPP_RT_ANC_CACHE_DIR')
    
    ADL_ASC_TEMPLATES = path.join(ANC_SCRIPTS_PATH,'asc_templates')

    # Ordered list of required algorithms (to be passed in)
    algList = ['VCM']
    #algList = ['VCM','AOT']
    
    # Collection shortnames of the required GridIP static ancillary datasets
    # FIXME : Poll ADL/cfg/ProEdrViirsCM_CFG.xml for this information

    # Create a list of algorithm module "pointers"
    algorithms = []
    for alg in algList :
        algName = Algorithms.modules[alg]
        algorithms.append(getattr(Algorithms,algName))

    # Obtain the required GridIP collection shortnames for each algorithm
    collectionShortNames = []
    for alg in algorithms :
        for shortName in alg.GridIP_collectionShortNames :
            LOG.info("Adding %s to the list of required collection short names..." \
                    %(shortName))
            collectionShortNames.append(shortName)

    # Remove duplicate shortNames
    collectionShortNames = list(set(collectionShortNames))
    LOG.info("collectionShortNames = %r" %(collectionShortNames))

    # Create a dict of GridIP class instances, which will handle ingest and granulation.
    GridIP_objects = {}
    for shortName in collectionShortNames :
        className = GridIP.classNames[shortName]
        GridIP_objects[shortName] = getattr(GridIP,className)(inDir=inDir)
        LOG.info("GridIP_objects[%s].blobDatasetName = %r" % (shortName,GridIP_objects[shortName].blobDatasetName))
        if (np.shape(GridIP_objects[shortName].collectionShortName) != () ):
            LOG.info("    GridIP_objects[%s].collectionShortName = %r" % (shortName,GridIP_objects[shortName].collectionShortName))
            LOG.info("    GridIP_objects[%s].xmlName = %r" % (shortName,GridIP_objects[shortName].xmlName))
            GridIP_objects[shortName].collectionShortName = shortName
            GridIP_objects[shortName].xmlName = GridIP_objects[shortName].xmlName[shortName]
            LOG.info("New GridIP_objects[%s].collectionShortName = %r" % (shortName,GridIP_objects[shortName].collectionShortName))
            LOG.info("New GridIP_objects[%s].xmlName = %r" % (shortName,GridIP_objects[shortName].xmlName))

    #sys.exit(0)

    # Loop through the required GridIP datasets and create the blobs.
    granIdKey = lambda x: (x['N_Granule_ID'])
    for dicts in sorted(geoDicts,key=granIdKey):
        for shortName in collectionShortNames :
        
            LOG.info("Processing dataset %s for %s" % (GridIP_objects[shortName].blobDatasetName,shortName))

            # Set the geolocation information in this ancillary object for the current granule...
            GridIP_objects[shortName].setGeolocationInfo(dicts)

            LOG.debug("min,max,range of latitude: %f %f %f" % (\
                    GridIP_objects[shortName].latMin,\
                    GridIP_objects[shortName].latMax,\
                    GridIP_objects[shortName].latRange))
            LOG.debug("min,max,range of longitude: %f %f %f" % (\
                    GridIP_objects[shortName].lonMin,\
                    GridIP_objects[shortName].lonMax,\
                    GridIP_objects[shortName].lonRange))

            LOG.debug("latitude corners: %r" % (GridIP_objects[shortName].latCrnList))
            LOG.debug("longitude corners: %r" % (GridIP_objects[shortName].lonCrnList))

            # Subset the gridded data for this ancillary object to cover the required lat/lon range.
            GridIP_objects[shortName].subset()

            # Granulate the gridded data in this ancillary object for the current granule...
            GridIP_objects[shortName].granulate(GridIP_objects)

            # Shipout the granulated data in this ancillary object to a blob/asc pair.
            GridIP_objects[shortName].shipOutToFile()


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
    CSPP_RT_ANC_CACHE_DIR = os.getenv('CSPP_RT_ANC_CACHE_DIR')

    ANC_SCRIPTS_PATH = path.join(CSPP_RT_HOME,'viirs')
    ADL_ASC_TEMPLATES = path.join(ANC_SCRIPTS_PATH,'asc_templates')

    #print "CSPP_RT_HOME =          %r" % (CSPP_RT_HOME)
    #print "CSPP_RT_ANC_CACHE_DIR = %r" % (CSPP_RT_ANC_CACHE_DIR)
    #print "ANC_SCRIPTS_PATH =      %r" % (ANC_SCRIPTS_PATH)
    #print "ADL_ASC_TEMPLATES =     %r" % (ADL_ASC_TEMPLATES)
    #print "ADL_HOME =              %r" % (ADL_HOME)

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

    auxillaryPaths = ['luts/viirs',
                      'luts/viirs',
                      'luts/viirs',
                      'luts/viirs',
                      'luts/viirs',
                      'NAAPS-ANC-Int',
                      'luts/viirs',
                      'luts/viirs',
                      'luts/viirs',
                      'luts/viirs']


    auxillarySourceFiles = []
    # Number of characters in a URID, plus the trailing "."...
    charsInUrid = 32+1

    for templatePath,blobTempFileName in zip(auxillaryPaths,auxillaryBlobTemplateFile) :
        blobTempFileName = path.join(CSPP_RT_ANC_CACHE_DIR,templatePath,blobTempFileName)
        if path.islink(blobTempFileName) :
            LOG.info("%s is a link, resolving auxillary filename..." %(blobTempFileName))
            auxillarySourceFile = path.basename(os.readlink(blobTempFileName))[charsInUrid:]
            auxillarySourceFiles.append(auxillarySourceFile)
        else :
            auxillarySourceFile = path.basename(blobTempFileName)
            auxillarySourceFiles.append(auxillarySourceFile)

        LOG.info("Auxillary filename : %s" %(auxillarySourceFile))

    for shortName,auxillarySourceFile in zip(auxillaryCollShortNames,auxillarySourceFiles) :
        LOG.info("%s --> %s" %(shortName,auxillarySourceFile))


    for shortName,ascTemplateFileName,blobTempFileName,templatePath,auxillarySourceFile in zip(auxillaryCollShortNames,auxillaryAscTemplateFile,auxillaryBlobTemplateFile,auxillaryPaths,auxillarySourceFiles):

        #LOG.info("Creating new %s asc file from template %s" % (shortName,ascTemplateFileName))

        # Create a new URID to be used in making the asc filenames

        URID_dict = _getURID()

        URID = URID_dict['URID']
        creationDate_nousecStr = URID_dict['creationDate_nousecStr']
        creationDateStr = URID_dict['creationDateStr']

        # The names for the new asc and blob files
        ascFileName = path.join(inDir,URID+'.asc')
        blobFileName = path.join(inDir,string.replace(blobTempFileName,'template',URID))

        # Make a new asc file from the template, and substitute for the various tags

        ascTemplateFileName = path.join(ADL_ASC_TEMPLATES,ascTemplateFileName)

        LOG.info("Creating new asc file %s from template %s" % \
                (path.basename(ascFileName),path.basename(ascTemplateFileName)))

        try:
            ascTemplateFile = open(ascTemplateFileName,"rt") # Open template file for reading
            ascFile = open(ascFileName,"wt") # create a new text file
        except Exception, err :
            LOG.error("%s, aborting." % (err))
            sys.exit(1)

        LOG.debug("Template file %s is %r with mode %s" %(ascTemplateFileName,'not open' if ascTemplateFile.closed else 'open',ascTemplateFile.mode))

        LOG.debug("New file %s is %r with mode %s" %(ascFileName,'not open' if ascFile.closed else 'open',ascFile.mode))

        for line in ascTemplateFile.readlines():
           line = line.replace("CSPP_URID",URID)
           line = line.replace("CSPP_CREATIONDATETIME_NOUSEC",creationDate_nousecStr)
           line = line.replace("CSPP_AUX_BLOB_FULLPATH",path.basename(blobFileName))
           line = line.replace("CSPP_CREATIONDATETIME",creationDateStr)
           line = line.replace("CSPP_AUX_SOURCE_FILE",auxillarySourceFile)
           ascFile.write(line) 

        ascFile.close()
        ascTemplateFile.close()

        # Create a link between the binary template file and working directory

        blobTempFileName = path.join(CSPP_RT_ANC_CACHE_DIR,templatePath,blobTempFileName)
        LOG.info("Creating the link %s -> %s" %(blobFileName,blobTempFileName))

        if not path.exists(blobFileName):
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


def _granulate_ANC(inDir,geoDicts):
    '''Granulates the input gridded blob files into the required ANC granulated datasets.'''

    import ANC
    import Algorithms
    global ancEndian 

    CSPP_RT_HOME = os.getenv('CSPP_RT_HOME')
    ANC_SCRIPTS_PATH = path.join(CSPP_RT_HOME,'viirs')
    CSPP_RT_ANC_CACHE_DIR = os.getenv('CSPP_RT_ANC_CACHE_DIR')
    
    ADL_ASC_TEMPLATES = path.join(ANC_SCRIPTS_PATH,'asc_templates')

    # Ordered list of required algorithms (to be passed in)
    algList = ['VCM']
    #algList = ['VCM','AOT']
    
    # Download the required NCEP grib files
    LOG.info("Downloading NCEP GRIB ancillary into cache...")
    gribFiles = ANC.retrieve_NCEP_grib_files(geoDicts)
    LOG.debug('dynamic ancillary GRIB files: %s' % repr(gribFiles))
    if (gribFiles == []) :
        LOG.error('Failed to find or retrieve any GRIB files, aborting.')
        sys.exit(1)

    # Transcode the NCEP GRIB files into ADL NCEP-ANC-Int
    gridBlobFiles = ANC.create_NCEP_grid_blobs(gribFiles)
    if (gridBlobFiles == []) :
        LOG.error('Failed to convert NCEP GRIB  files to blobs, aborting.')
        sys.exit(1)

    # Create a list of algorithm module "pointers"
    algorithms = []
    for alg in algList :
        algName = Algorithms.modules[alg]
        algorithms.append(getattr(Algorithms,algName))

    # Obtain the required ANC collection shortnames for each algorithm
    collectionShortNames = []
    for alg in algorithms :
        for shortName in alg.ANC_collectionShortNames :
            LOG.info("Adding %s to the list of required collection short names..." \
                    %(shortName))
            collectionShortNames.append(shortName)

    # Remove duplicate shortNames
    collectionShortNames = list(set(collectionShortNames))
    LOG.info("collectionShortNames = %r" %(collectionShortNames))

    # Create a dict of ANC class instances, which will handle ingest and 
    # granulation
    ANC_objects = {}
    for shortName in collectionShortNames :
        className = ANC.classNames[shortName]
        ANC_objects[shortName] = getattr(ANC,className)(inDir=inDir)
        LOG.info("ANC_objects[%s].blobDatasetName = %r" % (shortName,ANC_objects[shortName].blobDatasetName))

    # Open the NCEP gridded blob file
    # FIXME : Should be using two NCEP blob files, and averaging
    ncepXmlFile = path.join(ADL_HOME,'xml/ANC/NCEP_ANC_Int.xml')
    gridBlobFile = gridBlobFiles[0]

    if path.exists(ncepXmlFile):
        LOG.debug("We are using for %s: %s,%s" %('NCEP-ANC-Int',ncepXmlFile,gridBlobFile))
    
    endian = ancEndian
    ncepBlobObj = adl_blob.map(ncepXmlFile,gridBlobFile, endian=endian)
    ncepBlobArrObj = ncepBlobObj.as_arrays()
    LOG.debug("%s...\n%r" % (gridBlobFile,ncepBlobArrObj._fields))
    
    # Ingest the ANC gridded data and copy to the gridData attribute of the ANC objects
    for shortName in collectionShortNames :
        if ANC_objects[shortName].sourceType == 'NCEP_ANC_Int' :
            ANC_objects[shortName].sourceList = gridBlobFiles
            ANC_objects[shortName].ingest(ancBlob=ncepBlobArrObj)
        else :
            ANC_objects[shortName].ingest()
        LOG.info("Ingesting ANC_objects gridded  %s" % (shortName))

    # Loop through the required ANC datasets and create the blobs.
    granIdKey = lambda x: (x['N_Granule_ID'])
    for dicts in sorted(geoDicts,key=granIdKey):
        for shortName in collectionShortNames :
        
            LOG.info("Processing dataset %s for %s" % (ANC_objects[shortName].blobDatasetName,shortName))

            # Set the geolocation information in this ancillary object for the current granule...
            ANC_objects[shortName].setGeolocationInfo(dicts)

            # Granulate the gridded data in this ancillary object for the current granule...
            ANC_objects[shortName].granulate()

            # Shipout the granulated data in this ancillary object to a blob/asc pair.
            ANC_objects[shortName].shipOutToFile()


# Look through new log files for completed messages
# based on CrIS SDR original
#def check_log_files(work_dir, pid, xml):
    #"""
    #Find the log file
    #Look for success
    #Return True if found
    #Display log message and hint if problem occurred
    #"""
    ## retrieve exe name and log path from lw file
    #logpat = path.join(work_dir, "log", "*%d*.log" % pid)

    #n_err=0
    #for logFile in glob(logpat):
        #LOG.debug( "Checking Log file " +logFile +" for errors.")
        #n_err += adl_log.scan_log_file(COMMON_LOG_CHECK_TABLE, logFile)

    #if n_err == 0 :
        #LOG.debug(" Log: "+ logFile + " "+ xml + " Completed successfully" )
        #return True
    #else :
        #LOG.error( " Log: "+ logFile +" " + xml + " Completed unsuccessfully" )
        #return False


#def run_xml_files(work_dir, xml_files_to_process, setup_only=False, **additional_env):
    """Run each VIIRS EDR MASKS XML input in sequence
    return the list of granule IDs which crashed, and list of granule ids which did not create output
    """
    #crashed_runs = set()
    #no_output_runs = set()
    #geo_problem_runs = set()
    #bad_log_runs = set()
    #first = True

    # obtain pre-existing granule list
    #modGeoTCPattern = path.join(work_dir, 'GMTCO*.h5')
    #cmaskPattern = path.join(work_dir, 'IICMO*.h5')
    #activeFiresPattern = path.join(work_dir, 'AVAFO*.h5')
    # TODO : Enable for VIIRS AOT
    #aotIpPattern = path.join(work_dir, 'IVAOT*.h5')
    #aotEdrPattern = path.join(work_dir, 'VAOOO*.h5')
    #suspMatEdrPattern = path.join(work_dir, 'VSUMO*.h5')

    # prior_granules dicts contain (N_GranuleID,HDF5File) key,value pairs.
    # *ID dicts contain (HDF5File,N_GranuleID) key,value pairs.
    
    # Get the N_GranuleID and filename of existing cmask files in the work_dir ? ...
    #cmask_prior_granules, cmask_ID = h5_xdr_inventory(cmaskPattern, CM_GRANULE_ID_ATTR_PATH)
    #LOG.debug('Existing IICMO granules... %s' % (repr(cmask_prior_granules)))
    #cmask_prior_granules = set(cmask_prior_granules.keys())
    #LOG.debug(' Set of existing IICMO granules... %s' % (repr(cmask_prior_granules)))

    ## Get the (N_GranuleID,hdfFileName) pairs for the existing Active fires files
    #activeFires_prior_granules, afires_ID = h5_xdr_inventory(activeFiresPattern, AF_GRANULE_ID_ATTR_PATH)
    #activeFires_prior_granules = set(activeFires_prior_granules.keys())

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

    #for granule_id, xml in xml_files_to_process:

        #t1 = time()
        
        #cmd = [ADL_VIIRS_MASKS_EDR, xml]
        ## TODO : Enable for VIIRS AOT
        ##cmd = [ADL_VIIRS_AEROSOL_EDR, xml]
        
        #if setup_only:
            #print ' '.join(cmd)
        #else:
            #LOG.debug('executing "%s"' % ' '.join(cmd))
            #LOG.debug('additional environment variables: %s' % repr(additional_env))
            #try:
                #pid = sh(cmd, env=env(**additional_env), cwd=work_dir)
                #LOG.debug("%r ran as pid %d" % (cmd, pid))
                #if not check_log_files(work_dir, pid, xml):
                    #bad_log_runs.add(granule_id)

            #except CalledProcessError as oops:
                #LOG.debug(traceback.format_exc())
                #LOG.error('ProEdrViirsMasksController.exe failed on %r: %r. Continuing...' % (xml, oops))
                #crashed_runs.add(granule_id)
            #first = False

            ## check new IICMO output granules
            #cmask_new_granules, cmask_ID = h5_xdr_inventory(cmaskPattern, CM_GRANULE_ID_ATTR_PATH, state=cmask_ID)
            #LOG.debug('new IICMO granules after this run: %s' % (repr(cmask_new_granules)))
            #if granule_id not in cmask_new_granules:
                #LOG.warning('no IICMO HDF5 output for %s' % (granule_id))
                #no_output_runs.add(granule_id)
            #else:
                #filename = cmask_new_granules[granule_id]

            ## check new AVAFO output granules
            #afires_new_granules, afires_ID = h5_xdr_inventory(activeFiresPattern, AF_GRANULE_ID_ATTR_PATH, state=afires_ID)
            #LOG.debug('new AVAFO granules after this run: %s' % (repr(afires_new_granules)))
            #if granule_id not in afires_new_granules:
                #LOG.warning('no AVAFO HDF5 output for %s' % (granule_id))
                #no_output_runs.add(granule_id)
            #else:
                #filename = afires_new_granules[granule_id]

            ## TODO : Enable for VIIRS AOT
            ## check new IVAOT output granules
            ##aotIp_new_granules, aotIp_ID = h5_xdr_inventory(aotIpPattern, AOT_IP_GRANULE_ID_ATTR_PATH, state=aotIp_ID)
            ##LOG.debug('new IVAOT granules after this run: %s' % repr(aotIp_new_granules))
            ##if granule_id not in aotIp_new_granules:
                ##LOG.warning('no IVAOT HDF5 output for %s' % granule_id)
                ##no_output_runs.add(granule_id)
            ##else:
                ##filename = aotIp_new_granules[granule_id]

            ## TODO : Enable for VIIRS AOT
            ## check new VAOOO output granules
            ##aotEdr_new_granules, aotEdr_ID = h5_xdr_inventory(aotEdrPattern, AOT_EDR_GRANULE_ID_ATTR_PATH, state=aotEdr_ID)
            ##LOG.debug('new VAOOO granules after this run: %s' % repr(aotEdr_new_granules))
            ##if granule_id not in aotEdr_new_granules:
                ##LOG.warning('no VAOOO HDF5 output for %s' % granule_id)
                ##no_output_runs.add(granule_id)
            ##else:
                ##filename = aotEdr_new_granules[granule_id]

            ## TODO : Enable for VIIRS AOT
            ## check new VSUMO output granules
            ##suspMatEdr_new_granules, suspMatEdr_ID = h5_xdr_inventory(suspMatEdrPattern, SUSMAT_EDR_GRANULE_ID_ATTR_PATH, state=suspMatEdr_ID)
            ##LOG.debug('new VSUMO granules after this run: %s' % repr(suspMatEdr_new_granules))
            ##if granule_id not in suspMatEdr_new_granules:
                ##LOG.warning('no VSUMO HDF5 output for %s' % granule_id)
                ##no_output_runs.add(granule_id)
            ##else:
                ##filename = suspMatEdr_new_granules[granule_id]

        #t2 = time()
        #LOG.info ( "Controller ran in %f seconds." % (t2-t1))

    #LOG.debug("cmask_ID.values() = %r" % (cmask_ID.values()))
    #LOG.debug("set(cmask_ID.values()) = %r" % (set(cmask_ID.values())))
    #LOG.debug("cmask_prior_granules = %r" % (cmask_prior_granules))
    #cmask_granules_made = set(cmask_ID.values()) - cmask_prior_granules

    #LOG.debug("afires_ID.values() = %r" % (afires_ID.values()))
    #LOG.debug("set(afires_ID.values()) = %r" % (set(afires_ID.values())))
    #LOG.debug("activeFires_prior_granules = %r" % (activeFires_prior_granules))
    #activeFires_granules_made = set(afires_ID.values()) - activeFires_prior_granules

    ## TODO : Enable for VIIRS AOT
    ##LOG.debug("aotIp_ID.values() = \n%r" % (aotIp_ID.values()))
    ##LOG.debug("set(aotIp_ID.values()) = \n%r" % (set(aotIp_ID.values())))
    ##LOG.debug("aerosolIP_prior_granules = \n%r" % (aerosolIP_prior_granules))
    ##aerosolIP_granules_made = set(aotIp_ID.values()) - aerosolIP_prior_granules

    ## TODO : Enable for VIIRS AOT
    ##LOG.debug("aotEdr_ID.values() = \n%r" % (aotEdr_ID.values()))
    ##LOG.debug("set(aotEdr_ID.values()) = \n%r" % (set(aotEdr_ID.values())))
    ##LOG.debug("aerosolEDR_prior_granules = \n%r" % (aerosolEDR_prior_granules))
    ##aerosolEDR_granules_made = set(aotEdr_ID.values()) - aerosolEDR_prior_granules

    ## TODO : Enable for VIIRS AOT
    ##LOG.debug("suspMatEdr_ID.values() = \n%r" % (suspMatEdr_ID.values()))
    ##LOG.debug("set(suspMatEdr_ID.values()) = \n%r" % (set(suspMatEdr_ID.values())))
    ##LOG.debug("suspMatEDR_prior_granules = \n%r" % (suspMatEDR_prior_granules))
    ##suspMatEDR_granules_made = set(suspMatEdr_ID.values()) - suspMatEDR_prior_granules


    #LOG.info('cmask granules created: %s' %( ', '.join(list(cmask_granules_made))))
    #LOG.info('activeFires granules created: %s' % (', '.join(list(activeFires_granules_made))))
    ##LOG.info('aerosolIP granules created: %s' % ', '.join(list(aerosolIP_granules_made)))
    ##LOG.info('aerosolEDR granules created: %s' % ', '.join(list(aerosolEDR_granules_made)))
    ##LOG.info('suspMatEDR granules created: %s' % ', '.join(list(suspMatEDR_granules_made)))

    #if no_output_runs:
        #LOG.info('granules that failed to generate output: %s' % (', '.join(no_output_runs)))
    #if geo_problem_runs:
        #LOG.warning('granules which had no N_Geo_Ref: %s' % ', '.join(geo_problem_runs))
    #if crashed_runs:
        #LOG.warning('granules that crashed ADL: %s' % (', '.join(crashed_runs)))
    #if bad_log_runs:
        #LOG.warning('granules that produced logs indicating problems: %s' % (', '.join(bad_log_runs)))
    #if not cmask_granules_made:
        #LOG.warning('no HDF5 SDRs were created')

    #return crashed_runs, no_output_runs, geo_problem_runs, bad_log_runs


def _cleanup(work_dir, xml_glob, log_dir_glob, *more_dirs):
    """upon successful run, clean out work directory"""

    LOG.info("cleaning up work directory...")
    for fn in glob(path.join(work_dir, '????????-?????-????????-????????.*')):
        LOG.debug('removing %s' % (fn))
        os.unlink(fn)

    LOG.info("Removing task xml files...")
    for fn in glob(path.join(work_dir, xml_glob)):
        LOG.debug('removing task file %s' % (fn))
        os.unlink(fn)

    LOG.info("Removing log directories %s ..."%(log_dir_glob))
    for dirname in glob(path.join(work_dir,log_dir_glob)):
        LOG.debug('removing logs in %s' % (dirname))
        try :
            rmtree(dirname, ignore_errors=False)
        except Exception, err:
            LOG.warn( "%s" % (str(err)))

    # FIXME : Cannot seem to remove the log directory, complains that the dir is not empty.
    LOG.info("Removing other directories ...")
    for dirname in more_dirs:
        fullDirName = path.join(work_dir,dirname)
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

    # Determine the correct input file path and glob

    if not options.skipSdrUnpack :
        input_dir,inputGlobs = _create_input_file_globs(options.inputFiles)
        if inputGlobs['GEO'] is None or inputGlobs['MOD'] is None or inputGlobs['IMG'] is None :
            LOG.error("No input files found matching %s, aborting..."%(options.inputFiles))
            sys.exit(1)

    # Set the work directory

    work_dir = path.abspath(options.work_dir)
    LOG.debug('Setting the work directory to %r' % (work_dir))

    # Check the environment variables, and whether we can write to the working directory
    check_env(work_dir)
    
    # create work directory
    if not path.isdir(work_dir):
        LOG.info('creating directory %s' % (work_dir))
        os.makedirs(work_dir)
    log_dir = path.join(work_dir, 'log')
    if not path.isdir(log_dir):
        LOG.debug('creating directory %s' % (log_dir))
        os.makedirs(log_dir)
    #perf_dir = path.join(work_dir, 'perf')
    #if not path.isdir(perf_dir):
        #LOG.debug('creating directory %s' % (perf_dir))
        #os.makedirs(perf_dir)
    #anc_dir = path.join(work_dir, 'linked_data')
    #if not path.isdir(anc_dir):
        #LOG.debug('creating directory %s' % (anc_dir))
        #os.makedirs(anc_dir)

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
    anc_granules_to_process = None
    granules_to_process = None

    for geoType in geolocationShortNames :
        LOG.info("Searching for VIIRS geolocation %s for ancillary granule IDs..." % (geoType))
        anc_granules_to_process = sorted(list(sift_metadata_for_viirs_sdr(geoType,crossGran=None,work_dir=work_dir)))

        if anc_granules_to_process :
            LOG.info("Searching for VIIRS geolocation %s for product granule IDs..." % (geoType))
            granules_to_process = sorted(list(sift_metadata_for_viirs_sdr(geoType,crossGran=1,   work_dir=work_dir)))
        
        if granules_to_process :
            break
        else :
            LOG.info("\tNo granules for VIIRS geolocation %s" % (geoType))

    LOG.info('Finished sifting through metadata to find VIIRS SDR processing candidates')

    if granules_to_process :
        granIdKey = lambda x: (x['N_Granule_ID'])

        ancGransID = ', '.join(x['N_Granule_ID'] for x in sorted(anc_granules_to_process,key=granIdKey))
        LOG.info("Ancillary granule IDs: %s"%(ancGransID))

        prodGransID = ', '.join(x['N_Granule_ID'] for x in sorted(granules_to_process,key=granIdKey))
        LOG.info("Product granule IDs: %s"%(', '.join(x['N_Granule_ID'] for x in granules_to_process)))
    else :
        LOG.error("Found no granules to process!")
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

    LOG.info("%13s%28s%29s" % ('N_Granule_ID','ObservedStartTime','ObservedEndTime'))
    for dicts in granules_to_process :
        LOG.info("%15s%30s%30s"%(dicts['N_Granule_ID'],dicts['ObservedStartTime'],dicts['ObservedEndTime']))

    # Expand any user specifiers in the various paths

    CSPP_RT_ANC_PATH = os.getenv('CSPP_RT_ANC_PATH')
    CSPP_RT_ANC_CACHE_DIR = os.getenv('CSPP_RT_ANC_CACHE_DIR')
    CSPP_RT_ANC_TILE_PATH = os.getenv('CSPP_RT_ANC_TILE_PATH')
    LD_LIBRARY_PATH = os.getenv('LD_LIBRARY_PATH')
    DSTATICDATA = os.getenv('DSTATICDATA')

    LOG.debug("\nCSPP_RT_ANC_PATH:       %s" % (CSPP_RT_ANC_PATH))
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

        # Granulate the VIIRS ANC data

        _granulate_ANC(work_dir,anc_granules_to_process)

        # Granulate the VIIRS GridIP data

        _granulate_GridIP(work_dir,anc_granules_to_process)

        t2 = time()
        LOG.info( "Generation of VIIRS Masks ancillary data took %f seconds." % (t2-t1))

    else :

        LOG.info('Skipping retrieval and granulation of ancillary data.')


    # Specify the algorithm we want to run via the Lw XML file.

    if not options.skipAlgorithm :

        import Algorithms

        # Ordered list of required algorithms (to be passed in)
        algList = ['VCM']
        #algList = ['VCM','AOT']

        # Create a list of algorithm module "pointers"
        algorithms = []
        for alg in algList :
            algName = Algorithms.modules[alg]
            algorithms.append(getattr(Algorithms,algName))

        for alg in algorithms :

            # build XML configuration files for jobs that can be run
            LOG.info("Building %s XML files for %d granules" % \
                    (alg.AlgorithmString,len(granules_to_process)))

            # Generate the VIIRS Cloud Mask/Active Fires LW xml files for each N_Granule_ID...
            xml_files_to_process = alg.generate_viirs_edr_xml(work_dir, granules_to_process)
        
            # Generate the VIIRS Cloud Mask/Active Fires LW xml files for each N_Granule_ID...
            #xml_files_to_process = generate_viirs_masks_edr_xml(work_dir, granules_to_process)

            # TODO : Enable for VIIRS AOT
            # Generate the VIIRS Aerosol Optical Thickness LW xml files for each N_Granule_ID...
            #xml_files_to_process = generate_viirs_aerosol_edr_xml(work_dir, granules_to_process)

            LOG.info('%d granules to process: %s' % (len(xml_files_to_process), ''.join(name+' -> '+xmlfile+'\n' for (name,xmlfile) in xml_files_to_process)))

            LOG.info("Running VIIRS EDR Masks on XML files")
            crashed_runs, no_output_runs, geo_problem_runs, bad_log_runs = \
                    alg.run_xml_files(work_dir, xml_files_to_process, \
                    setup_only = False, WORK_DIR = work_dir, \
                    LINKED_ANCILLARY = work_dir, ADL_HOME=ADL_HOME)

            ## considered a noncritical problem if there were any crashed runs, runs that produced no output,
            ## runs where Geo failed, or runs where ADL logs indicated a problem
            noncritical_problem = crashed_runs or no_output_runs or geo_problem_runs or bad_log_runs

            ## print final disposition message and get return code
            environment_error=False
            rc = get_return_code(unpacking_problems, len(xml_files_to_process), \
                    len(no_output_runs), noncritical_problem, environment_error)

        ## if no errors or only non-critical errors: clean up
        if rc == 0 and not options.cspp_debug:
            LOG.debug("Cleaning up workspace...")
            _cleanup(work_dir, 'edr_viirs_masks*.xml', 'ProEdrViirsMasksController.exe_*', \
                    log_dir, anc_dir, perf_dir)

            # TODO : Enable for VIIRS AOT
            #_cleanup(work_dir, 'edr_viirs_aerosol*.xml', 'ProEdrViirsAerosolController.exe_*', log_dir, anc_dir, perf_dir)

        return rc
    
    else :

        LOG.info('Skipping execution of VIIRS Masks Controller.')


if __name__=='__main__':
    sys.exit(main())  
