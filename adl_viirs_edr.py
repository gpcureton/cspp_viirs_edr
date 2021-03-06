#!/usr/bin/env python
# encoding: utf-8
"""
adl_viirs_edr.py

Purpose: Run one or more of the VIIRS EDR Controllers using ADL.

Input:
    * One or more HDF5 VIIRS SDR input files, aggregated or single-granule.
      You need one granule before and one granule after a given granule to 
      process.
    * Static IGBP,LWM,NDVI files.
    * A work directory, typically, empty, in which to unpack the granules and 
      generate the output. If the work directory specified does not exist, it 
      will be created.

Output:
    * ADL VIIRS Cloud Mask Intermediate Product (IP) blob files, and associated 
      HDF5 output files
    * ADL VIIRS Active Fires Intermediate Product (IP) blob files, and 
      associated HDF5 output files 
    * Associated .asc metadata files.
    * Unless --debug is specified, all blob and asc files, and the log, perf 
      and ancillary directories will be deleted.

Details:
    * If you input a series of granules, the software will scan the work 
      directory.
      Thus, for N input SDR granules, you may get up to N-2 or less output EDR 
      granules if 
      they are all contiguous, depending on the algorithm.
    * It is ambiguous to provide several copies of the same granule in the work 
      directory; this will result in an error abort.

Preconditions:
    * Requires ADL_HOME, CSPP_RT_ANC_PATH, CSPP_RT_ANC_CACHE_DIR environment 
      variables are set.
    * Requires that any needed LD_LIBRARY_PATH is set.
    * Requires that DSTATICDATA is set.

Optional:
    * Environment variables CSPP_RT_ANC_PATH, CSPP_RT_ANC_CACHE_DIR, if static 
      ancillary data is not placed within ADL_HOME.

Minimum commandline:

    python adl_viirs_edr.py  --input_files=INPUTFILES

where...

    INPUTFILES: The fully qualified path to the input files. May be a directory 
    or a file glob.


Created by Geoff Cureton <geoff.cureton@ssec.wisc.edu> on 2011-09-30.
Copyright (c) 2011-2013 University of Wisconsin Regents. All rights reserved.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
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
from shutil import rmtree,copyfile
from glob import glob
from time import time
from datetime import datetime,timedelta

import numpy as np
from numpy import ma
import copy

import tables as pytables
from tables import exceptions as pyEx

from multiprocessing import Pool, Lock, Value, cpu_count

# skim and convert routines for reading .asc metadata fields of interest
import adl_blob2 as adl_blob
from adl_asc import skim_dir, contiguous_granule_groups,eliminate_duplicates, \
        granule_groups_contain, effective_anc_contains, \
        _is_contiguous, corresponding_asc_path, RDR_REQUIRED_KEYS, \
        POLARWANDER_REQUIRED_KEYS
from adl_common import anc_files_needed, link_ancillary_to_work_dir, unpack, \
        env, h5_xdr_inventory, get_return_code, check_env
from adl_common import check_and_convert_path,check_existing_env_var
from adl_common import ADL_HOME, ADL_VARS, CSPP_RT_HOME, CSPP_RT_ANC_PATH, \
         CSPP_RT_ANC_CACHE_DIR
from adl_post_process import repack_products, aggregate_products
import adl_log

# every module should have a LOG object
sourcename= file_Id.split(" ")
LOG = logging.getLogger(sourcename[1])
from adl_common import configure_logging
from adl_common import _test_logging as test_logging

import Algorithms

# we're not likely to succeed in processing using geolocation smaller than this 
# many bytes
MINIMUM_SDR_BLOB_SIZE = 81000000

# maximum delay between granule end time and next granule start time to 
# consider them contiguous
MAX_CONTIGUOUS_DELTA=timedelta(seconds = 5)


###################################################
#                  Global Data                    #
###################################################

environ['TZ'] = 'UTC'
hexPat = '[\\dA-Fa-f]'

def set_sdr_endian(inputEndianness) :
    '''
    Set the global sdr endianness variable.
    '''
    global sdrEndian 
    if inputEndianness=='big' :
        sdrEndian = adl_blob.BIG_ENDIAN
    elif inputEndianness=='little' :
        sdrEndian = adl_blob.LITTLE_ENDIAN
    else :
        LOG.error('Invalid value for the VIIRS SDR endianness : %s ' 
                % (inputEndianness))


def set_anc_endian(inputEndianness) :
    '''
    Set the global ancillary endianness variable.
    '''
    global ancEndian 
    # FIXME: Use adl_blob here...    
    if inputEndianness=='big' :
        ancEndian = adl_blob.BIG_ENDIAN
    elif inputEndianness=='little' :
        ancEndian = adl_blob.LITTLE_ENDIAN
    else :
        LOG.error('Invalid value for the VIIRS ancillary endianness : %s ' 
                % (inputEndianness))

    LOG.info("ancEndian is {}".format(ancEndian))


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

    inputGlob = None

    charsToKill = string.ascii_letters + string.digits + "."

    if (input_files is None):
        # Input file glob is of form "/path/to/files"
        LOG.debug('Path1')
        inputGlob = '*.h5'

    elif path.isfile(input_path) :
        # Input file glob is of form "/path/to/files/GMTCO_npp_d_t_e_b_c_cspp_dev.h5" 
        LOG.debug('Path2')
        fileGlob = string.rstrip(string.lstrip(string.split(input_files,"b")[0],charsToKill),charsToKill)
        LOG.debug("fileGlob = %s" %(fileGlob))
        inputGlob = "*%s*.h5" %(fileGlob)
        LOG.debug("Initial inputGlob = %s" %(inputGlob))
        while (string.find(inputGlob,"**")!= -1): 
            inputGlob = string.replace(inputGlob,"**","*")
            LOG.debug("New inputGlob = %s" %(inputGlob))

    elif ("*" in input_files):
        # Input file glob is of form "/path/to/files/*"
        LOG.debug('Path3')
        fileGlob = string.rstrip(string.lstrip(string.split(input_files,"b")[0],charsToKill),charsToKill)
        inputGlob = "*%s*.h5" %(fileGlob)
        LOG.debug("Initial inputGlob = %s" %(inputGlob))
        while (string.find(inputGlob,"**")!= -1): 
            inputGlob = string.replace(inputGlob,"**","*")
            LOG.debug("New inputGlob = %s" %(inputGlob))

    return input_dir,inputGlob


def _get_alg_cross_granules(algList,noAlgChain):
    # Determine the number of VIIRS SDR cross granules required to 
    # process the algorithm chain.
    cumulativeCrossGranules = {}
    for alg in algList :
        if not noAlgChain : 
            crossSum = Algorithms.crossGranules[alg]
            for preReq in Algorithms.prerequisites[alg]:
                if preReq is None :
                    pass
                else :
                    crossSum += Algorithms.crossGranules[preReq]
            cumulativeCrossGranules[alg] = crossSum
        else :
            cumulativeCrossGranules[alg] = Algorithms.crossGranules[alg]

        LOG.info("We require {} cross granules for {}".format(cumulativeCrossGranules[alg],alg))

    return cumulativeCrossGranules

def _get_geo_prefixes(algorithms):
    '''
    Determine the correct geolocation file HDF5 prefixes for unpacking
    '''
    # Determine what geolocation types are required for each algorithm
    requiredGeoShortname = []
    requiredGeoPrefix = []
    for alg in algorithms :
        for shortName in alg.GEO_collectionShortNames :
            LOG.info("Algorithm '%s' requires geolocation type %r (%s*.h5)" % \
                    (alg.AlgorithmName,shortName,Algorithms.geo_hdf5_prefix[shortName]))
            requiredGeoShortname.append(shortName)

    requiredGeoShortname = list(set(requiredGeoShortname))

    for shortName in requiredGeoShortname :
        requiredGeoPrefix.append(Algorithms.geo_hdf5_prefix[shortName])

    LOG.info('Required geolocation shortnames: %r' % (requiredGeoShortname))
    LOG.info('Required geolocation prefixes: %r' % (requiredGeoPrefix))

    return requiredGeoShortname,requiredGeoPrefix


def _get_radio_prefixes(algorithms):
    '''
    Determine the correct radiometric file HDF5 prefixes for unpacking
    '''
    # Determine what radiometric types are required for each algorithm
    requiredSdrShortname = []
    requiredSdrPrefix = []
    for alg in algorithms :
        for shortName in alg.SDR_collectionShortNames :
            LOG.info("Algorithm '%s' requires radiometric type %r (%s*.h5)" % \
                    (alg.AlgorithmName,shortName,Algorithms.sdr_hdf5_prefix[shortName]))
            requiredSdrShortname.append(shortName)

    requiredSdrShortname = list(set(requiredSdrShortname))
    requiredSdrShortname.sort()

    for shortName in requiredSdrShortname :
        requiredSdrPrefix.append(Algorithms.sdr_hdf5_prefix[shortName])

    LOG.info('Required radiometric shortnames: %r' % (requiredSdrShortname))
    LOG.info('Required radiometric prefixes: %r' % (requiredSdrPrefix))

    return requiredSdrShortname,requiredSdrPrefix


def _unpack_sdr(work_dir,input_dir,inputGlob):
    '''
    Unpack HDF5 VIIRS SDRs from the input directory to the work directory
    '''
    unpacking_problems = 0

    fileGlob = path.join(input_dir,inputGlob)
    h5_names = glob(fileGlob)
    h5_names.sort()

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

    LOG.debug("VIIRS SDR unpacking problems = %d" % (unpacking_problems))

    return unpacking_problems


def _get_sdr_blobs(inputFiles,work_dir,skipSdrUnpack,
        requiredShortname,requiredPrefix):
    """
    Unpack HDF5 SDR files from the inputFiles into the working directory, or 
    link the blob/asc pairs from inputFiles to the working directory.
    """

    if skipSdrUnpack :
        # Assume there are blobs at the input directory...
        unpacking_problems = 0

        # Determine the correct input file path and glob
        input_dir,_ = _create_input_file_globs(inputFiles)
        LOG.info("input_dir = {}".format(input_dir))

        # Link the VIIRS SDR blob/asc pairs in the input directory 
        # to the work directory

        for shortName in requiredShortname :

            # A key for sorting lists of granule dictionaries according to N_Granule_ID
            granIdKey = lambda x: (x['N_Granule_ID'])
            
            sdr_granules_to_process = sorted(list(sift_metadata_for_viirs_sdr(
                shortName,crossGran=None,work_dir=input_dir)))

            for dicts in sorted(sdr_granules_to_process,key=granIdKey) :
                if 'BlobPath' not in dicts.keys():
                    dicts['BlobPath'] = None 
                else:
                    pass 

                asc_file = dicts['_filename']
                blob_file = dicts['BlobPath']
                blob_file = path.join(input_dir,
                        path.basename(dicts['BlobPath'])) \
                                if dicts['BlobPath'] is not None else None

                LOG.debug("asc_file = {}".format(asc_file))
                LOG.debug("blob_file = {}".format(blob_file))

                local_asc_file = path.join(work_dir,
                        path.basename(asc_file))
                local_blob_file = path.join(work_dir,
                    path.basename(blob_file)) \
                            if blob_file is not None else None

                LOG.debug("local_asc_file = {}".format(local_asc_file))
                LOG.debug("local_blob_file = {}".format(local_blob_file))

                # Linking in the asc file
                if asc_file is not None and local_asc_file is not None:
                    if not path.exists(local_asc_file):
                        LOG.debug("Creating the link {} -> {}"
                                .format(local_asc_file,asc_file))
                        os.symlink(asc_file, local_asc_file)
                    else:
                        LOG.debug('{} already exists; continuing'
                                .format(local_asc_file))

                    try:
                        LOG.debug('testing {}'.format(local_asc_file))
                        s = os.stat(local_asc_file)
                    except OSError as oops:
                        LOG.error("link at {} is broken".format(local_asc_file))
                        raise

                # Linking in the blob file
                if blob_file is not None and local_blob_file is not None:

                    if not path.exists(local_blob_file):
                        LOG.debug("Creating the link {} -> {}"
                                .format(local_blob_file,blob_file))
                        os.symlink(blob_file, local_blob_file)
                    else:
                        LOG.debug('{} already exists; continuing'
                                .format(local_blob_file))

                    try:
                        LOG.debug('testing {}'.format(local_blob_file))
                        s = os.stat(local_blob_file)
                    except OSError as oops:
                        LOG.error("link at {} is broken"
                                .format(local_blob_file))
                        raise

    else:

        # Determine the correct input file path and glob
        input_dir,inputGlob = _create_input_file_globs(inputFiles)

        if inputGlob is None :
            LOG.error("No input files found matching %s, aborting..."%(inputFiles))
            sys.exit(1)

        # Unpack HDF5 VIIRS SDRs in the input directory to the work directory
        unpacking_problems = 0
        for prefix in requiredPrefix:
            LOG.info("prefix = {}".format(prefix))
            LOG.info("inputGlob = {}".format(inputGlob))
            fileGlob = "%s%s" % (prefix,inputGlob)
            LOG.info("Unpacking files matching {}".format(fileGlob))
            unpacking_problems += _unpack_sdr(work_dir,input_dir,fileGlob)

    return unpacking_problems


def _skim_viirs_sdr(collectionShortName,work_dir):
    """
    Skim for VIIRS SDR data meeting minimum requirements
    """
    for info in skim_dir(work_dir, N_Collection_Short_Name=collectionShortName):
        print info
        path = info['BlobPath']
        if not path.isfile(path) or os.stat(path).st_size < MINIMUM_SDR_BLOB_SIZE:
            LOG.info('Ignoring %s, invalid blob. For direct broadcast this can happen at start or end of pass.' % (path))
        else:
            yield info


def __is_contiguous(grana, granb, tolerance=MAX_CONTIGUOUS_DELTA):
    """
    This is a custom version of adl_asc._is_contiguous(), which keys off of 
    'StartTime' and 'EndTime', rather than 'ObservedStartTime' and 
    'ObservedEndTime' as is done here.
    """
    delta = granb['ObservedStartTime'] - grana['ObservedEndTime']
    LOG.debug('delta of %s between %r and %r' % (delta, grana['N_Granule_ID'], 
        granb['N_Granule_ID']))
    if (delta >= tolerance) :
        LOG.info( "Granule gap of %s sec end %s to start %s  " % (str(delta),
            grana['ObservedStartTime'].isoformat(),
            granb['ObservedEndTime'].isoformat()))
    return (delta < tolerance)


def _contiguous_granule_groups(granules, tolerance=MAX_CONTIGUOUS_DELTA, larger_granules_preferred=False):
    """
    given a sequence of granule dictionaries, yield a sequence of contiguous 
    granule groups as tuples tolerance, if provided, is a datetime.timedelta 
    object representing max differance between endtime and starttime

    This is a custom version of adl_asc.contiguous_granule_groups(), which keys 
    off of 'StartTime', rather than 'ObservedStartTime' as is done here.
    """
    
    # sort granules into start time order and eliminate exact duplicates
    # FUTURE: is lex-compare sufficient for A2/A1/etc
    #start_time_key = lambda x: (x['StartTime'], x.get('N_Granule_Version', None))
    start_time_key = lambda x: (x['ObservedStartTime'], 
            x.get('N_Granule_Version', None))
    granlist = eliminate_duplicates(sorted(granules, key = start_time_key),
            larger_granules_preferred=larger_granules_preferred)
    granset = set(x['N_Granule_ID'] for x in granlist)

    # it's ambiguous if we have a work directory with multiple different blobs for any given granule
    if len(granlist) != len(granset):
        LOG.error('Aborting. Multiple granule copies in work directory: %d files in directory representing only %d granules' % (len(granlist), len(granset)))
        LOG.error('To solve this, only unpack granules once in given work directory')
        return

    # pair off in start-time order
    # build a set containing contiguous granules
    # when we find a break in sequence, yield the set as an ordered tuple and start over
    #
    # If there is only one granule in the list, yield it immediately
    seq = {}
    if len(granset)==1 :

        a = granlist[0]
        seq[a['URID']] = a
        LOG.debug('contiguous sequence has %d granules' % (len(seq)))
        yield tuple(sorted(seq.values(), key=start_time_key))
        seq.clear()

    else :

        for a,b in zip(granlist[:-1], granlist[1:]):
            if a['N_Granule_ID']==b['N_Granule_ID']:
                LOG.error('Granule %r has been unpacked to this directory multiple times!' % (a['N_Granule_ID']))
                return
            if __is_contiguous(a, b, tolerance):
                seq[a['URID']] = a
                seq[b['URID']] = b
            else:
                LOG.debug('contiguous sequence has %d granules' % (len(seq)))
                yield tuple(sorted(seq.values(), key=start_time_key))
                seq.clear()
        # leftovers! yum!
        if seq:
            LOG.debug('Leftover contiguous sequence has %d granules' % (len(seq)))
            yield tuple(sorted(seq.values(), key=start_time_key))


def sift_metadata_for_viirs_sdr(collectionShortName, crossGran=None, work_dir='.'):
    """
    Search through the ASC metadata in a directory, grouping in StartTime order.
    Look for back-to-back granules and break into contiguous sequences.
    Yield sequence of granules we can process.
    """
    LOG.debug('Collecting information for VIIRS SDRs')

    work_dir = path.abspath(work_dir)

    geoGroupList = list(_contiguous_granule_groups(skim_dir(work_dir, required_keys=RDR_REQUIRED_KEYS, N_Collection_Short_Name=collectionShortName)))

    LOG.debug('geoGroupList : {}'.format(geoGroupList))

    if len(geoGroupList)==0:
        LOG.debug('No geoGroupList found...')
        return

    # Loop through the contigous granule groups 
    for group in _contiguous_granule_groups(skim_dir(work_dir, 
        N_Collection_Short_Name=collectionShortName)):
        ##- for VIIRS, we can process everything but the first and last granule
        ##- for CrIS, use [4:-4]
        LOG.debug('Contiguous granule group of length: {}'.format(len(group),))

        if not crossGran :
            startGran,endGran = None,None
        else :
            startGran,endGran = crossGran,-1*crossGran

        if startGran is not None:
            granIdx = startGran
        else :
            granIdx = None

        for gran in group[startGran:endGran]:
            if not granule_groups_contain(geoGroupList, gran):
                LOG.info("Insufficient VIIRS SDR coverage to process {} @ {} ({}) - skipping".format(gran['N_Granule_ID'], gran['StartTime'], gran['URID']))
                continue

            # If we have cross granules, add them to the dictionary
            if granIdx is not None:
                granIdx_prev = granIdx - 1
                granIdx_next = granIdx + 1
                granId_prev = group[granIdx_prev]['N_Granule_ID']
                granId = gran['N_Granule_ID']
                granId_next = group[granIdx_next]['N_Granule_ID']
                gran['N_Granule_ID_prev'] = granId_prev
                gran['N_Granule_ID_next'] = granId_next
                granIdx = granIdx + 1
                LOG.debug("Granule IDs are ({},{},{})".format(gran['N_Granule_ID_prev'],gran['N_Granule_ID'],gran['N_Granule_ID_next']))
            else:
                gran['N_Granule_ID_prev'] = startGran
                gran['N_Granule_ID_next'] = endGran

            LOG.debug('Processing opportunity: {} at {} with uuid {}'.format(gran['N_Granule_ID'], gran['StartTime'], gran['URID']))
            yield gran


def _strReplace(fileName,oldString,newString):
    """
    Check fileName for occurences of oldString, if found fileName is opened and oldString is 
    replaced with newString.
    """
    fileChanged=0
    try :
        with open(fileName) as thefile:
            content = thefile.read()                 # read entire file into memory
            replacedText = content.replace(oldString, newString)
        if replacedText != content:
            LOG.debug('Replacing occurence of "%s" in %s with "%s"' % (oldString,path.basename(fileName),newString))
            with open(fileName, 'w') as thefile:
                thefile.write(replacedText)
            fileChanged=1
        return fileChanged
    except Exception, err :
        LOG.error("{}.".format(err))


def _convert_datetime(s):
    "converter which takes strings from ASC and converts to computable datetime objects"
    pt = s.rfind('.')
    micro_s = s[pt+1:]
    micro_s += '0'*(6-len(micro_s))
    #when = dt.datetime.strptime(s[:pt], '%Y-%m-%d %H:%M:%S').replace(microsecond = int(micro_s))
    when = datetime.strptime(s[:pt], '%Y-%m-%d %H:%M:%S').replace(microsecond = int(micro_s))
    return when

def _convert_isodatetime(s):
    "converter which takes strings from ASC and converts to computable datetime objects"
    pt = s.rfind('.')
    micro_s = s[pt+1:]
    micro_s += '0'*(6-len(micro_s))
    #when = dt.datetime.strptime(s[:pt], '%Y-%m-%d %H:%M:%S').replace(microsecond = int(micro_s))
    when = datetime.strptime(s[:pt], '%Y-%m-%dT%H:%M:%S').replace(microsecond = int(micro_s))
    return when


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


def _fuse(*exps):
    "fuse regular expressions together into a single or-expression"
    return '|'.join(r'(?:%s)' % x for x in exps)


def _get_granule_ID(IET_StartTime,IET_EndTime):
    """
    Calculates the deterministic granule ID. From...
    ADL/CMN/Utilities/INF/util/gran/src/InfUtil_GranuleID.cpp
    """
    # 
    NPP_GRANULE_ID_BASETIME = int(os.environ.get('NPP_GRANULE_ID_BASETIME', 1698019234000000))
    granuleSize = 85350000      # microseconds

    # Subtract the spacecraft base time from the arbitrary time to obtain
    # an elapsed time. 
    elapsedTime = IET_StartTime - NPP_GRANULE_ID_BASETIME

    # Divide the elapsed time by the granule size to obtain the granule number; 
    # the integer division will give the desired floor value.
    granuleNumber = np.floor(elapsedTime / granuleSize)
    #granuleNumber = np.ceil(elapsedTime / granuleSize)

    # Multiply the granule number by the granule size, then add the spacecraft
    # base time to obtain the granule start boundary time. Add the granule
    # size to the granule start boundary time to obtain the granule end
    # boundary time.
    startBoundary = (granuleNumber * granuleSize) + NPP_GRANULE_ID_BASETIME
    endBoundary = startBoundary + granuleSize

    # assign the granule start and end boundary times to the class members 
    granuleStartTime = startBoundary
    granuleEndTime = endBoundary
    
    # multiply the granule number by the granule size
    # then divide by 10^5 to convert the microseconds to tenths of a second; 
    # the integer division will give the desired floor value.
    timeCode = int((granuleNumber * granuleSize) / 100000)

    N_Granule_ID = 'NPP{:0>12d}'.format(timeCode)

    return N_Granule_ID,granuleStartTime,granuleEndTime


def _get_bounding_granule_ID(IET_StartTime):
    """
    Given the IET start time of a granule, calculates the 
    previous and next N_Granule_ID.
    """
    granuleSize = 85350000.0    # microseconds

    granuleStartTime = IET_StartTime
    granuleEndTime = granuleStartTime + granuleSize

    prevGranuleStartTime = granuleStartTime - granuleSize
    nextGranuleStartTime = granuleStartTime + granuleSize

    prevGranuleEndTime = prevGranuleStartTime + granuleSize
    nextGranuleEndTime = nextGranuleStartTime + granuleSize

    prev_N_Granule_ID = _get_granule_ID(prevGranuleStartTime,prevGranuleEndTime)[0]
    next_N_Granule_ID = _get_granule_ID(nextGranuleStartTime,nextGranuleEndTime)[0]

    return prev_N_Granule_ID,next_N_Granule_ID


def _create_dummy_sdr(inDir,requiredGeoShortname,requiredSdrShortname,crossGranules):
    """
    Examine the unpacked geolocation and radiometric blob/asc pairs for a DB pass, and 
    generate enough bounding dummy blob/asc pairs to ensure there are at least as many
    VIIRS EDR outputs as geolocation inputs, for each algorithm.
    """

    global sdrEndian 

    from adl_asc import PAT_URID, PAT_GRANULE_ID, PAT_GRANULE_VERSION, PAT_COLLECTION, \
                        PAT_RANGEDATETIME, PAT_SOURCE, PAT_BLOBPATH, \
                        PAT_EFFECTIVEDATETIME, PAT_OBSERVEDDATETIME

    PAT_CREATEDATETIME   = r'\("CreationDateTime" DATETIME EQ "(?P<CreationDateTime>[- \d:.]+)"'

    patternDict = {}
    patternDict['GRANULE_ID']        = PAT_GRANULE_ID
    patternDict['GRANULE_VERSION']   = PAT_GRANULE_VERSION
    patternDict['COLLECTION']        = PAT_COLLECTION
    patternDict['RANGEDATETIME']     = PAT_RANGEDATETIME
    patternDict['SOURCE']            = PAT_SOURCE
    patternDict['URID']              = PAT_URID
    patternDict['BLOBPATH']          = PAT_BLOBPATH
    patternDict['EFFECTIVEDATETIME'] = PAT_EFFECTIVEDATETIME
    patternDict['OBSERVEDDATETIME']  = PAT_OBSERVEDDATETIME
    patternDict['CREATEDATETIME']    = PAT_CREATEDATETIME

    DCONFIG=ADL_VARS['DCONFIG']

    patternRe = {}
    for key in patternDict.keys():
        patternRe[key] = re.compile(patternDict[key], re.M)

    # A key for sorting lists of granule dictionaries according to N_Granule_ID
    granIdKey = lambda x: (x['N_Granule_ID'])

    # The time format required for the asc file time strings
    ascFileTimeFormatStr = '%Y-%m-%d %H:%M:%S.%f'

    # Make a list of the required filetypes....
    requiredShortnames = requiredGeoShortname + requiredSdrShortname
    LOG.info("Required types : {}".format(requiredShortnames))

    geo_sdr_Dicts = {}
    missingShortNames = []

    for shortName in requiredShortnames :
        LOG.info("Searching for candidate {} granules...".format(shortName))
        granuleList = sorted(list(sift_metadata_for_viirs_sdr(shortName,crossGran=None,work_dir=inDir)))

        #for key in dict.keys():
        for dict in granuleList:
            for key in ['_filename','BlobPath']:
                try:
                    LOG.debug("{} granuleList = {}:{}".format(shortName,key,dict[key]))
                except KeyError:
                    LOG.debug("{} granuleList = {}:NULL".format(shortName,key))

        if granuleList :
            granule_IDs = []
            sorted_Dicts = []
            for dicts in sorted(granuleList,key=granIdKey) :
                granule_IDs.append(dicts['N_Granule_ID'])
                sorted_Dicts.append(dicts)
            LOG.info("real granules = {}".format(granule_IDs))
            geo_sdr_Dicts[shortName] = {'granule_IDs':granule_IDs,'sorted_Dicts':sorted_Dicts}
        else :
            LOG.warn("Missing blob/asc files for {}...".format(shortName))
            missingShortNames.append(shortName)

    if missingShortNames != []:
        LOG.error("Missing blob/asc files for shortnames {}, aborting...".format(missingShortNames))
        sys.exit(1)
        
    granuleSize_small = 83625544 # microseconds
    granuleSize       = 85350000 # microseconds
    granuleSize_large = 85404800 # microseconds

    # Make a list of dummy granule IDs
    dummy_granule_IDs = [] 
    dummy_granule_dict = {} 

    for shortName in requiredShortnames :

        # Generate the dummy granules for the first real granule
        try :
            firstDict = geo_sdr_Dicts[shortName]['sorted_Dicts'][0]
        except KeyError:
            LOG.error("Missing blob/asc files for {}, aborting...".format(shortName))
            sys.exit(1)

        first_N_Collection_Short_Name = firstDict['N_Collection_Short_Name']
        first_blobDir = path.dirname(firstDict['_filename'])
        first_URID = firstDict['URID']
        first_BlobFile = path.join(first_blobDir,"{}.{}".format(first_URID,first_N_Collection_Short_Name))

        # If this is the geolocation, get the IET start time of the first scan...
        if shortName == requiredGeoShortname[0]:
            first_XmlFile = "{}.xml".format(string.replace(first_N_Collection_Short_Name,'-','_'))
            first_XmlFile = path.join(ADL_HOME,'xml/VIIRS',first_XmlFile)
            first_BlobObj = adl_blob.map(first_XmlFile,first_BlobFile,writable=False,endian=sdrEndian)
            first_scanStartTime = getattr(first_BlobObj,'scanStartTime')[0]

        first_ascFile = firstDict['_filename']
        LOG.debug("The first file for {} is {}".format(shortName,first_ascFile))

        first_ObservedStartTime = firstDict['ObservedStartTime']
        first_ObservedEndTime = firstDict['ObservedEndTime']
        first_StartTime = firstDict['StartTime']
        first_EndTime = firstDict['EndTime']

        for crossGranIdx in range(1,crossGranules+1):

            scanStartTime = first_scanStartTime - crossGranIdx * granuleSize # microseconds

            N_Granule_ID = _get_granule_ID(scanStartTime,scanStartTime+granuleSize_small)[0]

            ObservedStartTime = first_ObservedStartTime - timedelta(microseconds=crossGranIdx*granuleSize_large)
            ObservedEndTime = first_ObservedEndTime - timedelta(microseconds=crossGranIdx*granuleSize_large)

            StartTime = first_StartTime - timedelta(microseconds=crossGranIdx*granuleSize)
            EndTime = first_EndTime - timedelta(microseconds=crossGranIdx*granuleSize)

            URID_dict = _getURID()
            URID = URID_dict['URID']
            creationDate_nousecStr = URID_dict['creationDate_nousecStr']
            creationDateStr = URID_dict['creationDateStr']

            newAscFileName = path.join(first_blobDir,"{}.asc".format(URID))

            # Populate dictionary for the previous granule
            granDict = {}
            granDict['URID']              = '("{}" UR "{}"\n'.format(URID,creationDate_nousecStr)
            granDict['GRANULE_ID']        = '  ("N_Granule_ID" STRING EQ "{}")\n'.format(N_Granule_ID)
            granDict['GRANULE_VERSION']   = '  ("N_Granule_Version" STRING EQ "A1")\n'
            granDict['COLLECTION']        = '  ("N_Collection_Short_Name" STRING EQ "{}")\n'.format(first_N_Collection_Short_Name)
            granDict['RANGEDATETIME']     = '  ("RangeDateTime" DATETIMERANGE EQ "{}" "{}")\n'.format(StartTime.strftime(ascFileTimeFormatStr),EndTime.strftime(ascFileTimeFormatStr))
            granDict['SOURCE']            = None #'  ("N_Dataset_Source" STRING EQ "nfts")\n'
            granDict['URID']              = '("{}" UR "{}"\n'.format(URID,creationDate_nousecStr)
            granDict['BLOBPATH']          = '  ("{}/{}.{}" FILE\n'.format(first_blobDir,URID,first_N_Collection_Short_Name)
            granDict['EFFECTIVEDATETIME'] = None
            granDict['OBSERVEDDATETIME']  = '  ("ObservedDateTime" DATETIMERANGE EQ "{}" "{}")\n'.format(ObservedStartTime.strftime(ascFileTimeFormatStr),ObservedEndTime.strftime(ascFileTimeFormatStr))
            granDict['CREATEDATETIME']  = '  ("CreationDateTime" DATETIME EQ "{}")\n'.format(creationDateStr)

            try:
                ascFile = open(first_ascFile,"r") # Open template file for reading
            except Exception, err :
                LOG.error("{}, aborting.".format(err))

            try:
                newAscFile = open(newAscFileName,"w") # Open template file for reading
            except Exception, err :
                LOG.error("{}, aborting.".format(err))

            LOG.debug("Copying {} to {}".format(first_ascFile,newAscFileName))

            for line in ascFile.readlines():
                for key in patternRe.keys():
                    m = patternRe[key].search(line)
                    if m:
                        if granDict[key] is not None:
                            line = line.replace(line,granDict[key])
                            LOG.debug("{} --> {}".format(key,line[:-1]))

                newAscFile.write(line) 

            ascFile.close()
            newAscFile.close()

            # Copy the first blob file...
            newBlobFileName = "{}.{}".format(URID,first_N_Collection_Short_Name)
            LOG.debug("Copying {} to {}".format(first_BlobFile,path.join(first_blobDir,newBlobFileName)))
            try:
                copyfile(first_BlobFile,path.join(first_blobDir,newBlobFileName))
            except IOError:
                LOG.warning("Blob file {} does not exist to be copied...".format(first_BlobFile))

            dummy_granule_IDs.append(N_Granule_ID)

            try :
                dummy_granule_dict[N_Granule_ID][shortName] = None
            except :
                dummy_granule_dict[N_Granule_ID] = {shortName:None}

            dummy_granule_dict[N_Granule_ID][shortName] = URID


    for shortName in requiredShortnames :

        # Generate the dummy granules for the last real granule
        try:
            lastDict = geo_sdr_Dicts[shortName]['sorted_Dicts'][-1]
        except KeyError:
            LOG.error("Missing blob/asc files for {}, aborting...".format(shortName))
            sys.exit(1)

        last_N_Collection_Short_Name = lastDict['N_Collection_Short_Name']
        last_blobDir = path.dirname(lastDict['_filename'])
        last_URID = lastDict['URID']
        last_BlobFile = path.join(last_blobDir,"{}.{}".format(last_URID,last_N_Collection_Short_Name))

        # If this is the geolocation, get the IET start time of the first scan...
        if shortName == requiredGeoShortname[0]:
            last_XmlFile = "{}.xml".format(string.replace(last_N_Collection_Short_Name,'-','_'))
            last_XmlFile = path.join(ADL_HOME,'xml/VIIRS',last_XmlFile)
            last_BlobObj = adl_blob.map(last_XmlFile,last_BlobFile,writable=False,endian=sdrEndian)
            last_scanStartTime = getattr(last_BlobObj,'scanStartTime')[0]

        last_ascFile = lastDict['_filename']
        LOG.debug("The last file for {} is {}".format(shortName,last_ascFile))

        last_ObservedStartTime = lastDict['ObservedStartTime']
        last_ObservedEndTime = lastDict['ObservedEndTime']
        last_StartTime = lastDict['StartTime']
        last_EndTime = lastDict['EndTime']

        for crossGranIdx in range(1,crossGranules+1):

            scanStartTime = last_scanStartTime + crossGranIdx * granuleSize # microseconds

            N_Granule_ID = _get_granule_ID(scanStartTime,scanStartTime+granuleSize_small)[0]

            ObservedStartTime = last_ObservedStartTime + timedelta(microseconds=crossGranIdx*granuleSize_large)
            ObservedEndTime = last_ObservedEndTime + timedelta(microseconds=crossGranIdx*granuleSize_large)

            StartTime = last_StartTime + timedelta(microseconds=crossGranIdx*granuleSize)
            EndTime = last_EndTime + timedelta(microseconds=crossGranIdx*granuleSize)

            URID_dict = _getURID()
            URID = URID_dict['URID']
            creationDate_nousecStr = URID_dict['creationDate_nousecStr']
            creationDateStr = URID_dict['creationDateStr']

            newAscFileName = path.join(last_blobDir,"{}.asc".format(URID))

            # Populate dictionary for the previous granule
            granDict = {}
            granDict['URID']              = '("{}" UR "{}"\n'.format(URID,creationDate_nousecStr)
            granDict['GRANULE_ID']        = '  ("N_Granule_ID" STRING EQ "{}")\n'.format(N_Granule_ID)
            granDict['GRANULE_VERSION']   = '  ("N_Granule_Version" STRING EQ "A1")\n'
            granDict['COLLECTION']        = '  ("N_Collection_Short_Name" STRING EQ "{}")\n'.format(last_N_Collection_Short_Name)
            granDict['RANGEDATETIME']     = '  ("RangeDateTime" DATETIMERANGE EQ "{}" "{}")\n'.format(StartTime.strftime(ascFileTimeFormatStr),EndTime.strftime(ascFileTimeFormatStr))
            granDict['SOURCE']            = None #'  ("N_Dataset_Source" STRING EQ "nfts")\n'
            granDict['URID']              = '("{}" UR "{}"\n'.format(URID,creationDate_nousecStr)
            granDict['BLOBPATH']          = '  ("{}/{}.{}" FILE\n'.format(last_blobDir,URID,last_N_Collection_Short_Name)
            granDict['EFFECTIVEDATETIME'] = None
            granDict['OBSERVEDDATETIME']  = '  ("ObservedDateTime" DATETIMERANGE EQ "{}" "{}")\n'.format(ObservedStartTime.strftime(ascFileTimeFormatStr),ObservedEndTime.strftime(ascFileTimeFormatStr))
            granDict['CREATEDATETIME']  = '  ("CreationDateTime" DATETIME EQ "{}")\n'.format(creationDateStr)

            try:
                ascFile = open(last_ascFile,"r") # Open template file for reading
            except Exception, err :
                LOG.error("{}, aborting.".format(err))

            try:
                newAscFile = open(newAscFileName,"w") # Open template file for reading
            except Exception, err :
                LOG.error("{}, aborting.".format(err))

            LOG.debug("Copying {} to {}".format(last_ascFile,newAscFileName))

            for line in ascFile.readlines():
                for key in patternRe.keys():
                    m = patternRe[key].search(line)
                    if m:
                        if granDict[key] is not None:
                            line = line.replace(line,granDict[key])
                            LOG.debug("{} --> {}".format(key,line[:-1]))

                newAscFile.write(line) 

            ascFile.close()
            newAscFile.close()

            # Copy the last blob file...
            newBlobFileName = "{}.{}".format(URID,last_N_Collection_Short_Name)
            LOG.debug("Copying {} to {}".format(last_BlobFile,path.join(last_blobDir,newBlobFileName)))
            try:
                copyfile(last_BlobFile,path.join(last_blobDir,newBlobFileName))
            except IOError:
                LOG.warning("Blob file {} does not exist to be copied...".format(last_BlobFile))

            dummy_granule_IDs.append(N_Granule_ID)

            try :
                dummy_granule_dict[N_Granule_ID][shortName] = None
            except :
                dummy_granule_dict[N_Granule_ID] = {shortName:None}

            dummy_granule_dict[N_Granule_ID][shortName] = URID


    # Make a unique list of the collected dummy granule IDs
    dummy_granule_dict['N_Granule_ID'] = list(set(dummy_granule_IDs))
    
    # Check whether granule IDs make sense
    dummy_granule_IDs = dummy_granule_dict['N_Granule_ID']
    dummy_granule_IDs.sort()
    LOG.info("dummy granules = {}".format(dummy_granule_IDs))

    all_granules = dummy_granule_IDs + granule_IDs
    all_granules.sort()
    
    for granIdx in range(len(all_granules[:-1])): 
        thisGranNum=all_granules[granIdx][6:]
        nextGranNum=all_granules[granIdx+1][6:]
        granNumDiff=int(nextGranNum)-int(thisGranNum)
        try:
            assert granNumDiff==853 or granNumDiff==854
            LOG.info("{} and {} are {} deciseconds apart.".format(thisGranNum,nextGranNum,granNumDiff))
        except AssertionError:
            LOG.warn("{} and {} are {} deciseconds apart. Should be 853 or 854...".format(thisGranNum,nextGranNum,granNumDiff))

    return dummy_granule_dict


def _granulate_ANC(inDir,geoDicts,algList,dummy_granule_dict):
    '''Granulates the input gridded blob files into the required ANC granulated datasets.'''

    import ANC
    import Algorithms
    global sdrEndian 
    global ancEndian 

    ANC_SCRIPTS_PATH = path.join(CSPP_RT_HOME,'viirs')
    ADL_ASC_TEMPLATES = path.join(ANC_SCRIPTS_PATH,'asc_templates')

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

    # If there are ancillary types to process, proceed...
    if collectionShortNames != []:

        # Download the required NCEP grib files
        LOG.info("Downloading NCEP GRIB ancillary into cache...")
        gribFiles = ANC.retrieve_NCEP_grib_files(geoDicts)
        LOG.debug('dynamic ancillary GRIB files: %s' % repr(gribFiles))
        if (gribFiles == []) :
            LOG.error('Failed to find or retrieve any GRIB files, aborting.')
            sys.exit(1)

        # Transcode the NCEP GRIB files into ADL NCEP-ANC-Int
        ncepGridBlobFiles = []
        for gribFile in gribFiles:
            gridBlobFile = ANC.create_NCEP_grid_blobs(gribFile)
            ncepGridBlobFiles.append(gridBlobFile)

        if (ncepGridBlobFiles == []) :
            LOG.error('Failed to convert NCEP GRIB  files to blobs, aborting.')
            sys.exit(1)

        LOG.debug("NCEP ncepGridBlobFiles = %r" % (ncepGridBlobFiles))

        # Open the NCEP gridded blob file

        ncepXmlFile = path.join(ADL_HOME,'xml/ANC/NCEP_ANC_Int.xml')
        endian = adl_blob.LITTLE_ENDIAN
        ncepBlobObjs = []

        for gridBlobFile in ncepGridBlobFiles :
            timeObj = gridBlobFile[0]
            ncepBlobFile = gridBlobFile[1]
            ncepBlobObj = adl_blob.map(ncepXmlFile, ncepBlobFile, endian=endian)
            ncepBlobObjs.append([timeObj,ncepBlobObj])
            LOG.debug("%s...\n%r" % (ncepBlobFile,[field for field in ncepBlobObj._fields]))
        
        # Get the NAAPS AOT if required...
        if 'AOT' in algList :
            # Download the required NCEP grib files
            #LOG.info("Downloading NAAPS GRIB ancillary into cache...")
            #gribFiles = ANC.retrieve_NAAPS_grib_files(geoDicts)
            #LOG.debug('Dynamic ancillary GRIB files: %s' % repr(gribFiles))
            #if (gribFiles == []) :
                #LOG.error('Failed to find or retrieve any GRIB files, aborting.')
                #sys.exit(1)

            gribFiles = glob(path.join(CSPP_RT_ANC_CACHE_DIR,'NAAPS-ANC-Int/NAAPS.*.grib2'))
            LOG.info("We are using for NAAPS grib2 :%r" %(gribFiles))

            # Transcode the NAAPS GRIB files into ADL NAAPS-ANC-Int
            naapsGridBlobFiles = []
            for gribFile in gribFiles:
                gridBlobFile = ANC.create_NAAPS_grid_blobs(gribFile)
                naapsGridBlobFiles.append(gridBlobFile)

            if (naapsGridBlobFiles == []) :
                LOG.error('Failed to convert NAAPS GRIB  files to blobs, aborting.')
                sys.exit(1)

            LOG.debug("NAAPS naapsGridBlobFiles = %r" % (naapsGridBlobFiles))

            naapsXmlFile = path.join(ADL_HOME,'xml/ANC/NAAPS_ANC_Int.xml')
            endian = adl_blob.LITTLE_ENDIAN
            naapsBlobObjs = []

            for gridBlobFile in naapsGridBlobFiles :
                timeObj = gridBlobFile[0]
                naapsBlobFile = gridBlobFile[1]
                naapsBlobObj = adl_blob.map(naapsXmlFile, naapsBlobFile, endian=endian)
                naapsBlobObjs.append([timeObj,naapsBlobObj])
                LOG.debug("%s...\n%r" % (naapsBlobFile,[field for field in naapsBlobObj._fields]))

    else:
        LOG.info("No ancillary required, so no NCEP or NAAPS gridding...")


    # Create a dict of ANC class instances, which will handle ingest and granulation
    ANC_objects = {}
    for shortName in collectionShortNames :

        className = ANC.classNames[shortName]
        ANC_objects[shortName] = getattr(ANC,className)(inDir=inDir,sdrEndian=sdrEndian,ancEndian=ancEndian)
        LOG.debug("ANC_objects[%s].blobDatasetName = %r" % (shortName,ANC_objects[shortName].blobDatasetName))
        
        # Just in case the same ANC class handles more than one collection short name
        if (np.shape(ANC_objects[shortName].collectionShortName) != () ):
            LOG.debug("    ANC_objects[%s].collectionShortName = %r" % (shortName,ANC_objects[shortName].collectionShortName))
            LOG.debug("    ANC_objects[%s].xmlName = %r" % (shortName,ANC_objects[shortName].xmlName))
            ANC_objects[shortName].collectionShortName = shortName
            ANC_objects[shortName].xmlName = ANC_objects[shortName].xmlName[shortName]
            LOG.debug("New ANC_objects[%s].collectionShortName = %r" % (shortName,ANC_objects[shortName].collectionShortName))
            LOG.debug("New ANC_objects[%s].xmlName = %r" % (shortName,ANC_objects[shortName].xmlName))

    # Ingest the ANC gridded data and copy to the gridData attribute of the ANC objects
    for shortName in collectionShortNames :
        LOG.info("Ingesting gridded ANC_objects: %s" % (shortName))
        if ANC_objects[shortName].sourceType == 'NCEP_ANC_Int' :
            ANC_objects[shortName].sourceList = [files[1] for files in ncepGridBlobFiles]
            ANC_objects[shortName].ingest(ancBlob=ncepBlobObjs)
        elif ANC_objects[shortName].sourceType == 'NAAPS_ANC_Int' :
            ANC_objects[shortName].sourceList = [files[1] for files in naapsGridBlobFiles]
            ANC_objects[shortName].ingest(ancBlob=naapsBlobObjs)
        else :
            ANC_objects[shortName].ingest()

    # Loop through the required ANC datasets and create the blobs.
    granIdKey = lambda x: (x['N_Granule_ID'])
    for dicts in sorted(geoDicts,key=granIdKey):
        for shortName in collectionShortNames :
        
            LOG.info("Processing dataset %s for %s" % (ANC_objects[shortName].blobDatasetName,shortName))

            # Set the geolocation information in this ancillary object for the current granule...
            ANC_objects[shortName].setGeolocationInfo(dicts)

            # Granulate the gridded data in this ancillary object for the current granule...
            ANC_objects[shortName].granulate(ANC_objects)

            # Shipout the granulated data in this ancillary object to a blob/asc pair.
            URID = ANC_objects[shortName].shipOutToFile()

            # If this granule ID is in the list of dummy IDs, add this URID to the 
            # dummy_granule_dict dictionary.
            N_Granule_ID = dicts['N_Granule_ID']
            if N_Granule_ID in dummy_granule_dict.keys():
                try :
                    dummy_granule_dict[N_Granule_ID][shortName] = None
                except :
                    dummy_granule_dict[N_Granule_ID] = {shortName:None}

                dummy_granule_dict[N_Granule_ID][shortName] = URID
                
    return dummy_granule_dict


def _granulate_GridIP(inDir,geoDicts,algList,dummy_granule_dict):
    '''Granulates the input gridded static data into the required GridIP granulated datasets.'''

    import GridIP
    import Algorithms
    global sdrEndian 
    global ancEndian 

    ANC_SCRIPTS_PATH = path.join(CSPP_RT_HOME,'viirs')
    ADL_ASC_TEMPLATES = path.join(ANC_SCRIPTS_PATH,'asc_templates')

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
        GridIP_objects[shortName] = getattr(GridIP,className)(inDir=inDir,sdrEndian=sdrEndian)
        LOG.debug("GridIP_objects[%s].blobDatasetName = %r" % (shortName,GridIP_objects[shortName].blobDatasetName))
        
        # Just in case the same GridIP class handles more than one collection short name
        if (np.shape(GridIP_objects[shortName].collectionShortName) != () ):
            LOG.debug("    GridIP_objects[%s].collectionShortName = %r" % (shortName,GridIP_objects[shortName].collectionShortName))
            LOG.debug("    GridIP_objects[%s].xmlName = %r" % (shortName,GridIP_objects[shortName].xmlName))
            GridIP_objects[shortName].collectionShortName = shortName
            GridIP_objects[shortName].xmlName = GridIP_objects[shortName].xmlName[shortName]
            LOG.debug("New GridIP_objects[%s].collectionShortName = %r" % (shortName,GridIP_objects[shortName].collectionShortName))
            LOG.debug("New GridIP_objects[%s].xmlName = %r" % (shortName,GridIP_objects[shortName].xmlName))

    # Loop through the required GridIP datasets and create the blobs.
    granIdKey = lambda x: (x['N_Granule_ID'])
    for dicts in sorted(geoDicts,key=granIdKey):
        for shortName in collectionShortNames :
        
            LOG.info("Processing dataset %s for %s" % (GridIP_objects[shortName].blobDatasetName,shortName))

            # Set the geolocation information in this ancillary object for the current granule...
            GridIP_objects[shortName].setGeolocationInfo(dicts)

            # Subset the gridded data for this ancillary object to cover the required lat/lon range.
            GridIP_objects[shortName].subset()

            # Granulate the gridded data in this ancillary object for the current granule...
            GridIP_objects[shortName].granulate(GridIP_objects)

            # Shipout the granulated data in this ancillary object to a blob/asc pair.
            URID = GridIP_objects[shortName].shipOutToFile()

            # If this granule ID is in the list of dummy IDs, add this URID to the 
            # dummy_granule_dict dictionary.
            N_Granule_ID = dicts['N_Granule_ID']
            if N_Granule_ID in dummy_granule_dict.keys():
                try :
                    dummy_granule_dict[N_Granule_ID][shortName] = None
                except :
                    dummy_granule_dict[N_Granule_ID] = {shortName:None}

                dummy_granule_dict[N_Granule_ID][shortName] = URID
                
    return dummy_granule_dict


def __cleanup_dummy_files(work_dir, algList, noDummyGranules, dummy_granule_dict):
    '''
    Remove radiometric, geolocation and ancillary blob/asc pairs, and product HDF5
    files that correspond to the dummy values of N_Granule_ID.
    '''

    # Remove dummy SDR and ancillary files
    if not noDummyGranules:
        LOG.info("Removing dummy SDR and ancillary blob/asc file pairs...")
        for granID in dummy_granule_dict['N_Granule_ID']:
            for shortName in dummy_granule_dict[granID].keys():
                URID = dummy_granule_dict[granID][shortName]
                dummyGlob = "{}.*".format(URID)
                dummyFiles = glob(path.join(work_dir,dummyGlob))
                for files in dummyFiles:
                    if path.exists(files):
                        LOG.debug('Removing dummy {:15}:{:16} file -> {}'.format(granID,shortName,files))
                        os.unlink(files)


def __cleanup_dummy_h5_files(work_dir, algList, dummy_granule_dict):
    '''
    Remove product HDF5 files that correspond to the dummy values of N_Granule_ID.
    '''
    # Remove dummy HDF5 product files
    for alg in algList :
        LOG.info("Removing dummy {} HDF5 files...".format(alg))
        for prefix,nodeName in zip(Algorithms.edr_hdf5_prefix[alg],Algorithms.edr_hdf5_Gran_0[alg]):
            edr_glob = "{}*.h5".format(prefix)
            edr_glob = path.join(work_dir,edr_glob)
            edr_hdf5_files = glob(edr_glob)
            if edr_hdf5_files != []:
                for hdf5File in edr_hdf5_files:
                    try : 
                        hdf5Obj = pytables.openFile(hdf5File)
                        Gran_0 = hdf5Obj.getNode(nodeName)
                        thisGranID =  getattr(Gran_0.attrs,'N_Granule_ID')[0][0]
                        hdf5Obj.close()

                        try:
                            if thisGranID in dummy_granule_dict['N_Granule_ID']:
                                LOG.info('Removing dummy {1:5} file with granule ID {0:15}: {2:}'.format(thisGranID,prefix,hdf5File))
                                os.unlink(hdf5File)
                        except Exception, err :
                            LOG.warn("Problem removing HDF5 file: {}".format(hdf5File))
                            LOG.warn("{}".format(err))
                            LOG.debug(traceback.format_exc())

                    except Exception, err :
                        hdf5Obj.close()
                        LOG.warn("Problem retrieving granule ID for file {}: {}".format(hdf5File,err))
                        LOG.debug(traceback.format_exc())


def __cleanup(work_dir, dirs_to_remove):
    '''
    Remove radiometric, geolocation and ancillary blob/asc pairs, and product HDF5
    files that correspond to the dummy values of N_Granule_ID.
    '''
    # Remove SDR asc/blob file pairs
    LOG.info("Removing SDR blob/asc file pairs...")
    sdr_glob = path.join(work_dir,"*.VIIRS-[MI][1-9]*-SDR")
    geo_glob = path.join(work_dir,"*.VIIRS-[MI]*-GEO*")
    blobFiles = glob(sdr_glob) + glob(geo_glob)
    if blobFiles != [] :
        for blobFile in blobFiles:
            blobDir = path.dirname(blobFile)
            URID = string.split(path.basename(blobFile),".")[0]
            ascFile = path.join(blobDir,"{}.asc".format(URID))
            try:
                LOG.debug('Removing {}'.format(blobFile))
                os.unlink(blobFile)
            except Exception, err:
                LOG.warn( "{}".format(str(err)))
            try:
                LOG.debug('Removing {}'.format(ascFile))
                os.unlink(ascFile)
            except Exception, err:
                LOG.warn( "{}".format(str(err)))

    # Remove ANC and GridIP asc/blob file pairs
    LOG.info("Removing ANC and GridIP blob/asc file pairs...")
    anc_glob = path.join(work_dir,"*.VIIRS-ANC*")
    gridIP_glob = path.join(work_dir,"*.VIIRS-GridIP*")
    blobFiles = glob(anc_glob) + glob(gridIP_glob)
    if blobFiles != [] :
        for blobFile in blobFiles:
            blobDir = path.dirname(blobFile)
            URID = string.split(path.basename(blobFile),".")[0]
            ascFile = path.join(blobDir,"{}.asc".format(URID))
            try:
                LOG.debug('Removing {}'.format(blobFile))
                os.unlink(blobFile)
            except Exception, err:
                LOG.warn( "{}".format(str(err)))
            try:
                LOG.debug('Removing {}'.format(ascFile))
                os.unlink(ascFile)
            except Exception, err:
                LOG.warn( "{}".format(str(err)))

    # Remove all other asc/blob pairs (usually products).
    LOG.info("Remove all other asc/blob pairs (usually products)...")
    ascBlobFiles = glob(path.join(work_dir, '????????-?????-????????-????????.*'))
    if ascBlobFiles != [] :
        for ascBlobFile in ascBlobFiles:
            try:
                LOG.debug('Removing {}'.format(ascBlobFile))
                os.unlink(ascBlobFile)
            except Exception, err:
                LOG.warn( "{}".format(str(err)))

    # Remove log directory
    LOG.info("Removing other directories ...")
    for dirname in dirs_to_remove:
        fullDirName = path.join(work_dir,dirname)
        LOG.debug('Removing {}'.format(fullDirName))
        try :
            rmtree(fullDirName, ignore_errors=False)
        except Exception, err:
            LOG.warn( "{}".format(str(err)))


def _argparse():
    '''
    Method to encapsulate the option parsing and various setup tasks.
    '''

    import argparse

    endianChoices = ['little','big']
    algorithmChoices = ['VCM','AOT','SST','SRFREF','VI','ST','LST',
            'ATMOS','LAND','OCEAN','ALL','MPC']

    defaults = {'work_dir':'.',
                'skipSdrUnpack':False,
                'skipAuxLinking':False,
                'skipAncillary':False,
                'skipAlgorithm':False,
                'cspp_debug':False,
                'noAlgChain':False,
                'noDummyGranules':False,
                'processors':1,
                'sdr_Endianness':'little',
                'anc_Endianness':'little',
                'compress':False,
                'aggregate':False
                }

    description = '''Run one or more ADL VIIRS EDR Controllers.'''
    usage = "usage: %prog [mandatory args] [options]"
    version = __version__

    parser = argparse.ArgumentParser(
                                     #usage=usage,
                                     #version=version,
                                     description=description
                                     )

    # Mandatory arguments

    parser.add_argument('-i','--input_files',
                      action="store",
                      dest="inputFiles",
                      type=str,
                      required=True,
                      help='''The fully qualified path to the input files. May be
                            a directory or a file glob.'''
                      )

    parser.add_argument('--alg',
                      action="store",
                      dest="algorithm",
                      type=str,
                      required=True,
                      choices=algorithmChoices,
                      help='''The VIIRS algorithm to run.\n\n
                              Possible values are...
                              {}.
                           '''.format(algorithmChoices.__str__()[1:-1])
                      )

    # Optional arguments 

    parser.add_argument('-w','--work_directory',
                      action="store",
                      dest="work_dir",
                      type=str,
                      default=defaults['work_dir'],
                      help='''The directory which all activity will occur in, 
                      defaults to the current directory. 
                      [default: {}]'''.format(defaults['work_dir'])
                      )
    
    parser.add_argument('--skip_sdr_unpack',
                      action="store_true",
                      dest="skipSdrUnpack",
                      default=defaults['skipSdrUnpack'],
                      help='''Skip the unpacking of the VIIRS SDR HDF5 files. 
                      [default: {}]'''.format(defaults['skipSdrUnpack'])
                      )

    parser.add_argument('--skip_aux_linking',
                      action="store_true",
                      dest="skipAuxLinking",
                      default=defaults['skipAuxLinking'],
                      help='''Skip the the linking to auxillary files. 
                      [default: {}]'''.format(defaults['skipAuxLinking'])
                      )

    parser.add_argument('--skip_ancillary',
                      action="store_true",
                      dest="skipAncillary",
                      default=defaults['skipAncillary'],
                      help='''Skip the retrieval and granulation of ancillary 
                              data. 
                              [default: {}]'''.format(defaults['skipAncillary'])
                      )

    parser.add_argument('--skip_algorithm',
                      action="store_true",
                      dest="skipAlgorithm",
                      default=defaults['skipAlgorithm'],
                      help='''Skip running the VIIRS EDR algorithm(s). 
                      [default: {}]'''.format(defaults['skipAlgorithm'])
                      )

    parser.add_argument('--debug',
                      action="store_true",
                      dest="cspp_debug",
                      default=defaults['cspp_debug'],
                      help='''Enable debug mode on ADL and avoid cleaning 
                      workspace. 
                      [default: {}]'''.format(defaults['cspp_debug'])
                      )

    parser.add_argument('--no_chain',
                      action="store_true",
                      dest="noAlgChain",
                      default=defaults['noAlgChain'],
                      help='''Do not run prerequisite algorithms. 
                      [default: {}]'''.format(defaults['noAlgChain'])
                      )

    parser.add_argument('--no_dummy_granules',
                      action="store_true",
                      dest="noDummyGranules",
                      default=defaults['noDummyGranules'],
                      help='''Do not generate dummy SDR or ancillary 
                      cross granules. 
                      [default: {}]'''.format(defaults['noDummyGranules'])
                      )

    parser.add_argument('-p','--processors',
                      action="store",
                      dest="processors",
                      default=defaults['processors'],
                      type=int,
                      help='''Number of cpus to use for granule processing. 
                      [default: {}]'''.format(defaults['processors'])
                      )

    parser.add_argument('--sdr_endianness',
                      action="store",
                      dest="sdr_Endianness",
                      type=str,
                      default=defaults['sdr_Endianness'],
                      choices=endianChoices,
                      help='''The input VIIRS SDR endianness.\n\n
                              Possible values are...
                              {}. [default: '{}']
                           '''.format(endianChoices.__str__()[1:-1],
                               defaults['sdr_Endianness'])
                      )
    
    parser.add_argument('--anc_endianness',
                      action="store",
                      dest="anc_Endianness",
                      type=str,
                      default=defaults['anc_Endianness'],
                      choices=endianChoices,
                      help='''The input VIIRS ancillary endianness.\n\n
                              Possible values are...
                              {}. [default: '{}']
                           '''.format(endianChoices.__str__()[1:-1],
                               defaults['anc_Endianness'])
                      )

    parser.add_argument('-z', '--zip',
                      action="store_true",
                      dest="compress",
                      default=defaults['compress'],
                      help='''Enable product compression. 
                      [default: {}]'''.format(defaults['compress'])
                      )

    parser.add_argument('-a', '--aggregate',
                      action="store_true",
                      dest="aggregate",
                      default=defaults['aggregate'],
                      help='''Enable product nagg aggregation. 
                      [default: {}]'''.format(defaults['aggregate'])
                      )

    parser.add_argument('-v', '--verbose',
                      dest='verbosity',
                      action="count",
                      default=0,
                      help='''Each occurrence increases 
                      verbosity 1 level from INFO: -v=DEBUG'''
                      )

    #parser.add_argument('-V','--version',
                      #action='version',
                      #version='''{}\n
                                 #%(prog)s (myprog version 0.1)'''.format(version)
                      ##version='''%(prog)s (myprog version 0.1)'''
                      #)

    args = parser.parse_args()


    # Set the work directory
    work_dir = check_and_convert_path("WORK_DIR",args.work_dir)
    LOG.debug('Setting the work directory to %r' % (work_dir))

    # Set up the logging
    d = datetime.now()
    timestamp = d.isoformat()
    logname= "viirs_edr."+timestamp+".log"
    logfile= path.join(work_dir, logname )

    levels = [logging.ERROR, logging.WARN, logging.INFO, logging.DEBUG]
    level = levels[min(args.verbosity,3)]
    configure_logging(level,FILE=logfile)
    
    # create work directory
    if not path.isdir(work_dir):
        LOG.info('creating directory %s' % (work_dir))
        os.makedirs(work_dir)
    log_dir = path.join(work_dir, 'log')
    if not path.isdir(log_dir):
        LOG.debug('creating directory %s' % (log_dir))
        os.makedirs(log_dir)

    return args,work_dir ,log_dir



def main():

    options,work_dir,log_dir = _argparse()
    dpe_version = check_existing_env_var( 'DPE_VER',default_value=None )
    LOG.info("product version: DPE_VER for VIIRS EDR  %s"%dpe_version)

    LOG.info("Starting CSPP VIIRS EDR...")

    # Determine the ordered list of algs to satisfy the required 
    # algorithm's dependencies.
    algList = [options.algorithm]

    if not options.noAlgChain : 
        LOG.info("We ARE chaining algorithms...")
        thisAlg = algList[0]
        LOG.debug("thisAlg: %r" % (thisAlg))
        thisAlgPreReqs = copy.copy(Algorithms.prerequisites[thisAlg])
        LOG.debug("thisAlgPreReqs: %r" % (thisAlgPreReqs))
        thisAlgPreReqs.append(thisAlg)
        LOG.debug("thisAlgPreReqs: %r" % (thisAlgPreReqs))
        algList = thisAlgPreReqs
        LOG.debug("algList: %r" % (algList))
    else :
        LOG.info("We are NOT chaining algorithms...")

    LOG.info("List of required algorithm keys: {}".format(algList))

    # Create a list of algorithm class references
    algorithms = []
    for alg in algList :
        if alg in Algorithms.meta_algorithms:
            LOG.info("Skipping meta algorithm '{}' from list of algorithms pointers.".format(alg)) 
        else :
            algName = Algorithms.modules[alg]
            algorithms.append(getattr(Algorithms,algName))

    LOG.info("List of required algorithm class references: {}".format(algorithms))

    # Determine the number of VIIRS SDR cross granules required to 
    # process the algorithm chain.
    cumulativeCrossGranules = _get_alg_cross_granules(algList,options.noAlgChain)
    LOG.info("We require {} cross granules for {}".format(cumulativeCrossGranules[options.algorithm],\
                                                          options.algorithm))

    # Removing meta-algs from algList
    for metaAlg in Algorithms.meta_algorithms:
        try:
            algList.pop(algList.index(metaAlg))
            LOG.info("Removed meta algorithm '{}' from list of algorithms.".format(metaAlg)) 
        except:
            pass

    LOG.info("meta-updated list of required algorithm keys: {}".format(algList))

    # Determine what geolocation types are required for each algorithm, and 
    # get the required blob/asc pairs
    requiredGeoShortname,requiredGeoPrefix = _get_geo_prefixes(algorithms)
    geo_unpacking_problems = _get_sdr_blobs(options.inputFiles,work_dir,
            options.skipSdrUnpack,requiredGeoShortname,requiredGeoPrefix)

    # Determine what radiometric types are required for each algorithm, and 
    # get the required blob/asc pairs
    requiredSdrShortname,requiredSdrPrefix = _get_radio_prefixes(algorithms)
    radio_unpacking_problems = _get_sdr_blobs(options.inputFiles,work_dir,
            options.skipSdrUnpack,requiredSdrShortname,requiredSdrPrefix)

    unpacking_problems = geo_unpacking_problems + radio_unpacking_problems

    # Check geolocation SDR metadata for wrong N_Granule_Version, and fix...
    geo_blob_names = glob(path.join(work_dir,"*.*GEO-TC"))
    for geo_blob_name in geo_blob_names:
        ascName = corresponding_asc_path(geo_blob_name)
        fileChanged = _strReplace(ascName,
                '("N_Granule_Version" STRING EQ "A2")',
                '("N_Granule_Version" STRING EQ "A1")')
        if fileChanged :
            LOG.info("Fixed N_Granule_Version in %s metadata" 
                    % (path.basename(geo_blob_name)))


    # Check radiometric SDR metadata for wrong N_Granule_Version, and fix...
    sdr_blob_names = glob(path.join(work_dir,"*.*SDR"))
    for sdr_blob_name in sdr_blob_names:
        ascName = corresponding_asc_path(sdr_blob_name)
        fileChanged = _strReplace(ascName,
                '("N_Granule_Version" STRING EQ "A2")',
                '("N_Granule_Version" STRING EQ "A1")')
        if fileChanged :
            LOG.info("Fixed N_Granule_Version in %s metadata" 
                    % (path.basename(sdr_blob_name)))


    # Set the VIIRS SDR endianness from the input option...
    global sdrEndian
    set_sdr_endian(options.sdr_Endianness)

    # Create any required dummy geolocation and radiometric granules
    dummy_granule_dict = {}

    if not options.noDummyGranules:
        LOG.info("Creating dummy geolocation and radiometric granules...")
        dummy_granule_dict = _create_dummy_sdr(work_dir,requiredGeoShortname,
                requiredSdrShortname,cumulativeCrossGranules[options.algorithm])

        requiredShortnames = requiredGeoShortname + requiredSdrShortname

        LOG.info("Dummy granule IDs = {}".format(dummy_granule_dict['N_Granule_ID']))

        for granID in dummy_granule_dict['N_Granule_ID']:
            for shortName in requiredShortnames:
                LOG.debug("dummy_granule_dict[{}][{}] = {}".format(\
                        granID,shortName,dummy_granule_dict[granID][shortName]))


    # Read through ascii metadata and build up information table
    LOG.info('Sifting through geolocation metadata to find VIIRS SDR processing candidates...')
    anc_granules_to_process = None

    # Determine the candidate geolocation granules for which to generate ancillary.
    # Hint: ... Make ALL THE ANCILLARY! :-)
    for geoType in requiredGeoShortname :
        LOG.info("Searching for candidate %s geolocation granules for VIIRS ancillary..." % (geoType))
        anc_granules_to_process = sorted(list(sift_metadata_for_viirs_sdr(geoType,crossGran=None,work_dir=work_dir)))
        if anc_granules_to_process :
            break
        else :
            LOG.info("\tNo %s geolocation granules for VIIRS ancillary" % (geoType))

    # Determine the candidate geolocation granules for which we can generate VIIRS products.
    for alg in algorithms :
        for geoType in requiredGeoShortname :
            LOG.info("Searching for candidate %s geolocation granules for VIIRS %s ..." % (geoType,alg.AlgorithmName))
            crossGranules = cumulativeCrossGranules[alg.AlgorithmString]
            LOG.info("We need {} cross granules for {}".format(crossGranules,
                alg.AlgorithmString))
            alg.granules_to_process = sorted(list(sift_metadata_for_viirs_sdr(geoType, \
                    crossGran=crossGranules,work_dir=work_dir)))
            if alg.granules_to_process :
                break
            else :
                LOG.info("\tNo %s geolocation granule groups of length %d for VIIRS %s" % (geoType,
                    (2*crossGranules+1),alg.AlgorithmName))
        if not alg.granules_to_process :
            pass
        
    # A key for sorting lists of granule dictionaries according to N_Granule_ID
    granIdKey = lambda x: (x['N_Granule_ID'])

    # List ancillary processing candidate granule IDs
    if anc_granules_to_process :
        LOG.info("We have {} candidate granules for ancillary.".format(len(anc_granules_to_process)))
        LOG.info("{:^19}{:^30}{:^30}{:^20}{:^30}{:^30}{:^20}".format('N_Granule_ID',\
                                                         'ObservedStartTime',\
                                                         'ObservedEndTime',\
                                                         'Observed Duration',\
                                                         'StartTime',\
                                                         'EndTime',\
                                                         'Duration'
                                                         ))
        for dicts in sorted(anc_granules_to_process,key=granIdKey) :
            LOG.info("{:^19}{:^30}{:^30}{:^20}{:^30}{:^30}{:^20}".format(dicts['N_Granule_ID'], \
                                                      dicts['ObservedStartTime'].isoformat(), \
                                                      dicts['ObservedEndTime'].isoformat(), \
                                                      (dicts['ObservedEndTime']-dicts['ObservedStartTime']), \
                                                      dicts['StartTime'].isoformat(), \
                                                      dicts['EndTime'].isoformat(), \
                                                      (dicts['EndTime']-dicts['StartTime'])
                                                      ))

    # List algorithm processing candidate granule IDs
    for alg in algorithms :
        if alg.granules_to_process :
            LOG.info("")
            LOG.info("We have {} candidate granules for {}".format(len(alg.granules_to_process),alg.AlgorithmName))
            LOG.info("{:^19}{:^30}{:^30}{:^30}{:^30}{:^51}".format('N_Granule_ID','ObservedStartTime','ObservedEndTime','StartTime','EndTime','(prev,curr,next) N_Granule_ID'))
            for dicts in sorted(alg.granules_to_process,key=granIdKey) :
                LOG.info("{:^19}{:^30}{:^30}{:^30}{:^30}({:^17}{:^17}{:^17})".format(dicts['N_Granule_ID'], \
                                                          dicts['ObservedStartTime'].isoformat(), \
                                                          dicts['ObservedEndTime'].isoformat(), \
                                                          dicts['StartTime'].isoformat(), \
                                                          dicts['EndTime'].isoformat(), \
                                                          dicts['N_Granule_ID_prev'], \
                                                          dicts['N_Granule_ID'], \
                                                          dicts['N_Granule_ID_next']))
        else :
            LOG.info("We have %d candidate granules for %s" % (len(alg.granules_to_process),alg.AlgorithmName))
            num_xml_files_to_process = 0
            num_no_output_runs = 0
            noncritical_problem = False
            environment_error = False
            if not options.skipAlgorithm :
                return get_return_code(geo_unpacking_problems, num_xml_files_to_process, \
                        num_no_output_runs, noncritical_problem, environment_error)

    LOG.info("")
    LOG.info('Finished sifting through metadata to find VIIRS SDR processing candidates')

    
    # Set the VIIRS ancillary endianness from the input option...
    global ancEndian
    set_anc_endian(options.anc_Endianness)

    # Expand any user specifiers in the various paths
    LOG.info("CSPP_RT_ANC_PATH:       %s" % (CSPP_RT_ANC_PATH))
    LOG.info("CSPP_RT_ANC_CACHE_DIR:  %s" % (CSPP_RT_ANC_CACHE_DIR))
    #LOG.debug("CSPP_RT_ANC_TILE_PATH:  %s" % (CSPP_RT_ANC_TILE_PATH))
    #LOG.debug("LD_LIBRARY_PATH:     %s" % (LD_LIBRARY_PATH))
    #LOG.debug("DSTATICDATA:         %s" % (DSTATICDATA))

    # Retrieve and granulate the required ancillary data...
    if not options.skipAncillary :

        t1 = time()

        LOG.info('Retrieving and granulating ancillary data...')

        # Granulate the VIIRS ANC data
        dummy_granule_dict = _granulate_ANC(work_dir,anc_granules_to_process,algList,dummy_granule_dict)

        # Granulate the VIIRS GridIP data
        dummy_granule_dict = _granulate_GridIP(work_dir,anc_granules_to_process,algList,dummy_granule_dict)

        t2 = time()

        LOG.info( "Generation of VIIRS ancillary data took %f seconds." % (t2-t1))

    else :

        LOG.info('Skipping retrieval and granulation of ancillary data.')

    # List the SDR and ancillary dummy granules
    if not options.noDummyGranules:

        for granID in dummy_granule_dict['N_Granule_ID']:
            maxStrLen = len(max(dummy_granule_dict[granID].keys(), key=len))
            formatStr = "dummy_granule_dict[{{:15}}][{{:{0}}}] = {{:35}}".format(maxStrLen)
            for shortName in sorted(dummy_granule_dict[granID].keys()):
                LOG.info(formatStr.format(granID,shortName,dummy_granule_dict[granID][shortName]))



    # Link in auxillary files
    if not options.skipAuxLinking :
        LOG.info('Linking in the VIIRS EDR auxillary files...')

        ADL_VIIRS_ANC_GLOBS = []
        for alg in algorithms :
            for shortName in alg.AUX_collectionShortNames:
                aux_glob = string.replace(shortName,'-Int','')
                aux_glob = '*{}*'.format(aux_glob)
                ADL_VIIRS_ANC_GLOBS.append(aux_glob)
                LOG.info('Auxillary glob for {} = {}'.format(shortName,aux_glob))

        search_dirs = [CSPP_RT_ANC_CACHE_DIR if 
            CSPP_RT_ANC_CACHE_DIR is not None else anc_dir] \
                    + [CSPP_RT_ANC_PATH]

        LOG.info('Auxillary search directories = {}'.format(search_dirs))

        ancillary_files_neeeded = anc_files_needed(ADL_VIIRS_ANC_GLOBS, search_dirs)
        link_ancillary_to_work_dir(work_dir, ancillary_files_neeeded)

    else :
        LOG.info('Skipping linking in the VIIRS EDR auxillary files.')

    LOG.info("Number of processors is {}".format(options.processors))

    # Specify the algorithm we want to run via the Lw XML file.
    if not options.skipAlgorithm :

        for alg in algorithms :

            try :
                # build XML configuration files for jobs that can be run
                LOG.info("Building %s XML files for %d granules" % \
                        (alg.AlgorithmName,len(alg.granules_to_process)))

                #LOG.info("alg.granules_to_process: %r" %(alg.granules_to_process))

                # Generate the VIIRS LW xml files for each N_Granule_ID for this algorithm.
                xml_files_to_process = alg.generate_viirs_edr_xml(work_dir, alg.granules_to_process)
            
                LOG.info('%d granules to process: %s' % \
                        (len(xml_files_to_process), ''.join(name+' -> '+xmlfile+'\n' \
                        for (name,xmlfile) in xml_files_to_process)))


                LOG.info("Running VIIRS %s ..." % (alg.AlgorithmName))
                crashed_runs, no_output_runs, geo_problem_runs, bad_log_runs = \
                    alg.run_xml_files(work_dir, \
                                      xml_files_to_process, \
                                      nprocs = options.processors, \
                                      WORK_DIR = work_dir, \
                                      ADL_HOME = ADL_HOME, \
                                      COMPRESS=options.compress, \
                                      AGGREGATE=options.aggregate)

                LOG.debug('crashed_runs : {}'.format(crashed_runs))
                LOG.debug('no_output_runs : {}'.format(no_output_runs))
                LOG.debug('geo_problem_runs : {}'.format(geo_problem_runs))
                LOG.debug('bad_log_runs : {}'.format(bad_log_runs))

                ## considered a noncritical problem if there were any crashed runs, runs that produced no output,
                ## runs where Geo failed, or runs where ADL logs indicated a problem
                noncritical_problem = crashed_runs or no_output_runs or geo_problem_runs or bad_log_runs

                ## print final disposition message and get return code
                environment_error=False
                rc = get_return_code(unpacking_problems, len(xml_files_to_process), \
                        len(no_output_runs), noncritical_problem, environment_error)

                # If this alg failed, return error code and exit, preserving inputs and log files
                if rc != 0 :
                    LOG.warn("Non-zero error code %d for %s." % (rc, alg.AlgorithmName))

            except Exception:
                LOG.error(traceback.format_exc())
                rc = 1

            # if no errors or only non-critical errors: clean up
            LOG.info("{} return code : {}".format(alg.AlgorithmName,rc))


        if rc == 0 and not options.cspp_debug:
            for alg in algorithms :

                LOG.info("Cleaning up workspace for VIIRS {}...".format(alg.AlgorithmName))

                algorithmXmlGlob = '%s*.xml' % (alg.algorithmLWxml)
                algorithmLogGlob = '%s_*' % (alg.controllerName)

                alg.cleanup(work_dir, algorithmXmlGlob, algorithmLogGlob)
    
    else :

        LOG.info("Skipping execution of VIIRS %s ..." % (alg.AlgorithmName))


    # Aggregation...
    if options.aggregate is True:

        __cleanup_dummy_h5_files(work_dir, algList, dummy_granule_dict)
        
        for alg in algorithms :
            if alg.AlgorithmString in algList:
                LOG.info("Aggregating VIIRS {}...".format(alg.AlgorithmName))
                number_problems = aggregate_products(work_dir, alg.EDR_collectionShortNames)

    # General Cleanup...
    if not options.cspp_debug:

        # Remove dummy asc/blob pairs and HDF5 files
        __cleanup_dummy_files(work_dir, algList, options.noDummyGranules, dummy_granule_dict)
        
        # Remove dummy HDF5 files
        if options.aggregate is False:
            __cleanup_dummy_h5_files(work_dir, algList, dummy_granule_dict)

        # Remove what's left...
        __cleanup(work_dir, [log_dir])

    try :
        return rc
    except :
        return 0



if __name__=='__main__':
    sys.exit(main())  
