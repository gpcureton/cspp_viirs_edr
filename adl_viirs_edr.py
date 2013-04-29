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
      Thus, for N input SDR granules, you may get up to N-2 or less output EDR granules if 
      they are all contiguous, depending on the algorithm.
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

    python adl_viirs_edr.py  --input_files=INPUTFILES

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
from adl_asc import skim_dir, contiguous_granule_groups, granule_groups_contain, effective_anc_contains,_eliminate_duplicates,_is_contiguous, corresponding_asc_path, RDR_REQUIRED_KEYS, POLARWANDER_REQUIRED_KEYS

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

import Algorithms

# we're not likely to succeed in processing using geolocation smaller than this many bytes
MINIMUM_SDR_BLOB_SIZE = 81000000

# maximum delay between granule end time and next granule start time to consider them contiguous
MAX_CONTIGUOUS_DELTA=timedelta(seconds = 2)

CSPP_RT_ANC_HOME = path.abspath(os.getenv('CSPP_RT_ANC_HOME'))



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
        LOG.error('Invalid value for the VIIRS SDR endianness : %s ' % (inputEndianness))


def set_anc_endian(inputEndianness) :
    '''
    Set the global ancillary endianness variable.
    '''
    global ancEndian 
    if inputEndianness=='big' :
        ancEndian = adl_blob.BIG_ENDIAN
    elif inputEndianness=='little' :
        ancEndian = adl_blob.LITTLE_ENDIAN
    else :
        LOG.error('Invalid value for the VIIRS ancillary endianness : %s ' % (inputEndianness))


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

    if ((input_files is None) or (input_files == "*")):
        # Input file glob is of form "/path/to/files" or "/path/to/files/*"
        inputGlobs['GEO'] = 'GMTCO_npp*.h5'
        inputGlobs['MOD'] = 'SVM*_npp*.h5'
        inputGlobs['IMG'] = 'SVI*_npp*.h5'
    elif ((('GMTCO' in input_files) or ('SVM' in input_files) or ('SVI' in input_files)) and ('*' in input_files)) :
        # Input file glob is of form "/path/to/files/GMTCO*" or "/path/to/files/SVI*" or "/path/to/files/SVM*" 
        fileGlob = string.rstrip(string.lstrip(input_files,charsToKill),charsToKill)
        LOG.debug("fileGlob = %s" %(fileGlob))
        inputGlobs['GEO'] = "GMTCO%s.h5" %(fileGlob)
        inputGlobs['MOD'] = "SVM*%s.h5" %(fileGlob)
        inputGlobs['IMG'] = "SVI*%s.h5" %(fileGlob)
        for fileType in ['GEO','MOD','IMG']:
            inputGlobs[fileType] = string.replace(inputGlobs[fileType],"**","*")
    elif path.isfile(input_path) :
        # Input file glob is of form "/path/to/files/GMTCO_npp_d_t_e_b_c_cspp_dev.h5" 
        fileGlob = string.rstrip(string.lstrip(string.split(input_files,"b")[0],charsToKill),charsToKill)
        LOG.debug("fileGlob = %s" %(fileGlob))
        inputGlobs['GEO'] = "GMTCO%s*.h5" %(fileGlob)
        inputGlobs['MOD'] = "SVM*%s*.h5" %(fileGlob)
        inputGlobs['IMG'] = "SVI*%s*.h5" %(fileGlob)
        for fileType in ['GEO','MOD','IMG']:
            inputGlobs[fileType] = string.replace(inputGlobs[fileType],"**","*")

    return input_dir,inputGlobs


def _unpack_sdr(work_dir,input_dir,inputGlob):
    '''
    Unpack HDF5 VIIRS SDRs from the input directory to the work directory
    '''
    unpacking_problems = 0
    h5_names = glob(path.join(input_dir,inputGlob))

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
    #
    # If there is only one granule in the list, yield it immediately
    seq = {}
    if len(granset)==1 :

        a = granlist[0]
        seq[a['URID']] = a
        LOG.info('contiguous sequence has %d granules' % (len(seq)))
        yield tuple(sorted(seq.values(), key=start_time_key))
        seq.clear()

    else :

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
            LOG.info('Leftover contiguous sequence has %d granules' % (len(seq)))
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

    LOG.debug('geoGroupList : %r'%(geoGroupList))

    if len(geoGroupList)==0:
        LOG.debug('No geoGroupList found...')
        return

    # Loop through the contigous granule groups 
    for group in _contiguous_granule_groups(skim_dir(work_dir, N_Collection_Short_Name=collectionShortName)):
        ##- for VIIRS, we can process everything but the first and last granule
        ##- for CrIS, use [4:-4]
        LOG.info('Contiguous granule group of length: %r' % (len(group),))

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


def _strReplace(fileName,oldString,newString):
    """
    Check fileName for occurences of oldString, if found fileName is opened and oldString is 
    replaced with newString.
    """
    fileChanged=0
    with open(fileName) as thefile:
        content = thefile.read()                 # read entire file into memory
        replacedText = content.replace(oldString, newString)
    if replacedText != content:
        LOG.debug('Replacing occurence of "%s" in %s with "%s"' % (oldString,path.basename(fileName),newString))
        with open(fileName, 'w') as thefile:
            thefile.write(replacedText)
        fileChanged=1
    return fileChanged


def _granulate_ANC(inDir,geoDicts,algList):
    '''Granulates the input gridded blob files into the required ANC granulated datasets.'''

    import ANC
    import Algorithms
    global sdrEndian 
    global ancEndian 

    CSPP_RT_HOME = os.getenv('CSPP_RT_HOME')
    ANC_SCRIPTS_PATH = path.join(CSPP_RT_HOME,'viirs')
    CSPP_RT_ANC_CACHE_DIR = os.getenv('CSPP_RT_ANC_CACHE_DIR')
    
    ADL_ASC_TEMPLATES = path.join(ANC_SCRIPTS_PATH,'asc_templates')

    # Download the required NCEP grib files
    LOG.info("Downloading NCEP GRIB ancillary into cache...")
    gribFiles = ANC.retrieve_NCEP_grib_files(geoDicts)
    LOG.debug('dynamic ancillary GRIB files: %s' % repr(gribFiles))
    if (gribFiles == []) :
        LOG.error('Failed to find or retrieve any GRIB files, aborting.')
        sys.exit(1)

    # Transcode the NCEP GRIB files into ADL NCEP-ANC-Int
    gridBlobFiles = []
    for gribFile in gribFiles:
        gridBlobFile = ANC.create_NCEP_grid_blobs(gribFile)
        gridBlobFiles.append(gridBlobFile)

    if (gridBlobFiles == []) :
        LOG.error('Failed to convert NCEP GRIB  files to blobs, aborting.')
        sys.exit(1)

    # Open the NCEP gridded blob file

    ncepXmlFile = path.join(ADL_HOME,'xml/ANC/NCEP_ANC_Int.xml')
    endian = adl_blob.LITTLE_ENDIAN
    ncepBlobArrObjs = []

    for gridBlobFile in gridBlobFiles :
        timeObj = gridBlobFile[0]
        ncepBlobFile = gridBlobFile[1]
        ncepBlobObj = adl_blob.map(ncepXmlFile,ncepBlobFile, endian=endian)
        ncepBlobArrObj = ncepBlobObj.as_arrays()
        ncepBlobArrObjs.append([timeObj,ncepBlobArrObj])
        LOG.debug("%s...\n%r" % (ncepBlobFile,[field for field in ncepBlobArrObj._fields]))
    
    # Get the NAAPS AOT if required...
    if 'AOT' in algList :
        # Download the required NCEP grib files
        #LOG.info("Downloading NAAPS GRIB ancillary into cache...")
        #gribFiles = ANC.retrieve_NAAPS_grib_files(geoDicts)
        #LOG.debug('Dynamic ancillary GRIB files: %s' % repr(gribFiles))
        #if (gribFiles == []) :
            #LOG.error('Failed to find or retrieve any GRIB files, aborting.')
            #sys.exit(1)

        # Transcode the NAAPS GRIB files into ADL NAAPS-ANC-Int
        #gridBlobFiles = ANC.create_NAAPS_grid_blobs(gribFiles)
        #if (gridBlobFiles == []) :
            #LOG.error('Failed to convert NAAPS GRIB  files to blobs, aborting.')
            #sys.exit(1)

        # Link in the canned NAAPS gridded file
        NAAPSblobFile = path.join(CSPP_RT_ANC_CACHE_DIR,'NAAPS-ANC-Int','template.NAAPS-ANC-Int')
        #NAAPSblobFileLink = path.join(inDir,'template.NAAPS-ANC-Int')
        #os.symlink(NAAPSblobFile, NAAPSblobFileLink)

        # Open the NAAPS gridded blob file
        # FIXME : Should be using two NAAPS blob files, and averaging
        naapsXmlFile = path.join(ADL_HOME,'xml/ANC/NAAPS_ANC_Int.xml')
        #naapsGridBlobFile = glob(path.join(inDir,"*.NAAPS-ANC-Int"))[0]
        naapsGridBlobFile = NAAPSblobFile

        if path.exists(naapsXmlFile):
            LOG.info("We are using for %s: %s,%s" %('NAAPS-ANC-Int',naapsXmlFile,naapsGridBlobFile))
        
        # This is a BIG endian grid blob, usually will be little endian.
        naapsBlobObj = adl_blob.map(naapsXmlFile,naapsGridBlobFile, endian=adl_blob.BIG_ENDIAN)
        naapsBlobArrObj = naapsBlobObj.as_arrays()
        LOG.debug("%s...\n%r" % (naapsGridBlobFile,naapsBlobArrObj._fields))


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
            ANC_objects[shortName].sourceList = [files[1] for files  in gridBlobFiles]
            ANC_objects[shortName].ingest(ancBlob=ncepBlobArrObjs)
        elif ANC_objects[shortName].sourceType == 'NAAPS_ANC_Int' :
            ANC_objects[shortName].sourceList = [naapsGridBlobFile]
            ANC_objects[shortName].ingest(ancBlob=naapsBlobArrObj)
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
            ANC_objects[shortName].shipOutToFile()


def _granulate_GridIP(inDir,geoDicts,algList):
    '''Granulates the input gridded static data into the required GridIP granulated datasets.'''

    import GridIP
    import Algorithms
    global sdrEndian 
    global ancEndian 

    CSPP_RT_HOME = os.getenv('CSPP_RT_HOME')
    ANC_SCRIPTS_PATH = path.join(CSPP_RT_HOME,'viirs')
    CSPP_RT_ANC_CACHE_DIR = os.getenv('CSPP_RT_ANC_CACHE_DIR')
    
    ADL_ASC_TEMPLATES = path.join(ANC_SCRIPTS_PATH,'asc_templates')

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
            GridIP_objects[shortName].shipOutToFile()


def main():

    endianChoices = ['little','big']
    algorithmChoices = ['VCM','AOT','SST']

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

    mandatoryGroup.add_option('--alg',
                      action="store",
                      dest="algorithm",
                      type="choice",
                      #default='little',
                      choices=algorithmChoices,
                      help='''The VIIRS algorithm to run.\n\n
                              Possible values are...
                              %s.
                           ''' % (algorithmChoices.__str__()[1:-1]))

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

    optionalGroup.add_option('--no_chain',
                      action="store_true",
                      dest="noAlgChain",
                      default=False,
                      help="Do not run prerequisite algorithms.")

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
    mandatories = ['inputFiles','algorithm']
    mand_errors = ["Missing mandatory argument [-i input_files --input_files=input_files]",
                   "Missing mandatory argument [--alg=%r ]"%(algorithmChoices)]
    isMissingMand = False
    for m,m_err in zip(mandatories,mand_errors):
        if not options.__dict__[m]:
            isMissingMand = True
            parser.error(m_err)
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


    # Unpack HDF5 VIIRS geolocation SDRs in the input directory to the work directory
    geo_unpacking_problems = 0
    if not options.skipSdrUnpack :
        geo_unpacking_problems = _unpack_sdr(work_dir,input_dir,inputGlobs['GEO'])
    else :
        LOG.info('Skipping SDR GEO unpacking, assuming all VIIRS SDR blob and asc files are present.')


    # Ordered list of required algorithms (to be passed in)
    algList = [options.algorithm]

    # Determine the ordered list of algs to satisfy the required 
    # algorithm's dependencies.
    if not options.noAlgChain : 
        LOG.info("We ARE chaining algorithms...")
        thisAlg = algList[0]
        thisAlgPreReqs = copy.copy(Algorithms.prerequisites[thisAlg])
        thisAlgPreReqs.append(thisAlg)
        algList = thisAlgPreReqs
    else :
        LOG.info("We are NOT chaining algorithms...")

    LOG.info("Required algorithms: %r" % (algList))

    # Determine the number of VIIRS SDR cross granules required to 
    # process the algorithm chain.
    cumulativeCrossGranules = {}
    for alg in algList :
        if not options.noAlgChain : 
            crossSum = Algorithms.crossGranules[alg]
            for preReq in Algorithms.prerequisites[alg]:
                if preReq is None :
                    pass
                else :
                    crossSum += Algorithms.crossGranules[preReq]
            cumulativeCrossGranules[alg] = crossSum
        else :
            cumulativeCrossGranules[alg] = Algorithms.crossGranules[alg]

    for alg in algList :
        LOG.info("We require %d cross granules for %s" % (cumulativeCrossGranules[alg],alg))
    LOG.info("")

    # Create a list of algorithm module "pointers"
    algorithms = []
    for alg in algList :
        algName = Algorithms.modules[alg]
        algorithms.append(getattr(Algorithms,algName))

    # Read through ascii metadata and build up information table
    LOG.info('Sifting through geolocation metadata to find VIIRS SDR processing candidates...')
    #geolocationShortNames = ['VIIRS-MOD-RGEO-TC','VIIRS-MOD-RGEO','VIIRS-MOD-GEO-TC','VIIRS-MOD-GEO']
    geolocationShortNames = ['VIIRS-MOD-GEO-TC']
    anc_granules_to_process = None

    # Determine the candidate geolocation granules for which to generate ancillary.
    # Hint: ... Make ALL THE ANCILLARY! :-)
    for geoType in geolocationShortNames :
        LOG.info("Searching for candidate %s geolocation granules for VIIRS ancillary..." % (geoType))
        anc_granules_to_process = sorted(list(sift_metadata_for_viirs_sdr(geoType,crossGran=None,work_dir=work_dir)))
        if anc_granules_to_process :
            break
        else :
            LOG.info("\tNo %s geolocation granules for VIIRS ancillary" % (geoType))

    # Check geolocation metadata for wrong N_Granule_Version, and fix...
    if anc_granules_to_process :
        for grans in anc_granules_to_process:
            fileChanged=0
            if grans["N_Granule_Version"] == "A2" :
                LOG.info("In %s (%s): N_Granule_Version=%s, fixing..." % \
                        (path.basename(grans["_filename"]),grans["N_Granule_ID"],grans["N_Granule_Version"]))
                fileChanged = _strReplace(grans["_filename"],\
                        '("N_Granule_Version" STRING EQ "A2")','("N_Granule_Version" STRING EQ "A1")')

    # Determine the candidate geolocation granules for which we can generate VIIRS products.
    for alg in algorithms :
        for geoType in geolocationShortNames :
            LOG.info("Searching for candidate %s geolocation granules for VIIRS %s ..." % (geoType,alg.AlgorithmName))
            crossGranules = cumulativeCrossGranules[alg.AlgorithmString]
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

    if anc_granules_to_process :
        LOG.info("We have %d candidate granules for ancillary." % (len(anc_granules_to_process)))
        LOG.info("  %13s%28s%29s" % ('N_Granule_ID','ObservedStartTime','ObservedEndTime'))
        for dicts in sorted(anc_granules_to_process,key=granIdKey) :
            LOG.info("  %15s%30s%30s"%(dicts['N_Granule_ID'],dicts['ObservedStartTime'],dicts['ObservedEndTime']))

    for alg in algorithms :
        if alg.granules_to_process :
            LOG.info("")
            LOG.info("We have %d candidate granules for %s" % (len(alg.granules_to_process),alg.AlgorithmName))
            LOG.info("  %13s%28s%29s" % ('N_Granule_ID','ObservedStartTime','ObservedEndTime'))
            for dicts in sorted(alg.granules_to_process,key=granIdKey) :
                LOG.info("  %15s%30s%30s"%(dicts['N_Granule_ID'],dicts['ObservedStartTime'],dicts['ObservedEndTime']))
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

    
    # Set the VIIRS SDR endianness from the input option...
    global sdrEndian
    set_sdr_endian(options.sdr_Endianness)

    # Set the VIIRS ancillary endianness from the input option...
    global ancEndian
    set_anc_endian(options.anc_Endianness)

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

    # Retrieve and granulate the required ancillary data...
    if not options.skipAncillary :

        t1 = time()

        LOG.info('Retrieving and granulating ancillary data...')

        # Granulate the VIIRS ANC data
        _granulate_ANC(work_dir,anc_granules_to_process,algList)

        # Granulate the VIIRS GridIP data
        _granulate_GridIP(work_dir,anc_granules_to_process,algList)

        t2 = time()

        LOG.info( "Generation of VIIRS ancillary data took %f seconds." % (t2-t1))

    else :

        LOG.info('Skipping retrieval and granulation of ancillary data.')


    # Unpack HDF5 VIIRS radiometric SDRs in the input directory to the work directory
    mod_unpacking_problems = 0
    img_unpacking_problems = 0
    unpacking_problems = 0
    if not options.skipSdrUnpack :
        mod_unpacking_problems = _unpack_sdr(work_dir,input_dir,inputGlobs['MOD'])
        img_unpacking_problems = _unpack_sdr(work_dir,input_dir,inputGlobs['IMG'])

        unpacking_problems = geo_unpacking_problems + mod_unpacking_problems + img_unpacking_problems
        LOG.debug("Total VIIRS SDR unpacking problems = %d" % (unpacking_problems))
    else :
        LOG.info('Skipping SDR unpacking, assuming all VIIRS SDR blob and asc files are present.')
        unpacking_problems = geo_unpacking_problems + mod_unpacking_problems + img_unpacking_problems

    # Check radiometric SDR metadata for wrong N_Granule_Version, and fix...
    sdr_blob_names = glob(path.join(work_dir,"*.*SDR"))
    for sdr_blob_name in sdr_blob_names:
        ascName = corresponding_asc_path(sdr_blob_name)
        fileChanged = _strReplace(ascName,'("N_Granule_Version" STRING EQ "A2")','("N_Granule_Version" STRING EQ "A1")')
        if fileChanged :
            LOG.info("Fixed N_Granule_Version in %s metadata" % (path.basename(sdr_blob_name)))

    # Link in auxillary files
    if not options.skipAuxLinking :
        LOG.info('Linking in the VIIRS EDR auxillary files...')

        # Create a list of algorithm module "pointers"
        algs_for_Aux = []
        for alg in algList :
            algName = Algorithms.modules[alg]
            algs_for_Aux.append(getattr(Algorithms,algName))

        for alg in algs_for_Aux :
            alg.setupAuxillaryFiles(alg, work_dir)

        del(algs_for_Aux)

    else :
        LOG.info('Skipping linking in the VIIRS EDR auxillary files.')


    # Specify the algorithm we want to run via the Lw XML file.
    if not options.skipAlgorithm :

        # Create a list of algorithm module "pointers"
        algorithms = []
        for alg in algList :
            algName = Algorithms.modules[alg]
            algorithms.append(getattr(Algorithms,algName))

        for alg in algorithms :

            # build XML configuration files for jobs that can be run
            LOG.info("Building %s XML files for %d granules" % \
                    (alg.AlgorithmString,len(alg.granules_to_process)))

            # Generate the VIIRS LW xml files for each N_Granule_ID for this algorithm.
            xml_files_to_process = alg.generate_viirs_edr_xml(work_dir, alg.granules_to_process)
        
            LOG.info('%d granules to process: %s' % \
                    (len(xml_files_to_process), ''.join(name+' -> '+xmlfile+'\n' \
                    for (name,xmlfile) in xml_files_to_process)))

            LOG.info("Running VIIRS %s ..." % (alg.AlgorithmName))
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

            # If this alg failed, return error code and exit
            if rc != 0 :
                LOG.warn("Non-zero error code %d for %s, aborting." % (rc, alg.AlgorithmName))
                return rc

        ## if no errors or only non-critical errors: clean up
        LOG.info("Return code : %d" % (rc))
        if rc == 0 and not options.cspp_debug:
            LOG.info("Cleaning up workspace...")

            for alg in algorithms :

                algorithmXmlGlob = '%s*.xml' % (alg.algorithmLWxml)
                algorithmLogGlob = '%s_*' % (alg.controllerBinary)

                alg.cleanup(work_dir, algorithmXmlGlob, algorithmLogGlob, log_dir)

        return rc
    
    else :

        LOG.info("Skipping execution of VIIRS %s ..." % (alg.AlgorithmName))


if __name__=='__main__':
    sys.exit(main())  
