#!/usr/bin/env python
# encoding: utf-8
"""
Utils.py

Various methods that are used by other methods in the ANC module.

Created by Geoff Cureton on 2013-03-04.
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
import uuid
from datetime import datetime,timedelta

import adl_blob
from adl_common import ADL_HOME, CSPP_RT_ANC_PATH, CSPP_RT_ANC_CACHE_DIR, COMMON_LOG_CHECK_TABLE

# every module should have a LOG object
try :
    sourcename= file_Id.split(" ")
    LOG = logging.getLogger(sourcename[1])
except :
    LOG = logging.getLogger('Utils')


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


def getURID() :
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


def getAscLine(fileObj,searchString):
    ''' Parses a file and searches for a string in each line, returning 
        the line if the string is found.'''

    dataStr = ''
    try :
        while True :
            line = fileObj.readline()

            if searchString in line : 
                dataStr = "%s" % (string.replace(line,'\n',''));
                break

        fileObj.seek(0)

    except Exception, err:
        LOG.error('Exception: %r' % (err))
        fileObj.close()

    return dataStr


def getAscStructs(fileObj,searchString,linesOfContext):
    ''' Parses a file and searches for a string in each line, returning 
        the line (and a given number of lines of context) if the string 
        is found.'''

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
        LOG.error('Exception: %r' % (err))
        fileObj.close()
        return -1

    dataStr=''
    dataStr = "%s" % ("\n").join(['%s' % (str(lines)) for lines in dataList])

    return dataStr


def shipOutToFile(GridIPobj):
    '''
    Generate a blob/asc file pair from the input ancillary data object.
    '''

    # Set some environment variables and paths
    CSPP_RT_HOME = os.getenv('CSPP_RT_HOME')
    ANC_SCRIPTS_PATH = path.join(CSPP_RT_HOME,'viirs')
    CSPP_RT_ANC_CACHE_DIR = os.getenv('CSPP_RT_ANC_CACHE_DIR')
    ADL_ASC_TEMPLATES = path.join(ANC_SCRIPTS_PATH,'asc_templates')

    # Create new GridIP ancillary blob, and copy granulated data to it

    endian = GridIPobj.ancEndian
    xmlName = path.join(ADL_HOME,'xml/VIIRS',GridIPobj.xmlName)

    # Create a new URID to be used in making the asc filenames

    URID_dict = getURID()

    URID = URID_dict['URID']
    creationDate_nousecStr = URID_dict['creationDate_nousecStr']
    creationDateStr = URID_dict['creationDateStr']

    # Create a new directory in the input directory for the new ancillary
    # asc and blob files

    blobDir = GridIPobj.inDir

    ascFileName = path.join(blobDir,URID+'.asc')
    blobName = path.join(blobDir,URID+'.'+GridIPobj.collectionShortName)

    LOG.debug("ascFileName : %s" % (ascFileName))
    LOG.debug("blobName : %s" % (blobName))

    # Create a new ancillary blob, and copy the data to it.
    newGridIPblobObj = adl_blob.create(xmlName, blobName, endian=endian, overwrite=True)
    newGridIPblobArrObj = newGridIPblobObj.as_arrays()

    blobData = getattr(newGridIPblobArrObj,'data')
    blobData[:,:] = GridIPobj.data[:,:]

    # Make a new GridIP asc file from the template, and substitute for the various tags

    ascTemplateFileName = path.join(ADL_ASC_TEMPLATES,"VIIRS-ANC_Template.asc")

    LOG.debug("Creating new asc file\n%s\nfrom template\n%s" % (ascFileName,ascTemplateFileName))
    
    ANC_fileList = GridIPobj.sourceList
    for idx in range(len(ANC_fileList)) :
        ANC_fileList[idx] = path.basename(ANC_fileList[idx])
    ANC_fileList.sort()
    ancGroupRecipe = '    ("N_Anc_Filename" STRING EQ "%s")'
    ancFileStr = "%s" % ("\n").join([ancGroupRecipe % (str(files)) for files in ANC_fileList])

    LOG.debug("RangeDateTimeStr = %s\n" % (GridIPobj.RangeDateTimeStr))
    LOG.debug("GRingLatitudeStr = \n%s\n" % (GridIPobj.GRingLatitudeStr))
    LOG.debug("GRingLongitudeStr = \n%s\n" % (GridIPobj.GRingLongitudeStr))

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
       line = line.replace("CSPP_ANC_BLOB_FULLPATH",path.basename(blobName))
       line = line.replace("CSPP_ANC_COLLECTION_SHORT_NAME",GridIPobj.collectionShortName)
       line = line.replace("CSPP_GRANULE_ID",GridIPobj.geoDict['N_Granule_ID'])
       line = line.replace("CSPP_CREATIONDATETIME",creationDateStr)
       line = line.replace("  CSPP_RANGE_DATE_TIME",GridIPobj.RangeDateTimeStr)
       line = line.replace("  CSPP_GRINGLATITUDE",GridIPobj.GRingLatitudeStr)
       line = line.replace("  CSPP_GRINGLONGITUDE",GridIPobj.GRingLongitudeStr)
       line = line.replace("    CSPP_ANC_SOURCE_FILES",ancFileStr)
       ascFile.write(line) 

    ascFile.close()
    ascTemplateFile.close()
