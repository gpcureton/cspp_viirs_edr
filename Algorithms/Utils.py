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
from glob import glob
import uuid
from datetime import datetime

from adl_common import ADL_HOME, CSPP_RT_HOME, CSPP_RT_ANC_PATH, CSPP_RT_ANC_HOME, CSPP_RT_ANC_CACHE_DIR, COMMON_LOG_CHECK_TABLE


# log file scanning
import adl_log

# every module should have a LOG object
try :
    sourcename= file_Id.split(" ")
    LOG = logging.getLogger(sourcename[1])
except :
    LOG = logging.getLogger('Utils')


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
    logpat = path.join(work_dir, "log", "*%d*.log" % pid)

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


def _setupAuxillaryFiles(Alg_objects,workDir):
    '''
    Create asc files for the various auxillary files (tunable parameters etc...) 
    in the workDir, and create links to the binary auxillary files in workDir.
    '''


    ANC_SCRIPTS_PATH = path.join(CSPP_RT_HOME,'viirs')
    ADL_ASC_TEMPLATES = path.join(ANC_SCRIPTS_PATH,'asc_templates')


    auxillaryCollShortNames = Alg_objects.AUX_collectionShortNames

    auxillaryAscTemplateFile = Alg_objects.AUX_ascTemplateFile

    auxillaryBlobTemplateFile = Alg_objects.AUX_blobTemplateFile

    auxillaryPaths = Alg_objects.AUX_Paths


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

        URID_dict = getURID()

        URID = URID_dict['URID']
        creationDate_nousecStr = URID_dict['creationDate_nousecStr']
        creationDateStr = URID_dict['creationDateStr']

        # The names for the new asc and blob files
        ascFileName = path.join(workDir,URID+'.asc')
        blobFileName = path.join(workDir,string.replace(blobTempFileName,'template',URID))

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

