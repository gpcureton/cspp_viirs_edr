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
from subprocess import CalledProcessError, call
from glob import glob
from time import time
from shutil import rmtree

from Utils import check_log_files

# skim and convert routines for reading .asc metadata fields of interest
#import adl_asc
#from adl_asc import skim_dir, contiguous_granule_groups, granule_groups_contain, effective_anc_contains,_eliminate_duplicates,_is_contiguous, RDR_REQUIRED_KEYS, POLARWANDER_REQUIRED_KEYS
from adl_common import sh, unpack, env, h5_xdr_inventory
from adl_common import ADL_HOME, CSPP_RT_ANC_PATH, CSPP_RT_ANC_CACHE_DIR, COMMON_LOG_CHECK_TABLE

# log file scanning
import adl_log

# every module should have a LOG object
try :
    sourcename= file_Id.split(" ")
    LOG = logging.getLogger(sourcename[1])
except :
    LOG = logging.getLogger('AerosolOpticalThicknessIP')

AlgorithmString = 'AOT'

ANC_collectionShortNames = [
                           'VIIRS-ANC-Preci-Wtr-Mod-Gran',
                           'VIIRS-ANC-Temp-Surf2M-Mod-Gran',
                           'VIIRS-ANC-Wind-Speed-Mod-Gran',
                           'VIIRS-ANC-Wind-Direction-Mod-Gran',
                           'VIIRS-ANC-Tot-Col-Mod-Gran',
                           'VIIRS-ANC-Press-Surf-Mod-Gran',
                           'VIIRS-ANC-Optical-Depth-Mod-Gran'
                          ]

GridIP_collectionShortNames = []

controllerBinary = 'ProEdrViirsAerosolController.exe'
ADL_VIIRS_AEROSOL_EDR=path.abspath(path.join(ADL_HOME, 'bin', 'ProEdrViirsAerosolController.exe'))

algorithmLWxml = 'edr_viirs_aerosol'

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


def generate_viirs_edr_xml(work_dir, granule_seq):
    "generate XML files for VIIRS Masks EDR granule generation"
    to_process = []
    for gran in granule_seq:
        name = gran['N_Granule_ID']
        fnxml = '%s_%s.xml' % (algorithmLWxml,name)
        LOG.debug('writing XML file %r' % (fnxml))
        fpxml = file(path.join(work_dir, fnxml), 'wt')
        fpxml.write(xmlTemplate % gran)
        to_process.append([name,fnxml])
    return to_process


def cleanup(work_dir, xml_glob, log_dir_glob, *more_dirs):
    """upon successful run, clean out work directory"""

    LOG.info("Cleaning up work directory...")
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

    LOG.info("Removing other directories ...")
    for dirname in more_dirs:
        fullDirName = path.join(work_dir,dirname)
        LOG.debug('removing %s' % (fullDirName))
        try :
            rmtree(fullDirName, ignore_errors=False)
        except Exception, err:
            LOG.warn( "%s" % (str(err)))


def run_xml_files(work_dir, xml_files_to_process, setup_only=False, **additional_env):
    """Run each VIIRS EDR AEROSOL XML input in sequence.
       Return the list of granule IDs which crashed, 
       and list of granule IDs which did not create output.
    """
    crashed_runs = set()
    no_output_runs = set()
    geo_problem_runs = set()
    bad_log_runs = set()
    first = True

    # obtain pre-existing granule list
    aotIpPattern = path.join(work_dir, 'IVAOT*.h5')
    aotEdrPattern = path.join(work_dir, 'VAOOO*.h5')
    suspMatEdrPattern = path.join(work_dir, 'VSUMO*.h5')

    # prior_granules dicts contain (N_GranuleID,HDF5File) key,value pairs.
    # *ID dicts contain (HDF5File,N_GranuleID) key,value pairs.

    # Get the (N_GranuleID,hdfFileName) pairs for the existing Aerosol IP files
    aerosolIP_prior_granules, aotIp_ID = h5_xdr_inventory(aotIpPattern, AOT_IP_GRANULE_ID_ATTR_PATH)
    LOG.debug('Existing IVAOT granules... %s' % (repr(cmask_prior_granules)))

    aerosolIP_prior_granules = set(aerosolIP_prior_granules.keys())
    LOG.debug('Set of existing IVAOT granules... %s' % (repr(aerosolIP_prior_granules)))

    # Get the (N_GranuleID,hdfFileName) pairs for the existing Aerosol EDR files
    aerosolEDR_prior_granules, aotEdr_ID = h5_xdr_inventory(aotEdrPattern, AOT_EDR_GRANULE_ID_ATTR_PATH)
    LOG.debug('Existing VAOOO granules... %s' % (repr(cmask_prior_granules)))

    aerosolEDR_prior_granules = set(aerosolEDR_prior_granules.keys())
    LOG.debug('Set of existing VAOOO granules... %s' % (repr(aerosolEDR_prior_granules)))
    
    # Get the (N_GranuleID,hdfFileName) pairs for the existing Suspended Matter EDR files
    suspMatEDR_prior_granules, suspMatEdr_ID = h5_xdr_inventory(suspMatEdrPattern, SUSMAT_EDR_GRANULE_ID_ATTR_PATH)
    LOG.debug('Existing VSUMO granules... %s' % (repr(cmask_prior_granules)))

    suspMatEDR_prior_granules = set(suspMatEDR_prior_granules.keys())
    LOG.debug('Set of existing VSUMO granules... %s' % (repr(suspMatEDR_prior_granules)))


    for granule_id, xml in xml_files_to_process:

        t1 = time()
        
        cmd = [ADL_VIIRS_AEROSOL_EDR, xml]
        
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
                LOG.error('%s failed on %r: %r. Continuing...' % (controllerBinary, xml, oops))
                crashed_runs.add(granule_id)
            first = False

            # check new IVAOT output granules
            aotIp_new_granules, aotIp_ID = h5_xdr_inventory(aotIpPattern, AOT_IP_GRANULE_ID_ATTR_PATH, state=aotIp_ID)
            LOG.debug('new IVAOT granules after this run: %s' % repr(aotIp_new_granules))
            if granule_id not in aotIp_new_granules:
                LOG.warning('no IVAOT HDF5 output for %s' % granule_id)
                no_output_runs.add(granule_id)
            else:
                filename = aotIp_new_granules[granule_id]

            # check new VAOOO output granules
            aotEdr_new_granules, aotEdr_ID = h5_xdr_inventory(aotEdrPattern, AOT_EDR_GRANULE_ID_ATTR_PATH, state=aotEdr_ID)
            LOG.debug('new VAOOO granules after this run: %s' % repr(aotEdr_new_granules))
            if granule_id not in aotEdr_new_granules:
                LOG.warning('no VAOOO HDF5 output for %s' % granule_id)
                no_output_runs.add(granule_id)
            else:
                filename = aotEdr_new_granules[granule_id]

            # check new VSUMO output granules
            suspMatEdr_new_granules, suspMatEdr_ID = h5_xdr_inventory(suspMatEdrPattern, SUSMAT_EDR_GRANULE_ID_ATTR_PATH, state=suspMatEdr_ID)
            LOG.debug('new VSUMO granules after this run: %s' % repr(suspMatEdr_new_granules))
            if granule_id not in suspMatEdr_new_granules:
                LOG.warning('no VSUMO HDF5 output for %s' % granule_id)
                no_output_runs.add(granule_id)
            else:
                filename = suspMatEdr_new_granules[granule_id]

        t2 = time()
        LOG.info ( "Controller ran in %f seconds." % (t2-t1))


    LOG.debug("aotIp_ID.values() = \n%r" % (aotIp_ID.values()))
    LOG.debug("set(aotIp_ID.values()) = \n%r" % (set(aotIp_ID.values())))
    LOG.debug("aerosolIP_prior_granules = \n%r" % (aerosolIP_prior_granules))
    aerosolIP_granules_made = set(aotIp_ID.values()) - aerosolIP_prior_granules

    LOG.debug("aotEdr_ID.values() = \n%r" % (aotEdr_ID.values()))
    LOG.debug("set(aotEdr_ID.values()) = \n%r" % (set(aotEdr_ID.values())))
    LOG.debug("aerosolEDR_prior_granules = \n%r" % (aerosolEDR_prior_granules))
    aerosolEDR_granules_made = set(aotEdr_ID.values()) - aerosolEDR_prior_granules

    LOG.debug("suspMatEdr_ID.values() = \n%r" % (suspMatEdr_ID.values()))
    LOG.debug("set(suspMatEdr_ID.values()) = \n%r" % (set(suspMatEdr_ID.values())))
    LOG.debug("suspMatEDR_prior_granules = \n%r" % (suspMatEDR_prior_granules))
    suspMatEDR_granules_made = set(suspMatEdr_ID.values()) - suspMatEDR_prior_granules

    LOG.info('aerosolIP granules created: %s' % ', '.join(list(aerosolIP_granules_made)))
    LOG.info('aerosolEDR granules created: %s' % ', '.join(list(aerosolEDR_granules_made)))
    LOG.info('suspMatEDR granules created: %s' % ', '.join(list(suspMatEDR_granules_made)))


    if no_output_runs:
        LOG.info('Granules that failed to generate output: %s' % (', '.join(no_output_runs)))
    if geo_problem_runs:
        LOG.warning('Granules which had no N_Geo_Ref: %s' % ', '.join(geo_problem_runs))
    if crashed_runs:
        LOG.warning('Granules that crashed ADL: %s' % (', '.join(crashed_runs)))
    if bad_log_runs:
        LOG.warning('Granules that produced logs indicating problems: %s' % (', '.join(bad_log_runs)))
    if not aerosolIP_granules_made:
        LOG.warning('No Aerosol Optical Thickness IP HDF5 files were created')
    if not aerosolEDR_granules_made:
        LOG.warning('No Aerosol Optical Thickness EDR HDF5 files were created')
    if not suspMatEDR_granules_made:
        LOG.warning('No Suspended Matter EDR HDF5 files were created')

    return crashed_runs, no_output_runs, geo_problem_runs, bad_log_runs


def aggregate(fileGlob):
    '''
    Aggregate the files given by the input file glob.
    '''
    pass

