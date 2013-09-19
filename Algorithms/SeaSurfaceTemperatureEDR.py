#!/usr/bin/env python
# encoding: utf-8
"""
SeaSurfaceTemperatureEDR.py

 * DESCRIPTION:  Class containing data relevent to the VIIRS Sea Surface Temperature EDR. 

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

from Utils import check_log_files, _setupAuxillaryFiles

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
    LOG = logging.getLogger('SeaSurfaceTemperatureEDR')

AlgorithmString = 'SST'

AlgorithmName = 'Sea Surface Temperature EDR'

GEO_collectionShortNames = [
                            'VIIRS-MOD-GEO-TC'
                          ]

SDR_collectionShortNames = [
                              'VIIRS-M12-SDR',
                              'VIIRS-M15-SDR',
                              'VIIRS-M16-SDR'
                          ]

ANC_collectionShortNames = [
                           'VIIRS-ANC-Preci-Wtr-Mod-Gran',
                           'VIIRS-ANC-Temp-Surf2M-Mod-Gran',
                           'VIIRS-ANC-Wind-Speed-Mod-Gran',
                           'VIIRS-ANC-Wind-Direction-Mod-Gran',
                           'VIIRS-ANC-Tot-Col-Mod-Gran',
                           'VIIRS-ANC-Press-Surf-Mod-Gran',
                           'VIIRS-ANC-Optical-Depth-Mod-Gran'
                          ]
ANC_collectionShortNames = [
                             'VIIRS-ANC-Temp-Skin-Mod-Gran'
                           ]

GridIP_collectionShortNames = [
                                'VIIRS-I-Conc-IP'
                              ]

AUX_collectionShortNames = [
                            'VIIRS-SST-EDR-AC-Int',
                            'VIIRS-SST-EDR-DQTT-Int',
                            'VIIRS-SST-Coef-LUT'
                           ]

AUX_ascTemplateFile = [
                        'VIIRS-SST-EDR-AC-Int_Template.asc',
                        'VIIRS-SST-EDR-DQTT-Int_Template.asc',
                        'VIIRS-SST-Coef-LUT_Template.asc'
                      ]


AUX_blobTemplateFile = [
                        'template.VIIRS-SST-EDR-AC-Int',
                        'template.VIIRS-SST-EDR-DQTT-Int',
                        'template.VIIRS-SST-Coef-LUT'
                       ]

AUX_Paths = [
             'luts/viirs',
             'luts/viirs',
             'luts/viirs',
            ]

EDR_collectionShortNames = [
                           'VIIRS-SST-EDR'
                          ]

controllerBinary = 'ProEdrViirsSstController.exe'
ADL_VIIRS_SST_EDR=path.abspath(path.join(ADL_HOME, 'bin', controllerBinary))

algorithmLWxml = 'edr_viirs_sst'

# Attribute paths for SST EDR
attributePaths = {}
attributePaths['MOD_GEO_TC'] = 'Data_Products/VIIRS-MOD-GEO-TC/VIIRS-MOD-GEO-TC_Gran_0/N_Granule_ID'
attributePaths['CM_IP'] = 'Data_Products/VIIRS-CM-IP/VIIRS-CM-IP_Gran_0/N_Granule_ID'
attributePaths['AOT_IP'] = 'Data_Products/VIIRS-Aeros-Opt-Thick-IP/VIIRS-Aeros-Opt-Thick-IP_Gran_0/N_Granule_ID'
attributePaths['SST_EDR'] = 'Data_Products/VIIRS-SST-EDR/VIIRS-SST-EDR_Gran_0/N_Granule_ID'

MOD_GEO_TC_GRANULE_ID_ATTR_PATH = 'Data_Products/VIIRS-MOD-GEO-TC/VIIRS-MOD-GEO-TC_Gran_0/N_Granule_ID'
CM_GRANULE_ID_ATTR_PATH = 'Data_Products/VIIRS-CM-IP/VIIRS-CM-IP_Gran_0/N_Granule_ID'
AOT_IP_GRANULE_ID_ATTR_PATH = 'Data_Products/VIIRS-Aeros-Opt-Thick-IP/VIIRS-Aeros-Opt-Thick-IP_Gran_0/N_Granule_ID'
SST_GRANULE_ID_ATTR_PATH = 'Data_Products/VIIRS-SST-EDR/VIIRS-SST-EDR_Gran_0/N_Granule_ID'

# XML template for ProEdrViirsSstController.exe
# from ADL/cfg/dynamic/withMetadata/ProEdrViirsSstControllerLwFile.xml
xmlTemplate = """<InfTkConfig>
  <idpProcessName>ProEdrViirsSstController.exe</idpProcessName>
  <siSoftwareId></siSoftwareId>
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
  <actualScans>0</actualScans>
  <previousActualScans>0</previousActualScans>
  <nextActualScans>0</nextActualScans> 
  <usingMetadata>TRUE</usingMetadata>
  <configGuideName>ProEdrViirsSSTController_GuideList.cfg</configGuideName>

  <task>
    <taskType>EDR</taskType>
    <taskDetails1>%(N_Granule_ID)s</taskDetails1>
    <taskDetails2>%(N_Granule_Version)s</taskDetails2>
    <taskDetails3>NPP</taskDetails3>
    <taskDetails4>VIIRS</taskDetails4>
  </task>

</InfTkConfig>
"""


def setupAuxillaryFiles(Alg_objects,workDir):
    '''
    Call the generic Utils method to link in the auxillary files 
    specified in Alg_objects to the workDir directory.
    '''

    _setupAuxillaryFiles(Alg_objects,workDir)


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


def run_xml_files(work_dir, xml_files_to_process, setup_only=False, **additional_env):
    """Run each VIIRS Sea Surface Temperature EDR xml input in sequence.
       Return the list of granule IDs which crashed, 
       and list of granule IDs which did not create output.
    """
    crashed_runs = set()
    no_output_runs = set()
    geo_problem_runs = set()
    bad_log_runs = set()
    first = True

    # obtain pre-existing granule list
    modGeoTCPattern = path.join(work_dir, 'GMTCO*.h5')
    sstPattern = path.join(work_dir, 'VSSTO*.h5')

    # prior_granules dicts contain (N_GranuleID,HDF5File) key,value pairs.
    # *ID dicts contain (HDF5File,N_GranuleID) key,value pairs.
    
    # Get the N_GranuleID and filename of existing sst files in the work_dir ? ...
    sst_prior_granules, sst_ID = h5_xdr_inventory(sstPattern, SST_GRANULE_ID_ATTR_PATH)
    LOG.debug('Existing VSSTO granules... %s' % (repr(sst_prior_granules)))

    sst_prior_granules = set(sst_prior_granules.keys())
    LOG.debug('Set of existing VSSTO granules... %s' % (repr(sst_prior_granules)))


    for granule_id, xml in xml_files_to_process:

        t1 = time()
        
        cmd = [ADL_VIIRS_SST_EDR, xml]
        
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

            # check new VSSTO output granules
            sst_new_granules, sst_ID = h5_xdr_inventory(sstPattern, SST_GRANULE_ID_ATTR_PATH, state=sst_ID)
            LOG.debug('new VSSTO granules after this run: %s' % (repr(sst_new_granules)))
            if granule_id not in sst_new_granules:
                LOG.warning('no VSSTO HDF5 output for %s' % (granule_id))
                no_output_runs.add(granule_id)
            else:
                filename = sst_new_granules[granule_id]

        t2 = time()
        LOG.info ( "Controller ran in %f seconds." % (t2-t1))


    LOG.debug("sst_ID.values() = %r" % (sst_ID.values()))
    LOG.debug("set(sst_ID.values()) = %r" % (set(sst_ID.values())))
    LOG.debug("sst_prior_granules = %r" % (sst_prior_granules))
    sst_granules_made = set(sst_ID.values()) - sst_prior_granules

    LOG.info('Sea Surface Temperature granules created: %s' %( ', '.join(list(sst_granules_made))))


    if no_output_runs:
        LOG.info('Granules that failed to generate output: %s' % (', '.join(no_output_runs)))
    if geo_problem_runs:
        LOG.warning('Granules which had no N_Geo_Ref: %s' % ', '.join(geo_problem_runs))
    if crashed_runs:
        LOG.warning('Granules that crashed ADL: %s' % (', '.join(crashed_runs)))
    if bad_log_runs:
        LOG.warning('Granules that produced logs indicating problems: %s' % (', '.join(bad_log_runs)))
    if not sst_granules_made:
        LOG.warning('No Sea Surface Temperature EDR HDF5 files were created')

    return crashed_runs, no_output_runs, geo_problem_runs, bad_log_runs


def cleanup(work_dir, xml_glob, log_dir_glob, *more_dirs):
    """upon successful run, clean out work directory"""

    LOG.info("Cleaning up work directory...")

    # Remove asc/blob file pairs...
    ascBlobFiles = glob(path.join(work_dir, '????????-?????-????????-????????.*'))
    if ascBlobFiles != [] :
        for files in ascBlobFiles:
            LOG.debug('removing %s' % (files))
            os.unlink(files)

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


def aggregate(fileGlob):
    '''
    Aggregate the files given by the input file glob.
    '''
    pass

