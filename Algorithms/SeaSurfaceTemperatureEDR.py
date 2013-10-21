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
from time import time, sleep
from datetime import datetime, timedelta
from shutil import rmtree, move
import multiprocessing

from Utils import check_log_files, _setupAuxillaryFiles

from adl_common import sh, unpack, env, h5_xdr_inventory
from adl_common import ADL_HOME, CSPP_RT_HOME, CSPP_RT_ANC_PATH, CSPP_RT_ANC_HOME, CSPP_RT_ANC_CACHE_DIR, COMMON_LOG_CHECK_TABLE
from adl_post_process import repack_products, aggregate_products
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

controllerName = 'ProEdrViirsSstController'
controllerBinary = '{}.exe'.format(controllerName)
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
  <debugLogLevel>HIGH</debugLogLevel>
  <debugLevel>DBG_HIGH</debugLevel>
  <dbgDest>D_FILE</dbgDest>
  <enablePerf>FALSE</enablePerf>
  <perfPath>${WORK_DIR}/perf</perfPath>
  <dbgPath>${GRANULE_OUTPUT_DIR}/log</dbgPath>
  <initData>
     <domain>OPS</domain>
     <subDomain>SUBDOMAIN</subDomain>
     <startMode>INF_STARTMODE_COLD</startMode>
     <executionMode>INF_EXEMODE_PRIMARY</executionMode>
     <healthTimeoutPeriod>30</healthTimeoutPeriod>
  </initData>
  <lockinMem>FALSE</lockinMem>
  <rootDir>${GRANULE_OUTPUT_DIR}/log</rootDir>
  <inputPath>${WORK_DIR}</inputPath>
  <outputPath>${GRANULE_OUTPUT_DIR}</outputPath>
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


def move_products_to_work_directory(work_dir):
    """ Checks that proper products were produced.
    If product h5 file was produced blob and asc files are deleted
    If product was not produced files are left alone
    """
    dest = path.dirname(work_dir)

    h5Files          = glob(path.join(work_dir, "*.h5"))
    sstEdrBlobFiles  = glob(path.join(work_dir, "*.VIIRS-SST-EDR"))
    sstFSdrBlobFiles = glob(path.join(work_dir, "*.VIIRS-SST-FEDR"))
    ascFiles         = glob(path.join(work_dir, "*.asc"))

    files = h5Files + sstEdrBlobFiles + sstFSdrBlobFiles + ascFiles

    for file in files:
        try :
            LOG.debug( "Moving {} to {}".format(file,dest))
            move(file, dest)
        except Exception, err:
            LOG.warn( "%s" % (str(err)))

    return


def submit_granule(additional_env):
    "run a VIIRS EDR XML input in sequence"

    work_dir = additional_env['WORK_DIR']
    xml = additional_env['XML_FILE']
    granule_id = additional_env['N_Granule_ID']
    granule_output_dir = additional_env['GRANULE_OUTPUT_DIR']
    compress = additional_env['COMPRESS']

    # Pattern for expected output
    sstPattern = path.join(granule_output_dir, 'VSSTO*.h5')

    # prior_granules dicts contain (N_GranuleID,HDF5File) key,value pairs.
    # *ID dicts contain (HDF5File,N_GranuleID) key,value pairs.
    
    # Get the (N_GranuleID,hdfFileName) pairs for the existing Cloud Mask IP files
    sst_prior_granules, sst_ID = h5_xdr_inventory(sstPattern, SST_GRANULE_ID_ATTR_PATH)
    LOG.debug('Existing VSSTO granules in {}... {}'.format(granule_output_dir, sst_prior_granules))

    sst_prior_granules = set(sst_prior_granules.keys())
    LOG.debug('Set of existing VSSTO granules in {}... {}'.format(granule_output_dir, sst_prior_granules))

    # Specify the command line to execute.
    cmd = [ADL_VIIRS_SST_EDR, xml]
    #cmd = ['/bin/sleep','0.2']
    #cmd = ['/usr/bin/gdb', ADL_VIIRS_SST_EDR] # for debugging with gdb...

    LOG.info('executing "{}"'.format(' '.join(cmd)))
    LOG.debug('additional environment variables: {}'.format(additional_env))

    granule_diagnostic = {}
    granule_diagnostic['bad_log'] = False
    granule_diagnostic['crashed'] = False
    granule_diagnostic['geo_problem'] = False
    granule_diagnostic['N_Granule_ID'] = granule_id
    granule_diagnostic['no_output'] = []
    granule_diagnostic['output_file'] = []

    t1 = time()

    try:
        
        pid = sh(cmd, env=env(**additional_env), cwd=work_dir)

        LOG.info("{} ran as pid {}".format(cmd, pid))
        if not check_log_files(granule_output_dir, pid, xml):
            granule_diagnostic['bad_log'] = True

    except CalledProcessError as oops:
        pid = getattr(oops, 'pid', None)
        LOG.debug(traceback.format_exc())
        LOG.error('{} failed on {}: {}. Continuing...' % (cmd[0], xml, oops))
        granule_diagnostic['crashed'] = True

    t2 = time()

    LOG.info("{}({}) ran in {} seconds.".format(controllerName,granule_id,t2-t1))

    try :

        # check new SST output granules
        sst_new_granules, sst_ID = h5_xdr_inventory(sstPattern, SST_GRANULE_ID_ATTR_PATH, state=sst_ID)

        if granule_id not in sst_new_granules:
            LOG.warning('no VSSTO HDF5 output for {}'.format(granule_id))
            granule_diagnostic['no_output'].append(True)
            granule_diagnostic['output_file'].append(None)
        else :
            LOG.info('New VSSTO granule: {}'.format(repr(sst_new_granules)))
            sst_granules_made = set(sst_ID.values()) - sst_prior_granules
            LOG.info('{} granules created: {}'.format(AlgorithmName,', '.join(list(sst_granules_made))))
            granule_diagnostic['no_output'].append(False)
            granule_diagnostic['output_file'].append(path.basename(sst_new_granules[granule_id]))

    except Exception:
        LOG.warn(traceback.format_exc())

    if compress == "True":
        LOG.info('Compress products for {}'.format(granule_id))
        repack_products(granule_output_dir, EDR_collectionShortNames)
        LOG.info('Compress products for {} complete.'.format(granule_id))

    move_products_to_work_directory(granule_output_dir)

    return granule_diagnostic


def run_xml_files(work_dir, xml_files_to_process, nprocs=1, CLEANUP="True",COMPRESS=False,AGGREGATE=False, **additional_env):
    """Run each VIIRS Sea Surface Temperature EDR xml input in sequence.
       Return the list of granule IDs which crashed, 
       and list of granule IDs which did not create output.
    """

    total_granules = len(xml_files_to_process)
    LOG.info('{} granules to process'.format(total_granules))

    argument_dictionaries = []
    for granule_id, xml in xml_files_to_process:

        granule_output_dir = path.join(work_dir,"{}_{}".format(controllerName,granule_id))

        if not path.exists(granule_output_dir): os.mkdir(granule_output_dir)
        if not path.exists(path.join(granule_output_dir, "log")): os.mkdir(path.join(granule_output_dir, "log"))


        additional_envs = dict(
            COMPRESS=str(COMPRESS),
            N_Granule_ID=granule_id,
            XML_FILE=xml,
            GRANULE_OUTPUT_DIR=granule_output_dir,
            WORK_DIR=work_dir,
            ADL_HOME=ADL_HOME,
            CLEANUP=CLEANUP
        )

        argument_dictionaries.append(additional_envs)

    # Pattern for expected output in the root working directory
    sstPattern = path.join(work_dir, 'VSSTO*.h5')

    # prior_granules dicts contain (N_GranuleID,HDF5File) key,value pairs.
    # *ID dicts contain (HDF5File,N_GranuleID) key,value pairs.
    
    # Get the (N_GranuleID,hdfFileName) pairs for the existing Cloud Mask IP files
    sst_prior_granules, sst_ID = h5_xdr_inventory(sstPattern, SST_GRANULE_ID_ATTR_PATH)
    LOG.debug('Existing VSSTO granules in {}... {}'.format(granule_output_dir, sst_prior_granules))

    sst_prior_granules = set(sst_prior_granules.keys())
    LOG.debug('Set of existing VSSTO granules in {}... {}'.format(work_dir, sst_prior_granules))

    results = []

    if len(argument_dictionaries) > 0:

        # Create the multiprocessing infrastructure
        number_available = multiprocessing.cpu_count()

        if int(nprocs) > number_available:
            LOG.warning("More processors requested {} than available {}".format(nprocs, number_available))
            nprocs = number_available - 1

        pool = multiprocessing.Pool(int(nprocs))
        LOG.info('Creating pool supporting {} processes...\n'.format(nprocs))

        try:
            t1 = time()
            results = pool.map_async(submit_granule, argument_dictionaries).get(9999999)
            t2 = time()
            LOG.info ("Processed {} granules using {}/{} processes in {} seconds.\n".format(total_granules, \
                    nprocs, number_available, t2-t1))
        except KeyboardInterrupt:
            LOG.error("Got exception, stopping workers and exiting.")
            pool.terminate()
            pool.join()
            sys.exit(1)

    if AGGREGATE is True:
        LOG.info("Aggregating VIIRS {}...".format(AlgorithmName))
        number_problems = aggregate_products(work_dir, EDR_collectionShortNames)

    # check new IICMO output granules
    sst_new_granules, sst_ID = h5_xdr_inventory(sstPattern, SST_GRANULE_ID_ATTR_PATH, state=sst_ID)

    LOG.debug("sst_ID.values() = {}".format(sst_ID.values()))
    LOG.debug("set(sst_ID.values()) = {}".format(set(sst_ID.values())))
    LOG.debug("sst_prior_granules = {}".format(sst_prior_granules))
    sst_granules_made = set(sst_ID.values()) - sst_prior_granules

    LOG.info('{} granules created: {}'.format(AlgorithmName,', '.join(list(sst_granules_made))))

    crashed_runs = set()
    no_output_runs = set()
    geo_problem_runs = set()
    bad_log_runs = set()

    for dicts in results:
        LOG.debug("results[{}] : {}".format(dicts['N_Granule_ID'],dicts))
        if dicts['crashed']: crashed_runs.add(dicts['N_Granule_ID'])
        if dicts['no_output'][0]: no_output_runs.add(dicts['N_Granule_ID']) 
        if dicts['geo_problem']: geo_problem_runs.add(dicts['N_Granule_ID']) 
        if dicts['bad_log']: bad_log_runs.add(dicts['N_Granule_ID']) 

    if no_output_runs:
        LOG.warning('Granules that failed to generate output: %s' % (', '.join(no_output_runs)))
    if geo_problem_runs:
        LOG.warning('Granules which had no N_Geo_Ref: %s' % ', '.join(geo_problem_runs))
    if crashed_runs:
        LOG.warning('Granules that crashed ADL: %s' % (', '.join(crashed_runs)))
    if bad_log_runs:
        LOG.warning('Granules that produced logs indicating problems: %s' % (', '.join(bad_log_runs)))
    if not sst_granules_made:
        LOG.warning('No {} HDF5 files were created'.format(AlgorithmName))

    LOG.debug('no_output_runs : {}'.format(no_output_runs))
    LOG.debug('geo_problem_runs : {}'.format(geo_problem_runs))
    LOG.debug('crashed_runs : {}'.format(crashed_runs))
    LOG.debug('bad_log_runs : {}'.format(bad_log_runs))

    return crashed_runs, no_output_runs, geo_problem_runs, bad_log_runs


def cleanup(work_dir, xml_glob, log_dir_glob, *more_dirs):
    """upon successful run, clean out work directory"""

    LOG.info("Cleaning up work directory...")

    LOG.info("Removing task xml files...")
    for fn in glob(path.join(work_dir, xml_glob)):
        LOG.debug('removing task file {}'.format(fn))
        try :
            os.unlink(fn)
        except Exception, err:
            LOG.warn( "{}".format(str(err)))

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

