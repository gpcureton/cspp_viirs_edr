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
    LOG = logging.getLogger('AerosolOpticalThicknessIP')

AlgorithmString = 'AOT'

AlgorithmName = 'Aerosol Optical Thickness IP'

GEO_collectionShortNames = [
                            'VIIRS-MOD-GEO-TC'
                          ]

SDR_collectionShortNames = [
                              'VIIRS-M1-SDR',
                              'VIIRS-M2-SDR',
                              'VIIRS-M3-SDR',
                              'VIIRS-M4-SDR',
                              'VIIRS-M5-SDR',
                              'VIIRS-M6-SDR',
                              'VIIRS-M7-SDR',
                              'VIIRS-M8-SDR',
                              'VIIRS-M9-SDR',
                              'VIIRS-M10-SDR',
                              'VIIRS-M11-SDR',
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

GridIP_collectionShortNames = []

AUX_collectionShortNames = [
                           'VIIRS-Aeros-EDR-AC-Int',
                           'VIIRS-Aeros-EDR-DQTT-Int',
                           'AOT-ANC',
                           'VIIRS-AOT-LUT',
                           'VIIRS-AOT-Sunglint-LUT',
                           'VIIRS-SusMat-EDR-DQTT-Int'
                           ]

AUX_ascTemplateFile = [
                        'VIIRS-Aeros-EDR-AC-Int_Template.asc',
                        'VIIRS-Aeros-EDR-DQTT-Int_Template.asc',
                        'AOT-ANC_Template.asc',
                        'VIIRS-AOT-LUT_Template.asc',
                        'VIIRS-AOT-Sunglint-LUT_Template.asc',
                        'VIIRS-SusMat-EDR-DQTT-Int_Template.asc'
                      ]

AUX_blobTemplateFile = [
                         'template.VIIRS-Aeros-EDR-AC-Int',
                         'template.VIIRS-Aeros-EDR-DQTT-Int',
                         'template.AOT-ANC',
                         'template.VIIRS-AOT-LUT',
                         'template.VIIRS-AOT-Sunglint-LUT',
                         'template.VIIRS-SusMat-EDR-DQTT-Int'
                       ]

AUX_Paths = [
             'luts/viirs',
             'luts/viirs',
             'luts/viirs',
             'luts/viirs',
             'luts/viirs',
             'luts/viirs'
            ]

EDR_collectionShortNames = [
                           'VIIRS-Aeros-Opt-Thick-IP',
                           'VIIRS-Aeros-EDR-GEO',
                           'VIIRS-Aeros-EDR',
                           'VIIRS-Aeros-Modl-Info-IP',
                           'VIIRS-SusMat-EDR'
                          ]

controllerName = 'ProEdrViirsAerosolController'
controllerBinary = '{}.exe'.format(controllerName)
ADL_VIIRS_AEROSOL_EDR=path.abspath(path.join(ADL_HOME, 'bin', controllerBinary))

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
# from ADL/cfg/dynamic/withMetadata/ProEdrViirsAerosolControllerLwFile.xml
xmlTemplate = """<InfTkConfig>
  <idpProcessName>ProEdrViirsAerosolController.exe</idpProcessName>
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
  <actualScans>47</actualScans>
  <previousActualScans>48</previousActualScans>
  <nextActualScans>48</nextActualScans> 
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


def setupAuxillaryFiles(Alg_objects,workDir):
    '''
    Call the generic Utils method to link in the auxillary files 
    specified in Alg_objects to the workDir directory.
    '''

    _setupAuxillaryFiles(Alg_objects,workDir)


def generate_viirs_edr_xml(work_dir, granule_seq):
    "generate XML files for VIIRS Aerosol Optical Thickness IP granule generation"
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

    h5Files                  = glob(path.join(work_dir, "*.h5"))
    aerosEdrGeoBlobFiles     = glob(path.join(work_dir, "*.VIIRS-Aeros-EDR-GEO"))
    aerosEdrBlobFiles        = glob(path.join(work_dir, "*.VIIRS-Aeros-EDR"))
    aerosModlInfoIPBlobFiles = glob(path.join(work_dir, "*.VIIRS-Aeros-Modl-Info-IP"))
    aerosOptThickIPBlobFiles = glob(path.join(work_dir, "*.VIIRS-Aeros-Opt-Thick-IP"))
    susMatEdrBlobFiles       = glob(path.join(work_dir, "*.VIIRS-SusMat-EDR"))
    ascFiles                 = glob(path.join(work_dir, "*.asc"))

    files = h5Files + aerosEdrGeoBlobFiles + aerosEdrBlobFiles + aerosModlInfoIPBlobFiles + \
            aerosOptThickIPBlobFiles + susMatEdrBlobFiles + ascFiles

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

    # Pattern for expected output in the granule working directory
    aotIpPattern = path.join(granule_output_dir, 'IVAOT*.h5')
    aotEdrPattern = path.join(granule_output_dir, 'VAOOO*.h5')
    suspMatEdrPattern = path.join(granule_output_dir, 'VSUMO*.h5')

    # prior_granules dicts contain (N_GranuleID,HDF5File) key,value pairs.
    # *ID dicts contain (HDF5File,N_GranuleID) key,value pairs.
    
    # Get the (N_GranuleID,hdfFileName) pairs for the existing Aerosol IP files
    aerosolIP_prior_granules, aotIp_ID = h5_xdr_inventory(aotIpPattern, AOT_IP_GRANULE_ID_ATTR_PATH)
    LOG.debug('Existing IVAOT granules... {}'.format(aerosolIP_prior_granules))

    aerosolIP_prior_granules = set(aerosolIP_prior_granules.keys())
    LOG.debug('Set of existing IVAOT granules... {}'.format(aerosolIP_prior_granules))

    # Get the (N_GranuleID,hdfFileName) pairs for the existing Aerosol EDR files
    aerosolEDR_prior_granules, aotEdr_ID = h5_xdr_inventory(aotEdrPattern, AOT_EDR_GRANULE_ID_ATTR_PATH)
    LOG.debug('Existing VAOOO granules... {}'.format(aerosolEDR_prior_granules))

    aerosolEDR_prior_granules = set(aerosolEDR_prior_granules.keys())
    LOG.debug('Set of existing VAOOO granules... {}'.format(aerosolEDR_prior_granules))

    # Get the (N_GranuleID,hdfFileName) pairs for the existing Suspended Matter EDR files
    suspMatEDR_prior_granules, suspMatEdr_ID = h5_xdr_inventory(suspMatEdrPattern, SUSMAT_EDR_GRANULE_ID_ATTR_PATH)
    LOG.debug('Existing VSUMO granules... {}'.format(suspMatEDR_prior_granules))

    suspMatEDR_prior_granules = set(suspMatEDR_prior_granules.keys())
    LOG.debug('Set of existing VSUMO granules... {}'.format(suspMatEDR_prior_granules))


    # Specify the command line to execute.
    cmd = [ADL_VIIRS_AEROSOL_EDR, xml]
    #cmd = ['/bin/sleep','0.2']
    #cmd = ['/usr/bin/gdb', ADL_VIIRS_AEROSOL_EDR] # for debugging with gdb...

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
        LOG.info(traceback.format_exc())
        LOG.error('{} failed on {}: {}. Continuing...'.format(cmd[0], xml, oops))
        granule_diagnostic['crashed'] = True

    t2 = time()

    LOG.info("{}({}) ran in {} seconds.".format(controllerName,granule_id,t2-t1))

    try :

        # check new IVAOT output granules
        aotIp_new_granules, aotIp_ID = h5_xdr_inventory(aotIpPattern, AOT_IP_GRANULE_ID_ATTR_PATH, state=aotIp_ID)

        if granule_id not in aotIp_new_granules:
            LOG.warning('no IVAOT HDF5 output for {}'.format(granule_id))
            granule_diagnostic['no_output'].append(True)
            granule_diagnostic['output_file'].append(None)
        else :
            LOG.info('New IVAOT granule: {}'.format(aotIp_new_granules))
            aerosolIP_granules_made = set(aotIp_ID.values()) - aerosolIP_prior_granules
            LOG.info('{} granules created: {}'.format(AlgorithmName,', '.join(list(aerosolIP_granules_made))))
            granule_diagnostic['no_output'].append(False)
            granule_diagnostic['output_file'].append(path.basename(aotIp_new_granules[granule_id]))

        # check new VAOOO output granules
        aotEdr_new_granules, aotEdr_ID = h5_xdr_inventory(aotEdrPattern, AOT_EDR_GRANULE_ID_ATTR_PATH, state=aotEdr_ID)

        if granule_id not in aotEdr_new_granules:
            LOG.warning('no VAOOO HDF5 output for {}'.format(granule_id))
            granule_diagnostic['no_output'].append(True)
            granule_diagnostic['output_file'].append(None)
        else :
            LOG.info('New VAOOO granule: {}'.format(aotEdr_new_granules))
            aerosolEDR_granules_made = set(aotEdr_ID.values()) - aerosolEDR_prior_granules
            LOG.info('Aerosol Optical Thickness EDR granules created: {}'.format(', '.join(list(aerosolEDR_granules_made))))
            granule_diagnostic['no_output'].append(False)
            granule_diagnostic['output_file'].append(path.basename(aotEdr_new_granules[granule_id]))

        # check new VSUMO output granules
        suspMatEdr_new_granules, suspMatEdr_ID = h5_xdr_inventory(suspMatEdrPattern, SUSMAT_EDR_GRANULE_ID_ATTR_PATH, state=suspMatEdr_ID)

        if granule_id not in suspMatEdr_new_granules:
            LOG.warning('no VSUMO HDF5 output for {}'.format(granule_id))
            granule_diagnostic['no_output'].append(True)
            granule_diagnostic['output_file'].append(None)
        else :
            LOG.info('New VSUMO granule: {}'.format(suspMatEdr_new_granules))
            suspMatEDR_granules_made = set(suspMatEdr_ID.values()) - suspMatEDR_prior_granules
            LOG.info('Suspended Matter EDR granules created: {}'.format(', '.join(list(suspMatEDR_granules_made))))
            granule_diagnostic['no_output'].append(False)
            granule_diagnostic['output_file'].append(path.basename(suspMatEdr_new_granules[granule_id]))

    except Exception:
        LOG.warn(traceback.format_exc())

    if compress == "True":
        LOG.info('Compress products for {}'.format(granule_id))
        repack_products(granule_output_dir, EDR_collectionShortNames)
        LOG.info('Compress products for {} complete.'.format(granule_id))


    move_products_to_work_directory(granule_output_dir)

    return granule_diagnostic


def run_xml_files(work_dir, xml_files_to_process, nprocs=1, CLEANUP="True",COMPRESS=False,AGGREGATE=False, **additional_env):
    """Run each VIIRS Cloud Mask IP xml input in sequence.
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
    aotIpPattern = path.join(work_dir, 'IVAOT*.h5')
    aotEdrPattern = path.join(work_dir, 'VAOOO*.h5')
    suspMatEdrPattern = path.join(work_dir, 'VSUMO*.h5')

    # prior_granules dicts contain (N_GranuleID,HDF5File) key,value pairs.
    # *ID dicts contain (HDF5File,N_GranuleID) key,value pairs.
    
    # Get the (N_GranuleID,hdfFileName) pairs for the existing Aerosol IP files
    aerosolIP_prior_granules, aotIp_ID = h5_xdr_inventory(aotIpPattern, AOT_IP_GRANULE_ID_ATTR_PATH)
    LOG.debug('Existing IVAOT granules... {}'.format(aerosolIP_prior_granules))

    aerosolIP_prior_granules = set(aerosolIP_prior_granules.keys())
    LOG.debug('Set of existing IVAOT granules... {}'.format(aerosolIP_prior_granules))

    # Get the (N_GranuleID,hdfFileName) pairs for the existing Aerosol EDR files
    aerosolEDR_prior_granules, aotEdr_ID = h5_xdr_inventory(aotEdrPattern, AOT_EDR_GRANULE_ID_ATTR_PATH)
    LOG.debug('Existing VAOOO granules... {}'.format(aerosolEDR_prior_granules))

    aerosolEDR_prior_granules = set(aerosolEDR_prior_granules.keys())
    LOG.debug('Set of existing VAOOO granules... {}'.format(aerosolEDR_prior_granules))

    # Get the (N_GranuleID,hdfFileName) pairs for the existing Suspended Matter EDR files
    suspMatEDR_prior_granules, suspMatEdr_ID = h5_xdr_inventory(suspMatEdrPattern, SUSMAT_EDR_GRANULE_ID_ATTR_PATH)
    LOG.debug('Existing VSUMO granules... {}'.format(suspMatEDR_prior_granules))

    suspMatEDR_prior_granules = set(suspMatEDR_prior_granules.keys())
    LOG.debug('Set of existing VSUMO granules... {}'.format(suspMatEDR_prior_granules))

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


    # check new IVAOT output granules
    aotIp_new_granules, aotIp_ID = h5_xdr_inventory(aotIpPattern, AOT_IP_GRANULE_ID_ATTR_PATH, state=aotIp_ID)

    LOG.debug("aotIp_ID.values() = {}".format(aotIp_ID.values()))
    LOG.debug("set(aotIp_ID.values()) = {}".format(set(aotIp_ID.values())))
    LOG.debug("aerosolIP_prior_granules = {}".format(aerosolIP_prior_granules))
    aerosolIP_granules_made = set(aotIp_ID.values()) - aerosolIP_prior_granules

    LOG.info('{} granules created: {}'.format(AlgorithmName,', '.join(list(aerosolIP_granules_made))))

    # check new VAOOO output granules
    aotEdr_new_granules, aotEdr_ID = h5_xdr_inventory(aotEdrPattern, AOT_EDR_GRANULE_ID_ATTR_PATH, state=aotEdr_ID)

    LOG.debug("aotEdr_ID.values() = {}".format(aotEdr_ID.values()))
    LOG.debug("set(aotEdr_ID.values()) = {}".format(set(aotEdr_ID.values())))
    LOG.debug("aerosolEDR_prior_granules = {}".format(aerosolEDR_prior_granules))
    aerosolEDR_granules_made = set(aotEdr_ID.values()) - aerosolEDR_prior_granules

    LOG.info('Aerosol Optical Thickness EDR granules created: {}'.format(', '.join(list(aerosolEDR_granules_made))))

    # check new VSUMO output granules
    suspMatEdr_new_granules, suspMatEdr_ID = h5_xdr_inventory(suspMatEdrPattern, SUSMAT_EDR_GRANULE_ID_ATTR_PATH, state=suspMatEdr_ID)

    LOG.debug("suspMatEdr_ID.values() = {}".format(suspMatEdr_ID.values()))
    LOG.debug("set(suspMatEdr_ID.values()) = {}".format(set(suspMatEdr_ID.values())))
    LOG.debug("suspMatEDR_prior_granules = {}".format(suspMatEDR_prior_granules))
    suspMatEDR_granules_made = set(suspMatEdr_ID.values()) - suspMatEDR_prior_granules

    LOG.info('Suspended Matter EDR granules created: {}'.format(', '.join(list(suspMatEdr_new_granules))))

    crashed_runs = set()
    no_output_runs = set()
    geo_problem_runs = set()
    bad_log_runs = set()

    for dicts in results:
        LOG.debug("results[{}] : {}".format(dicts['N_Granule_ID'],dicts))
        if dicts['crashed']: crashed_runs.add(dicts['N_Granule_ID'])
        if (dicts['no_output'][0] and dicts['no_output'][1] and dicts['no_output'][2]): no_output_runs.add(dicts['N_Granule_ID']) 
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
    if not aerosolIP_granules_made:
        LOG.warning('No {} HDF5 files were created'.format(AlgorithmName))
    if not aerosolEDR_granules_made:
        LOG.warning('No Aerosol Optical Thickness EDR HDF5 files were created')
    if not suspMatEDR_granules_made:
        LOG.warning('No Suspended Matter EDR HDF5 files were created')

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

