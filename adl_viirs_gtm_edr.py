#!/usr/bin/env python
# encoding: utf-8
"""
$Id$

Purpose: Run the VIIRS Ground-Track Mercator (GTM) EDR products using Raytheon ADL 4.1

Input:
    One or more HDF5 VIIRS SDR input files with matching GEO data, aggregated or single-granule.
    A work directory, typically, empty, in which to unpack the granules and generate the output.
    If the work directory specified does not exist, it will be created.

Output:
    ADL VIIRS GTM EDR blob files, .asc metadata files, and HDF5 output granules will be created.

Details:
    If you input a series of granules, the software will scan the work directory.
    It is ambiguous to provide several copies of the same granule in the work directory;
    this will result in an error abort.

Preconditions:
    Requires ADL_HOME, CSPP_SDR_HOME, CSPP_RT_ANC_CACHE_DIR, CSPP_ANC_HOME  environment variables are set.
    Requires that any needed LD_LIBRARY_PATH is set.
    Requires that DSTATICDATA is set.


Copyright (c) 2013 University of Wisconsin Regents.
Licensed under GNU GPLv3.
"""

import os, sys, logging, glob, traceback,  time, shutil
import datetime as dt
from subprocess import CalledProcessError, call
from collections import namedtuple
from multiprocessing import Pool, Lock, Value, cpu_count

import h5py

# skim and convert routines for reading .asc metadata fields of interest
from adl_asc import skim_dir, skim_dir_collections, contiguous_granule_groups, RDR_REQUIRED_KEYS

import adl_log, adl_geo_ref
import adl_anc_retrieval

import xml.etree.ElementTree as ET
from adl_common import status_line, configure_logging, get_return_code, check_env, check_and_convert_path

# ancillary search and unpacker common routines
from adl_common import sh, anc_files_needed, link_ancillary_to_work_dir, unpack, env
from adl_common import unpack_h5s
from adl_common import COMMON_LOG_CHECK_TABLE,EXTERNAL_BINARY,CSPP_RT_ANC_CACHE_DIR,CSPP_RT_ANC_PATH,DDS_PRODUCT_FILE,ADL_HOME,CSPP_RT_ANC_TILE_PATH

from adl_post_process import repack_products, aggregate_products, add_geo_attribute_to_aggregates

# controller executables
ADL_VIIRS_IXX_GTM_EDR=os.path.join(ADL_HOME, 'bin', 'ProEdrViirsIChannelImagery.exe')
ADL_VIIRS_MXX_GTM_EDR=os.path.join(ADL_HOME, 'bin', 'ProEdrViirsMChannelImagery.exe')
ADL_VIIRS_NCC_GTM_EDR=os.path.join(ADL_HOME, 'bin', 'ProEdrViirsNccImagery.exe')

LOG = logging.getLogger('adl_viirs_gtm_edr')

# keys used in metadata dictionaries
K_FILENAME = '_asc_filename'

ANCILLARY_SUB_DIR="linked_data"


# and the patterns we're looking for
ADL_VIIRS_ANC_GLOBS = tuple()
#    '*VIIRS-SDR-DNB-F-PREDICTED-LUT*',   # adl 4.1
#    '*VIIRS-SDR-F-PREDICTED-LUT*',       # adl 4.1
##  ProSdrCmnGeo
#    '*CMNGEO-PARAM-LUT_npp*',
#    '*off_Planet-Eph-ANC*',
#    '*off_USNO-PolarWander*',
#    '*CmnGeo-SAA-AC_npp*',
##    '*Terrain-Eco-ANC-Tile*',
## RDR processing
#
## GEO processing
#    '*VIIRS-SDR-GEO-DNB-PARAM-LUT_npp*',
#    '*VIIRS-SDR-GEO-IMG-PARAM-LUT_npp*',
#    '*VIIRS-SDR-GEO-MOD-PARAM-LUT_npp*',
### CAL Processing
#    '*VIIRS-SDR-DNB-DN0-LUT_npp*',
#    '*VIIRS-SDR-DNB-RVF-LUT_npp*',
#    '*VIIRS-SDR-DG-ANOMALY-DN-LIMITS-LUT_npp*',
#    '*VIIRS-SDR-DNB-STRAY-LIGHT-LUT_npp*',
#    '*VIIRS-SDR-DNB-FRAME-TO-ZONE-LUT_npp*',
##    '*VIIRS-SDR-F-LUT_npp*',               # adl 4.0
#
#    '*VIIRS-SDR-GAIN-LUT_npp*',
#    '*VIIRS-SDR-HAM-ER-LUT*',
#    '*VIIRS-SDR-RTA-ER-LUT*',
#    '*VIIRS-SDR-OBC-ER-LUT_npp*',
#    '*VIIRS-SDR-OBC-RR-LUT_npp*',
#    '*VIIRS-SDR-EBBT-LUT_npp*',
#    '*VIIRS-SDR-TELE-COEFFS-LUT_npp*',
#    '*VIIRS-SDR-SOLAR-IRAD-LUT_npp*',
#    '*VIIRS-SDR-RSR-LUT_npp*',
#    '*VIIRS-SDR-OBS-TO-PIXELS-LUT_npp*',
#    '*VIIRS-SOLAR-DIFF-VOLT-LUT_npp*',
#    '*VIIRS-SDR-RADIOMETRIC-PARAM-LUT_npp*',
#    '*VIIRS-SDR-QA-LUT_npp*',
#    '*VIIRS-SDR-EMISSIVE-LUT_npp*',
#    '*VIIRS-SDR-REFLECTIVE-LUT_npp*',
#    '*VIIRS-SDR-RVF-LUT_npp*',
#    '*VIIRS-SDR-BB-TEMP-COEFFS-LUT_npp*',
#    '*VIIRS-SDR-DNB-C-COEFFS-LUT_npp*',
#    '*VIIRS-SDR-DELTA-C-LUT_npp*',
#    '*VIIRS-SDR-COEFF-A-LUT_npp*',
#    '*VIIRS-SDR-COEFF-B-LUT_npp*',
#        '*TLE-AUX*'
#                    )


# ADL_VIIRS_GEO_PRODUCT_SHORTNAMES = [
#     'VIIRS-MOD-GEO-TC',
#     'VIIRS-IMG-GEO-TC',
#     'VIIRS-DNB-GEO',
#    ]


# OPTIONAL_GEO_PRODUCTS=[
#     'VIIRS-IMG-GEO',
#     'VIIRS-MOD-GEO'
# ]


# ADL_VIIRS_SDR_PRODUCT_SHORTNAMES = [

#     'VIIRS-I1-SDR',
#     'VIIRS-I2-SDR',
#     'VIIRS-I3-SDR',
#     'VIIRS-I4-SDR',
#     'VIIRS-I5-SDR',
#     'VIIRS-M1-SDR',
#     'VIIRS-M2-SDR',
#     'VIIRS-M3-SDR',
#     'VIIRS-M4-SDR',
#     'VIIRS-M5-SDR',
#     'VIIRS-M6-SDR',
#     'VIIRS-M7-SDR',
#     'VIIRS-M8-SDR',
#     'VIIRS-M9-SDR',
#     'VIIRS-M10-SDR',
#     'VIIRS-M11-SDR',
#     'VIIRS-M12-SDR',
#     'VIIRS-M13-SDR',
#     'VIIRS-M14-SDR',
#     'VIIRS-M15-SDR',
#     'VIIRS-M16-SDR',
#     'VIIRS-DNB-SDR'

# ]

# ADL_VIIRS_SDR_intermediate_SHORTNAMES = [
#     'VIIRS-MOD-UNAGG-GEO',  # no nagg
#     'VIIRS-DualGain-Cal-IP',    # no nagg
#     'VIIRS-OBC-IP', # no nagg

#     'VIIRS-IMG-RGEO',
#     'VIIRS-MOD-RGEO',
#     'VIIRS-I1-FSDR',
#     'VIIRS-I2-FSDR',
#     'VIIRS-I3-FSDR',
#     'VIIRS-I4-FSDR',
#     'VIIRS-I5-FSDR',
#     'VIIRS-M1-FSDR',
#     'VIIRS-M2-FSDR',
#     'VIIRS-M3-FSDR',
#     'VIIRS-M4-FSDR',
#     'VIIRS-M5-FSDR',
#     'VIIRS-M6-FSDR',
#     'VIIRS-M7-FSDR',
#     'VIIRS-M8-FSDR',
#     'VIIRS-M9-FSDR',
#     'VIIRS-M10-FSDR',
#     'VIIRS-M11-FSDR',
#     'VIIRS-M12-FSDR',
#     'VIIRS-M14-FSDR',
#     'VIIRS-M15-FSDR',
#     'VIIRS-M16-FSDR',
#     'VIIRS-MOD-RGEO',
#     'VIIRS-MOD-RGEO-TC'
# ]


CHECK_REQUIRED_KEYS = ['N_Granule_ID', 'N_Collection_Short_Name']

OBSERVE_TIME = 'ObservedStartTime'
# table of NPP short names to data product ids

# PRODUCTID_2_SHORTNAME= dict()
# SHORTNAME_2_PRODUCTID = dict()


# WORK_DIR: directory that we unpack the input data into and accumulate final output to
# WORK_SUBDIR: output directory written to by each granule+kind task instance

XML_TMPL_VIIRS_MXX_GTM_EDR = """<InfTkConfig>
  <idpProcessName>ProEdrViirsMChannelImagery.exe</idpProcessName>
  <siSoftwareId></siSoftwareId>
  <isRestart>FALSE</isRestart>
  <useExtSiWriter>FALSE</useExtSiWriter>
  <debugLogLevel>NORMAL</debugLogLevel>
  <debugLevel>DBG_HIGH</debugLevel>
  <dbgDest>D_FILE</dbgDest>
  <enablePerf>FALSE</enablePerf>
  <perfPath>${WORK_DIR}</perfPath>
  <dbgPath>${WORK_DIR}</dbgPath>
  <initData>
     <domain>OPS</domain>
     <subDomain>SUBDOMAIN</subDomain>
     <startMode>INF_STARTMODE_COLD</startMode>
     <executionMode>INF_EXEMODE_PRIMARY</executionMode>
     <healthTimeoutPeriod>30</healthTimeoutPeriod>
  </initData>
  <lockinMem>FALSE</lockinMem>
  <rootDir>${WORK_SUBDIR}/log</rootDir>
  <inputPath>${WORK_DIR}:${WORK_SUBDIR}</inputPath>
  <outputPath>${WORK_SUBDIR}</outputPath>
  <dataStartIET>0</dataStartIET>
  <dataEndIET>0</dataEndIET>
  <actualScans>0</actualScans>
  <previousActualScans>0</previousActualScans>
  <nextActualScans>0</nextActualScans>
  <usingMetadata>TRUE</usingMetadata>
  <configGuideName>ProEdrViirsGtmMChannelImagery_GuideList.cfg</configGuideName>

  <task>
    <taskType>EDR</taskType>
    <taskDetails1>%(N_Granule_ID)s</taskDetails1>
    <taskDetails2>%(N_Granule_Version)s</taskDetails2>
    <taskDetails3>NPP</taskDetails3>
    <taskDetails4>VIIRS</taskDetails4>
  </task>

</InfTkConfig>
"""


XML_TMPL_VIIRS_IXX_GTM_EDR = """<InfTkConfig>
  <idpProcessName>ProEdrViirsIChannelImagery.exe</idpProcessName>
  <siSoftwareId></siSoftwareId>
  <isRestart>FALSE</isRestart>
  <useExtSiWriter>FALSE</useExtSiWriter>
  <debugLogLevel>NORMAL</debugLogLevel>
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
  <rootDir>${WORK_SUBDIR}/log</rootDir>
  <inputPath>${WORK_DIR}:${WORK_SUBDIR}</inputPath>
  <outputPath>${WORK_SUBDIR}</outputPath>
  <dataStartIET>0</dataStartIET>
  <dataEndIET>0</dataEndIET>
  <actualScans>0</actualScans>
  <previousActualScans>0</previousActualScans>
  <nextActualScans>0</nextActualScans>
  <usingMetadata>TRUE</usingMetadata>
  <configGuideName>ProEdrViirsGtmIChannelImagery_GuideList.cfg</configGuideName>

  <task>
    <taskType>EDR</taskType>
    <taskDetails1>%(N_Granule_ID)s</taskDetails1>
    <taskDetails2>%(N_Granule_Version)s</taskDetails2>
    <taskDetails3>NPP</taskDetails3>
    <taskDetails4>VIIRS</taskDetails4>
  </task>

</InfTkConfig>
"""


XML_TMPL_VIIRS_NCC_GTM_EDR = """<InfTkConfig>
  <idpProcessName>ProEdrViirsNccImagery.exe</idpProcessName>
  <siSoftwareId></siSoftwareId>
  <isRestart>FALSE</isRestart>
  <useExtSiWriter>FALSE</useExtSiWriter>
  <debugLogLevel>NORMAL</debugLogLevel>
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
  <rootDir>${WORK_SUBDIR}/log</rootDir>
  <inputPath>${WORK_DIR}:${WORK_SUBDIR}</inputPath>
  <outputPath>${WORK_SUBDIR}</outputPath>
  <dataStartIET>0</dataStartIET>
  <dataEndIET>0</dataEndIET>
  <actualScans>0</actualScans>
  <previousActualScans>0</previousActualScans>
  <nextActualScans>0</nextActualScans>
  <usingMetadata>TRUE</usingMetadata>
  <configGuideName>ProEdrViirsGtmNccImagery_GuideList.cfg</configGuideName>

  <task>
    <taskType>EDR</taskType>
    <taskDetails1>%(N_Granule_ID)s</taskDetails1>
    <taskDetails2>%(N_Granule_Version)s</taskDetails2>
    <taskDetails3>NPP</taskDetails3>
    <taskDetails4>VIIRS</taskDetails4>
  </task>

</InfTkConfig>
"""

# create a named tuple class holding the information we want for each group of SDRs
# geo_cn: geolocation collection name
# sdr_cn: sdr collection name
guidebook_info = namedtuple('guidebook_info', 'sdr_cns geo_cn template exe')

# note that for night time we only get M7,8,10,12,13,14,15,16, and I4,5
GTM_GUIDEBOOK = {
    'IXX': guidebook_info(set('VIIRS-I%d-SDR' % b for b in (4,5)), 'VIIRS-IMG-GEO', XML_TMPL_VIIRS_IXX_GTM_EDR, ADL_VIIRS_IXX_GTM_EDR),
    'MXX': guidebook_info(set('VIIRS-M%02d-SDR' % b for b in (7,8,10,12,13,14,15,16)), 'VIIRS-IMG-GEO', XML_TMPL_VIIRS_MXX_GTM_EDR, ADL_VIIRS_MXX_GTM_EDR),
    'NCC': guidebook_info(set('VIIRS-DNB-SDR'), 'VIIRS-DNB-GEO', XML_TMPL_VIIRS_NCC_GTM_EDR, ADL_VIIRS_NCC_GTM_EDR) # FIXME review this for accuracy
}


def _trim_geo_granules(gran_dict_seq):
    "sort granule sequence by (N_Granule_ID, N_Granule_Version); eliminate redundant old versions and output sequence"
    key = lambda g: (g['N_Granule_ID'], g['N_Granule_Version'])
    lst = list(gran_dict_seq)
    lst.sort(key=key)
    # go through sorted list, grabbing latest version of each granule
    dct = dict((g['N_Granule_ID'], g) for g in lst)
    return sorted(dct.values(), key=key)


def sift_metadata_for_viirs_gtm_edr(work_dir='.'):
    """
    search through the ASC metadata in a work directory, identifying processing candidate granules
    look for SDR + GEO and yield metadata dictionaries
    for each GEO granule available, sorted by N_Granule_Version
      look for corresponding SDR granules
      if non available, report error
      yield (kind, geo_granule_dictionary, sdr_N_Collection_Short_Names_seq)
    consults guidebook to know which collections are useful for a given kind of product
    """
    LOG.info("checking for VIIRS SDR and GEO blobs")

    # get {collection-name: list-of-dictionaries, ... } of all qualified metadata from unpacked blobs in work_dir
    meta = skim_dir_collections(work_dir, required_keys=CHECK_REQUIRED_KEYS)

    # FUTURE: reduce the depth of iteration
    # set of granules with available geolocation
    for kind, G in GTM_GUIDEBOOK.items():
        # list of available geo products for this group
        geo_granules = _trim_geo_granules( meta[G.geo_cn] )

        # check if we have at least one SDR collection and one GEO collection for this granule
        # if so, yield it
        for geo_granule in geo_granules:
            geo_gran_id = geo_granule['N_Granule_ID']
            geo_gran_ver = geo_granule['N_Granule_Version']

            for sdr_cn in [cn for cn in G.sdr_cns]:              
                sdr_grans = [g for g in meta[sdr_cn] if ((g['N_Granule_ID']==geo_gran_id) and (g['N_Granule_Version']==geo_gran_ver))]
                if not sdr_grans: 
                    LOG.warning('no SDR products found for %s:%s-v%d' % (kind, geo_gran_id, geo_gran_ver))
                else:
                    sdr_collections = [x['N_Collection_Short_Name'] for x in sdr_grans]
                    LOG.debug('found SDR collections %s for %s:%s-v%d' % (repr(sdr_collections), kind, geo_gran_id, geo_gran_ver))
                    yield (kind, geo_granule, sdr_collections)


def generate_gtm_edr_xml(kind, gran, work_dir):
    """
    writes XML file using template
    returns "granule_name.xml"
    :param work_dir: directory to write XML to
    :param kind: IXX, MXX, NCC string
    :param gran: granule dictionary as returned from sift_metadata_for_viirs_gtm_edr
    """
    name = gran['N_Granule_ID']
    xml_tmpl = GTM_GUIDEBOOK[kind].template
    fnxml = ('edr_viirs_gtm_%s_%s.xml' % (name, gran['N_Collection_Short_Name']))
    LOG.debug('writing XML file %r' % fnxml)
    with open(os.path.join(work_dir, fnxml), 'wt') as fpxml:
      fpxml.write(xml_tmpl % gran)
    status_line('Created ADL controller XML %s for %s:%s' % (fnxml, kind, name))
    return fnxml


# Look through new log files for completed messages
def check_logs_for_run(work_dir, pid, xml):
    """
    Find the log file 
    Look for success
    Return True if found
    Display log message and hint if problem occurred
    """
    # retrieve exe name and log path from lw file
    logDir = os.path.join(work_dir, "log")
    logExpression = "*" + str(pid) + "*.lo*"
    
    files = glob.glob(os.path.join(logDir, logExpression))
    status_line("Checking "+str(len(files))+" log files for errors"+logDir+" exp "+logExpression)
   
    n_err = 0
    err_files = set()
    for log_file in files:
        LOG.info("Checking Log file " + log_file + " for errors.")
        count = adl_log.scan_log_file(COMMON_LOG_CHECK_TABLE, log_file)
        n_err += count
        if count > 0:
            err_files.add(log_file)

    if n_err == 0:
        status_line("Processing of file: " + xml + " Completed successfully" )
        return True
    else:
        status_line("Processing of file: " + xml + " Completed unsuccessfully, Look at previous message" )
        LOG.debug("Log files with errors: " + ', '.join(err_files))
        return False


def transfer_gtm_edr_output(work_dir, work_subdir, kind, gran, sdr_collections):
    """
    examine work_subdir for products, based on sdr_collections that were available as input
    transfer products back from work_subdir to work_dir
    return product_filenames sequence, and transfer error sequence
    """
    products = []
    errors = []
    # FIXME: this is just a first wag at it and could stand to be more discriminating
    for h5path in glob.glob(os.path.join(work_subdir, '*.h5')):
        LOG.debug('transferring output %s' % h5path)
        h5filename = os.path.split(h5path)[-1]
        h5out = os.path.join(work_dir, h5filename)
        try:
            shutil.move(h5path, h5out)
            products.append(h5filename)
        except:
            errors.append('%s would not transfer' % h5path)
    return products, errors


task_input = namedtuple('task_input', 'kind granule sdr_collections work_dir env_dict')
task_output = namedtuple('task_output', 'kind granule_id product_filenames error_list')


def task_gtm_edr(task_in):
    """
    process a single task, returning a task_output tuple
    this is suitable for spinning off to a subprocess using multiprocessing.Pool
    """
    kind, gran, sdr_collections, work_dir, additional_env = task_in
    G = GTM_GUIDEBOOK[kind]

    granule_id = gran['N_Granule_ID']
    # list of error strings explaining what went wrong - we do this as a list of strings so that it can
    # 1. cross process boundary back to the master process
    # 2. properly be written to the screen in a final report
    errors = []

    # create granule_collection work-subdir which the controller will output to
    work_subdir = os.path.join(work_dir, 'GTM_%s_%s' % (kind, gran['N_Granule_ID']))  # e.g. GTM_IXX_NPP987654542516
    LOG.debug('granule %s:%s will run in %s' % (kind, gran['N_Granule_ID'], work_subdir))
    if os.path.isdir(work_subdir):
        LOG.error('directory %s already exists, re-use of work directories is discouraged' % work_subdir)
        return task_output(kind, granule_id, [], ['invalid work directory'])
    else:
        os.makedirs(os.path.join(work_subdir,'log'))

    # generate XML into work subdirectory
    xml_filename = generate_gtm_edr_xml(kind, gran, work_subdir)

    # run XML controller
    exe = G.exe
    cmd = [exe, xml_filename]
    local_env = {'WORK_SUBDIR': work_subdir}
    local_env.update(additional_env)

    status_line('Executing %s' % repr(cmd))
    LOG.debug('additional environment variables: %s' % repr(local_env))
    try:
        pid = sh(cmd, env=env(**local_env), cwd=work_subdir)
        LOG.debug("%r ran as pid %d" % (cmd, pid))
        ran_ok = check_logs_for_run(work_subdir, pid, xml_filename)

    except CalledProcessError as oops:
        pid = getattr(oops, 'pid', None)
        errors.append('process crashed')
        ran_ok = check_logs_for_run(work_subdir, pid, xml_filename)
        if not ran_ok:
            errors.append('log file problem')
        LOG.debug(traceback.format_exc())
        LOG.error('ProSdrViirsController.exe failed on %r: %r. Continuing...' % (xml_filename, oops))

    if not ran_ok:
        errors.append('logs were not error-free')

    # link the output from the work_subdir to the work_dir
    product_filenames, transfer_errors = transfer_gtm_edr_output(work_dir, work_subdir, kind, gran, sdr_collections)
    errors += list(transfer_errors)

    # if everything ran OK, clean up the intermediate stuff in our subdir
    if not errors: 
        LOG.debug('cleaning up %s, no errors' % work_subdir)
        shutil.rmtree(work_subdir)

    return task_output(kind, granule_id, product_filenames, errors)


def herd_viirs_gtm_edr_tasks(work_dir, nprocs=1, **additional_env):
    tasks = []
    # find all the things we want to do and build tasks for them
    for kind, geo_granule, sdr_collections in sift_metadata_for_viirs_gtm_edr(work_dir):
        tasks.append(task_input(kind, geo_granule, sdr_collections, work_dir, additional_env))
    parallel = Pool( int(nprocs) )
    try:
        results = parallel.map(task_gtm_edr, tasks)
    except (KeyboardInterrupt, SystemError) as ejectionseat:
        # note that we're depending on register_sigterm having been called for SystemError on SIGTERM
        LOG.warning('external termination detected, aborting subprocesses')
        parallel.terminate()
    # line up our tasks and run them with the processing pool
    # pick up task_output tuples 
    return results


def setup_directories(work_dir,anc_dir):
    """Create the working directory and a subdirectory for the logs
    :param work_dir: directory which we'll be creating work files in
    :param anc_dir: ancillary directory we'll be linking in
    """
    if not os.path.isdir(work_dir):
        LOG.info('Creating directory %s' % work_dir)
        os.makedirs(work_dir)

    log_dir = os.path.join(work_dir, 'log')
    if not os.path.isdir(log_dir):
        LOG.info('Creating log directory %s' % log_dir)
        os.makedirs(log_dir)

    if not os.path.exists( anc_dir ) :
        os.mkdir(anc_dir)


_registered_sigterm = False


def _sigterm(sig, frame):
    raise SystemError(sig)


def register_sigterm():
    global _registered_sigterm
    if not _registered_sigterm:
        import signal
        signal.signal(signal.SIGTERM, _sigterm)
        _registered_sigterm = True


def read_N_Geo_Ref(h5path):
    """
    Read the N_Geo_Ref attribute from the file if it exists, else return None
    :param h5path: hdf5 pathname
    :return: N_Geo_Ref filename, or None
    """
    h5 = h5py.File(h5path, 'r')
    zult = getattr(h5, 'N_GEO_Ref', None)
    h5.close()
    return zult


def input_list_including_geo(pathnames):
    """
    :param pathnames: sequence of pathnames we plan to unpack
    :return: set of pathnames, including geo filenames corresponding
    """
    zult = set()
    for path in pathnames:
        if not h5py.is_hdf5(path):
            LOG.warning('%s is not an HDF5 file, ignoring' % path)
            continue
        zult.add(path)
        geo = read_N_Geo_Ref(path)
        if geo is not None:
            dn, fn = os.path.split(path)
            LOG.info("adding %s as companion to %s" % (geo, fn))
            geo_path = os.path.join(dn, geo)
            if not h5py.is_hdf5(geo_path):
                LOG.warning('expected to find %s as an HDF5 file; may not be able to process %s' % (geo_path, fn))
                continue
            zult.add(geo_path)
    return zult


def viirs_gtm_edr(work_dir, h5_paths, nprocs = 1, compress=False, aggregate=False, allow_cache_update=True):
    """
    given a work directory and a series of hdf5 SDR and GEO paths
    make work directory
    unpack SDRs and GEOs into the workspace
    skim the metadata for tasks
    run the tasks, collecting task_output tuples
    clean up the work directory if errors weren't present
    report on the final outcome
    return a result code to pass back to the shell (0 for success, nonzero for error)
    :param work_dir: directory to work in
    :param h5_paths: SDR and GEO hdf5 path sequence
    :param nprocs: number of processors to use, default 1
    """
    check_env(work_dir)

    register_sigterm()

    anc_dir = os.path.join(work_dir, ANCILLARY_SUB_DIR)
    setup_directories(work_dir, anc_dir)

    status_line("Unpack the supplied inputs")
    h5_paths = list(input_list_including_geo(h5_paths))
    LOG.info('final list of inputs: %s' % repr(h5_paths))
    error_count = unpack_h5s(work_dir, h5_paths)
    # FIXME: this should discover the GEOs from the SDRs if only SDRs are provided
    # FIXME: compression
    # FIXME: aggregation
    # FIXME: cache update if we use any ancillary

    results = herd_viirs_gtm_edr_tasks(work_dir,
                                       nprocs=nprocs,
                                       LINKED_ANCILLARY=ANCILLARY_SUB_DIR,
                                       ADL_HOME=ADL_HOME,
                                       CSPP_RT_ANC_TILE_PATH=CSPP_RT_ANC_TILE_PATH
                                       )
    for products, errors in results:
        error_count += len(errors)
    return error_count




def main():
    """ Run Viirs GTM EDR processing
    """
    import argparse
    desc = """Build VIIRS GTM EDR work directory and run VIIRS GTM EDR."""
    parser = argparse.ArgumentParser(description = desc)
    parser.add_argument('-t', '--test',
                    action="store_true", default=False, help="run self-tests")
    parser.add_argument('-W', '--work-dir', metavar='work_dir', default='.',
                    help='work directory which all activity will occur in, defaults to current dir')

    parser.add_argument('-d', '--debug',
                    action="store_true", default=False, help="always retain intermediate files")

    parser.add_argument('-z', '--zip',
                    action="store_true", default=False, help="compress products with h5repack zip compression")

    parser.add_argument('-g', '--geo',
                    action="store_true", default=False, help="Retain terain un-correct GEO products")


    parser.add_argument('-a', '--aggregate',
                    action="store_true", default=False, help="aggregate products with nagg")

    parser.add_argument('-l', '--local',
                    action="store_true", default=False, help="disable download of remote ancillary data to cache")

    parser.add_argument('-p', '--processor',
                    type=int, default=1, help="Number of processors to use for band processing")

    parser.add_argument('-v', '--verbosity', action="count", default=0,
                    help='each occurrence increases verbosity 1 level through ERROR-WARNING-INFO-DEBUG')

    parser.add_argument('filenames', metavar='filename', type=str, nargs='+',
                   help='HDF5 VIIRS RDR file/s to process')

    args = parser.parse_args()

    levels = [logging.ERROR, logging.WARN, logging.INFO, logging.DEBUG]
    level = levels[args.verbosity if args.verbosity<4 else 3]

    do_cleanup = True

    work_dir = check_and_convert_path("WORK_DIR",args.work_dir)
    d = dt.date.today()
    timestamp = d.isoformat()
    log_name = "viirs_gtm_edr.%s.log" % timestamp
    logfile = os.path.join(work_dir, log_name)
    configure_logging(level, FILE=logfile)
    if args.debug == True:
        do_cleanup = False

    LOG.debug("Clean up: "+str(do_cleanup))
    LOG.info('CSPP execution work directory is %r' % work_dir)

    # if args.test:
    #     check_env(work_dir)
    #     grans = _test_sdr_granules(work_dir)
    #     if grans:
    #         LOG.debug('building XML files')
    #
    #     sys.exit(2)

    if not args:
        parser.print_help()
        return 1

    nprocs = args.processor
    if nprocs <= 0:
        nprocs = cpu_count()
        LOG.info('using nprocs=%d' % nprocs)

    rc = viirs_gtm_edr(work_dir, args.filenames, nprocs=nprocs,
                       compress=args.zip, aggregate=args.aggregate,
                       allow_cache_update=not args.local)

    if rc == 0 and not args.debug:
        os.remove(logfile)

    return rc


if __name__ == '__main__':
    sys.exit(main())
