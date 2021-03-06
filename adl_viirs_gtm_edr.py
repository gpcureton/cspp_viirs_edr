#!/usr/bin/env python
# encoding: utf-8
"""
$Id$

Purpose: Run the VIIRS Ground-Track Mercator (GTM) EDR products using Raytheon ADL 4.1

Input:
    A input directory with intermediate and end products from VIIRS SDR in BLOB+ASC form.
    This is typically generated in CSPP by using --edr parameter.
    A work directory, typically, empty, in which to unpack the granules and generate the output.
    If the work directory specified does not exist, it will be created.

Output:
    ADL VIIRS GTM EDR blob files, .asc metadata files, and HDF5 output granules will be created.

Details:
    If you input a series of granules, the software will scan the work directory.
    It is ambiguous to provide several copies of the same granule in the work directory;
    this will result in an error abort.

Preconditions:
    Requires ADL_HOME, CSPP_EDR_HOME, CSPP_RT_ANC_CACHE_DIR, CSPP_ANC_HOME  environment variables are set.
    Requires that any needed LD_LIBRARY_PATH is set.
    Requires that DSTATICDATA is set.


Copyright (c) 2013 University of Wisconsin Regents.
Licensed under GNU GPLv3.
"""

import os
import sys
import logging
import glob
import traceback
import shutil
import re
import datetime as dt
from subprocess import CalledProcessError
from collections import namedtuple
from multiprocessing import Pool, cpu_count
from pprint import pformat
import h5py

# skim and convert routines for reading .asc metadata fields of interest
from adl_asc import skim_dir, skim_dir_collections, contiguous_granule_groups, RDR_REQUIRED_KEYS, K_FILENAME

import adl_log
import adl_geo_ref
# import adl_anc_retrieval

from adl_common import status_line, configure_logging, get_return_code,  check_and_convert_path

# ancillary search and unpacker common routines
from adl_common import sh, link_ancillary_to_work_dir, env
from adl_common import COMMON_LOG_CHECK_TABLE, EXTERNAL_BINARY, CSPP_RT_ANC_CACHE_DIR, CSPP_RT_ANC_PATH, DDS_PRODUCT_FILE, ADL_HOME, CSPP_RT_ANC_TILE_PATH, ADL_VARS

from adl_post_process import repack_products, aggregate_products, add_geo_attribute_to_aggregates

# controller executables
ADL_VIIRS_IXX_GTM_EDR = os.path.join(ADL_HOME, 'bin', 'ProEdrViirsIChannelImagery.exe')
ADL_VIIRS_MXX_GTM_EDR = os.path.join(ADL_HOME, 'bin', 'ProEdrViirsMChannelImagery.exe')
ADL_VIIRS_NCC_GTM_EDR = os.path.join(ADL_HOME, 'bin', 'ProEdrViirsNccImagery.exe')

LOG = logging.getLogger('adl_viirs_gtm_edr')

ANCILLARY_SUB_DIR = "linked_data"

CHECK_REQUIRED_KEYS = ['N_Granule_ID', 'N_Collection_Short_Name']

# FUTURE: PRO_CROSSGRAN_FAIL may go away if we obtain -FSDR products from viirs_sdr.sh --edr
GTM_EDR_LOG_CHECK_TABLE = [('PRO_CROSSGRAN_FAIL', "Cross Granule dependency failure, more input may be needed?"),
                           ('INF_STATUSTYPE_TASK_INPUTNOTAVAIL', "Missing input?")]

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
  <inputPath>${WORK_DIR}:${WORK_SUBDIR}:${LINKED_ANCILLARY}:${INPUT_DIR}</inputPath>
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
  <inputPath>${WORK_DIR}:${WORK_SUBDIR}:${LINKED_ANCILLARY}:${INPUT_DIR}</inputPath>
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
  <inputPath>${WORK_DIR}:${WORK_SUBDIR}:${LINKED_ANCILLARY}:${INPUT_DIR}</inputPath>
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
guidebook_info = namedtuple('guidebook_info', 'sdr_cns edr_cns geo_cn template exe anc')

# note that for night time we only get M7,8,10,12,13,14,15,16, and I4,5
# guidebook tells us what to expect and how to deal with it
# Ref OAD for VIIRS GTM Imagery
# FUTURE: get ancillary requirements by reading ${ADL_HOME}/cfg/ProEdrViirs{IChannel,MChannel,Ncc}Imagery_CFG.xml
# FUTURE: check for FSDR blobs and use them if available?
# FUTURE: include GEO EDR output products e.g. , 'VIIRS-NCC-EDR-GEO'
# note that all of them need TLE and PolarWander anc!
GTM_GUIDEBOOK = {
    'IXX': guidebook_info(sdr_cns=['VIIRS-I%d-SDR' % b for b in (1, 2, 3, 4, 5)],
                          edr_cns=['VIIRS-I%d-IMG-EDR' % b for b in (1, 2, 3, 4, 5)],
                          geo_cn='VIIRS-IMG-GEO',
                          template=XML_TMPL_VIIRS_IXX_GTM_EDR,
                          exe=ADL_VIIRS_IXX_GTM_EDR,
                          anc=('CmnGeo-SAA-AC-Int',
                               'CMNGEO-PARAM-LUT',
                               'Planet-Eph-ANC')),
    'MXX': guidebook_info(sdr_cns=['VIIRS-M%02d-SDR' % b for b in (1, 4, 9, 14, 15, 16)],
                          edr_cns=['VIIRS-M%s-EDR' % q for q in ['1ST', '2ND', '3RD', '4TH', '5TH', '6TH']],
                          geo_cn='VIIRS-MOD-GEO',
                          template=XML_TMPL_VIIRS_MXX_GTM_EDR,
                          exe=ADL_VIIRS_MXX_GTM_EDR,
                          anc=('CmnGeo-SAA-AC-Int',
                               'CMNGEO-PARAM-LUT',
                               'Planet-Eph-ANC')),
    'NCC': guidebook_info(sdr_cns=['VIIRS-DNB-SDR'],
                          edr_cns=['VIIRS-NCC-EDR'],
                          geo_cn='VIIRS-DNB-GEO',
                          template=XML_TMPL_VIIRS_NCC_GTM_EDR,
                          exe=ADL_VIIRS_NCC_GTM_EDR,
                          anc=('CmnGeo-SAA-AC-Int',
                               'CMNGEO-PARAM-LUT',
                               'Planet-Eph-ANC',
                               'VIIRS-NCC-EDR-AC-Int',
                               'VIIRS-LUN-Phase-LUT',
                               'VIIRS-Sol-BRDF-LUT',
                               'VIIRS-Lun-BRDF-LUT',
                               'VIIRS-Ga-Val-Vs-Scene-Lun-Elev-LUT',
                               'VIIRS-Ga-Val-Vs-Scene-Sol-Elev-LUT'))
}


EDR_PRODUCTS_AGG = [
    cn for _,x in GTM_GUIDEBOOK.items() for cn in x.edr_cns
]

EDR_GEO_PRODUCTS_AGG = [
    'VIIRS-NCC-EDR-GEO',
    'VIIRS-IMG-GTM-EDR-GEO',
    'VIIRS-MOD-GTM-EDR-GEO'
]


def _geo_guide():
    """
    yield additions to GEO_GUIDE in adl_geo_ref.py
    This tells adl_geo_ref.geo_ref() how to match GEO with EDR files
    """
    for n, cn in enumerate(GTM_GUIDEBOOK['MXX'].edr_cns):
        yield ('VM%02d' % (n+1), 'GMGTO', r'Data_Products/%s/%s_Gran_0/N_Reference_ID' % (cn,cn))
    yield ('GMGTO', None, r'Data_Products/VIIRS-MOD-GTM-EDR-GEO/VIIRS-MOD-GTM-EDR-GEO_Gran_0/N_Reference_ID')
    for n, cn in enumerate(GTM_GUIDEBOOK['IXX'].edr_cns):
        yield ('VI%dBO' % (n+1), 'GIGTO', r'Data_Products/%s/%s_Gran_0/N_Reference_ID' % (cn,cn))
    yield ('GIGTO', None, r'Data_Products/VIIRS-IMG-GTM-EDR-GEO/VIIRS-IMG-GTM-EDR-GEO_Gran_0/N_Reference_ID')
    yield ('VNCCO', 'GNCCO', r'Data_Products/VIIRS-NCC-EDR/VIIRS-NCC-EDR_Gran_0/N_Reference_ID')
    yield ('GNCCO', None, r'Data_Products/VIIRS-NCC-EDR-GEO/VIIRS-NCC-EDR-GEO_Gran_0/N_Reference_ID')


GEO_GUIDE_PATCHED = False


def _patch_geo_guide():
    """
    FUTURE: integrate this into adl_geo_ref.py
    This adjusts for the difficulty of the Collection Name being 1-2-3-4-5-6 vs 1ST-2ND-3RD-4TH-5TH-6TH (thanks guys)
    on M collection. So we patch the lookup table and the filename regex for our use.
    :return:
    """
    global GEO_GUIDE_PATCHED
    if GEO_GUIDE_PATCHED:
        return
    LOG.debug('patching adl_geo_ref to permit operating on GTM output')
    new_guide_info = list(_geo_guide())
    LOG.debug(pformat(new_guide_info))
    adl_geo_ref.GEO_GUIDE += new_guide_info
    # previously kind would be "SVI", band would be 01-05, as an example; this is OK for SDR products
    # in this case we eat the whole prefix e.g. VI1BO as kind, and empty-string for band
    # we do this because the mapping between collection name and filename for VM* is not straightforward
    # VM03O => VIIRS-M3RD-EDR for instance
    # so we build an exhaustive table (_geo_guide) from GTM_GUIDEBOOK, push that into the adl_geo_ref table,
    # and patch the pathname regex to allow VM03O and VI3BO to scan successfully
    adl_geo_ref.RE_NPP = re.compile('(?P<kind>[A-Z0-9]+)(?P<band>[0-9]*)_(?P<sat>[A-Za-z0-9]+)_d(?P<date>\d+)'
                                    '_t(?P<start_time>\d+)_e(?P<end_time>\d+)_b(?P<orbit>\d+)_c(?P<created_time>\d+)'
                                    '_(?P<site>[a-zA-Z0-9]+)_(?P<domain>[a-zA-Z0-9]+)\.h5')
    GEO_GUIDE_PATCHED = True


def geo_ref(path):
    _patch_geo_guide()
    return adl_geo_ref.geo_ref(path)


def write_geo_ref(path):
    _patch_geo_guide()
    return adl_geo_ref.write_geo_ref(path)


def _trim_geo_granules(gran_dict_seq):
    """sort granule sequence by (N_Granule_ID, N_Granule_Version); eliminate redundant old versions and output sequence
    :param gran_dict_seq: sequence of granule metadata dictionaries
    """
    key = lambda g: (g['N_Granule_ID'], g['N_Granule_Version'])
    lst = list(gran_dict_seq)
    lst.sort(key=key)
    # go through sorted list, grabbing latest version of each granule
    dct = dict((g['N_Granule_ID'], g) for g in lst)
    return sorted(dct.values(), key=key)


def _crossgran_filter(geo_granules, n_crossgran=1):
    """
    given a sequence of geo granule metadata dictionaries,
    yield a sequence of geo granules which satisfy cross-granule dependencies
    in this case, we skip the first and last granule of each contiguous group
    this filter should be removed when we have cross-granule dependencies
    :param n_crossgran: number of cross-granules to check for, e.g. 1 implies +1/-1 granules are needed
    :param geo_granules: sequence of geo_granules to filter
    :return: filtered geo granules, in time order, eliminating granules not satisfying crossgran +/- n_crossgran
    """
    for group in contiguous_granule_groups(geo_granules):
        for geo_gran in group[n_crossgran:-n_crossgran]:
            yield geo_gran


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
        LOG.debug('before trimming, %d granule candidates for %s' % (len(meta[G.geo_cn]), G.geo_cn))
        # list of available geo products for this group
        sdr2edr = dict(zip(G.sdr_cns, G.edr_cns))  # FUTURE: This could be prebuilt and merged into guidebook
        geo_granules = _trim_geo_granules(meta[G.geo_cn])
        LOG.debug('found %d granule candidates' % len(geo_granules))

        # filter cross-granule, remove this if we can eliminate cross-granule input dependencies or make them optional
        geo_granules = list(_crossgran_filter(geo_granules))
        LOG.debug('after cross-granule dependencies, down to %d granules' % len(geo_granules))

        # check if we have at least one SDR collection that should produce output for each granule
        # if so, yield it
        LOG.debug('granules to check for SDR data: %s' % repr([x['N_Granule_ID'] for x in geo_granules]))
        for geo_granule in geo_granules:  # for each granule we have valid geo for
            geo_gran_id = geo_granule['N_Granule_ID']
            geo_gran_ver = geo_granule['N_Granule_Version']
            LOG.debug('checking available SDR collections for %s-v%s' % (geo_gran_id, geo_gran_ver))

            sdr_collections = []
            edr_collections = []
            # FUTURE: this will not properly identify cross-granule SDR problems, e.g. GEO is all present but I4 is not
            for sdr_cn in G.sdr_cns:  # see which SDR collections are available
                for g in meta[sdr_cn]:
                    if (g['N_Granule_ID'] == geo_gran_id) and (g['N_Granule_Version'] == geo_gran_ver):
                        sdr_collections.append(g['N_Collection_Short_Name'])
                        edr_collections.append(sdr2edr[g['N_Collection_Short_Name']])

            if not sdr_collections:
                LOG.warning('no SDR products found for %s:%s-v%s' % (kind, geo_gran_id, geo_gran_ver))
            else:
                LOG.debug(
                    'found SDR collections %s for %s:%s-v%s' % (repr(sdr_collections), kind, geo_gran_id, geo_gran_ver))
                yield (kind, geo_granule, sdr_collections, edr_collections, G.anc)


def anc_files_for_collections(collection_names, search_path_list, granules, optional=False):
    """
    find all the ancillary files that we need as a sequence to feed to link_ancillary_to_work_dir

    FUTURE: promote this into adl_asc and replace anc_files_needed with it, this should be more robust

    :param collection_names: sequence of collection names to search through metadata dictionaries for
    :param search_path_list: sequences of directories to look in
    :param granules: granules we want to ensure effectivity for these collections (not implemented!)
    :param optional: whether or not to raise an error if nothing is found for these collections
    :return:
    """
    for cn in collection_names:
        nfound = 0
        for search_dir in search_path_list:
            if not search_dir:
                continue
            LOG.debug("searching %s for %r" % (search_dir, cn))
            for nfo in skim_dir(search_dir, required_keys=['N_Collection_Short_Name'], N_Collection_Short_Name=cn):
                asc_path = nfo[K_FILENAME]
                blob_path = os.path.join(os.path.split(asc_path)[0], os.path.split(nfo['BlobPath'])[-1])
                nfound += 1
                LOG.debug('for ancillary %s, found %r / %r' % (cn, asc_path, blob_path))
                if os.path.isfile(asc_path) and os.path.isfile(blob_path):
                    yield asc_path
                    yield blob_path
                else:
                    if not os.path.isfile(asc_path):
                        LOG.error('could not find ancillary metadata %s' % asc_path)
                    if not os.path.isfile(blob_path):
                        LOG.error('could not find ancillary binary large object %s' % blob_path)
        if nfound == 0 and not optional:
            err = 'no ancillary data found for pattern %r in directories %s' % (cn, repr(search_path_list))
            LOG.error(err)
            raise EnvironmentError(err)
        else:
            LOG.debug('found %d files for %r' % (nfound, cn))


def populate_static_ancillary_links(anc_dir, ancillary_cns, geo_granules):
    """
    search static ancillary for LUTs and other required collections
    also, where possibly verify time range of ancillary to ensure we're covered
    FUTURE: promote nearly common routine shared with adl_atms_sdr.py up into adl_common.py

    :param anc_dir: ancillary directory to link into
    :param ancillary_cns: sequence of collection names we're looking for
    :param geo_granules: granule metadata dictionary list (for things like time ranges)
    """
    if CSPP_RT_ANC_CACHE_DIR is not None:
        search_dirs = [CSPP_RT_ANC_CACHE_DIR]
    else:
        search_dirs = []
    search_dirs += list(CSPP_RT_ANC_PATH.split(':'))
    LOG.debug('searching %s for static ancillary %s' % (repr(search_dirs), repr(ancillary_cns)))
    # convert collection names to filename globs
    link_ancillary_to_work_dir(anc_dir, anc_files_for_collections(ancillary_cns, search_dirs, geo_granules))


# def populate_dynamic_ancillary_links(anc_dir, work_dir, granules_to_process, allow_cache_update=True):
#     """
#     search for dynamic ancillary data (polar wander and TLE) for the granules we intend to process
#     if it's not in the cache already, download it
#     softlink those into the workspace
#     :param anc_dir: ancillary directory to link files into
#     :param work_dir: where to download ancillary files in the case that cache is unavailable
#     :param granules_to_process: list of granule metadata dictionaries, provide time ranges etc
#     :param allow_cache_update: whether or not to allow the helper scripts to download from Internet
#     FUTURE: originally in adl_atms_sdr.py, consider promoting common routine to adl_common
#     """
#     LOG.info("downloading TLE and PolarWander ancillary into cache and linking into workspace")
#     polar_files = adl_anc_retrieval.service_remote_ancillary(work_dir, granules_to_process,
#                                                              adl_anc_retrieval.kPOLAR,
#                                                              allow_cache_update=allow_cache_update)
#     tle_files = adl_anc_retrieval.service_remote_ancillary(work_dir, granules_to_process,
#                                                            adl_anc_retrieval.kTLE,
#                                                            allow_cache_update=allow_cache_update)
#     all_dyn_anc = list(polar_files) + list(tle_files)
#     LOG.debug('dynamic ancillary files: %s' % repr(all_dyn_anc))
#     LOG.info("Link the required ancillary data into the workspace")
#     link_ancillary_to_work_dir(anc_dir, all_dyn_anc)


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
    fnxml = ('edr_viirs_gtm_%s_%s.xml' % (kind, name))
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
    status_line("Checking " + str(len(files)) + " log files for errors" + logDir + " exp " + logExpression)

    n_err = 0
    err_files = set()
    for log_file in files:
        LOG.info("Checking Log file " + log_file + " for errors.")
        count = adl_log.scan_log_file(COMMON_LOG_CHECK_TABLE + GTM_EDR_LOG_CHECK_TABLE, log_file)
        n_err += count
        if count > 0:
            err_files.add(log_file)

    if n_err == 0:
        status_line("Processing of file: " + xml + " Completed successfully")
        return True
    else:
        status_line("Processing of file: " + xml + " Completed unsuccessfully, Look at previous message")
        LOG.debug("Log files with errors: " + ', '.join(err_files))
        return False


def transfer_gtm_edr_output(work_dir, work_subdir, kind, gran, sdr_cns, edr_cns):
    """
    examine work_subdir for products, based on sdr_collections that were available as input
    transfer products back from work_subdir to work_dir
    return product_filenames sequence, and transfer error sequence
    """
    products = []
    errors = []
    # FUTURE: this is a reasonable wag at it;
    # it should use EDR CNs and N_GEO_Ref preferably,
    # or h5 prefix in guidebook
    for h5path in glob.glob(os.path.join(work_subdir, '*.h5')):
        LOG.debug('transferring output %s' % h5path)
        h5filename = os.path.split(h5path)[-1]
        h5out = os.path.join(work_dir, h5filename)
        try:
            shutil.move(h5path, h5out)
            products.append(h5filename)
        except IOError:
            errors.append('%s would not transfer' % h5path)
    return products, errors


task_input = namedtuple('task_input', 'kind granule sdr_cns edr_cns work_dir env_dict cleanup aggregate compress')
task_output = namedtuple('task_output', 'kind granule_id product_filenames error_list')


def task_gtm_edr(task_in):
    """
    process a single task, returning a task_output tuple
    expect up to ngranules * (MXX,IXX,NCC) tasks
    this is suitable for spinning off to a subprocess using multiprocessing.Pool
    named tuples of primitive types serialize well across process boundaries
    """
    kind, gran, sdr_cns, edr_cns, work_dir, additional_env = task_in[:6]
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
        os.makedirs(os.path.join(work_subdir, 'log'))

    # generate XML into work subdirectory
    xml_filename = generate_gtm_edr_xml(kind, gran, work_subdir)

    # run XML controller
    exe = G.exe
    cmd = [exe, xml_filename]
    local_env = {'WORK_SUBDIR': work_subdir, 'WORK_DIR': work_dir}
    local_env.update(additional_env)

    status_line('Executing %s' % repr(cmd))
    LOG.debug('additional environment variables: %s' % repr(local_env))
    try:
        pid = sh(cmd, env=env(**local_env), cwd=work_subdir)
        LOG.debug("%r ran as pid %d" % (cmd, pid))
        ran_ok = check_logs_for_run(work_dir, pid, xml_filename)
        if not ran_ok:
            errors.append('logs were not error-free')

    except CalledProcessError as oops:
        pid = getattr(oops, 'pid', None)
        errors.append('process indicated failure or crash')
        ran_ok = check_logs_for_run(work_dir, pid, xml_filename)
        if not ran_ok:
            errors.append('log file problem')
        LOG.debug(traceback.format_exc())
        LOG.error('%s failed on %r: %r. Continuing...' % (exe, xml_filename, oops))

    # set the N_GEO_Ref attributes
    # FUTURE: something more elegant than a glob would be nice, but this will do
    LOG.debug('writing N_GEO_Ref attributes')
    for fn in glob.glob(os.path.join(work_subdir, 'V*h5')):
        LOG.debug('setting N_GEO_Ref for %s' % fn)
        write_geo_ref(fn)
        # SIDE EFFECT: we've already patched adl_geo_ref for our products.
        # now we can be confident that it'll work for the aggregate products as well.

    if task_in.compress and not task_in.aggregate:
        # then we can compress in-task, which is potentially parallel
        try:
            LOG.debug('compressing unaggregated output files')
            repack_products(work_subdir, EDR_PRODUCTS_AGG)
            repack_products(work_subdir, EDR_GEO_PRODUCTS_AGG)
        except KeyError as impossible_case:
            from adl_post_process import SHORTNAME_2_PRODUCTID
            LOG.error(repr(list(SHORTNAME_2_PRODUCTID.keys())))
            raise

    # link the output from the work_subdir to the work_dir
    product_filenames, transfer_errors = transfer_gtm_edr_output(work_dir, work_subdir, kind, gran, sdr_cns, edr_cns)
    errors += list(transfer_errors)

    # if everything ran OK, clean up the intermediate stuff in our subdir
    if not errors and task_in.cleanup:
        LOG.debug('cleaning up %s, no errors' % work_subdir)
        shutil.rmtree(work_subdir)

    return task_output(kind, granule_id, product_filenames, errors)


def herd_viirs_gtm_edr_tasks(work_dir, anc_dir, input_dir, nprocs=1, allow_cache_update=True,
                             cleanup=True, aggregate=False, compress=False, **additional_env):
    """
    skim work directory ASC metadata and decide which granules to process, and which controller to use
    generate task objects for processing candidates
    populate static and dynamic ancillary data with softlinks and/or data downloads
    execute tasks linearly (nprocs=1) or in parallel using multiprocessing Pool
    return result tuples indicating output files and any errors encountered

    :param input_dir: where to obtain inputs in ASC+BLOB form
    :param work_dir: outer work area in which we create GTM_* sub-workdirs
    :param anc_dir: linked_data ancillary directory, populated here and shared between tasks
    :param nprocs: number of processes to run in parallel for tasking
    :param allow_cache_update: whether or not to permit ancillary helper scripts to download from Internet
    :param cleanup: whether to clean up work directory afterward or not
    :param aggregate: whether to aggregate the output data into single HDF5 file
    :param compress: whether to run h5repack on output files or not
    :param additional_env: environment variables to pass on to controllers run by tasks
    :return: array of result tuples
    """
    tasks = []
    LOG.debug('sifting unpacked metadata for candidate granules to process')
    all_info = list(sift_metadata_for_viirs_gtm_edr(input_dir))

    all_anc_cns = set()
    all_geo_grans = []
    for kind, geo_granule, sdr_cns, edr_cns, anc_cns in all_info:
        all_anc_cns.update(set(anc_cns))
        all_geo_grans.append(geo_granule)

    # track down the set of all ancillary we'll be needing for these granules
    # FUTURE: consider doing this at a finer granularity
    LOG.debug('expect to need these static ancillary collections: %s' % repr(all_anc_cns))
    populate_static_ancillary_links(anc_dir, all_anc_cns, all_geo_grans)
    # LOG.debug('fetching dynamic ancillary for %d granules' % len(all_geo_grans))
    # populate_dynamic_ancillary_links(anc_dir, work_dir, all_geo_grans, allow_cache_update)

    LOG.debug('building task list for parallel processing')
    for kind, geo_granule, sdr_cns, edr_cns, anc_cns in all_info:
        tasks.append(task_input(kind, geo_granule, sdr_cns, edr_cns, work_dir, additional_env,
                                cleanup, aggregate, compress))

    if not tasks:
        LOG.error('no VIIRS SDR data was found to process - was VIIRS SDR run with --edr in order to preserve IP blobs?')

    LOG.debug('task list:')
    LOG.debug(pformat(tasks))

    if nprocs == 1:
        LOG.debug('running %d tasks without parallelism' % len(tasks))
        results = map(task_gtm_edr, tasks)
    else:
        LOG.debug('creating process pool size %d for %d tasks' % (int(nprocs), len(tasks)))
        parallel = Pool(int(nprocs))
        try:
            results = parallel.map(task_gtm_edr, tasks)
        except (KeyboardInterrupt, SystemError):
            # note that we're depending on register_sigterm having been called for SystemError on SIGTERM
            LOG.warning('external termination detected, aborting subprocesses')
            parallel.terminate()
            parallel.join()
            raise

    if aggregate:
        LOG.debug('aggregating collected products')
        aggregate_products(work_dir, EDR_GEO_PRODUCTS_AGG)
        aggregate_products(work_dir, EDR_PRODUCTS_AGG)
        _patch_geo_guide()
        add_geo_attribute_to_aggregates(work_dir, EDR_PRODUCTS_AGG)
        if compress:
            LOG.debug('compressing aggregated products')
            repack_products(work_dir, EDR_PRODUCTS_AGG)
            repack_products(work_dir, EDR_GEO_PRODUCTS_AGG)

    return results


def setup_directories(work_dir, anc_dir):
    """Create the working directory and a subdirectory for the logs
    :param work_dir: directory which we'll be creating work files in
    :param anc_dir: ancillary directory we'll be linking in
    """
    if not os.path.isdir(work_dir):
        LOG.info('Creating directory %s' % work_dir)
        os.makedirs(work_dir)

    log_dir = os.path.join(work_dir, 'log')   # FUTURE D-R-Y: this is also done in setup_directories
    if not os.path.isdir(log_dir):
        LOG.info('Creating log directory %s' % log_dir)
        os.makedirs(log_dir)

    if not os.path.exists(anc_dir):
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


def viirs_gtm_edr(work_dir, input_dir, nprocs=1, allow_cache_update=True,
                  cleanup=True, aggregate=False, compress=False):
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
    :param input_dir: directory containing intermediate and end-products from VIIRS SDR
    :param nprocs: number of processors to use, default 1
    :param compress: whether or not to apply h5repack to output HDF5 files
    :param aggregate: whether or not to use nagg to aggregate output HDF5 files
    :param cleanup: whether or not to clean up the work directory after finishing, if no errors occurred
    :param allow_cache_update: whether or not to allow dynamic ancillary to be downloaded
    """


    anc_dir = os.path.join(work_dir, ANCILLARY_SUB_DIR)
    setup_directories(work_dir, anc_dir)

    if not os.path.isdir(input_dir):
        LOG.error('%s is not a valid directory; please run VIIRS SDR with --edr' % input_dir)
        return 1
    input_dir = os.path.abspath(input_dir)

    # FUTURE: migrate directory environment settings to inner routine for better consistency
    results = herd_viirs_gtm_edr_tasks(work_dir, anc_dir, input_dir,
                                       nprocs=nprocs,
                                       allow_cache_update=allow_cache_update,
                                       cleanup=cleanup,
                                       aggregate=aggregate,
                                       compress=compress,
                                       LINKED_ANCILLARY=anc_dir,
                                       INPUT_DIR=input_dir,
                                       ADL_HOME=ADL_HOME,
                                       CSPP_RT_ANC_TILE_PATH=CSPP_RT_ANC_TILE_PATH,
                                       DCONFIG=ADL_VARS['DCONFIG'])
    LOG.debug(repr(results))
    error_count = 0
    LOG.info('Product list')

    if not results:
        LOG.error('no GTM products were created')
        error_count = 1

    for kind, granule, products, errors in results:
        LOG.info('%s:%s => %r' % (kind, granule, products))
        if errors:
            LOG.warning(' Problems reported: ' + ', '.join(errors))
            error_count += len(errors)
        else:
            LOG.info('  No problems reported')

    if error_count:
        status_line('%d problems were reported - please review above listing.')
    else:
        status_line('Completed successfully!')

    if cleanup and not error_count:
        LOG.debug('cleaning out ancillary and log directories')
        shutil.rmtree(anc_dir, ignore_errors=True)
        shutil.rmtree(os.path.join(work_dir, 'log'), ignore_errors=True) # FUTURE D-R-Y
    return error_count


def main():
    """ Run Viirs GTM EDR processing
    """
    import argparse

    desc = """Build VIIRS GTM EDR work directory and run VIIRS GTM EDR."""
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-t', '--test',
                        action="store_true", default=False, help="run self-tests")
    parser.add_argument('-W', '--work-dir', metavar='work_dir', default='.',
                        help='work directory which all activity will occur in, defaults to current dir')

    parser.add_argument('-d', '--debug',
                        action="store_true", default=False, help="always retain intermediate files")

    parser.add_argument('-z', '--zip',
                        action="store_true", default=False, help="compress products with h5repack zip compression")

    # parser.add_argument('-a', '--aggregate',
    #                     action="store_true", default=False, help="aggregate products with nagg")

    parser.add_argument('-l', '--local',
                        action="store_true", default=False, help="disable download of remote ancillary data to cache")

    parser.add_argument('-p', '--processors',
                        type=int, default=1, help="Number of processors to use for band processing")

    parser.add_argument('-v', '--verbosity', action="count", default=0,
                        help='each occurrence increases verbosity 1 level through ERROR-WARNING-INFO-DEBUG')

    parser.add_argument('input_dir', metavar='filename', type=str, nargs='+',
                        help='VIIRS SDR workspace directory to source SDR and intermediate products from')

    args = parser.parse_args()

    levels = [logging.ERROR, logging.WARN, logging.INFO, logging.DEBUG]
    level = levels[args.verbosity if args.verbosity < 4 else 3]

    do_cleanup = True

    work_dir = check_and_convert_path("WORK_DIR", args.work_dir)
    d = dt.date.today()
    timestamp = d.isoformat()
    log_name = "viirs_gtm_edr.%s.log" % timestamp
    logfile = os.path.join(work_dir, log_name)
    configure_logging(level, FILE=logfile)
    if args.debug:
        do_cleanup = False

    LOG.debug("Clean up: " + str(do_cleanup))
    LOG.info('CSPP execution work directory is %r' % work_dir)

    if 'TEST_GEO_REF' in os.environ:
        for fn in args.input_dir:
            print '%s => %s' % (fn, geo_ref(fn))
            return 0

    if 'TEST_GUIDE' in os.environ:
        from pprint import pprint
        pprint(list(_geo_guide()))
        return 0


    # if args.test:
    #     check_env(work_dir)
    #     grans = _test_sdr_granules(work_dir)
    #     if grans:
    #         LOG.debug('building XML files')
    #
    #     sys.exit(2)

    if not args or len(args.input_dir)>1:
        parser.print_help()
        return 1

    # register_sigterm()

    num_procs = args.processors
    if num_procs <= 0:
        num_procs = cpu_count()
        LOG.info('using %d cores' % num_procs)

    rc = viirs_gtm_edr(work_dir, args.input_dir[0], nprocs=num_procs, allow_cache_update=not args.local,
                       compress=args.zip, aggregate=False, cleanup=do_cleanup)

    if rc == 0 and not args.debug:
        os.remove(logfile)

    return rc


if __name__ == '__main__':
    sys.exit(main())
