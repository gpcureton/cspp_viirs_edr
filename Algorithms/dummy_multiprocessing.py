#!/usr/bin/env python
# encoding: utf-8
"""
dummy_multiprocessing.py

 * DESCRIPTION:  Class containing data relevent to the Dummy Multiprocessing Algorithm. 

Created by Geoff Cureton on 2013-08-22.
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
from random import randint, normalvariate
from shutil import rmtree, move
import numpy as np
import tables as pytables
import multiprocessing

from Utils import check_log_files, _setupAuxillaryFiles, getURID

# skim and convert routines for reading .asc metadata fields of interest
#import adl_asc
#from adl_asc import skim_dir, contiguous_granule_groups, granule_groups_contain, effective_anc_contains,_eliminate_duplicates,_is_contiguous, RDR_REQUIRED_KEYS, POLARWANDER_REQUIRED_KEYS
from adl_asc import skim_dir, granule_groups_contain, _eliminate_duplicates,_is_contiguous, skim, RE_LINE, convert, K_FILENAME
from adl_asc import PAT_URID, PAT_GRANULE_ID, PAT_GRANULE_VERSION, PAT_COLLECTION, PAT_RANGEDATETIME, PAT_SOURCE, PAT_BLOBPATH, PAT_EFFECTIVEDATETIME, PAT_OBSERVEDDATETIME
from adl_common import sh, unpack, env, h5_xdr_inventory
from adl_common import ADL_HOME, CSPP_RT_HOME, CSPP_RT_ANC_PATH, CSPP_RT_ANC_HOME, CSPP_RT_ANC_CACHE_DIR, COMMON_LOG_CHECK_TABLE

# N_GEO_Ref fix for ADL
from adl_geo_ref import write_geo_ref

# log file scanning
import adl_log

# every module should have a LOG object
try :
    sourcename= file_Id.split(" ")
    LOG = logging.getLogger(sourcename[1])
except :
    LOG = logging.getLogger('dummy_multiprocessing')

AlgorithmString = 'MPC'

AlgorithmName = 'Multiprocessing Test'

GEO_collectionShortNames = [
                            'VIIRS-MOD-GEO-TC'
                          ]

SDR_collectionShortNames = [
                              'VIIRS-M1-SDR'
                          ]


ANC_collectionShortNames = [
                           'VIIRS-ANC-Preci-Wtr-Mod-Gran'
                          ]

GridIP_collectionShortNames = [
                            'VIIRS-GridIP-VIIRS-Nbar-Ndvi-Mod-Gran'
                          ]

AUX_collectionShortNames = [
                            'VIIRS-CM-IP-AC-Int'
                           ]

AUX_ascTemplateFile = [
                        'VIIRS-CM-IP-AC-Int_Template.asc'
                      ]

AUX_blobTemplateFile = [
                         'template.VIIRS-CM-IP-AC-Int'
                       ]

AUX_Paths = [
             'luts/viirs'
            ]

controllerBinary = 'ProEdrViirsDummyController.exe'
ADL_VIIRS_DUMMY_EDR=path.abspath(path.join(ADL_HOME, 'bin', controllerBinary))

algorithmLWxml = 'edr_dummy_mpc'

# Attribute paths for MultiProcessing Dummy Product
attributePaths = {}
attributePaths['MOD_GEO_TC'] = 'Data_Products/VIIRS-MOD-GEO-TC/VIIRS-MOD-GEO-TC_Gran_0/N_Granule_ID'
attributePaths['DUMMY'] = 'Data_Products/VIIRS-DUMMY/VIIRS-DUMMY_Gran_0/N_Granule_ID'

MOD_GEO_TC_GRANULE_ID_ATTR_PATH = 'Data_Products/VIIRS-MOD-GEO-TC/VIIRS-MOD-GEO-TC_Gran_0/N_Granule_ID'
DUMMY_GRANULE_ID_ATTR_PATH = 'Data_Products/VIIRS-DUMMY/VIIRS-DUMMY_Gran_0/N_Granule_ID'

# maximum delay between granule end time and next granule start time to consider them contiguous
MAX_CONTIGUOUS_DELTA=timedelta(seconds = 2)

# XML template for ProEdrViirsMasksController.exe
# from ADL/cfg/dynamic/withMetadata/ProEdrViirsMasksControllerLwFile.xml
xmlTemplate = """<InfTkConfig>
  <idpProcessName>ProEdrViirsDummyController.exe</idpProcessName>
  <siSoftwareId />
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
  <outputPath>${GRANULE_OUTPUT_DIR}</outputPath>
  <dataStartIET>0000000000000000</dataStartIET>
  <dataEndIET>1111111111111111</dataEndIET>
  <actualScans>47</actualScans>
  <previousActualScans>48</previousActualScans>
  <nextActualScans>48</nextActualScans> 
  <usingMetadata>TRUE</usingMetadata>
  <configGuideName>ProEdrViirsDummyController_GuideList.cfg</configGuideName>

  <task>
    <taskType>EDR</taskType>
    <taskDetails1>%(N_Granule_ID)s</taskDetails1>
    <taskDetails2>%(N_Granule_Version)s</taskDetails2>
    <taskDetails3>NPP</taskDetails3>
    <taskDetails4>VIIRS</taskDetails4>
  </task>
</InfTkConfig>
"""


'''
def _contiguous_granule_groups(granules, tolerance=MAX_CONTIGUOUS_DELTA, larger_granules_preferred=False):
    """
    given a sequence of granule dictionaries, yield a sequence of contiguous granule groups as tuples
    tolerance, if provided, is a datetime.timedelta object representing max differance between endtime and starttime

    This is a custom version of adl_asc.contiguous_granule_groups(), which keys off of 'StartTime', rather than
    'ObservedStartTime' as is done here.
    """
    
    # sort granules into start time order and eliminate exact duplicates
    # FUTURE: is lex-compare sufficient for A2/A1/etc
    start_time_key = lambda x: (x['StartTime'], x.get('N_Granule_Version', None))
    #start_time_key = lambda x: (x['RangeDateTime'], x.get('N_Granule_Version', None))
    #sys.exit(0)
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
        LOG.debug('contiguous sequence has %d granules' % (len(seq)))
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
                LOG.debug('contiguous sequence has %d granules' % (len(seq)))
                yield tuple(sorted(seq.values(), key=start_time_key))
                seq.clear()
        # leftovers! yum!
        if seq:
            LOG.debug('Leftover contiguous sequence has %d granules' % (len(seq)))
            yield tuple(sorted(seq.values(), key=start_time_key))
'''

'''
def sift_metadata_for_viirs_ancillary(collectionShortName, crossGran=None, work_dir='.'):
    """
    Search through the ASC metadata in a directory, grouping in StartTime order.
    Look for back-to-back granules and break into contiguous sequences.
    Yield sequence of granules we can process.
    """
    LOG.debug('Collecting information for VIIRS SDRs')

    work_dir = path.abspath(work_dir)

    ancGroupList = list(_contiguous_granule_groups(skim_dir(work_dir, N_Collection_Short_Name=collectionShortName)))

    LOG.debug('ancGroupList : %r'%(ancGroupList))

    if len(ancGroupList)==0:
        LOG.debug('No ancGroupList found...')
        return

    # Loop through the contigous granule groups 
    for group in _contiguous_granule_groups(skim_dir(work_dir, N_Collection_Short_Name=collectionShortName)):
        ##- for VIIRS, we can process everything but the first and last granule
        ##- for CrIS, use [4:-4]
        LOG.debug('Contiguous granule group of length: %r' % (len(group),))

        if not crossGran :
            startGran,endGran = None,None
        else :
            startGran,endGran = crossGran,-1*crossGran

        for gran in group[startGran:endGran]:
            if not granule_groups_contain(ancGroupList, gran):
                LOG.warning("Insufficient VIIRS SDR coverage to process %s @ %s (%s) - skipping" % (gran['N_Granule_ID'], gran['StartTime'], gran['URID']))
                continue
                #pass
            LOG.debug('Processing opportunity: %r at %s with uuid %s' % (gran['N_Granule_ID'], gran['StartTime'], gran['URID']))
            yield gran
'''


def setupAuxillaryFiles(Alg_objects,workDir):
    '''
    Call the generic Utils method to link in the auxillary files 
    specified in Alg_objects to the workDir directory.
    '''

    _setupAuxillaryFiles(Alg_objects,workDir)


def generate_viirs_edr_xml(work_dir, granule_seq):
    "generate XML files for VIIRS dummy granule generation"
    to_process = []
    for gran in granule_seq:
        granId = gran['N_Granule_ID']
        xmlName = '%s_%s.xml' % (algorithmLWxml,granId)
        LOG.debug('writing XML file %r' % (xmlName))
        xml_fileObj = file(path.join(work_dir, xmlName), 'wt')
        xml_fileObj.write(xmlTemplate % gran)
        to_process.append([granId,xmlName])
    return to_process


def dummy_exe(additional_env):
    '''
    Dummy "executable"
    '''

    work_dir = additional_env['WORK_DIR']
    xml = additional_env['XML_FILE']
    granule_id = additional_env['N_Granule_ID']
    granule_output_dir = additional_env['GRANULE_OUTPUT_DIR']

    if granule_id in ['NPP000603201736', 'NPP000603197469']:
        # Simulate a failed granule
        pass
    else:
        # sleep() and dummy HDF5 file creation in lieu of an ADL algorithm...
        #sleep(normalvariate(2.5,0.01))
        sleep(2.)
        
        # Make the blob/asc output
        URID = getURID()['URID']
        ascFileName = '{}.asc'.format(path.join(granule_output_dir,URID))
        open(ascFileName, 'a').close()
        blobFileName = '{}.VIIRS-DUMMY'.format(path.join(granule_output_dir,URID))
        open(blobFileName, 'a').close()

        # Make the HDF5 output
        h5File = '%s/MPC_%s.h5'%(granule_output_dir,granule_id)
        LOG.info("Creating file {}...".format(h5File))
        fileObj = pytables.openFile(h5File, mode='w')
        fileObj.createGroup('/Data_Products/VIIRS-DUMMY','VIIRS-DUMMY_Gran_0',createparents=True)
        granID_array = np.array([[granule_id]])
        fileObj.setNodeAttr('/Data_Products/VIIRS-DUMMY/VIIRS-DUMMY_Gran_0','N_Granule_ID',granID_array)
        fileObj.close()

    # Create a log file corresponding to this granule
    logTime = datetime.utcnow()
    pid = randint(500,30000)
    logName = 'ProEdrViirsDummyController.exe_%s_%d' % (logTime.strftime("%Y%m%d_%H%M%S"),pid)
    os.mkdir(path.join(granule_output_dir,logName))
    open('%s.log'%(path.join(granule_output_dir,'log',logName)), 'a').close()

    return pid


def move_products_to_work_directory(work_dir):
    """ Checks that proper products were produced.
    If product h5 file was produced blob and asc files are deleted
    If product was not produced files are left alone
    """
    dest = path.dirname(work_dir)

    files = glob(path.join(work_dir, "*.h5"))
    for file in files:
        move(file, dest)

    return


def submit_granule(additional_env):
    "run a VIIRS EDR XML input in sequence"

    work_dir = additional_env['WORK_DIR']
    xml = additional_env['XML_FILE']
    granule_id = additional_env['N_Granule_ID']
    granule_output_dir = additional_env['GRANULE_OUTPUT_DIR']

    granule_diagnostic = {}
    granule_diagnostic['crashed'] = False
    granule_diagnostic['no_output'] = False
    granule_diagnostic['geo_problem'] = False
    granule_diagnostic['bad_log'] = False
    granule_diagnostic['N_Granule_ID'] = granule_id

    # Pattern for expected output
    dummyPattern = path.join(granule_output_dir, 'MPC*.h5')

    # prior_granules dicts contain (N_GranuleID,HDF5File) key,value pairs.
    # *ID dicts contain (HDF5File,N_GranuleID) key,value pairs.
    
    # Get the (N_GranuleID,hdfFileName) pairs for the existing dummy files
    dummy_prior_granules, dummy_ID = h5_xdr_inventory(dummyPattern, DUMMY_GRANULE_ID_ATTR_PATH)
    dummy_prior_granules = set(dummy_prior_granules.keys())
    LOG.debug('Set of existing dummy granules... {}'.format(repr(dummy_prior_granules)))

    # Specify the command line to execute.
    cmd = [ADL_VIIRS_DUMMY_EDR, xml]
    #cmd = ['/usr/bin/gdb', ADL_VIIRS_DUMMY_EDR] # for debugging with gdb...

    LOG.debug('executing "{}"'.format(' '.join(cmd)))
    LOG.debug('additional environment variables: {}'.format(additional_env))

    t1 = time()

    try:
        
        #pid = sh(cmd, env=env(**additional_env), cwd=work_dir)
        pid = dummy_exe(additional_env)

        LOG.debug("{} ran as pid {}".format(cmd, pid))
        if not check_log_files(granule_output_dir, pid, xml):
            granule_diagnostic['bad_log'] = True

    except CalledProcessError as oops:
        pid = getattr(oops, 'pid', None)
        LOG.debug(traceback.format_exc())
        LOG.error('{} failed on {}: {}. Continuing...' % (cmd[0], xml, oops))
        granule_diagnostic['crashed'] = True

    # check new MPC output granules
    dummy_new_granules, dummy_ID = h5_xdr_inventory(dummyPattern, DUMMY_GRANULE_ID_ATTR_PATH, state=dummy_ID)
    if granule_id not in dummy_new_granules:
        LOG.warning('no DUMMY HDF5 output for {}'.format(granule_id))
        granule_diagnostic['no_output'] = True
        granule_diagnostic['output_file'] = None
    else :
        LOG.info('New MPC granule: {}'.format(repr(dummy_new_granules)))
        dummy_granules_made = set(dummy_ID.values()) - dummy_prior_granules
        LOG.info('{} granules created: {}'.format(AlgorithmName,', '.join(list(dummy_granules_made))))
        granule_diagnostic['output_file'] = path.basename(dummy_new_granules[granule_id])

    t2 = time()

    LOG.info ("Controller ran in {} seconds.".format(t2-t1))

    move_products_to_work_directory(granule_output_dir)

    return granule_diagnostic


def run_xml_files(work_dir, xml_files_to_process, CLEANUP="True", nprocs="1", **additional_env):
    """Run each VIIRS dummy xml input in sequence.
       Return the list of granule IDs which crashed, 
       and list of granule IDs which did not create output.
    """

    total_granules = len(xml_files_to_process)
    LOG.info('{} granules to process'.format(total_granules))

    argument_dictionaries = []
    for granule_id, xml in xml_files_to_process:

        granule_output_dir = path.join(work_dir,"ProEdrViirsDummyController_%s" %(granule_id))

        try:
            os.mkdir(granule_output_dir)
            os.mkdir(os.path.join(granule_output_dir, "log"))

            additional_envs = dict(
                N_Granule_ID=granule_id,
                XML_FILE=xml,
                GRANULE_OUTPUT_DIR=granule_output_dir,
                WORK_DIR=work_dir,
                ADL_HOME=ADL_HOME,
                CLEANUP=CLEANUP
            )

            argument_dictionaries.append(additional_envs)

        except:
            LOG.info("Granule %s Skipped,  Already processed" % (granule_id))

    # Pattern for expected output in the root working directory
    dummyPattern = path.join(work_dir,'MPC*.h5')

    # prior_granules dicts contain (N_GranuleID,HDF5File) key,value pairs.
    # *ID dicts contain (HDF5File,N_GranuleID) key,value pairs.
    
    # Get the (N_GranuleID,hdfFileName) pairs for the existing dummy files
    dummy_prior_granules, dummy_ID = h5_xdr_inventory(dummyPattern, DUMMY_GRANULE_ID_ATTR_PATH)
    LOG.debug('Existing dummy granules... {}'.format(dummy_prior_granules))

    dummy_prior_granules = set(dummy_prior_granules.keys())
    LOG.debug('Set of existing dummy granules... {}'.format(dummy_prior_granules))

    results = []

    if len(argument_dictionaries) > 0:

        # Create the multiprocessing infrastructure
        number_available = multiprocessing.cpu_count()

        nprocs = 4 # FIXME: temporary
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

    # check new MPC output granules
    dummy_new_granules, dummy_ID = h5_xdr_inventory(dummyPattern, DUMMY_GRANULE_ID_ATTR_PATH, state=dummy_ID)

    LOG.debug("dummy_ID.values() = %r" % (dummy_ID.values()))
    LOG.debug("set(dummy_ID.values()) = %r" % (set(dummy_ID.values())))
    LOG.debug("dummy_prior_granules = %r" % (dummy_prior_granules))

    dummy_granules_made = set(dummy_ID.values()) - dummy_prior_granules
    LOG.info('{} granules created: {}'.format(AlgorithmName,', '.join(list(dummy_granules_made))))

    crashed_runs = set()
    no_output_runs = set()
    geo_problem_runs = set()
    bad_log_runs = set()

    for dicts in results:
        LOG.info("results[{}] : {}".format(dicts['N_Granule_ID'],dicts))
        if dicts['crashed']: crashed_runs.add(dicts['N_Granule_ID'])
        if dicts['no_output']: no_output_runs.add(dicts['N_Granule_ID']) 
        if dicts['geo_problem']: geo_problem_runs.add(dicts['N_Granule_ID']) 
        if dicts['bad_log']: bad_log_runs.add(dicts['N_Granule_ID']) 

    if no_output_runs:
        LOG.info('Granules that failed to generate output: %s' % (', '.join(no_output_runs)))
    if geo_problem_runs:
        LOG.warning('Granules which had no N_Geo_Ref: %s' % ', '.join(geo_problem_runs))
    if crashed_runs:
        LOG.warning('Granules that crashed ADL: %s' % (', '.join(crashed_runs)))
    if bad_log_runs:
        LOG.warning('Granules that produced logs indicating problems: %s' % (', '.join(bad_log_runs)))
    if not dummy_granules_made:
        LOG.warning('No {} HDF5 files were created'.format(AlgorithmName))

    return crashed_runs, no_output_runs, geo_problem_runs, bad_log_runs


def cleanup(work_dir, xml_glob, log_dir_glob, *more_dirs):
    """upon successful run, clean out work directory"""

    LOG.info("Cleaning up work directory...")

    LOG.info("Removing {} job directories...".format(AlgorithmName))
    job_dir_glob = path.join(work_dir,"ProEdrViirsDummyController_NPP*")
    job_dirs = glob(job_dir_glob)
    for dir in job_dirs:
        LOG.info('removing job directory {}'.format(dir))
        try :
            rmtree(dir, ignore_errors=False)
        except Exception, err:
            LOG.warn( "{}".format(str(err)))

    LOG.info("Removing task xml files...")
    for fn in glob(path.join(work_dir, xml_glob)):
        LOG.debug('removing task file {}'.format(fn))
        try :
            os.unlink(fn)
        except Exception, err:
            LOG.warn( "{}".format(str(err)))

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

