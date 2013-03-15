#!/usr/bin/env python
# encoding: utf-8
"""
$Id:$

Purpose: Run the VIIRS SDR using Raytheon ADL 3.1 or 4.0

Input:
    One or more HDF5 VIIRS RDR input files, aggregated or single-granule.
    You need one granule before and one granule after a given granule to process.
    A work directory, typically, empty, in which to unpack the granules and generate the output.
    If the work directory specified does not exist, it will be created.
    
Output:
    ADL VIIRS SDR blob files, .asc metadata files, and HDF5 output granules will be created.

Details:
    If you input a series of granules, the software will scan the work directory.
    Thus, for N input RDR granules, you may get up to N-2 output SDR granules if they are all contiguous.
    It is ambiguous to provide several copies of the same granule in the work directory;
    this will result in an error abort.
    The unpacker gives each unpacking of a given granule its own 
    
Preconditions:
    Requires ADL_HOME, CSPP_SDR_HOME, CSPP_RT_ANC_CACHE_DIR, CSPP_ANC_HOME  environment variables are set.
    Requires that any needed LD_LIBRARY_PATH is set.
    Requires that DSTATICDATA is set.
    

Copyright (c) 2011 University of Wisconsin Regents. 
Licensed under GNU GPLv3.
"""

import os, sys, logging, glob, traceback,  time
import datetime as dt
from subprocess import CalledProcessError, call

# skim and convert routines for reading .asc metadata fields of interest
from adl_asc import skim_dir, contiguous_granule_groups,  RDR_REQUIRED_KEYS

import adl_log, adl_geo_ref
import  adl_anc_retrieval 

import xml.etree.ElementTree as ET
from adl_common import status_line , configure_logging,get_return_code,check_env,check_and_convert_path
import shutil
# maximum delay between granule end time and next granule start time to consider them contiguous
MAX_CONTIGUOUS_SDR_DELTA=dt.timedelta(seconds = 600)

# ancillary search and unpacker common routines
from adl_common import sh, anc_files_needed, link_ancillary_to_work_dir,  unpack, env
from adl_common import  unpack_h5s
from adl_common import  COMMON_LOG_CHECK_TABLE, EXTERNAL_BINARY,CSPP_RT_ANC_CACHE_DIR,CSPP_RT_ANC_PATH,DDS_PRODUCT_FILE,ADL_HOME,CSPP_RT_ANC_TILE_PATH

from adl_post_process import repack_products,aggregate_products,add_geo_attribute_to_aggregates

ADL_VIIRS_SDR=os.path.join(ADL_HOME, 'bin', 'ProSdrViirsController.exe')



# every module should have a LOG object
LOG = logging.getLogger('adl_viirs_sdr')

# keys used in metadata dictionaries
K_FILENAME = '_asc_filename'

ANCILLARY_SUB_DIR="linked_data"

optimize=True

# locations of executables in ADL 3.1



# and the patterns we're looking for
ADL_VIIRS_ANC_GLOBS =  (
    '*VIIRS-SDR-DNB-F-PREDICTED-LUT*',   # adl 4.1
    '*VIIRS-SDR-F-PREDICTED-LUT*',       # adl 4.1
#  ProSdrCmnGeo   
    '*CMNGEO-PARAM-LUT_npp*',
    '*off_Planet-Eph-ANC*',
    '*off_USNO-PolarWander*', 
    '*CmnGeo-SAA-AC_npp*',
#    '*Terrain-Eco-ANC-Tile*',
# RDR processing

# GEO processing
    '*VIIRS-SDR-GEO-DNB-PARAM-LUT_npp*',
    '*VIIRS-SDR-GEO-IMG-PARAM-LUT_npp*',   
    '*VIIRS-SDR-GEO-MOD-PARAM-LUT_npp*',
## CAL Processing
    '*VIIRS-SDR-DNB-DN0-LUT_npp*',
    '*VIIRS-SDR-DNB-RVF-LUT_npp*',
    '*VIIRS-SDR-DG-ANOMALY-DN-LIMITS-LUT_npp*',
    '*VIIRS-SDR-DNB-STRAY-LIGHT-LUT_npp*',
    '*VIIRS-SDR-DNB-FRAME-TO-ZONE-LUT_npp*',
#    '*VIIRS-SDR-F-LUT_npp*',               # adl 4.0

    '*VIIRS-SDR-GAIN-LUT_npp*',
    '*VIIRS-SDR-HAM-ER-LUT*',
    '*VIIRS-SDR-RTA-ER-LUT*',
    '*VIIRS-SDR-OBC-ER-LUT_npp*',
    '*VIIRS-SDR-OBC-RR-LUT_npp*',
    '*VIIRS-SDR-EBBT-LUT_npp*',
    '*VIIRS-SDR-TELE-COEFFS-LUT_npp*',
    '*VIIRS-SDR-SOLAR-IRAD-LUT_npp*',
    '*VIIRS-SDR-RSR-LUT_npp*',
    '*VIIRS-SDR-OBS-TO-PIXELS-LUT_npp*',
    '*VIIRS-SOLAR-DIFF-VOLT-LUT_npp*',
    '*VIIRS-SDR-RADIOMETRIC-PARAM-LUT_npp*',
    '*VIIRS-SDR-QA-LUT_npp*',
    '*VIIRS-SDR-EMISSIVE-LUT_npp*',
    '*VIIRS-SDR-REFLECTIVE-LUT_npp*',
    '*VIIRS-SDR-RVF-LUT_npp*',
    '*VIIRS-SDR-BB-TEMP-COEFFS-LUT_npp*',
    '*VIIRS-SDR-DNB-C-COEFFS-LUT_npp*',
    '*VIIRS-SDR-DELTA-C-LUT_npp*',
    '*VIIRS-SDR-COEFF-A-LUT_npp*',
    '*VIIRS-SDR-COEFF-B-LUT_npp*',     
        '*TLE-AUX*'                
                    )


ADL_VIIRS_GEO_PRODUCT_SHORTNAMES = [
    'VIIRS-MOD-GEO-TC', 
    'VIIRS-IMG-GEO-TC',
    'VIIRS-DNB-GEO',
   ]


OPTIONAL_GEO_PRODUCTS=[
    'VIIRS-IMG-GEO',
    'VIIRS-MOD-GEO'
]

ADL_VIIRS_SDR_PRODUCT_SHORTNAMES = [

    'VIIRS-I1-SDR',
    'VIIRS-I2-SDR',
    'VIIRS-I3-SDR',
    'VIIRS-I4-SDR',
    'VIIRS-I5-SDR',
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
    'VIIRS-M13-SDR',
    'VIIRS-M14-SDR',
    'VIIRS-M15-SDR',
    'VIIRS-M16-SDR',
    'VIIRS-DNB-SDR'

]

ADL_VIIRS_SDR_intermediate_SHORTNAMES = [
    'VIIRS-MOD-UNAGG-GEO',  # no nagg
    'VIIRS-DualGain-Cal-IP',    # no nagg
    'VIIRS-OBC-IP', # no nagg

    
    'VIIRS-IMG-RGEO',
    'VIIRS-MOD-RGEO',
    'VIIRS-I1-FSDR',
    'VIIRS-I2-FSDR',
    'VIIRS-I3-FSDR',
    'VIIRS-I4-FSDR',
    'VIIRS-I5-FSDR',
    'VIIRS-M1-FSDR',
    'VIIRS-M2-FSDR',
    'VIIRS-M3-FSDR',
    'VIIRS-M4-FSDR',
    'VIIRS-M5-FSDR',
    'VIIRS-M6-FSDR',
    'VIIRS-M7-FSDR',
    'VIIRS-M8-FSDR',
    'VIIRS-M9-FSDR',
    'VIIRS-M10-FSDR',
    'VIIRS-M11-FSDR',
    'VIIRS-M12-FSDR',
    'VIIRS-M14-FSDR',
    'VIIRS-M15-FSDR',
    'VIIRS-M16-FSDR',
    'VIIRS-MOD-RGEO',
    'VIIRS-MOD-RGEO-TC'
]



CHECK_REQUIRED_KEYS = ['N_Granule_ID', 'N_Collection_Short_Name']

OBSERVE_TIME='ObservedStartTime'
# table of NPP short names to data product ids

PRODUCTID_2_SHORTNAME= dict()
SHORTNAME_2_PRODUCTID = dict()


  
        
# XML template for ProSdrViirsController.exe
XML_TMPL_VIIRS_SDR = """<InfTkConfig>
  <idpProcessName>ProSdrViirsController.exe</idpProcessName>
  <siSoftwareId></siSoftwareId>
  <isRestart>FALSE</isRestart>
  <useExtSiWriter>FALSE</useExtSiWriter>
  <debugLogLevel>NORMAL</debugLogLevel>
  <debugLevel>DBG_HIGH</debugLevel>
  <enablePerf>FALSE</enablePerf>
  <perfPath>${WORK_DIR}/log</perfPath>
  <dbgPath>${WORK_DIR}/log</dbgPath>
  <initData>
     <domain>OPS</domain>
     <subDomain>SUBDOMAIN</subDomain>
     <startMode>INF_STARTMODE_COLD</startMode>
     <executionMode>INF_EXEMODE_PRIMARY</executionMode>
     <healthTimeoutPeriod>30</healthTimeoutPeriod>
  </initData>
  <lockinMem>FALSE</lockinMem>
  <rootDir>${WORK_DIR}/log</rootDir>
  <inputPath>${WORK_DIR}:${LINKED_ANCILLARY}:${CSPP_RT_ANC_TILE_PATH}</inputPath>
  <outputPath>${WORK_DIR}</outputPath>
  <dataStartIET>0000000000000000</dataStartIET>
  <dataEndIET>1111111111111111</dataEndIET>
  <actualScans>48</actualScans>
  <previousActualScans>48</previousActualScans>
  <nextActualScans>48</nextActualScans> 
  <usingMetadata>TRUE</usingMetadata>
  <configGuideName>ProSdrViirs_GuideList.cfg</configGuideName>

  <task>
    <taskType>SDR</taskType>
    <taskDetails1>%(N_Granule_ID)s</taskDetails1>
    <taskDetails2>%(N_Granule_Version)s</taskDetails2>
    <taskDetails3>NPP</taskDetails3>
    <taskDetails4>VIIRS</taskDetails4>
  </task>

  <task>
    <taskType>SHUTDOWN</taskType>
    <taskDetails1></taskDetails1>
    <taskDetails2></taskDetails2>
    <taskDetails3></taskDetails3>
    <taskDetails4></taskDetails4>
  </task>
</InfTkConfig>
"""


def _skim_dir_collections(*args, **kwargs):
    """skim and return dictionary keyed on N_Collection_Short_Name, 
    each entry being a list of skimmed info dictionaries
    """
    from collections import defaultdict
    zult = defaultdict(list)
    for file_info in skim_dir(*args, **kwargs):
        key = file_info.get('N_Collection_Short_Name', '-unknown-collection-')
        zult[key].append(file_info)
    return zult


    
def sift_metadata_for_viirs_sdr(work_dir='.'):
    """
    search through the ASC metadata in a directory,
    grouping in StartTime order
    look for back-to-back granules and break into contiguous sequences
    check that S/C diary RDRs are available for granules of interest
    yield sequence of granules we can process
    """
    LOG.info('Collecting information for S/C diary RDRs')
    
    diaries = list(contiguous_granule_groups(skim_dir(work_dir,required_keys=RDR_REQUIRED_KEYS, N_Collection_Short_Name='SPACECRAFT-DIARY-RDR'),larger_granules_preferred=True))

    LOG.debug('sifting science RDR data for processing opportunities')
    Viirs_Science_RDRs = list(skim_dir(work_dir, N_Collection_Short_Name="VIIRS-SCIENCE-RDR"))
    
    status_line("Total Viirs Science RDRs: "+ str(len(Viirs_Science_RDRs)))


    start_time_key = lambda x: (x['StartTime'], x.get('N_Granule_Version', None))
 
    for group in sorted(Viirs_Science_RDRs, key = start_time_key):
            yield group           

def generate_sdr_xml(work_dir,gran):
        name = gran['N_Granule_ID']
        fnxml = 'sdr_viirs_%s.xml' % name
        LOG.debug('writing XML file %r' % fnxml)
        fpxml = file(os.path.join(work_dir, fnxml), 'wt')
        fpxml.write(XML_TMPL_VIIRS_SDR % gran)
        status_line('Created ADL controller XML %s for %s' % (fnxml,gran['N_Granule_ID']))

        return fnxml
#Description: Required input not available for Shortname

# table of ADL LOG error messages and the associated hint to correcting the problem
viirs_sdr_log_check_table = [
                            ("Verified RDR has invalid mode for geolocation/calibration","Likely bad data at start or end of pass "),\
                            ("PRO_FAIL Unable to read","Check CSPP/static and CSPP/cache directories"), \
                            ("Required input not available","Missing or out of date ancillary input, Check effectivity for listed Shortname."), \
                            ("PRO_FAIL runAlgorithm()","Algorithm failed"),("Completed unsuccessfully","Algorithm failed"), \
                            ("arbitrary time is invalid","Problem with input RDR,check NPP_GRANULE_ID_BASETIME"), \
                            ("Error retrieving data for USNO-POLARWANDER-UT1","POLAR WANDER file needs update,check NPP_GRANULE_ID_BASETIME"), \
                            ("Algorithm failed","Controller did not run, check log") \

    ] 
    
    


# Look through new log files for completed messages
def checkADLLogForSuccess(work_dir,pid,xml, remove_list) :
    """
    Find the log file 
    Look for success
    Return True if found
    Display log message and hint if problem occurred
    """
    start_time= time.time()
    LOG.debug("<< checkADLLogForSuccess  >>")
    
    # retrieve exe name and log path from lw file
    logDir = os.path.join(work_dir,"log")
    logExpression="*"+ str ( pid )+"*.lo*"
    
    
    files = glob.glob(os.path.join(logDir, logExpression))
    status_line("Checking "+str(len(files))+" log files for errors")
   
    n_err=0
    for logFile in files :
        count=0
        LOG.info( "Checking Log file " +logFile +" for errors.")
        
#        count = adl_log.scan_log_file(COMMON_LOG_CHECK_TABLE, logFile)            

        count = adl_log.scan_log_file(COMMON_LOG_CHECK_TABLE+viirs_sdr_log_check_table, logFile)            
        remove_list.append(logFile)
            
        n_err += count
       
    end_time=time.time() 
    elapsed_time = end_time-start_time
    LOG.debug("<< Check Log Elapsed time: %f >>"%(elapsed_time))
                    
         
        
        
    if n_err == 0 :
        status_line("Processing of file: "+ xml + " Completed successfully" )
        return True
    else :
        status_line("Processing of file: "+ xml + " Completed unsuccessfully, Look at previous message" )
        LOG.debug("Log file: "+logFile)
        return False


def is_granule_on_list(work_dir, wanted, found_granule_seq):
    """verify that granules matching an incoming sequence and specification were created
    returns sequence of granules that need to be converted
    example: blob_products_to_h5(work_dir, my_source_granules, adl_asc.skim_dir(work_dir, N_Collection_Short_Name='ATMS-SDR'))
    """
    found = dict((g['N_Granule_ID'],g) for g in found_granule_seq)

    name = wanted['N_Granule_ID']
             
    if name in found:
        LOG.debug('found granule for %s' % name)
        it=found.get(name)
        try :
            return it
        except KeyError: 
            LOG.error("No blob file for "+name)

    return None









def check_for_products(work_dir, gran_id, remove_list,product_dictionaries) :
    """ Checks that proper products were produced.
    If product h5 file was produced blob and asc files are deleted
    If product was not produced files are left alone
    """
    # look for all products that should have been produced
    # convert the products to H5
  
    global optimize

    
    LOG.info("Check for Granule: "+gran_id+" products.")
    start_time= time.time()
    LOG.debug("<< check_for_products  >>")



    problem=True
    total =0;
    good=0;
    for short_name in  sorted ( ADL_VIIRS_SDR_PRODUCT_SHORTNAMES) + sorted(  ADL_VIIRS_GEO_PRODUCT_SHORTNAMES )  :
        total = total + 1;
        
# must require granule id but I do not know it is used.  
        
        if optimize == True :
            gran_list = product_dictionaries[short_name]
        else :
            gran_list = list(skim_dir(work_dir,required_keys=CHECK_REQUIRED_KEYS,N_Collection_Short_Name=short_name,N_Granule_ID=gran_id))
 
        

        if len(gran_list) < 1 :
            try :
                LOG.debug("Problem : "+short_name+" No product produced.")
                SHORTNAME_2_PRODUCTID[short_name]
                problem = True
            except KeyError:
                LOG.warn("Short name not in product map:" + short_name)
        else :
            LOG.debug("Granule list length: "+ str(len( gran_list )))
            for it in  gran_list :
                try :
                    sdr_name=SHORTNAME_2_PRODUCTID[short_name]
                    dObj=it[OBSERVE_TIME ]
                
                    time_id_str  = dObj.strftime("_npp_d%Y%m%d_t%H%M%S")
             
                    #  Name Example: SVI04_npp_d20111121_t1805421_e1807062_b00346_c20120127203200212753_cspp_dev.h5
                    sdr_name=sdr_name+time_id_str+"*.h5"  
                    files = glob.glob(os.path.join(work_dir, sdr_name))
                    
                    
                    
                    if  len( files ) == 1 :
                        fullname= files[0]
                        LOG.debug("Product: "+fullname+" produced.")
                        
                        if fullname.find( "_d1958") != -1 :
                            LOG.warn("DELETE File with bad year: "+ fullname)
                            if os.path.exists(fullname)   :
                                os.remove(fullname)
                            problem = True
                        else :  
                            problem = False;
                            good = good + 1;
                            
                        " Night passes do not create a blob file for all asc files."
                            
                        ascname="noproperty"
                        blobname="noproperty"
                          
                        try :
                            blobname=it['BlobPath']
                        except KeyError:  
                            LOG.debug("Key error blob property")
                            LOG.debug("Blob: "+blobname )                         
                                    
                        try :
                            ascname=it['_filename']
                        except KeyError:  
                            LOG.debug("Key error on asc property")
                            LOG.debug("Blob: "+blobname )
                            b=blobname.split( ".")    
                            ascname=b[0]+".asc";
                                           
                        if  os.path.exists(ascname) :
                            remove_list.append( ascname )
                        if  os.path.exists(blobname) :
                            remove_list.append( blobname )
                    else :
                        LOG.error("H5 output: "+sdr_name+ " is missing")
                        a_exists=False;
                        b_exists=False;
                        if '_filename' in it.keys() :
                            ascname=it['_filename']
                            a_exists=os.path.exists(ascname)
                        if 'BlobPath' in it.keys() :		
                            blobname=it['BlobPath']
                            b_exists=os.path.exists(blobname)
                        LOG.error("Exists? "+str(a_exists)+" "+ascname)
                        LOG.error("Exists? "+str(b_exists)+" "+blobname)
                        
                except KeyError:  
                    LOG.error("Key error on blob property "+'ObservedStartTime' +" or "+short_name)
               
    LOG.info( str(good)+" out "+str(total)+" products produced for granule "+gran_id)
    
       
    end_time=time.time() 
    elapsed_time = end_time-start_time
    LOG.debug("<< Check Products Elapsed time: %f >>"%(elapsed_time))
                    

    return problem

##        skim_dir(work_dir,required_keys=CHECK_REQUIRED_KEYS,N_Collection_Short_Name=short_name,N_Granule_ID=gran_id)

def _skim_dir_collections(*args, **kwargs):
    """skim and return dictionary keyed on N_Collection_Short_Name, 
    each entry being a list of skimmed info dictionaries
    """
    from collections import defaultdict
    zult = defaultdict(list)
    for file_info in skim_dir(*args, **kwargs):
        key = file_info.get('N_Collection_Short_Name', '-unknown-collection-')
        zult[key].append(file_info)
    return zult

def run_xml_file(work_dir, xml_file, remove_list ,gran,setup_only=False, **additional_env):
    "run each VIIRS SDR XML input in sequence"
    
    start_time= time.time()
    LOG.debug("<< run_xml_file  >>")
    
    error_files=[]
    ran_ok=False
    cmd = [ADL_VIIRS_SDR, xml_file]
    if setup_only:
        print ' '.join(cmd)
    else:
        status_line('Executing %r with WORK_DIR=%r ' % (cmd, work_dir))

        LOG.debug('additional environment variables: %s' % repr(additional_env))

        try:
                
            pid = sh(cmd,env=env(**additional_env), cwd=work_dir)
            LOG.debug("%r ran as pid %d" % (cmd, pid))
            ran_ok=checkADLLogForSuccess(work_dir,pid,xml_file, remove_list)
        except CalledProcessError as oops:
            pid = getattr(oops,'pid', None)
            ran_ok=checkADLLogForSuccess(work_dir,pid,xml_file, remove_list)

            LOG.debug(traceback.format_exc())
            LOG.error('ProSdrViirsController.exe failed on %r: %r. Continuing...' % (xml_file, oops))
            error_files.append(xml_file)
            
 
    if error_files:
        LOG.warning('Had problems running these XML files: %r' % (error_files,))
        
        
       
    end_time=time.time() 
    elapsed_time = end_time-start_time
    LOG.debug("<< Run VIIRS SDR Elapsed time: %f >>"%(elapsed_time))
                    

    return ran_ok 
    
def build_product_name_maps():
    """ Read ADL short name to product map
    and create dictionaries for name convertion
    """

    tree= ET.parse(DDS_PRODUCT_FILE)
    p = tree.find("group/group") 
    
    links = list(p.iter("config"))   # Returns list of all links

    for i in links:
        productid=i.find('configValue').text
        shortname=i.find('name').text
        shortname=shortname.replace("_NPP","")
        
        PRODUCTID_2_SHORTNAME[productid] = shortname
        SHORTNAME_2_PRODUCTID[shortname]= productid


def add_geo_attribute_to_h5(work_dir, gran_id,product_dictionaries) :
    """ Adds GEO_REF property to produced granules
    """
     
    start_time= time.time()
    LOG.debug("<< add_geo_attribute_to_h5  >>")
    
   
    added=0

    LOG.debug("Add GEO_REF property to : "+gran_id+" products.")

    for short_name in  sorted ( ADL_VIIRS_SDR_PRODUCT_SHORTNAMES) :

        
# must require granule id but I do not know it is used.  
#        gran_list = list(skim_dir(work_dir,required_keys=CHECK_REQUIRED_KEYS,N_Collection_Short_Name=short_name,N_Granule_ID=gran_id))
       
        if optimize == True :
            gran_list = product_dictionaries[short_name]
        else :
            gran_list = list(skim_dir(work_dir,required_keys=CHECK_REQUIRED_KEYS,N_Collection_Short_Name=short_name,N_Granule_ID=gran_id))
        

        if len(gran_list) == 1 :
            LOG.debug("Granule list length: "+ str(len( gran_list )))
            for it in  gran_list :
                try :
                    sdr_name=SHORTNAME_2_PRODUCTID[short_name]
                    dObj=it[OBSERVE_TIME ]
                
                    time_id_str  = dObj.strftime("_npp_d%Y%m%d_t%H%M%S")
             
                    #  Name Example: SVI04_npp_d20111121_t1805421_e1807062_b00346_c20120127203200212753_cspp_dev.h5
                    sdr_name=sdr_name+time_id_str+"*.h5"  
                    files = glob.glob(os.path.join(work_dir, sdr_name))
                     
                    if  len( files ) == 1 :
                        fullname= files[0]
                        LOG.debug("GEO Product: "+fullname+" produced.") 
                        LOG.debug("Oberserved : " +time_id_str+" short_name: " + short_name) 
                        
                        adl_geo_ref.write_geo_ref(fullname)
                        added = added + 1                       
                except KeyError:  
                    LOG.debug("Key error on blob property "+'ObservedStartTime' +" or "+short_name)
               
    if  added > 0  :
            LOG.info("Added N_GEO_Ref properties to "+str(added)+" files.")
    
       
    end_time=time.time() 
    elapsed_time = end_time-start_time
    LOG.debug("<< Add GEO Elapsed time: %f >>"%(elapsed_time))
                    






 
def check_for_intermediate(work_dir, gran, remove_list,product_dictionaries) :
    """ Checks that proper products were produced.
    If product h5 file was produced blob and asc files are deleted
    If product was not produced files are left alone
    """
    # look for all products that should have been produced
    # convert the products to H5
   
     
    start_time= time.time()
    LOG.debug("<< check_for_intermediate  >>")
    
   
    
    LOG.info("Check for Granule: "+gran['N_Granule_ID']+" intermediate products.")

    gran_id=gran['N_Granule_ID']
    problem=True
    for short_name in  sorted ( ADL_VIIRS_SDR_intermediate_SHORTNAMES)  :
# must require granule id but I do not know it is used.
 
#        gran_list = list(skim_dir(work_dir,required_keys=CHECK_REQUIRED_KEYS,N_Collection_Short_Name=short_name,N_Granule_ID=gran_id))
        
        if optimize == True :
            gran_list = product_dictionaries[short_name]
        else :
            gran_list = list(skim_dir(work_dir,required_keys=CHECK_REQUIRED_KEYS,N_Collection_Short_Name=short_name,N_Granule_ID=gran_id))
        
        
        
        
        
        if len( gran_list )  < 1 :
            LOG.debug("No intermidate, Granule: "+gran_id+" Product: "+short_name)
        
        for it in  gran_list : 
                " Night passes do not create a blob file for all asc files."
                LOG.debug("Intermediate "+it['_filename']+" short "+short_name)
                ascname="noproperty"
                blobname="noproperty"
                        
                try :
                        blobname=it['BlobPath']   
                except KeyError:  
                        LOG.debug("Key error on blob property, Granule: "+gran_id+" Product: "+short_name)

                try :
                        ascname=it['_filename']
                except KeyError:  
                        LOG.debug("Key error on asc property")
                        LOG.debug("Blob: "+blobname )
                        b=blobname.split( ".")    
                        ascname=b[0]+".asc";


                               
                if  os.path.exists(ascname) :
                    LOG.debug("Append:" + ascname)
                    remove_list.append( ascname )
                            
                if  os.path.exists(blobname) :
                    LOG.debug("Append:" + blobname)
                    remove_list.append( blobname )
                  
                try :    
                    # now check for any h5 that might have been produced
                    sdr_name=SHORTNAME_2_PRODUCTID[short_name]
                    dObj=it[OBSERVE_TIME ]
                
                    time_id_str  = dObj.strftime("_npp_d%Y%m%d_t%H%M%S")
             
                    #  Name Example: SVI04_npp_d20111121_t1805421_e1807062_b00346_c20120127203200212753_cspp_dev.h5
                    sdr_name=sdr_name+time_id_str+"*.h5"  
                    files = glob.glob(os.path.join(work_dir, sdr_name))
                    for file in files:
                        if os.path.exists( file ) :
                            remove_list.append( file )
                        
                
                except KeyError:

                    i=1
                    
   
    end_time=time.time() 
    elapsed_time = end_time-start_time
    LOG.debug("<< Check Temporary Elapsed time: %f >>"%(elapsed_time))
                    
    return problem

def remove_inputs(work_dir) :
    viirs_inputs = "VIIRS-SCIENCE-RDR,SPACECRAFT-DIARY-RDR" 
    for input_name in viirs_inputs:
        pattern="*"+input_name
        for fn in glob.glob(os.path.join(work_dir, pattern)):
            b=fn.split( ".")            
            blob_id=b[0]+".asc"; 
            if os.path.exists(blob_id) :          
                os.remove(blob_id)
                
            blob_id=b[0]+".asc_not"; 
            if os.path.exists(blob_id) :          
                os.remove(blob_id)
            if os.path.exists(fn ) :
                os.remove(fn)
                
def remove_duplicate_asc(work_dir) :
    """
    remove any asc files that were ignored by the system
    """
    LOG.debug("Remove extra asc: "+work_dir)
    pattern="*.asc"
    for fn in glob.glob(os.path.join(work_dir, pattern)):
        LOG.debug("file: "+fn)
        if os.path.exists(fn ) :
            LOG.debug("Remove Dup: "+fn)
            os.remove(fn)               
                
    
def setup_directories(work_dir,anc_dir):
    # create work directory
    
    "Create the working directory and a subdirectory for the logs"
    
    if not os.path.isdir(work_dir):
        LOG.info('Creating directory %s' % work_dir)
        os.makedirs(work_dir)

        
    log_dir = os.path.join(work_dir, 'log')
    if not os.path.isdir(log_dir):
        LOG.info('Creating log directory %s' % log_dir)
        os.makedirs(log_dir)        

    if not os.path.exists( anc_dir ) :
        os.mkdir(anc_dir)

        

def find_granules_and_build_xml(work_dir):
    # read through ascii metadata and build up information table
    LOG.info('sifting through metadata to find VIIRS SDR processing candidates')
     
    start_time= time.time()
    LOG.debug("<< find_granules_and_build_xml  >>")
    
   
    granules_to_process=[]
    for gran in sift_metadata_for_viirs_sdr(work_dir) :
        granules_to_process.append( gran )
     
    LOG.debug(', '.join(x['N_Granule_ID'] for x in granules_to_process))
    if not granules_to_process:
        LOG.error("Found no granules to process!")
        return []
    else :
        LOG.info("Found %d granules to process."%(len(granules_to_process)))
        
       
    end_time=time.time() 
    elapsed_time = end_time-start_time
    LOG.debug("<< Find Granules Elapsed time: %f >>"%(elapsed_time))
                    

    
    return granules_to_process
       

def stage_the_ancillary_data() :
    """
    Stage all the ancillary data needed for the granules
    """
                 

def viirs_sdr(work_dir, h5_names, setup_only=False, cleanup=True, allow_cache_update=True,trim_granules=False):
    "process VIIRS RDR data in HDF5 format (any aggregation) to SDRs in an arbitrary work directory"

  
    global optimize


    ' Toggle to delete all files except h5 after executuion'
    cleanup_after_running = cleanup
    anc_dir=os.path.join(work_dir,ANCILLARY_SUB_DIR)
     
    "Create the working directory and a subdirectory for the logs" 
    setup_directories(work_dir, anc_dir)
    
    status_line("Unpack the supplied inputs")
    num_unpacking_problems =  unpack_h5s(work_dir, h5_names) 
     
    status_line("Search through the inputs for legal granule combinations")
    granules_to_process = find_granules_and_build_xml(work_dir)
    
    # used to create H5 file names
    build_product_name_maps()
  
    ####   EVALUATE AND RETRIEVE ANCILLARY DATA ##################
    status_line("Link the required ancillary data into the workspace")
        


    search_dirs = [CSPP_RT_ANC_CACHE_DIR if CSPP_RT_ANC_CACHE_DIR is not None else anc_dir] + [CSPP_RT_ANC_PATH]  
    problems_detected=0

    print search_dirs

    good_granules=0
    processed=0
    problem=False
    environment_error=False
    try:

        if len( granules_to_process ) > 0 :
            # get list of dynamic ancillary files.  Servicing may pull files from remote server.
            dynamic_ancillary_file_list =list()
            
            if allow_cache_update == True :
                LOG.info("Updating ancillary cache")
                dynamic_ancillary_file_list = adl_anc_retrieval.service_remote_ancillary(work_dir,granules_to_process,adl_anc_retrieval.kPOLAR)
                dynamic_ancillary_file_list2 = adl_anc_retrieval.service_remote_ancillary(work_dir,granules_to_process,adl_anc_retrieval.kTLE)
                for src_path in dynamic_ancillary_file_list2:
                    dynamic_ancillary_file_list.append( src_path )
        
            # get the static ancillary files needed
            ancillary_files_neeeded=anc_files_needed(ADL_VIIRS_ANC_GLOBS, search_dirs, granules_to_process)
        
            # create list of all ancillary needed
            for src_path in ancillary_files_neeeded:
                dynamic_ancillary_file_list.append( src_path )
        
            # link all the ancillary files to the ancillary directory.
            link_ancillary_to_work_dir(anc_dir, dynamic_ancillary_file_list)
        
        ##########  RUN THE VIIRS SDR  ##################3
        files_to_remove=[]
        file_that_will_be_removed  = []

        total_start_time= time.time()
        
        if trim_granules == True :
            # remove first and last granule
            item = granules_to_process.pop() # last item
            item = granules_to_process.pop(0) # first item


        for gran in granules_to_process:
            
            start_time= time.time()
            LOG.debug("<< RUN A GRANULE  >>")
           
            
            " For VIIRS files from previous granule need to stay around for next"
            " However we can remove then after they have been used to reduce parse time"
            " For ADL and error checking"
            processed = processed + 1
            if  cleanup_after_running == True :
                for file_to_remove in files_to_remove:
                    LOG.debug("Remove: "+ file_to_remove)
                    if os.path.exists(file_to_remove)   :
                        os.remove(file_to_remove)
                    else :
                        LOG.debug("Unable to remove:"+file_to_remove)
                
           
            files_to_remove            = file_that_will_be_removed
            file_that_will_be_removed  = []
            
            LOG.debug("Process: "+ str( gran))
            viirs_sdr_xml=generate_sdr_xml(work_dir,gran)
            file_that_will_be_removed.append(viirs_sdr_xml)    
 
            ran_ok = run_xml_file(work_dir, viirs_sdr_xml,file_that_will_be_removed,gran, setup_only = setup_only ,WORK_DIR=work_dir,LINKED_ANCILLARY=ANCILLARY_SUB_DIR,ADL_HOME=ADL_HOME,CSPP_RT_ANC_TILE_PATH=CSPP_RT_ANC_TILE_PATH)
            if ran_ok == False :
                LOG.info("Log indicates some problems.")
   
            LOG.debug("Checking that output granule blobs exist for granule: "+ gran['N_Granule_ID'])
            ################   Check For Products/ Errors  #####################3
            gran_id=gran['N_Granule_ID']
            
            
            product_dictionaries = None
            
            
            if optimize == True :
                product_dictionaries = _skim_dir_collections(work_dir, required_keys=CHECK_REQUIRED_KEYS,N_Granule_ID=gran_id)
            

            problem = check_for_products(work_dir, gran_id,file_that_will_be_removed,product_dictionaries)
            ###############   PATCH THE OUTPUTS ####################3  

            add_geo_attribute_to_h5(work_dir, gran_id,product_dictionaries)     

            check_for_intermediate(work_dir,gran,file_that_will_be_removed,product_dictionaries)
            
            
            stat=""
            if problem == False :
                good_granules += 1
#                status_line( 'Processing of %s Completed without problems'%viirs_sdr_xml)
                stat='Granule %s Completed without problems. '%gran_id
                
            if ran_ok == False or problem == True  :   
#               status_line( 'Processing of %s Completed with problems'%viirs_sdr_xml)
                stat='Granule %s Completed with problems. '%gran_id

                " Do not clean up files "
                problems_detected+=1
                
            status_line('%s %d out of %d granules processed, %d successfully'%(stat,processed,len(granules_to_process),good_granules))
            
            end_time=time.time() 
            elapsed_time = end_time-start_time
            LOG.debug("<< GRANULE Elapsed time: %f >>"%(elapsed_time))
           
            
        total_end_time=time.time() 
        elapsed_time = total_end_time-total_start_time
        LOG.debug("<< TOTAL Elapsed time: %f >>"%(elapsed_time))
            
           
                    
        if len( granules_to_process ) > 0 and cleanup_after_running == True :
            for file in files_to_remove + file_that_will_be_removed:      
                if os.path.exists( file)  :
                    LOG.debug( "Remove: "+file )
                    os.remove( file )
                else :
                    LOG.debug("Unable to remove:"+file)
        
        status_line('Final %d out of %d granules processed, %d successfully'%(processed,len(granules_to_process),good_granules))
        if cleanup_after_running == True  :
            remove_inputs(work_dir)
            remove_duplicate_asc(work_dir)
        
    except EnvironmentError:
        environment_error=True

    if cleanup_after_running == True  :
        dir_to_go=os.path.join(work_dir,"log")
        shutil.rmtree(dir_to_go)
        dir_to_go=os.path.join(work_dir,ANCILLARY_SUB_DIR)
        shutil.rmtree(dir_to_go)

    gran_failure =  len(granules_to_process) - good_granules
    noncritical_problem = 0

    rc = get_return_code(num_unpacking_problems, len(granules_to_process), gran_failure, noncritical_problem, environment_error)   
        
    return rc  



def _test_anc_names():
    "list ancillary files that would be processed"
    from pprint import pprint
    print "ancillary files:"
    pprint(list(anc_files_needed(ADL_VIIRS_ANC_GLOBS, ADL_CSPP_RT_ANC_PATH, [])))


def _test_sdr_granules(work_dir):
    "list granules we'd generate XML for"
    granules_to_process = list(sift_metadata_for_viirs_sdr(work_dir))
    from pprint import pprint
    pprint([x['N_Granule_ID'] for x in granules_to_process])
    return granules_to_process


def viirs_sdr_run(work_dir, h5_names,compress=False,aggregate=False,all_geo=False, setup_only=False, cleanup=True, allow_cache_update=True,trim_granules=False) :
    
    check_env(work_dir)
    
    if all_geo == True :
        ADL_VIIRS_GEO_PRODUCT_SHORTNAMES.extend(OPTIONAL_GEO_PRODUCTS)
    else :
        ADL_VIIRS_SDR_intermediate_SHORTNAMES.extend(OPTIONAL_GEO_PRODUCTS)
 
    rc = viirs_sdr( work_dir, h5_names ,cleanup=cleanup)

    if compress == True :
        repack_products(work_dir, ADL_VIIRS_SDR_PRODUCT_SHORTNAMES )
        repack_products(work_dir,  ADL_VIIRS_GEO_PRODUCT_SHORTNAMES )

    if aggregate == True :
        aggregate_products(work_dir,  ADL_VIIRS_GEO_PRODUCT_SHORTNAMES)
        aggregate_products(work_dir, ADL_VIIRS_SDR_PRODUCT_SHORTNAMES)
        add_geo_attribute_to_aggregates(work_dir, ADL_VIIRS_SDR_PRODUCT_SHORTNAMES)
    
    return rc


def main():
    """ Run Viirs SDR proccessing
    """   
    

    import argparse
    desc = """Build VIIRS SDR ProSdrViirsDBController work directory and run VIIRS SDR."""
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

    
    parser.add_argument('-v', '--verbosity', action="count", default=0,
                    help='each occurrence increases verbosity 1 level through ERROR-WARNING-INFO-DEBUG')
    
    parser.add_argument('filenames', metavar='filename', type=str, nargs='+',
                   help='HDF5 VIIRS RDR file/s to process')

    args = parser.parse_args()

    levels = [logging.ERROR, logging.WARN, logging.INFO, logging.DEBUG]
    level = levels[args.verbosity if args.verbosity<4 else 3]

    docleanup = True
    
    
   
    work_dir= check_and_convert_path("WORK_DIR",args.work_dir)
    d = dt.date.today()
    timestamp = d.isoformat()
    logname= "viirs_sdr."+timestamp+".log"
    logfile= os.path.join(work_dir, logname )
    configure_logging(level,FILE=logfile)
    if  args.debug == True:
        docleanup = False
        
    compress=False
    if args.zip == True :
        compress=True
        
    aggregate=False
    if args.aggregate == True :
        if len ( args.filenames ) > 1 :
            LOG.error( "Only one input allowed with aggregate option");
            sys.exit(2)
        aggregate=True

        
  
        
    
    LOG.debug("Clean up: "+str(docleanup))


    

    LOG.info('CSPP execution work directory is %r' % work_dir)
    
    if args.test:
        check_env(work_dir)
        grans = _test_sdr_granules(work_dir)
        if grans:
            LOG.debug('building XML files')
   
        sys.exit(2)

    if not args: 
        parser.print_help()
        return 1



    rc = viirs_sdr_run (work_dir, args.filenames,compress=compress,aggregate=aggregate,all_geo=args.geo, setup_only=False, cleanup=docleanup, allow_cache_update=not args.local)


    if rc == 0 and not args.debug :
        os.remove(logfile)

    return rc


if __name__=='__main__':
    sys.exit(main())


