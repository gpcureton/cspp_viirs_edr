#!/usr/bin/env python
# encoding: utf-8
"""
OpticalDepth.py

 * DESCRIPTION:  Class to granulate the ViirsAncOpticalDepth data product 

Created by Geoff Cureton on 2009-04-04.
Copyright (c) 2011 University of Wisconsin SSEC. All rights reserved.
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
import re
import uuid
from glob import glob
from time import time
from datetime import datetime,timedelta

from scipy import round_

import numpy as np
from numpy import ma
import copy
from bisect import bisect_left,bisect_right

import ctypes
from numpy.ctypeslib import ndpointer

import ViirsData

from NAAPStoBlob import NAAPSclass

# skim and convert routines for reading .asc metadata fields of interest
import adl_blob
import adl_asc
from adl_asc import skim_dir, contiguous_granule_groups, granule_groups_contain, effective_anc_contains,_eliminate_duplicates,_is_contiguous, RDR_REQUIRED_KEYS, POLARWANDER_REQUIRED_KEYS
from adl_common import ADL_HOME, CSPP_RT_ANC_PATH, CSPP_RT_ANC_CACHE_DIR, COMMON_LOG_CHECK_TABLE

# every module should have a LOG object
try :
    sourcename= file_Id.split(" ")
    LOG = logging.getLogger(sourcename[1])
except :
    LOG = logging.getLogger('OpticalDepth')

from Utils import getURID, getAscLine, getAscStructs, findDatelineCrossings, shipOutToFile

class OpticalDepth() :

    def __init__(self,inDir=None, sdrEndian=None, ancEndian=None):
        self.collectionShortName = 'VIIRS-ANC-Optical-Depth-Mod-Gran'
        self.xmlName = 'VIIRS_ANC_OPTICAL_DEPTH_MOD_GRAN.xml'
        self.blobDatasetName = 'aotGrid'
        self.dataType = 'float32'
        self.sourceType = 'NAAPS_ANC_Int'
        self.sourceList = ['']
        self.trimObj = ViirsData.ViirsTrimTable()

        if inDir is None :
            self.inDir = path.abspath(path.curdir)
        else :
            self.inDir = inDir

        if sdrEndian is None :
            self.sdrEndian = adl_blob.LITTLE_ENDIAN
        else :
            self.sdrEndian = sdrEndian

        if ancEndian is None :
            self.ancEndian = adl_blob.LITTLE_ENDIAN
        else :
            self.ancEndian = ancEndian


    def ingest(self,ancBlob=None):
        '''
        Ingest the ancillary dataset.
        '''
        self.gridData = getattr(ancBlob,self.blobDatasetName).astype(self.dataType)

        ###################################
        #-- From adl_viirs_edr_masks.py --#
        ###################################

        # Transcode the NAAPS GRIB files into NAAPS global grid blob files
        #gridBlobFiles = create_NAAPS_gridBlobs(naapsFiles)
        #gridBlobFiles = glob(path.join(work_dir,'*.NAAPS-ANC-Int'))
        #LOG.debug("gridBlobFiles: %r" % (gridBlobFiles)

        # Granulate the global grid NAAPS blob files
        #granBlobFiles = granulate_NAAPS_gridBlobs(work_dir,anc_granules_to_process,gridBlobFiles)
        #LOG.debug("granBlobFiles: %r" % (granBlobFiles)


    def setGeolocationInfo(self,dicts):
        '''
        Populate this class instance with the geolocation data for a single granule
        '''
        # Set some environment variables and paths
        CSPP_RT_HOME = os.getenv('CSPP_RT_HOME')
        ANC_SCRIPTS_PATH = path.join(CSPP_RT_HOME,'viirs')
        CSPP_RT_ANC_CACHE_DIR = os.getenv('CSPP_RT_ANC_CACHE_DIR')
        csppPython = os.getenv('PY')
    
        ADL_ASC_TEMPLATES = path.join(ANC_SCRIPTS_PATH,'asc_templates')

        # Collect some data from the geolocation dictionary
        self.geoDict = dicts
        URID = dicts['URID']
        geo_Collection_ShortName = dicts['N_Collection_Short_Name']
        N_Granule_ID = dicts['N_Granule_ID']
        ObservedStartTimeObj = dicts['ObservedStartTime']
        geoFiles = glob('%s/%s*' % (self.inDir,URID))
        geoFiles.sort()

        LOG.debug("\n###########################")
        LOG.debug("  Geolocation Information  ")
        LOG.debug("###########################")
        LOG.debug("N_Granule_ID : %r" % (N_Granule_ID))
        LOG.debug("ObservedStartTime : %s" % (ObservedStartTimeObj.__str__()))
        LOG.debug("N_Collection_Short_Name : %s" %(geo_Collection_ShortName))
        LOG.debug("URID : %r" % (URID))
        LOG.debug("geoFiles : %r" % (geoFiles))
        LOG.debug("###########################\n")

        # Do we have terrain corrected geolocation?

        terrainCorrectedGeo = True if 'GEO-TC' in geo_Collection_ShortName else False

        # Do we have long or short style geolocation field names?

        if (geo_Collection_ShortName=='VIIRS-MOD-GEO-TC' or geo_Collection_ShortName=='VIIRS-MOD-RGEO') :
            longFormGeoNames = True
            LOG.debug("We have long form geolocation names")
        elif (geo_Collection_ShortName=='VIIRS-MOD-GEO' or geo_Collection_ShortName=='VIIRS-MOD-RGEO-TC') :
            LOG.debug("We have short form geolocation names")
            longFormGeoNames = False
        else :
            LOG.error("Invalid geolocation shortname: %s",geo_Collection_ShortName)
            return -1

        # Get the geolocation xml file

        geoXmlFile = "%s.xml" % (string.replace(geo_Collection_ShortName,'-','_'))
        geoXmlFile = path.join(ADL_HOME,'xml/VIIRS',geoXmlFile)
        if path.exists(geoXmlFile):
            LOG.debug("We are using for %s: %s,%s" %(geo_Collection_ShortName,geoXmlFile,geoFiles[0]))

        # Open the geolocation blob and get the latitude and longitude

        endian = self.sdrEndian

        geoBlobObj = adl_blob.map(geoXmlFile,geoFiles[0], endian=endian)
        geoBlobArrObj = geoBlobObj.as_arrays()

        # Get scan_mode to find any bad scans

        scanMode = geoBlobArrObj.scan_mode[:]
        badScanIdx = np.where(scanMode==254)[0]
        LOG.debug("Bad Scans: %r" % (badScanIdx))

        # Detemine the min, max and range of the latitude and longitude, 
        # taking care to exclude any fill values.

        if longFormGeoNames :
            latitude = getattr(geoBlobArrObj,'latitude').astype('float')
            longitude = getattr(geoBlobArrObj,'longitude').astype('float')
        else :
            latitude = getattr(geoBlobArrObj,'lat').astype('float')
            longitude = getattr(geoBlobArrObj,'lon').astype('float')

        latitude = ma.masked_less(latitude,-800.)
        latMin,latMax = np.min(latitude),np.max(latitude)
        latRange = latMax-latMin

        longitude = ma.masked_less(longitude,-800.)
        lonMin,lonMax = np.min(longitude),np.max(longitude)
        lonRange = lonMax-lonMin

        LOG.debug("min,max,range of latitide: %f %f %f" % (latMin,latMax,latRange))
        LOG.debug("min,max,range of longitude: %f %f %f" % (lonMin,lonMax,lonRange))

        # Determine the latitude and longitude fill masks, so we can restore the 
        # fill values after we have scaled...

        latMask = latitude.mask
        lonMask = longitude.mask

        # Check if the geolocation is in radians, convert to degrees
        if 'RGEO' in geo_Collection_ShortName :
            LOG.debug("Geolocation is in radians, convert to degrees...")
            latitude = np.degrees(latitude)
            longitude = np.degrees(longitude)
        
            latMin,latMax = np.min(latitude),np.max(latitude)
            latRange = latMax-latMin
            lonMin,lonMax = np.min(longitude),np.max(longitude)
            lonRange = lonMax-lonMin

            LOG.debug("New min,max,range of latitude: %f %f %f" % (latMin,latMax,latRange))
            LOG.debug("New min,max,range of longitude: %f %f %f" % (lonMin,lonMax,lonRange))

        # Restore fill values to masked pixels in geolocation

        geoFillValue = self.trimObj.sdrTypeFill['VDNE_FLOAT64_FILL'][latitude.dtype.name]
        latitude = ma.array(latitude,mask=latMask,fill_value=geoFillValue)
        self.latitude = latitude.filled()

        geoFillValue = self.trimObj.sdrTypeFill['VDNE_FLOAT64_FILL'][longitude.dtype.name]
        longitude = ma.array(longitude,mask=lonMask,fill_value=geoFillValue)
        self.longitude = longitude.filled()

        # Record the corners, taking care to exclude any bad scans...
        nDetectors = 16
        firstGoodScan = np.where(scanMode<=2)[0][0]
        lastGoodScan = np.where(scanMode<=2)[0][-1]
        firstGoodRow = firstGoodScan * nDetectors
        lastGoodRow = lastGoodScan * nDetectors + nDetectors - 1

        latCrnList = [latitude[firstGoodRow,0],latitude[firstGoodRow,-1],latitude[lastGoodRow,0],latitude[lastGoodRow,-1]]
        lonCrnList = [longitude[firstGoodRow,0],longitude[firstGoodRow,-1],longitude[lastGoodRow,0],longitude[lastGoodRow,-1]]

        # Check for dateline/pole crossings
        num180Crossings = findDatelineCrossings(latCrnList,lonCrnList)
        LOG.info("We have %d dateline crossings."%(num180Crossings))

        # Copy the geolocation information to the class object
        self.latMin    = latMin
        self.latMax    = latMax
        self.latRange  = latRange
        self.lonMin    = lonMin
        self.lonMax    = lonMax
        self.lonRange  = lonRange
        self.scanMode  = scanMode
        self.latCrnList  = latCrnList
        self.lonCrnList  = lonCrnList
        self.num180Crossings  = num180Crossings

        # Parse the geolocation asc file to get struct information which will be 
        # written to the ancillary asc files

        geoAscFileName = path.join(self.inDir,URID+".asc")
        LOG.debug("\nOpening %s..." % (geoAscFileName))

        geoAscFile = open(geoAscFileName,'rt')

        #RangeDateTimeStr =  _getAscLine(geoAscFile,"RangeDateTime")
        self.RangeDateTimeStr =  getAscLine(geoAscFile,"ObservedDateTime")
        self.RangeDateTimeStr =  string.replace(self.RangeDateTimeStr,"ObservedDateTime","RangeDateTime")
        self.GRingLatitudeStr =  getAscStructs(geoAscFile,"GRingLatitude",12)
        self.GRingLongitudeStr = getAscStructs(geoAscFile,"GRingLongitude",12)

        geoAscFile.close()


    def _grid2Gran_bilinearInterp(self,dataLat, dataLon, gridData, gridLat, gridLon):
        '''Granulates a gridded dataset using an input geolocation'''

        nData = np.int64(dataLat.size)
        gridRows = np.int32(gridLat.shape[0])
        gridCols = np.int32(gridLat.shape[1])

        data = np.ones(np.shape(dataLat),dtype=np.float64)* -999.9
        dataIdx  = np.ones(np.shape(dataLat),dtype=np.int64) * -99999

        CSPP_RT_HOME = os.getenv('CSPP_RT_HOME')
        ANC_SCRIPTS_PATH = path.join(CSPP_RT_HOME,'viirs')

        libFile = path.join(ANC_SCRIPTS_PATH,'libgriddingAndGranulation.so.1.0.1')
        LOG.debug("Gridding and granulation library file: %s" % (libFile))
        lib = ctypes.cdll.LoadLibrary(libFile)
        grid2gran_bilinearInterp = lib.grid2gran_bilinearInterp
        grid2gran_bilinearInterp.restype = None
        grid2gran_bilinearInterp.argtypes = [
                ndpointer(ctypes.c_double,ndim=1,shape=(nData),flags='C_CONTIGUOUS'),
                ndpointer(ctypes.c_double,ndim=1,shape=(nData),flags='C_CONTIGUOUS'),
                ndpointer(ctypes.c_double,ndim=1,shape=(nData),flags='C_CONTIGUOUS'),
                ctypes.c_int64,
                ndpointer(ctypes.c_double,ndim=2,shape=(gridRows,gridCols),flags='C_CONTIGUOUS'),
                ndpointer(ctypes.c_double,ndim=2,shape=(gridRows,gridCols),flags='C_CONTIGUOUS'),
                ndpointer(ctypes.c_double,ndim=2,shape=(gridRows,gridCols),flags='C_CONTIGUOUS'),
                ndpointer(ctypes.c_int64,ndim=1,shape=(nData),flags='C_CONTIGUOUS'),
                ctypes.c_int32,
                ctypes.c_int32
                ]

        '''
        int snapGrid_ctypes(double *lat, 
                        double *lon, 
                        double *data, 
                        long nData, 
                        double *gridLat,
                        double *gridLon,
                        double *gridData,
                        long *gridDataIdx,
                        int nGridRows,
                        int nGridCols
                        )
        '''

        LOG.debug("Calling C routine grid2gran_bilinearInterp()...")

        retVal = grid2gran_bilinearInterp(dataLat,
                           dataLon,
                           data,
                           nData,
                           gridLat,
                           gridLon,
                           gridData,
                           dataIdx,
                           gridRows,
                           gridCols)

        LOG.debug("Returning from C routine grid2gran_bilinearInterp()")

        return data,dataIdx


    def granulate(self,ANC_objects):
        '''
        Granulate the ancillary dataset.
        '''
        LOG.info("Granulating %s ..." % (self.collectionShortName))

        degInc = 0.5

        lats = np.arange(361.)*degInc - 90.
        lons = np.arange(720.)*degInc - 180.
        latitude = self.latitude
        longitude = self.longitude

        # Flip so that lats are (-90 ... 90)
        gridData = self.gridData[::-1,:]

        if self.num180Crossings != 2 :

            gridData = np.roll(gridData,360)
            gridLon,gridLat = np.meshgrid(lons,lats)

            LOG.info("start,end NAAPS Grid Latitude values : %f,%f"%(gridLat[0,0],gridLat[-1,0]))
            LOG.info("start,end NAAPS Grid Longitude values : %f,%f"%(gridLon[0,0],gridLon[0,-1]))

        else :

            negLonIdx = np.where(lons<0)
            lons[negLonIdx] += 360.
            lons = np.roll(lons,360)
            gridLon,gridLat = np.meshgrid(lons,lats)

            longitudeNegIdx = np.where(longitude < 0.)
            longitude[longitudeNegIdx] += 360.

            LOG.info("start,end NAAPS Grid Latitude values : %f,%f"%(gridLat[0,0],gridLat[-1,0]))
            LOG.info("start,end NAAPS Grid Longitude values : %f,%f"%(gridLon[0,0],gridLon[0,-1]))


        LOG.debug("min of gridData  = %r"%(np.min(gridData)))
        LOG.debug("max of gridData  = %r"%(np.max(gridData)))

        t1 = time()
        data,dataIdx = self._grid2Gran_bilinearInterp(np.ravel(latitude),
                                  np.ravel(longitude),
                                  gridData.astype(np.float64),
                                  gridLat,
                                  gridLon)
        t2 = time()
        elapsedTime = t2-t1
        LOG.info("Granulation took %f seconds for %d points" % (elapsedTime,latitude.size))

        data = data.reshape(latitude.shape)
        dataIdx = dataIdx.reshape(latitude.shape)

        LOG.debug("Shape of granulated %s data is %s" % (self.collectionShortName,np.shape(data)))
        LOG.debug("Shape of granulated %s dataIdx is %s" % (self.collectionShortName,np.shape(dataIdx)))

        # Moderate resolution trim table arrays. These are 
        # bool arrays, and the trim pixels are set to True.
        modTrimMask = self.trimObj.createModTrimArray(nscans=48,trimType=bool)

        # Fill the required pixel trim rows in the granulated NCEP data with 
        # the ONBOARD_PT_FILL value for the correct data type

        fillValue = self.trimObj.sdrTypeFill['ONBOARD_PT_FILL'][data.dtype.name]        
        data = ma.array(data,mask=modTrimMask,fill_value=fillValue)
        self.data = data.filled()


    def shipOutToFile(self):
        ''' Pass the current class instance to this Utils method to generate 
            a blob/asc file pair from the input ancillary data object.'''

        # Set some environment variables and paths
        CSPP_RT_HOME = os.getenv('CSPP_RT_HOME')
        ANC_SCRIPTS_PATH = path.join(CSPP_RT_HOME,'viirs')
        CSPP_RT_ANC_CACHE_DIR = os.getenv('CSPP_RT_ANC_CACHE_DIR')
        ADL_ASC_TEMPLATES = path.join(ANC_SCRIPTS_PATH,'asc_templates')

        # Create new ANC ancillary blob, and copy granulated data to it

        endian = self.ancEndian
        xmlName = path.join(ADL_HOME,'xml/VIIRS',self.xmlName)

        # Create a new URID to be used in making the asc filenames

        URID_dict = getURID()

        URID = URID_dict['URID']
        creationDate_nousecStr = URID_dict['creationDate_nousecStr']
        creationDateStr = URID_dict['creationDateStr']

        # Create a new directory in the input directory for the new ancillary
        # asc and blob files

        blobDir = self.inDir

        ascFileName = path.join(blobDir,URID+'.asc')
        blobName = path.join(blobDir,URID+'.'+self.collectionShortName)

        LOG.debug("ascFileName : %s" % (ascFileName))
        LOG.debug("blobName : %s" % (blobName))

        # Create a new ancillary blob, and copy the data to it.
        newANCblobObj = adl_blob.create(xmlName, blobName, endian=endian, overwrite=True)
        newANCblobArrObj = newANCblobObj.as_arrays()

        # FIXME : Possible datasets are 'faot550' or 'aotSlant550', the latter of which requires the sensor
        #         zenith angle from the geolocation.
        blobData = getattr(newANCblobArrObj,'faot550')
        blobData[:,:] = self.data[:,:]

        # Make a new ANC asc file from the template, and substitute for the various tags

        ascTemplateFileName = path.join(ADL_ASC_TEMPLATES,"VIIRS-ANC_Template.asc")

        LOG.debug("Creating new asc file\n%s\nfrom template\n%s" % (ascFileName,ascTemplateFileName))
        
        ANC_fileList = self.sourceList
        for idx in range(len(ANC_fileList)) :
            ANC_fileList[idx] = path.basename(ANC_fileList[idx])
        ANC_fileList.sort()
        ancGroupRecipe = '    ("N_Anc_Filename" STRING EQ "%s")'
        ancFileStr = "%s" % ("\n").join([ancGroupRecipe % (str(files)) for files in ANC_fileList])

        LOG.debug("RangeDateTimeStr = %s\n" % (self.RangeDateTimeStr))
        LOG.debug("GRingLatitudeStr = \n%s\n" % (self.GRingLatitudeStr))
        LOG.debug("GRingLongitudeStr = \n%s\n" % (self.GRingLongitudeStr))

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
           line = line.replace("CSPP_ANC_COLLECTION_SHORT_NAME",self.collectionShortName)
           line = line.replace("CSPP_GRANULE_ID",self.geoDict['N_Granule_ID'])
           line = line.replace("CSPP_CREATIONDATETIME",creationDateStr)
           line = line.replace("  CSPP_RANGE_DATE_TIME",self.RangeDateTimeStr)
           line = line.replace("  CSPP_GRINGLATITUDE",self.GRingLatitudeStr)
           line = line.replace("  CSPP_GRINGLONGITUDE",self.GRingLongitudeStr)
           line = line.replace("    CSPP_ANC_SOURCE_FILES",ancFileStr)
           ascFile.write(line) 

        ascFile.close()
        ascTemplateFile.close()



    def retrieve_NAAPS_grib_files(geoDicts):
        ''' Download the NAAPS GRIB files which cover the dates of the geolocation files.'''

        CSPP_RT_HOME = os.getenv('CSPP_RT_HOME')
        ANC_SCRIPTS_PATH = path.join(CSPP_RT_HOME,'viirs')
        CSPP_RT_ANC_CACHE_DIR = os.getenv('CSPP_RT_ANC_CACHE_DIR')

        # FIXME : Fix rounding up of the seconds if the decisecond>=9.5 
        gribFiles = []
        for geoDict in geoDicts:
            #timeObj = geoDict['StartTime']
            timeObj = geoDict['ObservedStartTime']
            dateStamp = timeObj.strftime("%Y%m%d")
            seconds = repr(int(round(timeObj.second + float(timeObj.microsecond)/1000000.)))
            deciSeconds = int(round(float(timeObj.microsecond)/100000.))
            deciSeconds = repr(0 if deciSeconds > 9 else deciSeconds)
            startTimeStamp = "%s%s" % (timeObj.strftime("%H%M%S"),deciSeconds)

            #timeObj = geoDict['EndTime']
            timeObj = geoDict['ObservedEndTime']
            seconds = repr(int(round(timeObj.second + float(timeObj.microsecond)/1000000.)))
            deciSeconds = int(round(float(timeObj.microsecond)/100000.))
            deciSeconds = repr(0 if deciSeconds > 9 else deciSeconds)
            endTimeStamp = "%s%s" % (timeObj.strftime("%H%M%S"),deciSeconds)

            timeObj = geoDict['UnpackTime']
            unpackTimeStamp = timeObj.strftime("%Y%m%d%H%M%S%f")

            granuleName = "GMODO_npp_d%s_t%s_e%s_b00014_c%s.h5" % (dateStamp,startTimeStamp,endTimeStamp,unpackTimeStamp)

            try :
                LOG.info('Retrieving NAAPS files for %s ...' % (granuleName))
                cmdStr = '%s/cspp_retrieve_gdas_gfs.csh %s' % (ANC_SCRIPTS_PATH,granuleName)
                LOG.info('\t%s' % (cmdStr))
                args = shlex.split(cmdStr)
                LOG.debug('\t%s' % (repr(args)))

                procRetVal = 0
                procObj = subprocess.Popen(args,bufsize=0, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                procObj.wait()
                procRetVal = procObj.returncode

                procOutput = procObj.stdout.readlines()
                #procOutput = procObj.stdout.read()

                # FIXME : How to get this output to have linebreaks when using readlines()
                #LOG.debug(procOutput)
                
                for lines in procOutput:
                    if "GDAS/GFS file" in lines :
                        lines = string.replace(lines,'GDAS/GFS file 1: ','')
                        lines = string.replace(lines,'GDAS/GFS file 2: ','')
                        lines = string.replace(lines,'\n','')
                        gribFiles.append(lines)

                # TODO : On error, jump to a cleanup routine
                if not (procRetVal == 0) :
                    LOG.error('Retrieval of ancillary files failed for %s.' % (granuleName))
                    #sys.exit(procRetVal)

            except Exception, err:
                LOG.warn( "%s" % (str(err)))

        # Get a unique list of grib files that were fetched
        gribFiles.sort()
        gribFiles = dict(map(lambda i: (i,1),gribFiles)).keys()
        gribFiles.sort()

        for gribFile in gribFiles :
            LOG.info('Retrieved gribfile: %r' % (gribFile))

        return gribFiles


    def create_NAAPS_gridBlobs(gribFiles):
        '''Converts NAAPS GRIB files into NAAPS blobs'''

        from copy import deepcopy

        blobFiles = []

        CSPP_RT_HOME = os.getenv('CSPP_RT_HOME')
        ANC_SCRIPTS_PATH = path.join(CSPP_RT_HOME,'viirs')
        CSPP_RT_ANC_CACHE_DIR = os.getenv('CSPP_RT_ANC_CACHE_DIR')
        csppPython = os.getenv('PY')
        

        '''
        for files in gribFiles :
            gribPath = path.dirname(files)
            gribFile = path.basename(files)
            gribBlob = "%s_blob.le" % (gribFile)
            gribBlob = path.join(gribPath,gribBlob)
            LOG.debug("Creating NAAPS GRIB blob %s" % (gribBlob))

            if not path.exists(gribBlob):
                try :
                    LOG.info('Transcoding %s to %s ...' % (files,gribBlob))
                    NAAPSxml = path.join(ADL_HOME,'xml/ANC/NAAPS_ANC_Int.xml')

                    # Create the grib object and populate with the grib file data
                    NAAPSobj = NAAPSclass(gribFile=files)

                    # Write the contents of the NAAPSobj object to an ADL blob file
                    endian = adl_blob.LITTLE_ENDIAN
                    procRetVal = NAAPSclass.NAAPSgribToBlob_interpNew(NAAPSobj,NAAPSxml,gribBlob,endian=endian)

                    #blobFiles.append(gribBlob)
                    if not (procRetVal == 0) :
                        LOG.error('Transcoding of ancillary files failed for %s.' % (files))
                        sys.exit(procRetVal)
                    else :
                        LOG.info('Finished creating NAAPS GRIB blob %s' % (gribBlob))
                        blobFiles.append(gribBlob)

                except Exception, err:
                    LOG.warn( "%s" % (str(err)))
            else :
                LOG.info('Gridded NAAPS blob files %s exists, skipping.' % (gribBlob))
                blobFiles.append(gribBlob)
        '''

        blobFiles = []

        LOG.info('Returning NAAPS GRIB blob file names %r' % (blobFiles))
        return blobFiles


    def granulate_NAAPS_gridBlobs(inDir,geoDicts, gridBlobFiles):
        '''Granulates the input gridded blob files into the required NAAPS granulated datasets.'''

        global ancEndian 

        CSPP_RT_HOME = os.getenv('CSPP_RT_HOME')
        ANC_SCRIPTS_PATH = path.join(CSPP_RT_HOME,'masks')
        CSPP_RT_ANC_CACHE_DIR = os.getenv('CSPP_RT_ANC_CACHE_DIR')
        csppPython = os.getenv('PY')
        
        ADL_ASC_TEMPLATES = path.join(ANC_SCRIPTS_PATH,'asc_templates')
        
        # Collection shortnames of the required NAAPS ancillary datasets
        # FIXME : Poll ADL/cfg/ProEdrViirsCM_CFG.xml for this information

        masksCollShortNames = [
                               'VIIRS-ANC-Optical-Depth-Mod-Gran'
                              ]

        # Dictionary relating the required NAAPS collection short names and the 
        # NAAPS gridded blob dataset names

        NAAPS_shortNameToBlobName = {}
        NAAPS_shortNameToBlobName['VIIRS-ANC-Optical-Depth-Mod-Gran'] = 'aotGrid'

        NAAPS_shortNameToXmlName = {}
        NAAPS_shortNameToXmlName['VIIRS-ANC-Optical-Depth-Mod-Gran'] = 'VIIRS_ANC_OPTICAL_DEPTH_MOD_GRAN.xml'

        # Moderate resolution trim table arrays. These are 
        # bool arrays, and the trim pixels are set to True.

        trimObj = ViirsData.ViirsTrimTable()
        modTrimMask = trimObj.createModTrimArray(nscans=48,trimType=bool)

        # Open the NAAPS gridded blob file

        naapsXmlFile = path.join(ADL_HOME,'xml/ANC/NAAPS_ANC_Int.xml')
        # FIXME : Should be using two NAAPS blob files, and averaging
        gridBlobFile = gridBlobFiles[0]

        if path.exists(naapsXmlFile):
            LOG.debug("We are using for %s: %s,%s" %('NAAPS-ANC-Int',naapsXmlFile,gridBlobFile))
        
        endian = ancEndian

        naapsBlobObj = adl_blob.map(naapsXmlFile,gridBlobFile, endian=adl_blob.BIG_ENDIAN)
        naapsBlobArrObj = naapsBlobObj.as_arrays()

        LOG.debug("%s...\n%r" % (gridBlobFile,naapsBlobArrObj._fields))

        # Save the global grids of the required datasets into a dictionary...

        NAAPS_globalGridData = {}

        # Precipitable water
        blobDsetName = NAAPS_shortNameToBlobName['VIIRS-ANC-Optical-Depth-Mod-Gran']
        NAAPS_globalGridData['VIIRS-ANC-Optical-Depth-Mod-Gran'] = \
                getattr(naapsBlobArrObj,blobDsetName).astype('float')
        LOG.debug("Shape of dataset %s is %s" % (blobDsetName,np.shape(NAAPS_globalGridData['VIIRS-ANC-Optical-Depth-Mod-Gran'])))


        # Contruct a default 0.5 degree grid...

        degInc = 0.5
        grids = np.mgrid[-90.:90.+degInc:degInc,-180.:180.:degInc]
        gridLat,gridLon = grids[0],grids[1]

        # Loop through the geolocation files and granulate...

        for dicts in geoDicts :

            URID = dicts['URID']
            geo_Collection_ShortName = dicts['N_Collection_Short_Name']
            N_Granule_ID = dicts['N_Granule_ID']
            ObservedStartTimeObj = dicts['ObservedStartTime']
            geoFiles = glob('%s/%s*' % (inDir,URID))
            geoFiles.sort()

            LOG.debug("\n###########################")
            LOG.debug("  Geolocation Information  ")
            LOG.debug("###########################")
            LOG.debug("N_Granule_ID : %r" % (N_Granule_ID))
            LOG.debug("ObservedStartTime : %s" % (ObservedStartTimeObj.__str__()))
            LOG.debug("N_Collection_Short_Name : %s" %(geo_Collection_ShortName))
            LOG.debug("URID : %r" % (URID))
            LOG.debug("geoFiles : %r" % (geoFiles))
            LOG.debug("###########################\n")

            # Do we have terrain corrected geolocation?

            terrainCorrectedGeo = True if 'GEO-TC' in geo_Collection_ShortName else False

            # Do we have long or short style geolocation field names?

            if (geo_Collection_ShortName=='VIIRS-MOD-GEO-TC' or geo_Collection_ShortName=='VIIRS-MOD-RGEO') :
                longFormGeoNames = True
                LOG.debug("We have long form geolocation names")
            elif (geo_Collection_ShortName=='VIIRS-MOD-GEO' or geo_Collection_ShortName=='VIIRS-MOD-RGEO-TC') :
                LOG.debug("We have short form geolocation names")
                longFormGeoNames = False
            else :
                LOG.error("Invalid geolocation shortname: %s"% (geo_Collection_ShortName))
                return -1

            # Get the geolocation xml file

            geoXmlFile = "%s.xml" % (string.replace(geo_Collection_ShortName,'-','_'))
            geoXmlFile = path.join(ADL_HOME,'xml/VIIRS',geoXmlFile)
            if path.exists(geoXmlFile):
                LOG.debug("We are using for %s: %s,%s" %(geo_Collection_ShortName,geoXmlFile,geoFiles[0]))

            # Open the geolocation blob and get the latitude and longitude

            endian=sdrEndian

            geoBlobObj = adl_blob.map(geoXmlFile,geoFiles[0], endian=endian)
            geoBlobArrObj = geoBlobObj.as_arrays()

            # If we have the terrain corrected geolocation, get the terrain height

            if terrainCorrectedGeo :
                terrainHeight = geoBlobArrObj.height[:,:]

            # Get scan_mode to find any bad scans

            scanMode = geoBlobArrObj.scan_mode[:]
            badScanIdx = np.where(scanMode==254)[0]
            LOG.debug("Bad Scans: %r" % (badScanIdx))

            # Detemine the min, max and range of the latitude and longitude, 
            # taking care to exclude any fill values.

            if longFormGeoNames :
                latitude = getattr(geoBlobArrObj,'latitude').astype('float')
                longitude = getattr(geoBlobArrObj,'longitude').astype('float')
            else :
                latitude = getattr(geoBlobArrObj,'lat').astype('float')
                longitude = getattr(geoBlobArrObj,'lon').astype('float')
            
            latitude = ma.masked_less(latitude,-800.)
            latMin,latMax = np.min(latitude),np.max(latitude)
            latRange = latMax-latMin

            longitude = ma.masked_less(longitude,-800.)
            lonMin,lonMax = np.min(longitude),np.max(longitude)
            lonRange = lonMax-lonMin

            LOG.debug("min,max,range of latitide: %f %f %f" % (latMin,latMax,latRange))
            LOG.debug("min,max,range of longitude: %f %f %f" % (lonMin,lonMax,lonRange))

            # Determine the latitude and longitude fill masks, so we can restore the 
            # fill values after we have scaled...

            latMask = latitude.mask
            lonMask = longitude.mask

            # Check if the geolocation is in radians, convert to degrees
            # FIXME : This information is likely conveyed by whether the 
            # FIXME :     geolocation short-name is *-GEO-TC (degrees) or
            # FIXME :     *-RGEO_TC (radians).
            if (lonRange < 2.*np.pi) :
                LOG.debug("Geolocation is in radians, convert to degrees...")
                latitude = np.degrees(latitude)
                longitude = np.degrees(longitude)
            
                latMin,latMax = np.min(latitude),np.max(latitude)
                latRange = latMax-latMin
                lonMin,lonMax = np.min(longitude),np.max(longitude)
                lonRange = lonMax-lonMin

                LOG.debug("New min,max,range of latitide: %f %f %f" % (latMin,latMax,latRange))
                LOG.debug("New min,max,range of longitude: %f %f %f" % (lonMin,lonMax,lonRange))

            # Restore fill values to masked pixels in geolocation

            geoFillValue = trimObj.sdrTypeFill['VDNE_FLOAT64_FILL'][latitude.dtype.name]
            latitude = ma.array(latitude,mask=latMask,fill_value=geoFillValue)
            latitude = latitude.filled()

            geoFillValue = trimObj.sdrTypeFill['VDNE_FLOAT64_FILL'][longitude.dtype.name]
            longitude = ma.array(longitude,mask=lonMask,fill_value=geoFillValue)
            longitude = longitude.filled()

            # Parse the geolocation asc file to get struct information which will be 
            # written to the ancillary asc files

            geoAscFileName = path.join(inDir,URID+".asc")
            LOG.debug("\nOpening %s..." % (geoAscFileName))

            geoAscFile = open(geoAscFileName,'rt')

            #RangeDateTimeStr =  _getAscLine(geoAscFile,"RangeDateTime")
            RangeDateTimeStr =  _getAscLine(geoAscFile,"ObservedDateTime")
            RangeDateTimeStr =  string.replace(RangeDateTimeStr,"ObservedDateTime","RangeDateTime")
            GRingLatitudeStr =  _getAscStructs(geoAscFile,"GRingLatitude",12)
            GRingLongitudeStr =  _getAscStructs(geoAscFile,"GRingLongitude",12)

            geoAscFile.close()

            # Loop through the required NAAPS datasets and create the blobs.
            # FIXME : Handle pathological geolocation cases

            firstGranule = True

            for dSet in masksCollShortNames :
            
                LOG.debug("\n>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n")
                LOG.debug("Processing dataset %s for %s" % (NAAPS_shortNameToBlobName[dSet],dSet))

                # FIXME : Account for dateline and pole crossings...

                # Massage the NAAPS data array a bit...
                NAAPS_anc = np.array(NAAPS_globalGridData[dSet])[::-1,:]
                NAAPS_anc = np.roll(NAAPS_anc,360)

                if (firstGranule) :

                    LOG.debug("\nGranulating %s ..." % (dSet))
                    LOG.debug("latitide,longitude shapes: %s, %s" %(str(latitude.shape) , str(longitude.shape)))
                    LOG.debug("NAAPS_anc.shape = %s" % (str(NAAPS_anc.shape)))
                    LOG.debug("gridLat.shape = %s" % (str(gridLat.shape)))
                    LOG.debug("gridLon.shape = %s" % (str(gridLon.shape)))

                    LOG.debug("min of NAAPS_anc  = "%(np.min(NAAPS_anc)))
                    LOG.debug("max of NAAPS_anc  = "%(np.max(NAAPS_anc)))

                    data,dataIdx = _grid2Gran(np.ravel(latitude),
                                              np.ravel(longitude),
                                              NAAPS_anc.astype(np.float64),
                                              gridLat,
                                              gridLon)

                    data = data.reshape(latitude.shape)
                    dataIdx = dataIdx.reshape(latitude.shape)
                    firstGranule = False
                    LOG.debug("Shape of first granulated %s data is %s" % (dSet,np.shape(data)))
                    LOG.debug("Shape of first granulated %s dataIdx is %s" % (dSet,np.shape(dataIdx)))

                else :

                    LOG.debug("Granulating %s using existing data indices." % (dSet))
                    NAAPS_anc = np.ravel(NAAPS_anc)
                    data = np.ravel(NAAPS_anc)[np.ravel(dataIdx)]
                    data = data.reshape(latitude.shape)
                    LOG.debug("Shape of subsequent granulated %s is %s" % (dSet,np.shape(data)))

                # Fill the required pixel trim rows in the granulated NAAPS data with 
                # the ONBOARD_PT_FILL value for the correct data type

                fillValue = trimObj.sdrTypeFill['ONBOARD_PT_FILL'][data.dtype.name]        
                data = ma.array(data,mask=modTrimMask,fill_value=fillValue)
                data = data.filled()

                # Create new NAAPS ancillary blob, and copy granulated data to it

                endian = ancEndian
                xmlName = path.join(ADL_HOME,'xml/VIIRS',NAAPS_shortNameToXmlName[dSet])

                # Create a new URID to be used in making the asc filenames

                URID_dict = _getURID()

                URID = URID_dict['URID']
                creationDate_nousecStr = URID_dict['creationDate_nousecStr']
                creationDateStr = URID_dict['creationDateStr']

                # Create a new directory in the input directory for the new ancillary
                # asc and blob files

                blobDir = inDir

                ascFileName = path.join(blobDir,URID+'.asc')
                blobName = path.join(blobDir,URID+'.'+dSet)

                LOG.debug("ascFileName : %s" % (ascFileName))
                LOG.debug("blobName : %s" % (blobName))

                # Create a new ancillary blob, and copy the data to it.
                newNAAPSblobObj = adl_blob.create(xmlName, blobName, endian=endian, overwrite=True)
                newNAAPSblobArrObj = newNAAPSblobObj.as_arrays()

                blobData = getattr(newNAAPSblobArrObj,'faot550')
                blobData[:,:] = data[:,:]

                # Make a new NAAPS asc file from the template, and substitute for the various tags

                ascTemplateFileName = path.join(ADL_ASC_TEMPLATES,"VIIRS-ANC_Template.asc")

                LOG.debug("Creating new asc file\n%s\nfrom template\n%s" % (ascFileName,ascTemplateFileName))
                
                ANC_fileList = gridBlobFiles
                for idx in range(len(ANC_fileList)) :
                    ANC_fileList[idx] = path.basename(ANC_fileList[idx])
                ANC_fileList.sort()
                ancGroupRecipe = '    ("N_Anc_Filename" STRING EQ "%s")'
                ancFileStr = "%s" % ("\n").join([ancGroupRecipe % (str(files)) for files in ANC_fileList])

                LOG.debug("RangeDateTimeStr = %s\n" % (RangeDateTimeStr))
                LOG.debug("GRingLatitudeStr = \n%s\n" % (GRingLatitudeStr))
                LOG.debug("GRingLongitudeStr = \n%s\n" % (GRingLongitudeStr))


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
                   line = line.replace("CSPP_ANC_BLOB_FULLPATH",blobName)
                   line = line.replace("CSPP_ANC_COLLECTION_SHORT_NAME",dSet)
                   line = line.replace("CSPP_GRANULE_ID",N_Granule_ID)
                   line = line.replace("CSPP_CREATIONDATETIME",creationDateStr)
                   line = line.replace("  CSPP_RANGE_DATE_TIME",RangeDateTimeStr)
                   line = line.replace("  CSPP_GRINGLATITUDE",GRingLatitudeStr)
                   line = line.replace("  CSPP_GRINGLONGITUDE",GRingLongitudeStr)
                   line = line.replace("    CSPP_ANC_SOURCE_FILES",ancFileStr)
                   ascFile.write(line) 

                ascFile.close()
                ascTemplateFile.close()

