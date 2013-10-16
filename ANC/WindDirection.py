#!/usr/bin/env python
# encoding: utf-8
"""
WindDirection.py

 * DESCRIPTION:  Class to granulate the ViirsAncWindDirection data product 

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

from NCEPtoBlob import NCEPclass

# skim and convert routines for reading .asc metadata fields of interest
import adl_blob
import adl_asc
from adl_asc import skim_dir, contiguous_granule_groups, granule_groups_contain, effective_anc_contains,_eliminate_duplicates,_is_contiguous, RDR_REQUIRED_KEYS, POLARWANDER_REQUIRED_KEYS
from adl_common import ADL_HOME, CSPP_RT_HOME, CSPP_RT_ANC_PATH, CSPP_RT_ANC_CACHE_DIR, COMMON_LOG_CHECK_TABLE

# every module should have a LOG object
try :
    sourcename= file_Id.split(" ")
    LOG = logging.getLogger(sourcename[1])
except :
    LOG = logging.getLogger('WindDirection')

from Utils import getURID, getAscLine, getAscStructs, findDatelineCrossings, shipOutToFile

class WindDirection() :

    def __init__(self,inDir=None, sdrEndian=None, ancEndian=None):
        self.collectionShortName = 'VIIRS-ANC-Wind-Direction-Mod-Gran'
        self.xmlName = 'VIIRS_ANC_WIND_DIRECTION_MOD_GRAN.xml'
        self.blobDatasetName = ['uComponentOfWind', 'vComponentOfWind']
        self.dataType = 'float32'
        self.sourceType = 'NCEP_ANC_Int'
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

        dates = []
        ncepBlobFiles = []
        for gridBlobStruct in ancBlob:
            timeObj = gridBlobStruct[0]
            ncepBlobFile = gridBlobStruct[1]
            LOG.debug("VIIRS-ANC-Wind-Direction-Mod-Gran %s --> %s" % \
                    (ncepBlobFile,timeObj.strftime("%Y-%m-%d %H:%M:%S:%f")))
            dates.append(timeObj)
            ncepBlobFiles.append(ncepBlobFile)

        self.date_0 = dates[0]
        self.date_1 = dates[1]

        LOG.debug("Minimum NCEP date is: %s" %(self.date_0.strftime("%Y-%m-%d %H:%M:%S:%f")))
        LOG.debug("Maximum NCEP date is: %s" %(self.date_1.strftime("%Y-%m-%d %H:%M:%S:%f")))

        ncepBlobFile_0 = ncepBlobFiles[0]
        ncepBlobFile_1 = ncepBlobFiles[1]

        dsetName = self.blobDatasetName[0]
        self.uWind_0 = getattr(ncepBlobFile_0,dsetName).astype(self.dataType)
        self.uWind_1 = getattr(ncepBlobFile_1,dsetName).astype(self.dataType)

        dsetName = self.blobDatasetName[1]
        self.vWind_0 = getattr(ncepBlobFile_0,dsetName).astype(self.dataType)
        self.vWind_1 = getattr(ncepBlobFile_1,dsetName).astype(self.dataType)


    def setGeolocationInfo(self,dicts):
        '''
        Populate this class instance with the geolocation data for a single granule
        '''
        # Set some environment variables and paths
        ANC_SCRIPTS_PATH = path.join(CSPP_RT_HOME,'viirs')
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

        timeDelta = (self.date_1 - self.date_0).total_seconds()
        LOG.debug("timeDelta is %r seconds" %(timeDelta))

        timePrime = (ObservedStartTimeObj - self.date_0).total_seconds()
        LOG.debug("timePrime is %r seconds (%f percent along time interval)" % \
                (timePrime,(timePrime/timeDelta)*100.))

        delta_uWind = self.uWind_1 - self.uWind_0
        self.uWind = (delta_uWind/timeDelta) * timePrime + self.uWind_0

        uWind_0_avg = np.average(self.uWind_0)
        uWind_1_avg = np.average(self.uWind_1)
        uWind_avg = np.average(self.uWind)
        LOG.debug("average(uWind_0) = %f" %(np.average(self.uWind_0)))
        LOG.debug("average(uWind_1) = %f" %(np.average(self.uWind_1)))
        LOG.debug("average(uWind) = %f" %(np.average(self.uWind)))

        delta_vWind = self.vWind_1 - self.vWind_0
        self.vWind = (delta_vWind/timeDelta) * timePrime + self.vWind_0

        vWind_0_avg = np.average(self.vWind_0)
        vWind_1_avg = np.average(self.vWind_1)
        vWind_avg = np.average(self.vWind)
        LOG.debug("average(vWind_0) = %f" %(np.average(self.vWind_0)))
        LOG.debug("average(vWind_1) = %f" %(np.average(self.vWind_1)))
        LOG.debug("average(vWind) = %f" %(np.average(self.vWind)))

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
            LOG.warning("Geolocation is in radians, convert to degrees...")
            latitude = np.degrees(latitude)
            longitude = np.degrees(longitude)
        
            latMin,latMax = np.min(latitude),np.max(latitude)
            latRange = latMax-latMin
            lonMin,lonMax = np.min(longitude),np.max(longitude)
            lonRange = lonMax-lonMin

            LOG.debug("New min,max,range of latitide: %f %f %f" % (latMin,latMax,latRange))
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
        LOG.debug("We have %d dateline crossings."%(num180Crossings))

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
        uWind = self.uWind[::-1,:]
        vWind = self.vWind[::-1,:]

        if self.num180Crossings != 2 :

            uWind = np.roll(uWind,360)
            vWind = np.roll(vWind,360)
            gridLon,gridLat = np.meshgrid(lons,lats)

            LOG.debug("start,end NCEP Grid Latitude values : %f,%f"%(gridLat[0,0],gridLat[-1,0]))
            LOG.debug("start,end NCEP Grid Longitude values : %f,%f"%(gridLon[0,0],gridLon[0,-1]))

        else :

            negLonIdx = np.where(lons<0)
            lons[negLonIdx] += 360.
            lons = np.roll(lons,360)
            gridLon,gridLat = np.meshgrid(lons,lats)

            longitudeNegIdx = np.where(longitude < 0.)
            longitude[longitudeNegIdx] += 360.

            LOG.debug("start,end NCEP Grid Latitude values : %f,%f"%(gridLat[0,0],gridLat[-1,0]))
            LOG.debug("start,end NCEP Grid Longitude values : %f,%f"%(gridLon[0,0],gridLon[0,-1]))


        LOG.debug("min of uWind  = %r"%(np.min(uWind)))
        LOG.debug("max of uWind  = %r"%(np.max(uWind)))
        LOG.debug("min of vWind  = %r"%(np.min(vWind)))
        LOG.debug("max of vWind  = %r"%(np.max(vWind)))

        # Separately granulate the U and V wind magnitudes
        t1 = time()
        uData,dataIdx = self._grid2Gran_bilinearInterp(np.ravel(latitude),
                                  np.ravel(longitude),
                                  uWind.astype(np.float64),
                                  gridLat,
                                  gridLon)

        vData,dataIdx = self._grid2Gran_bilinearInterp(np.ravel(latitude),
                                  np.ravel(longitude),
                                  vWind.astype(np.float64),
                                  gridLat,
                                  gridLon)
        t2 = time()
        elapsedTime = t2-t1
        LOG.info("Granulation took %f seconds for %d points" % (elapsedTime,latitude.size))

        uData = uData.reshape(latitude.shape)
        vData = vData.reshape(latitude.shape)
        dataIdx = dataIdx.reshape(latitude.shape)

        LOG.debug("Shape of granulated %s uData is %s" % (self.collectionShortName,np.shape(uData)))
        LOG.debug("Shape of granulated %s vData is %s" % (self.collectionShortName,np.shape(vData)))
        LOG.debug("Shape of granulated %s dataIdx is %s" % (self.collectionShortName,np.shape(dataIdx)))

        # Combine the granulated U and V wind magnitudes into a direction...
        windDirection_rad  = np.pi/2. - np.arctan2(vData, uData)
        negWindDirIdx = np.where(windDirection_rad < 0.)
        windDirection_rad[negWindDirIdx] = windDirection_rad[negWindDirIdx] + 2.*np.pi
        data = np.degrees(windDirection_rad)

        # Moderate resolution trim table arrays. These are 
        # bool arrays, and the trim pixels are set to True.
        modTrimMask = self.trimObj.createModTrimArray(nscans=48,trimType=bool)

        # Fill the required pixel trim rows in the granulated NCEP data with 
        # the ONBOARD_PT_FILL value for the correct data type

        fillValue = self.trimObj.sdrTypeFill['ONBOARD_PT_FILL'][self.dataType]        
        data = ma.array(data,mask=modTrimMask,fill_value=fillValue)
        self.data = data.filled()


    def shipOutToFile(self):
        ''' Pass the current class instance to this Utils method to generate 
            a blob/asc file pair from the input ancillary data object.'''

        return shipOutToFile(self)
