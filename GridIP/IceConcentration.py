#!/usr/bin/env python
# encoding: utf-8
"""
IceConcentration.py

 * DESCRIPTION:  Class to granulate the Ice Concentration data product 

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
import re
import uuid
import shlex, subprocess
from subprocess import CalledProcessError, call
from glob import glob
from time import time
from datetime import datetime,timedelta

from scipy import round_
from scipy import interpolate

import numpy as np
from numpy import ma
import copy
from bisect import bisect_left,bisect_right

import ctypes
from numpy.ctypeslib import ndpointer

import pygrib

import ViirsData

# skim and convert routines for reading .asc metadata fields of interest
import adl_blob
import adl_asc
from adl_asc import skim_dir, contiguous_granule_groups, granule_groups_contain, effective_anc_contains,_eliminate_duplicates,_is_contiguous, RDR_REQUIRED_KEYS, POLARWANDER_REQUIRED_KEYS
from adl_common import ADL_HOME, CSPP_RT_HOME, CSPP_RT_ANC_PATH, CSPP_RT_ANC_CACHE_DIR, COMMON_LOG_CHECK_TABLE, env

# every module should have a LOG object
try :
    sourcename= file_Id.split(" ")
    LOG = logging.getLogger(sourcename[1])
except :
    LOG = logging.getLogger('IceConcentration')

from Utils import getURID, getAscLine, getAscStructs, findDatelineCrossings, shipOutToFile
from Utils import index, find_lt, find_le, find_gt, find_ge
#from Utils import plotArr

class IceConcentration() :

    def __init__(self,inDir=None, sdrEndian=None, ancEndian=None):
        self.collectionShortName = 'VIIRS-I-Conc-IP'
        self.xmlName = 'VIIRS_I_CONC_IP.xml'
        self.blobDatasetName = ['iceFraction','iceConcWeights']
        self.dataType = 'float32'
        self.sourceType = 'MMAB'
        self.sourceList = ['']
        self.trimObj = ViirsData.ViirsTrimTable()

        self.GridIP_collectionShortNames = [
                            'VIIRS-GridIP-VIIRS-Lwm-Mod-Gran'
                          ]

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
        geoAscFileName = dicts['_filename']
        geoBlobFileName = string.replace(geoAscFileName,'asc',geo_Collection_ShortName)

        LOG.debug("\n###########################")
        LOG.debug("  Geolocation Information  ")
        LOG.debug("###########################")
        LOG.debug("N_Granule_ID : %r" % (N_Granule_ID))
        LOG.debug("ObservedStartTime : %s" % (ObservedStartTimeObj.__str__()))
        LOG.debug("N_Collection_Short_Name : %s" %(geo_Collection_ShortName))
        LOG.debug("URID : %r" % (URID))
        LOG.debug("geoAscFileName : %r" % (geoAscFileName))
        LOG.debug("geoBlobFileName : %r" % (geoAscFileName))
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
            LOG.debug("We are using for %s: %s,%s" %(geo_Collection_ShortName,geoXmlFile,geoBlobFileName))

        # Open the geolocation blob and get the latitude and longitude

        endian = self.sdrEndian

        #geoBlobObj = adl_blob.map(geoXmlFile,geoFiles[0], endian=endian)
        geoBlobObj = adl_blob.map(geoXmlFile,geoBlobFileName, endian=endian)
        geoBlobArrObj = geoBlobObj.as_arrays()

        # Get scan_mode to find any bad scans

        #scanMode = geoBlobArrObj.scan_mode[:]
        scanMode = getattr(geoBlobArrObj,'scan_mode').astype('uint8')
        LOG.debug("Scan Mode = %r" % (scanMode))

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
            LOG.info("Geolocation is in radians, convert to degrees...")
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
        latitude = latitude.filled()

        geoFillValue = self.trimObj.sdrTypeFill['VDNE_FLOAT64_FILL'][longitude.dtype.name]
        longitude = ma.array(longitude,mask=lonMask,fill_value=geoFillValue)
        longitude = longitude.filled()

        # Shift the longitudes to be between -180 and 180 degrees
        if lonMax > 180. :
            LOG.debug("\nFinal min,max,range of longitude: %f %f %f" % (lonMin,lonMax,lonRange))
            # Scale to restore -ve longitudess, not necessarily # FIXME
            dateLineIdx = np.where(longitude>180.)
            LOG.debug("dateLineIdx = %r" % (dateLineIdx))
            longitude[dateLineIdx] -= 360.
            lonMax = np.max(ma.array(longitude,mask=lonMask))
            lonMin = np.min(ma.array(longitude,mask=lonMask))
            lonRange = lonMax-lonMin
            LOG.debug("\nFinal min,max,range of longitude: %f %f %f" % (lonMin,lonMax,lonRange))

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
        self.latitude  = latitude
        self.longitude = longitude
        self.scanMode  = scanMode
        self.latCrnList  = latCrnList
        self.lonCrnList  = lonCrnList
        self.num180Crossings  = num180Crossings

        # Parse the geolocation asc file to get struct information which will be 
        # written to the ancillary asc files

        LOG.debug("geolocation asc filename : %s"%(geoAscFileName))

        LOG.debug("\nOpening %s..." % (geoAscFileName))

        geoAscFile = open(geoAscFileName,'rt')

        #########
        #self.RangeDateTimeStr =  getAscLine(geoAscFile,"ObservedDateTime")
        #self.RangeDateTimeStr =  string.replace(self.RangeDateTimeStr,"ObservedDateTime","RangeDateTime")

        #self.North_Bounding_Coordinate_Str = getAscLine(geoAscFile,"North_Bounding_Coordinate")
        #self.South_Bounding_Coordinate_Str = getAscLine(geoAscFile,"South_Bounding_Coordinate")
        #self.East_Bounding_Coordinate_Str  = getAscLine(geoAscFile,"East_Bounding_Coordinate")
        #self.West_Bounding_Coordinate_Str  = getAscLine(geoAscFile,"West_Bounding_Coordinate")

        #self.GRingLatitudeStr =  getAscStructs(geoAscFile,"GRingLatitude",12)
        #self.GRingLongitudeStr = getAscStructs(geoAscFile,"GRingLongitude",12)
        ##########

        self.ascFileDict = {}
        self.ascFileDict["ObservedDateTime"]               = getAscLine(geoAscFile,"ObservedDateTime")
        self.ascFileDict["RangeDateTime"]                  = getAscLine(geoAscFile,"RangeDateTime")
        self.ascFileDict["N_Granule_ID"]                   = getAscLine(geoAscFile,"N_Granule_ID")
        self.ascFileDict["N_Granule_Version"]              = getAscLine(geoAscFile,"N_Granule_Version")
        self.ascFileDict["BeginningOrbitNumber"]           = getAscLine(geoAscFile,"BeginningOrbitNumber")

        self.ascFileDict["Ascending/Descending_Indicator"] = getAscLine(geoAscFile,"Ascending/Descending_Indicator")
        self.ascFileDict["North_Bounding_Coordinate"]      = getAscLine(geoAscFile,"North_Bounding_Coordinate")
        self.ascFileDict["South_Bounding_Coordinate"]      = getAscLine(geoAscFile,"South_Bounding_Coordinate")
        self.ascFileDict["East_Bounding_Coordinate"]       = getAscLine(geoAscFile,"East_Bounding_Coordinate")
        self.ascFileDict["West_Bounding_Coordinate"]       = getAscLine(geoAscFile,"West_Bounding_Coordinate")

        self.ascFileDict["N_Nadir_Latitude_Max"]           = getAscLine(geoAscFile,"N_Nadir_Latitude_Max")
        self.ascFileDict["N_Nadir_Latitude_Min"]           = getAscLine(geoAscFile,"N_Nadir_Latitude_Min")
        self.ascFileDict["N_Nadir_Longitude_Max"]          = getAscLine(geoAscFile,"N_Nadir_Longitude_Max")
        self.ascFileDict["N_Nadir_Longitude_Min"]          = getAscLine(geoAscFile,"N_Nadir_Longitude_Min")
        
        self.ascFileDict["GRingLatitude"] =  getAscStructs(geoAscFile,"GRingLatitude",12)
        self.ascFileDict["GRingLongitude"] =  getAscStructs(geoAscFile,"GRingLongitude",12)

        #self.ascFileDict["N_Satellite/Local_Azimuth_Angle_Max"] =  getAscLine(geoAscFile,"N_Satellite/Local_Azimuth_Angle_Max")
        #self.ascFileDict["N_Satellite/Local_Azimuth_Angle_Min"] =  getAscLine(geoAscFile,"N_Satellite/Local_Azimuth_Angle_Min")
        #self.ascFileDict["N_Satellite/Local_Zenith_Angle_Max"] =  getAscLine(geoAscFile,"N_Satellite/Local_Zenith_Angle_Max")
        #self.ascFileDict["N_Satellite/Local_Zenith_Angle_Min"] =  getAscLine(geoAscFile,"N_Satellite/Local_Zenith_Angle_Min")
        #self.ascFileDict["N_Solar_Azimuth_Angle_Max"] =  getAscLine(geoAscFile,"N_Solar_Azimuth_Angle_Max")
        #self.ascFileDict["N_Solar_Azimuth_Angle_Min"] =  getAscLine(geoAscFile,"N_Solar_Azimuth_Angle_Min")
        #self.ascFileDict["N_Solar_Zenith_Angle_Max"] =  getAscLine(geoAscFile,"N_Solar_Zenith_Angle_Max")
        #self.ascFileDict["N_Solar_Zenith_Angle_Min"] =  getAscLine(geoAscFile,"N_Solar_Zenith_Angle_Min")

        geoAscFile.close()

        self.ascFileKeys = ["ObservedDateTime", "RangeDateTime", \
                            "North_Bounding_Coordinate", "South_Bounding_Coordinate", \
                            "East_Bounding_Coordinate", "West_Bounding_Coordinate", \
                            "N_Granule_Version", "N_Granule_ID", \
                            "BeginningOrbitNumber", \
                            "Ascending/Descending_Indicator", \
                            "N_Nadir_Latitude_Max", "N_Nadir_Latitude_Min", \
                            "N_Nadir_Longitude_Max", "N_Nadir_Longitude_Min", \
                            "GRingLatitude", "GRingLongitude"
                            ]         

        for key in self.ascFileKeys:
            self.ascFileDict[key] = ["CSPP_%s"%(key),self.ascFileDict[key]]
            LOG.debug(">>>> ASCFILE : %s --> %s , %s"%(key,self.ascFileDict[key][0],self.ascFileDict[key][1]))


    def __retrieve_MMAB_files(self):
        ''' Download the MMAB Ice Concentration files which cover the dates of the geolocation files.'''

        ANC_SCRIPTS_PATH = path.join(CSPP_RT_HOME,'viirs')

        mmabFiles = []

        geoDict = self.geoDict

        timeObj = geoDict['ObservedStartTime']
        dateStamp = timeObj.strftime("%Y%j")

        try :
            LOG.info('Retrieving MMAB files for %s ...' % (dateStamp))
            cmdStr = '%s/get_anc_cspp_icec.csh %s' % (ANC_SCRIPTS_PATH,dateStamp)
            LOG.info('\t%s' % (cmdStr))
            args = shlex.split(cmdStr)

            procRetVal = 0
            procObj = subprocess.Popen(args,env=env(JPSS_LOCAL_ANC_DIR=CSPP_RT_ANC_CACHE_DIR,CSPP_RT_HOME=CSPP_RT_HOME),bufsize=0, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            procObj.wait()
            procRetVal = procObj.returncode

            procOutput = procObj.stdout.readlines()
            
            for lines in procOutput:
                if CSPP_RT_ANC_CACHE_DIR in lines.replace("//","/") :
                    lines = string.replace(lines,'\n','')
                    mmabFiles.append(lines)

            # TODO : On error, jump to a cleanup routine
            if not (procRetVal == 0) :
                LOG.error('Retrieval of MMAB files failed for %r.' % (dateStamp))
                #sys.exit(procRetVal)

        except Exception, err:
            LOG.exception( "%s" % (str(err)))

        if (mmabFiles == []) :
            LOG.error('Failed to find or retrieve any MMAB files for date %r, aborting...'%(timeObj))
            sys.exit(1)

        # Uniqify the list of MMAB files
        mmabFiles = list(set(mmabFiles))
        mmabFiles.sort()

        for mmabFile in mmabFiles :
            LOG.info('Retrieved MMAB file: %r' % (mmabFile))

        self.MMAB_fileNames = mmabFiles


    def subset(self):
        '''Subsets/ingests the MMAB Ice Concentration files.'''

        # Retrieve any required MMAB files.
        self.__retrieve_MMAB_files()

        # If more than one MMAB file, take the first
        MMAB_fileName = self.MMAB_fileNames[0]

        LOG.info("Opening the MMAB file %s" % (MMAB_fileName))
        try :
            engObj = pygrib.open(MMAB_fileName)
        except Exception, err :
            LOG.exception("%s"%(err))
            LOG.exception("Problem opening MMAB file (%s), aborting."%(MMAB_fileName))
            sys.exit(1)

        try :
            msgNum = engObj.messages
            if (msgNum == 0):
                LOG.error("Grib messages must start at 1: Incorrect message number (%s) for %s, aborting"%\
                        (msgNum,MMAB_fileName))
                engObj.close()
                sys.exit(1)
            else :
                message = engObj.message(msgNum)
                engObj.close()
        except Exception, err :
            LOG.exception("%s"%(err))
            engObj.close()
            sys.exit(1)

        #analDate = message.analDate
        #lats,lons = message.latlons()
        #keys = message.keys()
        #level = message['level']
        #typeOfLevel = message['typeOfLevel']
        self.distinctLatitudes = message['distinctLatitudes']
        self.distinctLongitudes = message['distinctLongitudes']
        self.gridData = message['values']
        LOG.debug("Shape of gridded %s data is %s" % (self.collectionShortName,np.shape(self.gridData)))


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


    def _grid2Gran(self, dataLat, dataLon, gridData, gridLat, gridLon):
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
        grid2gran = lib.grid2gran_nearest
        grid2gran.restype = None
        grid2gran.argtypes = [
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

        LOG.debug("Calling C routine grid2gran()...")

        retVal = grid2gran(dataLat,
                           dataLon,
                           data,
                           nData,
                           gridLat,
                           gridLon,
                           gridData,
                           dataIdx,
                           gridRows,
                           gridCols)

        LOG.debug("Returning from C routine grid2gran()")


        return data,dataIdx


    def granulate(self,GridIP_objects):
        '''Granulates the MMAB Snow and Ice files.'''

        import GridIP

        # Get the geolocation information
        geoDict = self.geoDict
        inDir = self.inDir

        # Obtain the required GridIP collection shortnames for this algorithm
        collectionShortNames = []
        for shortName in self.GridIP_collectionShortNames :
            LOG.info("Adding %s to the list of required collection short names..." \
                    %(shortName))
            collectionShortNames.append(shortName)

        # Create a dict of GridIP class instances, which will handle ingest and 
        # granulation
        for shortName in collectionShortNames :
            className = GridIP.classNames[shortName]
            GridIP_objects[shortName] = getattr(GridIP,className)(inDir=inDir,sdrEndian=self.sdrEndian)
            LOG.info("GridIP_objects[%s].blobDatasetName = %r" % (shortName,GridIP_objects[shortName].blobDatasetName))

        # Loop through the required GridIP datasets and create the blobs.
        for shortName in collectionShortNames :
        
            LOG.info("Processing dataset %s for %s" % (GridIP_objects[shortName].blobDatasetName,shortName))

            # Set the geolocation information in this ancillary object for the current granule...
            GridIP_objects[shortName].setGeolocationInfo(geoDict)

            LOG.debug("min,max,range of latitude: %f %f %f" % (\
                    GridIP_objects[shortName].latMin,\
                    GridIP_objects[shortName].latMax,\
                    GridIP_objects[shortName].latRange))
            LOG.debug("min,max,range of longitude: %f %f %f" % (\
                    GridIP_objects[shortName].lonMin,\
                    GridIP_objects[shortName].lonMax,\
                    GridIP_objects[shortName].lonRange))

            LOG.debug("latitude corners: %r" % (GridIP_objects[shortName].latCrnList))
            LOG.debug("longitude corners: %r" % (GridIP_objects[shortName].lonCrnList))

            # Subset the gridded data for this ancillary object to cover the required lat/lon range.
            GridIP_objects[shortName].subset()

            # Interpolate latitude and longitude to make IMG resolution geolocation
            modRows,modCols = 768,3200
            imgRows,imgCols = 1536,6400

            x = np.linspace(0,modRows-1,modRows)
            y = np.linspace(0,modCols-1,modCols)
            xnew = np.linspace(0,modRows-1,imgRows)
            ynew = np.linspace(0,modCols-1,imgCols)

            z = GridIP_objects[shortName].latitude
            fRect = interpolate.RectBivariateSpline( x, y, z)
            latitude_img = fRect(xnew,ynew)
            GridIP_objects[shortName].latitude = latitude_img

            z = GridIP_objects[shortName].longitude
            fRect = interpolate.RectBivariateSpline( x, y, z)
            longitude_img = fRect(xnew,ynew)
            GridIP_objects[shortName].longitude = longitude_img

            # Granulate the gridded data in this ancillary object for the current granule...
            GridIP_objects[shortName].granulate(GridIP_objects)


        LOG.info("Granulating %s ..." % (self.collectionShortName))

        # Generate a land mask
        LSM = GridIP_objects['VIIRS-GridIP-VIIRS-Lwm-Mod-Gran'].data
        DEM_LAND = GridIP_objects['VIIRS-GridIP-VIIRS-Lwm-Mod-Gran'].DEM_dict['DEM_LAND']
        DEM_COASTLINE = GridIP_objects['VIIRS-GridIP-VIIRS-Lwm-Mod-Gran'].DEM_dict['DEM_COASTLINE']
        LSM_LandMask = (ma.masked_equal(LSM,DEM_LAND) * ma.masked_equal(LSM,DEM_COASTLINE)).mask

        # Massage the gridded ice concentration to make it easier to granulate
        iceConc_Land    = ma.masked_equal(self.gridData,1.57).mask
        iceConc_BadData = ma.masked_equal(self.gridData,1.66).mask
        iceConc_Weather = ma.masked_equal(self.gridData,1.77).mask
        iceConc_Coast   = ma.masked_equal(self.gridData,1.95).mask
        iceConc_NoData  = ma.masked_equal(self.gridData,2.24).mask
        iceConc_genMask  = ma.masked_greater(self.gridData,1.).mask

        #totalMask = iceConc_Land * iceConc_BadData * iceConc_Weather * iceConc_Coast * iceConc_NoData
        totalMask = iceConc_genMask

        self.gridData = ma.array(self.gridData,mask=totalMask,fill_value=0.)
        self.gridData = self.gridData.filled()


        degInc = 0.5

        #lats = np.arange(361.)*degInc - 90.
        lats = self.distinctLatitudes
        lons = np.arange(720.)*degInc - 180.
        latitude = latitude_img
        longitude = longitude_img

        gridData = self.gridData

        if self.num180Crossings != 2 :

            gridData = np.roll(gridData,360)
            gridLon,gridLat = np.meshgrid(lons,lats)

            LOG.debug("start,end MMAB Grid Latitude values : %f,%f"%(gridLat[0,0],gridLat[-1,0]))
            LOG.debug("start,end MMAB Grid Longitude values : %f,%f"%(gridLon[0,0],gridLon[0,-1]))

        else :

            negLonIdx = np.where(lons<0)
            lons[negLonIdx] += 360.
            lons = np.roll(lons,360)
            gridLon,gridLat = np.meshgrid(lons,lats)

            longitudeNegIdx = np.where(longitude < 0.)
            longitude[longitudeNegIdx] += 360.

            LOG.debug("start,end MMAB Grid Latitude values : %f,%f"%(gridLat[0,0],gridLat[-1,0]))
            LOG.debug("start,end MMAB Grid Longitude values : %f,%f"%(gridLon[0,0],gridLon[0,-1]))


        LOG.debug("min of gridData  = %r"%(np.min(gridData)))
        LOG.debug("max of gridData  = %r"%(np.max(gridData)))

        t1 = time()
        #data,dataIdx = self._grid2Gran(np.ravel(latitude),
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

        #plotArr(data,'VIIRS-I-Conc-IP_%s.png'%(geoDict['N_Granule_ID']))

        LOG.debug("Shape of granulated %s data is %s" % (self.collectionShortName,np.shape(data)))
        LOG.debug("Shape of granulated %s dataIdx is %s" % (self.collectionShortName,np.shape(dataIdx)))

        # Apply the land mask to the granulation ice concentration...
        fillValue = self.trimObj.sdrTypeFill['NA_FILL'][self.dataType]        
        data = ma.array(data,mask=LSM_LandMask,fill_value=fillValue)
        data = data.filled()

        # Construct the ice concentration weights
        weights = np.ones(data.shape,dtype=self.dataType) * 0.1
        weights = ma.array(weights,mask=LSM_LandMask,fill_value=0.)
        weights = weights.filled()

        # Moderate resolution trim table arrays. These are 
        # bool arrays, and the trim pixels are set to True.
        imgTrimMask = self.trimObj.createImgTrimArray(nscans=48,trimType=bool)

        # Fill the required pixel trim rows in the granulated MMAB data with 
        # the ONBOARD_PT_FILL value for the correct data type

        fillValue = self.trimObj.sdrTypeFill['ONBOARD_PT_FILL'][self.dataType]

        data = ma.array(data,mask=imgTrimMask,fill_value=fillValue)
        self.data = data.filled()

        weights = ma.array(weights,mask=imgTrimMask,fill_value=fillValue)
        self.weights = weights.filled()


    def __shipOutToFile(self):
        '''
        Generate a blob/asc file pair from the input ancillary data object.
        This is a custom method for the ice concentration, which has more than 
        one blob dataset.
        '''

        # Set some environment variables and paths
        ANC_SCRIPTS_PATH = path.join(CSPP_RT_HOME,'viirs')
        ADL_ASC_TEMPLATES = path.join(ANC_SCRIPTS_PATH,'asc_templates')

        # Create new GridIP ancillary blob, and copy granulated data to it

        endian = self.ancEndian
        if endian is adl_blob.LITTLE_ENDIAN :
            endianString = "LE"
        else :
            endianString = "BE"

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
        newGridIPblobObj = adl_blob.create(xmlName, blobName, endian=endian, overwrite=True)
        newGridIPblobArrObj = newGridIPblobObj.as_arrays()

        blobData = getattr(newGridIPblobArrObj,self.blobDatasetName[0])
        blobData[:,:] = self.data[:,:]
        blobData = getattr(newGridIPblobArrObj,self.blobDatasetName[1])
        blobData[:,:] = self.weights[:,:]

        # Make a new GridIP asc file from the template, and substitute for the various tags

        #ascTemplateFileName = path.join(ADL_ASC_TEMPLATES,"VIIRS-GridIP-VIIRS_Template.asc")
        ascTemplateFileName = path.join(ADL_ASC_TEMPLATES,"VIIRS-I-Conc-IP_Template.asc")

        LOG.debug("Creating new asc file\n%s\nfrom template\n%s" % (ascFileName,ascTemplateFileName))
        
        ANC_fileList = self.sourceList
        for idx in range(len(ANC_fileList)) :
            ANC_fileList[idx] = path.basename(ANC_fileList[idx])
        ANC_fileList.sort()
        ancGroupRecipe = '    ("N_Anc_Filename" STRING EQ "%s")'
        ancFileStr = "%s" % ("\n").join([ancGroupRecipe % (str(files)) for files in ANC_fileList])

        #LOG.debug("RangeDateTimeStr = %s\n" % (self.RangeDateTimeStr))
        #LOG.debug("GRingLatitudeStr = \n%s\n" % (self.GRingLatitudeStr))
        #LOG.debug("GRingLongitudeStr = \n%s\n" % (self.GRingLongitudeStr))

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
           line = line.replace("CSPP_ANC_ENDIANNESS",endianString)
           line = line.replace("CSPP_CREATIONDATETIME",creationDateStr)

           #line = line.replace("CSPP_GRANULE_ID",self.geoDict['N_Granule_ID'])
           #line = line.replace("CSPP_CREATIONDATETIME",creationDateStr)
           #line = line.replace("  CSPP_RANGE_DATE_TIME",self.RangeDateTimeStr)
           #line = line.replace("  CSPP_GRINGLATITUDE",self.GRingLatitudeStr)
           #line = line.replace("  CSPP_GRINGLONGITUDE",self.GRingLongitudeStr)
           #line = line.replace("CSPP_NORTH_BOUNDING_COORD",self.North_Bounding_Coordinate_Str)
           #line = line.replace("CSPP_SOUTH_BOUNDING_COORD",self.South_Bounding_Coordinate_Str)
           #line = line.replace("CSPP_EAST_BOUNDING_COORD",self.East_Bounding_Coordinate_Str)
           #line = line.replace("CSPP_WEST_BOUNDING_COORD",self.West_Bounding_Coordinate_Str)

           for key in self.ascFileKeys:
               line = line.replace(self.ascFileDict[key][0],self.ascFileDict[key][1])

           ascFile.write(line) 

        ascFile.close()
        ascTemplateFile.close()

        return URID


    def shipOutToFile(self):
        ''' Pass the current class instance to this Utils method to generate 
            a blob/asc file pair from the input ancillary data object.'''

        #shipOutToFile(self)
        return self.__shipOutToFile()
