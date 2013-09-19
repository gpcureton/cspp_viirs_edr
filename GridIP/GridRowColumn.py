#!/usr/bin/env python
# encoding: utf-8
"""
GridRowColumn.py

 * DESCRIPTION:  Class to prepare the Grid-Row-Column data product 

Created by Geoff Cureton on 2013-03-14.
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

import tables as pytables
from tables import exceptions as pyEx

import ViirsData

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
    LOG = logging.getLogger('GridRowColumn')

from Utils import getURID, getAscLine, getAscStructs, findDatelineCrossings, shipOutToFile
from Utils import index, find_lt, find_le, find_gt, find_ge

class GridRowColumn() :

    def __init__(self,inDir=None, sdrEndian=None, ancEndian=None):
        self.collectionShortName = ['VIIRS-MOD-GRC-TC', 'VIIRS-MOD-GRC']
        self.xmlName = {'VIIRS-MOD-GRC-TC':'VIIRS_MOD_GRC_TC.xml', 'VIIRS-MOD-GRC':'VIIRS_MOD_GRC.xml'}
        self.blobDatasetName = 'gridRowCol'
        self.dataType = 'int64'
        self.sourceType = ''
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


    def setGeolocationInfo(self,dicts):
        '''
        Populate this class instance with the geolocation data for a single granule
        '''
        # Set some environment variables and paths
        CSPP_RT_HOME = os.getenv('CSPP_RT_HOME')
        ANC_SCRIPTS_PATH = path.join(CSPP_RT_HOME,'viirs')
        CSPP_RT_ANC_CACHE_DIR = os.getenv('CSPP_RT_ANC_CACHE_DIR')
    
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
        LOG.info("We have %d dateline crossings."%(num180Crossings))
        
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

        self.ObservedDateTimeStr =  getAscLine(geoAscFile,"ObservedDateTime")
        self.RangeDateTimeStr =  self.ObservedDateTimeStr
        self.RangeDateTimeStr =  string.replace(self.RangeDateTimeStr,"ObservedDateTime","RangeDateTime")

        self.GRingLatitudeStr =  getAscStructs(geoAscFile,"GRingLatitude",12)
        self.GRingLongitudeStr = getAscStructs(geoAscFile,"GRingLongitude",12)

        self.BeginningOrbitNumber  = getAscLine(geoAscFile,"BeginningOrbitNumber")
        self.N_Nadir_Latitude_Max  = getAscLine(geoAscFile,"N_Nadir_Latitude_Max")
        self.N_Nadir_Latitude_Min  = getAscLine(geoAscFile,"N_Nadir_Latitude_Min")
        self.N_Nadir_Longitude_Max = getAscLine(geoAscFile,"N_Nadir_Longitude_Max")
        self.N_Nadir_Longitude_Min = getAscLine(geoAscFile,"N_Nadir_Longitude_Min")

        self.North_Bounding_Coordinate_Str = getAscLine(geoAscFile,"North_Bounding_Coordinate")
        self.South_Bounding_Coordinate_Str = getAscLine(geoAscFile,"South_Bounding_Coordinate")
        self.East_Bounding_Coordinate_Str  = getAscLine(geoAscFile,"East_Bounding_Coordinate")
        self.West_Bounding_Coordinate_Str  = getAscLine(geoAscFile,"West_Bounding_Coordinate")

        self.N_Day_Night_Flag  = getAscLine(geoAscFile,"N_Day_Night_Flag")
        self.Ascending_Descending_Indicator  = getAscLine(geoAscFile,"Ascending/Descending_Indicator")

        geoAscFile.close()


    def subset(self):
        pass


    def granulate(self,GridIP_objects):
        ''' This method has been retasked to set the collection short name and 
            xml name to the correct value rather than a list.
        '''
        #self.collectionShortName = 
        #self.xmlName = 


    def shipOutToFile(self):
        ''' Pass the current class instance to this Utils method to generate 
            a blob/asc file pair from the input ancillary data object.
        '''

        # Set some environment variables and paths
        CSPP_RT_HOME = os.getenv('CSPP_RT_HOME')
        ANC_SCRIPTS_PATH = path.join(CSPP_RT_HOME,'viirs')
        CSPP_RT_ANC_CACHE_DIR = os.getenv('CSPP_RT_ANC_CACHE_DIR')
        ADL_ASC_TEMPLATES = path.join(ANC_SCRIPTS_PATH,'asc_templates')

        # Create new GridIP ancillary blob, and copy granulated data to it

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

        newGridIPblobObj = adl_blob.create(xmlName, blobName, endian=endian, overwrite=True)

        # Make a new GridIP asc file from the template, and substitute for the various tags

        templateName = "%s_Template.asc" % (self.collectionShortName)
        ascTemplateFileName = path.join(ADL_ASC_TEMPLATES,templateName)

        LOG.info("ascTemplateFileName = %s" % (ascTemplateFileName))

        LOG.info("Creating new asc file %s from template %s" % (ascFileName,ascTemplateFileName))
        
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
           line = line.replace("CSPP_BEGINNING_ORBIT_NUMBER",self.BeginningOrbitNumber)
           line = line.replace("CSPP_CREATIONDATETIME",creationDateStr)
           line = line.replace("  CSPP_OBSERVED_DATE_TIME",self.ObservedDateTimeStr)
           line = line.replace("  CSPP_RANGE_DATE_TIME",self.RangeDateTimeStr)
           line = line.replace("  CSPP_GRINGLATITUDE",self.GRingLatitudeStr)
           line = line.replace("  CSPP_GRINGLONGITUDE",self.GRingLongitudeStr)

           line = line.replace("  CSPP_NADIR_LATITUDE_MAX",self.N_Nadir_Latitude_Max)
           line = line.replace("  CSPP_NADIR_LATITUDE_MIN",self.N_Nadir_Latitude_Min)
           line = line.replace("  CSPP_NADIR_LONGITUDE_MAX",self.N_Nadir_Longitude_Max)
           line = line.replace("  CSPP_NADIR_LONGITUDE_MIN",self.N_Nadir_Longitude_Min)

           line = line.replace("  CSPP_NORTH_BOUNDING_COORD",self.North_Bounding_Coordinate_Str)
           line = line.replace("  CSPP_SOUTH_BOUNDING_COORD",self.South_Bounding_Coordinate_Str)
           line = line.replace("  CSPP_EAST_BOUNDING_COORD",self.East_Bounding_Coordinate_Str)
           line = line.replace("  CSPP_WEST_BOUNDING_COORD",self.West_Bounding_Coordinate_Str)

           line = line.replace("  CSPP_DAY_NIGHT_FLAG",self.N_Day_Night_Flag)
           line = line.replace("  CSPP_ASCENDING_DESCENDING_INDICATOR",self.Ascending_Descending_Indicator)
           ascFile.write(line) 

        ascFile.close()
        ascTemplateFile.close()

        return URID
