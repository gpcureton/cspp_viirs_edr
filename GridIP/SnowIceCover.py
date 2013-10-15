#!/usr/bin/env python
# encoding: utf-8
"""
SnowIceCover.py

 * DESCRIPTION:  Class to granulate the SnowIceCover data product 

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

import numpy as np
from numpy import ma
import copy
from bisect import bisect_left,bisect_right

import ctypes
from numpy.ctypeslib import ndpointer

from HDF4File import HDF4File

import ViirsData

# skim and convert routines for reading .asc metadata fields of interest
import adl_blob
import adl_asc
from adl_asc import skim_dir, contiguous_granule_groups, granule_groups_contain, effective_anc_contains,_eliminate_duplicates,_is_contiguous, RDR_REQUIRED_KEYS, POLARWANDER_REQUIRED_KEYS
from adl_common import ADL_HOME, CSPP_RT_HOME, CSPP_RT_ANC_PATH, CSPP_RT_ANC_CACHE_DIR, COMMON_LOG_CHECK_TABLE,env

# every module should have a LOG object
try :
    sourcename= file_Id.split(" ")
    LOG = logging.getLogger(sourcename[1])
except :
    LOG = logging.getLogger('SnowIceCover')

from Utils import getURID, getAscLine, getAscStructs, findDatelineCrossings, shipOutToFile, plotArr
from Utils import index, find_lt, find_le, find_gt, find_ge

class SnowIceCover() :

    def __init__(self,inDir=None, sdrEndian=None, ancEndian=None):
        self.collectionShortName = 'VIIRS-GridIP-VIIRS-Snow-Ice-Cover-Mod-Gran'
        self.xmlName = 'VIIRS_GRIDIP_VIIRS_SNOW_ICE_COVER_MOD_GRAN.xml'
        self.blobDatasetName = 'snowIceCover'
        self.dataType = 'uint8'
        self.sourceType = 'NISE'
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

        self.RangeDateTimeStr =  getAscLine(geoAscFile,"ObservedDateTime")
        self.RangeDateTimeStr =  string.replace(self.RangeDateTimeStr,"ObservedDateTime","RangeDateTime")

        self.GRingLatitudeStr =  getAscStructs(geoAscFile,"GRingLatitude",12)
        self.GRingLongitudeStr = getAscStructs(geoAscFile,"GRingLongitude",12)

        self.North_Bounding_Coordinate_Str = getAscLine(geoAscFile,"North_Bounding_Coordinate")
        self.South_Bounding_Coordinate_Str = getAscLine(geoAscFile,"South_Bounding_Coordinate")
        self.East_Bounding_Coordinate_Str  = getAscLine(geoAscFile,"East_Bounding_Coordinate")
        self.West_Bounding_Coordinate_Str  = getAscLine(geoAscFile,"West_Bounding_Coordinate")

        geoAscFile.close()


    def __retrieve_NISE_files(self):
        ''' Download the NISE Snow/Ice files which cover the dates of the geolocation files.'''

        ANC_SCRIPTS_PATH = path.join(CSPP_RT_HOME,'viirs')

        niseFiles = []

        geoDict = self.geoDict

        timeObj = geoDict['ObservedStartTime']
        dateStamp = timeObj.strftime("%Y%j")

        try :
            LOG.info('Retrieving NISE files for %s ...' % (dateStamp))
            cmdStr = '%s/get_anc_cspp_nise.csh %s' % (ANC_SCRIPTS_PATH,dateStamp)
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
                    niseFiles.append(lines)

            # TODO : On error, jump to a cleanup routine
            if not (procRetVal == 0) :
                LOG.error('Retrieval of NISE files failed for %r.' % (dateStamp))
                #sys.exit(procRetVal)

        except Exception, err:
            LOG.exception( "%s" % (str(err)))

        if (niseFiles == []) :
            LOG.error('Failed to find or retrieve any NISE files for date %r, aborting...'%(timeObj))
            sys.exit(1)

        # Uniqify the list of NISE files
        niseFiles = list(set(niseFiles))
        niseFiles.sort()

        for niseFile in niseFiles :
            LOG.info('Retrieved NISE file: %r' % (niseFile))

        self.NISE_fileNames = niseFiles


    def subset(self):
        '''Subsets the NISE Snow and Ice files.'''

        # Retrieve any required NISE files.
        self.__retrieve_NISE_files()

        # If more than one NISE file, take the first
        NISE_fileName = self.NISE_fileNames[0]

        LOG.debug("Opening the NISE file %s" % (NISE_fileName))
        try :
            fileObj = HDF4File(NISE_fileName)
        except Exception, err :
            LOG.exception("%s"%(err))
            LOG.exception("Problem opening NISE file (%s), aborting."%(NISE_fileName))
            sys.exit(1)

        try :

            northDsetName = "Northern Hemisphere/Data Fields/Extent"
            southDsetName = "Southern Hemisphere/Data Fields/Extent"

            LOG.debug("Retrieving NISE HDF4 path '%s'" % (northDsetName))
            self.nHemi = fileObj.get_dataset(northDsetName)

            LOG.debug("Retrieving NISE HDF4 path '%s'" % (southDsetName))
            self.sHemi = fileObj.get_dataset(southDsetName)

        except Exception, err :

            LOG.debug("EXCEPTION: %s" % (err))
            sys.exit(1)


    def granulate(self,GridIP_objects):
        '''Granulates the NISE Snow and Ice files.'''

        #import GridIP

        # Get the geolocation information
        geoDict = self.geoDict
        inDir = self.inDir

        LOG.info("Granulating %s ..." % (self.collectionShortName))

        LSM = GridIP_objects['VIIRS-GridIP-VIIRS-Lwm-Mod-Gran'].data

        latitude = self.latitude
        longitude = self.longitude

        # If we have a dateline crossing, remove the longitude discontinuity
        # by adding 360 degrees to the negative longitudes.
        if self.num180Crossings == 2 :
            longitudeNegIdx = np.where(longitude < 0.)
            longitude[longitudeNegIdx] += 360.

        latMin = self.latMin
        latMax = self.latMax
        lonMin = self.lonMin
        lonMax = self.lonMax

        LOG.debug("latMin = %d" % (latMin))
        LOG.debug("latMax = %d" % (latMax))
        LOG.debug("lonMin = %d" % (lonMin))
        LOG.debug("lonMax = %d" % (lonMax))

        nHemi = self.nHemi.astype('uint8')
        sHemi = self.sHemi.astype('uint8')

        xsize, ysize = nHemi.shape    

        nHemi = nHemi.flatten()
        sHemi = sHemi.flatten()
        
        pi = np.pi

        rg = 6371.228 / 25.067525

        r0 = (ysize-1) / 2.0
        s0 = (xsize-1) / 2.0

        gridrows, gridcols = latitude.shape

        posLatMask = ma.masked_less(latitude,0.).mask
        negLatMask = ma.masked_greater_equal(latitude,0.).mask

        if posLatMask.shape == () :
            posLatMask = np.zeros(latitude.shape,dtype='bool')

        if negLatMask.shape == () :
            negLatMask = np.zeros(latitude.shape,dtype='bool')

        # For positive latitudes...
        
        phi = np.radians(ma.array(latitude,mask=posLatMask,fill_value=0))
        lam = np.radians(ma.array(longitude,mask=posLatMask,fill_value=0))

        rho =  2. * rg * np.sin((pi / 4.0) - (phi / 2.0))

        r = r0 + rho * np.sin(lam)
        s = s0 + rho * np.cos(lam)          

        i_pos = round_(r, 0).astype('int') - 1
        j_pos = round_(s, 0).astype('int') - 1


        # For negative latitudes...

        phi = np.radians(ma.array(latitude,mask=negLatMask,fill_value=0))
        lam = np.radians(ma.array(longitude,mask=negLatMask,fill_value=0))

        rho = 2. * rg * np.cos((pi / 4.0) - (phi / 2.0))
        r = r0 + rho * np.sin(lam)
        s = s0 - rho * np.cos(lam)

        i_neg = round_(r, 0).astype('int') - 1
        j_neg = round_(s, 0).astype('int') - 1

        ###
        # Combine the +ve and -ve latitudes
        ###

        posIdx = np.where(posLatMask==False)
        negIdx = np.where(negLatMask==False)

        # convert i_pos and j_pos to an index into the raveled gridded data.
        k_pos = (j_pos * xsize) + i_pos
        k_neg = (j_neg * xsize) + i_neg

        # convert i_neg and j_neg to an index into the raveled sHemi..

        ###
        # Obtain the raw NISE snow/ice values
        ###
        nise_val = 999 * np.ones(latitude.shape,dtype='int')

        nise_val[posIdx] = nHemi[k_pos[posIdx]]
        nise_val[negIdx] = sHemi[k_neg[negIdx]]

        # From ADL/include/ProEdrViirsCMIPGbl.h ...

        CM_LAND = 1
        CM_COASTAL = 5
        CM_SEA_WATER = 3
        CM_IN_WATER = 2
        CM_SNOW = 1
        CM_NO_SNOW = 0

        # Permanent Ice
        permIceMask = ma.masked_equal(nise_val,101).mask

        # Dry Snow
        drySnowMask = ma.masked_equal(nise_val,103).mask

        # Wet Snow
        wetSnowMask = ma.masked_equal(nise_val,104).mask

        # Land and Coastal
        LSM_Land_Coastal = (ma.masked_equal(LSM,CM_LAND) * ma.masked_equal(LSM,CM_COASTAL)).mask

        # Sea Water and Inland Water
        LSM_SeaWater_InlandWater = (ma.masked_equal(LSM,CM_SEA_WATER) * ma.masked_equal(LSM,CM_IN_WATER)).mask

        # Transient sea ice > 25% concentration, AND (Land OR Coastal)
        transSeaIceCoastalMask = ma.masked_inside(nise_val,(25+1),(101-1)).mask
        
        # Transient sea ice > 25% concentration, AND Permanent Ice AND (Sea OR Inland Water)
        transSeaIceOceanMask = ma.masked_inside(nise_val,(25+1),(102-1)).mask

        snowMask = permIceMask + drySnowMask + wetSnowMask \
                + (transSeaIceCoastalMask * LSM_Land_Coastal) \
                + (transSeaIceOceanMask * LSM_SeaWater_InlandWater)

        snowIceMask = np.ones(nise_val.shape,dtype='int') * CM_NO_SNOW
        snowIceMask = ma.array(snowIceMask,mask=snowMask,fill_value=CM_SNOW)
        snowIceMask = snowIceMask.filled()
        snowIceMask = snowIceMask.astype(self.dataType)

        # Explicitly restore geolocation fill to the granulated data...
        fillMask = ma.masked_less(self.latitude,-800.).mask
        fillValue = self.trimObj.sdrTypeFill['MISS_FILL'][self.dataType]        
        data = ma.array(snowIceMask,mask=fillMask,fill_value=fillValue)
        data = data.filled()

        # Moderate resolution trim table arrays. These are 
        # bool arrays, and the trim pixels are set to True.
        modTrimMask = self.trimObj.createModTrimArray(nscans=48,trimType=bool)

        # Fill the required pixel trim rows in the granulated NISE data with 
        # the ONBOARD_PT_FILL value for the correct data type

        fillValue = self.trimObj.sdrTypeFill['ONBOARD_PT_FILL'][self.dataType]        
        data = ma.array(data,mask=modTrimMask,fill_value=fillValue)
        self.data = data.filled()


    def shipOutToFile(self):
        ''' Pass the current class instance to this Utils method to generate 
            a blob/asc file pair from the input ancillary data object.'''

        return shipOutToFile(self)
