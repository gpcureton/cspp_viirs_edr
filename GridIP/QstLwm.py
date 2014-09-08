#!/usr/bin/env python
# encoding: utf-8
"""
QstLwm.py

 * DESCRIPTION:  Class to granulate the QstLwm data product 

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

# skim and convert routines for reading .asc metadata fields of interest
import adl_blob2 as adl_blob
import adl_asc
from adl_asc import skim_dir, contiguous_granule_groups, granule_groups_contain, effective_anc_contains
from adl_asc import eliminate_duplicates,_is_contiguous, RDR_REQUIRED_KEYS, POLARWANDER_REQUIRED_KEYS
from adl_common import ADL_HOME, CSPP_RT_HOME

# every module should have a LOG object
try :
    sourcename= file_Id.split(" ")
    LOG = logging.getLogger(sourcename[1])
except :
    LOG = logging.getLogger('QstLwm')

from Utils import getURID, getAscLine, getAscStructs, findDatelineCrossings, shipOutToFile
from Utils import index, find_lt, find_le, find_gt, find_ge


class QstLwm() :

    def __init__(self,inDir=None, sdrEndian=None, ancEndian=None):
        self.collectionShortName = 'VIIRS-GridIP-VIIRS-Qst-Lwm-Mod-Gran'
        self.xmlName = 'VIIRS_GRIDIP_VIIRS_QST_LWM_MOD_GRAN.xml'
        self.blobDatasetName = 'qstlwm'
        self.dataType = 'uint8'
        self.sourceType = 'DEM-IGBP'
        self.sourceList = ['']
        self.trimObj = ViirsData.ViirsTrimTable()

        self.GridIP_collectionShortNames = [
                            'VIIRS-GridIP-VIIRS-Qst-Mod-Gran',
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

        # QST Land Water Mask type value enumerations to be used with the qstlwm field
        self.QSTLWM_list = [
            'EVERGREEN_NEEDLELEAF_FOREST', 'EVERGREEN_BROADLEAF_FORESTS', 'DECIDUOUS_NEEDLELEAF_FORESTS', \
            'DECIDUOUS_BROADLEAF_FORESTS', 'MIXED_FORESTS', 'CLOSED_SHRUBLANDS', 'OPEN_SHRUBLANDS', \
            'WOODY_SAVANNAS', 'SAVANNAS', 'GRASSLANDS', 'PERMANENT_WETLANDS', 'CROPLANDS', \
            'URBAN_AND_BUILDUP_LANDS', 'CROPLAND_NATURAL_VEGETATION', 'SNOW_AND_ICE', \
            'BARREN_OR_SPARSELY_VEGETATED', \
            'OCEAN_SEA', 'INLAND_WATER', 'COASTAL_WATER', \
            'UNCLASSIFIED_LAND', 'FILL_VALUE'
        ]

        self.QSTLWM_dict = {
            'EVERGREEN_NEEDLELEAF_FOREST'  : 1,
            'EVERGREEN_BROADLEAF_FORESTS'  : 2,
            'DECIDUOUS_NEEDLELEAF_FORESTS' : 3,
            'DECIDUOUS_BROADLEAF_FORESTS'  : 4,
            'MIXED_FORESTS'                : 5,
            'CLOSED_SHRUBLANDS'            : 6,
            'OPEN_SHRUBLANDS'              : 7,
            'WOODY_SAVANNAS'               : 8,
            'SAVANNAS'                     : 9,
            'GRASSLANDS'                   : 10,
            'PERMANENT_WETLANDS'           : 11,
            'CROPLANDS'                    : 12,
            'URBAN_AND_BUILDUP_LANDS'      : 13, 
            'CROPLAND_NATURAL_VEGETATION'  : 14,
            'SNOW_AND_ICE'                 : 15,
            'BARREN_OR_SPARSELY_VEGETATED' : 16,
            'OCEAN_SEA'                    : 17,
            'INLAND_WATER'                 : 18,
            'COASTAL_WATER'                : 19,
            'UNCLASSIFIED_LAND'            : 20,
            'FILL_VALUE'                   : 255
        }

        self.DEM_to_QSTLWM = np.array([
                self.QSTLWM_dict['OCEAN_SEA'],                    # DEM 0 : QSTLWM 17
                self.QSTLWM_dict['EVERGREEN_NEEDLELEAF_FOREST'],  # DEM 1 : QSTLWM  1
                self.QSTLWM_dict['COASTAL_WATER'],                # DEM 2 : QSTLWM 19
                self.QSTLWM_dict['INLAND_WATER'],                 # DEM 3 : QSTLWM 18
                self.QSTLWM_dict['SAVANNAS'],                     # DEM 4 : QSTLWM  9
                self.QSTLWM_dict['OCEAN_SEA'],                    # DEM 5 : QSTLWM 17
                self.QSTLWM_dict['OCEAN_SEA'],                    # DEM 6 : QSTLWM 17
                self.QSTLWM_dict['OCEAN_SEA']                     # DEM 7 : QSTLWM 17
            ],dtype=(self.dataType))


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

        LOG.debug("###########################")
        LOG.debug("  Geolocation Information  ")
        LOG.debug("###########################")
        LOG.debug("N_Granule_ID : %r" % (N_Granule_ID))
        LOG.debug("ObservedStartTime : %s" % (ObservedStartTimeObj.__str__()))
        LOG.debug("N_Collection_Short_Name : %s" %(geo_Collection_ShortName))
        LOG.debug("URID : %r" % (URID))
        LOG.debug("geoFiles : %r" % (geoFiles))
        LOG.debug("###########################")

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
            LOG.error("Invalid geolocation shortname: %s" %(geo_Collection_ShortName))
            return -1

        # Get the geolocation xml file

        geoXmlFile = "%s.xml" % (string.replace(geo_Collection_ShortName,'-','_'))
        geoXmlFile = path.join(ADL_HOME,'xml/VIIRS',geoXmlFile)
        if path.exists(geoXmlFile):
            LOG.debug("We are using for %s: %s,%s" %(geo_Collection_ShortName,geoXmlFile,geoFiles[0]))

        # Open the geolocation blob and get the latitude and longitude

        endian = self.sdrEndian

        geoBlobObj = adl_blob.map(geoXmlFile,geoFiles[0], endian=endian)

        # Get scan_mode to find any bad scans

        scanMode = geoBlobObj.scan_mode[:]
        badScanIdx = np.where(scanMode==254)[0]
        LOG.debug("Bad Scans: %r" % (badScanIdx))

        # Detemine the min, max and range of the latitude and longitude, 
        # taking care to exclude any fill values.

        if longFormGeoNames :
            if endian==adl_blob.BIG_ENDIAN:
                latitude = getattr(geoBlobObj,'latitude').byteswap()
                longitude = getattr(geoBlobObj,'longitude').byteswap()
                latitude = latitude.astype('float')
                longitude = longitude.astype('float')
            else:
                latitude = getattr(geoBlobObj,'latitude').astype('float')
                longitude = getattr(geoBlobObj,'longitude').astype('float')
        else :
            latitude = getattr(geoBlobObj,'lat').astype('float')
            longitude = getattr(geoBlobObj,'lon').astype('float')

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

        # Shift the longitudes to be between -180 and 180 degrees
        if lonMax > 180. :
            LOG.debug("\nFinal min,max,range of longitude: %f %f %f" % (lonMin,lonMax,lonRange))
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
        self.scanMode  = scanMode
        self.latitude  = latitude
        self.longitude = longitude
        self.latCrnList  = latCrnList
        self.lonCrnList  = lonCrnList
        self.num180Crossings  = num180Crossings

        # Parse the geolocation asc file to get struct information which will be 
        # written to the ancillary asc files

        geoAscFileName = path.join(self.inDir,URID+".asc")
        LOG.debug("\nOpening %s..." % (geoAscFileName))

        geoAscFile = open(geoAscFileName,'rt')

        self.ObservedDateTimeStr =  getAscLine(geoAscFile,"ObservedDateTime")
        #self.RangeDateTimeStr =  self.ObservedDateTimeStr
        self.RangeDateTimeStr =  getAscLine(geoAscFile,"ObservedDateTime")
        self.RangeDateTimeStr =  string.replace(self.RangeDateTimeStr,"ObservedDateTime","RangeDateTime")
        self.GRingLatitudeStr =  getAscStructs(geoAscFile,"GRingLatitude",12)
        self.GRingLongitudeStr = getAscStructs(geoAscFile,"GRingLongitude",12)

        self.North_Bounding_Coordinate_Str = getAscLine(geoAscFile,"North_Bounding_Coordinate")
        self.South_Bounding_Coordinate_Str = getAscLine(geoAscFile,"South_Bounding_Coordinate")
        self.East_Bounding_Coordinate_Str  = getAscLine(geoAscFile,"East_Bounding_Coordinate")
        self.West_Bounding_Coordinate_Str  = getAscLine(geoAscFile,"West_Bounding_Coordinate")

        geoAscFile.close()


    def subset(self):
        '''Subsets the IGBP and LWM datasets to cover the required geolocation range.'''
        pass


    def __QSTLWM(self,GridIP_objects):
        '''
        Combines the LWM and IGBP to make the QSTLWM
        '''
        
        '''
        LWM  = np.array([ 0, 1, 2, 3, 4, 5, 6, 7, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
            1, 1, 1, 1, 1, 1],dtype='uint8')
        IGBP = np.array([17,17,17,17,17,17,17,17, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,
            11,12,13,14,15,16],dtype='uint8')
        QSTLWM_exp = np.array([17,19,19,18, 9,17,17,17, 1, 2, 3, 4, 5, 6, 7, 8, 
            9,10,11,12,13,14,15,16],dtype='uint8')

        '''

        DEM_list = GridIP_objects['VIIRS-GridIP-VIIRS-Lwm-Mod-Gran'].DEM_list
        DEM_dict = GridIP_objects['VIIRS-GridIP-VIIRS-Lwm-Mod-Gran'].DEM_dict
        LWM = GridIP_objects['VIIRS-GridIP-VIIRS-Lwm-Mod-Gran'].data

        IGBP_dict = GridIP_objects['VIIRS-GridIP-VIIRS-Qst-Mod-Gran'].IGBP_dict
        IGBP = GridIP_objects['VIIRS-GridIP-VIIRS-Qst-Mod-Gran'].data

        QSTLWM_dict = self.QSTLWM_dict
        DEM_to_QSTLWM = self.DEM_to_QSTLWM
        QSTLWM_dtype = self.dataType

        LOG.debug("min of LWM  = %d" % (np.min(LWM)))
        LOG.debug("max of LWM  = %d" % (np.max(LWM)))
        LOG.debug("min of IGBP = %d" % (np.min(IGBP)))
        LOG.debug("max of IGBP = %d" % (np.max(IGBP)))

        # Construct the IGBP and LWM masks

        IGBP_mask = ma.masked_greater(IGBP,17).mask
        LWM_mask = ma.masked_greater(LWM,7).mask

        # Determine the combined mask

        totalMask = np.ravel(np.bitwise_or(IGBP_mask,LWM_mask))

        # Mask and compress the IGBP and LWM datasets

        IGBP_reduced = ma.array(np.ravel(IGBP),mask=totalMask).compressed()
        LWM_reduced = ma.array(np.ravel(LWM),mask=totalMask).compressed()

        DEM_LAND_mask = (LWM_reduced == DEM_dict['DEM_LAND'])
        DEM_water_mask = np.bitwise_not(DEM_LAND_mask)
        IGBP_waterBodies_mask = (IGBP_reduced == IGBP_dict['IGBP_WATERBODIES'])
        
        QSTLWM_coastalWater_mask = np.logical_and(DEM_LAND_mask,IGBP_waterBodies_mask)
        QSTLWM_IGBP_land_mask = np.invert(IGBP_waterBodies_mask)
        
        QSTLWM_IGBP_idx = np.where(QSTLWM_IGBP_land_mask)
        QSTLWM_coastalWater_idx = np.where(QSTLWM_coastalWater_mask)
        QSTLWM_LSM_idx = np.where(DEM_water_mask)
        
        QSTLWM_reduced = np.ones(IGBP_reduced.shape,dtype=self.dataType) * QSTLWM_dict['FILL_VALUE']
        
        QSTLWM_reduced[QSTLWM_IGBP_idx] = IGBP_reduced[QSTLWM_IGBP_idx]
        QSTLWM_reduced[QSTLWM_coastalWater_idx] = np.array([QSTLWM_dict['COASTAL_WATER']],dtype=QSTLWM_dtype)[0]
        QSTLWM_reduced[QSTLWM_LSM_idx] = DEM_to_QSTLWM[LWM_reduced[QSTLWM_LSM_idx]]

        # Restore to original shape...
        QSTLWM = np.ones(IGBP.shape,dtype=QSTLWM_dtype) * QSTLWM_dict['FILL_VALUE']
        QSTLWM = np.ravel(QSTLWM)
        validIdx = np.where(totalMask==False)
        QSTLWM[validIdx] = QSTLWM_reduced
        QSTLWM = QSTLWM.reshape(IGBP.shape)

        return QSTLWM


    def granulate(self,GridIP_objects):
        '''
        Granulates the IGBP and LWM datasets, and combines them create the QSTLWM data.
        '''

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
        #GridIP_objects = {}
        for shortName in collectionShortNames :
            className = GridIP.classNames[shortName]
            GridIP_objects[shortName] = getattr(GridIP,className)(inDir=inDir,sdrEndian=self.sdrEndian)
            LOG.debug("GridIP_objects[%s].blobDatasetName = %r" % (shortName,GridIP_objects[shortName].blobDatasetName))

        # Loop through the required GridIP datasets and create the blobs.
        for shortName in collectionShortNames :
        
            LOG.debug("Processing dataset %s for %s" % (GridIP_objects[shortName].blobDatasetName,shortName))

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

            # Granulate the gridded data in this ancillary object for the current granule...
            GridIP_objects[shortName].granulate(GridIP_objects)

        # Combine the LWM and IGBP to make the QSTLWM
        data = self.__QSTLWM(GridIP_objects)

        # Convert granulated data back to original type...
        data = data.astype(self.dataType)

        LOG.debug("Shape of granulated %s data is %s" % (self.collectionShortName,np.shape(data)))

        # Explicitly restore geolocation fill to the granulated data...
        fillMask = ma.masked_less(self.latitude,-800.).mask
        fillValue = self.trimObj.sdrTypeFill['MISS_FILL'][self.dataType]        
        data = ma.array(data,mask=fillMask,fill_value=fillValue)
        data = data.filled()

        # Moderate resolution trim table arrays. These are 
        # bool arrays, and the trim pixels are set to True.
        modTrimMask = self.trimObj.createModTrimArray(nscans=48,trimType=bool)

        # Fill the required pixel trim rows in the granulated GridIP data with 
        # the ONBOARD_PT_FILL value for the correct data type

        fillValue = self.trimObj.sdrTypeFill['ONBOARD_PT_FILL'][self.dataType]        
        data = ma.array(data,mask=modTrimMask,fill_value=fillValue)
        self.data = data.filled()


    def shipOutToFile(self):
        ''' Pass the current class instance to this Utils method to generate 
            a blob/asc file pair from the input ancillary data object.'''

        return shipOutToFile(self)
