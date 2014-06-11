#!/usr/bin/env python
# encoding: utf-8
"""
AnnMinMaxNdvi.py

 * DESCRIPTION:  Class to granulate the Annual Min/Max NDVI climatology 

Created by Geoff Cureton on 2014-05-13.
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
import h5py

import ViirsData

# skim and convert routines for reading .asc metadata fields of interest
import adl_blob
import adl_asc
from adl_asc import skim_dir, contiguous_granule_groups, \
        granule_groups_contain, effective_anc_contains,eliminate_duplicates, \
        _is_contiguous, RDR_REQUIRED_KEYS, POLARWANDER_REQUIRED_KEYS
from adl_common import ADL_HOME, CSPP_RT_HOME, CSPP_RT_ANC_HOME

# every module should have a LOG object
try :
    sourcename= file_Id.split(" ")
    LOG = logging.getLogger(sourcename[1])
except :
    LOG = logging.getLogger('AnnMinMaxNdvi')

from Utils import getURID, getAscLine, getAscStructs, \
        findDatelineCrossings, shipOutToFile
from Utils import index, find_lt, find_le, find_gt, find_ge


class AnnMinMaxNdvi() :

    def __init__(self,inDir=None, sdrEndian=None, ancEndian=None):
        self.collectionShortName = 'VIIRS-GridIP-VIIRS-Ann-Max-Min-Ndvi-Mod-Gran'
        self.xmlName = 'VIIRS_GRIDIP_VIIRS_ANN_MAX_MIN_NDVI_MOD_GRAN.xml'
        self.blobDatasetName = ['minNdvi','maxNdvi']
        #self.blobDatasetName = 'minNdvi'
        self.dataType = 'float32'
        self.sourceType = 'NDVI'
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
        if (geo_Collection_ShortName=='VIIRS-MOD-GEO-TC' 
                or geo_Collection_ShortName=='VIIRS-MOD-RGEO') :
            longFormGeoNames = True
            LOG.debug("We have long form geolocation names")
        elif (geo_Collection_ShortName=='VIIRS-MOD-GEO' 
                or geo_Collection_ShortName=='VIIRS-MOD-RGEO-TC') :
            LOG.debug("We have short form geolocation names")
            longFormGeoNames = False
        else :
            LOG.error("Invalid geolocation shortname: %s" %
                    (geo_Collection_ShortName))
            return -1

        # Get the geolocation xml file

        geoXmlFile = "%s.xml" % \
                (string.replace(geo_Collection_ShortName,'-','_'))
        geoXmlFile = path.join(ADL_HOME,'xml/VIIRS',geoXmlFile)
        if path.exists(geoXmlFile):
            LOG.debug("We are using for %s: %s,%s" %(geo_Collection_ShortName,
                geoXmlFile,geoFiles[0]))

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

        LOG.debug("min,max,range of latitide: %f %f %f" % (latMin,latMax,
            latRange))
        LOG.debug("min,max,range of longitude: %f %f %f" % (lonMin,lonMax,
            lonRange))

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

            LOG.debug("New min,max,range of latitude: %f %f %f" % (latMin,
                latMax,latRange))
            LOG.debug("New min,max,range of longitude: %f %f %f" % (lonMin,
                lonMax,lonRange))

        # Restore fill values to masked pixels in geolocation

        geoFillValue = self.trimObj. \
                sdrTypeFill['VDNE_FLOAT64_FILL'][latitude.dtype.name]
        latitude = ma.array(latitude,mask=latMask,fill_value=geoFillValue)
        self.latitude = latitude.filled()

        geoFillValue = self.trimObj. \
                sdrTypeFill['VDNE_FLOAT64_FILL'][longitude.dtype.name]
        longitude = ma.array(longitude,mask=lonMask,fill_value=geoFillValue)
        self.longitude = longitude.filled()

        # Shift the longitudes to be between -180 and 180 degrees
        if lonMax > 180. :
            LOG.debug("\nFinal min,max,range of longitude: %f %f %f" % (lonMin,
                lonMax,lonRange))
            dateLineIdx = np.where(longitude>180.)
            LOG.debug("dateLineIdx = %r" % (dateLineIdx))
            longitude[dateLineIdx] -= 360.
            lonMax = np.max(ma.array(longitude,mask=lonMask))
            lonMin = np.min(ma.array(longitude,mask=lonMask))
            lonRange = lonMax-lonMin
            LOG.debug("\nFinal min,max,range of longitude: %f %f %f" % (lonMin,
                lonMax,lonRange))

        # Record the corners, taking care to exclude any bad scans...
        nDetectors = 16
        firstGoodScan = np.where(scanMode<=2)[0][0]
        lastGoodScan = np.where(scanMode<=2)[0][-1]
        firstGoodRow = firstGoodScan * nDetectors
        lastGoodRow = lastGoodScan * nDetectors + nDetectors - 1

        latCrnList = [latitude[firstGoodRow,0],latitude[firstGoodRow,-1],
                latitude[lastGoodRow,0],latitude[lastGoodRow,-1]]
        lonCrnList = [longitude[firstGoodRow,0],longitude[firstGoodRow,-1],
                longitude[lastGoodRow,0],longitude[lastGoodRow,-1]]

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
        self.RangeDateTimeStr =  string.replace(self.RangeDateTimeStr,
                "ObservedDateTime","RangeDateTime")
        self.GRingLatitudeStr =  getAscStructs(geoAscFile,"GRingLatitude",12)
        self.GRingLongitudeStr = getAscStructs(geoAscFile,"GRingLongitude",12)

        self.North_Bounding_Coordinate_Str = getAscLine(geoAscFile,
                "North_Bounding_Coordinate")
        self.South_Bounding_Coordinate_Str = getAscLine(geoAscFile,
                "South_Bounding_Coordinate")
        self.East_Bounding_Coordinate_Str  = getAscLine(geoAscFile,
                "East_Bounding_Coordinate")
        self.West_Bounding_Coordinate_Str  = getAscLine(geoAscFile,
                "West_Bounding_Coordinate")

        geoAscFile.close()


    def subset(self):
        '''Subsets the global NDVI dataset to cover the required 
           geolocation range.'''

        geoDict = self.geoDict
        startTimeObj = geoDict['ObservedStartTime']
        LOG.debug("'ObservedStartTime' = %r" % (startTimeObj))

        # Get the subset of NDVI global dataset.

        NDVI_dLat = 60.*(1./3600.)
        NDVI_dLon = 60.*(1./3600.)

        NDVI_fileName = path.join(CSPP_RT_ANC_HOME,'NDVI/Ann-Max-Min-Ndvi.h5')
        self.sourceList.append(path.basename(NDVI_fileName))

        try :
            #NDVIobj = pytables.openFile(NDVI_fileName)
            NDVI_fileObj = h5py.File(NDVI_fileName,'r')
            #NDVI_node = NDVI_fileObj.getNode('/NDVI')
            min_ndvi_node = NDVI_fileObj['min_ndvi']
            max_ndvi_node = NDVI_fileObj['max_ndvi']
        except Exception, err :
            LOG.exception("%s"%(err))
            LOG.exception("Problem opening NDVI file (%s), aborting."
                    %(NDVI_fileName))
            sys.exit(1)

        try :
            #NDVI_gridLats = NDVIobj.getNode('/Latitude')[:]
            #NDVI_gridLons = NDVIobj.getNode('/Longitude')[:]
            NDVI_gridLats = NDVI_fileObj['Latitude'].value
            NDVI_gridLons = NDVI_fileObj['Longitude'].value

            LOG.debug("min,max NDVI Grid Latitude values : %f,%f"
                    %(NDVI_gridLats[0],NDVI_gridLats[-1]))
            LOG.debug("min,max NDVI Grid Longitude values : %f,%f"
                    %(NDVI_gridLons[0],NDVI_gridLons[-1]))

            #NDVI_scaleFactor = NDVI_node.attrs['scale_factor']
            #NDVI_offset = NDVI_node.attrs['add_offset']
            #NDVI_fillValue = NDVI_node.attrs['_FillValue']
            NDVI_scaleFactor = min_ndvi_node.attrs['scale_factor']
            NDVI_offset = min_ndvi_node.attrs['add_offset']
            NDVI_fillValue = int(min_ndvi_node.attrs['_FillValue'])


            latMin = self.latMin
            latMax = self.latMax
            lonMin = self.lonMin
            lonMax = self.lonMax

            NDVI_latMask = np.equal((NDVI_gridLats<(latMax+NDVI_dLat)),
                    (NDVI_gridLats>(latMin-NDVI_dLat)))
            NDVI_lonMask = np.equal((NDVI_gridLons<(lonMax+NDVI_dLon)),
                    (NDVI_gridLons>(lonMin-NDVI_dLon)))

            NDVI_latIdx = np.where(NDVI_latMask==True)[0]
            NDVI_lonIdx = np.where(NDVI_lonMask==True)[0]

            NDVI_latMinIdx = NDVI_latIdx[0]
            NDVI_latMaxIdx = NDVI_latIdx[-1]
            NDVI_lonMinIdx = NDVI_lonIdx[0]
            NDVI_lonMaxIdx = NDVI_lonIdx[-1]

            LOG.debug("NDVI_latMinIdx = %d" % (NDVI_latMinIdx))
            LOG.debug("NDVI_latMaxIdx = %d" % (NDVI_latMaxIdx))
            LOG.debug("NDVI_lonMinIdx = %d" % (NDVI_lonMinIdx))
            LOG.debug("NDVI_lonMaxIdx = %d" % (NDVI_lonMaxIdx))

            lat_subset = NDVI_gridLats[NDVI_latMinIdx:NDVI_latMaxIdx+1]
            self.gridLat = lat_subset

            if self.num180Crossings == 2 :

                # We have a dateline crossing, so subset the positude and negative
                # longitude grids and sandwich them together.
                posLonCrn = np.min(ma.masked_less_equal(
                    np.array(self.lonCrnList),0.))
                negLonCrn = np.max(ma.masked_outside(
                    np.array(self.lonCrnList),-800.,0.))
                posIdx = index(NDVI_gridLons,find_lt(NDVI_gridLons,posLonCrn))
                negIdx = index(NDVI_gridLons,find_gt(NDVI_gridLons,negLonCrn))

                posLons_subset = NDVI_gridLons[posIdx:]
                negLons_subset = NDVI_gridLons[:negIdx]
                lon_subset = np.concatenate((posLons_subset,negLons_subset))

                # Do the same with the minimum NDVI data
                posBlock = min_ndvi_node[NDVI_latMinIdx:NDVI_latMaxIdx+1,posIdx:]
                negBlock = min_ndvi_node[NDVI_latMinIdx:NDVI_latMaxIdx+1,:negIdx]
                min_NDVI_subset = np.concatenate((posBlock,negBlock),axis=1)

                # Do the same with the maximum NDVI data
                posBlock = max_ndvi_node[NDVI_latMinIdx:NDVI_latMaxIdx+1,posIdx:]
                negBlock = max_ndvi_node[NDVI_latMinIdx:NDVI_latMaxIdx+1,:negIdx]
                max_NDVI_subset = np.concatenate((posBlock,negBlock),axis=1)

            else :

                min_NDVI_subset = min_ndvi_node[NDVI_latMinIdx:NDVI_latMaxIdx+1,
                        NDVI_lonMinIdx:NDVI_lonMaxIdx+1]
                max_NDVI_subset = max_ndvi_node[NDVI_latMinIdx:NDVI_latMaxIdx+1,
                        NDVI_lonMinIdx:NDVI_lonMaxIdx+1]
                lon_subset = NDVI_gridLons[NDVI_lonMinIdx:NDVI_lonMaxIdx+1]

            self.gridLon = lon_subset

            # Unscale the NDVI data, and copy to the GridIP object
            self.gridData = []

            min_NDVI_mask = ma.masked_equal(min_NDVI_subset,NDVI_fillValue).mask
            min_NDVI_subset = min_NDVI_subset * NDVI_scaleFactor + NDVI_offset
            min_NDVI_subset = ma.array(min_NDVI_subset,mask=min_NDVI_mask,
                    fill_value=-999.9)
            min_NDVI_subset = min_NDVI_subset.filled()
            LOG.debug("min_NDVI_subset.shape = {}".format(min_NDVI_subset.shape))
            self.gridData.append(min_NDVI_subset.astype(self.dataType))

            max_NDVI_mask = ma.masked_equal(max_NDVI_subset,NDVI_fillValue).mask
            max_NDVI_subset = max_NDVI_subset * NDVI_scaleFactor + NDVI_offset
            max_NDVI_subset = ma.array(max_NDVI_subset,mask=max_NDVI_mask,
                    fill_value=-999.9)
            max_NDVI_subset = max_NDVI_subset.filled()
            LOG.debug("max_NDVI_subset.shape = {}".format(max_NDVI_subset.shape))
            self.gridData.append(max_NDVI_subset.astype(self.dataType))

            #NDVI_node.close()
            #NDVIobj.close()
            NDVI_fileObj.close()

        except Exception, err :

            LOG.debug("EXCEPTION: %s" % (err))
            LOG.debug(traceback.format_exc())
            #NDVI_node.close()
            #NDVIobj.close()
            NDVI_fileObj.close()

        LOG.debug("self.gridData.shape = {}".format(np.shape(self.gridData)))


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
        '''
        Granulate the GridIP NDVI dataset.
        '''

        # Generate the lat and lon grids, and flip them and the data 
        # over latitude
        gridLon,gridLat = np.meshgrid(self.gridLon,self.gridLat[::-1])
        gridData = self.gridData
        LOG.debug("gridData.shape = {}".format(np.shape(gridData)))
        gridData[0] = gridData[0][::-1,:]
        gridData[1] = gridData[1][::-1,:]
        LOG.debug("gridData.shape = {}".format(np.shape(gridData)))

        latitude = self.latitude
        longitude = self.longitude

        # If we have a dateline crossing, remove the longitude discontinuity
        # by adding 360 degrees to the negative longitudes.
        if self.num180Crossings == 2 :
            gridLonNegIdx = np.where(gridLon < 0.)
            gridLon[gridLonNegIdx] += 360.
            longitudeNegIdx = np.where(longitude < 0.)
            longitude[longitudeNegIdx] += 360.

        LOG.info("Granulating %s ..." % (self.collectionShortName))
        LOG.debug("latitide,longitude shapes: %s, %s"%(str(latitude.shape),
            str(longitude.shape)))
        LOG.debug("gridData.shape = {}".format(np.shape(gridData)))
        LOG.debug("gridLat.shape = %s" % (str(gridLat.shape)))
        LOG.debug("gridLon.shape = %s" % (str(gridLon.shape)))

        LOG.debug("min of gridData  = {}, {}".format(np.min(gridData[0]),
            np.min(gridData[1])))
        LOG.debug("max of gridData  = {}, {}".format(np.max(gridData[0]),
            np.max(gridData[1])))

        self.data = []

        for idx in range(2):
            t1 = time()
            data,dataIdx = self._grid2Gran(np.ravel(latitude),
                                      np.ravel(longitude),
                                      gridData[idx].astype(np.float64),
                                      gridLat.astype(np.float64),
                                      gridLon.astype(np.float64))
            t2 = time()
            elapsedTime = t2-t1
            LOG.info("Granulation took %f seconds for %d points" % 
                    (elapsedTime,latitude.size))

            data = data.reshape(latitude.shape)
            dataIdx = dataIdx.reshape(latitude.shape)

            # Convert granulated data back to original type...
            data = data.astype(self.dataType)

            LOG.debug("Shape of granulated %s data is %s" % 
                    (self.collectionShortName,np.shape(data)))
            LOG.debug("Shape of granulated %s dataIdx is %s" % 
                    (self.collectionShortName,np.shape(dataIdx)))

            # Explicitly restore geolocation fill to the granulated data...
            fillMask = ma.masked_less(self.latitude,-800.).mask
            fillValue = self.trimObj.sdrTypeFill['MISS_FILL'][self.dataType]        
            data = ma.array(data,mask=fillMask,fill_value=fillValue)
            data = data.filled()

            # What value are the bowtie deletion pixels
            ongroundPixelTrimValue = self.trimObj \
                    .sdrTypeFill['ONGROUND_PT_FILL'][data.dtype.name]
            LOG.debug("Onground Pixel Trim value is {}".
                    format(ongroundPixelTrimValue))
            onboardPixelTrimValue = self.trimObj \
                    .sdrTypeFill['ONBOARD_PT_FILL'][data.dtype.name]
            LOG.debug("Onboard Pixel Trim value is {}"
                    .format(onboardPixelTrimValue))

            onboardTrimMask = self.trimObj.createOnboardModTrimArray(nscans=48,
                    trimType=bool)
            ongroundTrimMask = self.trimObj.createModTrimArray(nscans=48,
                    trimType=bool)

            # FIXME: Get correct order of application of the on-board and 
            # on-ground pixel trim.

            # Apply the On-board pixel trim
            data = ma.array(data,mask=onboardTrimMask,
                    fill_value=ongroundPixelTrimValue)
            data = data.filled() # Substitute for the masked values with 
                                 # ongroundPixelTrimValue

            # Apply the On-ground pixel trim
            data = ma.array(data,mask=ongroundTrimMask,
                    fill_value=onboardPixelTrimValue)
            data = data.filled() # Substitute for the masked values with 
                                 # onboardPixelTrimValue

            self.data.append(data)


    def shipOutToFile(self):
        '''
        Generate a blob/asc file pair from the input ancillary data object.
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
        newGridIPblobObj = adl_blob.create(xmlName, blobName, endian=endian, 
                overwrite=True)
        newGridIPblobArrObj = newGridIPblobObj.as_arrays()

        for idx in range(2):
            blobData = getattr(newGridIPblobArrObj,self.blobDatasetName[idx])
            blobData[:,:] = self.data[idx][:,:]

        # Make a new GridIP asc file from the template, and substitute for the 
        # various tags

        ascTemplateFileName = path.join(ADL_ASC_TEMPLATES,
                "VIIRS-GridIP-VIIRS_Template.asc")

        LOG.debug("Creating new asc file\n%s\nfrom template\n%s" % 
                (ascFileName,ascTemplateFileName))
        
        ANC_fileList = self.sourceList
        for idx in range(len(ANC_fileList)) :
            ANC_fileList[idx] = path.basename(ANC_fileList[idx])
        ANC_fileList.sort()
        ancGroupRecipe = '    ("N_Anc_Filename" STRING EQ "%s")'
        ancFileStr = "%s" % ("\n").join([ancGroupRecipe % (str(files)) 
            for files in ANC_fileList])

        LOG.debug("RangeDateTimeStr = %s\n" % (self.RangeDateTimeStr))
        LOG.debug("GRingLatitudeStr = \n%s\n" % (self.GRingLatitudeStr))
        LOG.debug("GRingLongitudeStr = \n%s\n" % (self.GRingLongitudeStr))

        try:
            # Open template file for reading
            ascTemplateFile = open(ascTemplateFileName,"rt") 
            # create a new text file
            ascFile = open(ascFileName,"wt") 
        except Exception, err :
            LOG.error("%s, aborting." % (err))
            sys.exit(1)

        LOG.debug("Template file %s is %r with mode %s" %(ascTemplateFileName,
            'not open' if ascTemplateFile.closed else 'open',ascTemplateFile.mode))

        LOG.debug("New file %s is %r with mode %s" %(ascFileName,
            'not open' if ascFile.closed else 'open',ascFile.mode))

        for line in ascTemplateFile.readlines():
           line = line.replace("CSPP_URID",URID)
           line = line.replace("CSPP_CREATIONDATETIME_NOUSEC",
                   creationDate_nousecStr)
           line = line.replace("CSPP_ANC_BLOB_FULLPATH",path.basename(blobName))
           line = line.replace("CSPP_ANC_COLLECTION_SHORT_NAME",
                   self.collectionShortName)
           line = line.replace("CSPP_GRANULE_ID",self.geoDict['N_Granule_ID'])
           line = line.replace("CSPP_CREATIONDATETIME",creationDateStr)
           line = line.replace("  CSPP_RANGE_DATE_TIME",self.RangeDateTimeStr)
           line = line.replace("  CSPP_GRINGLATITUDE",self.GRingLatitudeStr)
           line = line.replace("  CSPP_GRINGLONGITUDE",self.GRingLongitudeStr)
           line = line.replace("CSPP_NORTH_BOUNDING_COORD",
                   self.North_Bounding_Coordinate_Str)
           line = line.replace("CSPP_SOUTH_BOUNDING_COORD",
                   self.South_Bounding_Coordinate_Str)
           line = line.replace("CSPP_EAST_BOUNDING_COORD",
                   self.East_Bounding_Coordinate_Str)
           line = line.replace("CSPP_WEST_BOUNDING_COORD",
                   self.West_Bounding_Coordinate_Str)
           line = line.replace("CSPP_ANC_ENDIANNESS",endianString)       
           #line = line.replace("    CSPP_ANC_SOURCE_FILES",ancFileStr)
           ascFile.write(line) 

        ascFile.close()
        ascTemplateFile.close()

        return URID
