#!/usr/bin/env python
# encoding: utf-8
"""
LandWaterMask.py

 * DESCRIPTION:  Class to granulate the DEM Land Water Mask data product 

Created by Geoff Cureton on 2013-03-05.
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
import adl_blob
import adl_asc
from adl_asc import skim_dir, contiguous_granule_groups, granule_groups_contain, effective_anc_contains,_eliminate_duplicates,_is_contiguous, RDR_REQUIRED_KEYS, POLARWANDER_REQUIRED_KEYS
from adl_common import ADL_HOME, CSPP_RT_ANC_PATH, CSPP_RT_ANC_CACHE_DIR, COMMON_LOG_CHECK_TABLE

# every module should have a LOG object
try :
    sourcename= file_Id.split(" ")
    LOG = logging.getLogger(sourcename[1])
except :
    LOG = logging.getLogger('LandWaterMask')

from Utils import getURID, getAscLine, getAscStructs, shipOutToFile
from Utils import index, find_lt, find_le, find_gt, find_ge

class LandWaterMask() :

    def __init__(self,inDir=None, sdrEndian=None, ancEndian=None):
        self.collectionShortName = 'VIIRS-GridIP-VIIRS-Lwm-Mod-Gran'
        self.xmlName = 'VIIRS_GRIDIP_VIIRS_LWM_MOD_GRAN.xml'
        self.blobDatasetName = 'landWaterMask'
        self.dataType = 'uint8'
        self.sourceType = 'DEM'
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

        # Digital Elevation Model (DEM) land sea mask types
        self.DEM_list = ['DEM_SHALLOW_OCEAN','DEM_LAND','DEM_COASTLINE','DEM_SHALLOW_INLAND_WATER','DEM_EPHEMERAL_WATER','DEM_DEEP_INLAND_WATER','DEM_MOD_CONT_OCEAN','DEM_DEEP_OCEAN']
        self.DEM_dict = {
            'DEM_SHALLOW_OCEAN'        : 0,
            'DEM_LAND'                 : 1,
            'DEM_COASTLINE'            : 2,
            'DEM_SHALLOW_INLAND_WATER' : 3,
            'DEM_EPHEMERAL_WATER'      : 4,
            'DEM_DEEP_INLAND_WATER'    : 5,
            'DEM_MOD_CONT_OCEAN'       : 6,
            'DEM_DEEP_OCEAN'           : 7
        }


    def subset(self,latMinList,latMaxList,lonMinList,lonMaxList,latCrnList,lonCrnList):
        '''Subsets the global elevation dataset to cover the required geolocation range.'''

        CSPP_RT_ANC_HOME = path.abspath(os.getenv('CSPP_RT_ANC_HOME'))
        
        DEM_dLat = 30.*(1./3600.)
        DEM_dLon = 30.*(1./3600.)

        DEM_gridLats = -1. * (np.arange(21600.) * DEM_dLat - 90.)
        DEM_gridLons = np.arange(43200.) * DEM_dLon - 180.

        # Get the subset of DEM global dataset.
        DEM_fileName = path.join(CSPP_RT_ANC_HOME,'LSM/dem30ARC_Global_LandWater_uncompressed.h5')

        dateLineCrossed,ascendingNode,descendingNode = isDatelineCrossed(latCrnList,lonCrnList)
        LOG.debug("dateLineCross is %r" % (dateLineCrossed))
        LOG.debug("ascendingNode is %r" % (ascendingNode))
        LOG.debug("descendingNode is %r\n" % (descendingNode))

        latMin = min(latMinList)
        latMax = max(latMaxList)
        lonMin = min(lonMinList)
        lonMax = max(lonMaxList)

        DEM_latMask = np.equal((DEM_gridLats<(latMax+DEM_dLat)),(DEM_gridLats>(latMin-DEM_dLat)))
        DEM_lonMask = np.equal((DEM_gridLons<(lonMax+DEM_dLon)),(DEM_gridLons>(lonMin-DEM_dLon)))

        DEM_latIdx = np.where(DEM_latMask==True)[0]
        DEM_lonIdx = np.where(DEM_lonMask==True)[0]

        DEM_latMinIdx = DEM_latIdx[0]
        DEM_latMaxIdx = DEM_latIdx[-1]
        DEM_lonMinIdx = DEM_lonIdx[0]
        DEM_lonMaxIdx = DEM_lonIdx[-1]

        LOG.debug("Opening DEM file %s" % (DEM_fileName))
        # TODO : Use original HDF4 file which contains elevation and LWM.
        DEMobj = pytables.openFile(DEM_fileName,'r')
        DEM_node = DEMobj.getNode('/demGRID/Data Fields/LandWater')

        try :

            lat_subset = DEM_gridLats[DEM_latMinIdx:DEM_latMaxIdx+1]

            if True in dateLineCrossed :
                posLonCrn = np.min(ma.masked_less_equal(np.array(lonCrnList),0.))
                negLonCrn = np.max(ma.masked_outside(np.array(lonCrnList),-800.,0.))
                posIdx = index(DEM_gridLons,find_lt(DEM_gridLons,posLonCrn))
                negIdx = index(DEM_gridLons,find_gt(DEM_gridLons,negLonCrn))

                posBlock = DEM_node[DEM_latMinIdx:DEM_latMaxIdx+1,posIdx:]
                negBlock = DEM_node[DEM_latMinIdx:DEM_latMaxIdx+1,:negIdx]

                DEM_subset = np.concatenate((posBlock,negBlock),axis=1)

                posLons_subset = DEM_gridLons[posIdx:]
                negLons_subset = DEM_gridLons[:negIdx]
                lon_subset = np.concatenate((posLons_subset,negLons_subset))

            else :

                DEM_subset = DEM_node[DEM_latMinIdx:DEM_latMaxIdx+1,DEM_lonMinIdx:DEM_lonMaxIdx+1]
                lon_subset = DEM_gridLons[DEM_lonMinIdx:DEM_lonMaxIdx+1]

            DEM_node.close()
            DEMobj.close()

        except Exception, err :

            LOG.debug("EXCEPTION: %s" % (err))

            DEM_node.close()
            DEMobj.close()

        return DEM_subset.astype('uint8'),lat_subset,lon_subset,DEM_fileName


    def granulate(self,geoDicts):
        '''Granulates the input DEM files.'''

        # Moderate resolution trim table arrays. These are 
        # bool arrays, and the trim pixels are set to True.

        trimObj = ViirsData.ViirsTrimTable()
        modTrimMask = trimObj.createModTrimArray(nscans=48,trimType=bool)

        # Get a bunch of information about the geolocation
        latitudeList,longitudeList,latMinList,latMaxList,lonMinList,lonMaxList,latCrnList,lonCrnList = _get_geo_Arrays(geoDicts)

        LOG.debug("\nGranules -->")
        LOG.debug("latMin : %r" % (latMinList))
        LOG.debug("latMax : %r" % (latMaxList))
        LOG.debug("lonMin : %r" % (lonMinList))
        LOG.debug("lonMax : %r" % (lonMaxList))
        LOG.debug("latCrnList : %r" % (latCrnList))
        LOG.debug("lonCrnList : %r\n" % (lonCrnList))

        DEM_subset,lat_subset,lon_subset,DEM_fileName = \
                _subset_DEM(latMinList,latMaxList,lonMinList,lonMaxList,latCrnList,lonCrnList)

        gridLon,gridLat = np.meshgrid(lon_subset,lat_subset[::-1])
        DEM_subset = DEM_subset[::-1,:]

        DEM_type = DEM_subset.dtype

        DEM_list = []

        for latitude,longitude,geoDict in zip(latitudeList,longitudeList,geoDicts):

            N_Granule_ID = geoDict['N_Granule_ID']

            LOG.debug("\nGranulating %s ..." % ('DEM'))
            LOG.debug("latitide,longitude shapes: %s, %s"%(str(latitude.shape) , str(longitude.shape)))
            LOG.debug("DEM_subset.shape = %s" % (str(DEM_subset.shape)))
            LOG.debug("gridLat.shape = %s" % (str(gridLat.shape)))
            LOG.debug("gridLon.shape = %s" % (str(gridLon.shape)))

            LOG.debug("min of DEM_subset  = %s"%(np.min(DEM_subset)))
            LOG.debug("max of DEM_subset  = %s"%(np.max(DEM_subset)))

            data,dataIdx = _grid2Gran(np.ravel(latitude),
                                      np.ravel(longitude),
                                      DEM_subset.astype(np.float64),
                                      gridLat.astype(np.float64),
                                      gridLon.astype(np.float64))

            data = data.reshape(latitude.shape)
            dataIdx = dataIdx.reshape(latitude.shape)
            LOG.debug("Shape of first granulated %s data is %s" % ('DEM',np.shape(data)))
            LOG.debug("Shape of first granulated %s dataIdx is %s" % ('DEM',np.shape(dataIdx)))

            # Convert granulated data back to original type...

            data = data.astype(DEM_type)

            # Fill the required pixel trim rows in the granulated NCEP data with 
            # the ONBOARD_PT_FILL value for the correct data type

            fillValue = trimObj.sdrTypeFill['ONBOARD_PT_FILL'][data.dtype.name]        
            data = ma.array(data,mask=modTrimMask,fill_value=fillValue)
            LOG.debug("min of DEM granule = %d"%(np.min(data)))
            LOG.debug("max of DEM granule = %d"%(np.max(data)))

            data = data.filled()

            DEM_list.append(data)

        return DEM_list


    def shipOutToFile(self):
        ''' Pass the current class instance to this Utils method to generate 
            a blob/asc file pair from the input ancillary data object.'''

        shipOutToFile(self)
