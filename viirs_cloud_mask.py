#!/usr/bin/env python
# encoding: utf-8
"""
viirs_cloud_mask.py

Purpose: Provide ingest functionality for the VIIRS cloud mask IP.

Created by Geoff Cureton on 2012-11-13.
Copyright (c) 2012 University of Wisconsin SSEC. All rights reserved.
"""

file_Date = '$Date$'
file_Revision = '$Revision$'
file_Author = '$Author$'
file_HeadURL = '$HeadURL$'
file_Id = '$Id$'

__author__ = 'G.P. Cureton <geoff.cureton@ssec.wisc.edu>'
__version__ = '$Id$'
__docformat__ = 'Epytext'


### System libraries
import time, string, sys
from glob import glob
from os import path,uname
import copy

import numpy as np
from numpy import ma as ma

import tables as pytables
from tables import exceptions as pyEx

### Local libraries
import viirs_edr_data as VD

class viirsCM:

    def __init__(self):
        '''
        This __init__ method sets static data for the class.
        '''
        self.ViirsCMbitMasks = VD.CloudMaskData.ViirsCMbitMasks

        self.ViirsCMbitShift = VD.CloudMaskData.ViirsCMbitShift

        self.ViirsCMbitMaskNames = VD.CloudMaskData.ViirsCMbitMaskNames

        # The possible values for the CM quantity being shown
        self.ViirsCMvalues = VD.CloudMaskData.ViirsCMvalues

        # The boundaries of the colourbar categories
        self.ViirsCMfillBoundaries = VD.CloudMaskData.ViirsCMfillBoundaries

        # Colours to fill in between contour boundaries
        self.ViirsCMfillColours = VD.CloudMaskData.ViirsCMfillColours

        # The tick labels of the colourbar categories
        self.ViirsCMtickNames = VD.CloudMaskData.ViirsCMtickNames


    def ingest(self,geoFile,cmFile,cmByte,cmBitFlag,shrink):
        '''
        This method ingests the lat, long and the specified byte from the
        VIIRS GMODO and IICMO files.
        '''

        # Get the input file names
        ViirsGeoFileName   = str(geoFile  )
        ViirsCMaskFileName = str(cmFile   )
        ViirsCMaskByte     = int(cmByte   )
        ViirsCMaskBitFlag  = int(cmBitFlag)
        shrinkFactor       = int(shrink   )

        self.ViirsGeoFileName = ViirsGeoFileName
        self.ViirsCMaskFileName = ViirsCMaskFileName
        self.ViirsCMaskByte = ViirsCMaskByte
        self.ViirsCMaskBitFlag = ViirsCMaskBitFlag
        self.shrinkFactor = shrinkFactor

        # Get the lat and long from the geolocation file
        try :
            if(pytables.isHDF5File(ViirsGeoFileName)) :
                try :
                    ViirsGeoFileObj = pytables.openFile(ViirsGeoFileName,mode='r')
                    print "Successfully opened geolocation file",ViirsGeoFileName
                except IOError :
                    print "error: Could not open geolocation file: ",ViirsGeoFileName
                    sys.exit(1)

                # Make a copy of the geolocation SDR, so we can close the file
                try :
                    print "Reading Latitude dataset..."
                    self.Lat = ViirsGeoFileObj.getNode('/All_Data/VIIRS-MOD-GEO_All/Latitude')[::shrinkFactor,::shrinkFactor]
                    print "done"
                    print "Reading Longitude dataset..."
                    self.Lon = ViirsGeoFileObj.getNode('/All_Data/VIIRS-MOD-GEO_All/Longitude')[::shrinkFactor,::shrinkFactor]
                    print "done"
                    print "Reading Solar Zenith Angle dataset..."
                    self.Sza = ViirsGeoFileObj.getNode('/All_Data/VIIRS-MOD-GEO_All/SolarZenithAngle')[::shrinkFactor,::shrinkFactor]
                    print "done"
                    print "Closing geolocation file"
                    ViirsGeoFileObj.close()
                except pyEx.NoSuchNodeError :
                    try :
                        print "Reading Latitude dataset..."
                        self.Lat = ViirsGeoFileObj.getNode('/All_Data/VIIRS-MOD-GEO-TC_All/Latitude')[::shrinkFactor,::shrinkFactor]
                        print "done"
                        print "Reading Longitude dataset..."
                        self.Lon = ViirsGeoFileObj.getNode('/All_Data/VIIRS-MOD-GEO-TC_All/Longitude')[::shrinkFactor,::shrinkFactor]
                        print "done"
                        print "Reading Solar Zenith Angle dataset..."
                        self.Sza = ViirsGeoFileObj.getNode('/All_Data/VIIRS-MOD-GEO-TC_All/SolarZenithAngle')[::shrinkFactor,::shrinkFactor]
                        print "done"
                        print "Closing geolocation file"
                        ViirsGeoFileObj.close()
                    except pyEx.NoSuchNodeError :
                        ViirsGeoFileObj.close()
                        errStr = "error: There is no corresponding Node in /All_Data/VIIRS-MOD-GEO-TC_All/ in %s" % (ViirsGeoFileName)
                        print errStr
                        sys.exit(1)
            else :
                print "error: %s is not a HDF5 file, aborting..."
                sys.exit(1)

        except pyEx.HDF5ExtError :
            print "error: Could not determine whether geolocation is a HDF5 file, aborting..."
            sys.exit(1)

        # Remove any unit-length dimensions from geolocation arrays
        self.Lon = np.squeeze(self.Lon)
        self.Lat = np.squeeze(self.Lat)
        self.Sza = np.squeeze(self.Sza)
        print "Shape of Latitude:",np.shape(self.Lat)
        print "Shape of Longitude:",np.shape(self.Lon)
        print "Shape of Solar Zenith Angle:",np.shape(self.Sza)

        # Get Viirs Cloud Mask file
        try :
            ViirsCMaskFileObj = pytables.openFile(ViirsCMaskFileName,'r')
            print "Successfully opened HDF5 Cloud Mask file ",ViirsCMaskFileName
        except :
            print "error: Could not open Cloud Mask file: ",ViirsCMaskFileName
            sys.exit(1)

        ViirsDataSetPath = '/All_Data/VIIRS-CM-IP_All/QF'+str(ViirsCMaskByte+1)+'_VIIRSCMIP'

        # Make a copy of the Cloud Mask SDS, so we can close the file
        try :
            self.ViirsCMaskSDS = ViirsCMaskFileObj.getNode(ViirsDataSetPath)[::shrinkFactor,::shrinkFactor]
            print "Type of self.ViirsCMaskSDS is: ",self.ViirsCMaskSDS.dtype.name

            # Get cloud mask quality
            if not (ViirsCMaskByte == 0) :
                ViirsDataSetPath = '/All_Data/VIIRS-CM-IP_All/QF1_VIIRSCMIP'
                self.ViirsCMquality = ViirsCMaskFileObj.getNode(ViirsDataSetPath)[::shrinkFactor,::shrinkFactor]
            else :
                self.ViirsCMquality = np.copy(self.ViirsCMaskSDS)
            print "Type of self.ViirsCMquality is: ",self.ViirsCMquality.dtype.name

            # Bitmask and shift VIIRS SDS
            self.ViirsCMaskSDS = np.bitwise_and(self.ViirsCMaskSDS,\
                self.ViirsCMbitMasks[ViirsCMaskByte][ViirsCMaskBitFlag])
            self.ViirsCMaskSDS = self.ViirsCMaskSDS \
            >> self.ViirsCMbitShift[ViirsCMaskByte][ViirsCMaskBitFlag]

            # Bitmask and shift CM quality
            self.ViirsCMquality = np.bitwise_and(self.ViirsCMquality,3) >> 0

            print "Closing Cloud Mask file"
            ViirsCMaskFileObj.close()

        except pyEx.NoSuchNodeError :
            print "error: There is no ",ViirsDataSetPath," dataset in",ViirsCMaskFileName
            ViirsCMaskFileObj.close()
            sys.exit(1)

        return 0
