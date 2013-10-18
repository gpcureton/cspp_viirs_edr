#!/usr/bin/env python
# encoding: utf-8
"""
viirs_vi_products.py

Purpose: Provide ingest functionality for VIIRS VI EDR.

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
from ViirsData import ViirsTrimTable
import viirs_edr_data as VD

class viirsVI:

    def __init__(self):
        '''
        This __init__ method sets static data for the class.
        '''
        pass

    def ingest(self,geoFile,viFile,SDS,shrink):
        '''
        This method ingests the lat, long and the specified SDS from the geolocation
        and VIIRS Vegetation Index EDR files.
        '''

        # Get the input file names
        ViirsGeoFileName       = str(geoFile )
        ViirsVIprodFileName    = str(viFile )
        ViirsVIprodSDStype     = str(SDS     )
        shrinkFactor           = int(shrink  )

        self.ViirsGeoFileName     = ViirsGeoFileName
        self.ViirsVIprodFileName  = ViirsVIprodFileName
        self.ViirsVIprodSDStype = SDS

        # Get the lat and long from the geolocation file
        try :
            if(pytables.isHDF5File(ViirsGeoFileName)) :
                try :
                    ViirsGeoFileObj = pytables.openFile(ViirsGeoFileName,mode='r')
                    print "Successfully opened geolocation file",ViirsGeoFileName
                except IOError :
                    print ">> error: Could not open geolocation file: ",ViirsGeoFileName
                    sys.exit(1)

                # Detemine the geolocation group name and related information
                group = ViirsGeoFileObj.getNode('/All_Data')
                geoGroupName = '/All_Data/'+group.__members__[0]
                group._g_close()

                # Make a copy of the geolocation SDR, so we can close the file
                try :
                    dataName = 'Latitude'
                    print "Reading %s dataset..." % (dataName)
                    geoNode = ViirsGeoFileObj.getNode('%s/%s' % (geoGroupName,dataName))
                    self.Lat = geoNode[::shrinkFactor,::shrinkFactor]
                    geoNode.close()
                    print "done"

                    dataName = 'Longitude'
                    print "Reading %s dataset..." % (dataName)
                    geoNode = ViirsGeoFileObj.getNode('%s/%s' % (geoGroupName,dataName))
                    self.Lon = geoNode[::shrinkFactor,::shrinkFactor]
                    geoNode.close()
                    print "done"
                    
                    #dataName = 'SolarZenithAngle'
                    #print "Reading %s dataset..." % (dataName)
                    #geoNode = ViirsGeoFileObj.getNode('%s/%s' % (geoGroupName,dataName))
                    #self.Sza = geoNode[::shrinkFactor,::shrinkFactor]
                    #geoNode.close()
                    #print "done"
                    
                    dataName = 'ModeGran'
                    print "Reading %s dataset..." % (dataName)
                    geoNode = ViirsGeoFileObj.getNode(geoGroupName+'/'+dataName)
                    self.ModeGran = geoNode.read()[0]
                    geoNode.close()
                    print "done"
                    
                    print "Closing geolocation file"
                    ViirsGeoFileObj.close()

                except pyEx.NoSuchNodeError :
                    print "error: There is no corresponding Node in /All_Data/VIIRS-MOD-GEO_All/ in ",ViirsGeoFileName
                    ViirsGeoFileObj.close()
                    sys.exit(1)

            else :
                print "error: %s is not a HDF5 file, aborting..." % (ViirsGeoFileName)
                sys.exit(1)

        except pyEx.HDF5ExtError :
            print "error: Could not determine whether geolocation isHDF4 or HDF5, aborting..."
            sys.exit(1)
        
        # Remove any unit-length dimensions from geolocation arrays
        self.Lon = np.squeeze(self.Lon)
        self.Lat = np.squeeze(self.Lat)
        print "Shape of Latitude:",np.shape(self.Lat)
        print "Shape of Longitude:",np.shape(self.Lon)
        #print "Shape of Solar Zenith Angle:",np.shape(self.Sza)

        # Open the VIIRS VI EDR file
        try :
            ViirsVIprodFileObj = pytables.openFile(ViirsVIprodFileName,mode='r')
            print "Successfully opened EDR file",ViirsVIprodFileName
        except IOError :
            print "Could not open EDR file: ",ViirsVIprodFileName
            sys.exit(1)

        # Get the VIIRS EDR dataset
        try :
            if self.ViirsVIprodSDStype == "NDVI" :
                ViirsDataSetPath = '/All_Data/VIIRS-VI-EDR_All/TOA_NDVI'
                self.ViirsVIprodSDS = ViirsVIprodFileObj.getNode(ViirsDataSetPath)[::shrinkFactor,::shrinkFactor]
                ViirsDataSetPath = '/All_Data/VIIRS-VI-EDR_All/TOA_NDVI_Factors'
                self.viFactors = ViirsVIprodFileObj.getNode(ViirsDataSetPath)[:]
            elif self.ViirsVIprodSDStype == "EVI" :
                ViirsDataSetPath = '/All_Data/VIIRS-VI-EDR_All/TOC_EVI'
                self.ViirsVIprodSDS = ViirsVIprodFileObj.getNode(ViirsDataSetPath)[::shrinkFactor,::shrinkFactor]
                ViirsDataSetPath = '/All_Data/VIIRS-VI-EDR_All/TOC_EVI_Factors'
                self.viFactors = ViirsVIprodFileObj.getNode(ViirsDataSetPath)[:]
            else :
                print 'Error: unknown product name "{}", aborting...'.format(self.ViirsVIprodSDStype)
                sys.exit(0)

            #ViirsDataSetPath = '/All_Data/VIIRS-VI-EDR_All/QF1_VIIRSVIEDR'
            #self.ViirsVI_QF1 =  ViirsVIprodFileObj.getNode(ViirsDataSetPath)[::shrinkFactor,::shrinkFactor]

            #ViirsDataSetPath = '/All_Data/VIIRS-VI-EDR_All/QF2_VIIRSVIEDR'
            #self.ViirsVI_QF2 =  ViirsVIprodFileObj.getNode(ViirsDataSetPath)[::shrinkFactor,::shrinkFactor]

            #ViirsDataSetPath = '/All_Data/VIIRS-VI-EDR_All/QF3_VIIRSVIEDR'
            #self.ViirsVI_QF3 =  ViirsVIprodFileObj.getNode(ViirsDataSetPath)[::shrinkFactor,::shrinkFactor]


            print "Closing VI Product file"
            ViirsVIprodFileObj.close()

        except pyEx.NoSuchNodeError :
            #print "error: There is no ",VIproduct.SDSname," dataset in",ViirsVIprodFileName
            print "error: There is no dataset in",ViirsVIprodFileName
            ViirsVIprodFileObj.close()
            sys.exit(1)

        return 0
