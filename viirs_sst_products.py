#!/usr/bin/env python
# encoding: utf-8
"""
viirs_sst_products.py

Purpose: Provide ingest functionality for VIIRS SST EDR.

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

class viirsSST:

    def __init__(self):
        '''
        This __init__ method sets static data for the class.
        '''
        pass

    def ingest(self,geoFile,sstFile,SDS,shrink):
        '''
        This method ingests the lat, long and the specified SDS from the geolocation
        and VIIRS Aerosol IP files.
        '''

        # Get the input file names
        ViirsGeoFileName       = str(geoFile )
        ViirsSSTprodFileName   = str(sstFile )
        ViirsSSTprodSDStype    = str(SDS     )
        shrinkFactor           = int(shrink  )

        self.ViirsGeoFileName       = ViirsGeoFileName
        self.ViirsSSTprodFileName  = ViirsSSTprodFileName
        self.ViirsSSTprodSDStype = SDS

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

        # Open the VIIRS SST EDR file
        try :
            ViirsSSTprodFileObj = pytables.openFile(ViirsSSTprodFileName,mode='r')
            print "Successfully opened EDR file",ViirsSSTprodFileName
        except IOError :
            print "Could not open EDR file: ",ViirsSSTprodFileName
            sys.exit(1)

        # Get the VIIRS EDR dataset
        try :
            #ViirsDataSetPath = '/All_Data/VIIRS-SST-EDR_All/BulkSST'
            ViirsDataSetPath = '/All_Data/VIIRS-SST-EDR_All/SkinSST'
            self.ViirsSSTprodSDS = ViirsSSTprodFileObj.getNode(ViirsDataSetPath)[::shrinkFactor,::shrinkFactor]
            #ViirsDataSetPath = '/All_Data/VIIRS-SST-EDR_All/BulkSSTFactors'
            ViirsDataSetPath = '/All_Data/VIIRS-SST-EDR_All/SkinSSTFactors'
            self.sstFactors = ViirsSSTprodFileObj.getNode(ViirsDataSetPath)[:]

            ViirsDataSetPath = '/All_Data/VIIRS-SST-EDR_All/QF1_VIIRSSSTEDR'
            self.ViirsSST_QF1 =  ViirsSSTprodFileObj.getNode(ViirsDataSetPath)[::shrinkFactor,::shrinkFactor]

            ViirsDataSetPath = '/All_Data/VIIRS-SST-EDR_All/QF2_VIIRSSSTEDR'
            self.ViirsSST_QF2 =  ViirsSSTprodFileObj.getNode(ViirsDataSetPath)[::shrinkFactor,::shrinkFactor]

            ViirsDataSetPath = '/All_Data/VIIRS-SST-EDR_All/QF3_VIIRSSSTEDR'
            self.ViirsSST_QF3 =  ViirsSSTprodFileObj.getNode(ViirsDataSetPath)[::shrinkFactor,::shrinkFactor]

            ViirsDataSetPath = '/All_Data/VIIRS-SST-EDR_All/QF4_VIIRSSSTEDR'
            self.ViirsSST_QF4 =  ViirsSSTprodFileObj.getNode(ViirsDataSetPath)[::shrinkFactor,::shrinkFactor]


            print "Closing SST Product file"
            ViirsSSTprodFileObj.close()

        except pyEx.NoSuchNodeError :
            #print "error: There is no ",SSTproduct.SDSname," dataset in",ViirsSSTprodFileName
            print "error: There is no dataset in",ViirsSSTprodFileName
            ViirsSSTprodFileObj.close()
            sys.exit(1)

        return 0
