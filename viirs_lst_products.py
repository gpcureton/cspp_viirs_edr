#!/usr/bin/env python
# encoding: utf-8
"""
viirs_lst_products.py

Purpose: Provide ingest functionality for VIIRS LST EDR.

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

class viirsLST:

    def __init__(self):
        '''
        This __init__ method sets static data for the class.
        '''
        pass

    def ingest(self,geoFile,lstFile,SDS,shrink):
        '''
        This method ingests the lat, long and the specified SDS from the geolocation
        and VIIRS Aerosol IP files.
        '''

        # Get the input file names
        ViirsGeoFileName       = str(geoFile )
        ViirsLSTprodFileName   = str(lstFile )
        ViirsLSTprodSDStype    = str(SDS     )
        shrinkFactor           = int(shrink  )

        self.ViirsGeoFileName       = ViirsGeoFileName
        self.ViirsLSTprodFileName  = ViirsLSTprodFileName
        self.ViirsLSTprodSDStype = SDS

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

        # Open the VIIRS LST EDR file
        try :
            ViirsLSTprodFileObj = pytables.openFile(ViirsLSTprodFileName,mode='r')
            print "Successfully opened EDR file",ViirsLSTprodFileName
        except IOError :
            print "Could not open EDR file: ",ViirsLSTprodFileName
            sys.exit(1)

        # Get the VIIRS EDR dataset
        try :
            print "Reading {}".format(0)
            ViirsDataSetPath = '/All_Data/VIIRS-LST-EDR_All/LandSurfaceTemperature'
            self.ViirsLSTprodSDS = ViirsLSTprodFileObj.getNode(ViirsDataSetPath)[::shrinkFactor,::shrinkFactor]
            print "Reading {}".format(1)
            ViirsDataSetPath = '/All_Data/VIIRS-LST-EDR_All/LSTFactors'
            self.lstFactors = ViirsLSTprodFileObj.getNode(ViirsDataSetPath)[:]

            ViirsDataSetPath = '/All_Data/VIIRS-LST-EDR_All/QF1_VIIRSLSTEDR'
            self.ViirsLST_QF1 =  ViirsLSTprodFileObj.getNode(ViirsDataSetPath)[::shrinkFactor,::shrinkFactor]

            ViirsDataSetPath = '/All_Data/VIIRS-LST-EDR_All/QF2_VIIRSLSTEDR'
            self.ViirsLST_QF2 =  ViirsLSTprodFileObj.getNode(ViirsDataSetPath)[::shrinkFactor,::shrinkFactor]

            ViirsDataSetPath = '/All_Data/VIIRS-LST-EDR_All/QF3_VIIRSLSTEDR'
            self.ViirsLST_QF3 =  ViirsLSTprodFileObj.getNode(ViirsDataSetPath)[::shrinkFactor,::shrinkFactor]

            print "Reading {}".format(2)

            print "Closing LST Product file"
            ViirsLSTprodFileObj.close()

        #except Exception, err :
        except pyEx.NoSuchNodeError :
            print traceback.format_exc()
            print "error: There is no dataset in",ViirsLSTprodFileName
            ViirsLSTprodFileObj.close()
            sys.exit(1)

        return 0
