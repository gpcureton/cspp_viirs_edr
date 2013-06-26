#!/usr/bin/env python
# encoding: utf-8
"""
viirs_cloud_products.py

Purpose: Provide ingest functionality for VIIRS cloud IP and EDR.

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

class viirsCld:

    def __init__(self):
        '''
        This __init__ method sets static data for the class.
        '''
        pass

    def ingestLite(self,geoFile,cldFile,SDS,shrink,scale='log'):
        '''
        This method ingests the lat, long and the specified SDS from the geolocation
        and VIIRS Cloud file.
        '''

        # Get the input file names
        ViirsGeoFileName   = str(geoFile)
        ViirsCProdFileName = str(cldFile)
        ViirsCProdSDStype  = str(SDS    )
        shrinkFactor       = int(shrink )
        scale              = str(scale  )

        self.ViirsGeoFileName    = ViirsGeoFileName
        self.ViirsCProdFileName  = ViirsCProdFileName
        self.ViirsCProdSDStype = SDS
        self.scale = scale

        # Get the various parameters related to the current product
        CloudProduct = VD.CloudProdData.CloudProduct[ViirsCProdSDStype]
        if isinstance(CloudProduct,VD.CloudProdData.CloudProd) :
            pass
        else :
            CloudProduct = VD.CloudProdData.CloudProduct[ViirsCProdSDStype][scale] # cot and eps

        # Get the lat and long from the geolocation file
        try :
            if(pytables.isHDF5File(ViirsGeoFileName)) :
                try :
                    ViirsGeoFileObj = pytables.openFile(ViirsGeoFileName,mode='r')
                    print "Successfully opened geolocation file",ViirsGeoFileName
                except IOError :
                    print "error: Could not open geolocation file: ",ViirsGeoFileName
                    sys.exit(1)

                # Detemine the sdr group name and related information
                group = ViirsGeoFileObj.getNode('/All_Data')
                geoGroupName = '/All_Data/'+group.__members__[0]
                group._g_close()

                # Make a copy of the geolocation SDR, so we can close the file
                try :
                    print "Reading Latitude dataset..."
                    dataName = 'Latitude'
                    geoNode = ViirsGeoFileObj.getNode('%s/%s' % (geoGroupName,dataName))
                    self.Lat = geoNode[::shrinkFactor,::shrinkFactor]
                    geoNode.close()
                    print "done"
                    print "Reading Longitude dataset..."
                    dataName = 'Longitude'
                    geoNode = ViirsGeoFileObj.getNode('%s/%s' % (geoGroupName,dataName))
                    self.Lon = geoNode[::shrinkFactor,::shrinkFactor]
                    geoNode.close()
                    print "done"
                    print "Reading SolarZenithAngle dataset..."
                    dataName = 'SolarZenithAngle'
                    geoNode = ViirsGeoFileObj.getNode('%s/%s' % (geoGroupName,dataName))
                    self.Sza = geoNode[::shrinkFactor,::shrinkFactor]
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
        self.Sza = np.squeeze(self.Sza)
        print "Shape of Latitude:",np.shape(self.Lat)
        print "Shape of Longitude:",np.shape(self.Lon)
        print "Shape of Solar Zenith Angle:",np.shape(self.Sza)

        # Open the VIIRS Cloud IP file
        try :
            ViirsCProdFileObj = pytables.openFile(ViirsCProdFileName,mode='r')
            print "Successfully opened IP file",ViirsCProdFileName
        except IOError :
            print "error: Could not open IP file: ",ViirsCProdFileName
            sys.exit(1)

        # Get the VIIRS IP dataset
        try :
            sdsName = CloudProduct.SDSname
            self.ViirsCProdSDS = ViirsCProdFileObj.getNode(sdsName)[::shrinkFactor,::shrinkFactor]
                #CloudProduct.SDSname)[::shrinkFactor,::shrinkFactor]
            print "Type of self.ViirsCProdSDS is: ",self.ViirsCProdSDS.dtype.name
            print "Shape of ViirsCProdSDS:",np.shape(self.ViirsCProdSDS)

            sdsName = CloudProduct.SDSphaseName
            self.ViirsCProdSDSphase = ViirsCProdFileObj.getNode(sdsName)[::shrinkFactor,::shrinkFactor]
                #CloudProduct.SDSphaseName)[::shrinkFactor,::shrinkFactor]
            print "Type of self.ViirsCProdSDSphase is: ",self.ViirsCProdSDSphase.dtype.name

            if self.ViirsCProdSDStype == 'cot' or self.ViirsCProdSDStype == 'eps':
                self.ViirsCProdSDSphase = np.bitwise_and(self.ViirsCProdSDSphase,224) >> 5

            if self.ViirsCProdSDStype == 'ctp' or self.ViirsCProdSDStype == 'ctt' or self.ViirsCProdSDStype == 'cth':
                self.ViirsCProdSDSphase = np.bitwise_and(self.ViirsCProdSDSphase,7)

            print "Type of self.ViirsCProdSDSphase is: ",self.ViirsCProdSDSphase.dtype.name
            print "Shape of ViirsCProdSDSphase:",np.shape(self.ViirsCProdSDSphase)

            print "Closing Cloud Product file"
            ViirsCProdFileObj.close()

        except pyEx.NoSuchNodeError :
            print "error: There is no ",sdsName," dataset in",ViirsCProdFileName
            ViirsCProdFileObj.close()
            sys.exit(1)

        return 0

    def ingest(self,geoFile,cmFile,cldFile,SDS,shrink,scale):
        '''
        This method ingests the lat, long and the specified SDS from the geolocation
        and VIIRS Cloud file.
        '''

        # Get the input file names
        ViirsGeoFileName   = str(geoFile)
        ViirsCMaskFileName = str(cmFile )
        ViirsCProdFileName = str(cldFile)
        ViirsCProdSDStype  = str(SDS    )
        shrinkFactor       = int(shrink )
        scale              = str(scale  )

        self.ViirsGeoFileName       = ViirsGeoFileName
        self.ViirsCProdFileName  = ViirsCProdFileName
        self.ViirsCProdSDStype = SDS
        self.scale = scale

        # Get the various parameters related to the current product
        CloudProduct = VD.CloudProdData.CloudProduct[ViirsCProdSDStype]
        if isinstance(CloudProduct,VD.CloudProdData.CloudProd) :
            pass
        else :
            CloudProduct = VD.CloudProdData.CloudProduct[ViirsCProdSDStype][scale] # cot and eps

        # Get the lat and long from the geolocation file
        try :
            if(pytables.isHDF5File(ViirsGeoFileName)) :
                try :
                    ViirsGeoFileObj = pytables.openFile(ViirsGeoFileName,mode='r')
                    print "Successfully opened geolocation file",ViirsGeoFileName
                except IOError :
                    print "Could not open geolocation file: ",ViirsGeoFileName
                    if __name__=='__main__':
                        sys.exit(1)
                    else :
                        return 1

                # Make a copy of the geolocation SDR, so we can close the file
                try :
                    self.Lat = np.copy(np.array(ViirsGeoFileObj.getNode('/All_Data/VIIRS-MOD-GEO_All/Latitude')))[::shrinkFactor,::shrinkFactor]
                    self.Lon = np.copy(np.array(ViirsGeoFileObj.getNode('/All_Data/VIIRS-MOD-GEO_All/Longitude')))[::shrinkFactor,::shrinkFactor]
                    print "Closing geolocation file"
                    ViirsGeoFileObj.close()
                except pyEx.NoSuchNodeError :
                    print "There is no corresponding Node in /All_Data/VIIRS-MOD-GEO_All/ in ",ViirsGeoFileName
                    ViirsGeoFileObj.close()
                    if __name__=='__main__':
                        sys.exit(1)
                    else :
                        return 1
            else :
                try :
                    ViirsGeoFileObj = SD( ViirsGeoFileName, SDC.READ )
                    print "Successfully opened geolocation file",ViirsGeoFileName
                except pyHDF.HDF4Error:
                    print "Could not open geolocation file: ",ViirsGeoFileName
                    if __name__=='__main__':
                        sys.exit(1)
                    else :
                        return 1

                # Make a copy of the geolocation SDR, so we can close the file
                try :
                    self.Lat = np.copy(ViirsGeoFileObj.select('Latitude')[::shrinkFactor,::shrinkFactor])
                    self.Lon = np.copy(ViirsGeoFileObj.select('Longitude')[::shrinkFactor,::shrinkFactor])
                    ViirsGeoFileObj.end()
                except pyHDF.HDF4Error:
                    print "Could not get geolocation dataset"
                    ViirsGeoFileObj.end()
                    if __name__=='__main__':
                        sys.exit(1)
                    else :
                        return 1

        except pyEx.HDF5ExtError :
            print "Could not determine whether geolocation isHDF4 or HDF5, aborting..."
            if __name__=='__main__':
                sys.exit(1)
            else :
                return 1

		# Remove any unit-length dimensions from geolocation arrays
		self.Lon = np.squeeze(self.Lon)
		self.Lat = np.squeeze(self.Lat)
		print "Shape of Latitude:",shape(self.Lat)
		print "Shape of Longitude:",shape(self.Lon)

        # Open Viirs Cloud Mask file
        try :
            ViirsCMaskFileObj = pytables.openFile(ViirsCMaskFileName,'r')
            print "Successfully opened HDF5 Cloud Mask file ",ViirsCMaskFileName
        except :
            print "Could not open Cloud Mask file: ",ViirsCMaskFileName
            if __name__=='__main__':
                sys.exit(1)
            else :
                return 1

        # Get the VIIRS VCM quantities (Phase and Quality)
        try :
            ViirsDataSetPath = '/All_Data/VIIRS-CM-IP_All/QF1_VIIRSCMIP'
            self.ViirsCMquality = np.copy(ViirsCMaskFileObj.getNode(ViirsDataSetPath))[::shrinkFactor,::shrinkFactor]
            self.ViirsCMquality = np.bitwise_and(self.ViirsCMquality,3) >> 0
            self.ViirsCMqualityMask = ma.masked_equal(self.ViirsCMquality,0)

            ViirsDataSetPath = '/All_Data/VIIRS-CM-IP_All/QF6_VIIRSCMIP'
            self.ViirsCMphase = np.copy(ViirsCMaskFileObj.getNode(ViirsDataSetPath))[::shrinkFactor,::shrinkFactor]
            self.ViirsCMphase = np.bitwise_and(self.ViirsCMphase,7) >> 0

            print "Closing Cloud Mask file"
            ViirsCMaskFileObj.close()

        except pyEx.NoSuchNodeError :
            print "There is no ",ViirsDataSetPath," SDR in",ViirsCMaskFileName
            ViirsCMaskFileObj.close()
            if __name__=='__main__':
                sys.exit(1)
            else :
                return 1

        # Open the VIIRS Cloud IP file
        try :
            ViirsCProdFileObj = pytables.openFile(ViirsCProdFileName,mode='r')
            print "Successfully opened IP file"
        except IOError :
            print "Could not open IP file: ",ViirsCProdFileName
            if __name__=='__main__':
                sys.exit(1)
            else :
                return 1

        # Get the VIIRS IP dataset
        try :
            self.ViirsCProdSDS = np.copy(ViirsCProdFileObj.getNode(CloudProduct.SDSname))[::shrinkFactor,::shrinkFactor]

            print "Closing Cloud Product file"
            ViirsCProdFileObj.close()

        except pyEx.NoSuchNodeError :
            print "There is no ",CloudProduct.SDSname," SDS in",ViirsCProdFileName
            ViirsCProdFileObj.close()
            if __name__=='__main__':
                sys.exit(1)
            else :
                return 1

        # Make logarithmic if appropriate
        if(CloudProduct.logScale):
            print "Log scaling the dataset..."
            self.ViirsCProdSDS = np.log10(self.ViirsCProdSDS)

