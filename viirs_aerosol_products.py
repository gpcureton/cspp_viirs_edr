#!/usr/bin/env python
# encoding: utf-8
"""
viirs_aerosol_products.py

Purpose: Provide ingest functionality for VIIRS aerosol IP and EDR.

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

from matplotlib import pyplot as ppl
from matplotlib.colors import ListedColormap
from matplotlib import cm as cm

from mpl_toolkits.basemap import Basemap,shiftgrid

import numpy as np
from numpy import ma as ma

import tables as pytables
from tables import exceptions as pyEx

### Local libraries
import ViirsData as VD

class viirsAero:

    def __init__(self):
        '''
        This __init__ method sets static data for the class.
        '''
        pass

    def ingest(self,geoFile,aeroFile,SDS,shrink,scale):
        '''
        This method ingests the lat, long and the specified SDS from the geolocation
        and VIIRS Aerosol IP files.
        '''

        # Get the input file names
        ViirsGeoFileName       = str(geoFile )
        ViirsAProdFileName     = str(aeroFile )
        ViirsAProdSDStype      = str(SDS     )
        shrinkFactor           = int(shrink  )
        scale                  = str(scale   )

        self.ViirsGeoFileName       = ViirsGeoFileName
        self.ViirsAProdFileName  = ViirsAProdFileName
        self.ViirsAProdSDStype = SDS
        self.scale = scale

        # Get the various parameters related to the current product
        AerosolProduct = VD.AerosolProdData.AerosolProduct[ViirsAProdSDStype]
        if isinstance(AerosolProduct,VD.AerosolProdData.AerosolProd) :
            pass
        else :
            AerosolProduct = VD.AerosolProdData.AerosolProduct[ViirsAProdSDStype][scale] # aot

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
                    
                    dataName = 'SolarZenithAngle'
                    print "Reading %s dataset..." % (dataName)
                    geoNode = ViirsGeoFileObj.getNode('%s/%s' % (geoGroupName,dataName))
                    self.Sza = geoNode[::shrinkFactor,::shrinkFactor]
                    geoNode.close()
                    print "done"
                    
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
        print "Shape of Solar Zenith Angle:",np.shape(self.Sza)

        # Open the VIIRS Aerosol IP file
        try :
            ViirsAProdFileObj = pytables.openFile(ViirsAProdFileName,mode='r')
            print "Successfully opened IP file",ViirsAProdFileName
        except IOError :
            print "Could not open IP file: ",ViirsAProdFileName
            sys.exit(1)

        # Get the VIIRS IP dataset
        try :
            self.ViirsAProdSDS = getobj(ViirsAProdFileObj,\
                AerosolProduct.SDSname)[::shrinkFactor,::shrinkFactor]
            print "Type of self.ViirsAProdSDS is: ",self.ViirsAProdSDS.dtype.name

            ViirsDataSetPath = '/All_Data/VIIRS-Aeros-Opt-Thick-IP_All/QF1'
            self.ViirsCMquality =  getobj(ViirsAProdFileObj,\
                ViirsDataSetPath)[::shrinkFactor,::shrinkFactor]
            self.ViirsCMquality = np.bitwise_and(self.ViirsCMquality,192) >> 6

            ViirsDataSetPath = '/All_Data/VIIRS-Aeros-Opt-Thick-IP_All/QF2'
            self.LandSeaMask =  getobj(ViirsAProdFileObj,\
                ViirsDataSetPath)[::shrinkFactor,::shrinkFactor]
            self.LandSeaMask = np.bitwise_and(self.LandSeaMask,112) >> 4

            ViirsDataSetPath = '/All_Data/VIIRS-Aeros-Opt-Thick-IP_All/QF3'
            self.ViirsAProdRet =  getobj(ViirsAProdFileObj,\
                ViirsDataSetPath)[::shrinkFactor,::shrinkFactor]
            self.ViirsAProdRet = np.bitwise_and(self.ViirsAProdRet,28) >> 2

            print "Closing Aerosol Product file"
            ViirsAProdFileObj.close()

        except pyEx.NoSuchNodeError :
            print "error: There is no ",AerosolProduct.SDSname," dataset in",ViirsAProdFileName
            ViirsAProdFileObj.close()
            sys.exit(1)

        return 0

    def ingestEDR(self,geoFile,aeroFile,shrink,scale):
        '''
        This method ingests the lat, long and the specified SDS from the geolocation
        and VIIRS Aerosol EDR files.
        '''

        # Get the input file names
        ViirsGeoFileName       = str(geoFile )
        ViirsAProdFileName     = str(aeroFile )
        #ViirsAProdSDStype      = str(SDS     )
        shrinkFactor           = int(shrink  )
        scale                  = str(scale   )

        self.ViirsGeoFileName       = ViirsGeoFileName
        self.ViirsAProdFileName  = ViirsAProdFileName
        #self.ViirsAProdSDStype = SDS
        self.scale = scale

        # Determine the correct fillValue
        trimObj = VD.ViirsTrimTable()
        eps = 1.e-6

        # Read in geolocation...
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
            print "Geolocation Group : %s " % (geoGroupName)
            isEdrGeo = ('VIIRS-Aeros-EDR-GEO_All' in geoGroupName)
            if not isEdrGeo :
                print ">> error: %s is not an EDR resolution aerosol geolocation file\n\taborting..."  % (ViirsGeoFileName)
                sys.exit(1)

            # Get the geolocation Latitude and Longitude nodes...
            try :
                dataName = 'Latitude'
                geoLatNode = ViirsGeoFileObj.getNode(geoGroupName+'/'+dataName)
                print "Shape of %s node is %s" % (geoGroupName+'/'+dataName,repr(geoLatNode.shape))
                dataName = 'Longitude'
                geoLonNode = ViirsGeoFileObj.getNode(geoGroupName+'/'+dataName)
                print "Shape of %s node is %s" % (geoGroupName+'/'+dataName,repr(geoLonNode.shape))
            except pyEx.NoSuchNodeError :
                print "\n>> error: No required node %s/%s in %s\n\taborting..." % (geoGroupName,dataName,ViirsGeoFileName)
                ViirsGeoFileObj.close()
                sys.exit(1)

            # Make a copy of the geolocation node data, so we can close the file
            try :
                dataName = 'Latitude'
                print "Reading %s dataset..." % (dataName)
                latArr = geoLatNode[:,:]
                geoLatNode.close()
                latArr = np.squeeze(latArr)
                print "done"
                dataName = 'Longitude'
                print "Reading %s dataset..." % (dataName)
                lonArr = geoLonNode[:,:]
                geoLonNode.close()
                lonArr = np.squeeze(lonArr)
                print "done"
                print "Closing geolocation file"
                ViirsGeoFileObj.close()
            except :
                print "\n>> error: Could not retrieve %/% node data in %s\n\taborting..." % (geoGroupName,dataName,ViirsGeoFileName)
                geoLatNode.close()
                geoLonNode.close()
                ViirsGeoFileObj.close()
                sys.exit(1)

        else :
            print "\n>> error: %s is not a HDF5 file,\n\taborting..." % (ViirsGeoFileName)
            sys.exit(1)
        
        # Read in dataSets...
        ViirsEDRFileName = ViirsAProdFileName
        if(pytables.isHDF5File(ViirsEDRFileName)) :
            try :
                ViirsEDRFileObj = pytables.openFile(ViirsEDRFileName,mode='r')
                print "Successfully opened edr file",ViirsEDRFileName
            except IOError :
                print ">> error: Could not open edr file: ",ViirsEDRFileName
                sys.exit(1)

            # Detemine the edr group name and related information
            group = ViirsEDRFileObj.getNode('/All_Data')
            edrGroupName = '/All_Data/'+group.__members__[0]
            group._g_close()
            print "Edr Group : %s " % (edrGroupName)
            isEdr = ('VIIRS-Aeros-EDR_All' in edrGroupName)
            if not isEdr :
                print ">> error: %s is not an EDR resolution aerosol file\n\taborting..."  % (ViirsEDRFileName)
                sys.exit(1)

            # Get the edr nodes...
            try :
                dataName = 'AerosolOpticalDepth_at_550nm'
                aot550Node = ViirsEDRFileObj.getNode(edrGroupName+'/'+dataName)
                print "Shape of %s node is %s" % (edrGroupName+'/'+dataName,repr(aot550Node.shape))
                dataName = 'AerosolOpticalDepthFactors'
                aotFactorsNode = ViirsEDRFileObj.getNode(edrGroupName+'/'+dataName)
                dataName = 'QF1_VIIRSAEROEDR'
                qf1Node = ViirsEDRFileObj.getNode(edrGroupName+'/'+dataName)
                print "Shape of %s node is %s" % (edrGroupName+'/'+dataName,repr(aotFactorsNode.shape))
            except pyEx.NoSuchNodeError :
                print "\n>> error: No required node %s/%s in %s\n\taborting..." % (edrGroupName,dataName,ViirsEDRFileName)
                ViirsEDRFileObj.close()
                sys.exit(1)

            # Make a copy of the edr node data, so we can close the file
            try :
                dataName = 'AerosolOpticalDepth_at_550nm'
                print "Reading %s dataset..." % (dataName)
                aot550 = aot550Node[:,:]
                aot550 = np.squeeze(aot550)
                aot550Node.close()
                print "done"
                dataName = 'AerosolOpticalDepthFactors'
                print "Reading %s dataset..." % (dataName)
                aotFactors = aotFactorsNode[:]
                aotFactors = np.squeeze(aotFactors)
                aotFactorsNode.close()
                print "done"
                dataName = 'QF1_VIIRSAEROEDR'
                print "Reading %s dataset..." % (dataName)
                qf1 = qf1Node[:]
                qf1 = np.squeeze(qf1Node)
                qf1Node.close()
                print "done"
                print "Closing edr file"
                ViirsEDRFileObj.close()

                print "Shape of aot550 is %s" % (repr(np.shape(aot550)))
                print "Shape of qf1 is %s" % (repr(np.shape(qf1)))
                print "Shape of latsArr is %s" % (repr(np.shape(latArr)))
                print "Shape of lonsArr is %s" % (repr(np.shape(lonArr)))

            except :
                print "\n>> error: Could not retrieve %/% node data in %s\n\taborting..." % (edrGroupName,dataName,ViirsEDRFileName)
                aot550Node.close()
                aotFactorsNode.close()
                qf1Node.close()
                ViirsEDRFileObj.close()
                sys.exit(1)

        else :
            print "\n>> error: %s is not a HDF5 file,\n\taborting..." % (ViirsEDRFileName)
            sys.exit(1)
        
        print "Creating some masks"
        try :
            aotArr  = aot550
            qf1Arr  = qf1

            # Determine masks for each fill type, for the AOT EDR
            aotFillMasks = {}
            for fillType in trimObj.sdrTypeFill.keys() :
                fillValue = trimObj.sdrTypeFill[fillType][aotArr.dtype.name]
                if 'float' in fillValue.__class__.__name__ :
                    aotFillMasks[fillType] = ma.masked_inside(aotArr,fillValue-eps,fillValue+eps).mask
                    if (aotFillMasks[fillType].__class__.__name__ != 'ndarray') :
                        aotFillMasks[fillType] = None
                elif 'int' in fillValue.__class__.__name__ :
                    aotFillMasks[fillType] = ma.masked_equal(aotArr,fillValue).mask
                    if (aotFillMasks[fillType].__class__.__name__ != 'ndarray') :
                        aotFillMasks[fillType] = None
                else :
                    print "Dataset was neither int not float... a worry"
                    pass

            # Construct the total mask from all of the various fill values
            totalMask = ma.array(np.zeros(aotArr.shape,dtype=np.bool))
            for fillType in trimObj.sdrTypeFill.keys() :
                if aotFillMasks[fillType] is not None :
                    totalMask = totalMask * ma.array(np.zeros(aotArr.shape,dtype=np.bool),\
                        mask=aotFillMasks[fillType])


            # Apply factors to those pixels which are not masked
            aotArr = aotArr * aotFactors[0] + aotFactors[1]
            # Extract the AOT Product Quality QF array
            qf1Arr = np.bitwise_and(qf1Arr,3) >> 0

        except :
            print ">> error: There was an exception..."
            sys.exit(1)

        self.Lat  = latArr
        self.Lon  = lonArr
        self.aotArr = aotArr
        self.qf1Arr = qf1Arr
            
        return 0
        #return lats,lons,data,lat_0,lon_0
