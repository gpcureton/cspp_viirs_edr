#!/usr/bin/env python
# encoding: utf-8
"""
ql_viirs_edr.py

Purpose: Create quicklook PNGs for VIIRS EDR products.

Minimum commandline...

export CSPP_HOME=/path/to/CSPP
source $CSPP_HOME/cspp_edr_env.sh
source $CSPP_HOME/common/cspp_common.sh

python ql_viirs_edr.py -g geofile.h5 -i ipfile.h5 -p CTP

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

import string, sys
from glob import glob
from os import path,uname
from time import time

import numpy as np
from  numpy import ma as ma
import scipy as scipy

import matplotlib
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap
from matplotlib.figure import Figure

matplotlib.use('Agg')
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

# This must come *after* the backend is specified.
import matplotlib.pyplot as ppl

from mpl_toolkits.basemap import Basemap

import optparse as optparse

#import VIIRS as VIIRS
import ViirsData
import viirs_cloud_mask as viirsCM
import viirs_cloud_products as viirsCld
import viirs_aerosol_products as viirsAero

import tables as pytables
from tables import exceptions as pyEx

###################################################
#                  Global Data                    #
###################################################

# VCM dataset selection
#global cmByte,cmBit 

def set_vcm_dset(newCmByte,newCmBit) :
    global cmByte,cmBit 
    cmByte,cmBit = newCmByte,newCmBit

###################################
#          Get File Lists         #
###################################

def granuleFiles(geoGlob,prodGlob):
    '''
    Returns sorted lists of the geolocation and product files.
    '''
    
    geoDir = path.dirname(path.abspath(path.expanduser(geoGlob)))
    prodDir = path.dirname(path.abspath(path.expanduser(prodGlob)))

    print "Initial geoGlob = ",geoGlob
    print "Initial prodGlob = ",prodGlob

    geoGlob = path.basename(path.abspath(path.expanduser(geoGlob)))
    prodGlob = path.basename(path.abspath(path.expanduser(prodGlob)))

    geoPrefix = string.split(geoGlob,'_')[0]

    print "geoDir = ",geoDir
    print "prodDir = ",prodDir
    print "geoGlob = ",geoGlob
    print "prodGlob = ",prodGlob
    print "geoPrefix = ",geoPrefix

    geoList_in = glob("%s/%s" % (geoDir,geoGlob))
    prodList_in = glob("%s/%s" % (prodDir,prodGlob))
    geoList_in.sort()
    prodList_in.sort()
    
    #print "prodList_in..."
    #for prodFile in prodList_in:
        #print prodFile

    geoList = []
    prodList = []
    #prodList = prodList_in
    for files in prodList_in :
        prod_arr = string.split(path.basename(files),"_")
        #print "prod_arr = ",prod_arr
        dateStamp = prod_arr[2]
        timeStamp = prod_arr[3]
        geoFileGlob="%s/%s*%s_%s*.h5" % (geoDir,geoPrefix,dateStamp,timeStamp)
        #print "dateStamp = ",dateStamp
        #print "timeStamp = ",timeStamp
        #print "geoFileGlob = ",geoFileGlob
        geoFile = glob("%s/%s*%s_%s*.h5" % (geoDir,geoPrefix,dateStamp,timeStamp))
        if (np.shape(geoFile)[0] != 0) :
            geoList.append(geoFile[0])
            prodList.append(files)
        else :
            #geoList.append(files)
            print " ... no match found for %s, appending %s" % ( geoFile, files)
            pass
    
    #for geoFile,prodFile in zip(geoList,prodList):
        #print geoFile,prodFile
    return geoList,prodList


###################################################
#              Granulation Functions              #
###################################################

def gran_VCM(geoList,cmList,shrink=1):
    '''
    Returns the granulated VCM
    '''

    global cmByte,cmBit 

    # Define the colourmap we want to use.
    print "Importing the CloudMaskProduct object..."
    CloudMaskProduct = ViirsData.CloudMaskData
    print "done"

    try :
        reload(viirsCM)
        reload(ViirsData)
        del(viirsCMobj)
        del(latArr)
        del(lonArr)
        del(cmArr)
        del(qualArr)
    except :
        pass

    print "Creating viirsCMobj..."
    viirsCMobj = viirsCM.viirsCM()
    print "done"

    # Determine the correct fillValue
    trimObj = ViirsData.ViirsTrimTable()
    eps = 1.e-6

    for grans in np.arange(len(geoList)):

        print "\nIngesting granule %d using cmByte=%d and cmBit=%d ..." % (grans,cmByte,cmBit)
        retArr = viirsCMobj.ingest(geoList[grans],cmList[grans],cmByte,cmBit,shrink)

        try :
            latArr  = np.vstack((latArr,viirsCMobj.Lat[:,:]))
            lonArr  = np.vstack((lonArr,viirsCMobj.Lon[:,:]))
        except NameError :
            latArr  = viirsCMobj.Lat[:,:]
            lonArr  = viirsCMobj.Lon[:,:]
        
        try :
            cmArr  = np.vstack((cmArr,viirsCMobj.ViirsCMaskSDS[:,:]))
            qualArr  = np.vstack((qualArr,viirsCMobj.ViirsCMquality[:,:]))
        except NameError :
            cmArr   = viirsCMobj.ViirsCMaskSDS[:,:]
            qualArr = viirsCMobj.ViirsCMquality[:,:]

        print "Intermediate latArr.shape = %s" % (str(latArr.shape))
        print "Intermediate lonArr.shape = %s" % (str(lonArr.shape))
        print "Intermediate cmArr.shape = %s" % (str(cmArr.shape))
        print "Intermediate qualArr.shape = %s" % (str(qualArr.shape))

    lat_0 = latArr[np.shape(latArr)[0]/2,np.shape(latArr)[1]/2]
    lon_0 = lonArr[np.shape(lonArr)[0]/2,np.shape(lonArr)[1]/2]

    try :
        # Determine masks for each fill type, for the VCM IP
        cmFillMasks = {}
        for fillType in trimObj.sdrTypeFill.keys() :
            fillValue = trimObj.sdrTypeFill[fillType][cmArr.dtype.name]
            if 'float' in fillValue.__class__.__name__ :
                cmFillMasks[fillType] = ma.masked_inside(cmArr,fillValue-eps,fillValue+eps).mask
                if (cmFillMasks[fillType].__class__.__name__ != 'ndarray') :
                    cmFillMasks[fillType] = None
            elif 'int' in fillValue.__class__.__name__ :
                cmFillMasks[fillType] = ma.masked_equal(cmArr,fillValue).mask
                if (cmFillMasks[fillType].__class__.__name__ != 'ndarray') :
                    cmFillMasks[fillType] = None
            else :
                print "Dataset was neither int not float... a worry"
                pass

        # Construct the total mask from all of the various fill values
        totalMask = ma.array(np.zeros(cmArr.shape,dtype=np.bool))
        for fillType in trimObj.sdrTypeFill.keys() :
            if cmFillMasks[fillType] is not None :
                totalMask = totalMask * ma.array(np.zeros(cmArr.shape,dtype=np.bool),\
                    mask=cmFillMasks[fillType])

        # Define the masks, and mask the main dataset
        ViirsCMqualityMask = ma.masked_equal(qualArr,0)
        ViirsCMclearMask = ma.masked_equal(cmArr,1)

        if cmByte==0 and cmBit==1 :
            totalMask = totalMask * ViirsCMqualityMask
        elif cmByte==5 and cmBit==0 :
            totalMask = totalMask * ViirsCMqualityMask * ViirsCMclearMask

        try :
            data = ma.array(cmArr,mask=totalMask.mask)
            lats = ma.array(latArr,mask=totalMask.mask)
            lons = ma.array(lonArr,mask=totalMask.mask)
        except ma.core.MaskError :
            print ">> error: Mask Error, probably mismatched geolocation and product array sizes, aborting..."
            sys.exit(1)

    except Exception, err :
        print ">> error: %s..." % (str(err))
        sys.exit(1)

    return lats,lons,data,lat_0,lon_0

def gran_COP(geoList,copList,dataSet,shrink=1):
    '''
    Returns the granulated COP dataset
    '''

    try :
        reload(viirsCld)
        reload(ViirsData)
        del(viirsCldObj)
        del(latArr)
        del(lonArr)
        del(copArr)
        del(phaseArr)
    except :
        pass

    print "Creating viirsCldObj..."
    reload(viirsCld)
    viirsCldObj = viirsCld.viirsCld()
    print "done"

    # Determine the correct fillValue
    trimObj = ViirsData.ViirsTrimTable()
    eps = 1.e-6

    for grans in np.arange(len(geoList)):
        print "Ingesting granule %d ..." % (grans)
        retArr = viirsCldObj.ingestLite(geoList[grans],copList[grans],dataSet,1)
        print "done"

        try :
            latArr  = viirsCldObj.Lat[:,:]
            lonArr  = viirsCldObj.Lon[:,:]
            
            lat_0 = latArr[np.shape(latArr)[0]/2,np.shape(latArr)[1]/2]
            lon_0 = lonArr[np.shape(lonArr)[0]/2,np.shape(lonArr)[1]/2]

            badGeo = False
            if not (-90. <= lat_0 <= 90.) :
                print "\n>> error: Latitude of granule midpoint (%f) does not satisfy (-90. <= lat_0 <= 90.)\nfor file %s\n\taborting..." % (lat_0,geoList[grans])
                badGeo = True
            if not (-180. <= lat_0 <= 180.) :
                print "\n>> error: Longitude of granule midpoint (%f) does not satisfy (-180. <= lon_0 <= 180.)\nfor file %s\n\taborting..." % (lon_0,geoList[grans])
                badGeo = True

            if badGeo :
                sys.exit(1)

            copArr   = viirsCldObj.ViirsCProdSDS[:,:]
            phaseArr = viirsCldObj.ViirsCProdSDSphase[:,:]


            # Determine masks for each fill type, for the COP IP
            copFillMasks = {}
            for fillType in trimObj.sdrTypeFill.keys() :
                fillValue = trimObj.sdrTypeFill[fillType][copArr.dtype.name]
                if 'float' in fillValue.__class__.__name__ :
                    copFillMasks[fillType] = ma.masked_inside(copArr,fillValue-eps,fillValue+eps).mask
                    if (copFillMasks[fillType].__class__.__name__ != 'ndarray') :
                        copFillMasks[fillType] = None
                elif 'int' in fillValue.__class__.__name__ :
                    copFillMasks[fillType] = ma.masked_equal(copArr,fillValue).mask
                    if (copFillMasks[fillType].__class__.__name__ != 'ndarray') :
                        copFillMasks[fillType] = None
                else :
                    print "Dataset was neither int not float... a worry"
                    pass

            # Construct the total mask from all of the various fill values
            totalMask = ma.array(np.zeros(copArr.shape,dtype=np.bool))
            for fillType in trimObj.sdrTypeFill.keys() :
                if copFillMasks[fillType] is not None :
                    totalMask = totalMask * ma.array(np.zeros(copArr.shape,dtype=np.bool),\
                        mask=copFillMasks[fillType])

            # Added the "bad data" from the COP phase quality flag.
            totalMask = totalMask * ma.masked_equal(phaseArr,0)

            try :
                data  = ma.array(copArr,  mask=totalMask.mask)
                phase = ma.array(phaseArr,mask=totalMask.mask)
                lats  = ma.array(latArr,  mask=totalMask.mask)
                lons  = ma.array(lonArr,  mask=totalMask.mask)
            except ma.core.MaskError :
                print ">> error: Mask Error, probably mismatched geolocation and product array sizes, aborting..."
                sys.exit(1)

        except :
            print ">> error: there was an exception..."
            sys.exit(1)

    return lats,lons,data,phase,lat_0,lon_0

def gran_COT_EDR(geoList,cotList,shrink=1):
    '''
    Returns the granulated COT EDR
    '''

    try :
        reload(viirsAero)
        reload(ViirsData)
        del(viirsAeroObj)
        del(latArr)
        del(lonArr)
        del(cotArr)
    except :
        pass

    print "Creating viirsCldObj..."
    reload(viirsCld)
    viirsCldObj = viirsCld.viirsCld()
    print "done"

    # Determine the correct fillValue
    trimObj = ViirsData.ViirsTrimTable()
    eps = 1.e-6

    for grans in np.arange(len(geoList)):
        print "Ingesting EDR granule %d ..." % (grans)

        # Read in geolocation...
        ViirsGeoFileName = geoList[grans]
        if(pytables.isHDF5File(ViirsGeoFileName)) :
            try :
                ViirsGeoFileObj = pytables.openFile(ViirsGeoFileName,mode='r')
                print "Successfully opened geolocation file",ViirsGeoFileName
            except IOError :
                print ">> error: Could not open geolocation file: ",ViirsGeoFileName
                sys.exit(1)

            # Detemine the geolocation group name and related information
            #group = getobj(ViirsGeoFileObj,'/All_Data')
            group = ViirsGeoFileObj.getNode('/All_Data')
            geoGroupName = '/All_Data/'+group.__members__[0]
            group._g_close()
            print "Geolocation Group : %s " % (geoGroupName)
            isEdrGeo = ('VIIRS-CLD-AGG-GEO_All' in geoGroupName)
            if not isEdrGeo :
                print ">> error: %s is not an EDR resolution cloud geolocation file\n\taborting..."  % (ViirsGeoFileName)
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
        ViirsEDRFileName = cotList[grans]
        if(pytables.isHDF5File(ViirsEDRFileName)) :
            try :
                ViirsEDRFileObj = pytables.openFile(ViirsEDRFileName,mode='r')
                print "Successfully opened edr file",ViirsEDRFileName
            except IOError :
                print ">> error: Could not open edr file: ",ViirsEDRFileName
                sys.exit(1)

            # Detemine the edr group name and related information
            #group = getobj(ViirsEDRFileObj,'/All_Data')
            group = ViirsEDRFileObj.getNode('/All_Data')
            edrGroupName = '/All_Data/'+group.__members__[0]
            group._g_close()
            print "Edr Group : %s " % (edrGroupName)
            isEdr = ('VIIRS-COT-EDR_All' in edrGroupName)
            if not isEdr :
                print ">> error: %s is not an EDR resolution cloud file\n\taborting..."  % (ViirsEDRFileName)
                sys.exit(1)

            # Get the edr nodes...
            try :
                dataName = 'AverageCloudOpticalThickness'
                cotNode = ViirsEDRFileObj.getNode(edrGroupName+'/'+dataName)
                print "Shape of %s node is %s" % (edrGroupName+'/'+dataName,repr(cotNode.shape))
                dataName = 'COTFactors'
                cotFactorsNode = ViirsEDRFileObj.getNode(edrGroupName+'/'+dataName)
                print "Shape of %s node is %s" % (edrGroupName+'/'+dataName,repr(cotFactorsNode.shape))
                dataName = 'QF3_VIIRSCOTAVGEDR'
                QF3_VIIRSCOTAVGEDR_Node = ViirsEDRFileObj.getNode(edrGroupName+'/'+dataName)
                print "Shape of %s node is %s" % (edrGroupName+'/'+dataName,repr(QF3_VIIRSCOTAVGEDR_Node.shape))
            except pyEx.NoSuchNodeError :
                print "\n>> error: No required node %s/%s in %s\n\taborting..." % (edrGroupName,dataName,ViirsEDRFileName)
                ViirsEDRFileObj.close()
                sys.exit(1)

            # Make a copy of the edr node data, so we can close the file
            try :
                dataName = 'AverageCloudOpticalThickness'
                print "Reading %s dataset..." % (dataName)
                cot = cotNode[:,:]
                cot = np.squeeze(cot)
                cotNode.close()
                print "done"
                dataName = 'COTFactors'
                print "Reading %s dataset..." % (dataName)
                cotFactors = cotFactorsNode[:]
                cotFactors = np.squeeze(cotFactors)
                cotFactorsNode.close()
                print "done"
                dataName = 'QF3_VIIRSCOTAVGEDR'
                print "Reading %s dataset..." % (dataName)
                cloudConf = np.bitwise_and(QF3_VIIRSCOTAVGEDR_Node[:,:],3) >> 0
                waterCldFrac = np.bitwise_and(QF3_VIIRSCOTAVGEDR_Node[:,:],12) >> 2
                mixedCldFrac = np.bitwise_and(QF3_VIIRSCOTAVGEDR_Node[:,:],192) >> 6
                QF3_VIIRSCOTAVGEDR_Node.close()
                print "done"
                print "Closing edr file"
                ViirsEDRFileObj.close()

                print "Shape of cot is %s" % (repr(np.shape(cot)))
                print "Shape of latsArr is %s" % (repr(np.shape(latArr)))
                print "Shape of lonsArr is %s" % (repr(np.shape(lonArr)))

            except :
                print "\n>> error: Could not retrieve %/% node data in %s\n\taborting..." % (edrGroupName,dataName,ViirsEDRFileName)
                cotNode.close()
                cotFactorsNode.close()
                QF3_VIIRSCOTAVGEDR_Node.close()
                ViirsEDRFileObj.close()
                sys.exit(1)

        else :
            print "\n>> error: %s is not a HDF5 file,\n\taborting..." % (ViirsEDRFileName)
            sys.exit(1)
        
        print "Creating some masks"
        try :
            
            lat_0 = latArr[np.shape(latArr)[0]/2,np.shape(latArr)[1]/2]
            lon_0 = lonArr[np.shape(lonArr)[0]/2,np.shape(lonArr)[1]/2]

            badGeo = False
            if not (-90. <= lat_0 <= 90.) :
                print "\n>> error: Latitude of granule midpoint (%f) does not satisfy (-90. <= lat_0 <= 90.)\nfor file %s\n\taborting..." % (lat_0,geoList[grans])
                badGeo = True
            if not (-180. <= lat_0 <= 180.) :
                print "\n>> error: Longitude of granule midpoint (%f) does not satisfy (-180. <= lon_0 <= 180.)\nfor file %s\n\taborting..." % (lon_0,geoList[grans])
                badGeo = True

            if badGeo :
                sys.exit(1)

            cotArr     = cot[:,:]

            cloudConfArr    = cloudConf[:,:]
            waterCldFracArr = waterCldFrac[:,:]
            mixedCldFracArr = mixedCldFrac[:,:]

            # Determine masks for each fill type, for the EPS EDR
            cotFillMasks = {}
            for fillType in trimObj.sdrTypeFill.keys() :
                fillValue = trimObj.sdrTypeFill[fillType][cotArr.dtype.name]
                if 'float' in fillValue.__class__.__name__ :
                    cotFillMasks[fillType] = ma.masked_inside(cotArr,fillValue-eps,fillValue+eps).mask
                    if (cotFillMasks[fillType].__class__.__name__ != 'ndarray') :
                        cotFillMasks[fillType] = None
                elif 'int' in fillValue.__class__.__name__ :
                    cotFillMasks[fillType] = ma.masked_equal(cotArr,fillValue).mask
                    if (cotFillMasks[fillType].__class__.__name__ != 'ndarray') :
                        cotFillMasks[fillType] = None
                else :
                    print "Dataset was neither int not float... a worry"
                    pass

            # Construct the total mask from all of the various fill values
            totalMask = ma.array(np.zeros(cotArr.shape,dtype=np.bool))
            for fillType in trimObj.sdrTypeFill.keys() :
                if cotFillMasks[fillType] is not None :
                    totalMask = totalMask * ma.array(np.zeros(cotArr.shape,dtype=np.bool),\
                        mask=cotFillMasks[fillType])

            try :
                data    = ma.array(cotArr,    mask=totalMask.mask)
                cloudConfArr    = ma.array(cloudConfArr   ,mask=totalMask.mask)
                waterCldFracArr = ma.array(waterCldFracArr,mask=totalMask.mask)
                mixedCldFracArr = ma.array(mixedCldFracArr,mask=totalMask.mask)
                lats    = ma.array(latArr,    mask=totalMask.mask)
                lons    = ma.array(lonArr,    mask=totalMask.mask)
            except ma.core.MaskError :
                print ">> error: Mask Error, probably mismatched geolocation and product array sizes, aborting..."
                sys.exit(1)

            data = data * cotFactors[0] + cotFactors[1]

            cldPhase = np.ones(data.shape,dtype=np.uint8) * trimObj.sdrTypeFill['NA_FILL']['uint8']

            print "cloudPhase.shape = ",cldPhase.shape
            print "cloudConfArr.shape = ",cloudConfArr.shape
            print "waterCldFracArr.shape = ",waterCldFracArr.shape
            print "mixedCldFracArr.shape = ",mixedCldFracArr.shape

            cloudFracMask = (cloudConfArr>0)
            waterFracMask = (waterCldFracArr>0)
            mixedFracMask = (mixedCldFracArr>0)

            waterMask = cloudFracMask * (waterFracMask + mixedFracMask)
            iceMask = cloudFracMask * np.logical_not(waterFracMask + mixedFracMask)

            print "waterMask.shape = ",waterMask.shape
            print "iceMask.shape = ",iceMask.shape

            print "Number of pixels with > 25% cloudy confidence: ",cloudFracMask.sum()
            print "Number of pixels with > 25% water cloud fraction: ",waterFracMask.sum()
            print "Number of pixels with > 25% mixed cloud fraction: ",mixedFracMask.sum()

            print "Number of water pixels: ",(cloudFracMask * (waterFracMask + mixedFracMask)).sum()
            print "Number of ice pixels: ",cloudFracMask.sum() - (waterFracMask + mixedFracMask).sum()

            print "Number of water pixels: ",waterMask.sum()
            print "Number of ice pixels: ",iceMask.sum()

            cldPhaseShape = cldPhase.shape
            cldPhase = np.ravel(cldPhase)

            wtrIdx = np.where(np.ravel(waterMask))
            iceIdx = np.where(np.ravel(iceMask))
            cldPhase[wtrIdx] = 3
            cldPhase[iceIdx] = 5

            cldPhase = np.reshape(cldPhase,cldPhaseShape)

        except Exception, err :
            print ">> error: %s..." % (str(err))
            sys.exit(1)

    print "lat_0,lon_0 = %f,%f" % (lat_0,lon_0)
    print "Shape of data is %s" % (repr(np.shape(data)))
    print "Shape of cldPhase is %s" % (repr(np.shape(cldPhase)))
    print "Shape of lats is %s" % (repr(np.shape(lats)))
    print "Shape of lons is %s" % (repr(np.shape(lons)))

    return lats,lons,data,cldPhase,lat_0,lon_0

def gran_EPS_EDR(geoList,epsList,shrink=1):
    '''
    Returns the granulated EPS EDR
    '''

    try :
        reload(viirsCld)
        reload(ViirsData)
        del(viirsCldObj)
        del(latArr)
        del(lonArr)
        del(epsArr)
    except :
        pass

    print "Creating viirsCldObj..."
    reload(viirsCld)
    viirsCldObj = viirsCld.viirsCld()
    print "done"

    # Determine the correct fillValue
    trimObj = ViirsData.ViirsTrimTable()
    eps = 1.e-6

    for grans in np.arange(len(geoList)):
        print "Ingesting EDR granule %d ..." % (grans)

        # Read in geolocation...
        ViirsGeoFileName = geoList[grans]
        if(pytables.isHDF5File(ViirsGeoFileName)) :
            try :
                ViirsGeoFileObj = pytables.openFile(ViirsGeoFileName,mode='r')
                print "Successfully opened geolocation file",ViirsGeoFileName
            except IOError :
                print ">> error: Could not open geolocation file: ",ViirsGeoFileName
                sys.exit(1)

            # Detemine the geolocation group name and related information
            #group = getobj(ViirsGeoFileObj,'/All_Data')
            group = ViirsGeoFileObj.getNode('/All_Data')
            geoGroupName = '/All_Data/'+group.__members__[0]
            group._g_close()
            print "Geolocation Group : %s " % (geoGroupName)
            isEdrGeo = ('VIIRS-CLD-AGG-GEO_All' in geoGroupName)
            if not isEdrGeo :
                print ">> error: %s is not an EDR resolution cloud geolocation file\n\taborting..."  % (ViirsGeoFileName)
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
        ViirsEDRFileName = epsList[grans]
        if(pytables.isHDF5File(ViirsEDRFileName)) :
            try :
                ViirsEDRFileObj = pytables.openFile(ViirsEDRFileName,mode='r')
                print "Successfully opened edr file",ViirsEDRFileName
            except IOError :
                print ">> error: Could not open edr file: ",ViirsEDRFileName
                sys.exit(1)

            # Detemine the edr group name and related information
            #group = getobj(ViirsEDRFileObj,'/All_Data')
            group = ViirsEDRFileObj.getNode('/All_Data')
            edrGroupName = '/All_Data/'+group.__members__[0]
            group._g_close()
            print "Edr Group : %s " % (edrGroupName)
            isEdr = ('VIIRS-CEPS-EDR_All' in edrGroupName)
            if not isEdr :
                print ">> error: %s is not an EDR resolution cloud file\n\taborting..."  % (ViirsEDRFileName)
                sys.exit(1)

            # Get the edr nodes...
            try :
                dataName = 'AverageCloudEffectiveParticleSize'
                epsNode = ViirsEDRFileObj.getNode(edrGroupName+'/'+dataName)
                print "Shape of %s node is %s" % (edrGroupName+'/'+dataName,repr(epsNode.shape))
                dataName = 'CEPSFactors'
                epsFactorsNode = ViirsEDRFileObj.getNode(edrGroupName+'/'+dataName)
                print "Shape of %s node is %s" % (edrGroupName+'/'+dataName,repr(epsFactorsNode.shape))
                dataName = 'QF3_VIIRSCEPSAVGEDR'
                QF3_VIIRSCEPSAVGEDR_Node = ViirsEDRFileObj.getNode(edrGroupName+'/'+dataName)
                print "Shape of %s node is %s" % (edrGroupName+'/'+dataName,repr(QF3_VIIRSCEPSAVGEDR_Node.shape))
            except pyEx.NoSuchNodeError :
                print "\n>> error: No required node %s/%s in %s\n\taborting..." % (edrGroupName,dataName,ViirsEDRFileName)
                ViirsEDRFileObj.close()
                sys.exit(1)

            # Make a copy of the edr node data, so we can close the file
            try :
                dataName = 'AverageCloudEffectiveParticleSize'
                print "Reading %s dataset..." % (dataName)
                eps = epsNode[:,:]
                eps = np.squeeze(eps)
                epsNode.close()
                print "done"
                dataName = 'CEPSFactors'
                print "Reading %s dataset..." % (dataName)
                epsFactors = epsFactorsNode[:]
                epsFactors = np.squeeze(epsFactors)
                epsFactorsNode.close()
                print "done"
                dataName = 'QF3_VIIRSCEPSAVGEDR'
                print "Reading %s dataset..." % (dataName)
                cloudConf = np.bitwise_and(QF3_VIIRSCEPSAVGEDR_Node[:,:],3) >> 0
                waterCldFrac = np.bitwise_and(QF3_VIIRSCEPSAVGEDR_Node[:,:],12) >> 2
                mixedCldFrac = np.bitwise_and(QF3_VIIRSCEPSAVGEDR_Node[:,:],192) >> 6
                QF3_VIIRSCEPSAVGEDR_Node.close()
                print "done"
                print "Closing edr file"
                ViirsEDRFileObj.close()

                print "Shape of eps is %s" % (repr(np.shape(eps)))
                print "Shape of latsArr is %s" % (repr(np.shape(latArr)))
                print "Shape of lonsArr is %s" % (repr(np.shape(lonArr)))

            except :
                print "\n>> error: Could not retrieve %/% node data in %s\n\taborting..." % (edrGroupName,dataName,ViirsEDRFileName)
                epsNode.close()
                epsFactorsNode.close()
                QF3_VIIRSCEPSAVGEDR_Node.close()
                ViirsEDRFileObj.close()
                sys.exit(1)

        else :
            print "\n>> error: %s is not a HDF5 file,\n\taborting..." % (ViirsEDRFileName)
            sys.exit(1)
        
        print "Creating some masks"
        try :
            
            lat_0 = latArr[np.shape(latArr)[0]/2,np.shape(latArr)[1]/2]
            lon_0 = lonArr[np.shape(lonArr)[0]/2,np.shape(lonArr)[1]/2]

            badGeo = False
            if not (-90. <= lat_0 <= 90.) :
                print "\n>> error: Latitude of granule midpoint (%f) does not satisfy (-90. <= lat_0 <= 90.)\nfor file %s\n\taborting..." % (lat_0,geoList[grans])
                badGeo = True
            if not (-180. <= lat_0 <= 180.) :
                print "\n>> error: Longitude of granule midpoint (%f) does not satisfy (-180. <= lon_0 <= 180.)\nfor file %s\n\taborting..." % (lon_0,geoList[grans])
                badGeo = True

            if badGeo :
                sys.exit(1)

            epsArr  = eps[:,:]

            cloudConfArr    = cloudConf[:,:]
            waterCldFracArr = waterCldFrac[:,:]
            mixedCldFracArr = mixedCldFrac[:,:]

            # Determine masks for each fill type, for the EPS EDR
            epsFillMasks = {}
            for fillType in trimObj.sdrTypeFill.keys() :
                fillValue = trimObj.sdrTypeFill[fillType][epsArr.dtype.name]
                if 'float' in fillValue.__class__.__name__ :
                    epsFillMasks[fillType] = ma.masked_inside(epsArr,fillValue-eps,fillValue+eps).mask
                    if (epsFillMasks[fillType].__class__.__name__ != 'ndarray') :
                        epsFillMasks[fillType] = None
                elif 'int' in fillValue.__class__.__name__ :
                    epsFillMasks[fillType] = ma.masked_equal(epsArr,fillValue).mask
                    if (epsFillMasks[fillType].__class__.__name__ != 'ndarray') :
                        epsFillMasks[fillType] = None
                else :
                    print "Dataset was neither int not float... a worry"
                    pass

            # Construct the total mask from all of the various fill values
            totalMask = ma.array(np.zeros(epsArr.shape,dtype=np.bool))
            for fillType in trimObj.sdrTypeFill.keys() :
                if epsFillMasks[fillType] is not None :
                    totalMask = totalMask * ma.array(np.zeros(epsArr.shape,dtype=np.bool),\
                        mask=epsFillMasks[fillType])

            try :
                data    = ma.array(epsArr,    mask=totalMask.mask)
                cloudConfArr    = ma.array(cloudConfArr   ,mask=totalMask.mask)
                waterCldFracArr = ma.array(waterCldFracArr,mask=totalMask.mask)
                mixedCldFracArr = ma.array(mixedCldFracArr,mask=totalMask.mask)
                lats    = ma.array(latArr,    mask=totalMask.mask)
                lons    = ma.array(lonArr,    mask=totalMask.mask)
            except ma.core.MaskError :
                print ">> error: Mask Error, probably mismatched geolocation and product array sizes, aborting..."
                sys.exit(1)

            data = data * epsFactors[0] + epsFactors[1]

            cldPhase = np.ones(data.shape,dtype=np.uint8) * trimObj.sdrTypeFill['NA_FILL']['uint8']

            cloudFracMask = (cloudConfArr>0)
            waterFracMask = (waterCldFracArr>0)
            mixedFracMask = (mixedCldFracArr>0)

            waterMask = cloudFracMask * (waterFracMask + mixedFracMask)
            iceMask = cloudFracMask * np.logical_not(waterFracMask + mixedFracMask)

            print "Number of pixels with > 25% cloudy confidence: ",cloudFracMask.sum()
            print "Number of pixels with > 25% water cloud fraction: ",waterFracMask.sum()
            print "Number of pixels with > 25% mixed cloud fraction: ",mixedFracMask.sum()

            print "Number of water pixels: ",(cloudFracMask * (waterFracMask + mixedFracMask)).sum()
            print "Number of ice pixels: ",cloudFracMask.sum() - (waterFracMask + mixedFracMask).sum()

            print "Number of water pixels: ",waterMask.sum()
            print "Number of ice pixels: ",iceMask.sum()

            cldPhaseShape = cldPhase.shape
            cldPhase = np.ravel(cldPhase)

            wtrIdx = np.where(np.ravel(waterMask))
            iceIdx = np.where(np.ravel(iceMask))
            cldPhase[wtrIdx] = 3
            cldPhase[iceIdx] = 5

            cldPhase = np.reshape(cldPhase,cldPhaseShape)

        except Exception, err :
            print ">> error: %s..." % (str(err))
            sys.exit(1)

    print "lat_0,lon_0 = %f,%f" % (lat_0,lon_0)
    print "Shape of data is %s" % (repr(np.shape(data)))
    print "Shape of cldPhase is %s" % (repr(np.shape(cldPhase)))
    print "Shape of lats is %s" % (repr(np.shape(lats)))
    print "Shape of lons is %s" % (repr(np.shape(lons)))

    return lats,lons,data,cldPhase,lat_0,lon_0

def gran_CTp(geoList,ctpList,dataSet,shrink=1):
    '''
    Returns the granulated CTp dataset
    '''

    try :
        reload(viirsCld)
        reload(ViirsData)
        del(viirsCldObj)
        del(latArr)
        del(lonArr)
        del(ctpArr)
        del(phaseArr)
    except :
        pass

    print "Creating viirsCldObj..."
    reload(viirsCld)
    viirsCldObj = viirsCld.viirsCld()
    print "done"

    # Determine the correct fillValue
    trimObj = ViirsData.ViirsTrimTable()
    eps = 1.e-6

    for grans in np.arange(len(geoList)):
        print "Ingesting granule %d ..." % (grans)
        retArr = viirsCldObj.ingestLite(geoList[grans],ctpList[grans],dataSet,1)
        print "done"

        try :
            latArr  = viirsCldObj.Lat[:,:]
            lonArr  = viirsCldObj.Lon[:,:]
            
            lat_0 = latArr[np.shape(latArr)[0]/2,np.shape(latArr)[1]/2]
            lon_0 = lonArr[np.shape(lonArr)[0]/2,np.shape(lonArr)[1]/2]

            badGeo = False
            if not (-90. <= lat_0 <= 90.) :
                print "\n>> error: Latitude of granule midpoint (%f) does not satisfy (-90. <= lat_0 <= 90.)\nfor file %s\n\taborting..." % (lat_0,geoList[grans])
                badGeo = True
            if not (-180. <= lat_0 <= 180.) :
                print "\n>> error: Longitude of granule midpoint (%f) does not satisfy (-180. <= lon_0 <= 180.)\nfor file %s\n\taborting..." % (lon_0,geoList[grans])
                badGeo = True

            if badGeo :
                sys.exit(1)

            ctpArr   = viirsCldObj.ViirsCProdSDS[:,:]
            phaseArr = viirsCldObj.ViirsCProdSDSphase[:,:]

            # Determine masks for each fill type, for the CTp IP
            ctpFillMasks = {}
            for fillType in trimObj.sdrTypeFill.keys() :
                fillValue = trimObj.sdrTypeFill[fillType][ctpArr.dtype.name]
                if 'float' in fillValue.__class__.__name__ :
                    ctpFillMasks[fillType] = ma.masked_inside(ctpArr,fillValue-eps,fillValue+eps).mask
                    if (ctpFillMasks[fillType].__class__.__name__ != 'ndarray') :
                        ctpFillMasks[fillType] = None
                elif 'int' in fillValue.__class__.__name__ :
                    ctpFillMasks[fillType] = ma.masked_equal(ctpArr,fillValue).mask
                    if (ctpFillMasks[fillType].__class__.__name__ != 'ndarray') :
                        ctpFillMasks[fillType] = None
                else :
                    print "Dataset was neither int not float... a worry"
                    pass

            # Construct the total mask from all of the various fill values
            totalMask = ma.array(np.zeros(ctpArr.shape,dtype=np.bool))
            for fillType in trimObj.sdrTypeFill.keys() :
                if ctpFillMasks[fillType] is not None :
                    totalMask = totalMask * ma.array(np.zeros(ctpArr.shape,dtype=np.bool),\
                        mask=ctpFillMasks[fillType])

            # Added the "bad data" from the CTp phase quality flag.
            totalMask = totalMask * ma.masked_equal(phaseArr,0)

            try :
                data  = ma.array(ctpArr,  mask=totalMask.mask)
                phase = ma.array(phaseArr,mask=totalMask.mask)
                lats  = ma.array(latArr,  mask=totalMask.mask)
                lons  = ma.array(lonArr,  mask=totalMask.mask)
            except ma.core.MaskError :
                print ">> error: Mask Error, probably mismatched geolocation and product array sizes, aborting..."
                sys.exit(1)

        except :
            print ">> error: There was an exception..."
            sys.exit(1)

    return lats,lons,data,phase,lat_0,lon_0

def gran_CTT_EDR(geoList,cttList,shrink=1):
    '''
    Returns the granulated EPS EDR
    '''

    try :
        reload(viirsCld)
        reload(ViirsData)
        del(viirsCldObj)
        del(latArr)
        del(lonArr)
        del(cttArr)
    except :
        pass

    print "Creating viirsCldObj..."
    reload(viirsCld)
    viirsCldObj = viirsCld.viirsCld()
    print "done"

    # Determine the correct fillValue
    trimObj = ViirsData.ViirsTrimTable()
    eps = 1.e-6

    for grans in np.arange(len(geoList)):
        print "Ingesting EDR granule %d ..." % (grans)

        # Read in geolocation...
        ViirsGeoFileName = geoList[grans]
        if(pytables.isHDF5File(ViirsGeoFileName)) :
            try :
                ViirsGeoFileObj = pytables.openFile(ViirsGeoFileName,mode='r')
                print "Successfully opened geolocation file",ViirsGeoFileName
            except IOError :
                print ">> error: Could not open geolocation file: ",ViirsGeoFileName
                sys.exit(1)

            # Detemine the geolocation group name and related information
            #group = getobj(ViirsGeoFileObj,'/All_Data')
            group = ViirsGeoFileObj.getNode('/All_Data')
            geoGroupName = '/All_Data/'+group.__members__[0]
            group._g_close()
            print "Geolocation Group : %s " % (geoGroupName)
            isEdrGeo = ('VIIRS-CLD-AGG-GEO_All' in geoGroupName)
            if not isEdrGeo :
                print ">> error: %s is not an EDR resolution cloud geolocation file\n\taborting..."  % (ViirsGeoFileName)
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
        ViirsEDRFileName = cttList[grans]
        if(pytables.isHDF5File(ViirsEDRFileName)) :
            try :
                ViirsEDRFileObj = pytables.openFile(ViirsEDRFileName,mode='r')
                print "Successfully opened edr file",ViirsEDRFileName
            except IOError :
                print ">> error: Could not open edr file: ",ViirsEDRFileName
                sys.exit(1)

            # Detemine the edr group name and related information
            #group = getobj(ViirsEDRFileObj,'/All_Data')
            group = ViirsEDRFileObj.getNode('/All_Data')
            edrGroupName = '/All_Data/'+group.__members__[0]
            group._g_close()
            print "Edr Group : %s " % (edrGroupName)
            isEdr = ('VIIRS-CTT-EDR_All' in edrGroupName)
            if not isEdr :
                print ">> error: %s is not an EDR resolution cloud file\n\taborting..."  % (ViirsEDRFileName)
                sys.exit(1)

            # Get the edr nodes...
            try :
                dataName = 'AverageCloudTopTemperature'
                cttNode = ViirsEDRFileObj.getNode(edrGroupName+'/'+dataName)
                print "Shape of %s node is %s" % (edrGroupName+'/'+dataName,repr(cttNode.shape))
                dataName = 'CTTFactors'
                cttFactorsNode = ViirsEDRFileObj.getNode(edrGroupName+'/'+dataName)
                print "Shape of %s node is %s" % (edrGroupName+'/'+dataName,repr(cttFactorsNode.shape))
                dataName = 'QF3_VIIRSCTTAVGEDR'
                QF3_VIIRSCTTAVGEDR_Node = ViirsEDRFileObj.getNode(edrGroupName+'/'+dataName)
                print "Shape of %s node is %s" % (edrGroupName+'/'+dataName,repr(QF3_VIIRSCTTAVGEDR_Node.shape))
            except pyEx.NoSuchNodeError :
                print "\n>> error: No required node %s/%s in %s\n\taborting..." % (edrGroupName,dataName,ViirsEDRFileName)
                ViirsEDRFileObj.close()
                sys.exit(1)

            # Make a copy of the edr node data, so we can close the file
            try :
                dataName = 'AverageCloudTopTemperature'
                print "Reading %s dataset..." % (dataName)
                ctt = cttNode[:,:]
                ctt = np.squeeze(ctt)
                cttNode.close()
                print "done"
                dataName = 'CTTFactors'
                print "Reading %s dataset..." % (dataName)
                cttFactors = cttFactorsNode[:]
                cttFactors = np.squeeze(cttFactors)
                cttFactorsNode.close()
                print "done"
                dataName = 'QF3_VIIRSCTTAVGEDR'
                print "Reading %s dataset..." % (dataName)
                cloudConf = np.bitwise_and(QF3_VIIRSCTTAVGEDR_Node[:,:],3) >> 0
                waterCldFrac = np.bitwise_and(QF3_VIIRSCTTAVGEDR_Node[:,:],12) >> 2
                mixedCldFrac = np.bitwise_and(QF3_VIIRSCTTAVGEDR_Node[:,:],192) >> 6
                QF3_VIIRSCTTAVGEDR_Node.close()
                print "done"
                print "Closing edr file"
                ViirsEDRFileObj.close()

                print "Shape of ctt is %s" % (repr(np.shape(ctt)))
                print "Shape of latsArr is %s" % (repr(np.shape(latArr)))
                print "Shape of lonsArr is %s" % (repr(np.shape(lonArr)))

            except :
                print "\n>> error: Could not retrieve %/% node data in %s\n\taborting..." % (edrGroupName,dataName,ViirsEDRFileName)
                cttNode.close()
                cttFactorsNode.close()
                waterCldFracNode.close()
                ViirsEDRFileObj.close()
                sys.exit(1)

        else :
            print "\n>> error: %s is not a HDF5 file,\n\taborting..." % (ViirsEDRFileName)
            sys.exit(1)
        
        print "Creating some masks"
        try :
            
            lat_0 = latArr[np.shape(latArr)[0]/2,np.shape(latArr)[1]/2]
            lon_0 = lonArr[np.shape(lonArr)[0]/2,np.shape(lonArr)[1]/2]

            badGeo = False
            if not (-90. <= lat_0 <= 90.) :
                print "\n>> error: Latitude of granule midpoint (%f) does not satisfy (-90. <= lat_0 <= 90.)\nfor file %s\n\taborting..." % (lat_0,geoList[grans])
                badGeo = True
            if not (-180. <= lat_0 <= 180.) :
                print "\n>> error: Longitude of granule midpoint (%f) does not satisfy (-180. <= lon_0 <= 180.)\nfor file %s\n\taborting..." % (lon_0,geoList[grans])
                badGeo = True

            if badGeo :
                sys.exit(1)

            cttArr  = ctt[:,:]

            cloudConfArr    = cloudConf[:,:]
            waterCldFracArr = waterCldFrac[:,:]
            mixedCldFracArr = mixedCldFrac[:,:]


            # Determine masks for each fill type, for the CTT EDR
            cttFillMasks = {}
            for fillType in trimObj.sdrTypeFill.keys() :
                fillValue = trimObj.sdrTypeFill[fillType][cttArr.dtype.name]
                if 'float' in fillValue.__class__.__name__ :
                    cttFillMasks[fillType] = ma.masked_inside(cttArr,fillValue-eps,fillValue+eps).mask
                    if (cttFillMasks[fillType].__class__.__name__ != 'ndarray') :
                        cttFillMasks[fillType] = None
                elif 'int' in fillValue.__class__.__name__ :
                    cttFillMasks[fillType] = ma.masked_equal(cttArr,fillValue).mask
                    if (cttFillMasks[fillType].__class__.__name__ != 'ndarray') :
                        cttFillMasks[fillType] = None
                else :
                    print "Dataset was neither int not float... a worry"
                    pass

            # Construct the total mask from all of the various fill values
            totalMask = ma.array(np.zeros(cttArr.shape,dtype=np.bool))
            for fillType in trimObj.sdrTypeFill.keys() :
                if cttFillMasks[fillType] is not None :
                    totalMask = totalMask * ma.array(np.zeros(cttArr.shape,dtype=np.bool),\
                        mask=cttFillMasks[fillType])

            try :
                data = ma.array(cttArr,mask=totalMask.mask)
                cloudConfArr    = ma.array(cloudConfArr   ,mask=totalMask.mask)
                waterCldFracArr = ma.array(waterCldFracArr,mask=totalMask.mask)
                mixedCldFracArr = ma.array(mixedCldFracArr,mask=totalMask.mask)
                lats = ma.array(latArr,mask=totalMask.mask)
                lons = ma.array(lonArr,mask=totalMask.mask)
            except ma.core.MaskError :
                print ">> error: Mask Error, probably mismatched geolocation and product array sizes, aborting..."
                sys.exit(1)

            data = data * cttFactors[0] + cttFactors[1]

            cldPhase = np.ones(data.shape,dtype=np.uint8) * trimObj.sdrTypeFill['NA_FILL']['uint8']

            cloudFracMask = (cloudConfArr>0)
            waterFracMask = (waterCldFracArr>0)
            mixedFracMask = (mixedCldFracArr>0)

            waterMask = cloudFracMask * (waterFracMask + mixedFracMask)
            iceMask = cloudFracMask * np.logical_not(waterFracMask + mixedFracMask)

            print "Number of pixels with > 25% cloudy confidence: ",cloudFracMask.sum()
            print "Number of pixels with > 25% water cloud fraction: ",waterFracMask.sum()
            print "Number of pixels with > 25% mixed cloud fraction: ",mixedFracMask.sum()

            print "Number of water pixels: ",(cloudFracMask * (waterFracMask + mixedFracMask)).sum()
            print "Number of ice pixels: ",cloudFracMask.sum() - (waterFracMask + mixedFracMask).sum()

            print "Number of water pixels: ",waterMask.sum()
            print "Number of ice pixels: ",iceMask.sum()

            cldPhaseShape = cldPhase.shape
            cldPhase = np.ravel(cldPhase)

            wtrIdx = np.where(np.ravel(waterMask))
            iceIdx = np.where(np.ravel(iceMask))
            cldPhase[wtrIdx] = 3
            cldPhase[iceIdx] = 5

        except Exception, err :
            print ">> error: %s..." % (str(err))
            sys.exit(1)

    print "lat_0,lon_0 = %f,%f" % (lat_0,lon_0)
    print "Shape of data is %s" % (repr(np.shape(data)))
    print "Shape of cldPhase is %s" % (repr(np.shape(cldPhase)))
    print "Shape of lats is %s" % (repr(np.shape(lats)))
    print "Shape of lons is %s" % (repr(np.shape(lons)))

    return lats,lons,data,cldPhase,lat_0,lon_0

def gran_CTP_EDR(geoList,ctpList,shrink=1):
    '''
    Returns the granulated EPS EDR
    '''

    try :
        reload(viirsCld)
        reload(ViirsData)
        del(viirsCldObj)
        del(latArr)
        del(lonArr)
        del(ctpArr)
    except :
        pass

    print "Creating viirsCldObj..."
    reload(viirsCld)
    viirsCldObj = viirsCld.viirsCld()
    print "done"

    # Determine the correct fillValue
    trimObj = ViirsData.ViirsTrimTable()
    eps = 1.e-6

    for grans in np.arange(len(geoList)):
        print "Ingesting EDR granule %d ..." % (grans)

        # Read in geolocation...
        ViirsGeoFileName = geoList[grans]
        if(pytables.isHDF5File(ViirsGeoFileName)) :
            try :
                ViirsGeoFileObj = pytables.openFile(ViirsGeoFileName,mode='r')
                print "Successfully opened geolocation file",ViirsGeoFileName
            except IOError :
                print ">> error: Could not open geolocation file: ",ViirsGeoFileName
                sys.exit(1)

            # Detemine the geolocation group name and related information
            #group = getobj(ViirsGeoFileObj,'/All_Data')
            group = ViirsGeoFileObj.getNode('/All_Data')
            geoGroupName = '/All_Data/'+group.__members__[0]
            group._g_close()
            print "Geolocation Group : %s " % (geoGroupName)
            isEdrGeo = ('VIIRS-CLD-AGG-GEO_All' in geoGroupName)
            if not isEdrGeo :
                print ">> error: %s is not an EDR resolution cloud geolocation file\n\taborting..."  % (ViirsGeoFileName)
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
        ViirsEDRFileName = ctpList[grans]
        if(pytables.isHDF5File(ViirsEDRFileName)) :
            try :
                ViirsEDRFileObj = pytables.openFile(ViirsEDRFileName,mode='r')
                print "Successfully opened edr file",ViirsEDRFileName
            except IOError :
                print ">> error: Could not open edr file: ",ViirsEDRFileName
                sys.exit(1)

            # Detemine the edr group name and related information
            #group = getobj(ViirsEDRFileObj,'/All_Data')
            group = ViirsEDRFileObj.getNode('/All_Data')
            edrGroupName = '/All_Data/'+group.__members__[0]
            group._g_close()
            print "Edr Group : %s " % (edrGroupName)
            isEdr = ('VIIRS-CTP-EDR_All' in edrGroupName)
            if not isEdr :
                print ">> error: %s is not an EDR resolution cloud file\n\taborting..."  % (ViirsEDRFileName)
                sys.exit(1)

            # Get the edr nodes...
            try :
                dataName = 'AverageCloudTopPressure'
                ctpNode = ViirsEDRFileObj.getNode(edrGroupName+'/'+dataName)
                print "Shape of %s node is %s" % (edrGroupName+'/'+dataName,repr(ctpNode.shape))
                dataName = 'CTPFactors'
                ctpFactorsNode = ViirsEDRFileObj.getNode(edrGroupName+'/'+dataName)
                print "Shape of %s node is %s" % (edrGroupName+'/'+dataName,repr(ctpFactorsNode.shape))
                dataName = 'QF3_VIIRSCTPAVGEDR'
                QF3_VIIRSCTPAVGEDR_Node = ViirsEDRFileObj.getNode(edrGroupName+'/'+dataName)
                print "Shape of %s node is %s" % (edrGroupName+'/'+dataName,repr(QF3_VIIRSCTPAVGEDR_Node.shape))
            except pyEx.NoSuchNodeError :
                print "\n>> error: No required node %s/%s in %s\n\taborting..." % (edrGroupName,dataName,ViirsEDRFileName)
                ViirsEDRFileObj.close()
                sys.exit(1)

            # Make a copy of the edr node data, so we can close the file
            try :
                dataName = 'AverageCloudTopPressure'
                print "Reading %s dataset..." % (dataName)
                ctp = ctpNode[:,:]
                ctp = np.squeeze(ctp)
                ctpNode.close()
                print "done"
                dataName = 'CTPFactors'
                print "Reading %s dataset..." % (dataName)
                ctpFactors = ctpFactorsNode[:]
                ctpFactors = np.squeeze(ctpFactors)
                ctpFactorsNode.close()
                print "done"
                dataName = 'QF3_VIIRSCTPAVGEDR'
                print "Reading %s dataset..." % (dataName)
                cloudConf = np.bitwise_and(QF3_VIIRSCTPAVGEDR_Node[:,:],3) >> 0
                waterCldFrac = np.bitwise_and(QF3_VIIRSCTPAVGEDR_Node[:,:],12) >> 2
                mixedCldFrac = np.bitwise_and(QF3_VIIRSCTPAVGEDR_Node[:,:],192) >> 6
                QF3_VIIRSCTPAVGEDR_Node.close()
                print "done"
                print "Closing edr file"
                ViirsEDRFileObj.close()

                print "Shape of ctp is %s" % (repr(np.shape(ctp)))
                print "Shape of latsArr is %s" % (repr(np.shape(latArr)))
                print "Shape of lonsArr is %s" % (repr(np.shape(lonArr)))

            except :
                print "\n>> error: Could not retrieve %/% node data in %s\n\taborting..." % (edrGroupName,dataName,ViirsEDRFileName)
                ctpNode.close()
                ctpFactorsNode.close()
                QF3_VIIRSCTPAVGEDR_Node.close()
                ViirsEDRFileObj.close()
                sys.exit(1)

        else :
            print "\n>> error: %s is not a HDF5 file,\n\taborting..." % (ViirsEDRFileName)
            sys.exit(1)
        
        print "Creating some masks"
        try :
            
            lat_0 = latArr[np.shape(latArr)[0]/2,np.shape(latArr)[1]/2]
            lon_0 = lonArr[np.shape(lonArr)[0]/2,np.shape(lonArr)[1]/2]

            badGeo = False
            if not (-90. <= lat_0 <= 90.) :
                print "\n>> error: Latitude of granule midpoint (%f) does not satisfy (-90. <= lat_0 <= 90.)\nfor file %s\n\taborting..." % (lat_0,geoList[grans])
                badGeo = True
            if not (-180. <= lat_0 <= 180.) :
                print "\n>> error: Longitude of granule midpoint (%f) does not satisfy (-180. <= lon_0 <= 180.)\nfor file %s\n\taborting..." % (lon_0,geoList[grans])
                badGeo = True

            if badGeo :
                sys.exit(1)

            ctpArr  = ctp[:,:]

            cloudConfArr    = cloudConf[:,:]
            waterCldFracArr = waterCldFrac[:,:]
            mixedCldFracArr = mixedCldFrac[:,:]

            # Determine masks for each fill type, for the CTT EDR
            ctpFillMasks = {}
            for fillType in trimObj.sdrTypeFill.keys() :
                fillValue = trimObj.sdrTypeFill[fillType][ctpArr.dtype.name]
                if 'float' in fillValue.__class__.__name__ :
                    ctpFillMasks[fillType] = ma.masked_inside(ctpArr,fillValue-eps,fillValue+eps).mask
                    if (ctpFillMasks[fillType].__class__.__name__ != 'ndarray') :
                        ctpFillMasks[fillType] = None
                elif 'int' in fillValue.__class__.__name__ :
                    ctpFillMasks[fillType] = ma.masked_equal(ctpArr,fillValue).mask
                    if (ctpFillMasks[fillType].__class__.__name__ != 'ndarray') :
                        ctpFillMasks[fillType] = None
                else :
                    print "Dataset was neither int not float... a worry"
                    pass

            # Construct the total mask from all of the various fill values
            totalMask = ma.array(np.zeros(ctpArr.shape,dtype=np.bool))
            for fillType in trimObj.sdrTypeFill.keys() :
                if ctpFillMasks[fillType] is not None :
                    totalMask = totalMask * ma.array(np.zeros(ctpArr.shape,dtype=np.bool),\
                        mask=ctpFillMasks[fillType])

            try :
                data = ma.array(ctpArr,mask=totalMask.mask)
                cloudConfArr    = ma.array(cloudConfArr   ,mask=totalMask.mask)
                waterCldFracArr = ma.array(waterCldFracArr,mask=totalMask.mask)
                mixedCldFracArr = ma.array(mixedCldFracArr,mask=totalMask.mask)
                lats = ma.array(latArr,mask=totalMask.mask)
                lons = ma.array(lonArr,mask=totalMask.mask)
            except ma.core.MaskError :
                print ">> error: Mask Error, probably mismatched geolocation and product array sizes, aborting..."
                sys.exit(1)

            data = data * ctpFactors[0] + ctpFactors[1]

            cldPhase = np.ones(data.shape,dtype=np.uint8) * trimObj.sdrTypeFill['NA_FILL']['uint8']

            cloudFracMask = (cloudConfArr>0)
            waterFracMask = (waterCldFracArr>0)
            mixedFracMask = (mixedCldFracArr>0)

            waterMask = cloudFracMask * (waterFracMask + mixedFracMask)
            iceMask = cloudFracMask * np.logical_not(waterFracMask + mixedFracMask)

            print "Number of pixels with > 25% cloudy confidence: ",cloudFracMask.sum()
            print "Number of pixels with > 25% water cloud fraction: ",waterFracMask.sum()
            print "Number of pixels with > 25% mixed cloud fraction: ",mixedFracMask.sum()

            print "Number of water pixels: ",(cloudFracMask * (waterFracMask + mixedFracMask)).sum()
            print "Number of ice pixels: ",cloudFracMask.sum() - (waterFracMask + mixedFracMask).sum()

            print "Number of water pixels: ",waterMask.sum()
            print "Number of ice pixels: ",iceMask.sum()

            cldPhaseShape = cldPhase.shape
            cldPhase = np.ravel(cldPhase)

            wtrIdx = np.where(np.ravel(waterMask))
            iceIdx = np.where(np.ravel(iceMask))
            cldPhase[wtrIdx] = 3
            cldPhase[iceIdx] = 5

            cldPhase = np.reshape(cldPhase,cldPhaseShape)

        except Exception, err :
            print ">> error: %s..." % (str(err))
            sys.exit(1)

    print "lat_0,lon_0 = %f,%f" % (lat_0,lon_0)
    print "Shape of data is %s" % (repr(np.shape(data)))
    print "Shape of cldPhase is %s" % (repr(np.shape(cldPhase)))
    print "Shape of lats is %s" % (repr(np.shape(lats)))
    print "Shape of lons is %s" % (repr(np.shape(lons)))

    return lats,lons,data,cldPhase,lat_0,lon_0

def gran_CTH_EDR(geoList,cthList,shrink=1):
    '''
    Returns the granulated EPS EDR
    '''

    try :
        reload(viirsCld)
        reload(ViirsData)
        del(viirsCldObj)
        del(latArr)
        del(lonArr)
        del(cthArr)
    except :
        pass

    print "Creating viirsCldObj..."
    reload(viirsCld)
    viirsCldObj = viirsCld.viirsCld()
    print "done"

    # Determine the correct fillValue
    trimObj = ViirsData.ViirsTrimTable()
    eps = 1.e-6

    for grans in np.arange(len(geoList)):
        print "Ingesting EDR granule %d ..." % (grans)

        # Read in geolocation...
        ViirsGeoFileName = geoList[grans]
        if(pytables.isHDF5File(ViirsGeoFileName)) :
            try :
                ViirsGeoFileObj = pytables.openFile(ViirsGeoFileName,mode='r')
                print "Successfully opened geolocation file",ViirsGeoFileName
            except IOError :
                print ">> error: Could not open geolocation file: ",ViirsGeoFileName
                sys.exit(1)

            # Detemine the geolocation group name and related information
            #group = getobj(ViirsGeoFileObj,'/All_Data')
            group = ViirsGeoFileObj.getNode('/All_Data')
            geoGroupName = '/All_Data/'+group.__members__[0]
            group._g_close()
            print "Geolocation Group : %s " % (geoGroupName)
            isEdrGeo = ('VIIRS-CLD-AGG-GEO_All' in geoGroupName)
            if not isEdrGeo :
                print ">> error: %s is not an EDR resolution cloud geolocation file\n\taborting..."  % (ViirsGeoFileName)
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
        ViirsEDRFileName = cthList[grans]
        if(pytables.isHDF5File(ViirsEDRFileName)) :
            try :
                ViirsEDRFileObj = pytables.openFile(ViirsEDRFileName,mode='r')
                print "Successfully opened edr file",ViirsEDRFileName
            except IOError :
                print ">> error: Could not open edr file: ",ViirsEDRFileName
                sys.exit(1)

            # Detemine the edr group name and related information
            #group = getobj(ViirsEDRFileObj,'/All_Data')
            group = ViirsEDRFileObj.getNode('/All_Data')
            edrGroupName = '/All_Data/'+group.__members__[0]
            group._g_close()
            print "Edr Group : %s " % (edrGroupName)
            isEdr = ('VIIRS-CTH-EDR_All' in edrGroupName)
            if not isEdr :
                print ">> error: %s is not an EDR resolution cloud file\n\taborting..."  % (ViirsEDRFileName)
                sys.exit(1)

            # Get the edr nodes...
            try :
                dataName = 'AverageCloudTopHeight'
                cthNode = ViirsEDRFileObj.getNode(edrGroupName+'/'+dataName)
                print "Shape of %s node is %s" % (edrGroupName+'/'+dataName,repr(cthNode.shape))
                dataName = 'CTHFactors'
                cthFactorsNode = ViirsEDRFileObj.getNode(edrGroupName+'/'+dataName)
                print "Shape of %s node is %s" % (edrGroupName+'/'+dataName,repr(cthFactorsNode.shape))
                dataName = 'QF3_VIIRSCTHAVGEDR'
                QF3_VIIRSCTHAVGEDR_Node = ViirsEDRFileObj.getNode(edrGroupName+'/'+dataName)
                print "Shape of %s node is %s" % (edrGroupName+'/'+dataName,repr(QF3_VIIRSCTHAVGEDR_Node.shape))
            except pyEx.NoSuchNodeError :
                print "\n>> error: No required node %s/%s in %s\n\taborting..." % (edrGroupName,dataName,ViirsEDRFileName)
                ViirsEDRFileObj.close()
                sys.exit(1)

            # Make a copy of the edr node data, so we can close the file
            try :
                dataName = 'AverageCloudTopHeight'
                print "Reading %s dataset..." % (dataName)
                cth = cthNode[:,:]
                cth = np.squeeze(cth)
                cthNode.close()
                print "done"
                dataName = 'CTHFactors'
                print "Reading %s dataset..." % (dataName)
                cthFactors = cthFactorsNode[:]
                cthFactors = np.squeeze(cthFactors)
                cthFactorsNode.close()
                print "done"
                dataName = 'QF3_VIIRSCTTAVGEDR'
                print "Reading %s dataset..." % (dataName)
                cloudConf = np.bitwise_and(QF3_VIIRSCTHAVGEDR_Node[:,:],3) >> 0
                waterCldFrac = np.bitwise_and(QF3_VIIRSCTHAVGEDR_Node[:,:],12) >> 2
                mixedCldFrac = np.bitwise_and(QF3_VIIRSCTHAVGEDR_Node[:,:],192) >> 6
                QF3_VIIRSCTHAVGEDR_Node.close()
                print "done"
                print "Closing edr file"
                ViirsEDRFileObj.close()

                print "Shape of cth is %s" % (repr(np.shape(cth)))
                print "Shape of latsArr is %s" % (repr(np.shape(latArr)))
                print "Shape of lonsArr is %s" % (repr(np.shape(lonArr)))

            except :
                print "\n>> error: Could not retrieve %/% node data in %s\n\taborting..." % (edrGroupName,dataName,ViirsEDRFileName)
                cthNode.close()
                cthFactorsNode.close()
                QF3_VIIRSCTHAVGEDR_Node.close()
                ViirsEDRFileObj.close()
                sys.exit(1)

        else :
            print "\n>> error: %s is not a HDF5 file,\n\taborting..." % (ViirsEDRFileName)
            sys.exit(1)
        
        print "Creating some masks"
        try :
            
            lat_0 = latArr[np.shape(latArr)[0]/2,np.shape(latArr)[1]/2]
            lon_0 = lonArr[np.shape(lonArr)[0]/2,np.shape(lonArr)[1]/2]

            badGeo = False
            if not (-90. <= lat_0 <= 90.) :
                print "\n>> error: Latitude of granule midpoint (%f) does not satisfy (-90. <= lat_0 <= 90.)\nfor file %s\n\taborting..." % (lat_0,geoList[grans])
                badGeo = True
            if not (-180. <= lat_0 <= 180.) :
                print "\n>> error: Longitude of granule midpoint (%f) does not satisfy (-180. <= lon_0 <= 180.)\nfor file %s\n\taborting..." % (lon_0,geoList[grans])
                badGeo = True

            if badGeo :
                sys.exit(1)

            cthArr  = cth[:,:]

            cloudConfArr    = cloudConf[:,:]
            waterCldFracArr = waterCldFrac[:,:]
            mixedCldFracArr = mixedCldFrac[:,:]

            # Determine masks for each fill type, for the CTH EDR
            cthFillMasks = {}
            for fillType in trimObj.sdrTypeFill.keys() :
                fillValue = trimObj.sdrTypeFill[fillType][cthArr.dtype.name]
                if 'float' in fillValue.__class__.__name__ :
                    cthFillMasks[fillType] = ma.masked_inside(cthArr,fillValue-eps,fillValue+eps).mask
                    if (cthFillMasks[fillType].__class__.__name__ != 'ndarray') :
                        cthFillMasks[fillType] = None
                elif 'int' in fillValue.__class__.__name__ :
                    cthFillMasks[fillType] = ma.masked_equal(cthArr,fillValue).mask
                    if (cthFillMasks[fillType].__class__.__name__ != 'ndarray') :
                        cthFillMasks[fillType] = None
                else :
                    print "Dataset was neither int not float... a worry"
                    pass

            # Construct the total mask from all of the various fill values
            totalMask = ma.array(np.zeros(cthArr.shape,dtype=np.bool))
            for fillType in trimObj.sdrTypeFill.keys() :
                if cthFillMasks[fillType] is not None :
                    totalMask = totalMask * ma.array(np.zeros(cthArr.shape,dtype=np.bool),\
                        mask=cthFillMasks[fillType])

            try :
                data = ma.array(cthArr,mask=totalMask.mask)
                cloudConfArr    = ma.array(cloudConfArr   ,mask=totalMask.mask)
                waterCldFracArr = ma.array(waterCldFracArr,mask=totalMask.mask)
                mixedCldFracArr = ma.array(mixedCldFracArr,mask=totalMask.mask)
                lats = ma.array(latArr,mask=totalMask.mask)
                lons = ma.array(lonArr,mask=totalMask.mask)
            except ma.core.MaskError :
                print ">> error: Mask Error, probably mismatched geolocation and product array sizes, aborting..."
                sys.exit(1)

            data = data * cthFactors[0] + cthFactors[1]

            cldPhase = np.ones(data.shape,dtype=np.uint8) * trimObj.sdrTypeFill['NA_FILL']['uint8']

            cloudFracMask = (cloudConfArr>0)
            waterFracMask = (waterCldFracArr>0)
            mixedFracMask = (mixedCldFracArr>0)

            waterMask = cloudFracMask * (waterFracMask + mixedFracMask)
            iceMask = cloudFracMask * np.logical_not(waterFracMask + mixedFracMask)

            print "Number of pixels with > 25% cloudy confidence: ",cloudFracMask.sum()
            print "Number of pixels with > 25% water cloud fraction: ",waterFracMask.sum()
            print "Number of pixels with > 25% mixed cloud fraction: ",mixedFracMask.sum()

            print "Number of water pixels: ",(cloudFracMask * (waterFracMask + mixedFracMask)).sum()
            print "Number of ice pixels: ",cloudFracMask.sum() - (waterFracMask + mixedFracMask).sum()

            print "Number of water pixels: ",waterMask.sum()
            print "Number of ice pixels: ",iceMask.sum()

            cldPhaseShape = cldPhase.shape
            cldPhase = np.ravel(cldPhase)

            wtrIdx = np.where(np.ravel(waterMask))
            iceIdx = np.where(np.ravel(iceMask))
            cldPhase[wtrIdx] = 3
            cldPhase[iceIdx] = 5

            cldPhase = np.reshape(cldPhase,cldPhaseShape)

        except Exception, err :
            print ">> error: %s..." % (str(err))
            sys.exit(1)

    print "lat_0,lon_0 = %f,%f" % (lat_0,lon_0)
    print "Shape of data is %s" % (repr(np.shape(data)))
    print "Shape of cldPhase is %s" % (repr(np.shape(cldPhase)))
    print "Shape of lats is %s" % (repr(np.shape(lats)))
    print "Shape of lons is %s" % (repr(np.shape(lons)))

    return lats,lons,data,cldPhase,lat_0,lon_0

def gran_AOT(geoList,aotList,shrink=1):
    '''
    Returns the granulated AOT
    '''

    try :
        reload(viirsAero)
        reload(ViirsData)
        del(viirsAeroObj)
        del(latArr)
        del(lonArr)
        del(aotArr)
        del(retArr)
        del(qualArr)
        del(lsmArr)
        del(newData)
        del(dataIdx)
    except :
        pass

    print "Creating viirsSdrObj..."
    reload(viirsAero)
    viirsAeroObj = viirsAero.viirsAero()
    print "done"

    # Determine the correct fillValue
    trimObj = ViirsData.ViirsTrimTable()
    onboard_pt_value = trimObj.sdrTypeFill['ONBOARD_PT_FILL']['float32']
    onground_pt_value = trimObj.sdrTypeFill['ONGROUND_PT_FILL']['float32']
    na_fill_value = trimObj.sdrTypeFill['NA_FILL']['float32']
    eps = 1.e-6
    print "onboard_pt_value = ",onboard_pt_value
    print "onground_pt_value = ",onground_pt_value
    print "na_fill_value = ",na_fill_value

    for grans in np.arange(len(geoList)):

        print "\nIngesting granule %d ..." % (grans)
        retList = viirsAeroObj.ingest(geoList[grans],aotList[grans],'aot',shrink,'linear')

        try :
            latArr  = np.vstack((latArr,viirsAeroObj.Lat[:,:]))
            lonArr  = np.vstack((lonArr,viirsAeroObj.Lon[:,:]))
            ModeGran = viirsAeroObj.ModeGran
            print "subsequent geo arrays..."
        except NameError :
            latArr  = viirsAeroObj.Lat[:,:]
            lonArr  = viirsAeroObj.Lon[:,:]
            ModeGran = viirsAeroObj.ModeGran
            print "first geo arrays..."

        try :
            aotArr  = np.vstack((aotArr ,viirsAeroObj.ViirsAProdSDS[:,:]))
            retArr  = np.vstack((retArr ,viirsAeroObj.ViirsAProdRet[:,:]))
            qualArr = np.vstack((qualArr,viirsAeroObj.ViirsCMquality[:,:]))
            lsmArr  = np.vstack((lsmArr ,viirsAeroObj.LandSeaMask[:,:]))
            print "subsequent aot arrays..."
        except NameError :
            aotArr  = viirsAeroObj.ViirsAProdSDS[:,:]
            retArr  = viirsAeroObj.ViirsAProdRet[:,:]
            qualArr = viirsAeroObj.ViirsCMquality[:,:]
            lsmArr  = viirsAeroObj.LandSeaMask[:,:]
            print "first aot arrays..."

        print "Intermediate aotArr.shape = %s" % (str(aotArr.shape))
        print "Intermediate retArr.shape = %s" % (str(retArr.shape))
        print "Intermediate qualArr.shape = %s" % (str(qualArr.shape))
        print "Intermediate lsmArr.shape = %s" % (str(lsmArr.shape))

    lat_0 = latArr[np.shape(latArr)[0]/2,np.shape(latArr)[1]/2]
    lon_0 = lonArr[np.shape(lonArr)[0]/2,np.shape(lonArr)[1]/2]

    print "lat_0,lon_0 = ",lat_0,lon_0

    try :
        # Determine masks for each fill type, for the VCM IP
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

        # Define the onboard and onground pixel trim and N/A masks
        onboard_pt_mask = ma.masked_inside(aotArr,onboard_pt_value-eps,onboard_pt_value+eps)
        onground_pt_mask = ma.masked_inside(aotArr,onground_pt_value-eps,onground_pt_value+eps)
        na_fill_mask = ma.masked_inside(aotArr,na_fill_value-eps,na_fill_value+eps)

        # Define the product CM quality and Aerosol retrieval type masks
        ViirsCMqualityMask = ma.masked_equal(qualArr,0)
        ViirsAProdRetMask  = ma.masked_not_equal(retArr,0)
        missingMask        = ma.masked_less(aotArr,-0.)

        # Define the land and water masks
        ViirsLandMask      = ma.masked_greater(lsmArr,1)
        ViirsWaterMask     = ma.masked_less(lsmArr,2)

        # Define the total mask
        totalMask = onboard_pt_mask * onground_pt_mask * \
                    ViirsCMqualityMask * ViirsAProdRetMask * \
                    missingMask * na_fill_mask

        try :
            data = ma.array(aotArr,mask=totalMask.mask)
            lats = ma.array(latArr,mask=totalMask.mask)
            lons = ma.array(lonArr,mask=totalMask.mask)
        except ma.core.MaskError :
            print ">> error: Mask Error, probably mismatched geolocation and product array sizes, aborting..."
            sys.exit(1)

    except Exception, err :
        print ">> error: %s..." % (str(err))
        sys.exit(1)

    print "gran_AOT ModeGran = ",ModeGran

    return lats,lons,data,lat_0,lon_0,ModeGran

def gran_AOT_EDR(geoList,aotList,shrink=1):
    '''
    Returns the granulated AOT EDR
    '''

    try :
        reload(viirsAero)
        reload(ViirsData)
        del(viirsAeroObj)
        del(latArr)
        del(lonArr)
        del(aotArr)
        del(retArr)
        del(qualArr)
        del(lsmArr)
        del(newData)
        del(dataIdx)
    except :
        pass

    print "Creating viirsSdrObj..."
    reload(viirsAero)
    viirsAeroObj = viirsAero.viirsAero()
    print "done"

    # Determine the correct fillValue
    trimObj = ViirsData.ViirsTrimTable()
    eps = 1.e-6

    for grans in np.arange(len(geoList)):
        print "Ingesting EDR granule %d ..." % (grans)

        # Read in geolocation...
        ViirsGeoFileName = geoList[grans]
        if(pytables.isHDF5File(ViirsGeoFileName)) :
            try :
                ViirsGeoFileObj = pytables.openFile(ViirsGeoFileName,mode='r')
                print "Successfully opened geolocation file",ViirsGeoFileName
            except IOError :
                print ">> error: Could not open geolocation file: ",ViirsGeoFileName
                sys.exit(1)

            # Detemine the geolocation group name and related information
            #group = getobj(ViirsGeoFileObj,'/All_Data')
            group = ViirsGeoFileObj.getNode('/All_Data')
            geoGroupName = '/All_Data/'+group.__members__[0]
            group._g_close()
            print "Geolocation Group : %s " % (geoGroupName)
            isEdrGeo = ('VIIRS-Aeros-EDR-GEO_All' in geoGroupName)
            if not isEdrGeo :
                print ">> error: %s is not an EDR resolution aerosol geolocation file\n\taborting..."  % (ViirsGeoFileName)
                sys.exit(1)

            # Determine if this is a day/night/both granule
            try:
                dataNode = ViirsGeoFileObj.getNode('/Data_Products/VIIRS-Aeros-EDR-GEO/VIIRS-Aeros-EDR-GEO_Gran_0')
                dayNightFlag = np.squeeze(dataNode.attrs['N_Day_Night_Flag'])
            except Exception, err :
                print ">> error: %s..." % (str(err))
                dataNode.close()

            # Get the geolocation Latitude and Longitude nodes...
            try :
                dataName = 'Latitude'
                geoLatNode = ViirsGeoFileObj.getNode(geoGroupName+'/'+dataName)
                print "Shape of %s node is %s" % (geoGroupName+'/'+dataName,repr(geoLatNode.shape))
                dataName = 'Longitude'
                geoLonNode = ViirsGeoFileObj.getNode(geoGroupName+'/'+dataName)
                print "Shape of %s node is %s" % (geoGroupName+'/'+dataName,repr(geoLonNode.shape))
                dataName = 'SolarZenithAngle'
                geoSzaNode = ViirsGeoFileObj.getNode(geoGroupName+'/'+dataName)
                print "Shape of %s node is %s" % (geoGroupName+'/'+dataName,repr(geoSzaNode.shape))
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
                dataName = 'SolarZenithAngle'
                print "Reading %s dataset..." % (dataName)
                szaArr = geoSzaNode[:,:]
                geoSzaNode.close()
                szaArr = np.squeeze(szaArr)
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
        
        # Try to determine if this is a day or night granule. Night will be determined
        # to be a SZA greater than 85 degrees.

        if vars().has_key('dayNightFlag'):
            if dayNightFlag=='Day':
                ModeGran = 1
            if dayNightFlag=='Night':
                ModeGran = 0
            if dayNightFlag=='Both':
                ModeGran = 2
            print "ModeGran = %d" % (ModeGran)
        else :
            szaMask = ma.masked_less(szaArr,85.).mask
            dayFraction = float(szaMask.sum())/float(szaArr.size)
            print "Day Fraction = %f" % (dayFraction)
            ModeGran = 2 # Default to a mixed day/night granule
            if dayFraction == 1.00 : ModeGran = 1
            if dayFraction == 0.00 : ModeGran = 0
            print "ModeGran = %d" % (ModeGran)

        # Read in dataSets...
        ViirsEDRFileName = aotList[grans]
        if(pytables.isHDF5File(ViirsEDRFileName)) :
            try :
                ViirsEDRFileObj = pytables.openFile(ViirsEDRFileName,mode='r')
                print "Successfully opened edr file",ViirsEDRFileName
            except IOError :
                print ">> error: Could not open edr file: ",ViirsEDRFileName
                sys.exit(1)

            # Detemine the edr group name and related information
            #group = getobj(ViirsEDRFileObj,'/All_Data')
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
                print "Closing edr file"
                ViirsEDRFileObj.close()

                print "Shape of aot550 is %s" % (repr(np.shape(aot550)))
                print "Shape of latsArr is %s" % (repr(np.shape(latArr)))
                print "Shape of lonsArr is %s" % (repr(np.shape(lonArr)))

            except :
                print "\n>> error: Could not retrieve %/% node data in %s\n\taborting..." % (edrGroupName,dataName,ViirsEDRFileName)
                aot550Node.close()
                aotFactorsNode.close()
                ViirsEDRFileObj.close()
                sys.exit(1)

        else :
            print "\n>> error: %s is not a HDF5 file,\n\taborting..." % (ViirsEDRFileName)
            sys.exit(1)
        
        print "Creating some masks"
        try :
            
            lat_0 = latArr[np.shape(latArr)[0]/2,np.shape(latArr)[1]/2]
            lon_0 = lonArr[np.shape(lonArr)[0]/2,np.shape(lonArr)[1]/2]

            badGeo = False
            if not (-90. <= lat_0 <= 90.) :
                print "\n>> error: Latitude of granule midpoint (%f) does not satisfy (-90. <= lat_0 <= 90.)\nfor file %s\n\taborting..." % (lat_0,geoList[grans])
                badGeo = True
            if not (-180. <= lat_0 <= 180.) :
                print "\n>> error: Longitude of granule midpoint (%f) does not satisfy (-180. <= lon_0 <= 180.)\nfor file %s\n\taborting..." % (lon_0,geoList[grans])
                badGeo = True

            if badGeo :
                sys.exit(1)

            aotArr  = aot550[:,:]

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

            try :
                data = ma.array(aotArr,mask=totalMask.mask)
                lats = ma.array(latArr,mask=totalMask.mask)
                lons = ma.array(lonArr,mask=totalMask.mask)
            except ma.core.MaskError :
                print ">> error: Mask Error, probably mismatched geolocation and product array sizes, aborting..."
                sys.exit(1)

            data = data * aotFactors[0] + aotFactors[1]

        except :
            print ">> error: There was an exception..."
            sys.exit(1)

    print "lat_0,lon_0 = %f,%f" % (lat_0,lon_0)
    print "Shape of data is %s" % (repr(np.shape(data)))
    print "Shape of lats is %s" % (repr(np.shape(lats)))
    print "Shape of lons is %s" % (repr(np.shape(lons)))

    return lats,lons,data,lat_0,lon_0,ModeGran


###################################################
#                 Plotting Functions              #
###################################################

def orthoPlot_VCM(gridLat,gridLon,gridData,lat_0=0.,lon_0=0.,pointSize=1.,scale=1.3,mapRes='c',\
    prodFileName='',outFileName='out.png',dpi=300,titleStr='VIIRS Cloud Mask'):
    '''
    Plots the VIIRS Cloud Mask on an orthographic projection
    '''

    # The min and max values of the dataset
    vmin = np.min(ViirsData.CloudMaskData.ViirsCMvalues[cmByte][cmBit])
    vmax = np.max(ViirsData.CloudMaskData.ViirsCMvalues[cmByte][cmBit])

    # If we have a zero size data array, make a dummy dataset
    # to span the allowed data range, which will be plotted with 
    # vanishing pointsize
    if (np.shape(gridLon)[0]==0) :
        print "We have no valid data, synthesising dummy data..."
        gridLat = np.array([lat_0,lat_0])
        gridLon = np.array([lon_0,lon_0])
        gridData = np.array([vmin,vmax])
        pointSize = 0.001

    # Setup plotting data

    figWidth = 5. # inches
    figHeight = 4. # inches

    cmap = ListedColormap(ViirsData.CloudMaskData.ViirsCMfillColours[cmByte][cmBit])
    numCats = np.array(ViirsData.CloudMaskData.ViirsCMfillColours[cmByte][cmBit]).size
    numBounds = numCats + 1
    print "Number of Categories: %d" %(numCats)
    print "Number of Boundaries: %d" %(numBounds)

    # The tick positions in colourbar ([0..1]) coordinates
    ViirsData.CloudMaskData.ViirsCMTickPos = np.arange(float(numBounds))/float(numCats)
    ViirsData.CloudMaskData.ViirsCMTickPos = ViirsData.CloudMaskData.ViirsCMTickPos[0 :-1] + \
        ViirsData.CloudMaskData.ViirsCMTickPos[1]/2.

    print "ViirsData.CloudMaskData.ViirsCMTickPos: %r" %(ViirsData.CloudMaskData.ViirsCMTickPos)

    # Create figure with default size, and create canvas to draw on
    fig = Figure(figsize=((figWidth,figHeight)))
    canvas = FigureCanvas(fig)

    # Create main axes instance, leaving room for colorbar at bottom,
    # and also get the Bbox of the axes instance
    ax_rect = [0.05, 0.18, 0.9, 0.75  ] # [left,bottom,width,height]
    ax = fig.add_axes(ax_rect)

    # Granule axis title
    ax_title = ppl.setp(ax,title=prodFileName)
    ppl.setp(ax_title,fontsize=6)
    ppl.setp(ax_title,family="monospace")

    # Create Basemap instance
    # set 'ax' keyword so pylab won't be imported.
    windowWidth = scale *(0.35*12000000.)
    windowHeight = scale *(0.50*9000000.)
    m = Basemap(width=windowWidth,height=windowHeight,projection='lcc',lon_0=lon_0,lat_0=lat_0,ax=ax,fix_aspect=True,resolution=mapRes)
    #m = Basemap(width=0.35*12000000.,height=0.65*9000000.,projection='lcc',lon_0=lon_0,lat_0=lat_0,ax=ax,fix_aspect=False,resolution=mapRes)
    #m = Basemap(width=0.75*12000000.,height=9000000.,projection='merc',lon_0=lon_0,lat_0=lat_0,ax=ax,fix_aspect=True,resolution=mapRes)
    #m = Basemap(projection='ortho',lon_0=lon_0,lat_0=lat_0,ax=ax,fix_aspect=True,resolution=mapRes)
    #m = Basemap(projection='ortho',lon_0=lon_0,lat_0=lat_0,\
        #ax=ax,fix_aspect=True,resolution=mapRes,\
        #llcrnrx = -1. * scale * 3200. * 750./2.,\
        #llcrnry = -1. * scale * 3200. * 750./2.,\
        #urcrnrx =       scale * 3200. * 750./2.,\
        #urcrnry =       scale * 3200. * 750./2.)

    x,y=m(np.array(gridLon),np.array(gridLat))

    # Some map style configuration stufff
    #m.drawlsmask(ax=ax,land_color='gray',ocean_color='black',lakes=True)
    m.drawmapboundary(ax=ax,linewidth=0.01,fill_color='grey')
    m.drawcoastlines(ax=ax,linewidth=0.3,color='black')
    #m.fillcontinents(ax=ax,color='gray',lake_color='black',zorder=0)
    #m.drawparallels(np.arange(-90.,120.,10.),linewidth=0.5,color='white')
    #m.drawmeridians(np.arange(0.,420.,10.),linewidth=0.5,color='white')
    #m.drawlsmask(ax=ax,land_color='grey',ocean_color='black',lakes=True)

    m.bluemarble()

    # Plot the granule data
    if cmByte==0 and cmBit==1 :
        vmin,vmax = -0.5,3.5 # FIXME : This is temporary
    elif cmByte==5 and cmBit==0 :
        vmin,vmax = -0.5,7.5 # FIXME : This is temporary

    #cs = m.scatter(x,y,s=pointSize,c=gridData,axes=ax,edgecolors='none',vmin=vmin,vmax=vmax,cmap=cmap)
    cs = m.pcolor(x,y,gridData,axes=ax,edgecolors='none',vmin=vmin,vmax=vmax,cmap=cmap,antialiased=False)

    # add a colorbar axis
    cax_rect = [0.05 , 0.05, 0.9 , 0.06 ] # [left,bottom,width,height]
    cax = fig.add_axes(cax_rect,frameon=False) # setup colorbar axes

    # Plot the colourbar
    cb = fig.colorbar(cs, cax=cax, orientation='horizontal')

    # Set the colourbar tick locations and ticklabels
    #tickPos = np.array([0,1,2,3])

    # Convert the tick positions to data coordinates
    tickPos = ViirsData.CloudMaskData.ViirsCMTickPos  * numCats - 0.5
    tickLabels = ViirsData.CloudMaskData.ViirsCMtickNames[cmByte][cmBit]

    #tickPos = [0.25,0.5,0.75]
    #tickLabels = ['0.25','0.5','0.75']
    #tickPos = [1.,2.,3.,4.]
    #tickLabels = ['1.','2.','3.','4.']

    print "tickPos: %r" %(tickPos)
    print "tickLabels: %r" %(tickLabels)

    # Old style...
    #ppl.setp(cax,xticks=tickPos)
    #ppl.setp(cax,xticklabels=tickLabels)
    # New style...
    cb.set_ticks(tickPos)
    cb.set_ticklabels(tickLabels)

    # Turn off the tickmarks on the colourbar
    ppl.setp(cax.get_xticklines(),visible=False)
    ppl.setp(cax.get_xticklabels(),fontsize=9)

    #
    # Add a small globe with the swath indicated on it #
    #

    # Create main axes instance, leaving room for colorbar at bottom,
    # and also get the Bbox of the axes instance
    glax_rect = [0.81, 0.75, 0.18, 0.20 ] # [left,bottom,width,height]
    glax = fig.add_axes(glax_rect)

    m_globe = Basemap(lat_0=0.,lon_0=0.,\
        ax=glax,resolution='c',area_thresh=10000.,projection='robin')

    # If we previously had a zero size data array, increase the pointSize
    # so the data points are visible on the global plot
    pointSize = 0.2
    if (np.shape(gridLon)[0]==2) :
        print "We have no valid data, synthesising dummy data..."
        pointSize = 5.

    x,y=m_globe(np.array(gridLon),np.array(gridLat))
    swath = np.zeros(np.shape(x),dtype=int)

    m_globe.drawcoastlines(ax=glax,linewidth=0.1)
    m_globe.fillcontinents(ax=glax,color='gray',zorder=0)
    m_globe.drawmapboundary(linewidth=0.1)

    p_globe = m_globe.scatter(x,y,s=pointSize,c="red",axes=glax,edgecolors='none')

    # Globe axis title
    glax_xlabel = ppl.setp(glax,xlabel=str(titleStr))
    ppl.setp(glax_xlabel,fontsize=6)

    # Redraw the figure
    canvas.draw()

    # save image 
    print "Writing file to ",outFileName
    canvas.print_figure(outFileName,dpi=dpi)

def orthoPlot_AOT(gridLat,gridLon,gridData,ModeGran, \
        vmin=-0.05,vmax=0.8,scale=1.3, \
        lat_0=0.,lon_0=0.,pointSize=1.,mapRes='c',cmap=None, \
        prodFileName='',outFileName='out.png',dpi=300,titleStr='VIIRS AOT'):
    '''
    Plots the VIIRS Aerosol Optical Thickness on an orthographic projection
    '''

    reload(ViirsData)

    # The plot range...
    print "vmin,vmax = ",vmin,vmax 

    # If we have a zero size data array, make a dummy dataset
    # to span the allowed data range, which will be plotted with 
    # vanishing pointsize
    if (np.shape(gridLon)[0]==0) :
        print "We have no valid data, synthesising dummy data..."
        gridLat = np.array([lat_0,lat_0])
        gridLon = np.array([lon_0,lon_0])
        gridData = np.array([0.,1.])
        pointSize = 0.001

    # Setup plotting data

    figWidth = 5. # inches
    figHeight = 4. # inches

    # Create figure with default size, and create canvas to draw on
    fig = Figure(figsize=((figWidth,figHeight)))
    canvas = FigureCanvas(fig)

    # Create main axes instance, leaving room for colorbar at bottom,
    # and also get the Bbox of the axes instance
    ax_rect = [0.05, 0.18, 0.9, 0.75  ] # [left,bottom,width,height]
    ax = fig.add_axes(ax_rect)

    # Granule axis title
    ax_title = ppl.setp(ax,title=prodFileName)
    ppl.setp(ax_title,fontsize=6)
    ppl.setp(ax_title,family="monospace")

    # Create Basemap instance
    # set 'ax' keyword so pylab won't be imported.
    windowWidth = scale *(0.35*12000000.)
    windowHeight = scale *(0.50*9000000.)
    m = Basemap(width=windowWidth,height=windowHeight,projection='lcc',lon_0=lon_0,lat_0=lat_0,ax=ax,fix_aspect=True,resolution=mapRes)
    #m = Basemap(width=0.35*12000000.,height=0.65*9000000.,projection='lcc',lon_0=lon_0,lat_0=lat_0,ax=ax,fix_aspect=False,resolution=mapRes)
    #m = Basemap(width=0.75*12000000.,height=9000000.,projection='merc',lon_0=lon_0,lat_0=lat_0,ax=ax,fix_aspect=True,resolution=mapRes)
    #m = Basemap(projection='ortho',lon_0=lon_0,lat_0=lat_0,ax=ax,fix_aspect=True,resolution=mapRes)
    #m = Basemap(projection='ortho',lon_0=lon_0,lat_0=lat_0,\
        #ax=ax,fix_aspect=True,resolution=mapRes,\
        #llcrnrx = -1. * scale * 3200. * 750./2.,\
        #llcrnry = -1. * scale * 3200. * 750./2.,\
        #urcrnrx =       scale * 3200. * 750./2.,\
        #urcrnry =       scale * 3200. * 750./2.)


    x,y=m(np.array(gridLon),np.array(gridLat))

    # Some map style configuration stufff
    #m.drawlsmask(ax=ax,land_color='gray',ocean_color='black',lakes=True)
    m.drawmapboundary(ax=ax,linewidth=0.01,fill_color='black')
    m.drawcoastlines(ax=ax,linewidth=0.3,color='white')
    m.fillcontinents(ax=ax,color='gray',lake_color='black',zorder=0)
    #m.drawparallels(np.arange(-90.,120.,30.),color='white')
    #m.drawmeridians(np.arange(0.,420.,60.),color='white')

    #m.bluemarble()

    # Plot the granule data
    print "shape of gridData is %s" % (repr(np.shape(gridData)))
    #cs = m.scatter(x,y,s=pointSize,c=gridData,axes=ax,edgecolors='none',vmin=vmin,vmax=vmax,cmap=cmap)
    cs = m.pcolor(x,y,gridData,axes=ax,edgecolors='none',vmin=vmin,vmax=vmax,cmap=cmap,antialiased=False)

    print "orthoPlot_AOT ModeGran = ",ModeGran
    if (ModeGran == 0) :
        print "Printing NIGHT text"
        fig.text(0.5, 0.555, 'NIGHT',fontsize=30, color='white',ha='center',va='center',alpha=0.6)

    # add a colorbar axis
    cax_rect = [0.05 , 0.05, 0.9 , 0.06 ] # [left,bottom,width,height]
    cax = fig.add_axes(cax_rect,frameon=False) # setup colorbar axes

    # Plot the colorbar.
    cb = fig.colorbar(cs, cax=cax, orientation='horizontal')
    ppl.setp(cax.get_xticklabels(),fontsize=9)

    # Colourbar title
    cax_title = ppl.setp(cax,title="AOT")
    ppl.setp(cax_title,fontsize=9)

    #
    # Add a small globe with the swath indicated on it #
    #

    # Create main axes instance, leaving room for colorbar at bottom,
    # and also get the Bbox of the axes instance
    glax_rect = [0.81, 0.75, 0.18, 0.20 ] # [left,bottom,width,height]
    glax = fig.add_axes(glax_rect)

    m_globe = Basemap(lat_0=0.,lon_0=0.,\
        ax=glax,resolution='c',area_thresh=10000.,projection='robin')

    # If we previously had a zero size data array, increase the pointSize
    # so the data points are visible on the global plot
    if (np.shape(gridLon)[0]==2) :
        pointSize = 5.

    x,y=m_globe(np.array(gridLon),np.array(gridLat))
    swath = np.zeros(np.shape(x),dtype=int)

    m_globe.drawcoastlines(ax=glax,linewidth=0.1)
    m_globe.fillcontinents(ax=glax,color='gray',zorder=0)
    m_globe.drawmapboundary(linewidth=0.1)

    p_globe = m_globe.scatter(x,y,s=pointSize,c="red",axes=glax,edgecolors='none')

    # Globe axis title
    glax_xlabel = ppl.setp(glax,xlabel=titleStr)
    ppl.setp(glax_xlabel,fontsize=6)

    # Redraw the figure
    canvas.draw()

    # save image 
    print "Writing file to ",outFileName
    canvas.print_figure(outFileName,dpi=dpi)


def orthoPlot_COP(gridLat,gridLon,gridData,gridPhase,dataSet, \
        lat_0=0.,lon_0=0.,\
        abScale='log',pointSize=1.,scale=1.3,mapRes='c',\
        prodFileName='',outFileName='out.png',dpi=300,titleStr='VIIRS COP'):
    '''
    Plots the VIIRS Cloud Optical Parameters (COP) on an orthographic projection
    '''

    # Setup plotting data
    reload(ViirsData)
    CloudProduct = ViirsData.CloudProdData.CloudProduct[dataSet][abScale]

    pointSize_water = pointSize
    pointSize_ice = pointSize

    cirrusMask     = ma.masked_equal(gridPhase,1) # Ice
    opaqueIceMask  = ma.masked_equal(gridPhase,2) # Ice
    waterMask      = ma.masked_equal(gridPhase,3) # Water
    mixedMask      = ma.masked_equal(gridPhase,4) # Water
    multiLayerMask = ma.masked_equal(gridPhase,5) # Ice

    totalMask = cirrusMask * opaqueIceMask * multiLayerMask # Masks the ice
    #gridData_water = ma.compressed(ma.array(gridData,mask=totalMask.mask))
    #gridLat_water  = ma.compressed(ma.array(gridLat,mask=totalMask.mask))
    #gridLon_water  = ma.compressed(ma.array(gridLon,mask=totalMask.mask))
    gridData_water = ma.array(gridData,mask=totalMask.mask)
    gridLat_water  = ma.array(gridLat,mask=totalMask.mask)
    gridLon_water  = ma.array(gridLon,mask=totalMask.mask)


    # If we have a zero size water array, make a dummy dataset
    # to span the allowed data range, which will be plotted with 
    # vanishing pointsize
    if (np.shape(gridLon_water)[0]==0) :
        print "We have no valid water data, synthesising dummy data..."
        gridLat_water = np.array([lat_0,lat_0])
        gridLon_water = np.array([lon_0,lon_0])
        gridData_water = np.array([CloudProduct.vmin_water,CloudProduct.vmax_water])
        pointSize_water = 0.001

    totalMask = waterMask * mixedMask # Masks the water
    #gridData_ice = ma.compressed(ma.array(gridData,mask=totalMask.mask))
    #gridLat_ice  = ma.compressed(ma.array(gridLat,mask=totalMask.mask))
    #gridLon_ice  = ma.compressed(ma.array(gridLon,mask=totalMask.mask))
    gridData_ice = ma.array(gridData,mask=totalMask.mask)
    gridLat_ice  = ma.array(gridLat,mask=totalMask.mask)
    gridLon_ice  = ma.array(gridLon,mask=totalMask.mask)

    # If we have a zero size ice array, make a dummy dataset
    # to span the allowed data range, which will be plotted with 
    # vanishing pointsize
    if (np.shape(gridLon_ice)[0]==0) :
        print "We have no valid ice data, synthesising dummy data..."
        gridLat_ice = np.array([lat_0,lat_0])
        gridLon_ice = np.array([lon_0,lon_0])
        gridData_ice = np.array([CloudProduct.vmin_ice,CloudProduct.vmax_ice])
        pointSize_ice = 0.001

    # Make logarithmic if appropriate
    if(CloudProduct.logScale):
        print "Log scaling the dataset..."
        gridData_water = np.log10(gridData_water)
        gridData_ice = np.log10(gridData_ice)
        print "Finished log scaling the dataset..."

    vmin_water = CloudProduct.vmin_water
    vmax_water = CloudProduct.vmax_water

    vmin_ice = CloudProduct.vmin_ice
    vmax_ice = CloudProduct.vmax_ice

    cmap_water = CloudProduct.cmap_water
    cmap_ice = CloudProduct.cmap_ice

    figWidth = 5. # inches
    figHeight = 4. # inches

    # Create figure with default size, and create canvas to draw on
    fig = Figure(figsize=((figWidth,figHeight)))
    canvas = FigureCanvas(fig)

    # Create main axes instance, leaving room for colorbar at bottom,
    # and also get the Bbox of the axes instance
    ax_rect = [0.05, 0.18, 0.9, 0.75  ] # [left,bottom,width,height]
    ax = fig.add_axes(ax_rect)

    # Granule axis title
    ax_title = ppl.setp(ax,title=prodFileName)
    ppl.setp(ax_title,fontsize=6)
    ppl.setp(ax_title,family="monospace")

    # Create Basemap instance
    # set 'ax' keyword so pylab won't be imported.
    #m = Basemap(projection='ortho',lon_0=lon_0,lat_0=lat_0,ax=ax,fix_aspect=True,resolution='c')
    m = Basemap(projection='ortho',lon_0=lon_0,lat_0=lat_0,\
        ax=ax,fix_aspect=True,resolution=mapRes,\
        llcrnrx = -1. * scale * 3200. * 750./2.,\
        llcrnry = -1. * scale * 3200. * 750./2.,\
        urcrnrx =       scale * 3200. * 750./2.,\
        urcrnry =       scale * 3200. * 750./2.)

    # Some map style configuration stufff
    #m.drawlsmask(ax=ax,land_color='gray',ocean_color='black',lakes=True)
    m.drawmapboundary(ax=ax,linewidth=0.01,fill_color='grey')
    m.drawcoastlines(ax=ax,linewidth=0.3,color='white')
    #m.fillcontinents(ax=ax,color='gray',lake_color='black',zorder=0)
    #m.drawparallels(np.arange(-90.,120.,30.),color='white')
    #m.drawmeridians(np.arange(0.,420.,60.),color='white')

    # Plot the granule water data
    x,y=m(np.array(gridLon_water),np.array(gridLat_water))
    #cs_water = m.scatter(x,y,s=pointSize_water,c=gridData_water,axes=ax,edgecolors='none',vmin=vmin_water,vmax=vmax_water,cmap=cmap_water)
    cs_water = m.pcolor(x,y,gridData_water,axes=ax,edgecolors='none',vmin=vmin_water,vmax=vmax_water,cmap=cmap_water,antialiased=False)

    # Plot the granule ice data
    x,y=m(np.array(gridLon_ice),np.array(gridLat_ice))
    #cs_ice = m.scatter(x,y,s=pointSize_ice,c=gridData_ice,axes=ax,edgecolors='none',vmin=vmin_ice,vmax=vmax_ice,cmap=cmap_ice)
    cs_ice = m.pcolor(x,y,gridData_ice,axes=ax,edgecolors='none',vmin=vmin_ice,vmax=vmax_ice,cmap=cmap_ice,antialiased=False)

    ### Colourbars
    cbar_WidthTotal = 0.9
    cbar_HeightTotal= 0.06
    cbarWidth = (cbar_WidthTotal - 0.05)/2.
    cbX0,cbY0,cbdX,cbdY  = 0.05,0.05,cbarWidth,cbar_HeightTotal
    rect_cbar_water = [cbX0,cbY0,cbdX,cbdY] # [left,bottom,width,height]
    cbX0,cbY0,cbdX,cbdY  = cbX0+cbarWidth+0.05,0.05,cbarWidth,cbar_HeightTotal
    rect_cbar_ice   = [cbX0,cbY0,cbdX,cbdY] # [left,bottom,width,height]

    # Add the water colorbar axis
    cax_water = fig.add_axes(rect_cbar_water,frameon=False) # setup colorbar axes
    if(CloudProduct.logScale):
        print "vmin_water,vmax_water = ",vmin_water,vmax_water
        cIndices = np.arange(-10,10)
        minIndex = cIndices[np.where(cIndices <= vmin_water)][-1]
        maxIndex = cIndices[np.where(cIndices >= vmax_water)][0]
        #minIndex,maxIndex = -1,2
        cbTicks = np.arange(minIndex,maxIndex+1,dtype=int)
        cbTickLabels = []
        for ticks in cbTicks :
            if ticks < 0 :
                cbTickLabels.append('{0:.{1}f}'.format(10**(ticks),np.abs(ticks)))
            else :
                cbTickLabels.append('{0:.{1}f}'.format(10**(ticks),0))
        cb_water = fig.colorbar(cs_water, cax=cax_water, orientation='horizontal',ticks=cbTicks)
        ppl.setp(cax_water,xticklabels=cbTickLabels)
    else :
        cb_water = fig.colorbar(cs_water, cax=cax_water, orientation='horizontal')
    ppl.setp(cax_water.get_xticklabels(),fontsize=6)
    cax_title = ppl.setp(cax_water,title=CloudProduct.cbarTitle_water)
    ppl.setp(cax_title,fontsize=9)

    # Add the ice colorbar axis
    cax_ice = fig.add_axes(rect_cbar_ice,frameon=False) # setup colorbar axes
    if(CloudProduct.logScale):
        cb_ice = fig.colorbar(cs_ice, cax=cax_ice, orientation='horizontal',ticks=cbTicks)
        ppl.setp(cax_ice,xticklabels=cbTickLabels)
    else :
        cb_ice = fig.colorbar(cs_ice, cax=cax_ice, orientation='horizontal')
    ppl.setp(cax_ice.get_xticklabels(),fontsize=6)
    cax_title = ppl.setp(cax_ice,title=CloudProduct.cbarTitle_ice)
    ppl.setp(cax_title,fontsize=9)

    #
    # Add a small globe with the swath indicated on it #
    #

    # Create main axes instance, leaving room for colorbar at bottom,
    # and also get the Bbox of the axes instance
    glax_rect = [0.81, 0.75, 0.18, 0.20 ] # [left,bottom,width,height]
    glax = fig.add_axes(glax_rect)

    m_globe = Basemap(lat_0=0.,lon_0=0.,\
        ax=glax,resolution='c',area_thresh=10000.,projection='robin')

    # If we previously had a zero size data array, increase the pointSize
    # so the data points are visible on the global plot
    if (np.shape(gridLon)[0]==0) :
        print "We have no valid data, synthesising dummy data..."
        gridLat = np.array([lat_0,lat_0])
        gridLon = np.array([lon_0,lon_0])
        pointSize = 5.

    x,y=m_globe(np.array(gridLon),np.array(gridLat))
    swath = np.zeros(np.shape(x),dtype=int)

    m_globe.drawcoastlines(ax=glax,linewidth=0.1)
    m_globe.fillcontinents(ax=glax,color='gray',zorder=0)
    m_globe.drawmapboundary(linewidth=0.1)

    p_globe = m_globe.scatter(x,y,s=pointSize,c="red",axes=glax,edgecolors='none')

    # Globe axis title
    glax_xlabel = ppl.setp(glax,xlabel=str(titleStr))
    ppl.setp(glax_xlabel,fontsize=6)

    # Redraw the figure
    canvas.draw()

    # save image 
    print "Writing file to ",outFileName
    canvas.print_figure(outFileName,dpi=dpi)

def orthoPlot_CTp(gridLat,gridLon,gridData,gridPhase,dataSet,lat_0=0.,lon_0=0.,\
        pointSize=1.,scale=1.3,mapRes='c',prodFileName='',outFileName='out.png',dpi=300,titleStr='VIIRS CTp'):
    '''
    Plots the VIIRS Cloud Top Parameters (CTp) on an orthographic projection
    '''

    # Setup plotting data
    reload(ViirsData)
    CloudProduct = ViirsData.CloudProdData.CloudProduct[dataSet]

    pointSize_water = pointSize
    pointSize_ice = pointSize

    badDataMask    = ma.masked_equal(gridPhase,0)
    waterMask      = ma.masked_equal(gridPhase,1) # Water
    iceMask        = ma.masked_equal(gridPhase,2) # Ice
    mixedMask      = ma.masked_equal(gridPhase,3) # Water

    totalMask = iceMask # Masks the ice
    #gridData_water = ma.compressed(ma.array(gridData,mask=totalMask.mask))
    #gridLat_water  = ma.compressed(ma.array(gridLat,mask=totalMask.mask))
    #gridLon_water  = ma.compressed(ma.array(gridLon,mask=totalMask.mask))
    gridData_water = ma.array(gridData,mask=totalMask.mask)
    gridLat_water  = ma.array(gridLat,mask=totalMask.mask)
    gridLon_water  = ma.array(gridLon,mask=totalMask.mask)

    # If we have a zero size water array, make a dummy dataset
    # to span the allowed data range, which will be plotted with 
    # vanishing pointsize
    if (np.shape(gridLon_water)[0]==0) :
        print "We have no valid water data, synthesising dummy data..."
        gridLat_water = np.array([lat_0,lat_0])
        gridLon_water = np.array([lon_0,lon_0])
        gridData_water = np.array([CloudProduct.vmin_water,CloudProduct.vmax_water])
        pointSize_water = 0.001

    totalMask = waterMask * mixedMask # Masks the water
    #gridData_ice = ma.compressed(ma.array(gridData,mask=totalMask.mask))
    #gridLat_ice  = ma.compressed(ma.array(gridLat,mask=totalMask.mask))
    #gridLon_ice  = ma.compressed(ma.array(gridLon,mask=totalMask.mask))
    gridData_ice = ma.array(gridData,mask=totalMask.mask)
    gridLat_ice  = ma.array(gridLat,mask=totalMask.mask)
    gridLon_ice  = ma.array(gridLon,mask=totalMask.mask)

    # If we have a zero size ice array, make a dummy dataset
    # to span the allowed data range, which will be plotted with 
    # vanishing pointsize
    if (np.shape(gridLon_ice)[0]==0) :
        print "We have no valid ice data, synthesising dummy data..."
        gridLat_ice = np.array([lat_0,lat_0])
        gridLon_ice = np.array([lon_0,lon_0])
        gridData_ice = np.array([CloudProduct.vmin_ice,CloudProduct.vmax_ice])
        pointSize_ice = 0.001

    # Make logarithmic if appropriate
    if(CloudProduct.logScale):
        print "Log scaling the dataset..."
        gridData_water = np.log10(gridData_water)
        gridData_ice = np.log10(gridData_ice)
        print "Finished log scaling the dataset..."

    vmin_water = CloudProduct.vmin_water
    vmax_water = CloudProduct.vmax_water

    vmin_ice = CloudProduct.vmin_ice
    vmax_ice = CloudProduct.vmax_ice

    cmap_water = CloudProduct.cmap_water
    cmap_ice = CloudProduct.cmap_ice

    figWidth = 5. # inches
    figHeight = 4. # inches

    # Create figure with default size, and create canvas to draw on
    fig = Figure(figsize=((figWidth,figHeight)))
    canvas = FigureCanvas(fig)

    # Create main axes instance, leaving room for colorbar at bottom,
    # and also get the Bbox of the axes instance
    ax_rect = [0.05, 0.18, 0.9, 0.75  ] # [left,bottom,width,height]
    ax = fig.add_axes(ax_rect)

    # Granule axis title
    ax_title = ppl.setp(ax,title=prodFileName)
    ppl.setp(ax_title,fontsize=6)
    ppl.setp(ax_title,family="monospace")

    # Create Basemap instance
    # set 'ax' keyword so pylab won't be imported.
    #m = Basemap(projection='ortho',lon_0=lon_0,lat_0=lat_0,ax=ax,fix_aspect=True,resolution='c')
    m = Basemap(projection='ortho',lon_0=lon_0,lat_0=lat_0,\
        ax=ax,fix_aspect=True,resolution=mapRes,\
        llcrnrx = -1. * scale * 3200. * 750./2.,\
        llcrnry = -1. * scale * 3200. * 750./2.,\
        urcrnrx =       scale * 3200. * 750./2.,\
        urcrnry =       scale * 3200. * 750./2.)

    # Some map style configuration stufff
    #m.drawlsmask(ax=ax,land_color='gray',ocean_color='black',lakes=True)
    m.drawmapboundary(ax=ax,linewidth=0.01,fill_color='grey')
    m.drawcoastlines(ax=ax,linewidth=0.3,color='white')
    #m.fillcontinents(ax=ax,color='gray',lake_color='black',zorder=0)
    #m.drawparallels(np.arange(-90.,120.,30.),color='white')
    #m.drawmeridians(np.arange(0.,420.,60.),color='white')

    # Plot the granule data

    x,y=m(np.array(gridLon_water),np.array(gridLat_water))
    #cs_water = m.scatter(x,y,s=pointSize_water,c=gridData_water,axes=ax,edgecolors='none',vmin=vmin_water,vmax=vmax_water,cmap=cmap_water)
    cs_water = m.pcolor(x,y,gridData_water,axes=ax,edgecolors='none',vmin=vmin_water,vmax=vmax_water,cmap=cmap_water,antialiased=False)

    x,y=m(np.array(gridLon_ice),np.array(gridLat_ice))
    #cs_ice = m.scatter(x,y,s=pointSize_ice,c=gridData_ice,axes=ax,edgecolors='none',vmin=vmin_ice,vmax=vmax_ice,cmap=cmap_ice)
    cs_ice = m.pcolor(x,y,gridData_ice,axes=ax,edgecolors='none',vmin=vmin_ice,vmax=vmax_ice,cmap=cmap_ice,antialiased=False)

    ### Colourbars

    cbar_WidthTotal = 0.9
    cbar_HeightTotal= 0.06
    cbarWidth = (cbar_WidthTotal - 0.05)/2.
    cbX0,cbY0,cbdX,cbdY  = 0.05,0.05,cbarWidth,cbar_HeightTotal
    rect_cbar_water = [cbX0,cbY0,cbdX,cbdY] # [left,bottom,width,height]
    cbX0,cbY0,cbdX,cbdY  = cbX0+cbarWidth+0.05,0.05,cbarWidth,cbar_HeightTotal
    rect_cbar_ice   = [cbX0,cbY0,cbdX,cbdY] # [left,bottom,width,height]

    # Add the water colorbar axis
    cax_water = fig.add_axes(rect_cbar_water,frameon=False) # setup colorbar axes
    cb_water = fig.colorbar(cs_water, cax=cax_water, orientation='horizontal')
    ppl.setp(cax_water.get_xticklabels(),fontsize=6)
    cax_title = ppl.setp(cax_water,title=CloudProduct.cbarTitle_water)
    ppl.setp(cax_title,fontsize=9)

    # Add the ice colorbar axis
    cax_ice = fig.add_axes(rect_cbar_ice,frameon=False) # setup colorbar axes
    cb_ice = fig.colorbar(cs_ice, cax=cax_ice, orientation='horizontal')
    ppl.setp(cax_ice.get_xticklabels(),fontsize=6)
    cax_title = ppl.setp(cax_ice,title=CloudProduct.cbarTitle_ice)
    ppl.setp(cax_title,fontsize=9)

    #
    # Add a small globe with the swath indicated on it #
    #

    # Create main axes instance, leaving room for colorbar at bottom,
    # and also get the Bbox of the axes instance
    glax_rect = [0.81, 0.75, 0.18, 0.20 ] # [left,bottom,width,height]
    glax = fig.add_axes(glax_rect)

    m_globe = Basemap(lat_0=0.,lon_0=0.,\
        ax=glax,resolution='c',area_thresh=10000.,projection='robin')

    # If we previously had a zero size data array, increase the pointSize
    # so the data points are visible on the global plot
    if (np.shape(gridLon)[0]==0) :
        print "We have no valid data, synthesising dummy data..."
        gridLat = np.array([lat_0,lat_0])
        gridLon = np.array([lon_0,lon_0])
        pointSize = 5.

    x,y=m_globe(np.array(gridLon),np.array(gridLat))
    swath = np.zeros(np.shape(x),dtype=int)

    m_globe.drawcoastlines(ax=glax,linewidth=0.1)
    m_globe.fillcontinents(ax=glax,color='gray',zorder=0)
    m_globe.drawmapboundary(linewidth=0.1)

    p_globe = m_globe.scatter(x,y,s=pointSize,c="red",axes=glax,edgecolors='none')

    # Globe axis title
    glax_xlabel = ppl.setp(glax,xlabel=titleStr)
    ppl.setp(glax_xlabel,fontsize=6)
    
    # Redraw the figure
    canvas.draw()

    # save image 
    print "Writing file to ",outFileName
    canvas.print_figure(outFileName,dpi=dpi)

###################################################
#                  Main Function                  #
###################################################

def main():

    prodChoices=['VCM','VCP','COT','COT_EDR','EPS','EPS_EDR','CTT','CTT_EDR','CTH','CTH_EDR','CTP','CTP_EDR','AOT','AOT_EDR','SDR']
    mapRes = ['c','l','i']

    description = \
    '''
    This is a brief description of %prog
    '''
    usage = "usage: %prog [mandatory args] [options]"
    version = version="%prog NCT3"
    parser = optparse.OptionParser(description=description,usage=usage,version=version)

    # Mandatory arguments
    mandatoryGroup = optparse.OptionGroup(parser, "Mandatory Arguments",
                        "At a minimum these arguments must be specified")

    mandatoryGroup.add_option('-g','--geo_file',
                      action="store",
                      dest="geoFile" ,
                      type="string",
                      help="The full path of the VIIRS geolocation file")
    mandatoryGroup.add_option('-i','--ip_file',
                      action="store",
                      dest="ipFile",
                      type="string",
                      help="The full path of the VIIRS IP/SDR file")
    mandatoryGroup.add_option('-p','--product',
                      action="store",
                      dest="ipProd",
                      type="choice",
                      choices=prodChoices,
                      help='''The VIIRS IP/SDR.\n\n
                           Possible values are...
                           %s
                           ''' % (prodChoices.__str__()[1:-1]))

    parser.add_option_group(mandatoryGroup)

    # Optional arguments
    optionalGroup = optparse.OptionGroup(parser, "Extra Options",
                        "These options may be used to customize plot characteristics.")

    optionalGroup.add_option('-r','--svn_revision',
                      action="store",
                      dest="svnRevision",
                      default=string.split(__version__," ")[2],
                      type="string",
                      help="The Subversion revision number/tag of this script")
    optionalGroup.add_option('--radiance',
                      action="store_true",
                      dest="isRadiance",
                      help="Show radiance for the VIIRS SDR.")
    optionalGroup.add_option('-R','--svn_repo_path',
                      action="store",
                      dest="svnRepoPath",
                      default="https://svn.ssec.wisc.edu/repos/geoffc/Python/VIIRS/"+path.basename(sys.argv[0]),
                      type="string",
                      help="The full Subversion repository path of this script [default is %default].")
    optionalGroup.add_option('--plotMin',
                      action="store",
                      type="float",
                      dest="plotMin",
                      help="Minimum value to plot.")
    optionalGroup.add_option('--plotMax',
                      action="store",
                      type="float",
                      dest="plotMax",
                      help="Maximum value to plot.")
    optionalGroup.add_option('-d','--dpi',
                      action="store",
                      dest="dpi",
                      default='200.',
                      type="float",
                      help="The resolution in dots per inch of the output png file. [default: %default]")
    optionalGroup.add_option('-s','--scale',
                      action="store",
                      dest="scale",
                      default='1.3',
                      type="float",
                      help="The extent of the plot viewport as a proportion of the VIIRS swath width (2400 km). [default: %default]")
    optionalGroup.add_option('-S','--stride',
                      action="store",
                      dest="stride",
                      #default='1',
                      type="int",
                      help="Sample every STRIDE pixels in the VIIRS IP/SDR product. [default: %default]")
    optionalGroup.add_option('-P','--pointSize',
                      action="store",
                      dest="pointSize",
                      #default='0.1',
                      type="float",
                      help="Size of the plot point used to represent each pixel. [default: %default]")
    optionalGroup.add_option('-m','--map_res',
                      action="store",
                      dest="mapRes",
                      default='c',
                      type="choice",
                      choices=mapRes,
                      help="The map coastline resolution. Possible values are 'c' (coarse),'l' (low) and 'i' (intermediate). [default: %default]")
    optionalGroup.add_option('-a','--map_annotation',
                      action="store",
                      dest="mapAnn",
                      #default='',
                      type="string",
                      help="The map legend describing the dataset being shown. [default: IPPROD]")
    optionalGroup.add_option('-o','--output_file',
                      action="store",
                      dest="outputFile",
                      default="out.png",
                      type="string",
                      help="The full path of the output png file. [default: %default]")

    parser.add_option_group(optionalGroup)

    # Parse the arguments from the command line
    (options, args) = parser.parse_args()

    # Copy the Product name to the map annotation
    if options.mapAnn == None :
        options.mapAnn = options.ipProd
    else :
        pass
        options.mapAnn = string.replace(options.mapAnn,'\\n','\n')

    # Check that all of the mandatory options are given. If one or more 
    # are missing, print error message and exit...
    mandatories = ['geoFile', 'ipFile','ipProd']
    mand_errors = ["Missing mandatory argument [-g GEOFILE | --geo_file=GEOFILE]",
                   "Missing mandatory argument [-i IPFILE | --ip_file=IPFILE]",
                   "Missing mandatory argument [-p IPPROD | --product=IPPROD]"
                  ]
    isMissingMand = False
    for m,m_err in zip(mandatories,mand_errors):
        if not options.__dict__[m]:
            isMissingMand = True
            print m_err
    if isMissingMand :
        parser.error("Incomplete mandatory arguments, aborting...")

    # Check that the input files actually exist
    if not glob(options.geoFile) :
        parser.error("Geolocation file \n\t%s\ndoes not exist, aborting..." % (options.geoFile))
    if not glob(options.ipFile) :
        parser.error("Product file \n\t%s\ndoes not exist, aborting..." % (options.ipFile))
        
    # We should have everything we need, run the program...

    #prodFileName='''%s\n%s %s''' % (path.basename(options.ipFile),options.svnRepoPath,str(options.svnRevision))
    prodFileName=''
    dataSet = string.lower((options.ipProd).lstrip())
    mapRes = str(options.mapRes)
    vmin = options.plotMin
    vmax = options.plotMax

    CloudData = ViirsData.CloudProdData.CloudProd()
    cloud_cmap = CloudData.cmap_ice_water
    CloudProdData = ViirsData.CloudProdData.CloudProduct

    # Some defaults plot values if the are not specified on the command line...

    stride_IP = 1
    pointSize_IP = 0.1

    stride_EDR = 1
    pointSize_EDR = 0.30

    stride_SDR_DNB = 7
    pointSize_SDR_DNB = 0.2

    stride_SDR_M = 1
    pointSize_SDR_M = 0.1

    stride_SDR_I = 3
    pointSize_SDR_I = 0.05

    global cmByte,cmBit 

    # Obtain matched lists of the geolocation and product files

    geoList,prodList = granuleFiles(options.geoFile,options.ipFile)

    # Call the granulation and plotting routine for the desired product...

    if 'VCM' in options.ipProd  or 'VCP' in options.ipProd:
        print "Calling VCM ingester..."
        stride = stride_IP if options.stride==None else options.stride
        if 'VCM' in options.ipProd :
            set_vcm_dset(0,1)
        if 'VCP' in options.ipProd :
            set_vcm_dset(5,0)

        lats,lons,vcmData,lat_0,lon_0 = gran_VCM(geoList,prodList,shrink=stride)

        print "Calling VCM plotter..."
        pointSize = pointSize_IP if options.pointSize==None else options.pointSize
        orthoPlot_VCM(lats,lons,vcmData,lat_0=lat_0,lon_0=lon_0,\
            pointSize=pointSize,scale=options.scale,mapRes=mapRes,\
            prodFileName=prodFileName,outFileName=options.outputFile,dpi=options.dpi,titleStr=options.mapAnn)

    if 'COT_EDR' in options.ipProd :
        print "Calling COT EDR ingester..."
        stride = stride_EDR if options.stride==None else options.stride
        lats,lons,cotData,cotPhase,lat_0,lon_0 = gran_COT_EDR([options.geoFile],[options.ipFile],shrink=stride)
        dset=string.split(dataSet,'_')[0]
        print "Calling COT plotter...",dset
        pointSize = pointSize_EDR if options.pointSize==None else options.pointSize
        orthoPlot_COP(lats,lons,cotData,cotPhase,dset,lat_0=lat_0,lon_0=lon_0,\
            pointSize=pointSize,scale=options.scale,mapRes=mapRes,\
            prodFileName=prodFileName,outFileName=options.outputFile,dpi=options.dpi,titleStr=options.mapAnn)

    elif 'EPS_EDR' in options.ipProd :
        print "Calling EPS EDR ingester..."
        stride = stride_EDR if options.stride==None else options.stride
        lats,lons,epsData,epsPhase,lat_0,lon_0 = gran_EPS_EDR([options.geoFile],[options.ipFile],shrink=stride)
        dset=string.split(dataSet,'_')[0]
        print "Calling EPS plotter...",dset
        pointSize = pointSize_EDR if options.pointSize==None else options.pointSize
        orthoPlot_COP(lats,lons,epsData,epsPhase,dset,lat_0=lat_0,lon_0=lon_0,\
            pointSize=pointSize,scale=options.scale,mapRes=mapRes,\
            prodFileName=prodFileName,outFileName=options.outputFile,dpi=options.dpi,titleStr=options.mapAnn)

    elif 'COT' in options.ipProd or 'EPS' in options.ipProd :
        print "Calling COP ingester..."
        stride = stride_IP if options.stride==None else options.stride
        lats,lons,copData,copPhase,lat_0,lon_0 = gran_COP([options.geoFile],[options.ipFile],\
            dataSet,shrink=stride)
        print "Calling COP plotter..."
        pointSize = pointSize_IP if options.pointSize==None else options.pointSize
        orthoPlot_COP(lats,lons,copData,copPhase,dataSet,lat_0=lat_0,lon_0=lon_0,\
            abScale='log',pointSize=pointSize,scale=options.scale,mapRes=mapRes,\
            prodFileName=prodFileName,outFileName=options.outputFile,dpi=options.dpi,titleStr=options.mapAnn)

    if 'CTT_EDR' in options.ipProd :
        print "Calling CTT EDR ingester..."
        stride = stride_EDR if options.stride==None else options.stride
        lats,lons,cttData,cttPhase,lat_0,lon_0 = gran_CTT_EDR([options.geoFile],[options.ipFile],shrink=stride)
        dset=string.split(dataSet,'_')[0]
        print "Calling CTT plotter...",dset
        pointSize = pointSize_EDR if options.pointSize==None else options.pointSize
        orthoPlot_CTp(lats,lons,cttData,cttPhase,dset,lat_0=lat_0,lon_0=lon_0,\
            pointSize=pointSize,scale=options.scale,mapRes=mapRes,\
            prodFileName=prodFileName,outFileName=options.outputFile,dpi=options.dpi,titleStr=options.mapAnn)

    elif 'CTP_EDR' in options.ipProd :
        print "Calling CTP EDR ingester..."
        stride = stride_EDR if options.stride==None else options.stride
        lats,lons,ctpData,ctpPhase,lat_0,lon_0 = gran_CTP_EDR([options.geoFile],[options.ipFile],shrink=stride)
        dset=string.split(dataSet,'_')[0]
        print "Calling CTP plotter...",dset
        pointSize = pointSize_EDR if options.pointSize==None else options.pointSize
        orthoPlot_CTp(lats,lons,ctpData,ctpPhase,dset,lat_0=lat_0,lon_0=lon_0,\
            pointSize=pointSize,scale=options.scale,mapRes=mapRes,\
            prodFileName=prodFileName,outFileName=options.outputFile,dpi=options.dpi,titleStr=options.mapAnn)

    elif 'CTH_EDR' in options.ipProd :
        print "Calling CTH EDR ingester..."
        stride = stride_EDR if options.stride==None else options.stride
        lats,lons,cthData,cthPhase,lat_0,lon_0 = gran_CTH_EDR([options.geoFile],[options.ipFile],shrink=stride)
        dset=string.split(dataSet,'_')[0]
        print "Calling CTH plotter...",dset
        pointSize = pointSize_EDR if options.pointSize==None else options.pointSize
        orthoPlot_CTp(lats,lons,cthData,cthPhase,dset,lat_0=lat_0,lon_0=lon_0,\
            pointSize=pointSize,scale=options.scale,mapRes=mapRes,\
            prodFileName=prodFileName,outFileName=options.outputFile,dpi=options.dpi,titleStr=options.mapAnn)

    elif 'CTT' in options.ipProd or 'CTP' in options.ipProd or 'CTH' in options.ipProd :
        print "Calling CTp ingester..."
        stride = stride_IP if options.stride==None else options.stride
        lats,lons,ctpData,ctpPhase,lat_0,lon_0 = gran_CTp([options.geoFile],[options.ipFile],\
            dataSet,shrink=stride)
        print "Calling CTp plotter..."
        pointSize = pointSize_IP if options.pointSize==None else options.pointSize
        orthoPlot_CTp(lats,lons,ctpData,ctpPhase,dataSet,lat_0=lat_0,lon_0=lon_0,\
            pointSize=pointSize,scale=options.scale,mapRes=mapRes,\
            prodFileName=prodFileName,outFileName=options.outputFile,dpi=options.dpi,titleStr=options.mapAnn)

    if 'AOT_EDR' in options.ipProd :
        print "Calling AOT EDR ingester..."
        stride = stride_EDR if options.stride==None else options.stride
        vmin = -0.05 if (vmin==None) else vmin
        vmax = 0.8 if (vmax==None) else vmax
        lats,lons,aotData,lat_0,lon_0,ModeGran = gran_AOT_EDR([options.geoFile],[options.ipFile],shrink=stride)
        print "Calling AOT plotter..."
        pointSize = pointSize_EDR if options.pointSize==None else options.pointSize
        orthoPlot_AOT(lats,lons,aotData,ModeGran,lat_0=lat_0,lon_0=lon_0,vmin=vmin,vmax=vmax,\
            pointSize=pointSize,scale=options.scale,mapRes=mapRes,cmap=cloud_cmap, \
            prodFileName=prodFileName,outFileName=options.outputFile,dpi=options.dpi,titleStr=options.mapAnn)
    elif 'AOT' in options.ipProd :
        print "Calling AOT ingester..."
        stride = stride_IP if options.stride==None else options.stride
        vmin = -0.05 if (vmin==None) else vmin
        vmax = 0.8 if (vmax==None) else vmax

        lats,lons,aotData,lat_0,lon_0,ModeGran = gran_AOT(geoList,prodList,shrink=stride)
        
        print "Calling AOT plotter..."
        pointSize = pointSize_IP if options.pointSize==None else options.pointSize
        orthoPlot_AOT(lats,lons,aotData,ModeGran,lat_0=lat_0,lon_0=lon_0,vmin=vmin,vmax=vmax, \
            pointSize=pointSize,scale=options.scale,mapRes=mapRes,cmap=cloud_cmap, \
            prodFileName=prodFileName,outFileName=options.outputFile,dpi=options.dpi,titleStr=options.mapAnn)

    print "Exiting..."
    sys.exit(0)

if __name__ == '__main__':
    main()


