#!/usr/bin/env python
# encoding: utf-8
"""
ql_viirs_edr.py

Purpose: Create quicklook PNGs for VIIRS EDR products.

Minimum commandline...

export CSPP_EDR_HOME=$(readlink -f /path/to/EDR)
source $CSPP_EDR_HOME/cspp_edr_env.sh

python ql_viirs_edr.py -g '/path/to/files/GMTCO*.h5' -i '/path/to/files/IICMO*.h5' -p VCM

  or

python ql_viirs_edr.py --geo_file=/path/to/files/GMTCO*.h5 --ip_file=/path/to/files/IICMO*.h5 -p VCM


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

import string, sys, traceback
from glob import glob
from os import path, uname
from time import time
import traceback

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

from ViirsData import ViirsTrimTable
import viirs_edr_data
import viirs_cloud_mask as viirsCM
import viirs_cloud_products as viirsCld
import viirs_aerosol_products as viirsAero
import viirs_sst_products as viirsSST
import viirs_vi_products as viirsVI
import viirs_st_products as viirsST
import viirs_lst_products as viirsLST

import tables as pytables
from tables import exceptions as pyEx

# every module should have a LOG object
# e.g. LOG.warning('my dog has fleas')
import logging
LOG = logging.getLogger(__file__)

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

    LOG.debug("Initial geoGlob = ",geoGlob)
    LOG.debug("Initial prodGlob = ",prodGlob)

    geoGlob = path.basename(path.abspath(path.expanduser(geoGlob)))
    prodGlob = path.basename(path.abspath(path.expanduser(prodGlob)))

    geoPrefix = string.split(geoGlob,'_')[0]

    LOG.debug("geoDir = ",geoDir)
    LOG.debug("prodDir = ",prodDir)
    LOG.debug("geoGlob = ",geoGlob)
    LOG.debug("prodGlob = ",prodGlob)
    LOG.debug("geoPrefix = ",geoPrefix)

    geoList_in = glob("%s/%s" % (geoDir,geoGlob))
    prodList_in = glob("%s/%s" % (prodDir,prodGlob))
    geoList_in.sort()
    prodList_in.sort()
    
    #LOG.debug("prodList_in...")
    #for prodFile in prodList_in:
        #LOG.debug(prodFile)

    geoList = []
    prodList = []
    #prodList = prodList_in
    for files in prodList_in :
        prod_arr = string.split(path.basename(files),"_")
        #LOG.debug("prod_arr = ",prod_arr
        dateStamp = prod_arr[2]
        timeStamp = prod_arr[3]
        geoFileGlob="%s/%s*%s_%s*.h5" % (geoDir,geoPrefix,dateStamp,timeStamp)
        #LOG.debug("dateStamp = ",dateStamp)
        #LOG.debug("timeStamp = ",timeStamp)
        #LOG.debug("geoFileGlob = ",geoFileGlob)
        geoFile = glob("%s/%s*%s_%s*.h5" % (geoDir,geoPrefix,dateStamp,timeStamp))
        if (np.shape(geoFile)[0] != 0) :
            geoList.append(geoFile[0])
            prodList.append(files)
        else :
            #geoList.append(files)
            LOG.debug(" ... no match found for {}, appending {}".format( geoFile, files))
            pass
    
    for geoFile,prodFile in zip(geoList,prodList):
        LOG.debug("{},{}".format(geoFile,prodFile))
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
    LOG.debug("Importing the CloudMaskProduct object...")
    CloudMaskProduct = viirs_edr_data.CloudMaskData
    LOG.debug("done")

    try :
        reload(viirsCM)
        reload(viirs_edr_data)
        del(viirsCMobj)
        del(latArr)
        del(lonArr)
        del(cmArr)
        del(qualArr)
    except :
        pass

    LOG.debug("Creating viirsCMobj...")
    viirsCMobj = viirsCM.viirsCM()
    LOG.debug("done")

    # Determine the correct fillValue
    trimObj = ViirsTrimTable()
    eps = 1.e-6

    for grans in np.arange(len(geoList)):

        LOG.debug("\nIngesting granule {} using cmByte={} and cmBit={} ...".format(grans,cmByte,cmBit))
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

        LOG.debug("Intermediate latArr.shape = {}".format(str(latArr.shape)))
        LOG.debug("Intermediate lonArr.shape = {}".format(str(lonArr.shape)))
        LOG.debug("Intermediate cmArr.shape = {}".format(str(cmArr.shape)))
        LOG.debug("Intermediate qualArr.shape = {}".format(str(qualArr.shape)))

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
                LOG.debug("Dataset was neither int not float... a worry")
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
            LOG.debug(">> error: Mask Error, probably mismatched geolocation and product array sizes, aborting...")
            sys.exit(1)

    except Exception, err :
        LOG.debug(">> error: {}...".format(str(err)))
        sys.exit(1)

    return lats,lons,data,lat_0,lon_0

def gran_COP(geoList,copList,dataSet,shrink=1):
    '''
    Returns the granulated COP dataset
    '''

    try :
        reload(viirsCld)
        reload(viirs_edr_data)
        del(viirsCldObj)
        del(latArr)
        del(lonArr)
        del(copArr)
        del(phaseArr)
    except :
        pass

    LOG.debug("Creating viirsCldObj...")
    reload(viirsCld)
    viirsCldObj = viirsCld.viirsCld()
    LOG.debug("done")

    # Determine the correct fillValue
    trimObj = ViirsTrimTable()
    eps = 1.e-6

    for grans in np.arange(len(geoList)):
        LOG.debug("Ingesting granule {} ...".format(grans))
        retArr = viirsCldObj.ingestLite(geoList[grans],copList[grans],dataSet,1)
        LOG.debug("done")

        try :
            latArr  = viirsCldObj.Lat[:,:]
            lonArr  = viirsCldObj.Lon[:,:]
            
            lat_0 = latArr[np.shape(latArr)[0]/2,np.shape(latArr)[1]/2]
            lon_0 = lonArr[np.shape(lonArr)[0]/2,np.shape(lonArr)[1]/2]

            badGeo = False
            if not (-90. <= lat_0 <= 90.) :
                LOG.debug("\n>> error: Latitude of granule midpoint ({}) does not satisfy (-90. <= lat_0 <= 90.)\nfor file {}\n\taborting...".format(lat_0,geoList[grans]))
                badGeo = True
            if not (-180. <= lat_0 <= 180.) :
                LOG.debug("\n>> error: Longitude of granule midpoint ({}) does not satisfy (-180. <= lon_0 <= 180.)\nfor file {}\n\taborting...".format(lon_0,geoList[grans]))
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
                    LOG.debug("Dataset was neither int not float... a worry")
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
                LOG.debug(">> error: Mask Error, probably mismatched geolocation and product array sizes, aborting...")
                sys.exit(1)

        except :
            LOG.debug(">> error: there was an exception...")
            sys.exit(1)

    return lats,lons,data,phase,lat_0,lon_0

def gran_COT_EDR(geoList,cotList,shrink=1):
    '''
    Returns the granulated COT EDR
    '''

    try :
        reload(viirsAero)
        reload(viirs_edr_data)
        del(viirsAeroObj)
        del(latArr)
        del(lonArr)
        del(cotArr)
    except :
        pass

    LOG.debug("Creating viirsCldObj...")
    reload(viirsCld)
    viirsCldObj = viirsCld.viirsCld()
    LOG.debug("done")

    # Determine the correct fillValue
    trimObj = ViirsTrimTable()
    eps = 1.e-6

    for grans in np.arange(len(geoList)):
        LOG.debug("Ingesting EDR granule {} ...".format(grans))

        # Read in geolocation...
        ViirsGeoFileName = geoList[grans]
        if(pytables.isHDF5File(ViirsGeoFileName)) :
            try :
                ViirsGeoFileObj = pytables.openFile(ViirsGeoFileName,mode='r')
                LOG.debug("Successfully opened geolocation file",ViirsGeoFileName)
            except IOError :
                LOG.debug(">> error: Could not open geolocation file: ",ViirsGeoFileName)
                sys.exit(1)

            # Detemine the geolocation group name and related information
            #group = getobj(ViirsGeoFileObj,'/All_Data')
            group = ViirsGeoFileObj.getNode('/All_Data')
            geoGroupName = '/All_Data/'+group.__members__[0]
            group._g_close()
            LOG.debug("Geolocation Group : {} ".format(geoGroupName))
            isEdrGeo = ('VIIRS-CLD-AGG-GEO_All' in geoGroupName)
            if not isEdrGeo :
                LOG.debug(">> error: {} is not an EDR resolution cloud geolocation file\n\taborting...".format(ViirsGeoFileName))
                sys.exit(1)

            # Get the geolocation Latitude and Longitude nodes...
            try :
                dataName = 'Latitude'
                geoLatNode = ViirsGeoFileObj.getNode(geoGroupName+'/'+dataName)
                LOG.debug("Shape of {} node is {}".format(geoGroupName+'/'+dataName,repr(geoLatNode.shape)))
                dataName = 'Longitude'
                geoLonNode = ViirsGeoFileObj.getNode(geoGroupName+'/'+dataName)
                LOG.debug("Shape of {} node is {}".format(geoGroupName+'/'+dataName,repr(geoLonNode.shape)))
            except pyEx.NoSuchNodeError :
                LOG.debug("\n>> error: No required node {}/{} in {}\n\taborting...".format(geoGroupName,dataName,ViirsGeoFileName))
                ViirsGeoFileObj.close()
                sys.exit(1)

            # Make a copy of the geolocation node data, so we can close the file
            try :
                dataName = 'Latitude'
                LOG.debug("Reading {} dataset...".format(dataName))
                latArr = geoLatNode[:,:]
                geoLatNode.close()
                latArr = np.squeeze(latArr)
                LOG.debug("done")
                dataName = 'Longitude'
                LOG.debug("Reading {} dataset...".format(dataName))
                lonArr = geoLonNode[:,:]
                geoLonNode.close()
                lonArr = np.squeeze(lonArr)
                LOG.debug("done")
                LOG.debug("Closing geolocation file")
                ViirsGeoFileObj.close()
            except :
                LOG.debug("\n>> error: Could not retrieve %/% node data in {}\n\taborting...".format(geoGroupName,dataName,ViirsGeoFileName))
                geoLatNode.close()
                geoLonNode.close()
                ViirsGeoFileObj.close()
                sys.exit(1)

        else :
            LOG.debug("\n>> error: {} is not a HDF5 file,\n\taborting...".format(ViirsGeoFileName))
            sys.exit(1)
        
        # Read in dataSets...
        ViirsEDRFileName = cotList[grans]
        if(pytables.isHDF5File(ViirsEDRFileName)) :
            try :
                ViirsEDRFileObj = pytables.openFile(ViirsEDRFileName,mode='r')
                LOG.debug("Successfully opened edr file",ViirsEDRFileName)
            except IOError :
                LOG.debug(">> error: Could not open edr file: ",ViirsEDRFileName)
                sys.exit(1)

            # Detemine the edr group name and related information
            #group = getobj(ViirsEDRFileObj,'/All_Data')
            group = ViirsEDRFileObj.getNode('/All_Data')
            edrGroupName = '/All_Data/'+group.__members__[0]
            group._g_close()
            LOG.debug("Edr Group : {} ".format(edrGroupName))
            isEdr = ('VIIRS-COT-EDR_All' in edrGroupName)
            if not isEdr :
                LOG.debug(">> error: {} is not an EDR resolution cloud file\n\taborting...".format(ViirsEDRFileName))
                sys.exit(1)

            # Get the edr nodes...
            try :
                dataName = 'AverageCloudOpticalThickness'
                cotNode = ViirsEDRFileObj.getNode(edrGroupName+'/'+dataName)
                LOG.debug("Shape of {} node is {}".format(edrGroupName+'/'+dataName,repr(cotNode.shape)))
                dataName = 'COTFactors'
                cotFactorsNode = ViirsEDRFileObj.getNode(edrGroupName+'/'+dataName)
                LOG.debug("Shape of {} node is {}".format(edrGroupName+'/'+dataName,repr(cotFactorsNode.shape)))
                dataName = 'QF3_VIIRSCOTAVGEDR'
                QF3_VIIRSCOTAVGEDR_Node = ViirsEDRFileObj.getNode(edrGroupName+'/'+dataName)
                LOG.debug("Shape of {} node is {}".format(edrGroupName+'/'+dataName,repr(QF3_VIIRSCOTAVGEDR_Node.shape)))
            except pyEx.NoSuchNodeError :
                LOG.debug("\n>> error: No required node {}/{} in {}\n\taborting...".format(edrGroupName,dataName,ViirsEDRFileName))
                ViirsEDRFileObj.close()
                sys.exit(1)

            # Make a copy of the edr node data, so we can close the file
            try :
                dataName = 'AverageCloudOpticalThickness'
                LOG.debug("Reading {} dataset...".format(dataName))
                cot = cotNode[:,:]
                cot = np.squeeze(cot)
                cotNode.close()
                LOG.debug("done")
                dataName = 'COTFactors'
                LOG.debug("Reading {} dataset...".format(dataName))
                cotFactors = cotFactorsNode[:]
                cotFactors = np.squeeze(cotFactors)
                cotFactorsNode.close()
                LOG.debug("done")
                dataName = 'QF3_VIIRSCOTAVGEDR'
                LOG.debug("Reading {} dataset...".format(dataName))
                cloudConf = np.bitwise_and(QF3_VIIRSCOTAVGEDR_Node[:,:],3) >> 0
                waterCldFrac = np.bitwise_and(QF3_VIIRSCOTAVGEDR_Node[:,:],12) >> 2
                mixedCldFrac = np.bitwise_and(QF3_VIIRSCOTAVGEDR_Node[:,:],192) >> 6
                QF3_VIIRSCOTAVGEDR_Node.close()
                LOG.debug("done")
                LOG.debug("Closing edr file")
                ViirsEDRFileObj.close()

                LOG.debug("Shape of cot is {}".format(repr(np.shape(cot))))
                LOG.debug("Shape of latsArr is {}".format(repr(np.shape(latArr))))
                LOG.debug("Shape of lonsArr is {}".format(repr(np.shape(lonArr))))

            except :
                LOG.debug("\n>> error: Could not retrieve %/% node data in {}\n\taborting...".format(edrGroupName,dataName,ViirsEDRFileName))
                cotNode.close()
                cotFactorsNode.close()
                QF3_VIIRSCOTAVGEDR_Node.close()
                ViirsEDRFileObj.close()
                sys.exit(1)

        else :
            LOG.debug("\n>> error: {} is not a HDF5 file,\n\taborting...".format(ViirsEDRFileName))
            sys.exit(1)
        
        LOG.debug("Creating some masks")
        try :
            
            lat_0 = latArr[np.shape(latArr)[0]/2,np.shape(latArr)[1]/2]
            lon_0 = lonArr[np.shape(lonArr)[0]/2,np.shape(lonArr)[1]/2]

            badGeo = False
            if not (-90. <= lat_0 <= 90.) :
                LOG.debug("\n>> error: Latitude of granule midpoint ({}) does not satisfy (-90. <= lat_0 <= 90.)\nfor file {}\n\taborting...".format(lat_0,geoList[grans]))
                badGeo = True
            if not (-180. <= lat_0 <= 180.) :
                LOG.debug("\n>> error: Longitude of granule midpoint ({}) does not satisfy (-180. <= lon_0 <= 180.)\nfor file {}\n\taborting...".format(lon_0,geoList[grans]))
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
                    LOG.debug("Dataset was neither int not float... a worry")
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
                LOG.debug(">> error: Mask Error, probably mismatched geolocation and product array sizes, aborting...")
                sys.exit(1)

            data = data * cotFactors[0] + cotFactors[1]

            cldPhase = np.ones(data.shape,dtype=np.uint8) * trimObj.sdrTypeFill['NA_FILL']['uint8']

            LOG.debug("cloudPhase.shape = ",cldPhase.shape)
            LOG.debug("cloudConfArr.shape = ",cloudConfArr.shape)
            LOG.debug("waterCldFracArr.shape = ",waterCldFracArr.shape)
            LOG.debug("mixedCldFracArr.shape = ",mixedCldFracArr.shape)

            cloudFracMask = (cloudConfArr>0)
            waterFracMask = (waterCldFracArr>0)
            mixedFracMask = (mixedCldFracArr>0)

            waterMask = cloudFracMask * (waterFracMask + mixedFracMask)
            iceMask = cloudFracMask * np.logical_not(waterFracMask + mixedFracMask)

            LOG.debug("waterMask.shape = ",waterMask.shape)
            LOG.debug("iceMask.shape = ",iceMask.shape)

            LOG.debug("Number of pixels with > 25% cloudy confidence: ",cloudFracMask.sum())
            LOG.debug("Number of pixels with > 25% water cloud fraction: ",waterFracMask.sum())
            LOG.debug("Number of pixels with > 25% mixed cloud fraction: ",mixedFracMask.sum())

            LOG.debug("Number of water pixels: ",(cloudFracMask * (waterFracMask + mixedFracMask)).sum())
            LOG.debug("Number of ice pixels: ",cloudFracMask.sum() - (waterFracMask + mixedFracMask).sum())

            LOG.debug("Number of water pixels: ",waterMask.sum())
            LOG.debug("Number of ice pixels: ",iceMask.sum())

            cldPhaseShape = cldPhase.shape
            cldPhase = np.ravel(cldPhase)

            wtrIdx = np.where(np.ravel(waterMask))
            iceIdx = np.where(np.ravel(iceMask))
            cldPhase[wtrIdx] = 3
            cldPhase[iceIdx] = 5

            cldPhase = np.reshape(cldPhase,cldPhaseShape)

        except Exception, err :
            LOG.debug(">> error: {}...".format(str(err)))
            sys.exit(1)

    LOG.debug("lat_0,lon_0 = {},{}".format(lat_0,lon_0))
    LOG.debug("Shape of data is {}".format(repr(np.shape(data))))
    LOG.debug("Shape of cldPhase is {}".format(repr(np.shape(cldPhase))))
    LOG.debug("Shape of lats is {}".format(repr(np.shape(lats))))
    LOG.debug("Shape of lons is {}".format(repr(np.shape(lons))))

    return lats,lons,data,cldPhase,lat_0,lon_0

def gran_EPS_EDR(geoList,epsList,shrink=1):
    '''
    Returns the granulated EPS EDR
    '''

    try :
        reload(viirsCld)
        reload(viirs_edr_data)
        del(viirsCldObj)
        del(latArr)
        del(lonArr)
        del(epsArr)
    except :
        pass

    LOG.debug("Creating viirsCldObj...")
    reload(viirsCld)
    viirsCldObj = viirsCld.viirsCld()
    LOG.debug("done")

    # Determine the correct fillValue
    trimObj = ViirsTrimTable()
    eps = 1.e-6

    for grans in np.arange(len(geoList)):
        LOG.debug("Ingesting EDR granule {} ...".format(grans))

        # Read in geolocation...
        ViirsGeoFileName = geoList[grans]
        if(pytables.isHDF5File(ViirsGeoFileName)) :
            try :
                ViirsGeoFileObj = pytables.openFile(ViirsGeoFileName,mode='r')
                LOG.debug("Successfully opened geolocation file",ViirsGeoFileName)
            except IOError :
                LOG.debug(">> error: Could not open geolocation file: ",ViirsGeoFileName)
                sys.exit(1)

            # Detemine the geolocation group name and related information
            #group = getobj(ViirsGeoFileObj,'/All_Data')
            group = ViirsGeoFileObj.getNode('/All_Data')
            geoGroupName = '/All_Data/'+group.__members__[0]
            group._g_close()
            LOG.debug("Geolocation Group : {} ".format(geoGroupName))
            isEdrGeo = ('VIIRS-CLD-AGG-GEO_All' in geoGroupName)
            if not isEdrGeo :
                LOG.debug(">> error: {} is not an EDR resolution cloud geolocation file\n\taborting...".format(ViirsGeoFileName))
                sys.exit(1)

            # Get the geolocation Latitude and Longitude nodes...
            try :
                dataName = 'Latitude'
                geoLatNode = ViirsGeoFileObj.getNode(geoGroupName+'/'+dataName)
                LOG.debug("Shape of {} node is {}".format(geoGroupName+'/'+dataName,repr(geoLatNode.shape)))
                dataName = 'Longitude'
                geoLonNode = ViirsGeoFileObj.getNode(geoGroupName+'/'+dataName)
                LOG.debug("Shape of {} node is {}".format(geoGroupName+'/'+dataName,repr(geoLonNode.shape)))
            except pyEx.NoSuchNodeError :
                LOG.debug("\n>> error: No required node {}/{} in {}\n\taborting...".format(geoGroupName,dataName,ViirsGeoFileName))
                ViirsGeoFileObj.close()
                sys.exit(1)

            # Make a copy of the geolocation node data, so we can close the file
            try :
                dataName = 'Latitude'
                LOG.debug("Reading {} dataset...".format(dataName))
                latArr = geoLatNode[:,:]
                geoLatNode.close()
                latArr = np.squeeze(latArr)
                LOG.debug("done")
                dataName = 'Longitude'
                LOG.debug("Reading {} dataset...".format(dataName))
                lonArr = geoLonNode[:,:]
                geoLonNode.close()
                lonArr = np.squeeze(lonArr)
                LOG.debug("done")
                LOG.debug("Closing geolocation file")
                ViirsGeoFileObj.close()
            except :
                LOG.debug("\n>> error: Could not retrieve %/% node data in {}\n\taborting...".format(geoGroupName,dataName,ViirsGeoFileName))
                geoLatNode.close()
                geoLonNode.close()
                ViirsGeoFileObj.close()
                sys.exit(1)

        else :
            LOG.debug("\n>> error: {} is not a HDF5 file,\n\taborting...".format(ViirsGeoFileName))
            sys.exit(1)
        
        # Read in dataSets...
        ViirsEDRFileName = epsList[grans]
        if(pytables.isHDF5File(ViirsEDRFileName)) :
            try :
                ViirsEDRFileObj = pytables.openFile(ViirsEDRFileName,mode='r')
                LOG.debug("Successfully opened edr file",ViirsEDRFileName)
            except IOError :
                LOG.debug(">> error: Could not open edr file: ",ViirsEDRFileName)
                sys.exit(1)

            # Detemine the edr group name and related information
            #group = getobj(ViirsEDRFileObj,'/All_Data')
            group = ViirsEDRFileObj.getNode('/All_Data')
            edrGroupName = '/All_Data/'+group.__members__[0]
            group._g_close()
            LOG.debug("Edr Group : {} ".format(edrGroupName))
            isEdr = ('VIIRS-CEPS-EDR_All' in edrGroupName)
            if not isEdr :
                LOG.debug(">> error: {} is not an EDR resolution cloud file\n\taborting...".format(ViirsEDRFileName))
                sys.exit(1)

            # Get the edr nodes...
            try :
                dataName = 'AverageCloudEffectiveParticleSize'
                epsNode = ViirsEDRFileObj.getNode(edrGroupName+'/'+dataName)
                LOG.debug("Shape of {} node is {}".format(edrGroupName+'/'+dataName,repr(epsNode.shape)))
                dataName = 'CEPSFactors'
                epsFactorsNode = ViirsEDRFileObj.getNode(edrGroupName+'/'+dataName)
                LOG.debug("Shape of {} node is {}".format(edrGroupName+'/'+dataName,repr(epsFactorsNode.shape)))
                dataName = 'QF3_VIIRSCEPSAVGEDR'
                QF3_VIIRSCEPSAVGEDR_Node = ViirsEDRFileObj.getNode(edrGroupName+'/'+dataName)
                LOG.debug("Shape of {} node is {}".format(edrGroupName+'/'+dataName,repr(QF3_VIIRSCEPSAVGEDR_Node.shape)))
            except pyEx.NoSuchNodeError :
                LOG.debug("\n>> error: No required node {}/{} in {}\n\taborting...".format(edrGroupName,dataName,ViirsEDRFileName))
                ViirsEDRFileObj.close()
                sys.exit(1)

            # Make a copy of the edr node data, so we can close the file
            try :
                dataName = 'AverageCloudEffectiveParticleSize'
                LOG.debug("Reading {} dataset...".format(dataName))
                eps = epsNode[:,:]
                eps = np.squeeze(eps)
                epsNode.close()
                LOG.debug("done")
                dataName = 'CEPSFactors'
                LOG.debug("Reading {} dataset...".format(dataName))
                epsFactors = epsFactorsNode[:]
                epsFactors = np.squeeze(epsFactors)
                epsFactorsNode.close()
                LOG.debug("done")
                dataName = 'QF3_VIIRSCEPSAVGEDR'
                LOG.debug("Reading {} dataset...".format(dataName))
                cloudConf = np.bitwise_and(QF3_VIIRSCEPSAVGEDR_Node[:,:],3) >> 0
                waterCldFrac = np.bitwise_and(QF3_VIIRSCEPSAVGEDR_Node[:,:],12) >> 2
                mixedCldFrac = np.bitwise_and(QF3_VIIRSCEPSAVGEDR_Node[:,:],192) >> 6
                QF3_VIIRSCEPSAVGEDR_Node.close()
                LOG.debug("done")
                LOG.debug("Closing edr file")
                ViirsEDRFileObj.close()

                LOG.debug("Shape of eps is {}".format(repr(np.shape(eps))))
                LOG.debug("Shape of latsArr is {}".format(repr(np.shape(latArr))))
                LOG.debug("Shape of lonsArr is {}".format(repr(np.shape(lonArr))))

            except :
                LOG.debug("\n>> error: Could not retrieve %/% node data in {}\n\taborting...".format(edrGroupName,dataName,ViirsEDRFileName))
                epsNode.close()
                epsFactorsNode.close()
                QF3_VIIRSCEPSAVGEDR_Node.close()
                ViirsEDRFileObj.close()
                sys.exit(1)

        else :
            LOG.debug("\n>> error: {} is not a HDF5 file,\n\taborting...".format(ViirsEDRFileName))
            sys.exit(1)
        
        LOG.debug("Creating some masks")
        try :
            
            lat_0 = latArr[np.shape(latArr)[0]/2,np.shape(latArr)[1]/2]
            lon_0 = lonArr[np.shape(lonArr)[0]/2,np.shape(lonArr)[1]/2]

            badGeo = False
            if not (-90. <= lat_0 <= 90.) :
                LOG.debug("\n>> error: Latitude of granule midpoint ({}) does not satisfy (-90. <= lat_0 <= 90.)\nfor file {}\n\taborting...".format(lat_0,geoList[grans]))
                badGeo = True
            if not (-180. <= lat_0 <= 180.) :
                LOG.debug("\n>> error: Longitude of granule midpoint ({}) does not satisfy (-180. <= lon_0 <= 180.)\nfor file {}\n\taborting...".format(lon_0,geoList[grans]))
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
                    LOG.debug("Dataset was neither int not float... a worry")
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
                LOG.debug(">> error: Mask Error, probably mismatched geolocation and product array sizes, aborting...")
                sys.exit(1)

            data = data * epsFactors[0] + epsFactors[1]

            cldPhase = np.ones(data.shape,dtype=np.uint8) * trimObj.sdrTypeFill['NA_FILL']['uint8']

            cloudFracMask = (cloudConfArr>0)
            waterFracMask = (waterCldFracArr>0)
            mixedFracMask = (mixedCldFracArr>0)

            waterMask = cloudFracMask * (waterFracMask + mixedFracMask)
            iceMask = cloudFracMask * np.logical_not(waterFracMask + mixedFracMask)

            LOG.debug("Number of pixels with > 25% cloudy confidence: ",cloudFracMask.sum())
            LOG.debug("Number of pixels with > 25% water cloud fraction: ",waterFracMask.sum())
            LOG.debug("Number of pixels with > 25% mixed cloud fraction: ",mixedFracMask.sum())

            LOG.debug("Number of water pixels: ",(cloudFracMask * (waterFracMask + mixedFracMask)).sum())
            LOG.debug("Number of ice pixels: ",cloudFracMask.sum() - (waterFracMask + mixedFracMask).sum())

            LOG.debug("Number of water pixels: ",waterMask.sum())
            LOG.debug("Number of ice pixels: ",iceMask.sum())

            cldPhaseShape = cldPhase.shape
            cldPhase = np.ravel(cldPhase)

            wtrIdx = np.where(np.ravel(waterMask))
            iceIdx = np.where(np.ravel(iceMask))
            cldPhase[wtrIdx] = 3
            cldPhase[iceIdx] = 5

            cldPhase = np.reshape(cldPhase,cldPhaseShape)

        except Exception, err :
            LOG.debug(">> error: {}...".format(str(err)))
            sys.exit(1)

    LOG.debug("lat_0,lon_0 = {},{}".format(lat_0,lon_0))
    LOG.debug("Shape of data is {}".format(repr(np.shape(data))))
    LOG.debug("Shape of cldPhase is {}".format(repr(np.shape(cldPhase))))
    LOG.debug("Shape of lats is {}".format(repr(np.shape(lats))))
    LOG.debug("Shape of lons is {}".format(repr(np.shape(lons))))

    return lats,lons,data,cldPhase,lat_0,lon_0

def gran_CTp(geoList,ctpList,dataSet,shrink=1):
    '''
    Returns the granulated CTp dataset
    '''

    try :
        reload(viirsCld)
        reload(viirs_edr_data)
        del(viirsCldObj)
        del(latArr)
        del(lonArr)
        del(ctpArr)
        del(phaseArr)
    except :
        pass

    LOG.debug("Creating viirsCldObj...")
    reload(viirsCld)
    viirsCldObj = viirsCld.viirsCld()
    LOG.debug("done")

    # Determine the correct fillValue
    trimObj = ViirsTrimTable()
    eps = 1.e-6

    for grans in np.arange(len(geoList)):
        LOG.debug("Ingesting granule {} ...".format(grans))
        retArr = viirsCldObj.ingestLite(geoList[grans],ctpList[grans],dataSet,1)
        LOG.debug("done")

        try :
            latArr  = viirsCldObj.Lat[:,:]
            lonArr  = viirsCldObj.Lon[:,:]
            
            lat_0 = latArr[np.shape(latArr)[0]/2,np.shape(latArr)[1]/2]
            lon_0 = lonArr[np.shape(lonArr)[0]/2,np.shape(lonArr)[1]/2]

            badGeo = False
            if not (-90. <= lat_0 <= 90.) :
                LOG.debug("\n>> error: Latitude of granule midpoint ({}) does not satisfy (-90. <= lat_0 <= 90.)\nfor file {}\n\taborting...".format(lat_0,geoList[grans]))
                badGeo = True
            if not (-180. <= lat_0 <= 180.) :
                LOG.debug("\n>> error: Longitude of granule midpoint ({}) does not satisfy (-180. <= lon_0 <= 180.)\nfor file {}\n\taborting...".format(lon_0,geoList[grans]))
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
                    LOG.debug("Dataset was neither int not float... a worry")
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
                LOG.debug(">> error: Mask Error, probably mismatched geolocation and product array sizes, aborting...")
                sys.exit(1)

        except :
            LOG.debug(">> error: There was an exception...")
            sys.exit(1)

    return lats,lons,data,phase,lat_0,lon_0

def gran_CTT_EDR(geoList,cttList,shrink=1):
    '''
    Returns the granulated EPS EDR
    '''

    try :
        reload(viirsCld)
        reload(viirs_edr_data)
        del(viirsCldObj)
        del(latArr)
        del(lonArr)
        del(cttArr)
    except :
        pass

    LOG.debug("Creating viirsCldObj...")
    reload(viirsCld)
    viirsCldObj = viirsCld.viirsCld()
    LOG.debug("done")

    # Determine the correct fillValue
    trimObj = ViirsTrimTable()
    eps = 1.e-6

    for grans in np.arange(len(geoList)):
        LOG.debug("Ingesting EDR granule {} ...".format(grans))

        # Read in geolocation...
        ViirsGeoFileName = geoList[grans]
        if(pytables.isHDF5File(ViirsGeoFileName)) :
            try :
                ViirsGeoFileObj = pytables.openFile(ViirsGeoFileName,mode='r')
                LOG.debug("Successfully opened geolocation file",ViirsGeoFileName)
            except IOError :
                LOG.debug(">> error: Could not open geolocation file: ",ViirsGeoFileName)
                sys.exit(1)

            # Detemine the geolocation group name and related information
            #group = getobj(ViirsGeoFileObj,'/All_Data')
            group = ViirsGeoFileObj.getNode('/All_Data')
            geoGroupName = '/All_Data/'+group.__members__[0]
            group._g_close()
            LOG.debug("Geolocation Group : {} ".format(geoGroupName))
            isEdrGeo = ('VIIRS-CLD-AGG-GEO_All' in geoGroupName)
            if not isEdrGeo :
                LOG.debug(">> error: {} is not an EDR resolution cloud geolocation file\n\taborting...".format(ViirsGeoFileName))
                sys.exit(1)

            # Get the geolocation Latitude and Longitude nodes...
            try :
                dataName = 'Latitude'
                geoLatNode = ViirsGeoFileObj.getNode(geoGroupName+'/'+dataName)
                LOG.debug("Shape of {} node is {}".format(geoGroupName+'/'+dataName,repr(geoLatNode.shape)))
                dataName = 'Longitude'
                geoLonNode = ViirsGeoFileObj.getNode(geoGroupName+'/'+dataName)
                LOG.debug("Shape of {} node is {}".format(geoGroupName+'/'+dataName,repr(geoLonNode.shape)))
            except pyEx.NoSuchNodeError :
                LOG.debug("\n>> error: No required node {}/{} in {}\n\taborting...".format(geoGroupName,dataName,ViirsGeoFileName))
                ViirsGeoFileObj.close()
                sys.exit(1)

            # Make a copy of the geolocation node data, so we can close the file
            try :
                dataName = 'Latitude'
                LOG.debug("Reading {} dataset...".format(dataName))
                latArr = geoLatNode[:,:]
                geoLatNode.close()
                latArr = np.squeeze(latArr)
                LOG.debug("done")
                dataName = 'Longitude'
                LOG.debug("Reading {} dataset...".format(dataName))
                lonArr = geoLonNode[:,:]
                geoLonNode.close()
                lonArr = np.squeeze(lonArr)
                LOG.debug("done")
                LOG.debug("Closing geolocation file")
                ViirsGeoFileObj.close()
            except :
                LOG.debug("\n>> error: Could not retrieve %/% node data in {}\n\taborting...".format(geoGroupName,dataName,ViirsGeoFileName))
                geoLatNode.close()
                geoLonNode.close()
                ViirsGeoFileObj.close()
                sys.exit(1)

        else :
            LOG.debug("\n>> error: {} is not a HDF5 file,\n\taborting...".format(ViirsGeoFileName))
            sys.exit(1)
        
        # Read in dataSets...
        ViirsEDRFileName = cttList[grans]
        if(pytables.isHDF5File(ViirsEDRFileName)) :
            try :
                ViirsEDRFileObj = pytables.openFile(ViirsEDRFileName,mode='r')
                LOG.debug("Successfully opened edr file",ViirsEDRFileName)
            except IOError :
                LOG.debug(">> error: Could not open edr file: ",ViirsEDRFileName)
                sys.exit(1)

            # Detemine the edr group name and related information
            #group = getobj(ViirsEDRFileObj,'/All_Data')
            group = ViirsEDRFileObj.getNode('/All_Data')
            edrGroupName = '/All_Data/'+group.__members__[0]
            group._g_close()
            LOG.debug("Edr Group : {} ".format(edrGroupName))
            isEdr = ('VIIRS-CTT-EDR_All' in edrGroupName)
            if not isEdr :
                LOG.debug(">> error: {} is not an EDR resolution cloud file\n\taborting...".format(ViirsEDRFileName))
                sys.exit(1)

            # Get the edr nodes...
            try :
                dataName = 'AverageCloudTopTemperature'
                cttNode = ViirsEDRFileObj.getNode(edrGroupName+'/'+dataName)
                LOG.debug("Shape of {} node is {}".format(edrGroupName+'/'+dataName,repr(cttNode.shape)))
                dataName = 'CTTFactors'
                cttFactorsNode = ViirsEDRFileObj.getNode(edrGroupName+'/'+dataName)
                LOG.debug("Shape of {} node is {}".format(edrGroupName+'/'+dataName,repr(cttFactorsNode.shape)))
                dataName = 'QF3_VIIRSCTTAVGEDR'
                QF3_VIIRSCTTAVGEDR_Node = ViirsEDRFileObj.getNode(edrGroupName+'/'+dataName)
                LOG.debug("Shape of {} node is {}".format(edrGroupName+'/'+dataName,repr(QF3_VIIRSCTTAVGEDR_Node.shape)))
            except pyEx.NoSuchNodeError :
                LOG.debug("\n>> error: No required node {}/{} in {}\n\taborting...".format(edrGroupName,dataName,ViirsEDRFileName))
                ViirsEDRFileObj.close()
                sys.exit(1)

            # Make a copy of the edr node data, so we can close the file
            try :
                dataName = 'AverageCloudTopTemperature'
                LOG.debug("Reading {} dataset...".format(dataName))
                ctt = cttNode[:,:]
                ctt = np.squeeze(ctt)
                cttNode.close()
                LOG.debug("done")
                dataName = 'CTTFactors'
                LOG.debug("Reading {} dataset...".format(dataName))
                cttFactors = cttFactorsNode[:]
                cttFactors = np.squeeze(cttFactors)
                cttFactorsNode.close()
                LOG.debug("done")
                dataName = 'QF3_VIIRSCTTAVGEDR'
                LOG.debug("Reading {} dataset...".format(dataName))
                cloudConf = np.bitwise_and(QF3_VIIRSCTTAVGEDR_Node[:,:],3) >> 0
                waterCldFrac = np.bitwise_and(QF3_VIIRSCTTAVGEDR_Node[:,:],12) >> 2
                mixedCldFrac = np.bitwise_and(QF3_VIIRSCTTAVGEDR_Node[:,:],192) >> 6
                QF3_VIIRSCTTAVGEDR_Node.close()
                LOG.debug("done")
                LOG.debug("Closing edr file")
                ViirsEDRFileObj.close()

                LOG.debug("Shape of ctt is {}".format(repr(np.shape(ctt))))
                LOG.debug("Shape of latsArr is {}".format(repr(np.shape(latArr))))
                LOG.debug("Shape of lonsArr is {}".format(repr(np.shape(lonArr))))

            except :
                LOG.debug("\n>> error: Could not retrieve %/% node data in {}\n\taborting...".format(edrGroupName,dataName,ViirsEDRFileName))
                cttNode.close()
                cttFactorsNode.close()
                waterCldFracNode.close()
                ViirsEDRFileObj.close()
                sys.exit(1)

        else :
            LOG.debug("\n>> error: {} is not a HDF5 file,\n\taborting...".format(ViirsEDRFileName))
            sys.exit(1)
        
        LOG.debug("Creating some masks")
        try :
            
            lat_0 = latArr[np.shape(latArr)[0]/2,np.shape(latArr)[1]/2]
            lon_0 = lonArr[np.shape(lonArr)[0]/2,np.shape(lonArr)[1]/2]

            badGeo = False
            if not (-90. <= lat_0 <= 90.) :
                LOG.debug("\n>> error: Latitude of granule midpoint ({}) does not satisfy (-90. <= lat_0 <= 90.)\nfor file {}\n\taborting...".format(lat_0,geoList[grans]))
                badGeo = True
            if not (-180. <= lat_0 <= 180.) :
                LOG.debug("\n>> error: Longitude of granule midpoint ({}) does not satisfy (-180. <= lon_0 <= 180.)\nfor file {}\n\taborting...".format(lon_0,geoList[grans]))
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
                    LOG.debug("Dataset was neither int not float... a worry")
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
                LOG.debug(">> error: Mask Error, probably mismatched geolocation and product array sizes, aborting...")
                sys.exit(1)

            data = data * cttFactors[0] + cttFactors[1]

            cldPhase = np.ones(data.shape,dtype=np.uint8) * trimObj.sdrTypeFill['NA_FILL']['uint8']

            cloudFracMask = (cloudConfArr>0)
            waterFracMask = (waterCldFracArr>0)
            mixedFracMask = (mixedCldFracArr>0)

            waterMask = cloudFracMask * (waterFracMask + mixedFracMask)
            iceMask = cloudFracMask * np.logical_not(waterFracMask + mixedFracMask)

            LOG.debug("Number of pixels with > 25% cloudy confidence: ",cloudFracMask.sum())
            LOG.debug("Number of pixels with > 25% water cloud fraction: ",waterFracMask.sum())
            LOG.debug("Number of pixels with > 25% mixed cloud fraction: ",mixedFracMask.sum())

            LOG.debug("Number of water pixels: ",(cloudFracMask * (waterFracMask + mixedFracMask)).sum())
            LOG.debug("Number of ice pixels: ",cloudFracMask.sum() - (waterFracMask + mixedFracMask).sum())

            LOG.debug("Number of water pixels: ",waterMask.sum())
            LOG.debug("Number of ice pixels: ",iceMask.sum())

            cldPhaseShape = cldPhase.shape
            cldPhase = np.ravel(cldPhase)

            wtrIdx = np.where(np.ravel(waterMask))
            iceIdx = np.where(np.ravel(iceMask))
            cldPhase[wtrIdx] = 3
            cldPhase[iceIdx] = 5

        except Exception, err :
            LOG.debug(">> error: {}...".format(str(err)))
            sys.exit(1)

    LOG.debug("lat_0,lon_0 = {},{}".format(lat_0,lon_0))
    LOG.debug("Shape of data is {}".format(repr(np.shape(data))))
    LOG.debug("Shape of cldPhase is {}".format(repr(np.shape(cldPhase))))
    LOG.debug("Shape of lats is {}".format(repr(np.shape(lats))))
    LOG.debug("Shape of lons is {}".format(repr(np.shape(lons))))

    return lats,lons,data,cldPhase,lat_0,lon_0

def gran_CTP_EDR(geoList,ctpList,shrink=1):
    '''
    Returns the granulated EPS EDR
    '''

    try :
        reload(viirsCld)
        reload(viirs_edr_data)
        del(viirsCldObj)
        del(latArr)
        del(lonArr)
        del(ctpArr)
    except :
        pass

    LOG.debug("Creating viirsCldObj...")
    reload(viirsCld)
    viirsCldObj = viirsCld.viirsCld()
    LOG.debug("done")

    # Determine the correct fillValue
    trimObj = ViirsTrimTable()
    eps = 1.e-6

    for grans in np.arange(len(geoList)):
        LOG.debug("Ingesting EDR granule {} ...".format(grans))

        # Read in geolocation...
        ViirsGeoFileName = geoList[grans]
        if(pytables.isHDF5File(ViirsGeoFileName)) :
            try :
                ViirsGeoFileObj = pytables.openFile(ViirsGeoFileName,mode='r')
                LOG.debug("Successfully opened geolocation file",ViirsGeoFileName)
            except IOError :
                LOG.debug(">> error: Could not open geolocation file: ",ViirsGeoFileName)
                sys.exit(1)

            # Detemine the geolocation group name and related information
            #group = getobj(ViirsGeoFileObj,'/All_Data')
            group = ViirsGeoFileObj.getNode('/All_Data')
            geoGroupName = '/All_Data/'+group.__members__[0]
            group._g_close()
            LOG.debug("Geolocation Group : {} ".format(geoGroupName))
            isEdrGeo = ('VIIRS-CLD-AGG-GEO_All' in geoGroupName)
            if not isEdrGeo :
                LOG.debug(">> error: {} is not an EDR resolution cloud geolocation file\n\taborting...".format(ViirsGeoFileName))
                sys.exit(1)

            # Get the geolocation Latitude and Longitude nodes...
            try :
                dataName = 'Latitude'
                geoLatNode = ViirsGeoFileObj.getNode(geoGroupName+'/'+dataName)
                LOG.debug("Shape of {} node is {}".format(geoGroupName+'/'+dataName,repr(geoLatNode.shape)))
                dataName = 'Longitude'
                geoLonNode = ViirsGeoFileObj.getNode(geoGroupName+'/'+dataName)
                LOG.debug("Shape of {} node is {}".format(geoGroupName+'/'+dataName,repr(geoLonNode.shape)))
            except pyEx.NoSuchNodeError :
                LOG.debug("\n>> error: No required node {}/{} in {}\n\taborting...".format(geoGroupName,dataName,ViirsGeoFileName))
                ViirsGeoFileObj.close()
                sys.exit(1)

            # Make a copy of the geolocation node data, so we can close the file
            try :
                dataName = 'Latitude'
                LOG.debug("Reading {} dataset...".format(dataName))
                latArr = geoLatNode[:,:]
                geoLatNode.close()
                latArr = np.squeeze(latArr)
                LOG.debug("done")
                dataName = 'Longitude'
                LOG.debug("Reading {} dataset...".format(dataName))
                lonArr = geoLonNode[:,:]
                geoLonNode.close()
                lonArr = np.squeeze(lonArr)
                LOG.debug("done")
                LOG.debug("Closing geolocation file")
                ViirsGeoFileObj.close()
            except :
                LOG.debug("\n>> error: Could not retrieve %/% node data in {}\n\taborting...".format(geoGroupName,dataName,ViirsGeoFileName))
                geoLatNode.close()
                geoLonNode.close()
                ViirsGeoFileObj.close()
                sys.exit(1)

        else :
            LOG.debug("\n>> error: {} is not a HDF5 file,\n\taborting...".format(ViirsGeoFileName))
            sys.exit(1)
        
        # Read in dataSets...
        ViirsEDRFileName = ctpList[grans]
        if(pytables.isHDF5File(ViirsEDRFileName)) :
            try :
                ViirsEDRFileObj = pytables.openFile(ViirsEDRFileName,mode='r')
                LOG.debug("Successfully opened edr file",ViirsEDRFileName)
            except IOError :
                LOG.debug(">> error: Could not open edr file: ",ViirsEDRFileName)
                sys.exit(1)

            # Detemine the edr group name and related information
            #group = getobj(ViirsEDRFileObj,'/All_Data')
            group = ViirsEDRFileObj.getNode('/All_Data')
            edrGroupName = '/All_Data/'+group.__members__[0]
            group._g_close()
            LOG.debug("Edr Group : {} ".format(edrGroupName))
            isEdr = ('VIIRS-CTP-EDR_All' in edrGroupName)
            if not isEdr :
                LOG.debug(">> error: {} is not an EDR resolution cloud file\n\taborting...".format(ViirsEDRFileName))
                sys.exit(1)

            # Get the edr nodes...
            try :
                dataName = 'AverageCloudTopPressure'
                ctpNode = ViirsEDRFileObj.getNode(edrGroupName+'/'+dataName)
                LOG.debug("Shape of {} node is {}".format(edrGroupName+'/'+dataName,repr(ctpNode.shape)))
                dataName = 'CTPFactors'
                ctpFactorsNode = ViirsEDRFileObj.getNode(edrGroupName+'/'+dataName)
                LOG.debug("Shape of {} node is {}".format(edrGroupName+'/'+dataName,repr(ctpFactorsNode.shape)))
                dataName = 'QF3_VIIRSCTPAVGEDR'
                QF3_VIIRSCTPAVGEDR_Node = ViirsEDRFileObj.getNode(edrGroupName+'/'+dataName)
                LOG.debug("Shape of {} node is {}".format(edrGroupName+'/'+dataName,repr(QF3_VIIRSCTPAVGEDR_Node.shape)))
            except pyEx.NoSuchNodeError :
                LOG.debug("\n>> error: No required node {}/{} in {}\n\taborting...".format(edrGroupName,dataName,ViirsEDRFileName))
                ViirsEDRFileObj.close()
                sys.exit(1)

            # Make a copy of the edr node data, so we can close the file
            try :
                dataName = 'AverageCloudTopPressure'
                LOG.debug("Reading {} dataset...".format(dataName))
                ctp = ctpNode[:,:]
                ctp = np.squeeze(ctp)
                ctpNode.close()
                LOG.debug("done")
                dataName = 'CTPFactors'
                LOG.debug("Reading {} dataset...".format(dataName))
                ctpFactors = ctpFactorsNode[:]
                ctpFactors = np.squeeze(ctpFactors)
                ctpFactorsNode.close()
                LOG.debug("done")
                dataName = 'QF3_VIIRSCTPAVGEDR'
                LOG.debug("Reading {} dataset...".format(dataName))
                cloudConf = np.bitwise_and(QF3_VIIRSCTPAVGEDR_Node[:,:],3) >> 0
                waterCldFrac = np.bitwise_and(QF3_VIIRSCTPAVGEDR_Node[:,:],12) >> 2
                mixedCldFrac = np.bitwise_and(QF3_VIIRSCTPAVGEDR_Node[:,:],192) >> 6
                QF3_VIIRSCTPAVGEDR_Node.close()
                LOG.debug("done")
                LOG.debug("Closing edr file")
                ViirsEDRFileObj.close()

                LOG.debug("Shape of ctp is {}".format(repr(np.shape(ctp))))
                LOG.debug("Shape of latsArr is {}".format(repr(np.shape(latArr))))
                LOG.debug("Shape of lonsArr is {}".format(repr(np.shape(lonArr))))

            except :
                LOG.debug("\n>> error: Could not retrieve %/% node data in {}\n\taborting...".format(edrGroupName,dataName,ViirsEDRFileName))
                ctpNode.close()
                ctpFactorsNode.close()
                QF3_VIIRSCTPAVGEDR_Node.close()
                ViirsEDRFileObj.close()
                sys.exit(1)

        else :
            LOG.debug("\n>> error: {} is not a HDF5 file,\n\taborting...".format(ViirsEDRFileName))
            sys.exit(1)
        
        LOG.debug("Creating some masks")
        try :
            
            lat_0 = latArr[np.shape(latArr)[0]/2,np.shape(latArr)[1]/2]
            lon_0 = lonArr[np.shape(lonArr)[0]/2,np.shape(lonArr)[1]/2]

            badGeo = False
            if not (-90. <= lat_0 <= 90.) :
                LOG.debug("\n>> error: Latitude of granule midpoint ({}) does not satisfy (-90. <= lat_0 <= 90.)\nfor file {}\n\taborting...".format(lat_0,geoList[grans]))
                badGeo = True
            if not (-180. <= lat_0 <= 180.) :
                LOG.debug("\n>> error: Longitude of granule midpoint ({}) does not satisfy (-180. <= lon_0 <= 180.)\nfor file {}\n\taborting...".format(lon_0,geoList[grans]))
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
                    LOG.debug("Dataset was neither int not float... a worry")
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
                LOG.debug(">> error: Mask Error, probably mismatched geolocation and product array sizes, aborting...")
                sys.exit(1)

            data = data * ctpFactors[0] + ctpFactors[1]

            cldPhase = np.ones(data.shape,dtype=np.uint8) * trimObj.sdrTypeFill['NA_FILL']['uint8']

            cloudFracMask = (cloudConfArr>0)
            waterFracMask = (waterCldFracArr>0)
            mixedFracMask = (mixedCldFracArr>0)

            waterMask = cloudFracMask * (waterFracMask + mixedFracMask)
            iceMask = cloudFracMask * np.logical_not(waterFracMask + mixedFracMask)

            LOG.debug("Number of pixels with > 25% cloudy confidence: ",cloudFracMask.sum())
            LOG.debug("Number of pixels with > 25% water cloud fraction: ",waterFracMask.sum())
            LOG.debug("Number of pixels with > 25% mixed cloud fraction: ",mixedFracMask.sum())

            LOG.debug("Number of water pixels: ",(cloudFracMask * (waterFracMask + mixedFracMask)).sum())
            LOG.debug("Number of ice pixels: ",cloudFracMask.sum() - (waterFracMask + mixedFracMask).sum())

            LOG.debug("Number of water pixels: ",waterMask.sum())
            LOG.debug("Number of ice pixels: ",iceMask.sum())

            cldPhaseShape = cldPhase.shape
            cldPhase = np.ravel(cldPhase)

            wtrIdx = np.where(np.ravel(waterMask))
            iceIdx = np.where(np.ravel(iceMask))
            cldPhase[wtrIdx] = 3
            cldPhase[iceIdx] = 5

            cldPhase = np.reshape(cldPhase,cldPhaseShape)

        except Exception, err :
            LOG.debug(">> error: {}...".format(str(err)))
            sys.exit(1)

    LOG.debug("lat_0,lon_0 = {},{}".format(lat_0,lon_0))
    LOG.debug("Shape of data is {}".format(repr(np.shape(data))))
    LOG.debug("Shape of cldPhase is {}".format(repr(np.shape(cldPhase))))
    LOG.debug("Shape of lats is {}".format(repr(np.shape(lats))))
    LOG.debug("Shape of lons is {}".format(repr(np.shape(lons))))

    return lats,lons,data,cldPhase,lat_0,lon_0

def gran_CTH_EDR(geoList,cthList,shrink=1):
    '''
    Returns the granulated EPS EDR
    '''

    try :
        reload(viirsCld)
        reload(viirs_edr_data)
        del(viirsCldObj)
        del(latArr)
        del(lonArr)
        del(cthArr)
    except :
        pass

    LOG.debug("Creating viirsCldObj...")
    reload(viirsCld)
    viirsCldObj = viirsCld.viirsCld()
    LOG.debug("done")

    # Determine the correct fillValue
    trimObj = ViirsTrimTable()
    eps = 1.e-6

    for grans in np.arange(len(geoList)):
        LOG.debug("Ingesting EDR granule {} ...".format(grans))

        # Read in geolocation...
        ViirsGeoFileName = geoList[grans]
        if(pytables.isHDF5File(ViirsGeoFileName)) :
            try :
                ViirsGeoFileObj = pytables.openFile(ViirsGeoFileName,mode='r')
                LOG.debug("Successfully opened geolocation file",ViirsGeoFileName)
            except IOError :
                LOG.debug(">> error: Could not open geolocation file: ",ViirsGeoFileName)
                sys.exit(1)

            # Detemine the geolocation group name and related information
            #group = getobj(ViirsGeoFileObj,'/All_Data')
            group = ViirsGeoFileObj.getNode('/All_Data')
            geoGroupName = '/All_Data/'+group.__members__[0]
            group._g_close()
            LOG.debug("Geolocation Group : {} ".format(geoGroupName))
            isEdrGeo = ('VIIRS-CLD-AGG-GEO_All' in geoGroupName)
            if not isEdrGeo :
                LOG.debug(">> error: {} is not an EDR resolution cloud geolocation file\n\taborting...".format(ViirsGeoFileName))
                sys.exit(1)

            # Get the geolocation Latitude and Longitude nodes...
            try :
                dataName = 'Latitude'
                geoLatNode = ViirsGeoFileObj.getNode(geoGroupName+'/'+dataName)
                LOG.debug("Shape of {} node is {}".format(geoGroupName+'/'+dataName,repr(geoLatNode.shape)))
                dataName = 'Longitude'
                geoLonNode = ViirsGeoFileObj.getNode(geoGroupName+'/'+dataName)
                LOG.debug("Shape of {} node is {}".format(geoGroupName+'/'+dataName,repr(geoLonNode.shape)))
            except pyEx.NoSuchNodeError :
                LOG.debug("\n>> error: No required node {}/{} in {}\n\taborting...".format(geoGroupName,dataName,ViirsGeoFileName))
                ViirsGeoFileObj.close()
                sys.exit(1)

            # Make a copy of the geolocation node data, so we can close the file
            try :
                dataName = 'Latitude'
                LOG.debug("Reading {} dataset...".format(dataName))
                latArr = geoLatNode[:,:]
                geoLatNode.close()
                latArr = np.squeeze(latArr)
                LOG.debug("done")
                dataName = 'Longitude'
                LOG.debug("Reading {} dataset...".format(dataName))
                lonArr = geoLonNode[:,:]
                geoLonNode.close()
                lonArr = np.squeeze(lonArr)
                LOG.debug("done")
                LOG.debug("Closing geolocation file")
                ViirsGeoFileObj.close()
            except :
                LOG.debug("\n>> error: Could not retrieve %/% node data in {}\n\taborting...".format(geoGroupName,dataName,ViirsGeoFileName))
                geoLatNode.close()
                geoLonNode.close()
                ViirsGeoFileObj.close()
                sys.exit(1)

        else :
            LOG.debug("\n>> error: {} is not a HDF5 file,\n\taborting...".format(ViirsGeoFileName))
            sys.exit(1)
        
        # Read in dataSets...
        ViirsEDRFileName = cthList[grans]
        if(pytables.isHDF5File(ViirsEDRFileName)) :
            try :
                ViirsEDRFileObj = pytables.openFile(ViirsEDRFileName,mode='r')
                LOG.debug("Successfully opened edr file",ViirsEDRFileName)
            except IOError :
                LOG.debug(">> error: Could not open edr file: ",ViirsEDRFileName)
                sys.exit(1)

            # Detemine the edr group name and related information
            #group = getobj(ViirsEDRFileObj,'/All_Data')
            group = ViirsEDRFileObj.getNode('/All_Data')
            edrGroupName = '/All_Data/'+group.__members__[0]
            group._g_close()
            LOG.debug("Edr Group : {} ".format(edrGroupName))
            isEdr = ('VIIRS-CTH-EDR_All' in edrGroupName)
            if not isEdr :
                LOG.debug(">> error: {} is not an EDR resolution cloud file\n\taborting...".format(ViirsEDRFileName))
                sys.exit(1)

            # Get the edr nodes...
            try :
                dataName = 'AverageCloudTopHeight'
                cthNode = ViirsEDRFileObj.getNode(edrGroupName+'/'+dataName)
                LOG.debug("Shape of {} node is {}".format(edrGroupName+'/'+dataName,repr(cthNode.shape)))
                dataName = 'CTHFactors'
                cthFactorsNode = ViirsEDRFileObj.getNode(edrGroupName+'/'+dataName)
                LOG.debug("Shape of {} node is {}".format(edrGroupName+'/'+dataName,repr(cthFactorsNode.shape)))
                dataName = 'QF3_VIIRSCTHAVGEDR'
                QF3_VIIRSCTHAVGEDR_Node = ViirsEDRFileObj.getNode(edrGroupName+'/'+dataName)
                LOG.debug("Shape of {} node is {}".format(edrGroupName+'/'+dataName,repr(QF3_VIIRSCTHAVGEDR_Node.shape)))
            except pyEx.NoSuchNodeError :
                LOG.debug("\n>> error: No required node {}/{} in {}\n\taborting...".format(edrGroupName,dataName,ViirsEDRFileName))
                ViirsEDRFileObj.close()
                sys.exit(1)

            # Make a copy of the edr node data, so we can close the file
            try :
                dataName = 'AverageCloudTopHeight'
                LOG.debug("Reading {} dataset...".format(dataName))
                cth = cthNode[:,:]
                cth = np.squeeze(cth)
                cthNode.close()
                LOG.debug("done")
                dataName = 'CTHFactors'
                LOG.debug("Reading {} dataset...".format(dataName))
                cthFactors = cthFactorsNode[:]
                cthFactors = np.squeeze(cthFactors)
                cthFactorsNode.close()
                LOG.debug("done")
                dataName = 'QF3_VIIRSCTTAVGEDR'
                LOG.debug("Reading {} dataset...".format(dataName))
                cloudConf = np.bitwise_and(QF3_VIIRSCTHAVGEDR_Node[:,:],3) >> 0
                waterCldFrac = np.bitwise_and(QF3_VIIRSCTHAVGEDR_Node[:,:],12) >> 2
                mixedCldFrac = np.bitwise_and(QF3_VIIRSCTHAVGEDR_Node[:,:],192) >> 6
                QF3_VIIRSCTHAVGEDR_Node.close()
                LOG.debug("done")
                LOG.debug("Closing edr file")
                ViirsEDRFileObj.close()

                LOG.debug("Shape of cth is {}".format(repr(np.shape(cth))))
                LOG.debug("Shape of latsArr is {}".format(repr(np.shape(latArr))))
                LOG.debug("Shape of lonsArr is {}".format(repr(np.shape(lonArr))))

            except :
                LOG.debug("\n>> error: Could not retrieve %/% node data in {}\n\taborting...".format(edrGroupName,dataName,ViirsEDRFileName))
                cthNode.close()
                cthFactorsNode.close()
                QF3_VIIRSCTHAVGEDR_Node.close()
                ViirsEDRFileObj.close()
                sys.exit(1)

        else :
            LOG.debug("\n>> error: {} is not a HDF5 file,\n\taborting...".format(ViirsEDRFileName))
            sys.exit(1)
        
        LOG.debug("Creating some masks")
        try :
            
            lat_0 = latArr[np.shape(latArr)[0]/2,np.shape(latArr)[1]/2]
            lon_0 = lonArr[np.shape(lonArr)[0]/2,np.shape(lonArr)[1]/2]

            badGeo = False
            if not (-90. <= lat_0 <= 90.) :
                LOG.debug("\n>> error: Latitude of granule midpoint ({}) does not satisfy (-90. <= lat_0 <= 90.)\nfor file {}\n\taborting...".format(lat_0,geoList[grans]))
                badGeo = True
            if not (-180. <= lat_0 <= 180.) :
                LOG.debug("\n>> error: Longitude of granule midpoint ({}) does not satisfy (-180. <= lon_0 <= 180.)\nfor file {}\n\taborting...".format(lon_0,geoList[grans]))
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
                    LOG.debug("Dataset was neither int not float... a worry")
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
                LOG.debug(">> error: Mask Error, probably mismatched geolocation and product array sizes, aborting...")
                sys.exit(1)

            data = data * cthFactors[0] + cthFactors[1]

            cldPhase = np.ones(data.shape,dtype=np.uint8) * trimObj.sdrTypeFill['NA_FILL']['uint8']

            cloudFracMask = (cloudConfArr>0)
            waterFracMask = (waterCldFracArr>0)
            mixedFracMask = (mixedCldFracArr>0)

            waterMask = cloudFracMask * (waterFracMask + mixedFracMask)
            iceMask = cloudFracMask * np.logical_not(waterFracMask + mixedFracMask)

            LOG.debug("Number of pixels with > 25% cloudy confidence: ",cloudFracMask.sum())
            LOG.debug("Number of pixels with > 25% water cloud fraction: ",waterFracMask.sum())
            LOG.debug("Number of pixels with > 25% mixed cloud fraction: ",mixedFracMask.sum())

            LOG.debug("Number of water pixels: ",(cloudFracMask * (waterFracMask + mixedFracMask)).sum())
            LOG.debug("Number of ice pixels: ",cloudFracMask.sum() - (waterFracMask + mixedFracMask).sum())

            LOG.debug("Number of water pixels: ",waterMask.sum())
            LOG.debug("Number of ice pixels: ",iceMask.sum())

            cldPhaseShape = cldPhase.shape
            cldPhase = np.ravel(cldPhase)

            wtrIdx = np.where(np.ravel(waterMask))
            iceIdx = np.where(np.ravel(iceMask))
            cldPhase[wtrIdx] = 3
            cldPhase[iceIdx] = 5

            cldPhase = np.reshape(cldPhase,cldPhaseShape)

        except Exception, err :
            LOG.debug(">> error: {}...".format(str(err)))
            sys.exit(1)

    LOG.debug("lat_0,lon_0 = {},{}".format(lat_0,lon_0))
    LOG.debug("Shape of data is {}".format(repr(np.shape(data))))
    LOG.debug("Shape of cldPhase is {}".format(repr(np.shape(cldPhase))))
    LOG.debug("Shape of lats is {}".format(repr(np.shape(lats))))
    LOG.debug("Shape of lons is {}".format(repr(np.shape(lons))))

    return lats,lons,data,cldPhase,lat_0,lon_0

def gran_AOT(geoList,aotList,shrink=1):
    '''
    Returns the granulated AOT
    '''

    try :
        reload(viirsAero)
        reload(viirs_edr_data)
        del(viirsAeroObj)
        del(latArr)
        del(lonArr)
        del(aotArr)
        del(retArr)
        del(qualArr)
        del(lsmArr)
    except :
        pass

    LOG.debug("Creating viirsAeroObj...")
    reload(viirsAero)
    viirsAeroObj = viirsAero.viirsAero()
    LOG.debug("done")

    # Determine the correct fillValue
    trimObj = ViirsTrimTable()
    eps = 1.e-6

    # Build up the swath...
    for grans in np.arange(len(geoList)):

        LOG.debug("\nIngesting granule {} ...".format(grans))
        retList = viirsAeroObj.ingest(geoList[grans],aotList[grans],'aot',shrink,'linear')

        try :
            latArr  = np.vstack((latArr,viirsAeroObj.Lat[:,:]))
            lonArr  = np.vstack((lonArr,viirsAeroObj.Lon[:,:]))
            ModeGran = viirsAeroObj.ModeGran
            LOG.debug("subsequent geo arrays...")
        except NameError :
            latArr  = viirsAeroObj.Lat[:,:]
            lonArr  = viirsAeroObj.Lon[:,:]
            ModeGran = viirsAeroObj.ModeGran
            LOG.debug("first geo arrays...")

        try :
            aotArr  = np.vstack((aotArr ,viirsAeroObj.ViirsAProdSDS[:,:]))
            retArr  = np.vstack((retArr ,viirsAeroObj.ViirsAProdRet[:,:]))
            qualArr = np.vstack((qualArr,viirsAeroObj.ViirsCMquality[:,:]))
            #lsmArr  = np.vstack((lsmArr ,viirsAeroObj.LandSeaMask[:,:]))
            LOG.debug("subsequent aot arrays...")
        except NameError :
            aotArr  = viirsAeroObj.ViirsAProdSDS[:,:]
            retArr  = viirsAeroObj.ViirsAProdRet[:,:]
            qualArr = viirsAeroObj.ViirsCMquality[:,:]
            #lsmArr  = viirsAeroObj.LandSeaMask[:,:]
            LOG.debug("first aot arrays...")

        LOG.debug("Intermediate aotArr.shape = {}".format(str(aotArr.shape)))
        LOG.debug("Intermediate retArr.shape = {}".format(str(retArr.shape)))
        LOG.debug("Intermediate qualArr.shape = {}".format(str(qualArr.shape)))
        #LOG.debug("Intermediate lsmArr.shape = {}".format(str(lsmArr.shape)))

    lat_0 = latArr[np.shape(latArr)[0]/2,np.shape(latArr)[1]/2]
    lon_0 = lonArr[np.shape(lonArr)[0]/2,np.shape(lonArr)[1]/2]

    LOG.debug("lat_0,lon_0 = ",lat_0,lon_0)

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
                LOG.debug("Dataset was neither int not float... a worry")
                pass

        # Construct the total mask from all of the various fill values
        fillMask = ma.array(np.zeros(aotArr.shape,dtype=np.bool))
        for fillType in trimObj.sdrTypeFill.keys() :
            if aotFillMasks[fillType] is not None :
                fillMask = fillMask * ma.array(np.zeros(aotArr.shape,dtype=np.bool),\
                    mask=aotFillMasks[fillType])

        # Define any masks based on the quality flags...
        ViirsCMqualityMask = ma.masked_equal(qualArr,0)     # VCM quality == poor
        ViirsAProdRetMask  = ma.masked_not_equal(retArr,0)  # Interp/NAAPS/Climo

        # Define the land and water masks
        #ViirsLandMask      = ma.masked_greater(lsmArr,1)
        #ViirsWaterMask     = ma.masked_less(lsmArr,2)

        # Define the total mask
        totalMask = fillMask * ViirsCMqualityMask * ViirsAProdRetMask

        try :
            data = ma.array(aotArr,mask=totalMask.mask)
            lats = ma.array(latArr,mask=totalMask.mask)
            lons = ma.array(lonArr,mask=totalMask.mask)
        except ma.core.MaskError :
            LOG.debug(">> error: Mask Error, probably mismatched geolocation and product array sizes, aborting...")
            sys.exit(1)

    except Exception, err :
        LOG.debug(">> error: {}...".format(str(err)))
        sys.exit(1)

    LOG.debug("gran_AOT ModeGran = ",ModeGran)

    return lats,lons,data,lat_0,lon_0,ModeGran


def gran_AOT_EDR(geoList,aotList,shrink=1):
    '''
    Returns the granulated AOT EDR
    '''

    try :
        reload(viirsAero)
        reload(viirs_edr_data)
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

    LOG.debug("Creating viirsAeroObj...")
    reload(viirsAero)
    viirsAeroObj = viirsAero.viirsAero()
    LOG.debug("done")

    # Determine the correct fillValue
    trimObj = ViirsTrimTable()
    eps = 1.e-6

    # Build up the swath...
    for grans in np.arange(len(geoList)):

        LOG.debug("Ingesting EDR granule {} ...".format(grans))

        # Read in geolocation...
        ViirsGeoFileName = geoList[grans]
        LOG.debug("Geo file name: {}".format(ViirsGeoFileName))

        try :
            ViirsGeoFileObj = pytables.openFile(ViirsGeoFileName,mode='r')
            LOG.debug("Successfully opened geolocation file",ViirsGeoFileName)
        except IOError :
            LOG.debug(">> error: Could not open geolocation file: ",ViirsGeoFileName)
            sys.exit(1)

        # Detemine the geolocation group name and related information
        #group = getobj(ViirsGeoFileObj,'/All_Data')
        group = ViirsGeoFileObj.getNode('/All_Data')
        geoGroupName = '/All_Data/'+group.__members__[0]
        group._g_close()
        LOG.debug("Geolocation Group : {} ".format(geoGroupName))
        isEdrGeo = ('VIIRS-Aeros-EDR-GEO_All' in geoGroupName)
        if not isEdrGeo :
            LOG.debug(">> error: {} is not an EDR resolution aerosol geolocation file\n\taborting...".format(ViirsGeoFileName))
            sys.exit(1)

        # Determine if this is a day/night/both granule
        try:
            dataNode = ViirsGeoFileObj.getNode('/Data_Products/VIIRS-Aeros-EDR-GEO/VIIRS-Aeros-EDR-GEO_Gran_0')
            dayNightFlag = np.squeeze(dataNode.attrs['N_Day_Night_Flag'])
        except Exception, err :
            LOG.debug(">> error: {}...".format(str(err)))
            dataNode.close()

        # Get the geolocation Latitude and Longitude nodes...
        try :
            dataName = 'Latitude'
            geoLatNode = ViirsGeoFileObj.getNode(geoGroupName+'/'+dataName)
            LOG.debug("Shape of {} node is {}".format(geoGroupName+'/'+dataName,repr(geoLatNode.shape)))
            dataName = 'Longitude'
            geoLonNode = ViirsGeoFileObj.getNode(geoGroupName+'/'+dataName)
            LOG.debug("Shape of {} node is {}".format(geoGroupName+'/'+dataName,repr(geoLonNode.shape)))
            dataName = 'SolarZenithAngle'
            geoSzaNode = ViirsGeoFileObj.getNode(geoGroupName+'/'+dataName)
            LOG.debug("Shape of {} node is {}".format(geoGroupName+'/'+dataName,repr(geoSzaNode.shape)))
        except pyEx.NoSuchNodeError :
            LOG.debug("\n>> error: No required node {}/{} in {}\n\taborting...".format(geoGroupName,dataName,ViirsGeoFileName))
            ViirsGeoFileObj.close()
            sys.exit(1)

        # Make a copy of the geolocation node data, so we can close the file
        try :
            dataName = 'Latitude'
            LOG.debug("Reading {} dataset...".format(dataName))
            latArr = np.vstack((latArr,geoLatNode[:,:]))
            geoLatNode.close()
            latArr = np.squeeze(latArr)
            LOG.debug("done")
            dataName = 'Longitude'
            LOG.debug("Reading {} dataset...".format(dataName))
            lonArr = np.vstack((lonArr,geoLonNode[:,:]))
            geoLonNode.close()
            lonArr = np.squeeze(lonArr)
            LOG.debug("done")
            dataName = 'SolarZenithAngle'
            LOG.debug("Reading {} dataset...".format(dataName))
            szaArr = np.vstack((szaArr,geoSzaNode[:,:]))
            geoSzaNode.close()
            szaArr = np.squeeze(szaArr)
            LOG.debug("done")
            LOG.debug("Closing geolocation file")
            ViirsGeoFileObj.close()
        except NameError :
            latArr = geoLatNode[:,:]
            lonArr = geoLonNode[:,:]
            szaArr = geoSzaNode[:,:]
            geoLatNode.close()
            geoLonNode.close()
            geoSzaNode.close()
            ViirsGeoFileObj.close()
        #except :
            #LOG.debug("\n>> error: Could not retrieve %/% node data in {}\n\taborting...".format(geoGroupName,dataName,ViirsGeoFileName))
            #geoLatNode.close()
            #geoLonNode.close()
            #geoSzaNode.close()
            #ViirsGeoFileObj.close()
            #sys.exit(1)

        
        # Try to determine if this is a day or night granule. Night will be determined
        # to be a SZA greater than 85 degrees.

        if vars().has_key('dayNightFlag'):
            if dayNightFlag=='Day':
                ModeGran = 1
            if dayNightFlag=='Night':
                ModeGran = 0
            if dayNightFlag=='Both':
                ModeGran = 2
            LOG.debug("ModeGran = {}".format(ModeGran))
        else :
            szaMask = ma.masked_less(szaArr,85.).mask
            dayFraction = float(szaMask.sum())/float(szaArr.size)
            LOG.debug("Day Fraction = {}".format(dayFraction))
            ModeGran = 2 # Default to a mixed day/night granule
            if dayFraction == 1.00 : ModeGran = 1
            if dayFraction == 0.00 : ModeGran = 0
            LOG.debug("ModeGran = {}".format(ModeGran))

        # Read in dataSets...
        ViirsEDRFileName = aotList[grans]

        try :
            ViirsEDRFileObj = pytables.openFile(ViirsEDRFileName,mode='r')
            LOG.debug("Successfully opened edr file",ViirsEDRFileName)
        except IOError :
            LOG.debug(">> error: Could not open edr file: ",ViirsEDRFileName)
            sys.exit(1)

        # Detemine the edr group name and related information
        #group = getobj(ViirsEDRFileObj,'/All_Data')
        group = ViirsEDRFileObj.getNode('/All_Data')
        edrGroupName = '/All_Data/'+group.__members__[0]
        group._g_close()
        LOG.debug("Edr Group : {} ".format(edrGroupName))
        isEdr = ('VIIRS-Aeros-EDR_All' in edrGroupName)
        if not isEdr :
            LOG.debug(">> error: {} is not an EDR resolution aerosol file\n\taborting...".format(ViirsEDRFileName))
            sys.exit(1)

        # Get the edr nodes...
        try :
            dataName = 'AerosolOpticalDepth_at_550nm'
            aot550Node = ViirsEDRFileObj.getNode(edrGroupName+'/'+dataName)
            LOG.debug("Shape of {} node is {}".format(edrGroupName+'/'+dataName,repr(aot550Node.shape)))
            dataName = 'AerosolOpticalDepthFactors'
            aotFactorsNode = ViirsEDRFileObj.getNode(edrGroupName+'/'+dataName)
            LOG.debug("Shape of {} node is {}".format(edrGroupName+'/'+dataName,repr(aotFactorsNode.shape)))
        except pyEx.NoSuchNodeError :
            LOG.debug("\n>> error: No required node {}/{} in {}\n\taborting...".format(edrGroupName,dataName,ViirsEDRFileName))
            ViirsEDRFileObj.close()
            sys.exit(1)

        # Make a copy of the edr node data, so we can close the file
        try :
            dataName = 'AerosolOpticalDepth_at_550nm'
            LOG.debug("Reading {} dataset...".format(dataName))
            aot550 = np.vstack((aot550,aot550Node[:,:]))
            aot550 = np.squeeze(aot550)
            aot550Node.close()
            LOG.debug("done")
            dataName = 'AerosolOpticalDepthFactors'
            LOG.debug("Reading {} dataset...".format(dataName))
            aotFactors = aotFactorsNode[:]
            aotFactors = np.squeeze(aotFactors)
            aotFactorsNode.close()
            LOG.debug("done")
            LOG.debug("Closing edr file")
            ViirsEDRFileObj.close()

            LOG.debug("Shape of aot550 is {}".format(repr(np.shape(aot550))))
            LOG.debug("Shape of latsArr is {}".format(repr(np.shape(latArr))))
            LOG.debug("Shape of lonsArr is {}".format(repr(np.shape(lonArr))))

        except NameError :
            aot550 = aot550Node[:,:]
            aotFactors = aotFactorsNode[:]
            aot550Node.close()
            aotFactorsNode.close()
            ViirsEDRFileObj.close()

        #except :
            #LOG.debug("\n>> error: Could not retrieve %/% node data in {}\n\taborting...".format(edrGroupName,dataName,ViirsEDRFileName))
            #aot550Node.close()
            #aotFactorsNode.close()
            #ViirsEDRFileObj.close()
            #sys.exit(1)
        
        LOG.debug("Creating some masks")
        try :
            
            lat_0 = latArr[np.shape(latArr)[0]/2,np.shape(latArr)[1]/2]
            lon_0 = lonArr[np.shape(lonArr)[0]/2,np.shape(lonArr)[1]/2]

            badGeo = False
            if not (-90. <= lat_0 <= 90.) :
                LOG.debug("\n>> error: Latitude of granule midpoint ({}) does not satisfy (-90. <= lat_0 <= 90.)\nfor file {}\n\taborting...".format(lat_0,geoList[grans]))
                badGeo = True
            if not (-180. <= lat_0 <= 180.) :
                LOG.debug("\n>> error: Longitude of granule midpoint ({}) does not satisfy (-180. <= lon_0 <= 180.)\nfor file {}\n\taborting...".format(lon_0,geoList[grans]))
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
                    LOG.debug("Dataset was neither int not float... a worry")
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
                LOG.debug(">> error: Mask Error, probably mismatched geolocation and product array sizes, aborting...")
                sys.exit(1)

            data = data * aotFactors[0] + aotFactors[1]

        except :
            LOG.debug(">> error: There was an exception...")
            sys.exit(1)

    LOG.debug("lat_0,lon_0 = {},{}".format(lat_0,lon_0))
    LOG.debug("Shape of data is {}".format(repr(np.shape(data))))
    LOG.debug("Shape of lats is {}".format(repr(np.shape(lats))))
    LOG.debug("Shape of lons is {}".format(repr(np.shape(lons))))

    return lats,lons,data,lat_0,lon_0,ModeGran


def gran_SST(geoList,sstList,shrink=1):
    '''
    Returns the granulated SST
    '''

    try :
        reload(viirsSST)
        reload(viirs_edr_data)
        del(viirsSSTObj)
        del(latArr)
        del(lonArr)
        del(sstArr)
        del(qf2Arr)
    except :
        pass

    LOG.debug("Creating viirsSSTObj...")
    reload(viirsSST)
    viirsSSTObj = viirsSST.viirsSST()
    LOG.debug("done")

    # Determine the correct fillValue
    trimObj = ViirsTrimTable()

    # Build up the swath...
    for grans in np.arange(len(geoList)):

        LOG.debug("\nIngesting granule {} ...".format(grans))
        retList = viirsSSTObj.ingest(geoList[grans],sstList[grans],'sst',shrink)

        try :
            latArr  = np.vstack((latArr,viirsSSTObj.Lat[:,:]))
            lonArr  = np.vstack((lonArr,viirsSSTObj.Lon[:,:]))
            ModeGran = viirsSSTObj.ModeGran
            LOG.debug("subsequent geo arrays...")
        except NameError :
            latArr  = viirsSSTObj.Lat[:,:]
            lonArr  = viirsSSTObj.Lon[:,:]
            ModeGran = viirsSSTObj.ModeGran
            LOG.debug("first geo arrays...")

        try :
            sstArr  = np.vstack((sstArr ,viirsSSTObj.ViirsSSTprodSDS[:,:]))
            qf1Arr = np.vstack((qf1Arr,viirsSSTObj.ViirsSST_QF1[:,:]))
            qf2Arr = np.vstack((qf2Arr,viirsSSTObj.ViirsSST_QF2[:,:]))
            qf3Arr = np.vstack((qf3Arr,viirsSSTObj.ViirsSST_QF3[:,:]))
            qf4Arr = np.vstack((qf4Arr,viirsSSTObj.ViirsSST_QF4[:,:]))
            LOG.debug("subsequent sst arrays...")
        except NameError :
            sstArr  = viirsSSTObj.ViirsSSTprodSDS[:,:]
            qf1Arr = viirsSSTObj.ViirsSST_QF1[:,:]
            qf2Arr = viirsSSTObj.ViirsSST_QF2[:,:]
            qf3Arr = viirsSSTObj.ViirsSST_QF3[:,:]
            qf4Arr = viirsSSTObj.ViirsSST_QF4[:,:]
            LOG.debug("first sst arrays...")

        LOG.debug("Intermediate sstArr.shape = {}".format(str(sstArr.shape)))
        LOG.debug("Intermediate qf1Arr.shape = {}".format(str(qf1Arr.shape)))
        LOG.debug("Intermediate qf2Arr.shape = {}".format(str(qf2Arr.shape)))
        LOG.debug("Intermediate qf3Arr.shape = {}".format(str(qf3Arr.shape)))
        LOG.debug("Intermediate qf4Arr.shape = {}".format(str(qf4Arr.shape)))

    lat_0 = latArr[np.shape(latArr)[0]/2,np.shape(latArr)[1]/2]
    lon_0 = lonArr[np.shape(lonArr)[0]/2,np.shape(lonArr)[1]/2]

    LOG.debug("lat_0,lon_0 = ",lat_0,lon_0)

    try :
        #Determine masks for each fill type, for the SST EDR
        sstFillMasks = {}
        for fillType in trimObj.sdrTypeFill.keys() :
            fillValue = trimObj.sdrTypeFill[fillType][sstArr.dtype.name]
            if 'float' in fillValue.__class__.__name__ :
                sstFillMasks[fillType] = ma.masked_inside(sstArr,fillValue-eps,fillValue+eps).mask
                if (sstFillMasks[fillType].__class__.__name__ != 'ndarray') :
                    sstFillMasks[fillType] = None
            elif 'int' in fillValue.__class__.__name__ :
                sstFillMasks[fillType] = ma.masked_equal(sstArr,fillValue).mask
                if (sstFillMasks[fillType].__class__.__name__ != 'ndarray') :
                    sstFillMasks[fillType] = None
            else :
                LOG.debug("Dataset was neither int not float... a worry")
                pass

        #Construct the total mask from all of the various fill values
        fillMask = ma.array(np.zeros(sstArr.shape,dtype=np.bool))
        for fillType in trimObj.sdrTypeFill.keys() :
            if sstFillMasks[fillType] is not None :
                fillMask = fillMask * ma.array(np.zeros(sstArr.shape,dtype=np.bool),\
                    mask=sstFillMasks[fillType])


        # Unscale the SST dataset
        sstArr =  sstArr * viirsSSTObj.sstFactors[0] + viirsSSTObj.sstFactors[1]

        # Skin SST quality mask
        SSTqualFlag = np.bitwise_and(qf1Arr,3) >> 0
        sstQualMask = ma.masked_equal(SSTqualFlag,0).mask

        # Skin SST valid range mask
        SSTvalidRangFlag = np.bitwise_and(qf3Arr,64) >> 6
        SSTvalidRangMask = ma.masked_equal(SSTvalidRangFlag,1).mask

        # Skin SST degraded mask
        SSTdegradedFlag = np.bitwise_and(qf4Arr,1) >> 0
        SSTdegradedMask = ma.masked_equal(SSTdegradedFlag,1).mask

        # Combine the fill mask and quality masks...
        totalMask = fillMask.mask + sstQualMask + SSTvalidRangMask + SSTdegradedMask

        try :
            data = ma.array(sstArr,mask=totalMask)
            lats = ma.array(latArr,mask=totalMask)
            lons = ma.array(lonArr,mask=totalMask)
        except ma.core.MaskError :
            LOG.debug(">> error: Mask Error, probably mismatched geolocation and product array sizes, aborting...")
            sys.exit(1)

    except Exception, err :
        LOG.debug(">> error: {}...".format(str(err)))
        sys.exit(1)

    LOG.debug("gran_SST ModeGran = ",ModeGran)

    return lats,lons,data,lat_0,lon_0,ModeGran


def gran_NDVI(geoList,ndviList,prodName='NDVI',shrink=1):
    '''
    Returns the granulated AOT
    '''

    try :
        reload(viirsVI)
        reload(viirs_edr_data)
        del(viirsVIObj)
        del(latArr)
        del(lonArr)
        del(ndviArr)
        #del(qf2Arr)
    except :
        pass

    LOG.debug("Creating viirsVIObj...")
    reload(viirsVI)
    viirsVIObj = viirsVI.viirsVI()
    LOG.debug("done")

    # Determine the correct fillValue
    trimObj = ViirsTrimTable()

    # Build up the swath...
    for grans in np.arange(len(geoList)):

        LOG.debug("\nIngesting granule {} ...".format(grans))
        retList = viirsVIObj.ingest(geoList[grans],ndviList[grans],prodName,shrink)

        try :
            latArr  = np.vstack((latArr,viirsVIObj.Lat[:,:]))
            lonArr  = np.vstack((lonArr,viirsVIObj.Lon[:,:]))
            ModeGran = viirsVIObj.ModeGran
            LOG.debug("subsequent geo arrays...")
        except NameError :
            latArr  = viirsVIObj.Lat[:,:]
            lonArr  = viirsVIObj.Lon[:,:]
            ModeGran = viirsVIObj.ModeGran
            LOG.debug("first geo arrays...")

        LOG.debug("Intermediate latArr = {}".format(str(latArr.shape)))
        LOG.debug("Intermediate lonArr = {}".format(str(lonArr.shape)))

        try :
            ndviArr  = np.vstack((ndviArr ,viirsVIObj.ViirsVIprodSDS[:,:]))
            #qf1Arr = np.vstack((qf1Arr,viirsVIObj.ViirsVI_QF1[:,:]))
            qf2Arr = np.vstack((qf2Arr,viirsVIObj.ViirsVI_QF2[:,:]))
            #qf3Arr = np.vstack((qf3Arr,viirsVIObj.ViirsVI_QF3[:,:]))
            LOG.debug("subsequent ndvi arrays...")
        except NameError :
            ndviArr  = viirsVIObj.ViirsVIprodSDS[:,:]
            #qf1Arr = viirsVIObj.ViirsVI_QF1[:,:]
            qf2Arr = viirsVIObj.ViirsVI_QF2[:,:]
            #qf3Arr = viirsVIObj.ViirsVI_QF3[:,:]
            LOG.debug("first ndvi arrays...")

        LOG.debug("Intermediate ndviArr.shape = {}".format(str(ndviArr.shape)))
        #LOG.debug("Intermediate qf1Arr.shape = {}".format(str(qf1Arr.shape)))
        LOG.debug("Intermediate qf2Arr.shape = {}".format(str(qf2Arr.shape)))
        #LOG.debug("Intermediate qf3Arr.shape = {}".format(str(qf3Arr.shape)))

    lat_0 = latArr[np.shape(latArr)[0]/2,np.shape(latArr)[1]/2]
    lon_0 = lonArr[np.shape(lonArr)[0]/2,np.shape(lonArr)[1]/2]

    LOG.debug("lat_0,lon_0 = ",lat_0,lon_0)

    try :
        #Determine masks for each fill type, for the NDVI EDR
        ndviFillMasks = {}
        for fillType in trimObj.sdrTypeFill.keys() :
            fillValue = trimObj.sdrTypeFill[fillType][ndviArr.dtype.name]
            if 'float' in fillValue.__class__.__name__ :
                ndviFillMasks[fillType] = ma.masked_inside(ndviArr,fillValue-eps,fillValue+eps).mask
                if (ndviFillMasks[fillType].__class__.__name__ != 'ndarray') :
                    ndviFillMasks[fillType] = None
            elif 'int' in fillValue.__class__.__name__ :
                ndviFillMasks[fillType] = ma.masked_equal(ndviArr,fillValue).mask
                if (ndviFillMasks[fillType].__class__.__name__ != 'ndarray') :
                    ndviFillMasks[fillType] = None
            else :
                LOG.debug("Dataset was neither int not float... a worry")
                pass

        #Construct the total mask from all of the various fill values
        fillMask = ma.array(np.zeros(ndviArr.shape,dtype=np.bool))
        for fillType in trimObj.sdrTypeFill.keys() :
            if ndviFillMasks[fillType] is not None :
                fillMask = fillMask * ma.array(np.zeros(ndviArr.shape,dtype=np.bool),\
                    mask=ndviFillMasks[fillType])


        # Unscale the NDVI dataset
        ndviArr =  ndviArr * viirsVIObj.viFactors[0] + viirsVIObj.viFactors[1]

        # Define some masks...
        #fillMask = ma.masked_less(ndviArr,-800.).mask
        
        VIlandWaterFlag = np.bitwise_and(qf2Arr[:,:],7) >> 0
        VIlandWaterMask = ma.masked_greater(VIlandWaterFlag,1).mask

        VIcldConfFlag = np.bitwise_and(qf2Arr[:,:],24) >> 3
        VIcldConfMask = ma.masked_not_equal(VIcldConfFlag,0).mask

        #NDVIqualFlag = np.bitwise_and(qf1Arr,3) >> 0
        #ndviQualMask = ma.masked_equal(NDVIqualFlag,0).mask

        # Combine the fill mask and quality masks...
        #totalMask = fillMask.mask
        totalMask = fillMask.mask + VIlandWaterMask + VIcldConfMask
        #totalMask = fillMask.mask + ndviQualMask
        #totalMask = np.zeros(ndviArr.shape,dtype=np.bool)

        try :
            data = ma.array(ndviArr,mask=totalMask)
            lats = ma.array(latArr,mask=totalMask)
            lons = ma.array(lonArr,mask=totalMask)
        except ma.core.MaskError :
            LOG.debug(">> error: Mask Error, probably mismatched geolocation and product array sizes, aborting...")
            traceback.print_exc(file=sys.stdout)
            sys.exit(1)

    except Exception, err :
        LOG.debug(">> error: {}...".format(str(err)))
        traceback.print_exc(file=sys.stdout)
        sys.exit(1)

    LOG.debug("gran_NDVI ModeGran = ",ModeGran)

    return lats,lons,data,lat_0,lon_0,ModeGran


def gran_ST(geoList,stList,shrink=1):
    '''
    Returns the granulated AOT
    '''

    try :
        reload(viirsST)
        reload(viirs_edr_data)
        del(viirsSTObj)
        del(latArr)
        del(lonArr)
        del(stArr)
        del(qf1Arr)
        del(qf2Arr)
    except :
        pass

    LOG.debug("Creating viirsSTObj...")
    reload(viirsST)
    viirsSTObj = viirsST.viirsST()
    LOG.debug("done")

    # Determine the correct fillValue
    trimObj = ViirsTrimTable()

    # Build up the swath...
    for grans in np.arange(len(geoList)):

        LOG.debug("\nIngesting granule {} ...".format(grans))
        retList = viirsSTObj.ingest(geoList[grans],stList[grans],'st',shrink)

        try :
            latArr  = np.vstack((latArr,viirsSTObj.Lat[:,:]))
            lonArr  = np.vstack((lonArr,viirsSTObj.Lon[:,:]))
            ModeGran = viirsSTObj.ModeGran
            LOG.debug("subsequent geo arrays...")
        except NameError :
            latArr  = viirsSTObj.Lat[:,:]
            lonArr  = viirsSTObj.Lon[:,:]
            ModeGran = viirsSTObj.ModeGran
            LOG.debug("first geo arrays...")

        try :
            stArr  = np.vstack((stArr ,viirsSTObj.ViirsSTprodSDS[:,:]))
            qf1Arr = np.vstack((qf1Arr,viirsSTObj.ViirsST_QF1[:,:]))
            qf2Arr = np.vstack((qf2Arr,viirsSTObj.ViirsST_QF2[:,:]))
            LOG.debug("subsequent st arrays...")
        except NameError :
            stArr  = viirsSTObj.ViirsSTprodSDS[:,:]
            qf1Arr = viirsSTObj.ViirsST_QF1[:,:]
            qf2Arr = viirsSTObj.ViirsST_QF2[:,:]
            LOG.debug("first st arrays...")

        LOG.debug("Intermediate stArr.shape = {}".format(str(stArr.shape)))
        LOG.debug("Intermediate qf1Arr.shape = {}".format(str(qf1Arr.shape)))
        LOG.debug("Intermediate qf2Arr.shape = {}".format(str(qf2Arr.shape)))

    lat_0 = latArr[np.shape(latArr)[0]/2,np.shape(latArr)[1]/2]
    lon_0 = lonArr[np.shape(lonArr)[0]/2,np.shape(lonArr)[1]/2]

    LOG.debug("lat_0,lon_0 = ",lat_0,lon_0)

    try :
        #Determine masks for each fill type, for the ST EDR
        stFillMasks = {}
        for fillType in trimObj.sdrTypeFill.keys() :
            fillValue = trimObj.sdrTypeFill[fillType][stArr.dtype.name]
            if 'float' in fillValue.__class__.__name__ :
                stFillMasks[fillType] = ma.masked_inside(stArr,fillValue-eps,fillValue+eps).mask
                if (stFillMasks[fillType].__class__.__name__ != 'ndarray') :
                    stFillMasks[fillType] = None
            elif 'int' in fillValue.__class__.__name__ :
                stFillMasks[fillType] = ma.masked_equal(stArr,fillValue).mask
                if (stFillMasks[fillType].__class__.__name__ != 'ndarray') :
                    stFillMasks[fillType] = None
            else :
                LOG.debug("Dataset was neither int not float... a worry")
                pass

        #Construct the total mask from all of the various fill values
        fillMask = ma.array(np.zeros(stArr.shape,dtype=np.bool))
        for fillType in trimObj.sdrTypeFill.keys() :
            if stFillMasks[fillType] is not None :
                fillMask = fillMask * ma.array(np.zeros(stArr.shape,dtype=np.bool),\
                    mask=stFillMasks[fillType])

        # Combine the fill mask and quality masks...
        totalMask = fillMask.mask

        try :
            data = ma.array(stArr,mask=totalMask)
            lats = ma.array(latArr,mask=totalMask)
            lons = ma.array(lonArr,mask=totalMask)
        except ma.core.MaskError :
            LOG.debug(">> error: Mask Error, probably mismatched geolocation and product array sizes, aborting...")
            sys.exit(1)

    except Exception, err :
        LOG.debug(">> error: {}...".format(str(err)))
        sys.exit(1)

    LOG.debug("gran_ST ModeGran = ",ModeGran)

    return lats,lons,data,lat_0,lon_0,ModeGran


    LOG.debug("lat_0,lon_0 = {},{}".format(lat_0,lon_0))
    LOG.debug("Shape of data is {}".format(repr(np.shape(data))))
    LOG.debug("Shape of lats is {}".format(repr(np.shape(lats))))
    LOG.debug("Shape of lons is {}".format(repr(np.shape(lons))))

    return lats,lons,data,lat_0,lon_0,ModeGran


def gran_LST(geoList,lstList,shrink=1):
    '''
    Returns the granulated LST
    '''

    try :
        reload(viirsLST)
        reload(viirs_edr_data)
        del(viirsLSTObj)
        del(latArr)
        del(lonArr)
        del(lstArr)
        del(qf2Arr)
    except :
        pass

    LOG.debug("Creating viirsLSTObj...")
    reload(viirsLST)
    viirsLSTObj = viirsLST.viirsLST()
    LOG.debug("done")

    # Determine the correct fillValue
    trimObj = ViirsTrimTable()

    # Build up the swath...
    for grans in np.arange(len(geoList)):

        LOG.debug("\nIngesting granule {} ...".format(grans))
        retList = viirsLSTObj.ingest(geoList[grans],lstList[grans],'lst',shrink)

        try :
            latArr  = np.vstack((latArr,viirsLSTObj.Lat[:,:]))
            lonArr  = np.vstack((lonArr,viirsLSTObj.Lon[:,:]))
            ModeGran = viirsLSTObj.ModeGran
            LOG.debug("subsequent geo arrays...")
        except NameError :
            latArr  = viirsLSTObj.Lat[:,:]
            lonArr  = viirsLSTObj.Lon[:,:]
            ModeGran = viirsLSTObj.ModeGran
            LOG.debug("first geo arrays...")

        try :
            lstArr  = np.vstack((lstArr ,viirsLSTObj.ViirsLSTprodSDS[:,:]))
            qf1Arr = np.vstack((qf1Arr,viirsLSTObj.ViirsLST_QF1[:,:]))
            qf2Arr = np.vstack((qf2Arr,viirsLSTObj.ViirsLST_QF2[:,:]))
            #qf3Arr = np.vstack((qf3Arr,viirsLSTObj.ViirsLST_QF3[:,:]))
            LOG.debug("subsequent lst arrays...")
        except NameError :
            lstArr  = viirsLSTObj.ViirsLSTprodSDS[:,:]
            qf1Arr = viirsLSTObj.ViirsLST_QF1[:,:]
            qf2Arr = viirsLSTObj.ViirsLST_QF2[:,:]
            #qf3Arr = viirsLSTObj.ViirsLST_QF3[:,:]
            LOG.debug("first lst arrays...")

        LOG.debug("Intermediate lstArr.shape = {}".format(str(lstArr.shape)))
        LOG.debug("Intermediate qf1Arr.shape = {}".format(str(qf1Arr.shape)))
        LOG.debug("Intermediate qf2Arr.shape = {}".format(str(qf2Arr.shape)))
        #LOG.debug("Intermediate qf3Arr.shape = {}".format(str(qf3Arr.shape)))

    lat_0 = latArr[np.shape(latArr)[0]/2,np.shape(latArr)[1]/2]
    lon_0 = lonArr[np.shape(lonArr)[0]/2,np.shape(lonArr)[1]/2]

    LOG.debug("lat_0,lon_0 = ",lat_0,lon_0)

    try :
        #Determine masks for each fill type, for the LST EDR
        lstFillMasks = {}
        for fillType in trimObj.sdrTypeFill.keys() :
            fillValue = trimObj.sdrTypeFill[fillType][lstArr.dtype.name]
            if 'float' in fillValue.__class__.__name__ :
                lstFillMasks[fillType] = ma.masked_inside(lstArr,fillValue-eps,fillValue+eps).mask
                if (lstFillMasks[fillType].__class__.__name__ != 'ndarray') :
                    lstFillMasks[fillType] = None
            elif 'int' in fillValue.__class__.__name__ :
                lstFillMasks[fillType] = ma.masked_equal(lstArr,fillValue).mask
                if (lstFillMasks[fillType].__class__.__name__ != 'ndarray') :
                    lstFillMasks[fillType] = None
            else :
                LOG.debug("Dataset was neither int not float... a worry")
                pass

        #Construct the total mask from all of the various fill values
        fillMask = ma.array(np.zeros(lstArr.shape,dtype=np.bool))
        for fillType in trimObj.sdrTypeFill.keys() :
            if lstFillMasks[fillType] is not None :
                fillMask = fillMask * ma.array(np.zeros(lstArr.shape,dtype=np.bool),\
                    mask=lstFillMasks[fillType])


        # Unscale the LST dataset
        lstArr =  lstArr * viirsLSTObj.lstFactors[0] + viirsLSTObj.lstFactors[1]

        # LST quality mask
        LSTqualFlag = np.bitwise_and(qf1Arr,3) >> 0
        lstQualMask = ma.masked_equal(LSTqualFlag,3).mask

        # LST valid range mask
        LSTvalidRangeFlag = np.bitwise_and(qf2Arr,2) >> 1
        LSTvalidRangeMask = ma.masked_equal(LSTvalidRangeFlag,1).mask

        # Combine the fill mask and quality masks...
        totalMask = fillMask.mask + lstQualMask + LSTvalidRangeMask

        try :
            data = ma.array(lstArr,mask=totalMask)
            lats = ma.array(latArr,mask=totalMask)
            lons = ma.array(lonArr,mask=totalMask)
        except ma.core.MaskError :
            LOG.debug(">> error: Mask Error, probably mismatched geolocation and product array sizes, aborting...")
            sys.exit(1)

    except Exception, err :
        LOG.debug(">> error: {}...".format(str(err)))
        sys.exit(1)

    LOG.debug("gran_LST ModeGran = ",ModeGran)

    return lats,lons,data,lat_0,lon_0,ModeGran


###################################################
#                 Plotting Functions              #
###################################################

def orthoPlot_VCM(gridLat,gridLon,gridData,lat_0=0.,lon_0=0.,pointSize=1.,scatterPlot=False,\
        scale=1.3,mapRes='c', prodFileName='',outFileName='VCM.png',dpi=300,titleStr='VIIRS Cloud Mask'):
    '''
    Plots the VIIRS Cloud Mask on an orthographic projection
    '''

    # The min and max values of the dataset
    vmin = np.min(viirs_edr_data.CloudMaskData.ViirsCMvalues[cmByte][cmBit])
    vmax = np.max(viirs_edr_data.CloudMaskData.ViirsCMvalues[cmByte][cmBit])

    # If we have a zero size data array, make a dummy dataset
    # to span the allowed data range, which will be plotted with 
    # vanishing pointsize
    if (np.shape(gridLon)[0]==0) :
        LOG.debug("We have no valid data, synthesising dummy data...")
        gridLat = np.array([lat_0,lat_0])
        gridLon = np.array([lon_0,lon_0])
        gridData = np.array([vmin,vmax])
        pointSize = 0.001

    # Setup plotting data

    figWidth = 5. # inches
    figHeight = 4. # inches

    cmap = ListedColormap(viirs_edr_data.CloudMaskData.ViirsCMfillColours[cmByte][cmBit])
    numCats = np.array(viirs_edr_data.CloudMaskData.ViirsCMfillColours[cmByte][cmBit]).size
    numBounds = numCats + 1
    LOG.debug("Number of Categories: {}".format(numCats))
    LOG.debug("Number of Boundaries: {}".format(numBounds))

    # The tick positions in colourbar ([0..1]) coordinates
    viirs_edr_data.CloudMaskData.ViirsCMTickPos = np.arange(float(numBounds))/float(numCats)
    viirs_edr_data.CloudMaskData.ViirsCMTickPos = viirs_edr_data.CloudMaskData.ViirsCMTickPos[0 :-1] + \
        viirs_edr_data.CloudMaskData.ViirsCMTickPos[1]/2.

    LOG.debug("viirs_edr_data.CloudMaskData.ViirsCMTickPos: {}".format(viirs_edr_data.CloudMaskData.ViirsCMTickPos))

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
    m.drawcoastlines(ax=ax,linewidth=0.3,color='black')
    m.fillcontinents(ax=ax,color='gray',lake_color='black',zorder=0)
    #m.drawparallels(np.arange(-90.,120.,10.),linewidth=0.5,color='white')
    #m.drawmeridians(np.arange(0.,420.,10.),linewidth=0.5,color='white')
    #m.drawlsmask(ax=ax,land_color='grey',ocean_color='black',lakes=True)

    #m.bluemarble()

    # Plot the granule data
    if cmByte==0 and cmBit==1 :
        vmin,vmax = -0.5,3.5 # FIXME : This is temporary
    elif cmByte==5 and cmBit==0 :
        vmin,vmax = -0.5,7.5 # FIXME : This is temporary

    if scatterPlot:
        cs = m.scatter(x,y,s=pointSize,c=gridData,axes=ax,edgecolors='none',vmin=vmin,vmax=vmax,cmap=cmap)
    else:
        cs = m.pcolor(x,y,gridData,axes=ax,edgecolors='none',vmin=vmin,vmax=vmax,cmap=cmap,antialiased=False)

    # add a colorbar axis
    cax_rect = [0.05 , 0.05, 0.9 , 0.06 ] # [left,bottom,width,height]
    cax = fig.add_axes(cax_rect,frameon=False) # setup colorbar axes

    # Plot the colourbar
    cb = fig.colorbar(cs, cax=cax, orientation='horizontal')

    # Set the colourbar tick locations and ticklabels
    #tickPos = np.array([0,1,2,3])

    # Convert the tick positions to data coordinates
    tickPos = viirs_edr_data.CloudMaskData.ViirsCMTickPos  * numCats - 0.5
    tickLabels = viirs_edr_data.CloudMaskData.ViirsCMtickNames[cmByte][cmBit]

    #tickPos = [0.25,0.5,0.75]
    #tickLabels = ['0.25','0.5','0.75']
    #tickPos = [1.,2.,3.,4.]
    #tickLabels = ['1.','2.','3.','4.']

    LOG.debug("tickPos: {}".format(tickPos))
    LOG.debug("tickLabels: {}".format(tickLabels))

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
        LOG.debug("We have no valid data, synthesising dummy data...")
        pointSize = 5.

    x,y=m_globe(np.array(gridLon),np.array(gridLat))
    swath = np.zeros(np.shape(x),dtype=int)

    m_globe.drawmapboundary(linewidth=0.1)
    m_globe.fillcontinents(ax=glax,color='gray',zorder=1)
    m_globe.drawcoastlines(ax=glax,linewidth=0.1,zorder=3)

    p_globe = m_globe.scatter(x,y,s=pointSize,c="red",axes=glax,edgecolors='none',zorder=2)

    # Globe axis title
    glax_xlabel = ppl.setp(glax,xlabel=str(titleStr))
    ppl.setp(glax_xlabel,fontsize=6)

    # Redraw the figure
    canvas.draw()

    # save image 
    LOG.debug("Writing file to ",outFileName)
    canvas.print_figure(outFileName,dpi=dpi)


def orthoPlot_AOT(gridLat,gridLon,gridData,ModeGran, \
        vmin=0.0,vmax=1.0,scale=1.3, \
        lat_0=0.,lon_0=0.,pointSize=1.,scatterPlot=False,mapRes='c',cmap=None, \
        prodFileName='',outFileName='AOT.png',dpi=300,titleStr='VIIRS AOT'):
    '''
    Plots the VIIRS Aerosol Optical Thickness on an orthographic projection
    '''

    reload(viirs_edr_data)

    # The plot range...
    LOG.debug("vmin,vmax = ",vmin,vmax)

    # If we have a zero size data array, make a dummy dataset
    # to span the allowed data range, which will be plotted with 
    # vanishing pointsize
    if (np.shape(gridLon)[0]==0) :
        LOG.debug("We have no valid data, synthesising dummy data...")
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
    if scatterPlot:
        cs = m.scatter(x,y,s=pointSize,c=gridData,axes=ax,edgecolors='none',vmin=vmin,vmax=vmax,cmap=cmap)
    else:
        cs = m.pcolor(x,y,gridData,axes=ax,edgecolors='none',vmin=vmin,vmax=vmax,cmap=cmap,antialiased=False)

    LOG.debug("orthoPlot_AOT ModeGran = ",ModeGran)
    if (ModeGran == 0) :
        LOG.debug("Printing NIGHT text")
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

    m_globe.drawmapboundary(linewidth=0.1)
    m_globe.fillcontinents(ax=glax,color='gray',zorder=1)
    m_globe.drawcoastlines(ax=glax,linewidth=0.1,zorder=3)

    p_globe = m_globe.scatter(x,y,s=pointSize,c="red",axes=glax,edgecolors='none',zorder=2)

    # Globe axis title
    glax_xlabel = ppl.setp(glax,xlabel=titleStr)
    ppl.setp(glax_xlabel,fontsize=6)

    # Redraw the figure
    canvas.draw()

    # save image 
    LOG.debug("Writing file to ",outFileName)
    canvas.print_figure(outFileName,dpi=dpi)


def orthoPlot_SST(gridLat,gridLon,gridData,ModeGran, \
        vmin=270.,vmax=320.,scale=1.3, \
        lat_0=0.,lon_0=0.,pointSize=1.,scatterPlot=False,mapRes='c',cmap=None, \
        prodFileName='',outFileName='VSSTO.png',dpi=300,titleStr='VIIRS SST EDR'):
    '''
    Plots the VIIRS Sea Surface Temperature on an orthographic projection
    '''

    reload(viirs_edr_data)
    SeaSurfaceTempProduct = viirs_edr_data.SeaSurfaceTempProdData.SeaSurfaceTempProd()
    cmap = SeaSurfaceTempProduct.cmap

    # The plot range...
    LOG.debug("vmin,vmax = ",vmin,vmax )

    # If we have a zero size data array, make a dummy dataset
    # to span the allowed data range, which will be plotted with 
    # vanishing pointsize
    if (np.shape(gridLon)[0]==0) :
        LOG.debug("We have no valid data, synthesising dummy data...")
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
    if scatterPlot:
        cs = m.scatter(x,y,s=pointSize,c=gridData,axes=ax,edgecolors='none',vmin=vmin,vmax=vmax,cmap=cmap)
    else:
        cs = m.pcolor(x,y,gridData,axes=ax,edgecolors='none',vmin=vmin,vmax=vmax,cmap=cmap,antialiased=False)

    LOG.debug("orthoPlot_SST ModeGran = ",ModeGran)
    #if (ModeGran == 0) :
        #LOG.debug("Printing NIGHT text")
        #fig.text(0.5, 0.555, 'NIGHT',fontsize=30, color='white',ha='center',va='center',alpha=0.6)

    # add a colorbar axis
    cax_rect = [0.05 , 0.05, 0.9 , 0.06 ] # [left,bottom,width,height]
    cax = fig.add_axes(cax_rect,frameon=False) # setup colorbar axes

    # Plot the colorbar.
    cb = fig.colorbar(cs, cax=cax, orientation='horizontal')
    ppl.setp(cax.get_xticklabels(),fontsize=9)

    # Colourbar title
    cax_title = ppl.setp(cax,title="Sea Surface Temperature ($\mathrm{K}$)")
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

    m_globe.drawmapboundary(linewidth=0.1)
    m_globe.fillcontinents(ax=glax,color='gray',zorder=1)
    m_globe.drawcoastlines(ax=glax,linewidth=0.1,zorder=3)

    p_globe = m_globe.scatter(x,y,s=pointSize,c="red",axes=glax,edgecolors='none',zorder=2)

    # Globe axis title
    glax_xlabel = ppl.setp(glax,xlabel=titleStr)
    ppl.setp(glax_xlabel,fontsize=6)

    # Redraw the figure
    canvas.draw()

    # save image 
    LOG.debug("Writing file to ",outFileName)
    canvas.print_figure(outFileName,dpi=dpi)


def orthoPlot_NDVI(gridLat,gridLon,gridData,ModeGran, \
        vmin=0.,vmax=1.,scale=1.3, \
        lat_0=0.,lon_0=0.,pointSize=1.,scatterPlot=False,mapRes='c',cmap=None, \
        prodFileName='',outFileName='VIVIO.png',dpi=300,titleStr='VIIRS VI EDR'):
    '''
    Plots the VIIRS Normalised Vegetation Index on an orthographic projection
    '''

    reload(viirs_edr_data)
    VegetationIndexProduct = viirs_edr_data.VegetationIndexProdData.VegetationIndexProd()
    cmap = VegetationIndexProduct.cmap

    # The plot range...
    LOG.debug("vmin,vmax = ",vmin,vmax )

    # If we have a zero size data array, make a dummy dataset
    # to span the allowed data range, which will be plotted with 
    # vanishing pointsize
    if (np.shape(gridLon)[0]==0) :
        LOG.debug("We have no valid data, synthesising dummy data...")
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
    if scatterPlot:
        cs = m.scatter(x,y,s=pointSize,c=gridData,axes=ax,edgecolors='none',vmin=vmin,vmax=vmax,cmap=cmap)
    else:
        cs = m.pcolor(x,y,gridData,axes=ax,edgecolors='none',vmin=vmin,vmax=vmax,cmap=cmap,antialiased=False)

    LOG.debug("orthoPlot_NDVI ModeGran = ",ModeGran)
    if (ModeGran == 0) :
        LOG.debug("Printing NIGHT text")
        fig.text(0.5, 0.555, 'NIGHT',fontsize=30, color='white',ha='center',va='center',alpha=0.6)

    # add a colorbar axis
    cax_rect = [0.05 , 0.05, 0.9 , 0.06 ] # [left,bottom,width,height]
    cax = fig.add_axes(cax_rect,frameon=False) # setup colorbar axes

    # Plot the colorbar.
    cb = fig.colorbar(cs, cax=cax, orientation='horizontal')
    ppl.setp(cax.get_xticklabels(),fontsize=9)

    # Colourbar title
    cax_title = ppl.setp(cax,title="Vegetation Index")
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

    m_globe.drawmapboundary(linewidth=0.1)
    m_globe.fillcontinents(ax=glax,color='gray',zorder=1)
    m_globe.drawcoastlines(ax=glax,linewidth=0.1,zorder=3)

    p_globe = m_globe.scatter(x,y,s=pointSize,c="red",axes=glax,edgecolors='none',zorder=2)

    # Globe axis title
    glax_xlabel = ppl.setp(glax,xlabel=titleStr)
    ppl.setp(glax_xlabel,fontsize=6)

    # Redraw the figure
    canvas.draw()

    # save image 
    LOG.debug("Writing file to ",outFileName)
    canvas.print_figure(outFileName,dpi=dpi)


def orthoPlot_LST(gridLat,gridLon,gridData,ModeGran, \
        vmin=270.,vmax=320.,scale=1.3, \
        lat_0=0.,lon_0=0.,pointSize=1.,scatterPlot=False,mapRes='c',cmap=None, \
        prodFileName='',outFileName='VSSTO.png',dpi=300,titleStr='VIIRS LST EDR'):
    '''
    Plots the VIIRS Land Surface Temperature on an orthographic projection
    '''

    reload(viirs_edr_data)
    LandSurfaceTempProduct = viirs_edr_data.LandSurfaceTempProdData.LandSurfaceTempProd()
    cmap = LandSurfaceTempProduct.cmap

    # The plot range...
    LOG.debug("vmin,vmax = ",vmin,vmax )

    # If we have a zero size data array, make a dummy dataset
    # to span the allowed data range, which will be plotted with 
    # vanishing pointsize
    if (np.shape(gridLon)[0]==0) :
        LOG.debug("We have no valid data, synthesising dummy data...")
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
    if scatterPlot:
        cs = m.scatter(x,y,s=pointSize,c=gridData,axes=ax,edgecolors='none',vmin=vmin,vmax=vmax,cmap=cmap)
    else:
        cs = m.pcolor(x,y,gridData,axes=ax,edgecolors='none',vmin=vmin,vmax=vmax,cmap=cmap,antialiased=False)

    LOG.debug("orthoPlot_LST ModeGran = ",ModeGran)
    #if (ModeGran == 0) :
        #LOG.debug("Printing NIGHT text")
        #fig.text(0.5, 0.555, 'NIGHT',fontsize=30, color='white',ha='center',va='center',alpha=0.6)

    # add a colorbar axis
    cax_rect = [0.05 , 0.05, 0.9 , 0.06 ] # [left,bottom,width,height]
    cax = fig.add_axes(cax_rect,frameon=False) # setup colorbar axes

    # Plot the colorbar.
    cb = fig.colorbar(cs, cax=cax, orientation='horizontal')
    ppl.setp(cax.get_xticklabels(),fontsize=9)

    # Colourbar title
    cax_title = ppl.setp(cax,title="Land Surface Temperature ($\mathrm{K}$)")
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

    m_globe.drawmapboundary(linewidth=0.1)
    m_globe.fillcontinents(ax=glax,color='gray',zorder=1)
    m_globe.drawcoastlines(ax=glax,linewidth=0.1,zorder=3)

    p_globe = m_globe.scatter(x,y,s=pointSize,c="red",axes=glax,edgecolors='none',zorder=2)

    # Globe axis title
    glax_xlabel = ppl.setp(glax,xlabel=titleStr)
    ppl.setp(glax_xlabel,fontsize=6)

    # Redraw the figure
    canvas.draw()

    # save image 
    LOG.debug("Writing file to ",outFileName)
    canvas.print_figure(outFileName,dpi=dpi)


def orthoPlot_COP(gridLat,gridLon,gridData,gridPhase,dataSet, \
        lat_0=0.,lon_0=0.,\
        abScale='log',pointSize=1.,scatterPlot=False,scale=1.3,mapRes='c',\
        prodFileName='',outFileName='out.png',dpi=300,titleStr='VIIRS COP'):
    '''
    Plots the VIIRS Cloud Optical Parameters (COP) on an orthographic projection
    '''

    # Setup plotting data
    reload(viirs_edr_data)
    CloudProduct = viirs_edr_data.CloudProdData.CloudProduct[dataSet][abScale]

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
        LOG.debug("We have no valid water data, synthesising dummy data...")
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
        LOG.debug("We have no valid ice data, synthesising dummy data...")
        gridLat_ice = np.array([lat_0,lat_0])
        gridLon_ice = np.array([lon_0,lon_0])
        gridData_ice = np.array([CloudProduct.vmin_ice,CloudProduct.vmax_ice])
        pointSize_ice = 0.001

    # Make logarithmic if appropriate
    if(CloudProduct.logScale):
        LOG.debug("Log scaling the dataset...")
        gridData_water = np.log10(gridData_water)
        gridData_ice = np.log10(gridData_ice)
        LOG.debug("Finished log scaling the dataset...")

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
    if scatterPlot:
        cs_water = m.scatter(x,y,s=pointSize_water,c=gridData_water,axes=ax,edgecolors='none',vmin=vmin_water,vmax=vmax_water,cmap=cmap_water)
    else:
        cs_water = m.pcolor(x,y,gridData_water,axes=ax,edgecolors='none',vmin=vmin_water,vmax=vmax_water,cmap=cmap_water,antialiased=False)

    # Plot the granule ice data
    x,y=m(np.array(gridLon_ice),np.array(gridLat_ice))
    if scatterPlot:
        cs_ice = m.scatter(x,y,s=pointSize_ice,c=gridData_ice,axes=ax,edgecolors='none',vmin=vmin_ice,vmax=vmax_ice,cmap=cmap_ice)
    else:
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
        LOG.debug("vmin_water,vmax_water = ",vmin_water,vmax_water)
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
        LOG.debug("We have no valid data, synthesising dummy data...")
        gridLat = np.array([lat_0,lat_0])
        gridLon = np.array([lon_0,lon_0])
        pointSize = 5.

    x,y=m_globe(np.array(gridLon),np.array(gridLat))
    swath = np.zeros(np.shape(x),dtype=int)

    m_globe.drawmapboundary(linewidth=0.1)
    m_globe.fillcontinents(ax=glax,color='gray',zorder=1)
    m_globe.drawcoastlines(ax=glax,linewidth=0.1,zorder=3)

    p_globe = m_globe.scatter(x,y,s=pointSize,c="red",axes=glax,edgecolors='none',zorder=2)

    # Globe axis title
    glax_xlabel = ppl.setp(glax,xlabel=str(titleStr))
    ppl.setp(glax_xlabel,fontsize=6)

    # Redraw the figure
    canvas.draw()

    # save image 
    LOG.debug("Writing file to ",outFileName)
    canvas.print_figure(outFileName,dpi=dpi)

def orthoPlot_CTp(gridLat,gridLon,gridData,gridPhase,dataSet,lat_0=0.,lon_0=0.,\
        pointSize=1.,scatterPlot=False,scale=1.3,mapRes='c',prodFileName='',outFileName='out.png',dpi=300,titleStr='VIIRS CTp'):
    '''
    Plots the VIIRS Cloud Top Parameters (CTp) on an orthographic projection
    '''

    # Setup plotting data
    reload(viirs_edr_data)
    CloudProduct = viirs_edr_data.CloudProdData.CloudProduct[dataSet]

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
        LOG.debug("We have no valid water data, synthesising dummy data...")
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
        LOG.debug("We have no valid ice data, synthesising dummy data...")
        gridLat_ice = np.array([lat_0,lat_0])
        gridLon_ice = np.array([lon_0,lon_0])
        gridData_ice = np.array([CloudProduct.vmin_ice,CloudProduct.vmax_ice])
        pointSize_ice = 0.001

    # Make logarithmic if appropriate
    if(CloudProduct.logScale):
        LOG.debug("Log scaling the dataset...")
        gridData_water = np.log10(gridData_water)
        gridData_ice = np.log10(gridData_ice)
        LOG.debug("Finished log scaling the dataset...")

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
    if scatterPlot:
        cs_water = m.scatter(x,y,s=pointSize_water,c=gridData_water,axes=ax,edgecolors='none',vmin=vmin_water,vmax=vmax_water,cmap=cmap_water)
    else:
        cs_water = m.pcolor(x,y,gridData_water,axes=ax,edgecolors='none',vmin=vmin_water,vmax=vmax_water,cmap=cmap_water,antialiased=False)

    x,y=m(np.array(gridLon_ice),np.array(gridLat_ice))
    if scatterPlot:
        cs_ice = m.scatter(x,y,s=pointSize_ice,c=gridData_ice,axes=ax,edgecolors='none',vmin=vmin_ice,vmax=vmax_ice,cmap=cmap_ice)
    else:
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
        LOG.debug("We have no valid data, synthesising dummy data...")
        gridLat = np.array([lat_0,lat_0])
        gridLon = np.array([lon_0,lon_0])
        pointSize = 5.

    x,y=m_globe(np.array(gridLon),np.array(gridLat))
    swath = np.zeros(np.shape(x),dtype=int)

    m_globe.drawmapboundary(linewidth=0.1)
    m_globe.fillcontinents(ax=glax,color='gray',zorder=1)
    m_globe.drawcoastlines(ax=glax,linewidth=0.1,zorder=3)

    p_globe = m_globe.scatter(x,y,s=pointSize,c="red",axes=glax,edgecolors='none',zorder=2)

    # Globe axis title
    glax_xlabel = ppl.setp(glax,xlabel=titleStr)
    ppl.setp(glax_xlabel,fontsize=6)
    
    # Redraw the figure
    canvas.draw()

    # save image 
    LOG.debug("Writing file to ",outFileName)
    canvas.print_figure(outFileName,dpi=dpi)

###################################################
#                  Main Function                  #
###################################################

def main():

    prodChoices=['VCM','VCP','AOT','AOT_EDR','SST_EDR','NDVI','EVI','LST']
    mapRes = ['c','l','i']

    description = \
    '''
    This is a brief description of %prog
    '''
    usage = "usage: %prog [mandatory args] [options]"
    version = version="%prog"
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

    #optionalGroup.add_option('-r','--svn_revision',
                      #action="store",
                      #dest="svnRevision",
                      #default=string.split(__version__," ")[2],
                      #type="string",
                      #help="The Subversion revision number/tag of this script")
    #optionalGroup.add_option('-R','--svn_repo_path',
                      #action="store",
                      #dest="svnRepoPath",
                      #default="https://svn.ssec.wisc.edu/repos/geoffc/Python/VIIRS/"+path.basename(sys.argv[0]),
                      #type="string",
                      #help="The full Subversion repository path of this script [default is %default].")
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
    optionalGroup.add_option('--lat_0',
                      action="store",
                      dest="lat_0",
                      #default='None',
                      type="float",
                      help="If given the plot is centered on the latitude lat_0. [default: %default]")
    optionalGroup.add_option('--lon_0',
                      action="store",
                      dest="lon_0",
                      #default='None',
                      type="float",
                      help="If given the plot is centered on the longitude lon_0. [default: %default]")
    optionalGroup.add_option('-S','--stride',
                      action="store",
                      dest="stride",
                      #default='1',
                      type="int",
                      help="Sample every STRIDE pixels in the VIIRS IP/SDR product. [default: %default]")
    optionalGroup.add_option('--scatter_plot',
                      action="store_true",
                      dest="doScatterPlot",
                      default=False,
                      help="Generate the plot using a scatterplot approach.")
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

    optionalGroup.add_option('-v', '--verbose',
                      dest='verbosity',
                      action="count",
                      default=2,
                      help="""each occurrence increases verbosity 1 level from 
                      ERROR: -v=WARNING -vv=INFO -vvv=DEBUG""")

    parser.add_option_group(optionalGroup)

    # Parse the arguments from the command line
    (options, args) = parser.parse_args()

    # Set up the logging
    console_logFormat = '%(asctime)s : (%(levelname)s):%(filename)s:%(funcName)s:%(lineno)d:  %(message)s'
    dateFormat = '%Y-%m-%d %H:%M:%S'
    levels = [logging.ERROR, logging.WARN, logging.INFO, logging.DEBUG]
    logging.basicConfig(level = levels[options.verbosity], 
            format = console_logFormat, 
            datefmt = dateFormat)

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
            LOG.debug("{}".format(m_err))
    if isMissingMand :
        parser.error("Incomplete mandatory arguments, aborting...")

    # Check that the input files actually exist
    if not glob(options.geoFile) :
        parser.error("Geolocation file \n\t%s\ndoes not exist, aborting..." % (options.geoFile))
    if not glob(options.ipFile) :
        parser.error("Product file \n\t%s\ndoes not exist, aborting..." % (options.ipFile))
        
    # We should have everything we need, run the program...

    #prodFileName='''%s\n%s %s''' % (path.basename(options.ipFile),options.svnRepoPath,str(options.svnRevision))
    #prodFileName='''%s''' % (path.basename(options.ipFile))
    prodFileName=''
    dataSet = string.lower((options.ipProd).lstrip())
    mapRes = str(options.mapRes)
    vmin = options.plotMin
    vmax = options.plotMax

    CloudData = viirs_edr_data.CloudProdData.CloudProd()
    cloud_cmap = CloudData.cmap_ice_water
    CloudProdData = viirs_edr_data.CloudProdData.CloudProduct

    # Some defaults plot values if the are not specified on the command line...

    doScatterPlot = options.doScatterPlot
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

    if options.ipProd == 'VCM' or options.ipProd == 'VCP' in options.ipProd:
        LOG.debug("Calling VCM ingester...")
        stride = stride_IP if options.stride==None else options.stride
        if 'VCM' in options.ipProd :
            set_vcm_dset(0,1)
        if 'VCP' in options.ipProd :
            set_vcm_dset(5,0)

        lats,lons,vcmData,gran_lat_0,gran_lon_0 = gran_VCM(geoList,prodList,shrink=stride)

        lat_0 = options.lat_0 if (options.lat_0 is not None) else gran_lat_0 
        lon_0 = options.lon_0 if (options.lon_0 is not None) else gran_lon_0 

        LOG.debug("Calling VCM plotter...")
        pointSize = pointSize_IP if options.pointSize==None else options.pointSize
        orthoPlot_VCM(lats,lons,vcmData,lat_0=lat_0,lon_0=lon_0,\
            pointSize=pointSize,scatterPlot=doScatterPlot,scale=options.scale,mapRes=mapRes,\
            prodFileName=prodFileName,outFileName=options.outputFile,dpi=options.dpi,titleStr=options.mapAnn)

    if options.ipProd == 'AOT_EDR' :
        LOG.debug("Calling AOT EDR ingester...")
        stride = stride_EDR if options.stride==None else options.stride
        vmin = 0.0 if (vmin==None) else vmin
        vmax = 1.0   if (vmax==None) else vmax

        lats,lons,aotData,gran_lat_0,gran_lon_0,ModeGran = gran_AOT_EDR(geoList,prodList,shrink=stride)
        
        lat_0 = options.lat_0 if (options.lat_0 is not None) else gran_lat_0 
        lon_0 = options.lon_0 if (options.lon_0 is not None) else gran_lon_0 

        LOG.debug("Calling AOT plotter...")
        pointSize = pointSize_EDR if options.pointSize==None else options.pointSize
        orthoPlot_AOT(lats,lons,aotData,ModeGran,lat_0=lat_0,lon_0=lon_0,vmin=vmin,vmax=vmax,\
            pointSize=pointSize,scatterPlot=doScatterPlot,scale=options.scale,mapRes=mapRes,cmap=cloud_cmap, \
            prodFileName=prodFileName,outFileName=options.outputFile,dpi=options.dpi,titleStr=options.mapAnn)

    elif options.ipProd == 'AOT' :
        LOG.debug("Calling AOT ingester...")
        stride = stride_IP if options.stride==None else options.stride
        vmin = 0.0 if (vmin==None) else vmin
        vmax = 1.0   if (vmax==None) else vmax

        lats,lons,aotData,gran_lat_0,gran_lon_0,ModeGran = gran_AOT(geoList,prodList,shrink=stride)
        
        lat_0 = options.lat_0 if (options.lat_0 is not None) else gran_lat_0 
        lon_0 = options.lon_0 if (options.lon_0 is not None) else gran_lon_0 

        LOG.debug("Calling AOT plotter...")
        pointSize = pointSize_IP if options.pointSize==None else options.pointSize
        orthoPlot_AOT(lats,lons,aotData,ModeGran,lat_0=lat_0,lon_0=lon_0,vmin=vmin,vmax=vmax, \
            pointSize=pointSize,scatterPlot=doScatterPlot,scale=options.scale,mapRes=mapRes,cmap=cloud_cmap, \
            prodFileName=prodFileName,outFileName=options.outputFile,dpi=options.dpi,titleStr=options.mapAnn)

    elif options.ipProd == 'SST_EDR' :
        LOG.debug("Calling SST_EDR ingester...")
        stride = stride_IP if options.stride==None else options.stride

        vmin = 271. if (vmin==None) else vmin
        vmax = 318. if (vmax==None) else vmax

        lats,lons,sstData,gran_lat_0,gran_lon_0,ModeGran = gran_SST(geoList,prodList,shrink=stride)
        
        lat_0 = options.lat_0 if (options.lat_0 is not None) else gran_lat_0 
        lon_0 = options.lon_0 if (options.lon_0 is not None) else gran_lon_0 

        LOG.debug("Calling SST plotter...")
        pointSize = pointSize_IP if options.pointSize==None else options.pointSize
        orthoPlot_SST(lats,lons,sstData,ModeGran,lat_0=lat_0,lon_0=lon_0,vmin=vmin,vmax=vmax, \
            pointSize=pointSize,scatterPlot=doScatterPlot,scale=options.scale,mapRes=mapRes,cmap=None, \
            prodFileName=prodFileName,outFileName=options.outputFile,dpi=options.dpi,titleStr=options.mapAnn)

    elif options.ipProd == 'NDVI'  or options.ipProd == 'EVI':
        LOG.debug("Calling NDVI ingester...")
        stride = stride_IP if options.stride==None else options.stride
        ipProd = options.ipProd

        vmin = 0.0 if (vmin==None) else vmin
        vmax = 1.0  if (vmax==None) else vmax

        lats,lons,ndviData,gran_lat_0,gran_lon_0,ModeGran = gran_NDVI(geoList,prodList,prodName=ipProd,shrink=stride)
        
        lat_0 = options.lat_0 if (options.lat_0 is not None) else gran_lat_0 
        lon_0 = options.lon_0 if (options.lon_0 is not None) else gran_lon_0 

        if ipProd == 'NDVI':
            mapAnn = 'TOA NDVI'
        if ipProd == 'EVI':
            mapAnn = 'TOC EVI'

        LOG.debug("Calling NDVI plotter...")
        pointSize = pointSize_IP if options.pointSize==None else options.pointSize
        orthoPlot_NDVI(lats,lons,ndviData,ModeGran,lat_0=lat_0,lon_0=lon_0,vmin=vmin,vmax=vmax, \
            pointSize=pointSize,scatterPlot=doScatterPlot,scale=options.scale,mapRes=mapRes,cmap=None, \
            prodFileName=prodFileName,outFileName=options.outputFile,dpi=options.dpi,titleStr=mapAnn)

    elif options.ipProd == 'ST':
        LOG.debug("Calling ST ingester...")
        stride = stride_IP if options.stride==None else options.stride

        vmin = 0 if (vmin==None) else vmin
        vmax = 17 if (vmax==None) else vmax

        lats,lons,sstData,gran_lat_0,gran_lon_0,ModeGran = gran_ST(geoList,prodList,shrink=stride)
        
        lat_0 = options.lat_0 if (options.lat_0 is not None) else gran_lat_0 
        lon_0 = options.lon_0 if (options.lon_0 is not None) else gran_lon_0 

        LOG.debug("Calling SST plotter...")
        pointSize = pointSize_IP if options.pointSize==None else options.pointSize
        orthoPlot_ST(lats,lons,sstData,ModeGran,lat_0=lat_0,lon_0=lon_0,vmin=vmin,vmax=vmax, \
            pointSize=pointSize,scatterPlot=doScatterPlot,scale=options.scale,mapRes=mapRes,cmap=None, \
            prodFileName=prodFileName,outFileName=options.outputFile,dpi=options.dpi,titleStr=options.mapAnn)

    elif options.ipProd == 'LST' :
        LOG.debug("Calling LST ingester...")
        stride = stride_IP if options.stride==None else options.stride

        vmin = 271. if (vmin==None) else vmin
        vmax = 318. if (vmax==None) else vmax

        lats,lons,lstData,gran_lat_0,gran_lon_0,ModeGran = gran_LST(geoList,prodList,shrink=stride)
        
        lat_0 = options.lat_0 if (options.lat_0 is not None) else gran_lat_0 
        lon_0 = options.lon_0 if (options.lon_0 is not None) else gran_lon_0 

        LOG.debug("Calling LST plotter...")
        pointSize = pointSize_IP if options.pointSize==None else options.pointSize
        orthoPlot_LST(lats,lons,lstData,ModeGran,lat_0=lat_0,lon_0=lon_0,vmin=vmin,vmax=vmax, \
            pointSize=pointSize,scatterPlot=doScatterPlot,scale=options.scale,mapRes=mapRes,cmap=None, \
            prodFileName=prodFileName,outFileName=options.outputFile,dpi=options.dpi,titleStr=options.mapAnn)


    #if 'COT_EDR' in options.ipProd :
        #LOG.debug("Calling COT EDR ingester...")
        #stride = stride_EDR if options.stride==None else options.stride
        #lats,lons,cotData,cotPhase,lat_0,lon_0 = gran_COT_EDR([options.geoFile],[options.ipFile],shrink=stride)
        #dset=string.split(dataSet,'_')[0]
        #LOG.debug("Calling COT plotter...",dset
        #pointSize = pointSize_EDR if options.pointSize==None else options.pointSize
        #orthoPlot_COP(lats,lons,cotData,cotPhase,dset,lat_0=lat_0,lon_0=lon_0,\
            #pointSize=pointSize,scatterPlot=doScatterPlot,scale=options.scale,mapRes=mapRes,\
            #prodFileName=prodFileName,outFileName=options.outputFile,dpi=options.dpi,titleStr=options.mapAnn)

    #elif 'EPS_EDR' in options.ipProd :
        #LOG.debug("Calling EPS EDR ingester...")
        #stride = stride_EDR if options.stride==None else options.stride
        #lats,lons,epsData,epsPhase,lat_0,lon_0 = gran_EPS_EDR([options.geoFile],[options.ipFile],shrink=stride)
        #dset=string.split(dataSet,'_')[0]
        #LOG.debug("Calling EPS plotter...",dset
        #pointSize = pointSize_EDR if options.pointSize==None else options.pointSize
        #orthoPlot_COP(lats,lons,epsData,epsPhase,dset,lat_0=lat_0,lon_0=lon_0,\
            #pointSize=pointSize,scatterPlot=doScatterPlot,scale=options.scale,mapRes=mapRes,\
            #prodFileName=prodFileName,outFileName=options.outputFile,dpi=options.dpi,titleStr=options.mapAnn)

    #elif 'COT' in options.ipProd or 'EPS' in options.ipProd :
        #LOG.debug("Calling COP ingester...")
        #stride = stride_IP if options.stride==None else options.stride
        #lats,lons,copData,copPhase,lat_0,lon_0 = gran_COP([options.geoFile],[options.ipFile],\
            #dataSet,shrink=stride)
        #LOG.debug("Calling COP plotter...")
        #pointSize = pointSize_IP if options.pointSize==None else options.pointSize
        #orthoPlot_COP(lats,lons,copData,copPhase,dataSet,lat_0=lat_0,lon_0=lon_0,\
            #abScale='log',pointSize=pointSize,scatterPlot=doScatterPlot,scale=options.scale,mapRes=mapRes,\
            #prodFileName=prodFileName,outFileName=options.outputFile,dpi=options.dpi,titleStr=options.mapAnn)

    #if 'CTT_EDR' in options.ipProd :
        #LOG.debug("Calling CTT EDR ingester...")
        #stride = stride_EDR if options.stride==None else options.stride
        #lats,lons,cttData,cttPhase,lat_0,lon_0 = gran_CTT_EDR([options.geoFile],[options.ipFile],shrink=stride)
        #dset=string.split(dataSet,'_')[0]
        #LOG.debug("Calling CTT plotter...",dset
        #pointSize = pointSize_EDR if options.pointSize==None else options.pointSize
        #orthoPlot_CTp(lats,lons,cttData,cttPhase,dset,lat_0=lat_0,lon_0=lon_0,\
            #pointSize=pointSize,scatterPlot=doScatterPlot,scale=options.scale,mapRes=mapRes,\
            #prodFileName=prodFileName,outFileName=options.outputFile,dpi=options.dpi,titleStr=options.mapAnn)

    #elif 'CTP_EDR' in options.ipProd :
        #LOG.debug("Calling CTP EDR ingester...")
        #stride = stride_EDR if options.stride==None else options.stride
        #lats,lons,ctpData,ctpPhase,lat_0,lon_0 = gran_CTP_EDR([options.geoFile],[options.ipFile],shrink=stride)
        #dset=string.split(dataSet,'_')[0]
        #LOG.debug("Calling CTP plotter...",dset
        #pointSize = pointSize_EDR if options.pointSize==None else options.pointSize
        #orthoPlot_CTp(lats,lons,ctpData,ctpPhase,dset,lat_0=lat_0,lon_0=lon_0,\
            #pointSize=pointSize,scatterPlot=doScatterPlot,scale=options.scale,mapRes=mapRes,\
            #prodFileName=prodFileName,outFileName=options.outputFile,dpi=options.dpi,titleStr=options.mapAnn)

    #elif 'CTH_EDR' in options.ipProd :
        #LOG.debug("Calling CTH EDR ingester...")
        #stride = stride_EDR if options.stride==None else options.stride
        #lats,lons,cthData,cthPhase,lat_0,lon_0 = gran_CTH_EDR([options.geoFile],[options.ipFile],shrink=stride)
        #dset=string.split(dataSet,'_')[0]
        #LOG.debug("Calling CTH plotter...",dset
        #pointSize = pointSize_EDR if options.pointSize==None else options.pointSize
        #orthoPlot_CTp(lats,lons,cthData,cthPhase,dset,lat_0=lat_0,lon_0=lon_0,\
            #pointSize=pointSize,scatterPlot=doScatterPlot,scale=options.scale,mapRes=mapRes,\
            #prodFileName=prodFileName,outFileName=options.outputFile,dpi=options.dpi,titleStr=options.mapAnn)

    #elif 'CTT' in options.ipProd or 'CTP' in options.ipProd or 'CTH' in options.ipProd :
        #LOG.debug("Calling CTp ingester...")
        #stride = stride_IP if options.stride==None else options.stride
        #lats,lons,ctpData,ctpPhase,lat_0,lon_0 = gran_CTp([options.geoFile],[options.ipFile],\
            #dataSet,shrink=stride)
        #LOG.debug("Calling CTp plotter...")
        #pointSize = pointSize_IP if options.pointSize==None else options.pointSize
        #orthoPlot_CTp(lats,lons,ctpData,ctpPhase,dataSet,lat_0=lat_0,lon_0=lon_0,\
            #pointSize=pointSize,scatterPlot=doScatterPlot,scale=options.scale,mapRes=mapRes,\
            #prodFileName=prodFileName,outFileName=options.outputFile,dpi=options.dpi,titleStr=options.mapAnn)

    LOG.debug("Exiting...")
    sys.exit(0)

if __name__ == '__main__':
    main()


