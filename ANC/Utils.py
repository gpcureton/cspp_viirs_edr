#!/usr/bin/env python
# encoding: utf-8
"""
Utils.py

Various methods that are used by other methods in the ANC module.

Created by Geoff Cureton on 2013-03-04.
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
import uuid
from datetime import datetime,timedelta

import numpy as np
from numpy import ma

import pygrib

import adl_blob
from adl_common import ADL_HOME, CSPP_RT_ANC_PATH, CSPP_RT_ANC_CACHE_DIR, COMMON_LOG_CHECK_TABLE

# Plotting stuff
import matplotlib
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap
from matplotlib.figure import Figure

matplotlib.use('Agg')
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

# This must come *after* the backend is specified.
import matplotlib.pyplot as ppl

# every module should have a LOG object
try :
    sourcename= file_Id.split(" ")
    LOG = logging.getLogger(sourcename[1])
except :
    LOG = logging.getLogger('Utils')


def getURID() :
    '''
    Create a new URID to be used in making the asc filenames
    '''
    
    URID_dict = {}

    URID_timeObj = datetime.utcnow()
    
    creationDateStr = URID_timeObj.strftime("%Y-%m-%d %H:%M:%S.%f")
    creationDate_nousecStr = URID_timeObj.strftime("%Y-%m-%d %H:%M:%S.000000")
    
    tv_sec = int(URID_timeObj.strftime("%s"))
    tv_usec = int(URID_timeObj.strftime("%f"))
    hostId_ = uuid.getnode()
    thisAddress = id(URID_timeObj)
    
    l = tv_sec + tv_usec + hostId_ + thisAddress
    
    URID = '-'.join( ('{0:08x}'.format(tv_sec)[:8],
                      '{0:05x}'.format(tv_usec)[:5],
                      '{0:08x}'.format(hostId_)[:8],
                      '{0:08x}'.format(l)[:8]) )
    
    URID_dict['creationDateStr'] = creationDateStr
    URID_dict['creationDate_nousecStr'] = creationDate_nousecStr
    URID_dict['tv_sec'] = tv_sec
    URID_dict['tv_usec'] = tv_usec
    URID_dict['hostId_'] = hostId_
    URID_dict['thisAddress'] = thisAddress
    URID_dict['URID'] = URID
    
    return URID_dict


def getAscLine(fileObj,searchString):
    ''' Parses a file and searches for a string in each line, returning 
        the line if the string is found.'''

    dataStr = ''
    try :
        while True :
            line = fileObj.readline()

            if searchString in line : 
                dataStr = "%s" % (string.replace(line,'\n',''));
                break

        fileObj.seek(0)

    except Exception, err:
        LOG.error('Exception: %r' % (err))
        fileObj.close()

    return dataStr


def getAscStructs(fileObj,searchString,linesOfContext):
    ''' Parses a file and searches for a string in each line, returning 
        the line (and a given number of lines of context) if the string 
        is found.'''

    dataList = []
    data_count = 0
    dataFound = False

    try :
        while True :
            line = fileObj.readline()

            if searchString in line : 
                dataFound = True

            if dataFound :
                dataStr = "%s" % (string.replace(line,'\n',''));
                dataList.append(dataStr)
                data_count += 1
            else :
                pass

            if (data_count == linesOfContext) :
                break

        fileObj.seek(0)

    except Exception, err:
        LOG.error('Exception: %r' % (err))
        fileObj.close()
        return -1

    dataStr=''
    dataStr = "%s" % ("\n").join(['%s' % (str(lines)) for lines in dataList])

    return dataStr


def findDatelineCrossings(latCrnList,lonCrnList):
    '''#----------------------------------------------------------------------------
    # Finds the places where the boundary points that will make up a polygon
    # cross the dateline.
    #
    # This method is heavily based on the AltNN NNfind_crossings() method
    #
    # NOTE:  This loop will find the place(s) where the boundary crosses 180
    # degrees longitude.  It will also record the index after the crossing
    # for the first two crossings.
    # 
    # NOTE:  Since the last point in the boundary is equal to the first point
    # in the boundary, there is no chance of a crossing between the last
    # and first points.
    #
    # initialize the number of crossings to zero
    # for loop over the boundary
    #    if the longitudes cross the 180 degree line, then
    #       increment the number of crossings
    #       if this is first crossing, then
    #          save the index after the crossing
    #       else if this is the second crossing
    #          save the index after the second crossing
    #       endif
    #    endif
    # end for loop
    #-------------------------------------------------------------------------'''

    status = 0
    numCrosses = 0

    # For an ascending granule, the corner points are numbered [0,1,3,2], from the southeast
    # corner moving anti-clockwise.

    LOG.debug("latCrnList = %r " % (latCrnList))
    LOG.debug("lonCrnList = %r " % (lonCrnList))

    for idx1,idx2 in zip([1,3,2],[0,1,3]):
        
        # Convert the longitudes to radians, and calculate the 
        # absolute difference
        lon1 = np.radians(lonCrnList[idx1])
        lon2 = np.radians(lonCrnList[idx2])
        lonDiff = np.fabs( lon1 - lon2 )
        
        if ( np.fabs(lonDiff) > np.pi ):

            # We have a crossing, incrememnt the number of crossings
            numCrosses += 1
            
            if(numCrosses == 1):

                # This was the first crossing
                cross1Idx_ = idx1

            elif(numCrosses == 2):

                # This was the second crossing
                cross2Idx_ = idx1

            else :

                # we should never get here
                return -1

    num180Crossings_ = numCrosses

    '''
    # now determine the minimum and maximum latitude
    maxLat_ = latCrnList[0]
    minLat_ = maxLat_

    for idx in [1,3,2]:
        if(latCrnList[idx] > maxLat_):
            # if current lat is bigger than maxLat_, make the current point the
            # maximum
            maxLat_ = latCrnList[idx]

        if(latCrnList[idx] < minLat_):
            # if current lat is smaller than minLat_, make the current point the
            # minimum
            minLat_ = latCrnList[idx]

    return num180Crossings_,minLat_,maxLat_
    '''

    return num180Crossings_


def shipOutToFile(ANCobj):
    '''
    Generate a blob/asc file pair from the input ancillary data object.
    '''

    # Set some environment variables and paths
    CSPP_RT_HOME = os.getenv('CSPP_RT_HOME')
    ANC_SCRIPTS_PATH = path.join(CSPP_RT_HOME,'viirs')
    CSPP_RT_ANC_CACHE_DIR = os.getenv('CSPP_RT_ANC_CACHE_DIR')
    ADL_ASC_TEMPLATES = path.join(ANC_SCRIPTS_PATH,'asc_templates')

    # Create new ANC ancillary blob, and copy granulated data to it

    endian = ANCobj.ancEndian
    if endian is adl_blob.LITTLE_ENDIAN :
        endianString = "LE"
    else :
        endianString = "BE"

    xmlName = path.join(ADL_HOME,'xml/VIIRS',ANCobj.xmlName)

    # Create a new URID to be used in making the asc filenames

    URID_dict = getURID()

    URID = URID_dict['URID']
    creationDate_nousecStr = URID_dict['creationDate_nousecStr']
    creationDateStr = URID_dict['creationDateStr']

    # Create a new directory in the input directory for the new ancillary
    # asc and blob files

    blobDir = ANCobj.inDir

    ascFileName = path.join(blobDir,URID+'.asc')
    blobName = path.join(blobDir,URID+'.'+ANCobj.collectionShortName)

    LOG.debug("ascFileName : %s" % (ascFileName))
    LOG.debug("blobName : %s" % (blobName))

    # Create a new ancillary blob, and copy the data to it.
    newANCblobObj = adl_blob.create(xmlName, blobName, endian=endian, overwrite=True)
    newANCblobArrObj = newANCblobObj.as_arrays()

    blobData = getattr(newANCblobArrObj,'data')
    blobData[:,:] = ANCobj.data[:,:]

    # Make a new ANC asc file from the template, and substitute for the various tags

    ascTemplateFileName = path.join(ADL_ASC_TEMPLATES,"VIIRS-ANC_Template.asc")

    LOG.debug("Creating new asc file\n%s\nfrom template\n%s" % (ascFileName,ascTemplateFileName))
    
    ANC_fileList = ANCobj.sourceList
    LOG.info("ANC_fileList = %r" % (ANC_fileList))

    for idx in range(len(ANC_fileList)) :
        ANC_fileList[idx] = path.basename(ANC_fileList[idx])
    ANC_fileList.sort()
    ancGroupRecipe = '    ("N_Anc_Filename" STRING EQ "%s")'
    ancFileStr = "%s" % ("\n").join([ancGroupRecipe % (str(files)) for files in ANC_fileList])

    LOG.debug("RangeDateTimeStr = %s\n" % (ANCobj.RangeDateTimeStr))
    LOG.debug("GRingLatitudeStr = \n%s\n" % (ANCobj.GRingLatitudeStr))
    LOG.debug("GRingLongitudeStr = \n%s\n" % (ANCobj.GRingLongitudeStr))

    try:
        ascTemplateFile = open(ascTemplateFileName,"rt") # Open template file for reading
        ascFile = open(ascFileName,"wt") # create a new text file
    except Exception, err :
        LOG.error("%s, aborting." % (err))
        sys.exit(1)

    LOG.debug("Template file %s is %r with mode %s" %(ascTemplateFileName,'not open' if ascTemplateFile.closed else 'open',ascTemplateFile.mode))

    LOG.debug("New file %s is %r with mode %s" %(ascFileName,'not open' if ascFile.closed else 'open',ascFile.mode))

    for line in ascTemplateFile.readlines():
       line = line.replace("CSPP_URID",URID)
       line = line.replace("CSPP_CREATIONDATETIME_NOUSEC",creationDate_nousecStr)
       line = line.replace("CSPP_ANC_BLOB_FULLPATH",path.basename(blobName))
       line = line.replace("CSPP_ANC_COLLECTION_SHORT_NAME",ANCobj.collectionShortName)
       line = line.replace("CSPP_GRANULE_ID",ANCobj.geoDict['N_Granule_ID'])
       line = line.replace("CSPP_CREATIONDATETIME",creationDateStr)
       line = line.replace("  CSPP_RANGE_DATE_TIME",ANCobj.RangeDateTimeStr)
       line = line.replace("  CSPP_GRINGLATITUDE",ANCobj.GRingLatitudeStr)
       line = line.replace("  CSPP_GRINGLONGITUDE",ANCobj.GRingLongitudeStr)
       line = line.replace("    CSPP_ANC_SOURCE_FILES",ancFileStr)
       line = line.replace("CSPP_ANC_ENDIANNESS",endianString)
       ascFile.write(line) 

    ascFile.close()
    ascTemplateFile.close()

    return URID


def retrieve_NCEP_grib_files(geoDicts):
    ''' Download the GRIB files which cover the dates of the geolocation files.'''

    import shlex, subprocess
    from subprocess import CalledProcessError, call

    CSPP_RT_HOME = os.getenv('CSPP_RT_HOME')
    ANC_SCRIPTS_PATH = path.join(CSPP_RT_HOME,'viirs')
    CSPP_RT_ANC_CACHE_DIR = os.getenv('CSPP_RT_ANC_CACHE_DIR')

    gribFiles = []

    for geoDict in geoDicts:
        timeObj = geoDict['ObservedStartTime']
        dateStamp = timeObj.strftime("%Y%m%d")
        seconds = repr(int(round(timeObj.second + float(timeObj.microsecond)/1000000.)))
        deciSeconds = int(round(float(timeObj.microsecond)/100000.))
        deciSeconds = repr(0 if deciSeconds > 9 else deciSeconds)
        startTimeStamp = "%s%s" % (timeObj.strftime("%H%M%S"),deciSeconds)

        timeObj = geoDict['ObservedEndTime']
        seconds = repr(int(round(timeObj.second + float(timeObj.microsecond)/1000000.)))
        deciSeconds = int(round(float(timeObj.microsecond)/100000.))
        deciSeconds = repr(0 if deciSeconds > 9 else deciSeconds)
        endTimeStamp = "%s%s" % (timeObj.strftime("%H%M%S"),deciSeconds)

        timeObj = geoDict['UnpackTime']
        unpackTimeStamp = timeObj.strftime("%Y%m%d%H%M%S%f")

        granuleName = "GMODO_npp_d%s_t%s_e%s_b00014_c%s.h5" % (dateStamp,startTimeStamp,endTimeStamp,unpackTimeStamp)

        try :
            LOG.info('Retrieving NCEP files for %s ...' % (granuleName))
            cmdStr = '%s/cspp_retrieve_gdas_gfs.csh %s' % (ANC_SCRIPTS_PATH,granuleName)
            LOG.debug('\t%s' % (cmdStr))
            args = shlex.split(cmdStr)

            procRetVal = 0
            procObj = subprocess.Popen(args,bufsize=0, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            procObj.wait()
            procRetVal = procObj.returncode

            procOutput = procObj.stdout.readlines()
            #procOutput = procObj.stdout.read()

            # FIXME : How to get this output to have linebreaks when using readlines()
            #LOG.debug(procOutput)
            
            for lines in procOutput:
                if "GDAS/GFS file" in lines :
                    lines = string.replace(lines,'GDAS/GFS file 1: ','')
                    lines = string.replace(lines,'GDAS/GFS file 2: ','')
                    lines = string.replace(lines,'\n','')
                    gribFiles.append(lines)

            # TODO : On error, jump to a cleanup routine
            if not (procRetVal == 0) :
                LOG.error('Retrieval of ancillary files failed for %s.' % (granuleName))
                #sys.exit(procRetVal)

        except Exception, err:
            LOG.warn( "%s" % (str(err)))

    # Uniqify the list of GRIB files
    gribFiles = list(set(gribFiles))
    gribFiles.sort()

    for gribFile in gribFiles :
        LOG.info('Retrieved GRIB file: %r' % (gribFile))

    return gribFiles


def create_NCEP_grid_blobs(gribFile):
    '''Converts NCEP GRIB files into NCEP blobs'''

    from NCEPtoBlob import NCEPclass
    from thermo import rh_to_mr
    rh_to_mr_vec = np.vectorize(rh_to_mr)
    from copy import deepcopy

    CSPP_RT_HOME = os.getenv('CSPP_RT_HOME')
    ANC_SCRIPTS_PATH = path.join(CSPP_RT_HOME,'viirs')
    CSPP_RT_ANC_CACHE_DIR = os.getenv('CSPP_RT_ANC_CACHE_DIR')
    csppPython = os.getenv('PY')
    
    # Get the valid time for the grib file...
    gribFileObj = pygrib.open(gribFile)
    msg = gribFileObj.select(name="Temperature")[0]
    validDate = msg.validDate
    gribFileObj.close()
    LOG.info('NCEP GRIB file %s has valid date %r' % \
            (path.basename(gribFile),validDate.strftime("%Y-%m-%d %H:%M:%S:%f")))

    gribPath = path.dirname(gribFile)
    gribFile = path.basename(gribFile)
    gribBlob = "%s_blob.le" % (gribFile)
    gribBlob = path.join(gribPath,gribBlob)
    gribFile = path.join(gribPath,gribFile)
    LOG.debug('Candidate grib blob file name is %s ...' % (gribBlob))

    if not path.exists(gribBlob):
        LOG.info('Grib blob file %s does not exist, creating...' % (path.basename(gribBlob)))
        try :
            LOG.info('Transcoding %s to %s ...' % \
                    (path.basename(gribFile),path.basename(gribBlob)))

            NCEPxml = path.join(ADL_HOME,'xml/ANC/NCEP_ANC_Int.xml')

            # Create the grib object and populate with the grib file data
            NCEPobj = NCEPclass(gribFile=gribFile)
            LOG.debug('Successfully created NCEPobj...')

            # Convert surface pressure from Pa to mb or hPa ...
            # Ref: ADL/CMN/Utilities/ING/MSD/NCEP/src/IngMsdNCEP_Converter.cpp
            # Ref: applyScalingFactor(currentBuffer, GRID_SIZE, 0.01);
            LOG.debug('Converting the surface pressure from Pa to mb...')
            NCEPobj.NCEPmessages['surfacePressure'].data /= 100.

            # Convert total column ozone from DU or kg m**-2 to Atm.cm ...
            # Ref: ADL/CMN/Utilities/ING/MSD/NCEP/src/IngMsdNCEP_Converter.cpp
            # Code: const float DOBSON_TO_ATMSCM_SCALING_FACTOR = .001;
            # Code: applyScalingFactor(currentBuffer, GRID_SIZE,DOBSON_TO_ATMSCM_SCALING_FACTOR);
            LOG.debug('Convert total column ozone from DU or kg m**-2 to Atm.cm ...')
            NCEPobj.NCEPmessages['totalColumnOzone'].data /= 1000.

            # Convert total precipitable water kg m^{-2} to cm ...
            # Ref: ADL/CMN/Utilities/ING/MSD/NCEP/src/IngMsdNCEP_Converter.cpp
            # Code: applyScalingFactor(currentBuffer, GRID_SIZE, .10);
            LOG.debug('Convert total precipitable water kg m^{-2} to cm ...')
            NCEPobj.NCEPmessages['totalPrecipitableWater'].data /= 10.

            # Convert specific humidity in kg.kg^{-1} to water vapor mixing ratio in g.kg^{-1}
            # Ref: ADL/CMN/Utilities/ING/MSD/NCEP/src/IngMsdNCEP_Converter.cpp
            # Code: void IngMsdNCEP_Converter::applyWaterVaporMixingRatio()
            # Code: destination[i] = 1000 * (destination[i]/ (1-destination[i]));
            LOG.debug('Convert specific humidity in kg.kg^{-1} to water vapor mixing ratio in g.kg^{-1}')
            moistureObj = NCEPobj.NCEPmessages['waterVaporMixingRatioLayers'].messageLevelData

            temperatureObj = NCEPobj.NCEPmessages['temperatureLayers'].messageLevelData

            # Compute the 100mb mixing ratio in g/kg
            LOG.debug('Compute the 100mb mixing ratio in g/kg')
            if moistureObj['100'].name == 'Specific humidity':
                specHumidity_100mb = moistureObj['100'].data
                mixingRatio_100mb = 1000. * specHumidity_100mb/(1. - specHumidity_100mb)
            elif  moistureObj['100'].name == 'Relative humidity':
                relativeHumidity_100mb = moistureObj['100'].data
                temperature_100mb = temperatureObj['100'].data
                mixingRatio_100mb = rh_to_mr_vec(relativeHumidity_100mb,100.,temperature_100mb)
            else :
                pass

            LOG.debug('Compute the pressure level mixing ratios in g/kg...')
            for level,levelIdx in NCEPobj.NCEP_LAYER_LEVELS.items() : 
                levelStrIdx = level[:-2]
                pressure = NCEPobj.NCEP_LAYER_VALUES[levelIdx]

                LOG.debug('pressure = %f mb'%(pressure))
                if pressure < 100. :
                    # Copy the 100mb message object to this pressure, and assign the correct
                    # mixing ratio
                    moistureObj[levelStrIdx] = deepcopy(moistureObj['100'])
                    moistureObj[levelStrIdx].level = levelStrIdx

                    # Compute the mixing ratio in g/kg
                    mixingRatio = np.maximum(mixingRatio_100mb,0.003) * ((pressure/100.)**3.)
                    mixingRatio = np.maximum(mixingRatio_100mb,0.003)
                else :
                    # Compute the mixing ratio in g/kg
                    if moistureObj[levelStrIdx].name == 'Specific humidity':
                        specHumidity = moistureObj[levelStrIdx].data 
                        mixingRatio = 1000. * specHumidity/(1. - specHumidity)
                    elif  moistureObj[levelStrIdx].name == 'Relative humidity':
                        relativeHumidity = moistureObj[levelStrIdx].data
                        temperature = temperatureObj[levelStrIdx].data
                        mixingRatio = rh_to_mr_vec(relativeHumidity,pressure,temperature)
                    else :
                        pass

                moistureObj[levelStrIdx].data = mixingRatio

            # Write the contents of the NCEPobj object to an ADL blob file
            LOG.debug('Writing the contents of the the NCEPobj object to ADL blob file %s'%(gribBlob))
            endian = adl_blob.LITTLE_ENDIAN
            procRetVal = NCEPclass.NCEPgribToBlob_interpNew(NCEPobj,NCEPxml,gribBlob,endian=endian)

            if not (procRetVal == 0) :
                LOG.error('Transcoding of ancillary files failed for %s.' % (gribFile))
                sys.exit(procRetVal)
            else :
                LOG.debug('Finished creating NCEP GRIB blob %s' % (gribBlob))
                if not os.path.exists(gribBlob) :
                    LOG.error("Blob file error %s "%gribBlob)

        except Exception, err:
            LOG.warn( "%s" % (str(err)))
    else :
        LOG.info('NCEP global GRIB blob file %s exists, skipping.' % (path.basename(gribBlob)))


    LOG.info('Returning GRIB blob file %s with valid date %r' % \
            (path.basename(gribBlob),validDate.strftime("%Y-%m-%d %H:%M:%S:%f")))

    return validDate, gribBlob


def create_NAAPS_grid_blobs(gribFile):
    '''Converts NAAPS GRIB files into NAAPS blobs'''

    from NAAPStoBlob import NAAPSclass
    from copy import deepcopy

    CSPP_RT_HOME = os.getenv('CSPP_RT_HOME')
    ANC_SCRIPTS_PATH = path.join(CSPP_RT_HOME,'viirs')
    CSPP_RT_ANC_CACHE_DIR = os.getenv('CSPP_RT_ANC_CACHE_DIR')
    csppPython = os.getenv('PY')

    # Get the valid time for the grib file...
    gribFileObj = pygrib.open(gribFile)
    msg = gribFileObj.message(1)
    validDate = msg.validDate
    gribFileObj.close()
    LOG.info('NAAPS GRIB file %s has valid date %r' % \
            (path.basename(gribFile),validDate.strftime("%Y-%m-%d %H:%M:%S:%f")))

    gribPath = path.dirname(gribFile)
    gribFile = path.basename(gribFile)
    gribBlob = "%s_blob.le" % (gribFile)
    gribBlob = path.join(gribPath,gribBlob)
    gribFile = path.join(gribPath,gribFile)
    LOG.debug('Candidate grib blob file name is %s ...' % (gribBlob))

    if not path.exists(gribBlob):
        try :
            LOG.info('Transcoding %s to %s ...' % \
                    (path.basename(gribFile),path.basename(gribBlob)))

            NAAPSxml = path.join(ADL_HOME,'xml/ANC/NAAPS_ANC_Int.xml')

            # Create the grib object and populate with the grib file data
            NAAPSobj = NAAPSclass(gribFile=gribFile)
            LOG.debug('Successfully created NAAPSobj...')

            # Write the contents of the NAAPSobj object to an ADL blob file
            endian = adl_blob.LITTLE_ENDIAN
            procRetVal = NAAPSclass.NAAPSgribToBlob_interpNew(NAAPSobj,NAAPSxml,gribBlob,endian=endian)

            if not (procRetVal == 0) :
                LOG.error('Transcoding of ancillary files failed for %s.' % (gribFile))
                sys.exit(procRetVal)
            else :
                LOG.info('Finished creating NAAPS GRIB blob %s' % (gribBlob))

        except Exception, err:
            LOG.warn( "NAAPS: %s" % (str(err)))
    else :
        LOG.info('NAAPS global GRIB blob file %s exists, skipping.' % (path.basename(gribBlob)))


    LOG.info('Returning GRIB blob file %s with valid date %r' % \
            (path.basename(gribBlob),validDate.strftime("%Y-%m-%d %H:%M:%S:%f")))

    return validDate, gribBlob


def plotArr(data,pngName):
    '''
    Plot the input array, with a colourbar.
    '''

    plotTitle =  string.replace(pngName,".png","")
    cbTitle   =  "Value"
    #vmin,vmax =  0,1
    vmin,vmax =  None,None

    # Create figure with default size, and create canvas to draw on
    scale=1.5
    fig = Figure(figsize=(scale*8,scale*3))
    canvas = FigureCanvas(fig)

    # Create main axes instance, leaving room for colorbar at bottom,
    # and also get the Bbox of the axes instance
    ax_rect = [0.05, 0.18, 0.9, 0.75  ] # [left,bottom,width,height]
    ax = fig.add_axes(ax_rect)

    # Granule axis title
    ax_title = ppl.setp(ax,title=plotTitle)
    ppl.setp(ax_title,fontsize=12)
    ppl.setp(ax_title,family="sans-serif")

    # Plot the data
    data = ma.masked_less(data,-800.)
    im = ax.imshow(data,axes=ax,interpolation='nearest',vmin=vmin,vmax=vmax)
    
    # add a colorbar axis
    cax_rect = [0.05 , 0.05, 0.9 , 0.10 ] # [left,bottom,width,height]
    cax = fig.add_axes(cax_rect,frameon=False) # setup colorbar axes

    # Plot the colorbar.
    cb = fig.colorbar(im, cax=cax, orientation='horizontal')
    ppl.setp(cax.get_xticklabels(),fontsize=9)

    # Colourbar title
    cax_title = ppl.setp(cax,title=cbTitle)
    ppl.setp(cax_title,fontsize=9)

    # Redraw the figure
    canvas.draw()

    # save image 
    canvas.print_figure(pngName,dpi=200)


