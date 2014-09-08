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

from subprocess import CalledProcessError

import numpy as np
from bisect import bisect_left,bisect_right

import adl_blob2 as adl_blob
from adl_common import sh, env
from adl_common import ADL_HOME, CSPP_RT_HOME, CSPP_RT_ANC_PATH, CSPP_RT_ANC_CACHE_DIR, COMMON_LOG_CHECK_TABLE

# every module should have a LOG object
try :
    sourcename= file_Id.split(" ")
    LOG = logging.getLogger(sourcename[1])
except :
    LOG = logging.getLogger('Utils')


def index(a, x):
    '''Locate the leftmost value exactly equal to x'''
    i = bisect_left(a, x)
    if i != len(a) and a[i] == x:
        return i
    raise ValueError


def find_lt(a, x):
    '''Find rightmost value less than x'''
    i = bisect_left(a, x)
    if i:
        return a[i-1]
    raise ValueError


def find_le(a, x):
    '''Find rightmost value less than or equal to x'''
    i = bisect_right(a, x)
    if i:
        return a[i-1]
    raise ValueError


def find_gt(a, x):
    '''Find leftmost value greater than x'''
    i = bisect_right(a, x)
    if i != len(a):
        return a[i]
    raise ValueError


def find_ge(a, x):
    '''Find leftmost item greater than or equal to x'''
    i = bisect_left(a, x)
    if i != len(a):
        return a[i]
    raise ValueError


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


def check_exe(exeName):
    ''' Check that a required executable is in the path...'''
    try:
        retVal = sh(['which',exeName])
        LOG.info("{} is in the PATH...".format(exeName))
    except CalledProcessError:
        LOG.error("Required executable {} is not in the path or is not installed, aborting.".format(exeName))
        sys.exit(1)


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


def shipOutToFile(GridIPobj):
    '''
    Generate a blob/asc file pair from the input ancillary data object.
    '''

    # Set some environment variables and paths
    ANC_SCRIPTS_PATH = path.join(CSPP_RT_HOME,'viirs')
    ADL_ASC_TEMPLATES = path.join(ANC_SCRIPTS_PATH,'asc_templates')

    # Create new GridIP ancillary blob, and copy granulated data to it

    endian = GridIPobj.ancEndian
    if endian is adl_blob.LITTLE_ENDIAN :
        endianString = "LE"
    else :
        endianString = "BE"

    xmlName = path.join(ADL_HOME,'xml/VIIRS',GridIPobj.xmlName)

    # Create a new URID to be used in making the asc filenames

    URID_dict = getURID()

    URID = URID_dict['URID']
    creationDate_nousecStr = URID_dict['creationDate_nousecStr']
    creationDateStr = URID_dict['creationDateStr']

    # Create a new directory in the input directory for the new ancillary
    # asc and blob files

    blobDir = GridIPobj.inDir

    ascFileName = path.join(blobDir,URID+'.asc')
    blobName = path.join(blobDir,URID+'.'+GridIPobj.collectionShortName)

    LOG.debug("ascFileName : %s" % (ascFileName))
    LOG.debug("blobName : %s" % (blobName))

    # Create a new ancillary blob, and copy the data to it.
    newGridIPblobObj = adl_blob.create(xmlName, blobName, endian=endian, overwrite=True)

    # TODO: This should be a loop, so we can cycle through any data and 
    #       quality flag arrays.
    blobData = getattr(newGridIPblobObj,GridIPobj.blobDatasetName)
    blobData[:,:] = GridIPobj.data[:,:]

    # Make a new GridIP asc file from the template, and substitute for the various tags

    ascTemplateFileName = path.join(ADL_ASC_TEMPLATES,"VIIRS-GridIP-VIIRS_Template.asc")

    LOG.debug("Creating new asc file\n%s\nfrom template\n%s" % (ascFileName,ascTemplateFileName))
    
    ANC_fileList = GridIPobj.sourceList
    for idx in range(len(ANC_fileList)) :
        ANC_fileList[idx] = path.basename(ANC_fileList[idx])
    ANC_fileList.sort()
    ancGroupRecipe = '    ("N_Anc_Filename" STRING EQ "%s")'
    ancFileStr = "%s" % ("\n").join([ancGroupRecipe % (str(files)) for files in ANC_fileList])

    LOG.debug("RangeDateTimeStr = %s\n" % (GridIPobj.RangeDateTimeStr))
    LOG.debug("GRingLatitudeStr = \n%s\n" % (GridIPobj.GRingLatitudeStr))
    LOG.debug("GRingLongitudeStr = \n%s\n" % (GridIPobj.GRingLongitudeStr))

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
       line = line.replace("CSPP_ANC_COLLECTION_SHORT_NAME",GridIPobj.collectionShortName)
       line = line.replace("CSPP_GRANULE_ID",GridIPobj.geoDict['N_Granule_ID'])
       line = line.replace("CSPP_CREATIONDATETIME",creationDateStr)
       line = line.replace("  CSPP_RANGE_DATE_TIME",GridIPobj.RangeDateTimeStr)
       line = line.replace("  CSPP_GRINGLATITUDE",GridIPobj.GRingLatitudeStr)
       line = line.replace("  CSPP_GRINGLONGITUDE",GridIPobj.GRingLongitudeStr)
       line = line.replace("CSPP_NORTH_BOUNDING_COORD",GridIPobj.North_Bounding_Coordinate_Str)
       line = line.replace("CSPP_SOUTH_BOUNDING_COORD",GridIPobj.South_Bounding_Coordinate_Str)
       line = line.replace("CSPP_EAST_BOUNDING_COORD",GridIPobj.East_Bounding_Coordinate_Str)
       line = line.replace("CSPP_WEST_BOUNDING_COORD",GridIPobj.West_Bounding_Coordinate_Str)
       line = line.replace("CSPP_ANC_ENDIANNESS",endianString)       
       #line = line.replace("    CSPP_ANC_SOURCE_FILES",ancFileStr)
       ascFile.write(line) 

    ascFile.close()
    ascTemplateFile.close()

    return URID


def plotArr(data,pngName,vmin=None,vmax=None):
    '''
    Plot the input array, with a colourbar.
    '''

    # Plotting stuff
    import matplotlib
    import matplotlib.cm as cm
    from matplotlib.colors import ListedColormap
    from matplotlib.figure import Figure

    matplotlib.use('Agg')
    from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

    # This must come *after* the backend is specified.
    import matplotlib.pyplot as ppl

    LOG.info("Plotting a GridIP dataset {}".format(pngName))

    plotTitle =  string.replace(pngName,".png","")
    cbTitle   =  "Value"
    #vmin,vmax =  0,1


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


