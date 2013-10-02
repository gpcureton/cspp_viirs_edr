#!/usr/bin/env python
# encoding: utf-8
"""
NAAPStoBlob.py

The purpose of this script is to ingest NAAPS grib files and convert
to NAAPS blob format required by Raytheon's Algorithm Development Library (ADL).

Minimum commandline(s)...

python NAAPStoBlob.py -g gribFile              # To audit a grib file...
python NAAPStoBlob.py -x xmlFile -g gribFile   # To transcode a grib file...

Created by Geoff Cureton <geoff.cureton@ssec.wisc.edu> on 2011-09-06.
Copyright (c) 2011-2013 University of Wisconsin Regents. All rights reserved.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
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
from copy import deepcopy
from os import path,uname
from time import time

import numpy as np
from numpy import ma as ma

from scipy import interpolate

import pygrib

import adl_blob as adl

import optparse as optparse

# every module should have a LOG object
# e.g. LOG.warning('my dog has fleas')
import logging
LOG = logging.getLogger('NAAPStoBlob')


class NAAPSclass(object):

    """ Static Attributes """

    NAAPS_LAYER_LEVELS = {
                          '10mb' : 0,    '20mb' : 1,   '30mb' : 2,   '50mb' : 3,   '70mb' : 4,  '100mb' : 5,
                         '150mb' : 6,   '200mb' : 7,  '250mb' : 8,  '300mb' : 9,  '350mb' : 10, '400mb' : 11, 
                         '450mb' : 12,  '500mb' : 13, '550mb' : 14, '600mb' : 15, '650mb' : 16, '700mb' : 17, 
                         '750mb' : 18,  '800mb' : 19, '850mb' : 20, '900mb' : 21, '925mb' : 22, '950mb' : 23, 
                         '975mb' : 24, '1000mb' : 25
                        }

    NAAPS_LAYER_VALUES = np.array([10.0, 20.0, 30.0, 50.0, 70.0, 100.0, 150.0, 200.0, 
                                  250.0, 300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 
                                  600.0, 650.0, 700.0, 750.0, 800.0, 850.0, 900.0, 
                                  925.0, 950.0, 975.0, 1000.0 ]);

    unitsDict = {'%' : '%', \
        'Dobson' : '$\mathrm{Dobson}$', \
      'Fraction' : '$\mathrm{Fraction}$', \
           'gpm' : '$\mathrm{gpm}$', \
       'J kg-1'  : '$\mathrm{J\, kg^{-1}}$', \
             'K' : '$\mathrm{K}$', \
      'kg kg-1'  : '$\mathrm{kg\, kg^{-1}}$', \
       'g kg-1'  : '$\mathrm{g\, kg^{-1}}$', \
       'kg m-2'  : '$\mathrm{kg\, m^{-2}}$', \
             'm' : '$\mathrm{m}$', \
        'm s-1'  : '$\mathrm{m\, s^{-1}}$', \
            'Pa' : '$\mathrm{Pa}$', \
       'Pa s-1'  : '$\mathrm{Pa\, s^{-1}}$', \
    'Proportion' : '$\mathrm{Proportion}$', \
          's-1'  : '$\mathrm{s^{-1}}$'}

    # The dictionary keys in each of these dicts are the names for the various datasets
    # used the ADL NAAPS blob.

    gfsMsgKeys = [
        'aotGrid'
    ]

    ancDataSetTitle = {}
    ancDataSetTitle['aotGrid'] = 'aotGrid'

    ancName = {}
    ancName['aotGrid'] = 'unknown'

    ancParameterName = {}
    ancParameterName['aotGrid'] = 192

    ancShortName = {}
    ancShortName['aotGrid'] = 'unknown'

    ancParameterUnits = {}
    ancParameterUnits['aotGrid'] = 192

    ancUnits = {}
    ancUnits['aotGrid'] = 'unknown'

    ancTypeOfLevel = {}
    ancTypeOfLevel['aotGrid'] = 'surface'

    ancLevel = {}
    ancLevel['aotGrid'] = 0


    def __init__(self,gribFile=None,isListing=False):
        if gribFile is not None :
            gribFile = path.abspath(path.expanduser(gribFile))
            print "The grib file is %s" % (gribFile)

            # Setup some class attributes
            self.gribFile = gribFile
            self.Latitude = None
            self.Longitude = None

            # A default 0.5 degree grid which should match the gfs grid
            self.blob_Dlat,self.blob_Dlon = 0.5,0.5
            self.blob_Nlats,self.blob_Nlons = 361,720
            grids = np.mgrid[-90.:90.+self.blob_Dlat:self.blob_Dlat,0.:360.:self.blob_Dlon]
            self.blobLatitude,self.blobLongitude = grids[0],grids[1]
    
            #self.blobLatitude = np.arange(self.blob_Nlats)*self.blob_Dlat - 90.
            #self.blobLongitude = np.arange(self.blob_Nlons)*self.blob_Dlon - 180.

            self.NAAPSmessages = {}

            returnObj = not isListing

            if not returnObj :
                print "Auditing %s\n" % (gribFile)

            self.__NAAPSaudit(gribFile,returnObj=returnObj,levelInhPa=None)

        else :
            pass

    class getNAAPSlevelMessage(object):
        '''
        getNAAPSlevelMessage

        Instantiations of the getNAAPSlevelMessage class serve as containers for multiple 
        getNAAPSmessage objects making up a profile, and have as attributes various properties 
        common to each pressure level. The getNAAPSmessage object for each level contains the array 
        data for that level, along with some related attributes.
        '''
        def __init__(self,fileObj,name,parameterName,shortName,units,parameterUnits,typeOfLevel,level,pressureUnits):
            self.name = name
            self.parameterName = parameterName
            self.shortName = shortName
            self.units = units
            self.parameterUnits = parameterUnits
            self.typeOfLevel = typeOfLevel
            self.level = level
            self.pressureUnits = pressureUnits

            self.dLon = fileObj.dLon

            self.messageLevelData = {}
            for message in fileObj.msgList :
                pressureLevel = str(message['level'])
                self.msgList = [message]
                self.messageLevelData[pressureLevel] = NAAPSclass.getNAAPSmessage(self,name,parameterName,shortName,units,parameterUnits,typeOfLevel,pressureLevel,pressureUnits)
                self.messageLevelData[pressureLevel].height = "%s %s" %(pressureLevel,self.pressureUnits)

    class getNAAPSmessage(object):
        '''
        getNAAPSmessage

        Instantiations of the getNAAPSmessage class contain the array data for a single GRIB message,
        along with some related attributes.
        '''
        def __init__(self,outerObj,name,parameterName,shortName,units,parameterUnits,typeOfLevel,level,pressureUnits):
            self.name = name
            self.parameterName = parameterName
            self.shortName = shortName
            self.units = units
            self.parameterUnits = parameterUnits
            self.typeOfLevel = typeOfLevel
            self.level = level
            self.pressureUnits = pressureUnits

            self.data = outerObj.msgList[0].values
            self.dataMin = outerObj.msgList[0].min
            self.dataMax = outerObj.msgList[0].max

    @staticmethod
    def NAAPSfileAudit(gribFile,levelInhPa=None):
        '''
        NAAPSfileAudit

        This static method lists a manifest of the supplied grib file's messages 
        which match the required ADL blob datasets.
        '''
        gribFile = path.abspath(path.expanduser(gribFile))
        NAAPSobj = NAAPSclass(gribFile=gribFile, isListing=True)

    def __NAAPSaudit(self,gribFile,returnObj=False,levelInhPa=None) :
        '''
        __NAAPSaudit

        This private method opens the supplied grib file, and populates a dictionary with getNAAPSlevelMessage 
        and getNAAPSmessage objects matching the required ADL blob datasets.

        If returnObj is False, just lists the mathing grib messages and exits. 
        '''

        gribFileObj = pygrib.open(gribFile)

        gfsMsgKeys = NAAPSclass.gfsMsgKeys

        gfs_Msg = {}
        plotData = {}

        # Populate the gfs_Msg dict
        for msgKey in gfsMsgKeys :

            # Retrieve some identifying information about this message type
            name = NAAPSclass.ancName[msgKey]
            shortName = NAAPSclass.ancShortName[msgKey]
            typeOfLevel = NAAPSclass.ancTypeOfLevel[msgKey]

            if (typeOfLevel == 'isobaricInhPa') :
                level = NAAPSclass.ancLevel[msgKey] if (levelInhPa==None) else levelInhPa
            else :
                level = NAAPSclass.ancLevel[msgKey]

            print "Ingesting(%s) : %s for typeOfLevel %s..." % (msgKey,repr(name),repr(typeOfLevel))
            print "... for level %s\n" % (repr(level))
            gfs_Msg[msgKey] = gribFileObj.select(name=name,typeOfLevel=typeOfLevel,level=level)

            # The one or more messages that conform to the choice of name, typeOfLevel and level
            msgList = gfs_Msg[msgKey]
            self.msgList = msgList
            numMessages = len(msgList)
            if (numMessages>0) :
                
                #msgListIdx = 0
                msgListIdx = -1

                name           = msgList[msgListIdx]['name']
                parameterName  = msgList[msgListIdx]['parameterName']
                shortName      = msgList[msgListIdx]['shortName']
                parameterUnits = msgList[msgListIdx]['parameterUnits']
                units          = msgList[msgListIdx]['units']
                typeOfLevel    = msgList[msgListIdx]['typeOfLevel']
                #level          = msgList[msgListIdx]['level']
                pressureUnits  = msgList[msgListIdx]['pressureUnits']

                #LOG.debug("message keys = %r" % (msgList[msgListIdx].keys()))

                latitudeOfFirstGridPoint = msgList[msgListIdx]['latitudeOfFirstGridPoint']
                longitudeOfFirstGridPoint = msgList[msgListIdx]['longitudeOfFirstGridPoint']
                iDirectionIncrementGiven = msgList[msgListIdx]['iDirectionIncrementGiven']
                jDirectionIncrementGiven = msgList[msgListIdx]['jDirectionIncrementGiven']
                ijDirectionIncrementGiven = msgList[msgListIdx]['ijDirectionIncrementGiven']
                latitudeOfLastGridPoint = msgList[msgListIdx]['latitudeOfLastGridPoint']
                longitudeOfLastGridPoint = msgList[msgListIdx]['longitudeOfLastGridPoint']
                iDirectionIncrement = msgList[msgListIdx]['iDirectionIncrement']
                jDirectionIncrement = msgList[msgListIdx]['jDirectionIncrement']
                iScansNegatively = msgList[msgListIdx]['iScansNegatively']
                iScansPositively = msgList[msgListIdx]['iScansPositively']
                #jScansNegatively = msgList[msgListIdx]['jScansNegatively']
                jScansPositively = msgList[msgListIdx]['jScansPositively']
                latitudeOfFirstGridPointInDegrees = msgList[msgListIdx]['latitudeOfFirstGridPointInDegrees']
                longitudeOfFirstGridPointInDegrees = msgList[msgListIdx]['longitudeOfFirstGridPointInDegrees']
                latitudeOfLastGridPointInDegrees = msgList[msgListIdx]['latitudeOfLastGridPointInDegrees']
                longitudeOfLastGridPointInDegrees = msgList[msgListIdx]['longitudeOfLastGridPointInDegrees']
                iDirectionIncrementInDegrees = msgList[msgListIdx]['iDirectionIncrementInDegrees']
                jDirectionIncrementInDegrees = msgList[msgListIdx]['jDirectionIncrementInDegrees']
                latLonValues = msgList[msgListIdx]['latLonValues']
                latitudes = msgList[msgListIdx]['latitudes']
                longitudes = msgList[msgListIdx]['longitudes']
                distinctLatitudes = msgList[msgListIdx]['distinctLatitudes']
                distinctLongitudes = msgList[msgListIdx]['distinctLongitudes']
                latlons_Latitudes,latlons_Longitudes = msgList[msgListIdx].latlons()

                #LOG.debug("latitudeOfFirstGridPoint = %r" % (latitudeOfFirstGridPoint))
                #LOG.debug("longitudeOfFirstGridPoint = %r" % (longitudeOfFirstGridPoint))
                #LOG.debug("iDirectionIncrementGiven = %r" % (iDirectionIncrementGiven))
                #LOG.debug("jDirectionIncrementGiven = %r" % (jDirectionIncrementGiven))
                #LOG.debug("ijDirectionIncrementGiven = %r" % (ijDirectionIncrementGiven))
                #LOG.debug("latitudeOfLastGridPoint = %r" % (latitudeOfLastGridPoint))
                #LOG.debug("longitudeOfLastGridPoint = %r" % (longitudeOfLastGridPoint))
                #LOG.debug("iDirectionIncrement = %r" % (iDirectionIncrement))
                #LOG.debug("jDirectionIncrement = %r" % (jDirectionIncrement))
                #LOG.debug("iScansNegatively = %r" % (iScansNegatively))
                #LOG.debug("iScansPositively = %r" % (iScansPositively))
                #LOG.debug("jScansNegatively = %r" % (jScansNegatively))
                #LOG.debug("jScansPositively = %r" % (jScansPositively))
                LOG.debug("latitudeOfFirstGridPointInDegrees = %r" % (latitudeOfFirstGridPointInDegrees))
                LOG.debug("latitudeOfLastGridPointInDegrees = %r" % ( latitudeOfLastGridPointInDegrees))
                LOG.debug("longitudeOfFirstGridPointInDegrees = %r" % (longitudeOfFirstGridPointInDegrees))
                LOG.debug("longitudeOfLastGridPointInDegrees = %r" % (longitudeOfLastGridPointInDegrees))
                LOG.debug("iDirectionIncrementInDegrees = %r" % (iDirectionIncrementInDegrees))
                LOG.debug("jDirectionIncrementInDegrees = %r" % (jDirectionIncrementInDegrees))
                LOG.debug("latLonValues = %r , %r" % (latLonValues.shape,latLonValues))
                LOG.debug("latitudes = %r , %r" % (latitudes.shape,latitudes))
                #LOG.debug("longitudes = %r , %r" % (longitudes.shape,longitudes))
                LOG.debug("distinctLatitudes = %r , %r" % (distinctLatitudes.shape,distinctLatitudes))
                LOG.debug("distinctLongitudes = %r , %r" % (distinctLongitudes.shape,distinctLongitudes))
                LOG.debug("latlons_Latitudes = %r , %r" % (latlons_Latitudes.shape,latlons_Latitudes[:,0]))
                LOG.debug("latlons_Longitudes = %r , %r" % (latlons_Longitudes.shape,latlons_Longitudes[0,:]))

                self.Nlats = msgList[msgListIdx]['Nj']
                self.Nlons = msgList[msgListIdx]['Ni']
                self.dLat =  msgList[msgListIdx]['jDirectionIncrementInDegrees']
                self.dLon =  msgList[msgListIdx]['iDirectionIncrementInDegrees']
                #self.Latitude, self.Longitude = msgList[msgListIdx].latlons()
                self.Latitude = distinctLatitudes \
                    if distinctLatitudes[0] < distinctLatitudes[-1] else distinctLatitudes[::-1]
                self.Longitude = distinctLongitudes

                LOG.debug("self.Nlats =%r" % ( self.Nlats))
                LOG.debug("self.Nlons =%r" % ( self.Nlons))
                LOG.debug("self.dLat = %r" % ( self.dLat ))
                LOG.debug("self.dLon = %r" % ( self.dLon ))
                LOG.debug("self.Latitude = %r , %r" % (self.Latitude.shape,self.Latitude))
                LOG.debug("self.Longitude = %r , %r" % (self.Longitude.shape,self.Longitude))


                LOG.debug("There are %d messages for %s with...\n\tname = \"%s\"\n\tparameterName = \"%s\"\n\tshortName = %s\n\tparameterUnits = %s\n\tunits = %s\n\ttypeOfLevel = %s\n\tpressureUnits = %s\n\tlevel = %s" \
                    % (len(msgList), msgKey, name, parameterName, shortName, 
                            parameterUnits, units, typeOfLevel, pressureUnits, level
                      ))

                # TODO : If we just want a list we should stop here...
                if not returnObj :
                    pass
                    #print "Auditing %s\n" % (gribFile)
                else :

                    # Determine the level index, name and height of the dataset
                    if (msgList[msgListIdx]['typeOfLevel'] == 'isobaricInhPa') :
                        level = level[(26-numMessages):] if (numMessages>1) else level

                        print "level = ",level
                        self.NAAPSmessages[msgKey] = self.getNAAPSlevelMessage(self,name,parameterName,shortName,units,parameterUnits,typeOfLevel,level,pressureUnits)
                        
                    else :

                        self.NAAPSmessages[msgKey] = self.getNAAPSmessage(self,name,parameterName,shortName,units,parameterUnits,typeOfLevel,level,pressureUnits)

                        if (typeOfLevel == 'heightAboveGround') :
                            print "msgList: ",msgList
                            height = str(msgList[msgListIdx]['level'])+'m'
                        elif (typeOfLevel == 'surface' or \
                              typeOfLevel == 'tropopause' or \
                              typeOfLevel == 'meanSea') :
                            height = str(msgList[msgListIdx]['typeOfLevel'])
                        elif (typeOfLevel == 'unknown' or typeOfLevel == 'entireAtmosphere') :
                            height = 'columnar'
                        else :
                            pass

                        self.NAAPSmessages[msgKey].height = height

            else :
                LOG.debug("There are %d messages for %s" % (len(msgList),msgKey))

            del(msgList)

            LOG.debug("###########################")

        gribFileObj.close()

        return 0


    def NAAPSgribToBlob_interp(gribObj,xmlFile,newNAAPSblob,endian=adl.LITTLE_ENDIAN):
        '''
        NAAPSgribToBlob_interp

        Method to copy the contents of the grib object gribObj to a newly
        created ADL blob file newNAAPSblob. If the grib lat/lon grid is coarser 
        than the required ADL grid, interpolation is performed.
        '''

        # Determine whether we need to do interpolation.
        needsInterp = False
        if (gribObj.dLon != gribObj.blob_Dlon) or (gribObj.dLat != gribObj.blob_Dlat) \
           or (gribObj.Nlons != gribObj.blob_Nlons) or (gribObj.Nlats != gribObj.blob_Nlats):
            needsInterp = True
            x = gribObj.Latitude[:,0]
            y = gribObj.Longitude[0,:]
            xnew = gribObj.blobLatitude[:,0]
            ynew = gribObj.blobLongitude[0,:]

        # Create a blank blob file from the xml file
        newBlobObj = adl.create(xmlFile,newNAAPSblob,endian=endian,overwrite=True)

        for msgKey in gribObj.gfsMsgKeys :
            print "\nCopying NAAPS[%s] to blob...\n" % (msgKey)
            msgObj = gribObj.NAAPSmessages[msgKey]
            blobArr = newBlobObj.__getattribute__(msgKey)

            try :
                if needsInterp :
                    z = gribObj.NAAPSmessages[msgKey].data
                    fRect = interpolate.RectBivariateSpline( x, y, z)
                    zNew = fRect(xnew,ynew)
                    for row in np.arange(gribObj.blob_Nlats):
                        blobArr[row][:] = zNew[gribObj.blob_Nlats-row-1,:]
                else :
                    for row in np.arange(gribObj.Nlats):
                        blobArr[row][:] = gribObj.NAAPSmessages[msgKey].data[row,:]

                print "Shape of %s blob array is %s" % (msgKey,repr(np.shape(blobArr)))

            except Exception, err:
                print "ERROR: %s" % (str(err))
                print "There was a problem assigning %s" % (msgKey)

    def NAAPSgribToBlob_interpNew(gribObj,xmlFile,newNAAPSblob,endian=adl.LITTLE_ENDIAN,reverseLat=False):
        '''
        NAAPSgribToBlob_interpNew

        Method to copy the contents of the grib object gribObj to a newly
        created ADL blob file newNAAPSblob. If the grib lat/lon grid is coarser 
        than the required ADL grid, interpolation is performed.

        This method employs proper numpy slicing of arrays rather than 
        "lists-of-lists" previously needed for arrays derived from ADL blobs.
        '''

        LOG.debug('Inside NAAPSgribToBlob_interpNew...')

        # Set latitude orientation in final output
        latDir = -1 if reverseLat else 1

        # Determine whether we need to do interpolation.
        needsInterp = False
        if (gribObj.dLon != gribObj.blob_Dlon) or (gribObj.dLat != gribObj.blob_Dlat) \
           or (gribObj.Nlons != gribObj.blob_Nlons) or (gribObj.Nlats != gribObj.blob_Nlats):
            needsInterp = True
            x = gribObj.Latitude#[:,0]
            y = gribObj.Longitude#[0,:]
            xnew = gribObj.blobLatitude[:,0]
            ynew = gribObj.blobLongitude[0,:]


        LOG.debug("Creating empty NAAPS blob file %s..." % (newNAAPSblob))
        newBlobObj = adl.create(xmlFile,newNAAPSblob,endian=endian,overwrite=True)
        newBlobArrObj = newBlobObj.as_arrays()

        for msgKey in gribObj.gfsMsgKeys :
            LOG.debug("Copying NAAPS[%s] to blob..." % (msgKey))
            msgObj = gribObj.NAAPSmessages[msgKey]
            blobArr = getattr(newBlobArrObj, msgKey)

            try :
                if needsInterp :
                    z = gribObj.NAAPSmessages[msgKey].data
                    fRect = interpolate.RectBivariateSpline( x, y, z)
                    zNew = fRect(xnew,ynew)
                    blobArr[:,:] = zNew[::latDir,:]
                else :
                    blobArr[:,:] = gribObj.NAAPSmessages[msgKey].data[::latDir,:]

                LOG.debug("Shape of %s blob array is %s" % (msgKey,repr(np.shape(blobArr))))

            except Exception, err:
                LOG.error("%s" % (str(err)))
                LOG.error("There was a problem assigning %s" % (msgKey))
                return -1

        return 0


###################################################
#                  Main Function                  #
###################################################

def main():

    #transChoices=['hdf5', 'binary']
    endianChoices=['big', 'little']

    description = '''Transcode an NAAPS GFS or GDAS grib file into the binary blob format required by Raytheon's Algorithm Development Library (ADL).'''
    usage = "usage: %prog [mandatory args] [options]"
    version = version="%prog "+__version__

    parser = optparse.OptionParser(description=description,usage=usage,version=version)

    # Mandatory arguments
    mandatoryGroup = optparse.OptionGroup(parser, "Mandatory Arguments",
                        "At a minimum these arguments must be specified...")

    mandatoryGroup.add_option('-g','--grib_file',
                      action="store",
                      dest="gribFile",
                      type="string",
                      help="The full path of the NAAPS grib file")

    parser.add_option_group(mandatoryGroup)

    # Optional arguments
    optionalGroup = optparse.OptionGroup(parser, "Extra Options",
                        "These options may be used to customize behaviour of this program.")

    optionalGroup.add_option('-x','--xml_file',
                      action="store",
                      dest="xmlFile" ,
                      type="string",
                      help="The full path of the ADL xml file describing the blob contents")

    optionalGroup.add_option('-l','--list',
                      action="store_true",
                      dest="isListing",
                      #default=False,
                      help="List the datasets contained in this grib file which match the required datasets for the ADL NAAPS blob file, and exit.")

    optionalGroup.add_option('-e','--endian',
                      action="store",
                      dest="endian",
                      default="little",
                      type="choice",
                      choices=endianChoices,
                      help="The endianess of the output blob file [default: %default]."+' Possible values are... %s' % (endianChoices.__str__()[0 :])
                      )

    optionalGroup.add_option('--reverse_latitude',
                      action="store_true",
                      dest="reverseLat",
                      default=False,
                      help="Reverse the latitude direction of the datasets.")

    optionalGroup.add_option('-o','--output_file',
                      action="store",
                      dest="outputFile",
                      default="NAAPS-ANC-Int_grib",
                      type="string",
                      help="The full path of the transcoded output file. [default: %default]")

    parser.add_option('-v', '--verbose',
                      dest='verbosity',
                      action="count",
                      default=0,
                      help='each occurrence increases verbosity 1 level from ERROR: -v=WARNING -vv=INFO -vvv=DEBUG')

    parser.add_option_group(optionalGroup)

    # Parse the arguments from the command line
    (options, args) = parser.parse_args()

    # Set up the logging levels
    levels = [logging.ERROR, logging.WARN, logging.INFO, logging.DEBUG]
    logging.basicConfig(level = levels[options.verbosity])

    # Check that all of the mandatory options are given. If one or more 
    # are missing, print error message and exit...
    mandatories = ['gribFile']
    mand_errors = [
                   "Missing mandatory argument [-grib GRIBFILE]"
                  ]
    isMissingMand = False
    for m,m_err in zip(mandatories,mand_errors):
        if not options.__dict__[m]:
            isMissingMand = True
            print m_err
    if isMissingMand :
        parser.error("Incomplete mandatory arguments, aborting...")

    # Expand the provided xml and grib file names into their fully qualified names
    if options.xmlFile :
        options.xmlFile = path.expanduser(options.xmlFile)
    if options.gribFile :
        options.gribFile = path.expanduser(options.gribFile)

    # Check that the grib file exists
    if not glob(options.gribFile) :
        parser.error("Grib file \n\t%s\ndoes not exist, aborting..." % (options.gribFile))
        
    if options.isListing :

        # Check that we don't have mutually exclusive options
        listError = False
        listErrorStr = "\n"
        if ( options.endian is not None ):
            listError = True
            print "We have endianness!"
            listErrorStr += "\toptions -l/--list and -e/--endian are mutually exclusive\n"
        if ( options.outputFile is not None ) :
            print "We have output file!"
            listError = True
            listErrorStr += "\toptions -l/--list and -o/--output_file  are mutually exclusive\n"

        if listError :
            parser.error(listErrorStr)

        NAAPSclass.NAAPSfileAudit(options.gribFile,levelInhPa=None)

    else : 

        # We want to transcode to blob, check that we have a valid xml file
        if options.xmlFile is None :
            parser.error("XML file required for transcoding to blob not provided, aborting...")

        if not glob(options.xmlFile) :
            parser.error("XML file \n\t%s\ndoes not exist, aborting..." % (options.xmlFile))
    
        # Create the grib object
        if options.outputFile is None :
            options.outputFile = "NAAPS-ANC-Int_grib"

        # Create the grib object and populate with the grib file data
        gribObj = NAAPSclass(gribFile = options.gribFile)

        # Write the contents of the gribObj object to an ADL blob file
        endian = adl.LITTLE_ENDIAN if options.endian=="little" else adl.BIG_ENDIAN
        #NAAPSclass.NAAPSgribToBlob_interp(gribObj,options.xmlFile,options.outputFile,endian=endian)
        NAAPSclass.NAAPSgribToBlob_interpNew(gribObj,options.xmlFile,options.outputFile,\
                endian=endian,reverseLat=options.reverseLat)

    LOG.info('Exiting...')
    return 0

if __name__=='__main__':
    sys.exit(main())

