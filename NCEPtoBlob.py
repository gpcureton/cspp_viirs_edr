#!/usr/bin/env python
# encoding: utf-8
"""
NCEPtoBlob.py

The purpose of this script is to ingest NCEP grib files and convert
to NCEP blob format required by Raytheon's Algorithm Development Library (ADL).

Minimum commandline(s)...

python NCEPtoBlob.py -g gribFile              # To audit a grib file...
python NCEPtoBlob.py -x xmlFile -g gribFile   # To transcode a grib file...

Created by Geoff Cureton <geoff.cureton@ssec.wisc.edu> on 2011-06-27.
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

import adl_blob
import adl_blob2
from thermo import rh_to_mr
rh_to_mr_vec = np.vectorize(rh_to_mr)

import optparse as optparse

# every module should have a LOG object
# e.g. LOG.warning('my dog has fleas')
import logging
LOG = logging.getLogger('NCEPtoBlob')


class NCEPclass(object):

    """ Static Attributes """

    NCEP_LAYER_LEVELS = {
                          '10mb' : 0,    '20mb' : 1,   '30mb' : 2,   '50mb' : 3,   '70mb' : 4,  '100mb' : 5,
                         '150mb' : 6,   '200mb' : 7,  '250mb' : 8,  '300mb' : 9,  '350mb' : 10, '400mb' : 11, 
                         '450mb' : 12,  '500mb' : 13, '550mb' : 14, '600mb' : 15, '650mb' : 16, '700mb' : 17, 
                         '750mb' : 18,  '800mb' : 19, '850mb' : 20, '900mb' : 21, '925mb' : 22, '950mb' : 23, 
                         '975mb' : 24, '1000mb' : 25
                        }

    NCEP_LAYER_VALUES = np.array([10.0, 20.0, 30.0, 50.0, 70.0, 100.0, 150.0, 200.0, 
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
    # used the ADL NCEP blob.

    gfsMsgKeys = [
        'geopotentialHeightLayers',
        'temperatureLayers',
        'waterVaporMixingRatioLayers',
        'pressureReducedToMSL',
        'uComponentOfWind',
        'vComponentOfWind',
        'surfacePressure',
        'skinTemperature',
        'surfaceTemperature',
        'totalPrecipitableWater',
        'surfaceGeopotentialHeight',
        'surfaceSpecificHumidity',
        'tropopauseGeopotentialHeight',
        'totalColumnOzone'
    ]

    ancDataSetTitle = {}
    ancDataSetTitle['geopotentialHeightLayers'] = 'geopotentialHeightLayers'
    ancDataSetTitle['temperatureLayers'] = 'temperatureLayers'
    ancDataSetTitle['waterVaporMixingRatioLayers'] = ['waterVaporMixingRatioLayers','waterVaporMixingRatioLayers']
    ancDataSetTitle['pressureReducedToMSL'] = 'pressureReducedToMSL'
    ancDataSetTitle['uComponentOfWind'] = 'uComponentOfWind'
    ancDataSetTitle['vComponentOfWind'] = 'vComponentOfWind'
    ancDataSetTitle['surfacePressure'] = 'surfacePressure'
    ancDataSetTitle['skinTemperature'] = 'skinTemperature'
    ancDataSetTitle['surfaceTemperature'] = 'surfaceTemperature'
    ancDataSetTitle['totalPrecipitableWater'] = 'totalPrecipitableWater'
    ancDataSetTitle['surfaceGeopotentialHeight'] = 'surfaceGeopotentialHeight'
    ancDataSetTitle['surfaceSpecificHumidity'] = 'surfaceSpecificHumidity'
    ancDataSetTitle['tropopauseGeopotentialHeight'] = 'tropopauseGeopotentialHeight'
    ancDataSetTitle['totalColumnOzone'] = ['totalColumnOzone','totalColumnOzone']

    ancName = {}
    ancName['geopotentialHeightLayers'] = 'Geopotential Height'
    ancName['temperatureLayers'] = 'Temperature'
    ancName['waterVaporMixingRatioLayers'] = ['Relative humidity','Specific humidity']
    ancName['pressureReducedToMSL'] = ['Pressure reduced to MSL','Mean sea level pressure']
    ancName['uComponentOfWind'] = '10 metre U wind component'
    ancName['vComponentOfWind'] = '10 metre V wind component'
    ancName['surfacePressure'] = 'Surface pressure'
    ancName['skinTemperature'] = 'Temperature'
    ancName['surfaceTemperature'] = '2 metre temperature'
    ancName['totalPrecipitableWater'] = 'Precipitable water'
    ancName['surfaceGeopotentialHeight'] = 'Orography'
    ancName['surfaceSpecificHumidity'] = 'Specific humidity'
    ancName['tropopauseGeopotentialHeight'] = 'Geopotential Height'
    ancName['totalColumnOzone'] = ['Total ozone','Total column ozone']

    ancParameterName = {}
    ancParameterName['geopotentialHeightLayers'] = ['Geopotential height','GH Geopotential height gpm']
    ancParameterName['temperatureLayers'] = ['Temperature','T Temperature K']
    ancParameterName['waterVaporMixingRatioLayers'] = ['Relative humidity','Specific humidity','R Relative humidity %']
    ancParameterName['pressureReducedToMSL'] = ['Pressure reduced to MSL','MSL Mean sea level pressure Pa']
    ancParameterName['uComponentOfWind'] = ['10 metre U wind component','U U-component of wind m s**-1']
    ancParameterName['vComponentOfWind'] = ['10 metre V wind component','V V-component of wind m s**-1']
    ancParameterName['surfacePressure'] = ['Surface pressure','P Pressure Pa']
    ancParameterName['skinTemperature'] = ['Temperature','T Temperature K']
    ancParameterName['surfaceTemperature'] = ['Temperature','T Temperature K']
    ancParameterName['totalPrecipitableWater'] = ['Precipitable water','None Precipitable water kg m**-2']
    ancParameterName['surfaceGeopotentialHeight'] = ['Geopotential height','GH Geopotential height gpm']
    ancParameterName['surfaceSpecificHumidity'] = ['Specific humidity','Q Specific humidity kg kg**-1']
    ancParameterName['tropopauseGeopotentialHeight'] = ['Geopotential height','GH Geopotential height gpm']
    ancParameterName['totalColumnOzone'] = ['Total ozone','TCO3 Total']

    ancShortName = {}
    ancShortName['geopotentialHeightLayers'] = 'gh'
    ancShortName['temperatureLayers'] = 't'
    ancShortName['waterVaporMixingRatioLayers'] = ['r','q']
    ancShortName['pressureReducedToMSL'] = ['prmsl','msl']
    ancShortName['uComponentOfWind'] = '10u'
    ancShortName['vComponentOfWind'] = '10v'
    ancShortName['surfacePressure'] = 'sp'
    ancShortName['skinTemperature'] = 't'
    ancShortName['surfaceTemperature'] = '2t'
    ancShortName['totalPrecipitableWater'] = 'pwat'
    ancShortName['surfaceGeopotentialHeight'] = 'orog'
    ancShortName['surfaceSpecificHumidity'] = 'qv_2m'
    ancShortName['tropopauseGeopotentialHeight'] = 'gh'
    ancShortName['totalColumnOzone'] = ['tozne','tco3']

    ancParameterUnits = {}
    ancParameterUnits['geopotentialHeightLayers'] = ['gpm','unknown']
    ancParameterUnits['temperatureLayers'] = ['K','unknown']
    ancParameterUnits['waterVaporMixingRatioLayers'] = ['%','kg kg-1','unknown']
    ancParameterUnits['pressureReducedToMSL'] = ['Pa','unknown']
    ancParameterUnits['uComponentOfWind'] = ['m s-1','unknown']
    ancParameterUnits['vComponentOfWind'] = ['m s-1','unknown']
    ancParameterUnits['surfacePressure'] = ['Pa','unknown']
    ancParameterUnits['skinTemperature'] = ['K','unkknown']
    ancParameterUnits['surfaceTemperature'] = ['K','unkknown']
    ancParameterUnits['totalPrecipitableWater'] = ['kg m-2','unknown']
    ancParameterUnits['surfaceGeopotentialHeight'] = ['gpm','unknown']
    ancParameterUnits['surfaceSpecificHumidity'] = ['kg kg-1','unknown']
    ancParameterUnits['tropopauseGeopotentialHeight'] = ['gpm','unknown']
    ancParameterUnits['totalColumnOzone'] = ['Dobson','column']

    ancUnits = {}
    ancUnits['geopotentialHeightLayers'] = 'gpm'
    ancUnits['temperatureLayers'] = 'K'
    ancUnits['waterVaporMixingRatioLayers'] = ['%','kg kg**-1']
    ancUnits['pressureReducedToMSL'] = 'Pa'
    ancUnits['uComponentOfWind'] = 'm s**-1'
    ancUnits['vComponentOfWind'] = 'm s**-1'
    ancUnits['surfacePressure'] = 'Pa'
    ancUnits['skinTemperature'] = 'K'
    ancUnits['surfaceTemperature'] = 'K'
    ancUnits['totalPrecipitableWater'] = 'kg m**-2'
    ancUnits['surfaceGeopotentialHeight'] = 'm'
    ancUnits['surfaceSpecificHumidity'] = 'kg kg**-1'
    ancUnits['tropopauseGeopotentialHeight'] = 'gpm'
    ancUnits['totalColumnOzone'] = ['Dobson','kg m**-2']

    ancTypeOfLevel = {}
    ancTypeOfLevel['geopotentialHeightLayers'] = 'isobaricInhPa'
    ancTypeOfLevel['temperatureLayers'] = 'isobaricInhPa'
    ancTypeOfLevel['waterVaporMixingRatioLayers'] = 'isobaricInhPa'
    ancTypeOfLevel['pressureReducedToMSL'] = 'meanSea'
    ancTypeOfLevel['uComponentOfWind'] = 'heightAboveGround'
    ancTypeOfLevel['vComponentOfWind'] = 'heightAboveGround'
    ancTypeOfLevel['surfacePressure'] = 'surface'
    ancTypeOfLevel['skinTemperature'] = 'surface'
    ancTypeOfLevel['surfaceTemperature'] = 'heightAboveGround'
    ancTypeOfLevel['totalPrecipitableWater'] = ['unknown','entireAtmosphere']
    ancTypeOfLevel['surfaceGeopotentialHeight'] = 'surface'
    ancTypeOfLevel['surfaceSpecificHumidity'] = 'heightAboveGround'
    ancTypeOfLevel['tropopauseGeopotentialHeight'] = 'tropopause'
    ancTypeOfLevel['totalColumnOzone'] = ['unknown','entireAtmosphere']

    ancLevel = {}
    ancLevel['geopotentialHeightLayers'] = NCEP_LAYER_VALUES
    ancLevel['temperatureLayers'] = NCEP_LAYER_VALUES
    ancLevel['waterVaporMixingRatioLayers'] = NCEP_LAYER_VALUES
    ancLevel['pressureReducedToMSL'] = 0
    ancLevel['uComponentOfWind'] = 10
    ancLevel['vComponentOfWind'] = 10
    ancLevel['surfacePressure'] = 0
    ancLevel['skinTemperature'] = 0
    ancLevel['surfaceTemperature'] = 2
    ancLevel['totalPrecipitableWater'] = 0
    ancLevel['surfaceGeopotentialHeight'] = 0
    ancLevel['surfaceSpecificHumidity'] = 2
    ancLevel['tropopauseGeopotentialHeight'] = 0
    ancLevel['totalColumnOzone'] = 0


    def __init__(self,gribFile=None,isListing=False):
        if gribFile is not None :
            gribFile = path.abspath(path.expanduser(gribFile))

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

            self.NCEPmessages = {}

            returnObj = not isListing

            if not returnObj :
                LOG.info("Auditing %s\n" % (gribFile))

            self.__NCEPaudit(gribFile,returnObj=returnObj,levelInhPa=None)

        else :
            pass

    class getNCEPlevelMessage(object):
        '''
        getNCEPlevelMessage

        Instantiations of the getNCEPlevelMessage class serve as containers for multiple 
        getNCEPmessage objects making up a profile, and have as attributes various properties 
        common to each pressure level. The getNCEPmessage object for each level contains the array 
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
                self.messageLevelData[pressureLevel] = NCEPclass.getNCEPmessage(self,name,parameterName,shortName,units,parameterUnits,typeOfLevel,pressureLevel,pressureUnits)
                self.messageLevelData[pressureLevel].height = "%s %s" %(pressureLevel,self.pressureUnits)

    class getNCEPmessage(object):
        '''
        getNCEPmessage

        Instantiations of the getNCEPmessage class contain the array data for a single GRIB message,
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
    def NCEPfileAudit(gribFile,levelInhPa=None):
        '''
        NCEPfileAudit

        This static method lists a manifest of the supplied grib file's messages 
        which match the required ADL blob datasets.
        '''
        gribFile = path.abspath(path.expanduser(gribFile))
        NCEPobj = NCEPclass(gribFile=gribFile, isListing=True)

    def __NCEPaudit(self,gribFile,returnObj=False,levelInhPa=None) :
        '''
        __NCEPaudit

        This private method opens the supplied grib file, and populates a dictionary with getNCEPlevelMessage 
        and getNCEPmessage objects matching the required ADL blob datasets.

        If returnObj is False, just lists the mathing grib messages and exits. 
        '''

        gribFileObj = pygrib.open(gribFile)

        gfsMsgKeys = NCEPclass.gfsMsgKeys

        gfs_Msg = {}
        plotData = {}

        # Populate the gfs_Msg dict
        for msgKey in gfsMsgKeys :

            # Retrieve some identifying information about this message type
            name = NCEPclass.ancName[msgKey]
            shortName = NCEPclass.ancShortName[msgKey]
            typeOfLevel = NCEPclass.ancTypeOfLevel[msgKey]

            if (typeOfLevel == 'isobaricInhPa') :
                level = NCEPclass.ancLevel[msgKey] if (levelInhPa==None) else levelInhPa
            else :
                level = NCEPclass.ancLevel[msgKey]

            LOG.debug("Ingesting(%s) : %s for typeOfLevel %s...\n... for level %s\n" % \
                    (msgKey,repr(name),repr(typeOfLevel),repr(level)))
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

                #sys.exit(0)

                LOG.debug("There are %d messages for %s with...\n\tname = \"%s\"\n\tparameterName = \"%s\"\n\tshortName = %s\n\tparameterUnits = %s\n\tunits = %s\n\ttypeOfLevel = %s\n\tpressureUnits = %s\n\tlevel = %s" \
                    % (len(msgList), msgKey, name, parameterName, shortName, 
                            parameterUnits, units, typeOfLevel, pressureUnits, level
                      ))

                # TODO : If we just want a list we should stop here...
                if not returnObj :
                    pass
                    #LOG.info("Auditing %s\n" % (gribFile))
                else :

                    # Determine the level index, name and height of the dataset
                    if (msgList[msgListIdx]['typeOfLevel'] == 'isobaricInhPa') :
                        level = level[(26-numMessages):] if (numMessages>1) else level

                        LOG.debug("level = %r" % (level))
                        self.NCEPmessages[msgKey] = self.getNCEPlevelMessage(self,name,parameterName,shortName,units,parameterUnits,typeOfLevel,level,pressureUnits)
                        
                    else :

                        self.NCEPmessages[msgKey] = self.getNCEPmessage(self,name,parameterName,shortName,units,parameterUnits,typeOfLevel,level,pressureUnits)

                        if (typeOfLevel == 'heightAboveGround') :
                            LOG.debug("msgList: " % (msgList))
                            height = str(msgList[msgListIdx]['level'])+'m'
                        elif (typeOfLevel == 'surface' or \
                              typeOfLevel == 'tropopause' or \
                              typeOfLevel == 'meanSea') :
                            height = str(msgList[msgListIdx]['typeOfLevel'])
                        elif (typeOfLevel == 'unknown' or typeOfLevel == 'entireAtmosphere') :
                            height = 'columnar'
                        else :
                            pass

                        self.NCEPmessages[msgKey].height = height

            else :
                LOG.debug("There are %d messages for %s" % (len(msgList),msgKey))

            del(msgList)

            LOG.debug("###########################")

        gribFileObj.close()

        return 0


    def NCEPgribToBlob(gribObj,xmlFile,newNCEPblob):
        '''
        NCEPgribToBlob

        Method to copy the contents of the grib object gribObj to a newly
        created ADL blob file newNCEPblob.
        '''

        # Create a blank blob file from the xml file
        #newBlobObj = adl_blob.create(xmlFile,newNCEPblob,endian=adl_blob.LITTLE_ENDIAN,overwrite=True)
        newBlobObj = adl_blob2.create(xmlFile,newNCEPblob,endian=adl_blob2.LITTLE_ENDIAN,overwrite=True)

        #newBlobArrObj = newBlobObj.as_arrays()

        for msgKey in gribObj.gfsMsgKeys :
            LOG.debug("\nCopying NCEP[%s] to blob...\n" % (msgKey))
            msgObj = gribObj.NCEPmessages[msgKey]
            #blobArr = getattr(newBlobArrObj, msgKey)
            blobArr = getattr(newBlobObj, msgKey)

            if msgObj.typeOfLevel=='isobaricInhPa' :
                for level,levelIdx in gribObj.NCEP_LAYER_LEVELS.items() : 
                    levelStrIdx = level[:-2]
                    try :
                        for row in np.arange(gribObj.Nlats):
                            blobArr[levelIdx][row][:] = gribObj.NCEPmessages[msgKey].messageLevelData[levelStrIdx].data[gribObj.Nlats-row-1,:]

                    except Exception, err:
                        LOG.error("%s" % (str(err)))
                        LOG.error("There was a problem assigning layer %s of %s" % (level,msgKey))
            else :
                try :
                    for row in np.arange(gribObj.Nlats):
                        blobArr[row][:] = gribObj.NCEPmessages[msgKey].data[gribObj.Nlats-row-1,:]
                    LOG.debug("Shape of %s blob array is %s" % (msgKey,repr(np.shape(blobArr))))

                except Exception, err:
                    LOG.error("%s" % (str(err)))
                    LOG.error("There was a problem assigning %s" % (msgKey))

    def NCEPgribToBlob_interp(gribObj,xmlFile,newNCEPblob,endian=adl_blob2.LITTLE_ENDIAN):
        '''
        NCEPgribToBlob_interp

        Method to copy the contents of the grib object gribObj to a newly
        created ADL blob file newNCEPblob. If the grib lat/lon grid is coarser 
        than the required ADL grid, interpolation is performed.
        '''

        # Determine whether we need to do interpolation.
        needsInterp = False
        if (gribObj.dLon != gribObj.blob_Dlon) or (gribObj.dLat != gribObj.blob_Dlat) \
           or (gribObj.Nlons != gribObj.blob_Nlons) or (gribObj.Nlats != gribObj.blob_Nlats):
            needsInterp = True
            x = gribObj.Latitude#[:,0]
            y = gribObj.Longitude#[0,:]
            xnew = gribObj.blobLatitude[:,0]
            ynew = gribObj.blobLongitude[0,:]

        # Create a blank blob file from the xml file
        #newBlobObj = adl_blob.create(xmlFile,newNCEPblob,endian=endian,overwrite=True)
        newBlobObj = adl_blob2.create(xmlFile,newNCEPblob,endian=endian,overwrite=True)

        #newBlobArrObj = newBlobObj.as_arrays()

        for msgKey in gribObj.gfsMsgKeys :
            LOG.debug("\nCopying NCEP[%s] to blob...\n" % (msgKey))
            msgObj = gribObj.NCEPmessages[msgKey]
            #blobArr = getattr(newBlobArrObj, msgKey)
            blobArr = getattr(newBlobObj, msgKey)

            if msgObj.typeOfLevel=='isobaricInhPa' :
                for level,levelIdx in gribObj.NCEP_LAYER_LEVELS.items() : 
                    levelStrIdx = level[:-2]
                    try :
                        if needsInterp :
                            z = gribObj.NCEPmessages[msgKey].messageLevelData[levelStrIdx].data
                            fRect = interpolate.RectBivariateSpline( x, y, z)
                            zNew = fRect(xnew,ynew)
                            for row in np.arange(gribObj.blob_Nlats):
                                blobArr[levelIdx][row][:] = zNew[gribObj.blob_Nlats-row-1,:]
                        else :
                            for row in np.arange(gribObj.blob_Nlats):
                                blobArr[levelIdx][row][:] = gribObj.NCEPmessages[msgKey].messageLevelData[levelStrIdx].data[gribObj.blob_Nlats-row-1,:]

                    except Exception, err:
                        LOG.error("%s" % (str(err)))
                        LOG.error("There was a problem assigning layer %s of %s" % (level,msgKey))
                pass
            else :
                try :
                    if needsInterp :
                        z = gribObj.NCEPmessages[msgKey].data
                        fRect = interpolate.RectBivariateSpline( x, y, z)
                        zNew = fRect(xnew,ynew)
                        for row in np.arange(gribObj.blob_Nlats):
                            blobArr[row][:] = zNew[gribObj.blob_Nlats-row-1,:]
                    else :
                        for row in np.arange(gribObj.Nlats):
                            blobArr[row][:] = gribObj.NCEPmessages[msgKey].data[gribObj.blob_Nlats-row-1,:]

                except Exception, err:
                    LOG.error("%s" % (str(err)))
                    LOG.error("There was a problem assigning %s" % (msgKey))

            LOG.debug("\tShape of %s blob array is %s" % (msgKey,repr(np.shape(blobArr))))


    def NCEPgribToBlob_interpNew(gribObj,xmlFile,newNCEPblob,endian=adl_blob2.LITTLE_ENDIAN,reverseLat=False):
        '''
        NCEPgribToBlob_interpNew

        Method to copy the contents of the grib object gribObj to a newly
        created ADL blob file newNCEPblob. If the grib lat/lon grid is coarser 
        than the required ADL grid, interpolation is performed.

        This method employs proper numpy slicing of arrays rather than 
        "lists-of-lists" previously needed for arrays derived from ADL blobs.
        '''

        # Set latitude orientation in final output
        latDir = 1 if reverseLat else -1

        # Determine whether we need to do interpolation.
        needsInterp = False
        if (gribObj.dLon != gribObj.blob_Dlon) or (gribObj.dLat != gribObj.blob_Dlat) \
           or (gribObj.Nlons != gribObj.blob_Nlons) or (gribObj.Nlats != gribObj.blob_Nlats):
            needsInterp = True
            x = gribObj.Latitude#[:,0]
            y = gribObj.Longitude#[0,:]
            xnew = gribObj.blobLatitude[:,0]
            ynew = gribObj.blobLongitude[0,:]


        LOG.debug("Creating empty NCEP blob file %s..." % (newNCEPblob))

        try :
            #newBlobObj = adl_blob.create(xmlFile,newNCEPblob,endian=endian,overwrite=True)
            newBlobObj = adl_blob2.create(xmlFile,newNCEPblob,endian=endian,overwrite=True)
        except Exception, err:
                    LOG.error("%s" % (str(err)))
                    LOG.error("Blob creation error %s" % (newNCEPblob))
                    return -1


        #newBlobArrObj = newBlobObj.as_arrays()

        for msgKey in gribObj.gfsMsgKeys :
            LOG.debug("Copying NCEP[%s] to blob...\n" % (msgKey))
            msgObj = gribObj.NCEPmessages[msgKey]
            #blobArr = getattr(newBlobArrObj, msgKey)
            blobArr = getattr(newBlobObj, msgKey)

            if msgObj.typeOfLevel=='isobaricInhPa' :
                for level,levelIdx in gribObj.NCEP_LAYER_LEVELS.items() : 
                    levelStrIdx = level[:-2]
                    try :
                        if needsInterp :
                            LOG.debug("%r" % (x))
                            z = gribObj.NCEPmessages[msgKey].messageLevelData[levelStrIdx].data
                            fRect = interpolate.RectBivariateSpline( x, y, z)
                            zNew = fRect(xnew,ynew)
                            blobArr[levelIdx][:,:] = zNew[::latDir,:]
                        else :
                            blobArr[levelIdx,:,:] = gribObj.NCEPmessages[msgKey].messageLevelData[levelStrIdx].data[::latDir,:]

                    except Exception, err:
                        LOG.error("%s" % (str(err)))
                        LOG.error("There was a problem assigning layer %s of %s" % (level,msgKey))
                        return -1

            else :
                try :
                    if needsInterp :
                        z = gribObj.NCEPmessages[msgKey].data
                        fRect = interpolate.RectBivariateSpline( x, y, z)
                        zNew = fRect(xnew,ynew)
                        blobArr[:,:] = zNew[::latDir,:]
                    else :
                        blobArr[:,:] = gribObj.NCEPmessages[msgKey].data[::latDir,:]

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

    description = '''Transcode an NCEP GFS or GDAS grib file into the binary 
    blob format required by Raytheon's Algorithm Development Library (ADL).'''
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
                      help="The full path of the NCEP grib file")

    parser.add_option_group(mandatoryGroup)

    # Optional arguments
    optionalGroup = optparse.OptionGroup(parser, "Extra Options",
                        "These options may be used to customize behaviour of this program.")

    optionalGroup.add_option('-x','--xml_file',
                      action="store",
                      dest="xmlFile" ,
                      type="string",
                      help="""The full path of the ADL xml file describing the 
                      blob contents""")

    optionalGroup.add_option('-l','--list',
                      action="store_true",
                      dest="isListing",
                      #default=False,
                      help="""List the datasets contained in this grib file which 
match the required datasets for the ADL NCEP blob file, and 
exit.""")

    optionalGroup.add_option('-e','--endian',
                      action="store",
                      dest="endian",
                      default="little",
                      type="choice",
                      choices=endianChoices,
                      help="""The endianess of the output blob file [default: %default]. 
Possible values are... {}""".format(endianChoices)
                      )

    optionalGroup.add_option('--reverse_latitude',
                      action="store_true",
                      dest="reverseLat",
                      default=False,
                      help="Reverse the latitude direction of the datasets.")

    optionalGroup.add_option('-o','--output_file',
                      action="store",
                      dest="outputFile",
                      default="NCEP-ANC-Int_grib",
                      type="string",
                      help="""The full path of the transcoded output file. 
                      [default: %default]""")

    parser.add_option('-v', '--verbose',
                      dest='verbosity',
                      action="count",
                      default=0,
                      help="""each occurrence increases verbosity 1 level 
                      from ERROR: -v=WARNING -vv=INFO -vvv=DEBUG""")

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
            LOG.error(m_err)
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
            LOG.debug("We have endianness!")
            listErrorStr += "\toptions -l/--list and -e/--endian are mutually exclusive\n"
        if ( options.outputFile is not None ) :
            LOG.debug("We have output file!")
            listError = True
            listErrorStr += "\toptions -l/--list and -o/--output_file  are mutually exclusive\n"

        if listError :
            parser.error(listErrorStr)

        NCEPclass.NCEPfileAudit(options.gribFile,levelInhPa=None)

    else : 

        # We want to transcode to blob, check that we have a valid xml file
        if options.xmlFile is None :
            parser.error("XML file required for transcoding to blob not provided, aborting...")

        if not glob(options.xmlFile) :
            parser.error("XML file \n\t%s\ndoes not exist, aborting..." % (options.xmlFile))
    
        # Create the grib object
        if options.outputFile is None :
            options.outputFile = "NCEP-ANC-Int_grib"

        # Create the grib object and populate with the grib file data
        gribObj = NCEPclass(gribFile = options.gribFile)

        # Scale/convert various datasets to match ADL spec...
        # TODO : This should be in a separate method...

        # Convert surface pressure from Pa to mb or hPa ...
        # Ref: ADL/CMN/Utilities/ING/MSD/NCEP/src/IngMsdNCEP_Converter.cpp
        # Ref: applyScalingFactor(currentBuffer, GRID_SIZE, 0.01);
        gribObj.NCEPmessages['surfacePressure'].data /= 100.

        # Convert total column ozone from DU or kg m**-2 to Atm.cm ...
        # Ref: ADL/CMN/Utilities/ING/MSD/NCEP/src/IngMsdNCEP_Converter.cpp
        # Code: const float DOBSON_TO_ATMSCM_SCALING_FACTOR = .001;
        # Code: applyScalingFactor(currentBuffer, GRID_SIZE,DOBSON_TO_ATMSCM_SCALING_FACTOR);
        gribObj.NCEPmessages['totalColumnOzone'].data /= 1000.

        # Convert total precipitable water kg m^{-2} to cm ...
        # Ref: ADL/CMN/Utilities/ING/MSD/NCEP/src/IngMsdNCEP_Converter.cpp
        # Code: applyScalingFactor(currentBuffer, GRID_SIZE, .10);
        gribObj.NCEPmessages['totalPrecipitableWater'].data /= 10.

        # Convert specific humidity in kg.kg^{-1} to water vapor mixing ratio in g.kg^{-1}
        # Ref: ADL/CMN/Utilities/ING/MSD/NCEP/src/IngMsdNCEP_Converter.cpp
        # Code: void IngMsdNCEP_Converter::applyWaterVaporMixingRatio()
        # Code: destination[i] = 1000 * (destination[i]/ (1-destination[i]));
        moistureObj = gribObj.NCEPmessages['waterVaporMixingRatioLayers'].messageLevelData
        temperatureObj = gribObj.NCEPmessages['temperatureLayers'].messageLevelData

        # Compute the 100mb mixing ratio in g/kg
        if moistureObj['100'].name == 'Specific humidity':
            specHumidity_100mb = moistureObj['100'].data
            mixingRatio_100mb = 1000. * specHumidity_100mb/(1. - specHumidity_100mb)
        elif  moistureObj['100'].name == 'Relative humidity':
            relativeHumidity_100mb = moistureObj['100'].data
            temperature_100mb = temperatureObj['100'].data
            mixingRatio_100mb = rh_to_mr_vec(relativeHumidity_100mb,100.,temperature_100mb)
        else :
            pass

        for level,levelIdx in gribObj.NCEP_LAYER_LEVELS.items() : 
            levelStrIdx = level[:-2]
            pressure = gribObj.NCEP_LAYER_VALUES[levelIdx]

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

        # Write the contents of the gribObj object to an ADL blob file
        endian = adl_blob.LITTLE_ENDIAN if options.endian=="little" else adl_blob.BIG_ENDIAN
        #NCEPclass.NCEPgribToBlob_interp(gribObj,options.xmlFile,options.outputFile,endian=endian)
        NCEPclass.NCEPgribToBlob_interpNew(gribObj,options.xmlFile,options.outputFile,\
                endian=endian,reverseLat=options.reverseLat)

    LOG.info('Exiting...')
    return 0

if __name__=='__main__':
    sys.exit(main())

