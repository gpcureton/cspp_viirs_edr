#!/usr/bin/env python
# encoding: utf-8
"""
CSPP_VIIRS_NCEP-ANC-Int_blobPlots.py

Purpose: Create PNGs from the NCEP global blob files.

Input:
    * Various inputs.

Output:
    * None

Details:
    * None

Preconditions:
    * None

Optional:
    * 

Minimum commandline...

    export CSPP_EDR_HOME=$(readlink -f /path/to/EDR)
    . $CSPP_EDR_HOME/cspp_edr_env.sh

    python CSPP_VIIRS_NCEP-ANC-Int_blobPlots.py -b '/path/to/blob/file' -x '/path/to/xml/file'


Created by Geoff Cureton <geoff.cureton@ssec.wisc.edu> on 2014-06-30.
Copyright (c) 2014 University of Wisconsin Regents. All rights reserved.

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

#############
import os, sys
from os import path, uname, mkdir
from glob import glob
import string, traceback
from time import time
import re

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

from mpl_toolkits.basemap import Basemap,addcyclic,shiftgrid

import optparse as optparse

#import adl_blob
import adl_blob2 as adl_blob
from ViirsData import ViirsTrimTable
import viirs_edr_data

# every module should have a LOG object
# e.g. LOG.warning('my dog has fleas')
import logging
LOG = logging.getLogger(__file__)


class ANCclass():

    def __init__(self):

        self.data_names = [
                               'geopotentialHeightLayers',
                               'temperatureLayers',
                               'waterVaporMixingRatioLayers',
                               'totalPrecipitableWater',
                               'surfaceTemperature',
                               'uComponentOfWind',
                               'vComponentOfWind',
                               'surfacePressure',
                               'pressureReducedToMSL',
                               'totalColumnOzone',
                               'surfaceGeopotentialHeight',
                               'tropopauseGeopotentialHeight',
                               'surfaceSpecificHumidity',
                               'skinTemperature'
                              ]


        self.ancDataSetTitle = {}
        self.ancDataSetTitle['geopotentialHeightLayers'] = 'geopotentialHeightLayers'
        self.ancDataSetTitle['temperatureLayers'] = 'temperatureLayers'
        self.ancDataSetTitle['waterVaporMixingRatioLayers'] = 'waterVaporMixingRatioLayers'
        self.ancDataSetTitle['totalPrecipitableWater'] = 'totalPrecipitableWater'
        self.ancDataSetTitle['surfaceTemperature'] = 'surfaceTemperature'
        self.ancDataSetTitle['uComponentOfWind'] = 'uComponentOfWind'
        self.ancDataSetTitle['vComponentOfWind'] = 'vComponentOfWind'
        self.ancDataSetTitle['surfacePressure'] = 'surfacePressure'
        self.ancDataSetTitle['pressureReducedToMSL'] = 'pressureReducedToMSL'
        self.ancDataSetTitle['totalColumnOzone'] = 'totalColumnOzone'
        self.ancDataSetTitle['surfaceGeopotentialHeight'] = 'surfaceGeopotentialHeight'
        self.ancDataSetTitle['tropopauseGeopotentialHeight'] = 'tropopauseGeopotentialHeight'
        self.ancDataSetTitle['surfaceSpecificHumidity'] = 'surfaceSpecificHumidity'
        self.ancDataSetTitle['skinTemperature'] = 'skinTemperature'

        self.NCEP_LAYER_LEVELS = {
                              '10mb' : 0,    '20mb' : 1,   '30mb' : 2,   '50mb' : 3,   '70mb' : 4,  '100mb' : 5,
                             '150mb' : 6,   '200mb' : 7,  '250mb' : 8,  '300mb' : 9,  '350mb' : 10, '400mb' : 11, 
                             '450mb' : 12,  '500mb' : 13, '550mb' : 14, '600mb' : 15, '650mb' : 16, '700mb' : 17, 
                             '750mb' : 18,  '800mb' : 19, '850mb' : 20, '900mb' : 21, '925mb' : 22, '950mb' : 23, 
                             '975mb' : 24, '1000mb' : 25
                            }

        self.NCEP_LAYER_VALUES = np.array([10.0, 20.0, 30.0, 50.0, 70.0, 100.0, 150.0, 200.0, 
                                      250.0, 300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 
                                      600.0, 650.0, 700.0, 750.0, 800.0, 850.0, 900.0, 
                                      925.0, 950.0, 975.0, 1000.0 ]);

        self.plotDescr = {}
        self.plotDescr['geopotentialHeightLayers'] = 'Geopotential Height (gpm)'
        self.plotDescr['temperatureLayers'] = r'Temperature (K)'
        self.plotDescr['waterVaporMixingRatioLayers'] = 'Water Vapor Mixing Ratio (g kg$^{-1}$)'
        self.plotDescr['totalPrecipitableWater']   = r'Total Precipitable Water (cm, $\times$10 kg m$^{-2}$)'
        self.plotDescr['surfaceTemperature'] = r'Air Temperature @ 2m (K)'
        self.plotDescr['uComponentOfWind']  = r'Surface Wind velocity (ms$^{-1}$)'
        self.plotDescr['vComponentOfWind']  = r'Surface Wind velocity (ms$^{-1}$)'
        self.plotDescr['surfacePressure']  = r'Surface Pressure (mb, $\times$10$^{2}$ Pa)'
        self.plotDescr['pressureReducedToMSL']  = r'Pressure (Pa)'
        self.plotDescr['totalColumnOzone']     = r'Total Column Ozone (Atm cm, $\times$10$^{3}$DU)'
        self.plotDescr['surfaceGeopotentialHeight'] = 'Surface Geopotential Height (gpm)'
        self.plotDescr['tropopauseGeopotentialHeight'] = 'Tropopause Geopotential Height (gpm)'
        self.plotDescr['surfaceSpecificHumidity'] = 'Specific Humidity (kg kg$^{-1}$)'
        self.plotDescr['skinTemperature'] = 'Sea Surface Skin Temperature (K)'

        self.plotLims = {}
        self.plotLims['geopotentialHeightLayers'] = [None,None]
        self.plotLims['temperatureLayers'] = [None,None]
        self.plotLims['waterVaporMixingRatioLayers'] = [None,None]
        self.plotLims['totalPrecipitableWater'] = [0.,10.]
        self.plotLims['surfaceTemperature'] = [275.,315.]
        self.plotLims['uComponentOfWind'] = [0.,20.]
        self.plotLims['vComponentOfWind'] = [0.,20.]
        self.plotLims['surfacePressure'] = [800.,1013.]
        self.plotLims['pressureReducedToMSL'] = [80000.,101300.]
        self.plotLims['totalColumnOzone'] = [0.2,0.40]
        self.plotLims['surfaceGeopotentialHeight'] = [None,None]
        self.plotLims['tropopauseGeopotentialHeight'] = [None,None]
        self.plotLims['surfaceSpecificHumidity'] = [0.002,0.04]
        self.plotLims['skinTemperature'] = [275.,315.]


    def plotAncData(self,gridLat,gridLon,gridData,plotData) :

        ancType   =  plotData['ancType']
        plotTitle =  plotData['plotTitle']
        blobName  =  plotData['NCEP_ANC_BlobFile']
        cbTitle   =  plotData['cbTitle']
        vmin,vmax =  plotData['plotLims'][0], plotData['plotLims'][1]
        dpi       =  plotData['dpi']
        prefix       =  plotData['prefix']

        # Create figure with default size, and create canvas to draw on
        scale=1.5
        fig = Figure(figsize=(scale*8,scale*5))
        canvas = FigureCanvas(fig)

        # Create main axes instance, leaving room for colorbar at bottom,
        # and also get the Bbox of the axes instance
        ax_rect = [0.05, 0.18, 0.9, 0.75  ] # [left,bottom,width,height]
        ax = fig.add_axes(ax_rect)

        # Granule axis title
        LOG.info("plotTitle = {}".format(plotTitle))
        ax_title = ppl.setp(ax,title=plotTitle)
        ppl.setp(ax_title,fontsize=10)
        ppl.setp(ax_title,family="sans-serif")

        # Create the basemap object
        llcrnrlon = self.lonMin
        llcrnrlat = self.latMin
        urcrnrlon = self.lonMax
        urcrnrlat = self.latMax
         
        lon_0 = 0.
        gridData, gridLon = shiftgrid(lon_0, gridData, gridLon[0], start=True)

        if llcrnrlon != None and llcrnrlat != None and \
                urcrnrlon != None and urcrnrlat != None:
            m = Basemap(projection='cyl',llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat,
                    urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,ax=ax,resolution='l')
        else:
            #m = Basemap(projection='cyl',lon_0=0.,ax=ax,resolution='l')
            m = Basemap(projection='cyl',lon_0=lon_0,ax=ax,resolution='l')
            #m = Basemap(projection='cyl',ax=ax,resolution='l')

        x,y = m(gridLon,gridLat)

        # Plot the data
        im = m.imshow(gridData[:,:],axes=ax,interpolation='nearest',vmin=vmin,vmax=vmax)
        
        m.drawmapboundary(ax=ax,linewidth=0.01,fill_color='grey')
        m.drawcoastlines(ax=ax,linewidth=0.5)

        # draw parallels
        delat = 30.
        circles = np.arange(-90.,90.+delat,delat)
        m.drawparallels(circles,ax=ax,labelstyle="+/-",labels=[1,0,0,0])

        # draw meridians
        delon = 60.
        meridians = np.arange(-180,180,delon)
        m.drawmeridians(meridians,ax=ax,labelstyle="+/-",labels=[0,0,0,1])

        # add a colorbar axis
        cax_rect = [0.05 , 0.05, 0.9 , 0.06 ] # [left,bottom,width,height]
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
        levelStr = "_"+string.replace(string.split(plotTitle,'@')[1],')','') if "@" in plotTitle else ""
        levelStr = string.replace(levelStr," ","")
        LOG.debug("LevelStr = {}".format(levelStr))
        LOG.debug("pngDir = {}".format(self.pngDir))
        pngFile = "{}/{}{}_NCEP-ANC-Int_{}{}.png".format(self.pngDir,
                prefix,path.basename(blobName),
                self.ancDataSetTitle[ancType],levelStr)
        LOG.info("Writing file to {}".format(pngFile))
        canvas.print_figure(pngFile,dpi=dpi)

        del(m)
        return 0


    def plot_ncep_data(self,NCEP_ANC_xmlFile,NCEP_ANC_BlobFile,endian,level='1000mb'):
        '''
        Loop through the NCEP datasets that have been selected for plotting,
        and call the plotting routine for each one.
        '''

        NCEP_ancObj = adl_blob.map(NCEP_ANC_xmlFile,NCEP_ANC_BlobFile,
                endian=endian)

        # Check whether we are plotting a pressure level or something else...
        for dataset in self.plot_datasets :
            LOG.info("Ingesting dataset {}".format(dataset))

            NCEP_anc = getattr(NCEP_ancObj,dataset)

            try :
                if 'Layers' in dataset:
                    LOG.debug("Reading a layer dataset")
                    layer = self.NCEP_LAYER_LEVELS[level]
                    NCEP_anc = NCEP_anc[layer,:,:]
                else:
                    LOG.debug("Reading a level dataset")
                    NCEP_anc = NCEP_anc[:,:]

                LOG.debug("Inital shape(NCEP_anc) = {}".format(np.shape(NCEP_anc)))

                if endian == adl_blob.BIG_ENDIAN: 
                    NCEP_anc = NCEP_anc[::-1,:] # For ADL BE blobs only

                #NCEP_anc = np.roll(NCEP_anc,360) # Roll flattened array
                #NCEP_anc = np.roll(NCEP_anc,0,axis=0) # Roll over rows
                #NCEP_anc = np.roll(NCEP_anc,360,axis=1) # Roll over cols

                NCEP_anc_min = np.min(NCEP_anc)
                NCEP_anc_max = np.max(NCEP_anc)
                LOG.debug("min(NCEP_anc) = {}".format(np.min(NCEP_anc_min)))
                LOG.debug("max(NCEP_anc) = {}".format(np.max(NCEP_anc_max)))

                ### A default 0.5 degree grid...
                degInc = 0.5
                grids = np.mgrid[-90.:90.+degInc:degInc,-180.:180.:degInc] # original
                gridLat,gridLon = grids[0],grids[1]

                LOG.debug("Final shape(NCEP_anc) = {}".format(np.shape(NCEP_anc)))
                LOG.debug("Final shape(gridLon) = {}".format(np.shape(gridLon)))
                LOG.debug("Final shape(gridLat) = {}".format(np.shape(gridLat)))

                plotData = {}
                plotData['ancType'] = dataset
                plotData['plotTitle'] = "{}\n{}".format(path.basename(NCEP_ANC_BlobFile),
                        plotData['ancType'])
                plotData['cbTitle'] = self.plotDescr[dataset]
                plotData['plotLims'] = self.plotLims[dataset]
                plotData['NCEP_ANC_BlobFile'] = NCEP_ANC_BlobFile
                plotData['dpi'] = self.dpi
                plotData['prefix'] = self.outputFilePrefix

                self.plotAncData(gridLat,gridLon,NCEP_anc,plotData)

            except Exception, err :
                LOG.debug(traceback.format_exc())

        del(NCEP_ancObj)

        return 0


def _argparse():
    '''
    Method to encapsulate the option parsing and various setup tasks.
    '''

    import argparse

    ANCobj = ANCclass()

    prodChoices = ANCobj.data_names
    endian_choices = ['big','little']

    defaults = {
                'plotProduct' : None,
                'endianness' : 'little',
                'plotPass' : False,
                'plotMin'  : None,
                'plotMax'  : None,
                'latMin'  : None,
                'latMax'  : None,
                'lonMin'  : None,
                'lonMax'  : None,
                'dpi'      : 100,
                'mapAnn'   : None,
                'pngDir'   : None,
                'outputFilePrefix' : ""
                }

    description = '''Boilerplate code which shows how to use argparse, and tries to exercise
most of the input types.'''
    usage = "usage: %prog [mandatory args] [options]"
    version = __version__

    parser = argparse.ArgumentParser()

    # Mandatory arguments

    parser.add_argument('-b','--blob_file',
                      action='store',
                      dest='blob_file',
                      type=str,
                      required=True,
                      help="""The fully qualified path to the input blob file 
                      directory."""
                      )

    parser.add_argument('-x','--xml_file',
                      action='store',
                      dest='xml_file',
                      type=str,
                      required=True,
                      help="""The fully qualified path to the input xml file 
                      directory."""
                      )

    # Optional arguments 

    parser.add_argument('-p','--plotProduct',
                      action="store",
                      dest="plotProduct",
                      default=defaults["plotProduct"],
                      type=str,
                      choices=prodChoices,
                      help='''Which VIIRS ANC product to plot. If not specified,
                      all available products will be plotted. Possible options 
                      are...\n
                              {}.
                           '''.format(prodChoices.__str__()[1:-1])
                      )

    parser.add_argument('-e','--endianness',
                      action="store",
                      dest="endianness",
                      default=defaults["endianness"],
                      type=str,
                      choices=endian_choices,
                      help='''The endianness of the NCEP-ANC-Int blob file.
                      Possible options are...\n
                              {}.
                           '''.format(endian_choices.__str__()[1:-1])
                      )

    parser.add_argument('--plotMin',
                      action="store",
                      dest="plotMin",
                      default=defaults["plotMin"],
                      type=float,
                      help="Minimum value to plot.".format(defaults["plotMin"])
                      )

    parser.add_argument('--plotMax',
                      action="store",
                      dest="plotMax",
                      default=defaults["plotMax"],
                      type=float,
                      help="Maximum value to plot.".format(defaults["plotMax"])
                      )

    parser.add_argument('--latMin',
                      action="store",
                      dest="latMin",
                      default=defaults["latMin"],
                      type=float,
                      help="Minimum latitude of plot.".format(defaults["latMin"])
                      )

    parser.add_argument('--latMax',
                      action="store",
                      dest="latMax",
                      default=defaults["latMax"],
                      type=float,
                      help="Maximum latitude of plot.".format(defaults["latMax"])
                      )

    parser.add_argument('--lonMin',
                      action="store",
                      dest="lonMin",
                      default=defaults["lonMin"],
                      type=float,
                      help="Minimum longitude of plot.".format(defaults["lonMin"])
                      )

    parser.add_argument('--lonMax',
                      action="store",
                      dest="lonMax",
                      default=defaults["lonMax"],
                      type=float,
                      help="Maximum longitude of plot.".format(defaults["lonMax"])
                      )

    parser.add_argument('-d','--dpi',
                      action="store",
                      dest="dpi",
                      default=defaults["dpi"],
                      type=int,
                      help="""The resolution in dots per inch of the output png file. 
                      [default: {}]""".format(defaults["dpi"])
                      )

    parser.add_argument('--mapAnn',
                      action="store",
                      dest="mapAnn",
                      default=defaults["mapAnn"],
                      type=str,
                      help="""The map legend describing the dataset being shown.
                      [default: {}]""".format(defaults["mapAnn"])
                      )

    parser.add_argument('--pngDir',
                      action="store",
                      dest="pngDir",
                      default=defaults["pngDir"],
                      type=str,
                      help="""The directory where png files will be written. 
                      [default: {}]""".format(defaults["pngDir"])
                      )

    parser.add_argument('-o','--output_file_prefix',
                      action="store",
                      dest="outputFilePrefix",
                      default=defaults["outputFilePrefix"],
                      type=str,
                      help="""String to prefix to the automatically generated 
                      png names, which are of the form 
                      <N_Collection_Short_Name>_<N_Granule_ID>.png. 
                      [default: {}]""".format(defaults["outputFilePrefix"])
                      )

    parser.add_argument("-v", "--verbose",
                      dest='verbosity',
                      action="count", 
                      default=0,
                      help='each occurrence increases verbosity 1 level from ERROR: -v=WARNING -vv=INFO -vvv=DEBUG')


    args = parser.parse_args()


    # Set up the logging
    console_logFormat = '%(asctime)s : (%(levelname)s):%(filename)s:%(funcName)s:%(lineno)d:  %(message)s'
    dateFormat = '%Y-%m-%d %H:%M:%S'
    levels = [logging.ERROR, logging.WARN, logging.INFO, logging.DEBUG]
    logging.basicConfig(level = levels[args.verbosity], 
            format = console_logFormat, 
            datefmt = dateFormat)


    return args,ANCobj


###################################################
#                  Main Function                  #
###################################################


def main():
    '''
    The main method.
    '''

    options,ANCobj = _argparse()

    blob_file = options.blob_file
    xml_file  = options.xml_file
    plotProduct  = options.plotProduct
    endianness  = options.endianness
    plotMin = options.plotMin
    plotMax = options.plotMax
    latMin = options.latMin
    latMax = options.latMax
    lonMin = options.lonMin
    lonMax = options.lonMax
    dpi = options.dpi
    mapAnn = options.mapAnn
    pngDir = options.pngDir
    outputFilePrefix = options.outputFilePrefix


    blob_file = path.abspath(path.expanduser(blob_file))
    xml_file = path.abspath(path.expanduser(xml_file))

    pngDir = '.' if (pngDir == None) else pngDir
    ANCobj.pngDir = path.abspath(path.expanduser(pngDir))

    ANCobj.plot_datasets = ANCobj.data_names if (plotProduct==None) else [plotProduct]

    # If there is only one dataset, and we have specified the plot limits...
    if len(ANCobj.plot_datasets)==1 :
        dataset = ANCobj.plot_datasets[0]
        if plotMin != None :
            ANCobj.plotLims[dataset][0] = plotMin
        if plotMax != None :
            ANCobj.plotLims[dataset][1] = plotMax

    # Set the endianness
    if endianness == 'big':
        endian = adl_blob.BIG_ENDIAN
    else:
        endian = adl_blob.LITTLE_ENDIAN

    # Set the lat and lon plot limits
    ANCobj.latMin = latMin
    ANCobj.latMax = latMax
    ANCobj.lonMin = lonMin
    ANCobj.lonMax = lonMax

    ANCobj.dpi = dpi
    ANCobj.outputFilePrefix = outputFilePrefix

    LOG.info("blob_file = {}".format(blob_file))
    LOG.info("xml_file = {}".format(xml_file))
    LOG.info("plotProduct = {}".format(ANCobj.plot_datasets))
    LOG.info("endianness = {}".format(endianness))
    LOG.info("plotMin = {}".format(plotMin))
    LOG.info("plotMax = {}".format(plotMax))
    LOG.info("latMin = {}".format(latMin))
    LOG.info("latMax = {}".format(latMax))
    LOG.info("lonMin = {}".format(lonMin))
    LOG.info("lonMax = {}".format(lonMax))
    LOG.info("dpi = {}".format(dpi))
    LOG.info("mapAnn = {}".format(mapAnn))
    LOG.info("pngDir = {}".format(pngDir))
    LOG.info("outputFilePrefix = {}".format(outputFilePrefix))

    try :

        ANCobj.plot_ncep_data(xml_file,blob_file,endian)
        #pass

    except Exception, err :
        LOG.debug(traceback.format_exc())


    return 0


if __name__=='__main__':
    sys.exit(main())  
