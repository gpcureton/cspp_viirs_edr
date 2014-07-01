#!/usr/bin/env python
# encoding: utf-8
"""
CSPP_VIIRS_ANC_blobPlots.py

Purpose: Create swath projection quicklook PNGs from the NCEP granulated blob files.

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

    python CSPP_VIIRS_ANC_blobPlots.py -b '/path/to/blobs' -x '/path/to/xml'


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

import optparse as optparse

import adl_blob
import adl_blob2
from ViirsData import ViirsTrimTable
import viirs_edr_data

# every module should have a LOG object
# e.g. LOG.warning('my dog has fleas')
import logging
LOG = logging.getLogger(__file__)


### Moderate and Imager resolution trim table arrays. These are 
### bool arrays, and the trim pixels are set to True.
trimObj = ViirsTrimTable()
modTrimMask = trimObj.createModTrimArray(nscans=48,trimType=bool)


def get_blob_asc_dict(xmlDir,blobDir,collShortNames):
    blob_asc_Files = {}
    
    blobDir = path.expanduser(blobDir)

    if path.isdir(blobDir):
        pass
    else :
        blobDir = path.dirname(blobDir)

    for shortName in collShortNames :
        fileGlob = path.join(path.expanduser(blobDir),'*.%s'%(shortName))

        blobFiles = glob(fileGlob)
        if blobFiles != []:
            blobFiles.sort()
            blobDict = {}
            LOG.info('{} --> '.format(shortName))
            for files in blobFiles :
                blobFile = path.basename(files)
                ascFile = string.replace(blobFile,shortName,'asc')
                ascFullPath = path.join(blobDir,'%s'%(ascFile))
                ascFileObj = open(ascFullPath,"r")
                for line in ascFileObj:
                    if re.search("N_Granule_ID", line):
                        granID = string.split(line,'"')[3]
                ascFileObj.close()
                blobDict[granID] = [blobFile,ascFile]
            blob_asc_Files[shortName] = blobDict

    return blob_asc_Files

#def NISE_globalPlot(NISE_fileName,pngDir=None):

    #from HDF4File import HDF4File

    #if pngDir is None :
        #pngDir = path.abspath(path.curdir)

    #latStart,latEnd = 2*3000,2*4500
    #lonStart,lonEnd = 2*3000,2*5000
    #print "NISE_latMinIdx = %d" % (latEnd)
    #print "NISE_latMaxIdx = %d" % (latStart)
    #print "NISE_lonMinIdx = %d" % (lonStart)
    #print "NISE_lonMaxIdx = %d" % (lonEnd)

    #try :
        #fileObj = HDF4File(NISE_fileName)
    #except Exception, err :
        #print "%s" % (err)
        #print "Problem opening NISE file (%s), aborting." % (NISE_fileName)
        #sys.exit(1)

    #try :

        #northDsetName = "Northern Hemisphere/Data Fields/Extent"
        #southDsetName = "Southern Hemisphere/Data Fields/Extent"

        #print "Retrieving NISE HDF4 path '%s'" % (northDsetName)
        #nHemi = fileObj.get_dataset(northDsetName)

        #print "Retrieving NISE HDF4 path '%s'" % (southDsetName)
        #sHemi = fileObj.get_dataset(southDsetName)

    #except Exception, err :

        #print "EXCEPTION: %s" % (err)
        #sys.exit(1)

    #for data,title in zip([nHemi,sHemi],['Northern','Sourthern']):

        #plotTitle =  "NISE %s Hemi : %s" %(title,path.basename(NISE_fileName))
        #cbTitle   =  "Snow / Ice"
        ##vmin,vmax =  0,1
        #vmin,vmax =  None,None

        ## Create figure with default size, and create canvas to draw on
        #scale=1.5
        #fig = Figure(figsize=(scale*8,scale*8))
        #canvas = FigureCanvas(fig)

        ## Create main axes instance, leaving room for colorbar at bottom,
        ## and also get the Bbox of the axes instance
        #ax_rect = [0.05, 0.18, 0.9, 0.75  ] # [left,bottom,width,height]
        #ax = fig.add_axes(ax_rect)

        ## Granule axis title
        #ax_title = ppl.setp(ax,title=plotTitle)
        #ppl.setp(ax_title,fontsize=12)
        #ppl.setp(ax_title,family="sans-serif")

        ## Plot the data
        #im = ax.imshow(data,axes=ax,interpolation='nearest',vmin=vmin,vmax=vmax)
        ##im = ax.imshow(data,axes=ax,interpolation='nearest')
        
        ## add a colorbar axis
        #cax_rect = [0.05 , 0.05, 0.9 , 0.06 ] # [left,bottom,width,height]
        #cax = fig.add_axes(cax_rect,frameon=False) # setup colorbar axes

        ## Plot the colorbar.
        #cb = fig.colorbar(im, cax=cax, orientation='horizontal')
        #ppl.setp(cax.get_xticklabels(),fontsize=9)

        ## Colourbar title
        #cax_title = ppl.setp(cax,title=cbTitle)
        #ppl.setp(cax_title,fontsize=9)

        ## Redraw the figure
        #canvas.draw()

        ## save image 
        #pngFile = "%s/%s_%s.png" % (pngDir,path.basename(NISE_fileName),title)
        #print "Writing file to ",pngFile
        #canvas.print_figure(pngFile,dpi=200)


class ANCclass():

    def __init__(self):


        self.collShortNames = [
                               'VIIRS-ANC-Preci-Wtr-Mod-Gran',
                               'VIIRS-ANC-Temp-Surf2M-Mod-Gran',
                               'VIIRS-ANC-Wind-Speed-Mod-Gran',
                               'VIIRS-ANC-Wind-Direction-Mod-Gran',
                               'VIIRS-ANC-Surf-Ht-Mod-Gran',
                               'VIIRS-ANC-Press-Surf-Mod-Gran',
                               'VIIRS-ANC-Tot-Col-Mod-Gran',
                               'VIIRS-ANC-Optical-Depth-Mod-Gran',
                               'VIIRS-ANC-Geopot-Ht-Lev-Mod-Gran',
                               'VIIRS-ANC-Sp-Humd-Surf-Mod-Gran',
                               'VIIRS-ANC-Temp-Skin-Mod-Gran'
                              ]


        self.xmlName = {}
        self.xmlName['VIIRS-ANC-Preci-Wtr-Mod-Gran'] = 'VIIRS_ANC_PRECI_WTR_MOD_GRAN.xml'
        self.xmlName['VIIRS-ANC-Temp-Surf2M-Mod-Gran'] = 'VIIRS_ANC_TEMP_SURF2M_MOD_GRAN.xml'
        self.xmlName['VIIRS-ANC-Wind-Speed-Mod-Gran'] = 'VIIRS_ANC_WIND_SPEED_MOD_GRAN.xml'
        self.xmlName['VIIRS-ANC-Wind-Direction-Mod-Gran'] = 'VIIRS_ANC_WIND_DIRECTION_MOD_GRAN.xml'
        self.xmlName['VIIRS-ANC-Surf-Ht-Mod-Gran'] = 'VIIRS_ANC_SURF_HT_MOD_GRAN.xml'
        self.xmlName['VIIRS-ANC-Press-Surf-Mod-Gran'] = 'VIIRS_ANC_PRESS_SURF_MOD_GRAN.xml'
        self.xmlName['VIIRS-ANC-Tot-Col-Mod-Gran'] = 'VIIRS_ANC_TOT_COL_MOD_GRAN.xml'
        self.xmlName['VIIRS-ANC-Optical-Depth-Mod-Gran'] = 'VIIRS_ANC_OPTICAL_DEPTH_MOD_GRAN.xml'
        self.xmlName['VIIRS-ANC-Geopot-Ht-Lev-Mod-Gran'] = 'VIIRS_ANC_GEOPOT_HT_LEV_MOD_GRAN.xml'
        self.xmlName['VIIRS-ANC-Sp-Humd-Surf-Mod-Gran'] = 'VIIRS_ANC_SP_HUMD_SURF_MOD_GRAN.xml'
        self.xmlName['VIIRS-ANC-Temp-Skin-Mod-Gran'] = 'VIIRS_ANC_TEMP_SKIN_MOD_GRAN.xml'

        self.plotDescr = {}
        self.plotDescr['VIIRS-ANC-Preci-Wtr-Mod-Gran']   = r'Total Precipitable Water (cm, $\times$10 kg m$^{-2}$)'
        #self.plotDescr['VIIRS-ANC-Preci-Wtr-Mod-Gran']   = r'Total Precipitable Water (kg m$^{-2}$)'
        self.plotDescr['VIIRS-ANC-Temp-Surf2M-Mod-Gran'] = r'Air Temperature @ 2m (K)'
        self.plotDescr['VIIRS-ANC-Wind-Speed-Mod-Gran']  = r'Surface Wind velocity (ms$^{-1}$)'
        self.plotDescr['VIIRS-ANC-Wind-Direction-Mod-Gran'] = r'Surface Wind Direction (degrees)'
        self.plotDescr['VIIRS-ANC-Surf-Ht-Mod-Gran']     = r'GMTCO Surface Height (m)'
        self.plotDescr['VIIRS-ANC-Press-Surf-Mod-Gran']  = r'Surface Pressure (mb, $\times$10$^{2}$ Pa)'
        self.plotDescr['VIIRS-ANC-Tot-Col-Mod-Gran']     = r'Total Column Ozone (Atm cm, $\times$10$^{3}$DU)'
        self.plotDescr['VIIRS-ANC-Optical-Depth-Mod-Gran'] = 'NAAPS Total Column Aerosol Optical Depth'
        self.plotDescr['VIIRS-ANC-Geopot-Ht-Lev-Mod-Gran'] = 'Surface Geopotential Height (gpm)'
        self.plotDescr['VIIRS-ANC-Sp-Humd-Surf-Mod-Gran'] = 'Surface Specific Humidity'
        self.plotDescr['VIIRS-ANC-Temp-Skin-Mod-Gran'] = 'Sea Surface Skin Temperature (K)'

        self.plotLims = {}
        self.plotLims['VIIRS-ANC-Preci-Wtr-Mod-Gran']    = [0.,10.]
        self.plotLims['VIIRS-ANC-Temp-Surf2M-Mod-Gran']  = [275.,315.]
        self.plotLims['VIIRS-ANC-Wind-Speed-Mod-Gran']   = [0.,20.]
        self.plotLims['VIIRS-ANC-Wind-Direction-Mod-Gran'] =[0.,360.]
        self.plotLims['VIIRS-ANC-Surf-Ht-Mod-Gran']      = [None, None]
        self.plotLims['VIIRS-ANC-Press-Surf-Mod-Gran']  = [800.,1013.]
        self.plotLims['VIIRS-ANC-Tot-Col-Mod-Gran']     = [0.2,0.40]
        self.plotLims['VIIRS-ANC-Optical-Depth-Mod-Gran'] = [None, None]
        self.plotLims['VIIRS-ANC-Geopot-Ht-Lev-Mod-Gran'] = [None, None]
        self.plotLims['VIIRS-ANC-Sp-Humd-Surf-Mod-Gran'] = [0.002,0.04]
        self.plotLims['VIIRS-ANC-Temp-Skin-Mod-Gran'] = [275.,315.]

        self.dataName = {}
        self.dataName['VIIRS-ANC-Preci-Wtr-Mod-Gran'] = 'data'
        self.dataName['VIIRS-ANC-Temp-Surf2M-Mod-Gran'] = 'data'
        self.dataName['VIIRS-ANC-Wind-Speed-Mod-Gran'] = 'data'
        self.dataName['VIIRS-ANC-Wind-Direction-Mod-Gran'] = 'data'
        self.dataName['VIIRS-ANC-Surf-Ht-Mod-Gran'] = 'data'
        self.dataName['VIIRS-ANC-Press-Surf-Mod-Gran'] = 'data'
        self.dataName['VIIRS-ANC-Tot-Col-Mod-Gran'] = 'data'
        self.dataName['VIIRS-ANC-Optical-Depth-Mod-Gran'] = 'faot550'
        self.dataName['VIIRS-ANC-Geopot-Ht-Lev-Mod-Gran'] = 'data'
        self.dataName['VIIRS-ANC-Sp-Humd-Surf-Mod-Gran'] = 'data'
        self.dataName['VIIRS-ANC-Temp-Skin-Mod-Gran'] = 'data'


    def set_blob_dict(self,xmlDir,blobPath,shortName):

        self.xmlDir = xmlDir
        self.blobPath = path.expanduser(blobPath)
        self.blob_dict = get_blob_asc_dict(self.xmlDir,self.blobPath,shortName)


    def plot_ANC_pass(self,pngDir=None,endian=adl_blob.LITTLE_ENDIAN):

        if pngDir is None :
            pngDir = path.abspath(path.curdir)

        dpi = self.dpi

        xmlDir = self.xmlDir
        blobPath = self.blobPath

        xmlName = self.xmlName
        plotDescr = self.plotDescr
        plotLims = self.plotLims

        blobDict = self.blob_dict
        collShortNames = blobDict.keys()

        for shortName in collShortNames :

            dataName = self.dataName[shortName]
            LOG.info("dataName = {}".format(dataName))

            xmlFile = path.join(xmlDir,xmlName[shortName])

            granID_list = blobDict[shortName].keys()
            granID_list.sort()

            for granID in granID_list :

                LOG.info('{} --> {}'.format(shortName, granID))
                
                blobFile = blobDict[shortName][granID][0]
                blobFile = path.join(blobPath,'%s'%(blobFile))
                blobObj = adl_blob.map(xmlFile,blobFile,endian=endian)
                blobArrObj = blobObj.as_arrays()

                dataGranule = getattr(blobArrObj,dataName)
                
                LOG.info("{} is of kind {}".format(shortName,dataGranule.dtype.kind))
                if (dataGranule.dtype.kind == 'O') :
                    dataGranule = data.astype(np.float)

                # Concatenate the granules.
                try :
                    data = np.vstack((data,dataGranule))
                    LOG.info("data shape = {}".format(data.shape))
                except :
                    data = dataGranule[:,:]
                    LOG.info("data shape = {}".format(data.shape))

            LOG.info("Final data shape = {}".format(data.shape))

            # Assuming this is a descending granule, flip it...
            data = data[::-1,::-1]

            # What value are the bowtie deletion pixels
            ongroundPixelTrimValue = trimObj.sdrTypeFill['ONGROUND_PT_FILL'][data.dtype.name]
            print "Onground Pixel Trim value is {}".format(ongroundPixelTrimValue)

            # Create onboard and onground pixel trim mask arrays, for the total number of
            # scans in the pass...
            numGranules = len(blobDict[shortName].keys())
            numScans = numGranules * 48
            ongroundTrimMask = trimObj.createModTrimArray(nscans=numScans,trimType=bool)

            # Apply the On-ground pixel trim
            data = ma.array(data,mask=ongroundTrimMask,fill_value=ongroundPixelTrimValue)
            data = data.filled() # Substitute for the masked values with onboardPixelTrimValue

            vmin,vmax = plotLims[shortName]
            
            plotTitle = '{} : {}'.format(shortName,granID)
            cbTitle = plotDescr[shortName]

            # Create figure with default size, and create canvas to draw on
            passRows = float(data.shape[0])
            passCols = float(data.shape[1])
            aspect = passRows/passCols
            scale = 5.
            fig = Figure(figsize=(scale*1.,scale*aspect))
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
            im = ax.imshow(ma.masked_less(data,-800.),axes=ax,interpolation='nearest',vmin=vmin,vmax=vmax)
            #im = ax.imshow(ma.masked_less(data,-800.),axes=ax,interpolation='nearest')
            
            ppl.setp(ax.get_xticklabels(), visible=False)
            ppl.setp(ax.get_yticklabels(), visible=False)
            ppl.setp(ax.get_xticklines(),visible=False)
            ppl.setp(ax.get_yticklines(),visible=False)

            # add a colorbar axis
            cax_rect = [0.05 , 0.05, 0.9 , 0.10 ] # [left,bottom,width,height]
            cax = fig.add_axes(cax_rect,frameon=False) # setup colorbar axes

            # Plot the colorbar.
            cb = fig.colorbar(im, cax=cax, orientation='horizontal')
            ppl.setp(cax.get_xticklabels(),fontsize=9)

            # Colourbar title
            cax_title = ppl.setp(cax,title=cbTitle)
            ppl.setp(cax_title,fontsize=10)

            # Redraw the figure
            canvas.draw()

            # Save the figure to a png file...
            pngFile = path.join(pngDir,'%s_%s.png' % (shortName,granID))
            LOG.info("Writing to {} ...".format(pngFile))
            canvas.print_figure(pngFile,dpi=dpi)

            ppl.close('all')
            del(data)


    def plot_ANC_granules(self,pngDir=None,endian=adl_blob.LITTLE_ENDIAN):

        if pngDir is None :
            pngDir = path.abspath(path.curdir)

        xmlDir = self.xmlDir
        blobPath = self.blobPath

        xmlName = self.xmlName
        plotDescr = self.plotDescr
        plotLims = self.plotLims

        blobDict = self.blob_dict
        #collShortNames = blobDict.keys()
        collShortNames = self.collShortNames

        for shortName in collShortNames :

            dataName = self.dataName[shortName]

            for granID in  blobDict[shortName].keys() :
                print '%s --> %s ' % (shortName, granID)
                xmlFile = path.join(xmlDir,xmlName[shortName])
                blobFile = blobDict[shortName][granID][0]
                blobFile = path.join(blobPath,'%s'%(blobFile))
                blobObj = adl_blob.map(xmlFile,blobFile,endian=endian)
                print "dataName = %s" % (dataName)
                blobArrObj = blobObj.as_arrays()
                data = getattr(blobArrObj,dataName)
                print "%s is of kind %r" % (shortName,data.dtype.kind)
                if (data.dtype.kind == 'O') :
                    data = data.astype(np.float)

                vmin,vmax = plotLims[shortName]
                plotTitle = '%s : %s' % (shortName,granID)
                cbTitle = plotDescr[shortName]

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
                im = ax.imshow(ma.masked_less(data,-800.),axes=ax,interpolation='nearest',vmin=vmin,vmax=vmax)
                #im = ax.imshow(ma.masked_less(data,-800.),axes=ax,interpolation='nearest')
                
                ppl.setp(ax.get_xticklabels(), visible=False)
                ppl.setp(ax.get_yticklabels(), visible=False)
                ppl.setp(ax.get_xticklines(),visible=False)
                ppl.setp(ax.get_yticklines(),visible=False)

                # add a colorbar axis
                cax_rect = [0.05 , 0.05, 0.9 , 0.10 ] # [left,bottom,width,height]
                cax = fig.add_axes(cax_rect,frameon=False) # setup colorbar axes

                # Plot the colorbar.
                cb = fig.colorbar(im, cax=cax, orientation='horizontal')
                ppl.setp(cax.get_xticklabels(),fontsize=9)

                # Colourbar title
                cax_title = ppl.setp(cax,title=cbTitle)
                ppl.setp(cax_title,fontsize=10)

                # Redraw the figure
                canvas.draw()

                # Save the figure to a png file...
                pngFile = path.join(pngDir,'%s_%s.png' % (shortName,granID))
                canvas.print_figure(pngFile,dpi=100)

                ppl.close('all')


def _argparse():
    '''
    Method to encapsulate the option parsing and various setup tasks.
    '''

    import argparse

    ANCobj = ANCclass()

    prodChoices = ANCobj.collShortNames

    defaults = {
                'plotProduct' : None,
                'plotPass' : False,
                'plotMin'  : None,
                'plotMax'  : None,
                'dpi'      : 200,
                'mapAnn'   : None,
                'pngDir'   : None,
                'outputFilePrefix' : None
                }

    description = '''Boilerplate code which shows how to use argparse, and tries to exercise
most of the input types.'''
    usage = "usage: %prog [mandatory args] [options]"
    version = __version__

    parser = argparse.ArgumentParser()

    # Mandatory arguments

    parser.add_argument('-b','--blob_dir',
                      action='store',
                      dest='blob_dir',
                      type=str,
                      required=True,
                      help="""The fully qualified path to the input blob file 
                      directory."""
                      )

    parser.add_argument('-x','--xml_dir',
                      action='store',
                      dest='xml_dir',
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

    parser.add_argument('--plotPass',
                      action="store_true",
                      dest="plotPass",
                      default=defaults["plotPass"],
                      help="""Concatenate the product blob files into a pass
                      [default: {}]""".format(defaults["plotPass"])
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

    blob_dir = options.blob_dir
    xml_dir  = options.xml_dir
    plotProduct  = options.plotProduct
    plotPass = options.plotPass
    plotMin = options.plotMin
    plotMax = options.plotMax
    dpi = options.dpi
    mapAnn = options.mapAnn
    pngDir = options.pngDir
    outputFilePrefix = options.outputFilePrefix


    blob_dir = path.abspath(path.expanduser(blob_dir))
    xml_dir = path.abspath(path.expanduser(xml_dir))

    pngDir = '.' if (pngDir == None) else pngDir
    pngDir = path.abspath(path.expanduser(pngDir))

    plotProduct = ANCobj.collShortNames if (plotProduct==None) else [plotProduct]

    # If there is only one dataset, and we have specified the plot limits...
    if len(plotProduct)==1 :
        dataset = plotProduct[0]
        if plotMin != None :
            ANCobj.plotLims[dataset][0] = plotMin
        if plotMax != None :
            ANCobj.plotLims[dataset][1] = plotMax

    ANCobj.dpi = dpi

    LOG.info("blob_dir = {}".format(blob_dir))
    LOG.info("xml_dir = {}".format(xml_dir))
    LOG.info("plotProduct = {}".format(plotProduct))
    LOG.info("plotPass = {}".format(plotPass))
    LOG.info("plotMin = {}".format(plotMin))
    LOG.info("plotMax = {}".format(plotMax))
    LOG.info("dpi = {}".format(dpi))
    LOG.info("mapAnn = {}".format(mapAnn))
    LOG.info("pngDir = {}".format(pngDir))
    LOG.info("outputFilePrefix = {}".format(outputFilePrefix))

    try :
        ANCobj.set_blob_dict(xml_dir,blob_dir,plotProduct)

        if plotPass :
            LOG.info("Plotting a pass")
            ANCobj.plot_ANC_pass(pngDir=pngDir,endian=adl_blob.LITTLE_ENDIAN)
        else:
            LOG.info("Plotting a single granule")

    except Exception, err:
        traceback.print_exc(file=sys.stdout)



    return 0


if __name__=='__main__':
    sys.exit(main())  
