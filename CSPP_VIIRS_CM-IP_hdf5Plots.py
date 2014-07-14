#!/usr/bin/env python
# encoding: utf-8
"""
CSPP_VIIRS_CM-IP_hdf5Plots.py

Purpose: Create swath projection quicklook PNGs from the VIIRS SST EDR HDF5 files.
         Images can be created for the EDR product or the associated quality flags.

Minimum commandline...

export CSPP_EDR_HOME=$(readlink -f /path/to/EDR)
. ${CSPP_EDR_HOME}/cspp_edr_runtime.sh

python CSPP_VIIRS_CM-IP_hdf5Plots.py -i '/path/to/files/IICMO*.h5'

  or

python CSPP_VIIRS_CM-IP_hdf5Plots.py --input_files=/path/to/files/IICMO*.h5


Created by Geoff Cureton on 2013-06-04.
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

#############

import os, sys
from os import path, uname, mkdir
from glob import glob
import string, logging, traceback
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

import optparse as optparse

from ViirsData import ViirsTrimTable
import viirs_edr_data

import tables as pytables
from tables import exceptions as pyEx


# every module should have a LOG object
# e.g. LOG.warning('my dog has fleas')
import logging
LOG = logging.getLogger(__file__)

dpi=200

### Moderate and Imager resolution trim table arrays. These are 
### bool arrays, and the trim pixels are set to True.
trimObj = ViirsTrimTable()
modTrimMask = trimObj.createModTrimArray(nscans=48,trimType=bool)


def get_hdf5_dict(hdf5Path,filePrefix):
    shortNameDict = {}
    
    hdf5Path = path.abspath(path.expanduser(hdf5Path))
    LOG.info("hdf5Path = %s" % (hdf5Path))

    hdf5Dir = path.dirname(hdf5Path)
    hdf5Glob = path.basename(hdf5Path)

    LOG.info("hdf5Dir = %s" % (hdf5Dir))
    LOG.info("hdf5Glob = %s" % (hdf5Glob))

    if (hdf5Glob == '' or hdf5Glob == '*'):
        LOG.info("prefix = %s" % (filePrefix))
        hdf5Glob = path.join(hdf5Dir,'%s_*.h5'%(filePrefix))
    else :
        hdf5Glob = path.join(hdf5Dir,'%s'%(hdf5Glob))

    LOG.info("Final hdf5Glob = %s" % (hdf5Glob))

    hdf5Files = glob(hdf5Glob)
    if hdf5Files != []:
        hdf5Files.sort()
        granIdDict = {}
        for files in hdf5Files :

            # Open the hdf5 file
            fileObj = pytables.openFile(files)
            
            # Get a "pointer" to the granules attribute group.
            VIIRS_CM_IP_Gran_0 = fileObj.getNode('/Data_Products/VIIRS-CM-IP/VIIRS-CM-IP_Gran_0')
            
            # Retrieve a few attributes
            granID =  getattr(VIIRS_CM_IP_Gran_0.attrs,'N_Granule_ID')[0][0]
            LOG.info('N_Granule_ID = %s' % (granID))
            
            dayNightFlag =  getattr(VIIRS_CM_IP_Gran_0.attrs,'N_Day_Night_Flag')[0][0]
            LOG.info('N_Day_Night_Flag = %s' % (dayNightFlag))
            
            shortName = fileObj.getNodeAttr('/Data_Products/VIIRS-CM-IP','N_Collection_Short_Name')[0][0]
            LOG.info('N_Collection_Short_Name = %s' % (shortName))
            
            # Strip the path from the filename
            hdf5File = path.basename(files)

            # Add the granule information to the dictionary, keyed with the granule ID...
            granIdDict[granID] = [hdf5File,fileObj]

        shortNameDict[shortName] = granIdDict

    return shortNameDict


class VCMclass():

    def __init__(self,hdf5Dir):

        self.hdf5Dir = hdf5Dir

        self.collShortNames = [
                               'VIIRS-CM-IP'
                              ]

        self.plotDescr = {}
        self.plotDescr['VIIRS-CM-IP'] = ['Cloud Mask','Cloud Phase']

        self.plotLims = {}
        self.plotLims['VIIRS-CM-IP'] = [0, 3]

        self.dataName = {}
        self.dataName['VIIRS-CM-IP'] = ['/All_Data/VIIRS-CM-IP_All/QF1_VIIRSCMIP', \
                                        '/All_Data/VIIRS-CM-IP_All/QF6_VIIRSCMIP']

        self.byteIdx = {}
        self.byteIdx['VIIRS-CM-IP'] = [0,5]

        self.byteDsetIdx = {}
        self.byteDsetIdx['VIIRS-CM-IP'] = [1,0]

        self.hdf5_dict = get_hdf5_dict(hdf5Dir,'IICMO')


    def plot_VCM_granules(self,plotProd='IP',pngDir=None,pngPrefix=None,annotation='',dpi=300):

        if pngDir is None :
            pngDir = path.abspath(path.curdir)

        plotDescr = self.plotDescr
        plotLims = self.plotLims

        hdf5_dict = self.hdf5_dict
        collShortNames = hdf5_dict.keys()

        LOG.info('collShortNames = %r' % (collShortNames))

        CMD = viirs_edr_data.CloudMaskData

        for shortName in collShortNames :

            LOG.info('shortName = %s' % (shortName))

            if (plotProd == 'IP'):

                dataNames = self.dataName[shortName]
                byteIdx = self.byteIdx[shortName]
                byteDsetIdx = self.byteDsetIdx[shortName]
                plotDescrs = plotDescr[shortName]
                prodNames = ['CMask','CPhase']

            elif (plotProd == 'CMask'):

                dataNames = [self.dataName[shortName][0]]
                byteIdx = [self.byteIdx[shortName][0]]
                byteDsetIdx = [self.byteDsetIdx[shortName][0]]
                plotDescrs = [plotDescr[shortName][0]]
                prodNames = ['CMask']

            elif (plotProd == 'CPhase'):

                dataNames = [self.dataName[shortName][1]]
                byteIdx = [self.byteIdx[shortName][1]]
                byteDsetIdx = [self.byteDsetIdx[shortName][1]]
                plotDescrs = [plotDescr[shortName][1]]
                prodNames = ['CPhase']

            granID_list =  hdf5_dict[shortName].keys()
            granID_list.sort()

            for granID in granID_list :

                LOG.info('%s --> %s ' % (shortName, granID))
                
                hdf5Obj = hdf5_dict[shortName][granID][1]
                
                VIIRS_CM_IP_Gran_0 = hdf5Obj.getNode('/Data_Products/VIIRS-CM-IP/VIIRS-CM-IP_Gran_0')
                dayNightFlag =  getattr(VIIRS_CM_IP_Gran_0.attrs,'N_Day_Night_Flag')[0][0]
                orbitNumber =  getattr(VIIRS_CM_IP_Gran_0.attrs,'N_Beginning_Orbit_Number')[0][0]
                LOG.info('N_Day_Night_Flag = %s' % (dayNightFlag))
                LOG.info('N_Beginning_Orbit_Number = %s' % (orbitNumber))
                orient = -1 if dayNightFlag == 'Day' else 1

                for dataName,byte,dSet,plotDescr,prodName in zip(dataNames,byteIdx,byteDsetIdx,plotDescrs,prodNames):

                    data = hdf5Obj.getNode(dataName)[:,:]
                    
                    pixelTrimValue = trimObj.sdrTypeFill['ONBOARD_PT_FILL'][data.dtype.name]
                    LOG.info("pixelTrimValue is %r" % (pixelTrimValue))

                    # Bit mask and shift to get the Cloud Mask
                    bitMask  = CMD.ViirsCMbitMasks[byte][dSet]
                    bitShift = CMD.ViirsCMbitShift[byte][dSet]

                    vmin,vmax = CMD.ViirsCMvalues[byte][dSet][0], CMD.ViirsCMvalues[byte][dSet][-1]
                    data = np.bitwise_and(data ,bitMask) >> bitShift

                    # Apply the moderate pixel trim, so that we can properly mask them out at plot time.
                    data = ma.array(data,mask=modTrimMask,fill_value=trimObj.sdrTypeFill['ONBOARD_PT_FILL'][data.dtype.name])
                    data = data.filled()

                    cmap = ListedColormap(CMD.ViirsCMfillColours[byte][dSet])

                    numCats = np.array(CMD.ViirsCMfillColours[byte][dSet]).size
                    numBounds = numCats + 1

                    tickPos = np.arange(float(numBounds))/float(numCats)
                    tickPos = tickPos[0 :-1] + tickPos[1]/2.

                    LOG.info("numCats = {}".format(numCats))
                    LOG.info("tickPos = {}".format(tickPos))

                    plotTitle = '%s : %s %s' % (shortName,granID,annotation)
                    cbTitle = plotDescr

                    # Create figure with default size, and create canvas to draw on
                    scale=1.5
                    fig = Figure(figsize=(scale*8,scale*3))
                    canvas = FigureCanvas(fig)

                    # Create main axes instance, leaving room for colorbar at bottom,
                    # and also get the Bbox of the axes instance
                    ax_rect = [0.05, 0.18, 0.9, 0.75  ] # [left,bottom,width,height]
                    ax = fig.add_axes(ax_rect,axis_bgcolor='lightgray')

                    # Granule axis title
                    ax_title = ppl.setp(ax,title=plotTitle)
                    ppl.setp(ax_title,fontsize=12)
                    ppl.setp(ax_title,family="sans-serif")

                    # Plot the data
                    LOG.info("%s is of kind %r" % (shortName,data.dtype.kind))
                    #LOG.info(data)
                    if (data.dtype.kind =='i' or data.dtype.kind =='u'):
                        data = ma.masked_greater(data,200)
                    else:
                        data = ma.masked_less(data,-800.)

                    im = ax.imshow(ma.masked_less(data[::orient,::orient],0),interpolation='nearest',vmin=vmin,vmax=vmax,cmap=cmap)
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
                    ppl.setp(cax.get_xticklines(),visible=False)

                    # Set the colourbar tick locations and ticklabels
                    #ppl.setp(cb.ax,xticks=CMD.ViirsCMTickPos) # In colorbar axis coords (0..1)
                    cb.set_ticks(vmax*tickPos) # In data coords (0..3)
                    ppl.setp(cb.ax,xticklabels=CMD.ViirsCMtickNames[byte][dSet])

                    # Colourbar title
                    cax_title = ppl.setp(cax,title=cbTitle)
                    ppl.setp(cax_title,fontsize=10)

                    # Turn off the tickmarks on the colourbar
                    #ppl.setp(cb.ax.get_xticklines(),visible=False)
                    #ppl.setp(cb.ax.get_xticklabels(),fontsize=9)

                    # Redraw the figure
                    canvas.draw()

                    # Save the figure to a png file...
                    pngFile = path.join(pngDir,'%s%s_%s_%s.png' % (pngPrefix,shortName,granID,prodName))
                    canvas.print_figure(pngFile,dpi=dpi)
                    LOG.info("Writing to %s..." % (pngFile))

                    ppl.close('all')

                hdf5Obj.close()


    def plot_VCM_pass(self,plotProd='IP',pngDir=None,pngPrefix=None,
            annotation='',dpi=300):

        if pngDir is None :
            pngDir = path.abspath(path.curdir)

        plotDescr = self.plotDescr
        plotLims = self.plotLims

        hdf5_dict = self.hdf5_dict
        collShortNames = hdf5_dict.keys()

        LOG.info('collShortNames = %r' % (collShortNames))

        CMD = viirs_edr_data.CloudMaskData

        for shortName in collShortNames :

            LOG.info('shortName = %s' % (shortName))

            if (plotProd == 'IP'):

                dataNames = self.dataName[shortName]
                byteIdx = self.byteIdx[shortName]
                byteDsetIdx = self.byteDsetIdx[shortName]
                plotDescrs = plotDescr[shortName]
                prodNames = ['CMask','CPhase']

            elif (plotProd == 'CMask'):

                dataNames = [self.dataName[shortName][0]]
                byteIdx = [self.byteIdx[shortName][0]]
                byteDsetIdx = [self.byteDsetIdx[shortName][0]]
                plotDescrs = [plotDescr[shortName][0]]
                prodNames = ['CMask']

            elif (plotProd == 'CPhase'):

                dataNames = [self.dataName[shortName][1]]
                byteIdx = [self.byteIdx[shortName][1]]
                byteDsetIdx = [self.byteDsetIdx[shortName][1]]
                plotDescrs = [plotDescr[shortName][1]]
                prodNames = ['CPhase']

            granID_list =  hdf5_dict[shortName].keys()
            granID_list.sort()

            for dataName,byte,dSet,plotDescr,prodName in zip(dataNames,byteIdx,byteDsetIdx,plotDescrs,prodNames):

                # Read in the data from the granules and concatenate
                for granID in granID_list :

                    LOG.info('%s --> %s ' % (shortName, granID))
                    
                    hdf5Obj = hdf5_dict[shortName][granID][1]
                    
                    VIIRS_CM_IP_Gran_0 = hdf5Obj.getNode('/Data_Products/VIIRS-CM-IP/VIIRS-CM-IP_Gran_0')
                    dayNightFlag =  getattr(VIIRS_CM_IP_Gran_0.attrs,'N_Day_Night_Flag')[0][0]
                    orbitNumber =  getattr(VIIRS_CM_IP_Gran_0.attrs,'N_Beginning_Orbit_Number')[0][0]
                    LOG.info('N_Day_Night_Flag = %s' % (dayNightFlag))
                    LOG.info('N_Beginning_Orbit_Number = %s' % (orbitNumber))
                    orient = -1 if dayNightFlag == 'Day' else 1

                    dataGranule = hdf5Obj.getNode(dataName)[:,:]

                    # Concatenate the granules.
                    try :
                        data = np.vstack((data,dataGranule))
                        LOG.info("data shape = {}".format(data.shape))
                    except NameError :
                        data = dataGranule[:,:]
                        #data = hdf5Obj.getNode(dataName)[:,:]
                        LOG.info("data shape = {}".format(data.shape))

                LOG.info("Final data shape = {}".format(data.shape))

                # What value are the bowtie deletion pixels
                onboardPixelTrimValue = trimObj.sdrTypeFill['ONBOARD_PT_FILL'][data.dtype.name]
                LOG.info("Onboard Pixel Trim value is {}".format(onboardPixelTrimValue))
                ongroundPixelTrimValue = trimObj.sdrTypeFill['ONGROUND_PT_FILL'][data.dtype.name]
                LOG.info("Onground Pixel Trim value is {}".format(ongroundPixelTrimValue))

                # Create onboard and onground pixel trim mask arrays, for the total number of
                # scans in the pass...
                numGranules = len(granID_list)
                numScans = numGranules * 48
                onboardTrimMask = trimObj.createOnboardModTrimArray(nscans=numScans,trimType=bool)
                ongroundTrimMask = trimObj.createOngroundModTrimArray(nscans=numScans,trimType=bool)

                # Bit mask and shift to get the Cloud Mask
                bitMask  = CMD.ViirsCMbitMasks[byte][dSet]
                bitShift = CMD.ViirsCMbitShift[byte][dSet]

                vmin,vmax = CMD.ViirsCMvalues[byte][dSet][0], CMD.ViirsCMvalues[byte][dSet][-1]
                data = np.bitwise_and(data ,bitMask) >> bitShift

                # Apply the On-board pixel trim
                data = ma.array(data,mask=onboardTrimMask,fill_value=onboardPixelTrimValue)
                data = data.filled() # Substitute for the masked values with onboardPixelTrimValue

                # Apply the On-ground pixel trim
                data = ma.array(data,mask=ongroundTrimMask,fill_value=ongroundPixelTrimValue)
                data = data.filled() # Substitute for the masked values with ongroundPixelTrimValue

                # Flip the pass depending on whether this is an ascending or decending pass
                data = data[::orient,::orient]

                # Define the colourbar
                cmap = ListedColormap(CMD.ViirsCMfillColours[byte][dSet])

                numCats = np.array(CMD.ViirsCMfillColours[byte][dSet]).size
                numBounds = numCats + 1

                tickPos = np.arange(float(numBounds))/float(numCats)
                tickPos = tickPos[0 :-1] + tickPos[1]/2.

                LOG.info("numCats = {}".format(numCats))
                LOG.info("tickPos = {}".format(tickPos))

                # Set the plot and colourbar titles...
                plotTitle = '%s : orbit %s %s' % (shortName,orbitNumber,annotation)
                cbTitle = plotDescr

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
                ax = fig.add_axes(ax_rect,axis_bgcolor='lightgray')

                # Granule axis title
                ax_title = ppl.setp(ax,title=plotTitle)
                ppl.setp(ax_title,fontsize=12)
                ppl.setp(ax_title,family="sans-serif")

                # Remove the ticks and ticklabels on the main axis
                ppl.setp(ax.get_xticklabels(), visible=False)
                ppl.setp(ax.get_yticklabels(), visible=False)
                ppl.setp(ax.get_xticklines(),visible=False)
                ppl.setp(ax.get_yticklines(),visible=False)

                # Mask the data
                LOG.info("%s is of kind %r" % (shortName,data.dtype.kind))
                if (data.dtype.kind =='i' or data.dtype.kind =='u'):
                    data = ma.masked_greater(data,247)
                else:
                    data = ma.masked_less(data,-800.)

                # Plot the dataset on the main plotting axis
                im = ax.imshow(data,interpolation='nearest',vmin=vmin,vmax=vmax,
                        cmap=cmap)

                # add a colorbar axis
                cax_rect = [0.05 , 0.05, 0.9 , 0.08 ] # [left,bottom,width,height]
                cax = fig.add_axes(cax_rect,frameon=False) # setup colorbar axes

                # Plot the colorbar.
                cb = fig.colorbar(im, cax=cax, orientation='horizontal')
                ppl.setp(cax.get_xticklabels(),fontsize=9)
                ppl.setp(cax.get_xticklines(),visible=False)

                # Set the colourbar tick locations and ticklabels
                #ppl.setp(cb.ax,xticks=CMD.ViirsCMTickPos) # In colorbar axis coords (0..1)
                cb.set_ticks(vmax*tickPos) # In data coords (0..3)
                ppl.setp(cb.ax,xticklabels=CMD.ViirsCMtickNames[byte][dSet])

                # Colourbar title
                cax_title = ppl.setp(cax,title=cbTitle)
                ppl.setp(cax_title,fontsize=10)

                # Turn off the tickmarks on the colourbar
                ppl.setp(cb.ax.get_xticklines(),visible=False)
                ppl.setp(cb.ax.get_xticklabels(),fontsize=9)

                # Redraw the figure
                canvas.draw()

                # Save the figure to a png file...
                pngFile = path.join(pngDir,'%s%s_b%s_%s.png' % (pngPrefix,shortName,orbitNumber,prodName))
                canvas.print_figure(pngFile,dpi=dpi)
                LOG.info("Writing to %s..." % (pngFile))

                ppl.close('all')

                del(data)


    def plot_VCM_granule_tests(self,plotProd='QF',pngDir=None,pngPrefix=None,annotation='',dpi=300):

        if pngDir is None :
            pngDir = path.abspath(path.curdir)

        hdf5_dict = self.hdf5_dict
        collShortNames = hdf5_dict.keys()

        if (plotProd == 'QF'):
            byteList = [0,1,2,3,4,5]
        else :
            byteList = [int(plotProd.strip('QF'))-1]

        LOG.info('collShortNames = %r' % (collShortNames))

        CMD = viirs_edr_data.CloudMaskData

        for shortName in collShortNames :

            LOG.info('shortName = %s' % (shortName))

            dataName = self.dataName[shortName]

            granID_list =  hdf5_dict[shortName].keys()
            granID_list.sort()

            for granID in granID_list :

                LOG.info('%s --> %s ' % (shortName, granID))
                hdf5Obj = hdf5_dict[shortName][granID][1]

                VIIRS_CM_IP_Gran_0 = hdf5Obj.getNode('/Data_Products/VIIRS-CM-IP/VIIRS-CM-IP_Gran_0')
                dayNightFlag =  getattr(VIIRS_CM_IP_Gran_0.attrs,'N_Day_Night_Flag')[0][0]
                orbitNumber =  getattr(VIIRS_CM_IP_Gran_0.attrs,'N_Beginning_Orbit_Number')[0][0]
                LOG.info('N_Day_Night_Flag = %s' % (dayNightFlag))
                LOG.info('N_Beginning_Orbit_Number = %s' % (orbitNumber))
                orient = -1 if dayNightFlag == 'Day' else 1

                for byte in byteList :

                    LOG.info("")

                    plots = list(CMD.ViirsCMbitMaskNames[byte])
                    for item in plots:
                        if item == 'Spare':
                            plots.remove(item)
                    LOG.info("plots = ",plots)

                    if (plots == []):

                        LOG.info("There are no valid datasets in this byte array, skipping...")

                    else :

                        numPlots = len(plots)
                        numRows = np.ceil(float(numPlots)/2.)

                        figWidth = 15. # inches
                        figHeight = 3. * numPlots/2. # inches
                        fig = Figure(figsize=((figWidth,figHeight)))
                        canvas = FigureCanvas(fig)

                        Ax = []
                        Im = []
                        Txt = []
                        Cb = []

                        plotIdx = 1

                        for dSet in range(np.shape(CMD.ViirsCMbitMasks[byte])[0]):

                            LOG.info('\ndSet = %d' %(dSet))

                            dSetName = '/All_Data/VIIRS-CM-IP_All/QF%s_VIIRSCMIP' % (str(byte+1))
                            byteData = hdf5Obj.getNode(dSetName)[:,:]

                            plotTitle = '%s : %s %s (Byte %d)' % (shortName,granID,annotation,byte)
                            fig.text(0.5, 0.95, plotTitle, fontsize=16, color='black', ha='center', va='bottom', alpha=1.0)

                            if (CMD.ViirsCMbitMaskNames[byte][dSet] == 'Spare') :
                                LOG.info("Skipping dataset with byte = %d, dSet = %d" % (byte, dSet))
                            else :
                                LOG.info("byte = %d, dSet = %d" % (byte, dSet))

                                byteMask = CMD.ViirsCMbitMasks[byte][dSet]
                                byteShift = CMD.ViirsCMbitShift[byte][dSet]

                                LOG.info("byteMask = %d, byteShift = %d, dSetName = %s" % (byteMask, byteShift, dSetName))

                                data = np.bitwise_and(byteData,byteMask) >> byteShift
                                vmin,vmax = CMD.ViirsCMvalues[byte][dSet][0], CMD.ViirsCMvalues[byte][dSet][-1]
                                LOG.info("vmin = %d, vmax = %d" % (vmin, vmax))

                                cmap = ListedColormap(CMD.ViirsCMfillColours[byte][dSet])

                                numCats = np.array(CMD.ViirsCMfillColours[byte][dSet]).size
                                numBounds = numCats + 1

                                tickPos = np.arange(float(numBounds))/float(numCats)
                                tickPos = tickPos[0 :-1] + tickPos[1]/2.

                                LOG.info("numCats = {}".format(numCats))
                                LOG.info("tickPos = {}".format(tickPos))

                                titleStr = CMD.ViirsCMbitMaskNames[byte][dSet]
                                LOG.info("titleStr = %s" % (titleStr))

                                Ax.append(fig.add_subplot(numRows,2,plotIdx))
                                LOG.info("data.dtype.__str__() = %s" % (data.dtype.__str__()))
                                Im.append(Ax[dSet].imshow(data.astype('int')[::orient,::orient], vmin=vmin, vmax=vmax, interpolation='nearest',cmap=cmap))
                                Txt.append(Ax[dSet].set_title(titleStr))

                                ppl.setp(Txt[dSet],fontsize=10)
                                ppl.setp(Ax[dSet].get_xticklines(), visible=False)
                                ppl.setp(Ax[dSet].get_yticklines(), visible=False)
                                ppl.setp(Ax[dSet].get_xticklabels(), visible=False)
                                ppl.setp(Ax[dSet].get_yticklabels(), visible=False)

                                Cb.append(fig.colorbar(Im[dSet], orientation='horizontal', pad=0.05))

                                LOG.info("Cb byte = %d, dSet = %d" % (byte, dSet))
                                LOG.info("CMD.ViirsCMbitMaskNames[%d][%d] = %s" % \
                                        (byte,dSet,CMD.ViirsCMbitMaskNames[byte][dSet]))
                                LOG.info("CMD.ViirsCMfillColours[%d][%d] = %s" % \
                                        (byte,dSet,CMD.ViirsCMfillColours[byte][dSet]))

                                Cb[dSet].set_ticks(vmax*tickPos)
                                ppl.setp(Cb[dSet].ax,xticklabels=CMD.ViirsCMtickNames[byte][dSet])
                                ppl.setp(Cb[dSet].ax.get_xticklabels(),fontsize=6)
                                ppl.setp(Cb[dSet].ax.get_xticklines(),visible=False)

                                plotIdx += 1

                        pngFile = path.join(pngDir,'%s%s_%s_QF%s.png' % (pngPrefix,shortName,granID,str(byte+1)))
                        LOG.info("Writing to %s..." % (pngFile))

                        canvas.draw()
                        canvas.print_figure(pngFile,dpi=dpi)

                        ppl.close('all')

                hdf5Obj.close()


    def plot_VCM_pass_tests(self,plotProd='QF',pngDir=None,pngPrefix=None,annotation='',dpi=300):

        if pngDir is None :
            pngDir = path.abspath(path.curdir)

        hdf5_dict = self.hdf5_dict
        collShortNames = hdf5_dict.keys()

        if (plotProd == 'QF'):
            byteList = [0,1,2,3,4,5]
        else :
            byteList = [int(plotProd.strip('QF'))-1]

        LOG.info('collShortNames = %r' % (collShortNames))

        CMD = viirs_edr_data.CloudMaskData

        for shortName in collShortNames :

            LOG.info('shortName = %s' % (shortName))

            dataName = self.dataName[shortName]

            granID_list =  hdf5_dict[shortName].keys()
            granID_list.sort()


            for byte in byteList :

                LOG.info("")

                plots = list(CMD.ViirsCMbitMaskNames[byte])
                for item in plots:
                    if item == 'Spare':
                        plots.remove(item)
                LOG.info("plots = ",plots)


                if (plots == []):

                    LOG.info("There are no valid datasets in this byte array, skipping...")

                else :

                    for dSet in range(np.shape(CMD.ViirsCMbitMasks[byte])[0]):

                        LOG.info('\ndSet = %d' %(dSet))

                        for granID in granID_list :

                            LOG.info('%s --> %s ' % (shortName, granID))

                            hdf5Obj = hdf5_dict[shortName][granID][1]

                            VIIRS_CM_IP_Gran_0 = hdf5Obj.getNode('/Data_Products/VIIRS-CM-IP/VIIRS-CM-IP_Gran_0')
                            dayNightFlag =  getattr(VIIRS_CM_IP_Gran_0.attrs,'N_Day_Night_Flag')[0][0]
                            orbitNumber =  getattr(VIIRS_CM_IP_Gran_0.attrs,'N_Beginning_Orbit_Number')[0][0]
                            LOG.info('N_Day_Night_Flag = %s' % (dayNightFlag))
                            LOG.info('N_Beginning_Orbit_Number = %s' % (orbitNumber))
                            orient = -1 if dayNightFlag == 'Day' else 1

                            dSetName = '/All_Data/VIIRS-CM-IP_All/QF%s_VIIRSCMIP' % (str(byte+1))
                            byteDataGranule = hdf5Obj.getNode(dSetName)[:,:]

                            # Concatenate the granules.
                            try :
                                byteData = np.vstack((byteData,byteDataGranule))
                                LOG.info("byteData shape = %s\n" %(str(byteData.shape)))
                            except :
                                byteData = byteDataGranule[:,:]
                                LOG.info("byteData shape = %s\n" %(str(byteData.shape)))

                        # What value are the bowtie deletion pixels
                        ongroundPixelTrimValue = trimObj.sdrTypeFill['ONGROUND_PT_FILL'][byteData.dtype.name]
                        LOG.info("Onground Pixel Trim value is {}".format(ongroundPixelTrimValue))
                        onboardPixelTrimValue = trimObj.sdrTypeFill['ONBOARD_PT_FILL'][byteData.dtype.name]
                        LOG.info("Onboard Pixel Trim value is {}\n".format(onboardPixelTrimValue))

                        # Create onboard and onground pixel trim mask arrays, for the total number of
                        # scans in the pass...
                        numGranules = len(granID_list)
                        numScans = numGranules * 48
                        onboardTrimMask = trimObj.createOnboardModTrimArray(nscans=numScans,trimType=bool)
                        ongroundTrimMask = trimObj.createModTrimArray(nscans=numScans,trimType=bool)

                        # Apply the On-board pixel trim
                        byteData = ma.array(byteData,mask=onboardTrimMask,fill_value=ongroundPixelTrimValue)
                        #byteData = byteData.filled() # Substitute for the masked values with ongroundPixelTrimValue

                        # Apply the On-board pixel trim
                        byteData = ma.array(byteData,mask=ongroundTrimMask,fill_value=onboardPixelTrimValue)
                        #byteData = byteData.filled() # Substitute for the masked values with onboardPixelTrimValue

                        # Flip the pass depending on whether this is an ascending or decending pass
                        byteData = byteData[::orient,::orient]

                        if (CMD.ViirsCMbitMaskNames[byte][dSet] == 'Spare') :
                            LOG.info("Skipping dataset with byte = %d, dSet = %d" % (byte, dSet))
                        else :
                            LOG.info("byte = %d, dSet = %d" % (byte, dSet))

                            byteMask = CMD.ViirsCMbitMasks[byte][dSet]
                            byteShift = CMD.ViirsCMbitShift[byte][dSet]

                            LOG.info("byteMask = %d, byteShift = %d, dSetName = %s" % (byteMask, byteShift, dSetName))

                            data = np.bitwise_and(byteData,byteMask) >> byteShift
                            vmin,vmax = CMD.ViirsCMvalues[byte][dSet][0], CMD.ViirsCMvalues[byte][dSet][-1]
                            LOG.info("vmin = %d, vmax = %d" % (vmin, vmax))

                            cmap = ListedColormap(CMD.ViirsCMfillColours[byte][dSet])

                            numCats = np.array(CMD.ViirsCMfillColours[byte][dSet]).size
                            numBounds = numCats + 1

                            tickPos = np.arange(float(numBounds))/float(numCats)
                            tickPos = tickPos[0 :-1] + tickPos[1]/2.

                            LOG.info("numCats = {}".format(numCats))
                            LOG.info("tickPos = {}".format(tickPos))

                            # Set the plot and colourbar titles...
                            plotTitle = '%s : orbit %s %s (Byte %d)' % (shortName,orbitNumber,annotation,byte)
                            cbTitle = CMD.ViirsCMbitMaskNames[byte][dSet]

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
                            ax = fig.add_axes(ax_rect,axis_bgcolor='lightgray')

                            # Granule axis title
                            ax_title = ppl.setp(ax,title=plotTitle)
                            ppl.setp(ax_title,fontsize=12)
                            ppl.setp(ax_title,family="sans-serif")

                            # Remove the ticks and ticklabels on the main axis
                            ppl.setp(ax.get_xticklabels(), visible=False)
                            ppl.setp(ax.get_yticklabels(), visible=False)
                            ppl.setp(ax.get_xticklines(),visible=False)
                            ppl.setp(ax.get_yticklines(),visible=False)

                            # Mask the data
                            if (data.dtype.kind =='i' or data.dtype.kind =='u'):
                                LOG.info("%s is of kind %r" % (shortName,data.dtype.kind))
                                data = ma.masked_greater(data,247)
                            else:
                                LOG.info("%s is of kind %r" % (shortName,data.dtype.kind))
                                data = ma.masked_less(data,-800.)

                            # Plot the dataset on the main plotting axis
                            im = ax.imshow(data,interpolation='nearest',vmin=vmin,vmax=vmax,cmap=cmap)

                            # add a colorbar axis
                            cax_rect = [0.05 , 0.05, 0.9 , 0.08 ] # [left,bottom,width,height]
                            cax = fig.add_axes(cax_rect,frameon=False) # setup colorbar axes

                            # Plot the colorbar.
                            cb = fig.colorbar(im, cax=cax, orientation='horizontal')
                            ppl.setp(cax.get_xticklabels(),fontsize=9)
                            ppl.setp(cax.get_xticklines(),visible=False)

                            # Set the colourbar tick locations and ticklabels
                            #ppl.setp(cb.ax,xticks=CMD.ViirsCMTickPos) # In colorbar axis coords (0..1)
                            cb.set_ticks(vmax*tickPos) # In data coords (0..3)
                            ppl.setp(cb.ax,xticklabels=CMD.ViirsCMtickNames[byte][dSet])

                            # Colourbar title
                            cax_title = ppl.setp(cax,title=cbTitle)
                            ppl.setp(cax_title,fontsize=10)

                            # Turn off the tickmarks on the colourbar
                            ppl.setp(cb.ax.get_xticklines(),visible=False)
                            ppl.setp(cb.ax.get_xticklabels(),fontsize=9)

                            # Redraw the figure
                            canvas.draw()

                            # Save the figure to a png file...
                            pngFile = path.join(pngDir,'%s%s_b%s_QF%s_%d.png' % (pngPrefix,shortName,orbitNumber,str(byte+1),dSet))
                            canvas.print_figure(pngFile,dpi=dpi)
                            LOG.info("Writing to %s..." % (pngFile))

                            ppl.close('all')

                        del(byteData)



###################################################
#                  Main Function                  #
###################################################

def main():

    prodChoices=['IP','QF','CMask','CPhase','QF1','QF2','QF3','QF4','QF5','QF6']

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

    mandatoryGroup.add_option('-i','--input_files',
                      action="store",
                      dest="hdf5Files" ,
                      type="string",
                      help="The fully qualified path to the input IICMO HDF5 files. May be a directory or a file glob.")

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
    optionalGroup.add_option('--pass',
                      action="store_true",
                      dest="plotPass",
                      help="Concatenate the granules")
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
    optionalGroup.add_option('-a','--map_annotation',
                      action="store",
                      dest="mapAnn",
                      #default='',
                      type="string",
                      help="The map legend describing the dataset being shown. [default: IPPROD]")
    optionalGroup.add_option('-p','--product',
                      action="store",
                      dest="plotProduct",
                      type="choice",
                      choices=prodChoices,
                      help='''The VIIRS CM IP or QF datasets to plot.\n\n
                           Possible values are...
                           %s
                           ''' % (prodChoices.__str__()[1:-1]))
    optionalGroup.add_option('--png_dir',
                      action="store",
                      dest="pngDir" ,
                      type="string",
                      help="The directory where png files will be written.")

    optionalGroup.add_option('-o','--output_file_prefix',
                      action="store",
                      dest="outputFilePrefix",
                      default="",
                      type="string",
                      help="""String to prefix to the automatically generated 
                      png names, which are of the form 
                      <N_Collection_Short_Name>_<N_Granule_ID>_<dset>.png. 
                      [default: %default]""")

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

    # Check that all of the mandatory options are given. If one or more 
    # are missing, LOG.info(error message and exit...
    mandatories = ['hdf5Files']
    mand_errors = ["Missing mandatory argument [-i HDF5FILES | --input_files=HDF5FILES]"
                  ]
    isMissingMand = False
    for m,m_err in zip(mandatories,mand_errors):
        if not options.__dict__[m]:
            isMissingMand = True
            LOG.info(m_err)
    if isMissingMand :
        parser.error("Incomplete mandatory arguments, aborting...")

    vmin = options.plotMin
    vmax = options.plotMax

    hdf5Path = path.abspath(path.expanduser(options.hdf5Files))

    LOG.info("hdf5Path = %s" % (hdf5Path))
    pngDir = '.' if (options.pngDir is None) else options.pngDir
    pngDir = path.abspath(path.expanduser(pngDir))
    LOG.info("pngDir = %s" % (pngDir))
    if not path.isdir(pngDir):
        LOG.info("Output image directory %s does not exist, creating..." % (pngDir))
        try:
            mkdir(pngDir,0755)
        except Exception, err :
            LOG.info("%s" % (err))
            LOG.info("Creating directory %s failed, aborting..." % (pngDir))
            sys.exit(1)

    pngPrefix = options.outputFilePrefix
    dpi = options.dpi
    plotProduct = options.plotProduct
    plotPass = options.plotPass

    plotEDR = False
    plotQF = False

    if (plotProduct is None):
        plotEDR = True
        plotQF = True
        edrPlotProduct = 'IP'
        qfPlotProduct = 'QF'
    else :
        if ('IP' in plotProduct) \
           or ('CMask' in plotProduct) \
           or ('CPhase' in plotProduct) :
            plotEDR = True
            edrPlotProduct = plotProduct

        if ('QF' in plotProduct) :
            plotQF = True
            qfPlotProduct = plotProduct

    if plotEDR :
        try :
            VCMobj = VCMclass(hdf5Path)
            if plotPass :
                VCMobj.plot_VCM_pass(plotProd=edrPlotProduct,pngDir=pngDir,pngPrefix=pngPrefix,dpi=dpi)
            else:
                VCMobj.plot_VCM_granules(plotProd=edrPlotProduct,pngDir=pngDir,pngPrefix=pngPrefix,dpi=dpi)

            #pytables.file.close_open_files()
        except Exception, err:
            LOG.info(traceback.format_exc())            
            pytables.file.close_open_files()

    if plotQF :
        try :
            VCMobj = VCMclass(hdf5Path)
            if plotPass :
                VCMobj.plot_VCM_pass_tests(plotProd=qfPlotProduct,pngDir=pngDir,pngPrefix=pngPrefix,dpi=dpi)
            else:
                VCMobj.plot_VCM_granule_tests(plotProd=qfPlotProduct,pngDir=pngDir,pngPrefix=pngPrefix,dpi=dpi)

            pytables.file.close_open_files()
        except Exception, err:
            LOG.info(traceback.format_exc())            
            pytables.file.close_open_files()

    LOG.info("Exiting...")
    sys.exit(0)


if __name__ == '__main__':
    main()

