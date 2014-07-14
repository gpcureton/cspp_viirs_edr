#!/usr/bin/env python
# encoding: utf-8
"""
CSPP_VIIRS_LST_hdf5Plots.py

Purpose: Create swath projection quicklook PNGs from the VIIRS LST EDR HDF5 files.
         Images can be created for the EDR product or the associated quality flags.

Minimum commandline...

export CSPP_EDR_HOME=$(readlink -f /path/to/EDR)
source $CSPP_EDR_HOME/cspp_edr_env.sh

python CSPP_VIIRS_LST_hdf5Plots.py -i '/path/to/files/VLSTO*.h5'

  or

python CSPP_VIIRS_LST_hdf5Plots.py --input_files=/path/to/files/VLSTO*.h5


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
    print "hdf5Path = %s" % (hdf5Path)

    hdf5Dir = path.dirname(hdf5Path)
    hdf5Glob = path.basename(hdf5Path)

    print "hdf5Dir = %s" % (hdf5Dir)
    print "hdf5Glob = %s" % (hdf5Glob)

    if (hdf5Glob == '' or hdf5Glob == '*'):
        print "prefix = %s" % (filePrefix)
        hdf5Glob = path.join(hdf5Dir,'%s_*.h5'%(filePrefix))
    else :
        hdf5Glob = path.join(hdf5Dir,'%s'%(hdf5Glob))

    print "Final hdf5Glob = %s" % (hdf5Glob)

    hdf5Files = glob(hdf5Glob)
    if hdf5Files != []:
        hdf5Files.sort()
        granIdDict = {}
        for files in hdf5Files :

            # Open the hdf5 file
            fileObj = pytables.openFile(files)

            # Get a "pointer" to the granules attribute group.
            VIIRS_LST_EDR_Gran_0 = fileObj.getNode('/Data_Products/VIIRS-LST-EDR/VIIRS-LST-EDR_Gran_0')

            # Retrieve a few attributes
            granID =  getattr(VIIRS_LST_EDR_Gran_0.attrs,'N_Granule_ID')[0][0]
            print 'N_Granule_ID = %s' % (granID)

            dayNightFlag =  getattr(VIIRS_LST_EDR_Gran_0.attrs,'N_Day_Night_Flag')[0][0]
            print 'N_Day_Night_Flag = %s' % (dayNightFlag)

            shortName = fileObj.getNodeAttr('/Data_Products/VIIRS-LST-EDR','N_Collection_Short_Name')[0][0]
            print 'N_Collection_Short_Name = %s' % (shortName)

            # Strip the path from the filename
            hdf5File = path.basename(files)

            # Add the granule information to the dictionary, keyed with the granule ID...
            granIdDict[granID] = [hdf5File,fileObj]

        shortNameDict[shortName] = granIdDict

    return shortNameDict


class LSTclass():

    def __init__(self,hdf5Dir):

        self.hdf5Dir = hdf5Dir

        self.collShortNames = [
                               'VIIRS-LST-EDR',
                              ]

        self.plotDescr = {}
        self.plotDescr['VIIRS-LST-EDR'] = ['Land Surface Temperature (K)']

        self.plotLims = {}
        self.plotLims['VIIRS-LST-EDR'] = [None,None]

        self.dataName = {}
        self.dataName['VIIRS-LST-EDR'] = ['/All_Data/VIIRS-LST-EDR_All/LandSurfaceTemperature']

        self.dataFactors = {}
        self.dataFactors['VIIRS-LST-EDR'] = ['/All_Data/VIIRS-LST-EDR_All/LSTFactors']

        self.hdf5_dict = get_hdf5_dict(hdf5Dir,'VLSTO')


    def plot_LST_granules(self,plotProd='EDR',vmin=None,vmax=None,pngDir=None,pngPrefix=None,annotation='',dpi=300):

        if pngDir is None :
            pngDir = path.abspath(path.curdir)

        plotDescr = self.plotDescr
        plotLims = self.plotLims

        hdf5_dict = self.hdf5_dict
        collShortNames = hdf5_dict.keys()

        print 'collShortNames = %r' % (collShortNames)

        for shortName in collShortNames :

            print 'shortName = %s' % (shortName)

            if (plotProd == 'EDR'):

                dataNames = self.dataName[shortName]
                factorsNames = self.dataFactors[shortName]
                plotDescrs = plotDescr[shortName]
                prodNames = ['LST']


            granID_list =  hdf5_dict[shortName].keys()
            granID_list.sort()

            for granID in granID_list :

                print '%s --> %s ' % (shortName, granID)

                hdf5Obj = hdf5_dict[shortName][granID][1]

                VIIRS_LST_EDR_Gran_0 = hdf5Obj.getNode('/Data_Products/VIIRS-LST-EDR/VIIRS-LST-EDR_Gran_0')
                dayNightFlag =  getattr(VIIRS_LST_EDR_Gran_0.attrs,'N_Day_Night_Flag')[0][0]
                print 'N_Day_Night_Flag = %s' % (dayNightFlag)
                orient = -1 if dayNightFlag == 'Day' else 1

                for dataName,factorsName,plotDescr,prodName in zip(dataNames,factorsNames,plotDescrs,prodNames):

                    data = hdf5Obj.getNode(dataName)[:,:]
                    factors = hdf5Obj.getNode(factorsName)[:]
                    data = data*factors[0] + factors[1]

                    LSTqualFlag = hdf5Obj.getNode('/All_Data/VIIRS-LST-EDR_All/QF1_VIIRSLSTEDR')
                    LSTqualFlag = np.bitwise_and(LSTqualFlag,3) >> 0
                    LSTqualFlagMask = ma.masked_equal(LSTqualFlag,0).mask

                    pixelTrimValue = trimObj.sdrTypeFill['ONGROUND_PT_FILL'][data.dtype.name]
                    print "pixelTrimValue is %r" % (pixelTrimValue)

                    # Apply the moderate pixel trim, so that we can properly mask them out at plot time.
                    data = ma.array(data,mask=modTrimMask,fill_value=trimObj.sdrTypeFill['ONBOARD_PT_FILL'][data.dtype.name])
                    data = data.filled()

                    plotTitle = '%s : %s %s' % (shortName,granID,annotation)
                    cbTitle = plotDescr

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
                    print "%s is of kind %r" % (shortName,data.dtype.kind)
                    if (data.dtype.kind =='i' or data.dtype.kind =='u'):
                        fill_mask = ma.masked_greater(data,200).mask
                    else:
                        fill_mask = ma.masked_less(data,-800.).mask

                    # Construct the total mask

                    totalMask = LSTqualFlagMask + fill_mask

                    # Mask the aerosol so we only have the retrievals
                    data = ma.masked_array(data,mask=totalMask)
                    
                    im = ax.imshow(data[::orient,::orient],interpolation='nearest',vmin=vmin,vmax=vmax)
                    
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
                    ppl.setp(cax.get_xticklines(),visible=True)

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
                    print "Writing to %s..." % (pngFile)

                    ppl.close('all')

                hdf5Obj.close()


    def plot_LST_pass(self,plotProd='EDR',vmin=None,vmax=None,pngDir=None,
            pngPrefix=None,annotation='',dpi=300):

        if pngDir is None :
            pngDir = path.abspath(path.curdir)

        plotDescr = self.plotDescr
        plotLims = self.plotLims

        hdf5_dict = self.hdf5_dict
        collShortNames = hdf5_dict.keys()

        print 'collShortNames = %r' % (collShortNames)

        for shortName in collShortNames :

            print 'shortName = %s' % (shortName)

            if (plotProd == 'EDR'):

                dataNames = self.dataName[shortName]
                factorsNames = self.dataFactors[shortName]
                plotDescrs = plotDescr[shortName]
                prodNames = ['LST']


            granID_list =  hdf5_dict[shortName].keys()
            granID_list.sort()

            for dataName,factorsName,plotDescr,prodName in zip(dataNames,factorsNames,plotDescrs,prodNames):

                # Read in the data from the granules and concatenate
                for granID in granID_list :

                    print '%s --> %s ' % (shortName, granID)

                    hdf5Obj = hdf5_dict[shortName][granID][1]

                    VIIRS_LST_EDR_Gran_0 = hdf5Obj.getNode('/Data_Products/VIIRS-LST-EDR/VIIRS-LST-EDR_Gran_0')
                    dayNightFlag =  getattr(VIIRS_LST_EDR_Gran_0.attrs,'N_Day_Night_Flag')[0][0]
                    orbitNumber =  getattr(VIIRS_LST_EDR_Gran_0.attrs,'N_Beginning_Orbit_Number')[0][0]
                    print 'N_Day_Night_Flag = %s' % (dayNightFlag)
                    print 'N_Beginning_Orbit_Number = %s' % (orbitNumber)
                    orient = -1 if dayNightFlag == 'Day' else 1

                    dataGranule = hdf5Obj.getNode(dataName)[:,:]
                    factors = hdf5Obj.getNode(factorsName)[:]
                    qf1ArrGranule = hdf5Obj.getNode('/All_Data/VIIRS-LST-EDR_All/QF1_VIIRSLSTEDR')
                    qf2ArrGranule = hdf5Obj.getNode('/All_Data/VIIRS-LST-EDR_All/QF2_VIIRSLSTEDR')
                    qf3ArrGranule = hdf5Obj.getNode('/All_Data/VIIRS-LST-EDR_All/QF3_VIIRSLSTEDR')

                    # Concatenate the granules.
                    try :
                        data = np.vstack((data,dataGranule))
                        qf1Arr = np.vstack((qf1Arr,qf1ArrGranule))
                        qf2Arr = np.vstack((qf2Arr,qf1ArrGranule))
                        qf3Arr = np.vstack((qf3Arr,qf1ArrGranule))
                        print "data shape = {}".format(data.shape)
                        print "qf1Arr shape = {}".format(qf1Arr.shape)
                        print "qf2Arr shape = {}".format(qf2Arr.shape)
                        print "qf3Arr shape = {}\n".format(qf3Arr.shape)
                    except NameError :
                        data = dataGranule[:,:]
                        qf1Arr = qf1ArrGranule[:,:]
                        qf2Arr = qf2ArrGranule[:,:]
                        qf3Arr = qf3ArrGranule[:,:]
                        print "data shape = {}".format(data.shape)
                        print "qf1Arr shape = {}".format(qf1Arr.shape)
                        print "qf2Arr shape = {}".format(qf2Arr.shape)
                        print "qf3Arr shape = {}\n".format(qf3Arr.shape)

                print "Final data shape = {}".format(data.shape)
                print "Final qf1Arr shape = {}".format(qf1Arr.shape)
                print "Final qf2Arr shape = {}".format(qf2Arr.shape)
                print "Final qf3Arr shape = {}\n".format(qf3Arr.shape)

                try :
                    #Determine masks for each fill type, for the LST EDR
                    lstFillMasks = {}
                    for fillType in trimObj.sdrTypeFill.keys() :
                        fillValue = trimObj.sdrTypeFill[fillType][data.dtype.name]
                        if 'float' in fillValue.__class__.__name__ :
                            lstFillMasks[fillType] = ma.masked_inside(data,fillValue-eps,fillValue+eps).mask
                            if (lstFillMasks[fillType].__class__.__name__ != 'ndarray') :
                                lstFillMasks[fillType] = None
                        elif 'int' in fillValue.__class__.__name__ :
                            lstFillMasks[fillType] = ma.masked_equal(data,fillValue).mask
                            if (lstFillMasks[fillType].__class__.__name__ != 'ndarray') :
                                lstFillMasks[fillType] = None
                        else :
                            print "Dataset was neither int not float... a worry"
                            pass

                    #Construct the total mask from all of the various fill values
                    fillMask = ma.array(np.zeros(data.shape,dtype=np.bool))
                    for fillType in trimObj.sdrTypeFill.keys() :
                        if lstFillMasks[fillType] is not None :
                            fillMask = fillMask * ma.array(np.zeros(data.shape,dtype=np.bool),\
                                mask=lstFillMasks[fillType])


                    # Unscale the LST dataset
                    data = data * factors[0] + factors[1]

                    # LST quality mask
                    LSTqualFlag = np.bitwise_and(qf1Arr,3) >> 0
                    lstQualMask = ma.masked_equal(LSTqualFlag,3).mask

                    # LST valid range mask
                    LSTvalidRangeFlag = np.bitwise_and(qf2Arr,2) >> 1
                    lstValidRangeMask = ma.masked_equal(LSTvalidRangeFlag,1).mask

                    # LST water bodies LWM mask
                    #LSTwaterBackgroundFlag = np.bitwise_and(qf3Arr,7) >> 0
                    #lstWaterBackgroundMask = ma.masked_inside(LSTwaterBackgroundFlag,2,5).mask

                    # LST water bodies surface type mask
                    LSTwaterBodiesFlag = np.bitwise_and(qf3Arr,248) >> 3
                    lstWaterBodiesMask = ma.masked_equal(LSTwaterBodiesFlag,17).mask

                    # Combine the fill mask and quality masks...
                    totalMask = fillMask.mask + lstQualMask \
                            + lstWaterBodiesMask

                    try :
                        data = ma.array(data,mask=totalMask)
                    except ma.core.MaskError :
                        print ">> error: Mask Error, probably mismatched geolocation and product array sizes, aborting..."
                        sys.exit(1)

                except Exception, err :
                    print ">> error: %s..." % (str(err))
                    sys.exit(1)

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

                # Plot the data
                print "%s is of kind %r" % (shortName,data.dtype.kind)
                if (data.dtype.kind =='i' or data.dtype.kind =='u'):
                    fill_mask = ma.masked_greater(data,200).mask
                else:
                    fill_mask = ma.masked_less(data,-800.).mask

                # Construct the total mask

                totalMask = fill_mask# + LSTqualFlagMask

                # Mask the aerosol so we only have the retrievals
                data = ma.masked_array(data,mask=totalMask)
                
                # Flip the pass depending on whether this is an ascending or decending pass
                data = data[::orient,::orient]

                im = ax.imshow(data,interpolation='nearest',vmin=vmin,vmax=vmax)
                
                # add a colorbar axis
                cax_rect = [0.05 , 0.05, 0.9 , 0.08 ] # [left,bottom,width,height]
                cax = fig.add_axes(cax_rect,frameon=False) # setup colorbar axes

                # Plot the colorbar.
                cb = fig.colorbar(im, cax=cax, orientation='horizontal')
                ppl.setp(cax.get_xticklabels(),fontsize=9)
                ppl.setp(cax.get_xticklines(),visible=True)

                # Colourbar title
                cax_title = ppl.setp(cax,title=cbTitle)
                ppl.setp(cax_title,fontsize=10)

                # Turn off the tickmarks on the colourbar
                #ppl.setp(cb.ax.get_xticklines(),visible=False)
                #ppl.setp(cb.ax.get_xticklabels(),fontsize=9)

                # Redraw the figure
                canvas.draw()

                # Save the figure to a png file...
                pngFile = path.join(pngDir,'%s%s_b%s_%s.png' % (pngPrefix,shortName,orbitNumber,prodName))
                canvas.print_figure(pngFile,dpi=dpi)
                print "Writing to %s..." % (pngFile)

                ppl.close('all')

                del(data)
                del(qf1Arr)
                del(qf2Arr)
                del(qf3Arr)


    def plot_LST_tests(self,plotProd='QF',pngDir=None,pngPrefix=None,annotation='',dpi=300):

        if pngDir is None :
            pngDir = path.abspath(path.curdir)

        hdf5_dict = self.hdf5_dict
        collShortNames = hdf5_dict.keys()

        if (plotProd == 'QF'):
            byteList = [0,1,2,3]
        else :
            byteList = [int(plotProd.strip('QF'))]

        print 'collShortNames = %r' % (collShortNames)

        CMD = viirs_edr_data.SeaSurfaceTempProdData

        for shortName in collShortNames :

            print 'shortName = %s' % (shortName)

            granID_list =  hdf5_dict[shortName].keys()
            granID_list.sort()

            for granID in granID_list :

                print '%s --> %s ' % (shortName, granID)
                hdf5Obj = hdf5_dict[shortName][granID][1]

                VIIRS_LST_EDR_Gran_0 = hdf5Obj.getNode('/Data_Products/VIIRS-LST-EDR/VIIRS-LST-EDR_Gran_0')
                dayNightFlag =  getattr(VIIRS_LST_EDR_Gran_0.attrs,'N_Day_Night_Flag')[0][0]
                print 'N_Day_Night_Flag = %s' % (dayNightFlag)
                orient = -1 if dayNightFlag == 'Day' else 1

                for byte in byteList :

                    print ""

                    plots = list(CMD.ViirsLSTqualBitMaskNames[byte])
                    for item in plots:
                        if item == 'Spare':
                            plots.remove(item)
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

                    for dSet in range(np.shape(CMD.ViirsLSTqualBitMasks[byte])[0]):

                        print '\ndSet = %d' %(dSet)

                        dSetName = '/All_Data/VIIRS-LST-EDR_All/QF%s_VIIRSLSTEDR'%(str(byte+1))
                        byteData = hdf5Obj.getNode(dSetName)[:,:]

                        plotTitle = '%s : %s %s (Byte %d)' % (shortName,granID,annotation,byte)
                        fig.text(0.5, 0.95, plotTitle, fontsize=16, color='black', ha='center', va='bottom', alpha=1.0)

                        if (CMD.ViirsLSTqualBitMaskNames[byte][dSet] == 'Spare') :
                            print "Skipping dataset with byte = %d, dSet = %d" % (byte, dSet)
                        else :
                            print "byte = %d, dSet = %d" % (byte, dSet)

                            byteMask = CMD.ViirsLSTqualBitMasks[byte][dSet]
                            byteShift = CMD.ViirsLSTqualBitShift[byte][dSet]

                            print "byteMask = %d, byteShift = %d, dSetName = %s" % (byteMask, byteShift, dSetName)

                            data = np.bitwise_and(byteData,byteMask) >> byteShift
                            vmin,vmax = CMD.ViirsLSTqualvalues[byte][dSet][0], CMD.ViirsLSTqualvalues[byte][dSet][-1]
                            print "vmin = %d, vmax = %d" % (vmin, vmax)

                            cmap = ListedColormap(CMD.ViirsLSTqualFillColours[byte][dSet])

                            numCats = np.array(CMD.ViirsLSTqualFillColours[byte][dSet]).size
                            numBounds = numCats + 1

                            tickPos = np.arange(float(numBounds))/float(numCats)
                            tickPos = tickPos[0 :-1] + tickPos[1]/2.

                            print "numCats = ",numCats
                            print "tickPos = ",tickPos

                            titleStr = CMD.ViirsLSTqualBitMaskNames[byte][dSet]
                            print "titleStr = %s" % (titleStr)

                            Ax.append(fig.add_subplot(numRows,2,plotIdx))
                            print "data.dtype.__str__() = %s" % (data.dtype.__str__())
                            Im.append(Ax[dSet].imshow(data.astype('int')[::orient,::orient], vmin=vmin, vmax=vmax, interpolation='nearest',cmap=cmap))
                            Txt.append(Ax[dSet].set_title(titleStr))

                            ppl.setp(Txt[dSet],fontsize=10)
                            ppl.setp(Ax[dSet].get_xticklines(), visible=False)
                            ppl.setp(Ax[dSet].get_yticklines(), visible=False)
                            ppl.setp(Ax[dSet].get_xticklabels(), visible=False)
                            ppl.setp(Ax[dSet].get_yticklabels(), visible=False)

                            Cb.append(fig.colorbar(Im[dSet], orientation='horizonal', pad=0.05))

                            print "Cb byte = %d, dSet = %d" % (byte, dSet)
                            print "CMD.ViirsLSTqualTickNames[%d][%d] = %s" % \
                                    (byte,dSet,CMD.ViirsLSTqualTickNames[byte][dSet])
                            print "CMD.ViirsLSTqualFillColours[%d][%d] = %s" % \
                                    (byte,dSet,CMD.ViirsLSTqualFillColours[byte][dSet])

                            Cb[dSet].set_ticks(vmax*tickPos)
                            ppl.setp(Cb[dSet].ax,xticklabels=CMD.ViirsLSTqualTickNames[byte][dSet])
                            ppl.setp(Cb[dSet].ax.get_xticklabels(),fontsize=6)
                            ppl.setp(Cb[dSet].ax.get_xticklines(),visible=False)

                            plotIdx += 1

                    pngFile = path.join(pngDir,'%s%s_%s_QF%s.png' % (pngPrefix,shortName,granID,str(byte+1)))
                    print "Writing to %s..." % (pngFile)

                    canvas.draw()
                    canvas.print_figure(pngFile,dpi=dpi)

                    ppl.close('all')

                hdf5Obj.close()


    def plot_LST_pass_tests(self,plotProd='QF',pngDir=None,pngPrefix=None,annotation='',dpi=300):

        if pngDir is None :
            pngDir = path.abspath(path.curdir)

        hdf5_dict = self.hdf5_dict
        collShortNames = hdf5_dict.keys()

        if (plotProd == 'QF'):
            byteList = [0,1,2,3]
        else :
            byteList = [int(plotProd.strip('QF'))-1]

        print 'collShortNames = %r' % (collShortNames)

        CMD = viirs_edr_data.SeaSurfaceTempProdData

        for shortName in collShortNames :

            print 'shortName = %s' % (shortName)

            dataName = self.dataName[shortName]

            granID_list =  hdf5_dict[shortName].keys()
            granID_list.sort()


            for byte in byteList :

                print ""

                plots = list(CMD.ViirsLSTqualBitMaskNames[byte])
                for item in plots:
                    if item == 'Spare':
                        plots.remove(item)
                print "plots = ",plots


                if (plots == []):

                    print "There are no valid datasets in this byte array, skipping..."

                else :

                    for dSet in range(np.shape(CMD.ViirsLSTqualBitMasks[byte])[0]):

                        print '\ndSet = %d' %(dSet)

                        for granID in granID_list :

                            print '%s --> %s ' % (shortName, granID)

                            hdf5Obj = hdf5_dict[shortName][granID][1]

                            VIIRS_LST_EDR_Gran_0 = hdf5Obj.getNode('/Data_Products/VIIRS-LST-EDR/VIIRS-LST-EDR_Gran_0')
                            dayNightFlag =  getattr(VIIRS_LST_EDR_Gran_0.attrs,'N_Day_Night_Flag')[0][0]
                            orbitNumber =  getattr(VIIRS_LST_EDR_Gran_0.attrs,'N_Beginning_Orbit_Number')[0][0]
                            print 'N_Day_Night_Flag = %s' % (dayNightFlag)
                            print 'N_Beginning_Orbit_Number = %s' % (orbitNumber)
                            orient = -1 if dayNightFlag == 'Day' else 1

                            dSetName = '/All_Data/VIIRS-LST-EDR_All/QF%s_VIIRSLSTEDR'%(str(byte+1))
                            byteDataGranule = hdf5Obj.getNode(dSetName)[:,:]

                            # Concatenate the granules.
                            try :
                                byteData = np.vstack((byteData,byteDataGranule))
                                print "byteData shape = %s\n" %(str(byteData.shape))
                            except :
                                byteData = byteDataGranule[:,:]
                                print "byteData shape = %s\n" %(str(byteData.shape))

                        # What value are the bowtie deletion pixels
                        ongroundPixelTrimValue = trimObj.sdrTypeFill['ONGROUND_PT_FILL'][byteData.dtype.name]
                        print "Onground Pixel Trim value is {}".format(ongroundPixelTrimValue)
                        onboardPixelTrimValue = trimObj.sdrTypeFill['ONBOARD_PT_FILL'][byteData.dtype.name]
                        print "Onboard Pixel Trim value is {}\n".format(onboardPixelTrimValue)

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

                        if (CMD.ViirsLSTqualBitMaskNames[byte][dSet] == 'Spare') :
                            print "Skipping dataset with byte = %d, dSet = %d" % (byte, dSet)
                        else :
                            print "byte = %d, dSet = %d" % (byte, dSet)

                            byteMask = CMD.ViirsLSTqualBitMasks[byte][dSet]
                            byteShift = CMD.ViirsLSTqualBitShift[byte][dSet]

                            print "byteMask = %d, byteShift = %d, dSetName = %s" % (byteMask, byteShift, dSetName)

                            data = np.bitwise_and(byteData,byteMask) >> byteShift
                            vmin,vmax = CMD.ViirsLSTqualvalues[byte][dSet][0], CMD.ViirsLSTqualvalues[byte][dSet][-1]
                            print "vmin = %d, vmax = %d" % (vmin, vmax)

                            cmap = ListedColormap(CMD.ViirsLSTqualFillColours[byte][dSet])

                            numCats = np.array(CMD.ViirsLSTqualFillColours[byte][dSet]).size
                            numBounds = numCats + 1

                            tickPos = np.arange(float(numBounds))/float(numCats)
                            tickPos = tickPos[0 :-1] + tickPos[1]/2.

                            print "numCats = ",numCats
                            print "tickPos = ",tickPos

                            # Set the plot and colourbar titles...
                            plotTitle = '%s : orbit %s %s (Byte %d)' % (shortName,orbitNumber,annotation,byte)
                            cbTitle = CMD.ViirsLSTqualBitMaskNames[byte][dSet]

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
                                print "%s is of kind %r" % (shortName,data.dtype.kind)
                                data = ma.masked_greater(data,247)
                            else:
                                print "%s is of kind %r" % (shortName,data.dtype.kind)
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
                            ppl.setp(cb.ax,xticklabels=CMD.ViirsLSTqualTickNames[byte][dSet])

                            # Colourbar title
                            cax_title = ppl.setp(cax,title=cbTitle)
                            ppl.setp(cax_title,fontsize=10)

                            # Turn off the tickmarks on the colourbar
                            ppl.setp(cb.ax.get_xticklines(),visible=False)
                            ppl.setp(cb.ax.get_xticklabels(),fontsize=7)

                            # Redraw the figure
                            canvas.draw()

                            # Save the figure to a png file...
                            pngFile = path.join(pngDir,'%s%s_b%s_QF%s_%d.png' % (pngPrefix,shortName,orbitNumber,str(byte+1),dSet))
                            canvas.print_figure(pngFile,dpi=dpi)
                            print "Writing to %s..." % (pngFile)

                            ppl.close('all')

                        del(byteData)


###################################################
#                  Main Function                  #
###################################################

def main():

    prodChoices=['EDR','QF','LST','QF1','QF2','QF3']

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
                      help="The fully qualified path to the input VLSTO HDF5 files. May be a directory or a file glob.")

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
                      help='''The VIIRS LST EDR or QF datasets to plot.\n\n
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
                      help="""String to prefix to the automatically generated png names, which are of
the form <N_Collection_Short_Name>_<N_Granule_ID>_<dset>.png. [default: %default]""")

    parser.add_option_group(optionalGroup)

    # Parse the arguments from the command line
    (options, args) = parser.parse_args()


    # Check that all of the mandatory options are given. If one or more 
    # are missing, print error message and exit...
    mandatories = ['hdf5Files']
    mand_errors = ["Missing mandatory argument [-i HDF5FILES | --input_files=HDF5FILES]"
                  ]
    isMissingMand = False
    for m,m_err in zip(mandatories,mand_errors):
        if not options.__dict__[m]:
            isMissingMand = True
            print m_err
    if isMissingMand :
        parser.error("Incomplete mandatory arguments, aborting...")

    vmin = options.plotMin
    vmax = options.plotMax

    hdf5Path = path.abspath(path.expanduser(options.hdf5Files))

    print "hdf5Path = %s" % (hdf5Path)
    pngDir = '.' if (options.pngDir is None) else options.pngDir
    pngDir = path.abspath(path.expanduser(pngDir))
    print "pngDir = %s" % (pngDir)
    if not path.isdir(pngDir):
        print "Output image directory %s does not exist, creating..." % (pngDir)
        try:
            mkdir(pngDir,0755)
        except Exception, err :
            print "%s" % (err)
            print "Creating directory %s failed, aborting..." % (pngDir)
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
        edrPlotProduct = 'EDR'
        qfPlotProduct = 'QF'
    else :
        if ('EDR' in plotProduct) \
           or ('Skin' in plotProduct) \
           or ('Bulk' in plotProduct) :
            plotEDR = True
            edrPlotProduct = plotProduct

        if ('QF' in plotProduct) :
            plotQF = True
            qfPlotProduct = plotProduct

    if plotEDR :
        try :
            LSTobj = LSTclass(hdf5Path)
            if plotPass :
                LSTobj.plot_LST_pass(plotProd=edrPlotProduct,
                vmin=vmin,vmax=vmax,pngDir=pngDir,pngPrefix=pngPrefix,dpi=dpi)
            else:
                LSTobj.plot_LST_granules(plotProd=edrPlotProduct,vmin=vmin,vmax=vmax,pngDir=pngDir,pngPrefix=pngPrefix,dpi=dpi)

            pytables.file.close_open_files()
        except Exception, err:
            traceback.print_exc(file=sys.stdout)
            pytables.file.close_open_files()

    if plotQF :
        try :
            LSTobj = LSTclass(hdf5Path)
            if plotPass :
                LSTobj.plot_LST_pass_tests(plotProd=qfPlotProduct,pngDir=pngDir,pngPrefix=pngPrefix,dpi=dpi)
            else:
                LSTobj.plot_LST_tests(plotProd=qfPlotProduct,pngDir=pngDir,pngPrefix=pngPrefix,dpi=dpi)

            pytables.file.close_open_files()
        except Exception, err:
            traceback.print_exc(file=sys.stdout)
            pytables.file.close_open_files()

    print "Exiting..."
    sys.exit(0)


if __name__ == '__main__':
    main()
