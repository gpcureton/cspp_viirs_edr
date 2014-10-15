#!/usr/bin/env python
# encoding: utf-8
"""
CSPP_VIIRS_SDR_compare.py

Purpose: Create a histogram of the VIIRS SDR, for one two sets of one or more 
         SVM*.h5 granules.

Minimum commandline...

export CSPP_EDR_HOME=$(readlink -f /path/to/EDR)
. $CSPP_EDR_HOME/cspp_edr_env.sh

python CSPP_VIIRS_SDR_compare.py  --input_file_1=dir1/IVAOT_*.h5 \
        --input_file_2=dir2/IVAOT_*.h5


Created by Geoff Cureton on 2014-10-14.
Copyright (c) 2014 University of Wisconsin SSEC. All rights reserved.
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
from datetime import datetime

import numpy as np
from  numpy import ma as ma
import scipy as scipy

import matplotlib
import matplotlib.cm as cm
from matplotlib.colors import Colormap,normalize,LinearSegmentedColormap,\
        ListedColormap,LogNorm
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
import h5py


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
    LOG.debug('hdf5Path = {}'.format(hdf5Path))

    hdf5Dir = path.dirname(hdf5Path)
    hdf5Glob = path.basename(hdf5Path)

    LOG.debug('hdf5Dir = {}'.format(hdf5Dir))
    LOG.debug('hdf5Glob = {}'.format(hdf5Glob))

    if (hdf5Glob == '' or hdf5Glob == '*'):
        LOG.debug('prefix = {}'.format(filePrefix))
        hdf5Glob = path.join(hdf5Dir,'%s_*.h5'%(filePrefix))
    else :
        hdf5Glob = path.join(hdf5Dir,'%s'%(hdf5Glob))

    LOG.debug('Final hdf5Glob = {}'.format(hdf5Glob))

    hdf5Files = glob(hdf5Glob)
    LOG.debug('Final hdf5Files = {}'.format(hdf5Files))

    if hdf5Files != []:
        hdf5Files.sort()
        granIdDict = {}
        for files in hdf5Files :

            # Open the hdf5 file
            fileObj = h5py.File(files,'r')

            all_data_group = fileObj['/All_Data']
            leaf_name = all_data_group.items()[0][0]
            LOG.debug('All_Data leaf name = {}'.format(leaf_name))
            all_data_leaf = all_data_group.items()[0][1]

            data_products_group = fileObj['/Data_Products']
            leaf_name = data_products_group.items()[0][0]
            LOG.debug('Data_Products leaf name = {}'.format(leaf_name))
            data_products_leaf = data_products_group.items()[0][1]
            SDR_Gran_0_name = data_products_leaf.items()[1][0]

            # Get a "pointer" to the granules attribute group.
            node_path = path.join('/Data_Products',leaf_name,SDR_Gran_0_name)
            SDR_Gran_0 = fileObj[node_path]

            # Retrieve a few attributes
            granID = SDR_Gran_0.attrs['N_Granule_ID'][0][0]
            LOG.debug('N_Granule_ID = {}'.format(granID))

            dayNightFlag = SDR_Gran_0.attrs['N_Day_Night_Flag'][0][0]
            LOG.debug('N_Day_Night_Flag = {}'.format(dayNightFlag))

            group_name = path.join('/Data_Products',leaf_name)
            shortName = fileObj[group_name].attrs['N_Collection_Short_Name'][0][0]
            LOG.debug('N_Collection_Short_Name = {}'.format(shortName))

            # Strip the path from the filename
            hdf5File = path.basename(files)

            # Add the granule information to the dictionary, keyed with the granule ID...
            granIdDict[granID] = [hdf5File,fileObj]

        shortNameDict[shortName] = granIdDict

    return shortNameDict


class SDRclass():

    def __init__(self,hdf5Dir_1,hdf5Dir_2):

        self.hdf5Dir_1 = hdf5Dir_1
        self.hdf5Dir_2 = hdf5Dir_2

        self.hdf5_dict_1 = get_hdf5_dict(hdf5Dir_1,'SVM07')
        self.hdf5_dict_2 = get_hdf5_dict(hdf5Dir_2,'SVM07')

        self.collShortNames_1 = self.hdf5_dict_1.keys()
        self.collShortNames_2 = self.hdf5_dict_2.keys()

        LOG.info('collShortNames_1 = {}'.format(self.collShortNames_1))
        LOG.info('collShortNames_2 = {}'.format(self.collShortNames_2))


    def SDR_histogram_generate(self,plotProd='Radiance',
            histBins=20,histMin=None,histMax=None):

        hdf5_dict_1 = self.hdf5_dict_1
        hdf5_dict_2 = self.hdf5_dict_2
        collShortNames_1 = self.collShortNames_1
        collShortNames_2 = self.collShortNames_2

        LOG.info('collShortNames_1 = {}'.format(collShortNames_1))
        LOG.info('collShortNames_2 = {}'.format(collShortNames_2))

        for shortName in collShortNames_1 :

            LOG.info('shortName = {}'.format(shortName))
            LOG.info('plotProd = {}'.format(plotProd))

            granID_list_1 =  hdf5_dict_1[shortName].keys()
            granID_list_1.sort()
            granID_list_2 =  hdf5_dict_2[shortName].keys()
            granID_list_2.sort()

            # Create onboard and onground pixel trim mask arrays, for 
            # the total number of scans in the pass...
            numScans = 48
            onboardTrimMask = trimObj.createOnboardModTrimArray(
                    nscans=numScans,trimType=bool)
            ongroundTrimMask = trimObj.createOngroundModTrimArray(
                    nscans=numScans,trimType=bool)

            H_list = []

            # Read in the data from the granules and concatenate
            for granID_1,granID_2 in zip(granID_list_1,granID_list_2) :

                assert granID_1==granID_2,\
                        "Granule IDs ({},{}) for these granules are not equal"\
                        .format(granID_1,granID_2)

                hdf5Obj_1 = hdf5_dict_1[shortName][granID_1][1]
                hdf5Obj_2 = hdf5_dict_2[shortName][granID_2][1]
                
                # Get the dataset and factors names
                self.dataNames_1 = {}
                all_data_group = hdf5Obj_1['/All_Data']
                leaf_path = all_data_group.items()[0][0]
                LOG.debug('hdf5Obj_1 All_Data leaf name = {}'.format(leaf_path))
                self.dataNames_1[shortName] = ["/All_Data/{}/Reflectance".format(leaf_path), \
                                               "/All_Data/{}/Radiance".format(leaf_path)]

                self.dataNames_2 = {}
                all_data_group = hdf5Obj_2['/All_Data']
                leaf_path = all_data_group.items()[0][0]
                LOG.debug('hdf5Obj_2 All_Data leaf name = {}'.format(leaf_path))
                self.dataNames_2[shortName] = ["/All_Data/{}/Reflectance".format(leaf_path), \
                                               "/All_Data/{}/Radiance".format(leaf_path)]

                dataNames_1 = self.dataNames_1[shortName]
                dataNames_2 = self.dataNames_2[shortName]

                LOG.debug(dataNames_1)
                LOG.debug(dataNames_2)

                for dataNames in dataNames_1:
                    if plotProd in dataNames:
                        dataName_1 = dataNames
                        dataFactorsName_1 = "{}Factors".format(dataName_1)
                        break
                for dataNames in dataNames_2:
                    if plotProd in dataNames:
                        dataName_2 = dataNames
                        dataFactorsName_2 = "{}Factors".format(dataName_2)
                        break

                LOG.debug('dataName_1 = {}'.format(dataName_1))
                LOG.debug('dataFactorsName_1 = {}'.format(dataFactorsName_1))
                LOG.debug('dataName_2 = {}'.format(dataName_2))
                LOG.debug('dataFactorsName_2 = {}'.format(dataFactorsName_2))

                LOG.info('{} --> {},{}'.format(dataName_1,granID_1,granID_2))

                LOG.debug("dataName_1 is {}".format(dataName_1))
                LOG.debug("dataName_2 is {}".format(dataName_2))
                try:
                    data_1 = hdf5Obj_1[dataName_1].value[:,:]
                    data_2 = hdf5Obj_2[dataName_2].value[:,:]
                except Exception, err :
                    LOG.error("{}".format(err))
                    close_hdf5_files(self)
                    sys.exit(-1)

                LOG.debug("data_1 is {}".format(data_1))
                LOG.debug("data_2 is {}".format(data_2))

                # Construct fill masks to cover whatever isn't covered by
                # the bow-tie pixels.
                LOG.debug("{} is of kind {}".format(shortName,data_1.dtype.kind))
                if (data_1.dtype.kind =='i' or data_1.dtype.kind =='u'):
                    LOG.debug("Performing mask of integer types")
                    fill_mask_1 = ma.masked_greater(data_1,65528).mask
                    fill_mask_2 = ma.masked_greater(data_2,65528).mask
                else:
                    LOG.debug("Performing mask of float types")
                    fill_mask_1 = ma.masked_less(data_1,-800.).mask
                    fill_mask_2 = ma.masked_less(data_2,-800.).mask

                # Unscale the datasets
                try:
                    factors_1 = hdf5Obj_1[dataFactorsName_1].value[:]
                    LOG.debug("factors_1 is {}".format(factors_1))
                    LOG.debug("factors_1 is of kind {}".format(factors_1.dtype.kind))
                    data_1 = data_1.astype(factors_1.dtype)
                    data_1 = data_1 * factors_1[0] + factors_1[1]
                except Exception, err :
                    LOG.warn("{}".format(err))
                    LOG.warn("Unable to open {}".format(dataFactorsName_1))

                try:
                    factors_2 = hdf5Obj_2[dataFactorsName_2].value[:]
                    LOG.debug("factors_2 is {}".format(factors_2))
                    LOG.debug("factors_2 is of kind {}".format(factors_2.dtype.kind))
                    data_2 = data_2.astype(factors_2.dtype)
                    data_2 = data_2 * factors_2[0] + factors_2[1]
                except Exception, err :
                    LOG.warn("{}".format(err))
                    LOG.warn("Unable to open {}".format(dataFactorsName_2))

                LOG.debug("data_1 is {}".format(data_1))
                LOG.debug("data_2 is {}".format(data_2))

                # Construct the total masks
                totalMask_1 = fill_mask_1 + onboardTrimMask
                totalMask_2 = fill_mask_1 + onboardTrimMask
                totalMask = np.ravel(totalMask_1 + totalMask_2)

                LOG.debug("totalMask_1.sum() = {}".format(totalMask_1.sum()))
                LOG.debug("totalMask_2.sum() = {}".format(totalMask_2.sum()))
                LOG.debug("totalMask.sum() = {}".format(totalMask.sum()))

                # Flatten the datasets
                data_1 = np.ravel(data_1)
                data_2 = np.ravel(data_2)
                LOG.debug("ravelled data_1.shape is {}".format(data_1.shape))
                LOG.debug("ravelled data_2.shape is {}".format(data_2.shape))

                # Mask the SDR so we only have the radiometric values
                data_1 = ma.masked_array(data_1,mask=totalMask)
                data_2 = ma.masked_array(data_2,mask=totalMask)

                LOG.debug("data_1.mask.sum() = {}".format(data_1.mask.sum()))
                LOG.debug("data_2.mask.sum() = {}".format(data_2.mask.sum()))

                ## Compress the datasets
                data_1 = ma.compressed(data_1)
                data_2 = ma.compressed(data_2)
                LOG.debug("compressed data_1.shape is {}".format(data_1.shape))
                LOG.debug("compressed data_2.shape is {}".format(data_2.shape))

                LOG.debug("data_1 is {}".format(data_1))
                LOG.debug("data_2 is {}".format(data_2))

                ## Generate the histogram for this granule

                LOG.info("Creating histogram...")

                vmin = np.min(data_1) if (histMin == None) else histMin
                vmax = np.max(data_1) if (histMax == None) else histMax
                LOG.debug("vmin is {}".format(vmin))
                LOG.debug("vmax is {}".format(vmax))

                histRange = np.array([[vmin,vmax],[vmin,vmax]])

                H, xedges, yedges = np.histogram2d(data_2,data_1,
                        bins=histBins,range=histRange,normed=False)

                H_list.append(H)

        return H_list, xedges, yedges


def _histogramPlot(xedges, yedges,histogram, 
        vmin=None,vmax=None,histMin=None,histMax=None,scale=None,
        axis_label_1=None,axis_label_2=None,plot_title=r'',pngDpi=300, 
        cmap=None, pngName='SDR_hist.png'):

    figWidth = 5. # inches
    figHeight = 4.2 # inches

    # Create figure with default size, and create canvas to draw on
    fig = Figure(figsize=((figWidth,figHeight)))
    canvas = FigureCanvas(fig)

    ax_len = 0.80
    ax_x_len = ax_y_len = ax_len

    x0,y0 = 0.07,0.10
    x1,y1 = x0+ax_len,y0+ax_len

    cax_x_pad = 0.0
    cax_x_len = 0.05
    cax_y_len = ax_len

    ax_rect  = [x0, y0, ax_len , ax_len  ] # [left,bottom,width,height]
    cax_rect = [x1+cax_x_pad , y0, cax_x_len , cax_y_len ] # [left,bottom,width,height]

    LOG.debug("ax_rect = {}".format(ax_rect))
    LOG.debug("cax_rect = {}".format(cax_rect))

    timeString = 'Creation date: %s' %(datetime.strftime(datetime.utcnow(),"%Y-%m-%d %H:%M:%S Z"))
    fig.text(0.98, 0.01, timeString,fontsize=5, color='gray',ha='right',va='bottom',alpha=0.9)

    # Set the histogram ranges
    histRange = [histMin, histMax]
    LOG.debug("_histogramPlot Histogram range: {}".format(histRange))

    countsRange = [vmin, vmax]
    LOG.debug("_histogramPlot Counts range: {}".format(countsRange))

    # Create main axes instance, leaving room for colorbar at bottom,
    # and also get the Bbox of the axes instance
    ax = fig.add_axes(ax_rect)
    Nbins = len(xedges) - 1
    parity = np.linspace(histMin,histMax,Nbins)
    parLine = ax.plot(parity,parity,'--')
    ppl.setp(parLine,color='gray')

    # New style...
    #ax.set_ylim(0.8,-0.05)
    #ax.set_ylim(ax.get_ylim()[::-1])
    #print "ax.get_ylim() = ",ax.get_ylim()
    # Old Style
    #ax.set_ylim(-0.05,0.8,-1)
    LOG.debug("ax.get_ylim() = {}".format(ax.get_ylim()))
    
    ppl.setp(ax.get_xticklabels(),fontsize=6)
    ppl.setp(ax.get_yticklabels(),fontsize=6)
    ppl.setp(ax,xlabel=axis_label_1)
    ppl.setp(ax,ylabel=axis_label_2)
    ax_title = ppl.setp(ax,title=plot_title)

    # The extent values just change the axis limits, they do NOT
    # alter which part of the array gets plotted.
    extent = [yedges[0], yedges[-1], xedges[0], xedges[-1]]
    H = histogram
    H = H.astype(np.float)

    LOG.debug("Histogram.shape = {}".format(H.shape))
    LOG.debug("xedges.shape = {}".format(xedges.shape))
    LOG.debug("yedges.shape = {}".format(yedges.shape))
    LOG.debug("Histogram min,max = {},{}".format(np.min(H),np.max(H)))

    H /= np.max(H)
    LOG.debug("Scaled Histogram min,max = {},{}".format(np.min(H),np.max(H)))

    cs = ax.imshow(H[:,:], extent=extent, vmin=0.001, vmax = 1.,
            interpolation='nearest',origin='lower',norm=LogNorm(vmin=0.001, vmax=1.))

    # add a colorbar.
    cax = fig.add_axes(cax_rect,frameon=False) # setup colorbar axes

    t = [0.001,0.01,0.1,1.]
    cb = fig.colorbar(cs, cax=cax, ticks=t, format='$%.3f$', orientation='vertical')

    ppl.setp(cax.get_yticklabels(),fontsize=6)
    cax_title = ppl.setp(cax,title="counts/counts$_{max}$")
    ppl.setp(cax_title,fontsize=5)

    # Redraw the figure
    canvas.draw()

    # save image 
    LOG.info("Creating the image file {}".format((pngName)))
    canvas.print_figure(pngName,dpi=pngDpi)


def close_hdf5_files(SDRobj):
    # Close the open HDF5 files...
    for dicts in [SDRobj.hdf5_dict_1,SDRobj.hdf5_dict_2]:
        for granID in np.sort(dicts['VIIRS-M7-SDR'].keys()):

            try :
                h5File = dicts['VIIRS-M7-SDR'][granID][0]
                h5Obj  = dicts['VIIRS-M7-SDR'][granID][1]
                LOG.info('Closing file object for {}'.format(h5File))
                h5Obj.close()
            except Exception, err:
                traceback.print_exc(file=sys.stdout)


def _argparse():
    '''
    Method to encapsulate the option parsing and various setup tasks.
    '''

    import argparse

    prodChoices=['Reflectance','Radiance']

    defaults = {
                'axis_label_1'     : r'M07 Reflectance (original)',
                'axis_label_2'     : r'M07 Reflectance (new)',
                'plot_title'       : r'',
                'isLand'           : False,
                'isOcean'          : False,
                'colormap'         : None,
                'histMin'          : None,
                'histMax'          : None,
                'histBins'         : 50,
                'isLogPlot'        : False,
                'dpi'              : 200,
                'mapAnn'           : "",
                'pngDir'           : None,
                'outputFile' : "SDR_hist.png",
                'verbosity'        : 0
                }

    description = \
    '''
    This is a brief description of %prog
    '''
    usage = "usage: %prog [mandatory args] [options]"
    version = __version__

    parser = argparse.ArgumentParser()

    # Mandatory arguments

    parser.add_argument('--input_file_1',
                      action='store',
                      dest='input_file_1',
                      type=str,
                      required=True,
                      help='''The fully qualified path to the first set of input 
                              files. May be a directory or a file glob.'''
                      )

    parser.add_argument('--input_file_2',
                      action='store',
                      dest='input_file_2',
                      type=str,
                      required=True,
                      help='''The fully qualified path to the second set of input 
                              files. May be a directory or a file glob.'''
                      )

    # Optional arguments

    parser.add_argument('--axis_label_1',
                      action="store",
                      dest="axis_label_1" ,
                      default=defaults["axis_label_1"],
                      type=str,
                      help='''The label of the ordinate axis.'''
                      )

    parser.add_argument('--axis_label_2',
                      action="store",
                      dest="axis_label_2" ,
                      default=defaults["axis_label_2"],
                      type=str,
                      help='''The label of the abcissa axis.'''
                      )

    parser.add_argument('--plot_title',
                      action="store",
                      dest="plot_title" ,
                      default=defaults["plot_title"],
                      type=str,
                      help='''The plot title.'''
                      )

    parser.add_argument('--land',
                      action="store_true",
                      dest="isLand",
                      default=defaults["isLand"],
                      help='''Select the pixels which occur over land.'''
                      )

    parser.add_argument('--ocean',
                      action="store_true",
                      dest="isOcean",
                      default=defaults["isOcean"],
                      help='''Select the pixels which occur over ocean.'''
                      )

    parser.add_argument('--colormap',
                      action="store",
                      dest="colormap" ,
                      default=defaults["colormap"],
                      type=str,
                      help='''The color map used for the histogram.'''
                      )

    parser.add_argument('--histMin',
                      action="store",
                      type=float,
                      dest="histMin",
                      default=defaults["histMin"],
                      help='''Minimum histogram value.'''
                      )

    parser.add_argument('--histMax',
                      action="store",
                      type=float,
                      dest="histMax",
                      default=defaults["histMax"],
                      help='''Maximum value of histogram.'''
                      )

    parser.add_argument('--histBins',
                      action="store",
                      type=int,
                      dest="histBins",
                      default=defaults["histBins"],
                      help='''Number of histogram levels.'''
                      )

    parser.add_argument('--logplot',
                      action="store_true",
                      dest="isLogPlot",
                      default=defaults["isLogPlot"],
                      help='''Plot the data product on a logarithmic scale.'''
                      )

    parser.add_argument('-d','--dpi',
                      action="store",
                      dest="dpi",
                      default=defaults["dpi"],
                      type=float,
                      help='''The resolution in dots per inch of the output png file'''
    )

    parser.add_argument('-a','--map_annotation',
                      action="store",
                      dest="mapAnn",
                      default=defaults["mapAnn"],
                      type=str,
                      help='''The map legend describing the dataset being shown.'''
                      )

    parser.add_argument('-p','--product',
                      action="store",
                      dest="plotProduct",
                      type=str,
                      choices=prodChoices,
                      help='''The VIIRS SDR datasets to plot.\n\n
                           Possible values are...
                           {}
                           '''.format(prodChoices.__str__()[1:-1])
                           )

    parser.add_argument('--png_dir',
                      action="store",
                      dest="pngDir" ,
                      default=defaults["pngDir"],
                      type=str,
                      help='''The directory where png files will be written.'''
                      )

    parser.add_argument('-o','--output_file',
                      action="store",
                      dest="outputFile",
                      default=defaults["outputFile"],
                      type=str,
                      help='''Output file name.png. 
                      '''
                      )

    parser.add_argument("-v", "--verbose",
                      dest='verbosity',
                      action="count", 
                      default=defaults["verbosity"],
                      help='''each occurrence increases verbosity 1 level from 
                      ERROR: -v=WARNING -vv=INFO -vvv=DEBUG''')


    args = parser.parse_args()


    # Set up the logging
    console_logFormat = '%(asctime)s : (%(levelname)s):%(filename)s:%(funcName)s:%(lineno)d:  %(message)s'
    dateFormat = '%Y-%m-%d %H:%M:%S'
    levels = [logging.ERROR, logging.WARN, logging.INFO, logging.DEBUG]
    logging.basicConfig(level = levels[args.verbosity], 
            format = console_logFormat, 
            datefmt = dateFormat)


    return args


###################################################
#                  Main Function                  #
###################################################

def main():
    '''
    The main method.
    '''

    options = _argparse()

    input_file_1 = options.input_file_1
    input_file_2 = options.input_file_2

    axis_label_1     = options.axis_label_1
    axis_label_2     = options.axis_label_2
    plot_title       = options.plot_title
    plotLand         = options.isLand
    plotOcean        = options.isOcean
    colormap         = options.colormap
    histMin          = options.histMin
    histMax          = options.histMax
    histBins         = options.histBins
    isLogPlot        = options.isLogPlot
    dpi              = options.dpi
    mapAnn           = options.mapAnn
    plotProduct      = options.plotProduct
    pngDir           = options.pngDir
    outputFile       = options.outputFile
    verbosity        = options.verbosity

    LOG.debug("Input option 'plotLand'         = {}".format(plotLand))
    LOG.debug("Input option 'plotOcean'        = {}".format(plotOcean))
    LOG.debug("Input option 'colormap'         = {}".format(colormap))
    LOG.debug("Input option 'histMin'          = {}".format(histMin))
    LOG.debug("Input option 'histMax'          = {}".format(histMax))
    LOG.debug("Input option 'histBins'         = {}".format(histBins))
    LOG.debug("Input option 'isLogPlot'        = {}".format(isLogPlot))
    LOG.debug("Input option 'dpi'              = {}".format(dpi))
    LOG.debug("Input option 'mapAnn'           = {}".format(mapAnn))
    LOG.debug("Input option 'plotProduct'      = {}".format(plotProduct))
    LOG.debug("Input option 'pngDir'           = {}".format(pngDir))
    LOG.debug("Input option 'outputFile'       = {}".format(outputFile))
    LOG.debug("Input option 'verbosity'        = {}".format(verbosity))

    hdf5Path_1 = path.abspath(path.expanduser(input_file_1))
    hdf5Path_2 = path.abspath(path.expanduser(input_file_2))

    LOG.debug("hdf5Path_1 = {}".format(hdf5Path_1))
    LOG.debug("hdf5Path_2 = {}".format(hdf5Path_2))


    pngDir = '.' if (pngDir is None) else pngDir
    pngDir = path.abspath(path.expanduser(pngDir))
    LOG.info("pngDir = {}".format(pngDir))
    if not path.isdir(pngDir):
        LOG.info("Output image directory {} does not exist, creating...".format(pngDir))
        try:
            mkdir(pngDir,0755)
        except Exception, err :
            LOG.info("{}".format(err))
            LOG.info("Creating directory {} failed, aborting...".format(pngDir))
            sys.exit(1)

    plotEDR = False

    if (plotProduct is None):
        edrPlotProduct = 'Reflectance'
    else :
        edrPlotProduct = plotProduct

    try :
        SDRobj = SDRclass(hdf5Path_1,hdf5Path_2)


        H_list, xedges, yedges = SDRobj.SDR_histogram_generate(
                plotProd=edrPlotProduct,
                histBins=histBins,histMin=histMin,histMax=histMax)

        H = np.zeros(H_list[0].shape,dtype=H_list[0].dtype)
        for hists in H_list:
            H = H + hists

        histRange = [xedges[0], xedges[-1]]
        LOG.info("_histogramPlot Histogram range: {}".format(histRange))

        plot_options={}
        plot_options['histMin'] = xedges[0]
        plot_options['histMax'] = xedges[-1]
        plot_options['axis_label_1'] = axis_label_1
        plot_options['axis_label_2'] = axis_label_2
        plot_options['plot_title'] = plot_title
        plot_options['pngDpi'] = 300
        plot_options['cmap'] = None
        plot_options['pngName'] = outputFile

        _histogramPlot(xedges,yedges,H,**plot_options)

    except Exception, err:
        traceback.print_exc(file=sys.stdout)

    # Close the open HDF5 files...
    close_hdf5_files(SDRobj)

    print "Exiting..."
    sys.exit(0)


if __name__ == '__main__':
    main()

