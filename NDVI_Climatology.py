#!/usr/bin/env python
# encoding: utf-8
"""
NDVI_Climatology.py

Purpose: This script will compute the min and max global gridded NDVI from the 
         16-day NDVI 
         
         The Filled Normalized Difference Vegetative Index (NDVI)
         Product, which is computed from the (White-Sky) Filled Land Surface
         Albedo Map Product, is a global data set of spatially complete NDVI
         maps for 23 sixteen-day periods per year (001, 017, ... 353). There are
         two types of Filled NDVI Products: 1-minute Map Products and coarser
         resolution Statistical Products.

         Map Products, containing spatially complete NDVI data, are generated at
         1-minute resolution on an equal-angle grid.

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

Minimum commandline:

    python NDVI_Climatology.py  [mandatory options]


Created by Geoff Cureton <geoff.cureton@ssec.wisc.edu> on 2014-04-30.
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

__author__ = 'Geoff Cureton <geoff.cureton@ssec.wisc.edu>'
__version__ = '$Id$'
__docformat__ = 'Epytext'

import os
from os import path,uname,environ
import sys
import logging
import traceback
import string
import re
import uuid
from shutil import rmtree,copyfile
from glob import glob
import copy
from time import time
from datetime import datetime,timedelta

import numpy as np
from numpy import ma

import tables as pytables
import h5py

import matplotlib
import matplotlib.cm as cm
from matplotlib.colors import Colormap, normalize, LinearSegmentedColormap,ListedColormap

from matplotlib.figure import Figure

matplotlib.use('Agg')
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

# This must come *after* the backend is specified.
import matplotlib.pyplot as ppl
from mpl_toolkits.basemap import Basemap,addcyclic,shiftgrid

import viirs_edr_data

# every module should have a LOG object
LOG = logging.getLogger(__file__)


def _create_input_file_globs(inputFiles):
    '''
    Determine the correct input file path and globs
    '''
    input_path = path.abspath(path.expanduser(inputFiles))
    if path.isdir(input_path) :
        input_dir = input_path
        input_files = None
    else :
        input_dir = path.dirname(input_path)
        input_files = path.basename(input_path)

    LOG.debug("input_path = %s" %(input_path))
    LOG.debug("input_dir = %s" %(input_dir))
    LOG.debug("input_files = %s" %(input_files))

    inputGlob = None

    charsToKill = string.ascii_letters + string.digits + "."

    if (input_files is None):
        # Input file glob is of form "/path/to/files"
        LOG.debug('Path1')
        inputGlob = '*.h5'

    elif path.isfile(input_path) :
        # Input file glob is of form "/path/to/files/full_file_name.h5" 
        LOG.debug('Path2')
        fileGlob = string.rstrip(string.lstrip(string.split(input_files,"b")[0],
            charsToKill),charsToKill)
        LOG.debug("fileGlob = %s" %(fileGlob))
        inputGlob = "*%s*.h5" %(fileGlob)
        LOG.debug("Initial inputGlob = %s" %(inputGlob))
        while (string.find(inputGlob,"**")!= -1): 
            inputGlob = string.replace(inputGlob,"**","*")
            LOG.debug("New inputGlob = %s" %(inputGlob))

    elif ("*" in input_files):
        # Input file glob is of form "/path/to/files/something*else"
        LOG.debug('Path3')
        #fileGlob = string.rstrip(string.lstrip(string.split(input_files,"b")[0],
        #   charsToKill),charsToKill)
        fileGlob = input_files
        inputGlob = "{}".format(fileGlob)
        LOG.debug("Initial inputGlob = %s" %(inputGlob))
        while (string.find(inputGlob,"**")!= -1): 
            inputGlob = string.replace(inputGlob,"**","*")
            LOG.debug("New inputGlob = %s" %(inputGlob))

    return input_dir,inputGlob


def _create_ndvi_object_h5py(input_files,dset_dicts,output_file):
    ''' Create a new output file and populate the attributes. '''

    LOG.info("Opening the min/max NDVI file {}".format(path.basename(output_file)))

    f = h5py.File(output_file,'w')

    # Set some global attributes on the output file
    f.attrs['Author'] = __author__
    f.attrs['Source'] = file_HeadURL
    f.attrs['Version'] = __version__
    f.attrs['input 16-day NDVI files'] = np.array(input_files)

    # Copy the datasets and their associated attributes to the 
    # hdf5 file.
    for dset in ['Latitude','Longitude','min_ndvi','max_ndvi']:
        LOG.debug(dset)
        f[dset] = dset_dicts[dset]['data']
        attr_list = dset_dicts[dset].keys()
        attr_list.remove('data')
        LOG.debug(attr_list)
        for attr_key in attr_list:
            f[dset].attrs[attr_key] = dset_dicts[dset][attr_key]
    
    LOG.debug("Closing HDF5 file")

    f.close()


def _read_ndvi_object_h5py(input_file):
    ''' Create a new output file and populate the attributes. '''

    LOG.info("Opening the min/max NDVI file {}".format(path.basename(input_file)))

    f = h5py.File(input_file,'r')

    dset_dicts = {}

    for dset in ['Latitude','Longitude','min_ndvi','max_ndvi']:
        dset_dicts[dset] = {}

        # Get the list of attributes for this dataset
        attr_list = f[dset].attrs.keys()
        LOG.debug(attr_list)
        for attr_key in attr_list:
            dset_dicts[dset][attr_key] = f[dset].attrs[attr_key]

        # Get the array data for this dataset
        dset_dicts[dset]['data'] = f[dset].value

    for dset in dset_dicts.keys():
        LOG.debug(dset)
        for attr_key in dset_dicts[dset].keys():
            LOG.debug('\t{} = {}'.format(attr_key,dset_dicts[dset][attr_key]))

    f.close()

    return dset_dicts


def _plot_ndvi(dset_dicts,dset,title=None,png_name=None,dpi=200):

        Latitude  = dset_dicts['Latitude']['data']
        Longitude = dset_dicts['Longitude']['data']
        dataset = dset_dicts[dset]['data'][::-1,:]
        fill_value = dset_dicts[dset]['_FillValue']
        scale_factor = dset_dicts[dset]['scale_factor']
        offset = dset_dicts[dset]['add_offset']

        dataset_mask = ma.masked_equal(dataset,fill_value).mask
        dataset = scale_factor * dataset.astype('float32') + offset
        dataset = ma.array(dataset,mask=dataset_mask)

        # A default 0.5 degree grid...
        lon,lat = np.meshgrid(Longitude,Latitude)

        # Create figure with default size, and create canvas to draw on
        scale=1.5
        fig = Figure(figsize=(scale*8,scale*5))
        canvas = FigureCanvas(fig)

        # Create main axes instance, leaving room for colorbar at bottom,
        # and also get the Bbox of the axes instance
        ax_rect = [0.05, 0.18, 0.9, 0.75  ] # [left,bottom,width,height]
        ax = fig.add_axes(ax_rect)

        # Granule axis title
        ax_title = ppl.setp(ax,title=title)
        ppl.setp(ax_title,fontsize=12)
        ppl.setp(ax_title,family="sans-serif")

        # Create the basemap object
        m = Basemap(projection='cyl',lon_0=0.,ax=ax)
        x,y = m(lon,lat)

        # Get the colormap
        VegetationIndexProduct \
            = viirs_edr_data.VegetationIndexProdData.VegetationIndexProd()
        cmap = VegetationIndexProduct.cmap

        # Plot the data
        im = m.imshow(dataset,axes=ax,interpolation='nearest',
                vmin=-0.05,vmax=1.0,cmap=cmap)
        
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
        cax_title = ppl.setp(cax,title='NDVI')
        ppl.setp(cax_title,fontsize=9)

        # Redraw the figure
        canvas.draw()

        # save image
        if png_name == None:
            png_name = "{}.png".format(dset)
        LOG.info("Writing image file to {}".format(png_name))
        canvas.print_figure(png_name,dpi=dpi)

        del(m)
        return 0


def _argparse():
    '''
    Method to encapsulate the option parsing and various setup tasks.
    '''

    import argparse

    ndvi_choices=['min_ndvi','max_ndvi']

    defaults = {
                'input_file' : None,
                'plot_ndvi'  : False,
                'stride'     : 1,
                'output_file': 'ann_min_max_ndvi.h5',
                'ndvi_choice': 'min_ndvi',
                'dpi'        : 200.
                }

    description = '''This script will compute the min and max global gridded 
                     NDVI from the 16-day NDVI.'''
    usage = "usage: %prog [mandatory args] [options]"
    version = __version__

    parser = argparse.ArgumentParser()

    # Mandatory arguments

    parser.add_argument('-i','--input_file',
                      action='store',
                      dest='input_file',
                      type=str,
                      required=True,
                      help='''The fully qualified path to the input files. May be 
                              a directory or a file glob.'''
                      )

    # Optional arguments 

    parser.add_argument('--which_ndvi',
                      action="store",
                      dest="ndvi_choice",
                      default=defaults["ndvi_choice"],
                      type=str,
                      choices=ndvi_choices,
                      help='''Which Annual min/max NDVI dataset to plot. 
                              Possible options here are...\n
                              {}.
                           '''.format(ndvi_choices.__str__()[1:-1])
                      )

    parser.add_argument('--plot_ndvi',
                      action="store_true",
                      dest="plot_ndvi",
                      default=defaults["plot_ndvi"],
                      help='''Plot the annual min and max NDVI from a HDF5 
                              file generated previously 
                              [default: {}]'''.format(defaults["plot_ndvi"])
                      )

    parser.add_argument('--dpi',
                      action="store",
                      dest="dpi",
                      default=defaults["dpi"],
                      type=float,
                      help='''An example of an option to set a float variable 
                      [default: {}]'''.format(defaults["dpi"])
                      )

    parser.add_argument('--stride',
                      action="store",
                      dest="stride",
                      default=defaults["stride"],
                      type=int,
                      help='''An example of an option to set a int variable 
                      [default: {}]'''.format(defaults["stride"])
                      )

    parser.add_argument('--output_file',
                      action="store",
                      dest="output_file",
                      default=defaults["output_file"],
                      type=str,
                      help='''The filename of the output annual min/max NDVI HDF5 
                      file [default: {}]'''.format(defaults["output_file"])
                      )

    parser.add_argument("-v", "--verbose",
                      dest='verbosity',
                      action="count", 
                      default=0,
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


def main():
    '''
    The main method.
    '''

    options = _argparse()

    input_file = options.input_file
    plot_ndvi  = options.plot_ndvi
    stride = options.stride
    output_file = options.output_file
    dpi = options.dpi
    ndvi_choice  = options.ndvi_choice


    LOG.info("Input option 'input_file'  = {} ".format(input_file))
    LOG.info("Input option 'plot_ndvi'   = {} ".format(plot_ndvi))
    LOG.info("Input option 'stride'      = {} ".format(stride))
    LOG.info("Input option 'output_file' = {} ".format(output_file))
    LOG.info("Input option 'dpi'         = {} ".format(dpi))
    LOG.info("Input option 'ndvi_choice'   =  {} ".format(ndvi_choice))

    input_dir,input_glob = _create_input_file_globs(input_file)

    if input_glob is None :
        LOG.error("No input files found matching %s, aborting...".format(input_file))
        return 1

    LOG.info("Input directory is {}".format(input_dir))
    LOG.info("Input glob is {}".format(input_glob))

    # Get a list of input files...
    input_files = glob(path.join(input_dir,input_glob))
    input_files.sort()

    for idx in range(len(input_files)):
        input_files[idx] = path.basename(input_files[idx])

    LOG.info(input_files)

    # Plot the NDVI climatology and exit
    if plot_ndvi :
        for files in input_files:
            dset_dicts = _read_ndvi_object_h5py(files)
            
            plotTitle = {'min_ndvi':'Minimum Annual NDVI',
                    'max_ndvi':'Maximum Annual NDVI'}
            png_name = "{}_{}.png".format(string.split(files,'.h5')[:-1][0],
                    ndvi_choice)
            _plot_ndvi(dset_dicts,ndvi_choice,title=plotTitle[ndvi_choice],
                    png_name=png_name,dpi=dpi)

        return 0


    # Do some file operations with h5py

    fileObj = h5py.File(path.join(input_dir,input_files[0]),"r")

    dset_dicts = {}

    # FIXME: Ensure that attributes have correct type (e.g.: _FillValue)
    for dset in ['Latitude','Longitude','NDVI']:
        dset_dicts[dset] = {}
        obj = fileObj[dset]
        try:
            LOG.debug("Checking {} ...".format(dset))
            LOG.debug(obj.attrs.keys())
            LOG.debug(obj.attrs.items())
        except KeyError:
            pass
        for attr_key in obj.attrs.keys():
            dset_dicts[dset][attr_key] = obj.attrs[attr_key]

    fileObj.close()

    # Create some dictionaries

    dset_dicts['min_ndvi'] = {}
    dset_dicts['max_ndvi'] = {}

    for attr_key in dset_dicts['NDVI'].keys():
        dset_dicts['min_ndvi'][attr_key] = dset_dicts['NDVI'][attr_key]
        dset_dicts['max_ndvi'][attr_key] = dset_dicts['NDVI'][attr_key]

    for dset in dset_dicts.keys():
        LOG.debug(dset)
        for attr_key in dset_dicts[dset].keys():
            LOG.debug('\t{} = {}'.format(attr_key,dset_dicts[dset][attr_key]))

    # Create the man and max NDVI arrays

    for ndvi_file in input_files:
        LOG.info('Opening NDVI file {}...'.format(ndvi_file))
        fileObj = h5py.File(path.join(input_dir,ndvi_file),"r")
        ndvi = fileObj['NDVI'].value[::stride,::stride]
        LOG.debug('NDVI shape is = {}'.format(ndvi.shape))

        try:
            LOG.info("Checking current NDVI against previous")
            min_idx = ndvi < prev_ndvi
            max_idx = ndvi > prev_ndvi
            min_ndvi[min_idx] = ndvi[min_idx]
            max_ndvi[max_idx] = ndvi[max_idx]
            prev_ndvi = ndvi
        except Exception, err :
            #LOG.debug(traceback.format_exc())
            LOG.info("Initialsing minimum NDVI")
            min_ndvi = copy.copy(ndvi)
            max_ndvi = copy.copy(ndvi)
            prev_ndvi = ndvi

            Latitude = fileObj['Latitude'].value[::stride]
            Longitude = fileObj['Longitude'].value[::stride]
            LOG.debug('Latitude shape is = {}'.format(Latitude.shape))
            LOG.debug('Longitude shape is = {}'.format(Longitude.shape))

        fileObj.close()

    dset_dicts['Latitude']['data'] = Latitude
    dset_dicts['Longitude']['data'] = Longitude
    dset_dicts['min_ndvi']['data'] = min_ndvi
    dset_dicts['max_ndvi']['data'] = max_ndvi

    _create_ndvi_object_h5py(input_files,dset_dicts,output_file)


    return 0


if __name__=='__main__':
    sys.exit(main())  
