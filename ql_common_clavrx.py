#!/usr/bin/env python
# encoding: utf-8
"""
ql_common.py
$Id$

Purpose: Amend ADL HDF5 SDR products with N_GEO_Ref attributes,
based on guide-book information and ASC metadata track-back.

Created by rayg@ssec.wisc.edu Jan 2012.
Copyright (c) 2011 University of Wisconsin SSEC. All rights reserved.
Licensed under GNU GPLv3.
"""


import os, sys, logging
import numpy as np
import pyproj
import h5py 
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.image as mpimg
import matplotlib.cm as cm
from collections import namedtuple, OrderedDict as ordict
import pyresample as pr
from datetime import timedelta, datetime
from scipy import misc
from glob import glob


# in scripts/common/
from adl_geo_ref import RE_NPP

LOG = logging.getLogger(__name__)

EXTENTS_FACTOR = 0.55
DEFAULT_X_SIZE = 600
DEFAULT_Y_SIZE = 600

class evaluator(object):
    """
    evaluate expressions in format statements
    e.g. '%(4*2)d %(", ".join([c]*b))s' % evaluator(a=2,b=3,c='foo')
    """
    def __init__(self,**argd):
        vars(self).update(argd)
    def __getitem__(self, expr):
        return eval(expr, globals(), vars(self))

# monkeypatch h5py File objects to have a .var('/path/to/variable') routine
h5py.File.var = lambda h5,pth: reduce( lambda x,a: x[a] if a else x, pth.split('/'), h5)


MAX_CONTIGUOUS_DELTA_IN_FILENAME_TIMESTAMP=timedelta(seconds=5)


def _startend(date, start_time, end_time, **kwargs):
     s = datetime.strptime('%sT%s' % (date, start_time[:6]), '%Y%m%dT%H%M%S')
     e = datetime.strptime('%sT%s' % (date, end_time[:6]), '%Y%m%dT%H%M%S')
     if (e < s):
         e += timedelta(hours=24)
     return s,e


DTZ = timedelta(0)


def are_cdfcb_filenames_contiguous(filenames, groupdicts=None, tolerance=MAX_CONTIGUOUS_DELTA_IN_FILENAME_TIMESTAMP):
    """
    return sequence of booleans saying whether neighboring granules are contiguous or not
    """
    if not groupdicts:
        groupdicts = [RE_NPP.match(os.path.split(x)[-1]).groupdict() for x in filenames]
    time_ranges = [_startend(x['date'], x['start_time'], x['end_time']) for x in groupdicts]
    nfo = [(filename, start, end) for (filename, (start, end)) in zip(filenames, time_ranges)]

    for (na,sa,ea),(nb,sb,eb) in zip(nfo[:-1], nfo[1:]):
        dt = sb - ea
        yield (dt >= DTZ) and (dt <= tolerance)


def info(*sdr_filenames):
    "return a dictionary of information about one or more files; requires valid NPP filenames"
    nfos = [RE_NPP.match(os.path.split(fn)[-1]).groupdict() for fn in sdr_filenames]

    neighbors_are_contig = list(are_cdfcb_filenames_contiguous(sdr_filenames, nfos)) + [False]
    swath_breaks_after_filenames = set()
    for flows_into_successor, filename in zip(neighbors_are_contig, sdr_filenames):
        if not flows_into_successor:
            swath_breaks_after_filenames.add(filename)
    LOG.debug('swath breaks after %s' % repr(swath_breaks_after_filenames))

    return dict( date=nfos[0]['date'],
                 start_time = nfos[0]['start_time'],
                 end_time = nfos[-1]['end_time'],
                 orbit = nfos[0]['orbit'],
                 site = nfos[0]['site'],
                 domain = nfos[0]['domain'],
                 kind  = nfos[0]['kind'],
                 sat = nfos[0]['sat'],
                 break_after=swath_breaks_after_filenames)


def piecewise_scaled(data, scale, dtype=None):
    """given a scale variable that consists of one or more m,b pairs
    and a data variable that contains 1 or more scan-lines
    apply proper scaling parameter to blocks of scan-lines
    NPP stores individual granules with m,b scales; however
    nagg just concatenates variables together.
    """
    scale = scale.flatten()
    lines = data.shape[0]
    if dtype is None:
        dtype = data.dtype
    else:
        data = data.astype(dtype)
    group_count = scale.shape[0] / 2
    group = lines / group_count
    LOG.debug('group size for %r is %d lines' % (data.shape, group))
    zult = np.empty(data.shape, dtype=dtype)
    for n in range(group_count):
        m,b = scale[n*2:n*2+2]
        LOG.debug('scale group %d: %fx + %f' % (n,m,b))
        zult[n*group:(n+1)*group, :] = data[n*group:(n+1)*group, :] * m + b
    return zult


def expand_paths(paths, glob_pattern):
    "iterate through a series of paths, expanding directories using a glob pattern"
    for path in paths:
        if os.path.isdir(path):
            for fn in sorted(glob(os.path.join(path, glob_pattern))):
                yield fn
        elif os.path.exists(path):
            yield path
        else:
            raise ValueError('path %s goes does not refer to a directory or file' % path)


#
# Projection Math
#

# 
# FIXME: prove that this works across the dateline
# 

nmx = namedtuple('nmx', 'min median max')

def min_median_max(x):
    xf = x.flatten()
    xk = np.isfinite(xf)
    return nmx(np.min(xf[xk]), np.median(xf[xk]), np.max(xf[xk]))
    

# see http://www.remotesensing.org/geotiff/proj_list/
# http://matplotlib.github.com/basemap/users/mapsetup.html

def _tuple2args(parms):
    s = ' '.join( '+%s=%s' % (k,v) for (k,v) in parms )
    return s.encode('ascii')

def _polar_area(nmx_lon, nmx_lat, x_size, y_size):
    """
      +proj=stere +lat_ts=Latitude at natural origin 
              +lat_0=90
              +lon_0=Longitude at natural origin
              +k_0=Scale factor at natural origin (normally 1.0)
              +x_0=False Easting
              +y_0=False Northing
    """
    polar_lat = -90 if nmx_lat.median < 0 else 90
    parms = (   ('proj', 'stere'),
                ('lat_0', polar_lat),
                ('lat_ts', nmx_lat.median),
                ('lon_0', nmx_lon.median),
                ('ellps', 'WGS84'),
                ('units', 'm')
            )
    proj4_args = _tuple2args(parms)
    LOG.debug(proj4_args)
    
    # compute extents
    proj = pyproj.Proj(proj4_args)
    x0,y0 = proj(0, polar_lat)
    x1,y1 = proj(nmx_lon.min, nmx_lat.min)
    x2,y2 = proj(nmx_lon.max, nmx_lat.max)
    dx = max([(x2-x0) / EXTENTS_FACTOR, (x1-x0) / EXTENTS_FACTOR])
    dy = max([(y2-y0) / EXTENTS_FACTOR, (y1-y0) / EXTENTS_FACTOR])
    md = max(dx, dy)
    area_extent = (-md, -md, md, md)

    area_id = proj_id = 'polargrid'
    area_name = 'Polar Stereographic Grid'
    return pr.utils.get_area_def(area_id, area_name, proj_id, proj4_args, x_size, y_size, area_extent)

def _midlatitude_area(nmx_lon, nmx_lat, x_size, y_size):
    """
     +proj=lcc   +lat_1=Latitude of natural origin
             +lon_0=Longitude of natural origin
             +k_0=Scale factor at natural origin
             +x_0=False Origin Easting
             +y_0=False Origin Northing
    """
    dlat = max( [nmx_lat.max - nmx_lat.median, nmx_lat.median - nmx_lat.min] ) / 2.0
    lat_0 = nmx_lat.median
    parms = (   ('proj', 'stere'),
                ('lat_0', lat_0),
#                ('lat_1', lat_0 - dlat),
#                ('lat_2', lat_0 + dlat),
                ('lat_ts', lat_0),
                ('lon_0', nmx_lon.median),
                ('ellps', 'WGS84'),
                ('units', 'm')
            )
    proj4_args = _tuple2args(parms)
    LOG.debug(proj4_args)

    # compute extents
    proj = pyproj.Proj(proj4_args)
    x1,y1 = proj(nmx_lon.min, nmx_lat.min)
    x2,y2 = proj(nmx_lon.max, nmx_lat.max)
    dx = (x2-x1) * EXTENTS_FACTOR
    dy = (y2-y1) * EXTENTS_FACTOR
    md = max(dx, dy)
    area_extent = (-md, -md, md, md)

    area_id = proj_id = 'stereogrid'
    area_name = 'Stereoscopic Grid'
    return pr.utils.get_area_def(area_id, area_name, proj_id, proj4_args, x_size, y_size, area_extent)


_equatorial_area = _midlatitude_area

def _global_area(nmx_lon, nmx_lat, x_size, y_size):
    " returns area def for cylindrical projection (cyl)"
    parms = (   ('proj', 'eqc'),
            )
    proj4_args = _tuple2args(parms)
    LOG.debug(proj4_args)

    # compute extents (entire globe)
    proj = pyproj.Proj(proj4_args)
    lon_bbox = (-180, 180)
    lat_bbox = (-90, 90)
    x, y = proj(lon_bbox, lat_bbox)
    area_extent = (x[0], y[0], x[1], y[1])

    # get area def
    area_id = proj_id = 'pc_world'
    area_name = 'Plate Carree world map'
    return pr.utils.get_area_def(area_id, area_name, proj_id, proj4_args, x_size, y_size, area_extent)
    
# def _equatorial_area(nmx_lon, nmx_lat, x_size, y_size):
#     """
#      +proj=lcc   +lat_1=Latitude of natural origin
#              +lon_0=Longitude of natural origin
#              +k_0=Scale factor at natural origin
#              +x_0=False Origin Easting
#              +y_0=False Origin Northing
#     """
#     dy = max( [nmx_lat.max - nmx_lat.median, nmx_lat.median - nmx_lat.min] ) / 2.0
#     lat_0 = nmx_lat.median
#     parms = (   ('proj', 'lcc'),
#                 ('lat_0', lat_0),
#                 ('lat_1', lat_0 - dy),
#                 ('lat_2', lat_0 + dy),
#                 ('lon_0', nmx_lon.median),
#                 ('ellps', 'WGS84'),
#                 ('units', 'm')
#             )
#     proj4_args = _tuple2args(parms)
#     LOG.debug(proj4_args)
#     
#     # compute extents
#     proj = pyproj.Proj(proj4_args)
#     x1,y1 = proj(nmx_lon.min, nmx_lat.min)
#     x2,y2 = proj(nmx_lon.max, nmx_lat.max)
#     dx = (x2-x1) * EXTENTS_FACTOR
#     dy = (y2-y1) * EXTENTS_FACTOR
#     md = max(dx, dy)
#     area_extent = (-md, -md, md, md)
# 
#     area_id = proj_id = 'lccgrid'
#     area_name = 'Equatorial Grid'
#     return pr.utils.get_area_def(area_id, area_name, proj_id, proj4_args, x_size, y_size, area_extent)


def optimal_projection(nmx_lon, nmx_lat, size = (DEFAULT_X_SIZE, DEFAULT_Y_SIZE)):
    "create an area definition suitable for a given set of latitude and longitude arrays"
    x_size, y_size = size
    zone = abs(nmx_lat.median)
    lat_range = nmx_lat.max - nmx_lat.min
    LOG.debug('nmx_lat = %r nmx_lon = %r' % (nmx_lat, nmx_lon))
    if zone > 66.0:
        return _polar_area(nmx_lon, nmx_lat, x_size, y_size)
    if zone > 23.0:
        return _midlatitude_area(nmx_lon, nmx_lat, x_size, y_size)
    if lat_range > 120.0:
        return _global_area(nmx_lon, nmx_lat, x_size, y_size)
    return _equatorial_area(nmx_lon, nmx_lat, x_size, y_size)


# def centered_projection( center_lat, center_lon, radius, size = (DEFAULT_X_SIZE, DEFAULT_Y_SIZE) ):
#     "for a center latitude and longitude and radius in meters, return an area definition"
#     x_size, y_size = size
#     
#     if abs(lat) > 66.0: 
#         # build polar stereographic
#         
#     else:
#         # build lcc
#         dy = max( [nmx_lat.max - nmx_lat.median, nmx_lat.median - nmx_lat.min] ) / 2.0
#         lat_0 = center_lat
#         proj = pyproj.Proj('+proj=lcc +ellps=WGS84 +units=m')
#         x0,y0 = proj(lon, lat)
#         xlat, xlon = proj(x0+radius, y0+radius, inverse=True)
#         
#         parms = (   ('proj', 'lcc'),
#                     ('lat_0', lat_0),
#                     ('lat_1', lat_0 - dy),
#                     ('lat_2', lat_0 + dy),
#                     ('lon_0', nmx_lon.median),
#                     ('ellps', 'WGS84'),
#                     ('units', 'm')
#                 )
#         proj4_args = _tuple2args(parms)
#         LOG.debug(proj4_args)
#     
#         # compute extents
#         proj = pyproj.Proj(proj4_args)
#         x1,y1 = proj(nmx_lon.min, nmx_lat.min)
#         x2,y2 = proj(nmx_lon.max, nmx_lat.max)
#         dx = (x2-x1) * EXTENTS_FACTOR
#         dy = (y2-y1) * EXTENTS_FACTOR
#         md = max(dx, dy)
#         area_extent = (-md, -md, md, md)
#     
#         area_id = proj_id = 'lccgrid'
#         area_name = 'Midlatitude LCC Grid'
#         # area_extent = (-4000000,-3000000,4000000,3000000) # FIXME (xll, yll, xur, yur)
#         return pr.utils.get_area_def(area_id, area_name, proj_id, proj4_args, x_size, y_size, area_extent)

        
    
    

def _test_area():
    from pyresample import utils
    area_id = 'mygrid'
    area_name = 'My Grid'
    proj_id = 'mygrid'
    proj4_args = '+proj=laea +lat_0=45 +lon_0=-90 +a=6371228.0 +units=m'
    x_size = 2048
    y_size = 2048
    area_extent = (-4000000,-3000000,4000000,3000000)
    area_def = pr.utils.get_area_def(area_id, area_name, proj_id, proj4_args,
                                     x_size, y_size, area_extent)
    return area_def

def swath2grid(area, swath, lon, lat):
    "return gridded image using pyresample, along with basemap object"
# fixme - make sure we still handle non-masked arrays
    try:
        nodata_mask = np.ma.getmaskarray(swath)
        lon = np.ma.masked_array(lon, nodata_mask)
        lat = np.ma.masked_array(lat, nodata_mask)
    except:
        pass

    swath_def = pr.geometry.SwathDefinition(lons=lon, lats=lat)
    image = pr.kd_tree.resample_nearest(swath_def, swath, area, nprocs=2,
                                        radius_of_influence=20000, fill_value=None)
    return image
    

def plot_gridded(area, image, **kwargs):
    bmap = pr.plot.area_def2basemap(area)
    # bmng = bmap.bluemarble()
    col = bmap.imshow(image, origin='upper')
    return bmap


def plot(swath, lon, lat, size = (DEFAULT_X_SIZE, DEFAULT_Y_SIZE)):
    nmx_lon = min_median_max(lon)
    nmx_lat = min_median_max(lat)
    area_def = optimal_projection(lon, lat, size)
    image = swath2grid(area_def, swath, lon, lat)
    return plot_gridded(area_def, image)
    

#BORDER_COLOR = 'cadetblue' # '#20ff40'
BORDER_COLOR = 'grey' # '#666666'


def draw_coastlines_parallels_meridians(m, nmx_lon, nmx_lat, bg='white', **kwargs):
    m.drawcoastlines(color=BORDER_COLOR, linewidth=0.5)
    m.drawmapboundary(fill_color=bg)
    m.drawcountries(color=BORDER_COLOR, linewidth=0.5)
#    m.drawstates(color=BORDER_COLOR, linewidth=0.5)
    pstride = np.ceil((nmx_lat.max-nmx_lat.min)/20.0)*5.0
    parallels = np.arange(-90.,90,pstride)
    m.drawparallels(parallels, labels=[1,0,0,1], **kwargs)
    mstride = np.ceil((nmx_lon.max-nmx_lon.min)/20.0)*5.0
    meridians = np.arange(0.,360.,mstride)
    m.drawmeridians(meridians, labels=[1,0,0,1], **kwargs)


def pre_scale( swath , vmin=None,vmax=None) :
    n,m,x = min_median_max(swath.flatten())
    LOG.debug(" max: %f median: %f min: %f" % (x,m,n))

    max = float ( vmax or x )
    min = float ( vmin or n )
    range = max - min
    
    np.subtract(swath,min,out=swath)
    np.divide(swath, range, out=swath)
    np.sqrt(swath, out=swath)
    np.multiply(swath, (range), out=swath)
    np.add(swath,min, out=swath)

    n,m,x = min_median_max(swath.flatten())
    LOG.debug(" max: %f median: %f min: %f" % (x,m,n))

    return swath

def _is_empty(arr):
    if np.prod(arr.shape) == 0:
        LOG.info('array has no dimension')
        return True
    if not hasattr(arr, 'mask'):
        LOG.warning('masked arrays should be used for quicklooks')
        n,m,x = min_median_max(arr)
        # look for "standard-until-it's-not" missing value, FIXME this is a bad assumption
        if n == x == -999.0:
            return True
        # test finite-ness, i.e. look for nans
        if np.sum(np.isfinite(arr).flatten()) == len(arr.flatten()):
            return True
        return False
    return np.sum(arr.mask.flatten()) == len(arr.mask.flatten())



def raw_image(pngname, swath, label=None, dpi=175, scale='black2white',sqrt_enhance=False,vmin=None,vmax=None):
    """ map VIIRS image to base map
        generate raw image and remapped image
    """
  
#    if sqrt_enhance == True :
#	LOG.debug("Enhance")
#	swath = pre_scale(swath)

    n,m,x = min_median_max(swath)
    pngraw = pngname.replace(".png",".raw.png" ) 

    if _is_empty(swath):
        LOG.warn('No data for "%s"' % pngraw)
        return 1

    if	scale=='black2white' :
        cmap=cm.get_cmap('gray')
    else :
        LOG.debug("Reverse scale")
        cmap=cm.get_cmap('gray_r')

    mpl.rc('font', size = 6)

    fig = plt.figure()
    
    plt.imshow(swath,vmin=vmin,vmax=vmax,cmap=cmap)
    if label:
        plt.title(label)

    fig.savefig(pngraw, dpi=dpi, format="png")
    
    return fig



def map_image(pngname, swath, lon, lat, size = (DEFAULT_X_SIZE, DEFAULT_Y_SIZE), label=None, dpi=175, scale='black2white',
              sqrt_enhance=False, vmin=None, vmax=None, area_thresh=500000, swath_lengths=None, cmap=None):
    """ map VIIRS image to base map
        generate raw image and remapped image
    """

    if _is_empty(swath):
        LOG.warn('No data for "%s"' % pngname)
        return None, None

    n,m,x = min_median_max(swath)
    #LOG.debug(" max: "+str((swath.data).max())+" min: "+str((swath.data).min()))

    if n==x or np.isnan(n) or np.isnan(x):
        LOG.warning('no data in array, min-max check found flat or empty field with minimum %s' % n)
        return None, None

    if sqrt_enhance:
        LOG.debug("Enhance")
        swath = pre_scale(swath,vmax=vmax,vmin=vmin)

    if cmap is not None:
        pass
    elif scale == 'black2white' :
        cmap = cm.get_cmap('gray')
    else:
        LOG.debug("Reverse scale")
        cmap = cm.get_cmap('gray_r')

    mpl.rc('font', size=6)

    #fig = plt.figure()
    #pngraw = pngname.replace(".png",".raw.png" ) 
    #
    #plt.imshow(swath,vmin=vmin,vmax=vmax,cmap=cmap)
    #if label:
    #    plt.title(label)
    #
    #fig.savefig(pngraw, dpi=dpi, format="png")
    fig = plt.figure()
    
    nmx_lon = min_median_max(lon)
    nmx_lat = min_median_max(lat)
    area_def = optimal_projection(nmx_lon, nmx_lat, size)
    bmap = pr.plot.area_def2basemap(area_def, resolution='i',area_thresh=area_thresh)

    if not swath_lengths:
        image = swath2grid(area_def, swath, lon, lat)
        widget = bmap.imshow(image, origin='upper',interpolation='none',vmin=vmin,vmax=vmax,cmap=cmap)
    else:
        start=0
        for lines in swath_lengths:
            stop = start + lines
            image = swath2grid(area_def, swath[start:stop,:], lon[start:stop,:], lat[start:stop,:])
            widget = bmap.imshow(image, origin='upper', interpolation='none', vmin=vmin, vmax=vmax, cmap=cmap)
            start += lines

    draw_coastlines_parallels_meridians(bmap, nmx_lon, nmx_lat, color=BORDER_COLOR)
    if label:
        plt.title(label)
    if pngname and not sqrt_enhance:
        plt.colorbar()

    fig.savefig(pngname, dpi=dpi, format="png")
    
    return fig, bmap


def quicklook(pngname, swath, lon, lat, size = (DEFAULT_X_SIZE, DEFAULT_Y_SIZE), label=None):
    nmx_lon = min_median_max(lon)
    nmx_lat = min_median_max(lat)
    area_def = optimal_projection(nmx_lon, nmx_lat, size)
    image = swath2grid(area_def, swath, lon, lat)    
    pr.plot.save_quicklook(pngname, area_def, image, label=label, backend='Agg',
                            coast_res = 'i', num_meridians=10, num_parallels=10)

def map_contours(pngname, swath, lon, lat, size = (DEFAULT_X_SIZE, DEFAULT_Y_SIZE), label=None, dpi=175,
                 swath_lengths=None, font_size=None, vmin=None, vmax=None):

    return plot_data(pngname, swath, lon, lat, size=size, label=label, dpi=dpi,
                 swath_lengths=swath_lengths, font_size=font_size, vmin=vmin, vmax=vmax, plot_type='contour')
    
def map_scatter(pngname, swath, lon, lat, size = (DEFAULT_X_SIZE, DEFAULT_Y_SIZE), label=None, dpi=175,
                 swath_lengths=None, font_size=None, vmin=None, vmax=None):

    return plot_data(pngname, swath, lon, lat, size=size, label=label, dpi=dpi,
                 swath_lengths=swath_lengths, font_size=font_size, vmin=vmin, vmax=vmax, plot_type='scatter')
    

    
def plot_data(pngname, swath, lon, lat, size = (DEFAULT_X_SIZE, DEFAULT_Y_SIZE), label=None, dpi=175,
                 swath_lengths=None, font_size=None, vmin=None, vmax=None, plot_type=None):
    """plots with the specified plot type ('contour' or 'scatter'), or if not specified chooses the best for
    the dataset.
    """
    
#   Take missing data out of the picture    
    try:
        nodata_mask = np.ma.getmaskarray(swath)
        lon = np.ma.masked_array(lon, nodata_mask)
        lat = np.ma.masked_array(lat, nodata_mask)
    except:
        pass

    fig = plt.figure()
    nmx_lon = min_median_max(lon)
    nmx_lat = min_median_max(lat)
    if not vmin:
        vmin = np.nanmin(swath)
    if not vmax:
        vmax = np.nanmax(swath)
    area_def = optimal_projection(nmx_lon, nmx_lat, size)
    bmap = pr.plot.area_def2basemap(area_def, resolution='i')

    # If plot type not specified, choose best type for the projection. Avoid these combinations:
    # 'cyl' projection with contour (observed smearing in polar region in global image), and
    # smaller datasets with scatter (dots become sparse).
    if not plot_type:
        if bmap.projection == 'cyl':
            plot_type = 'scatter'
        else:
            plot_type = 'contour'
        LOG.debug("Chose plot type %s based on projection %s" % (plot_type, bmap.projection))
        
    x,y = bmap(lon, lat)
    if swath_lengths is None:
        if plot_type == 'scatter':
            bmap.scatter(x.flatten(), y.flatten(), c=swath.squeeze().flatten(), marker=',', faceted=False, s=2,
                         vmin=vmin, vmax=vmax)
        elif plot_type == 'contour':
            bmap.contourf(x, y, swath.squeeze(), 100, vmin=vmin, vmax=vmax)
        else:
            raise ValueError('Plot type not supported: %s' % plot_type)
        
    else:
        line0 = 0
        for lines in swath_lengths:
            LOG.debug("plotting swath of %d lines" % lines)
            if plot_type == 'scatter':
                bmap.scatter(x[line0:line0+lines, :].flatten(), y[line0:line0+lines, :].flatten(),
                             c=swath.squeeze()[line0:line0+lines, :].flatten(), marker=',', faceted=False, s=2,
                             vmin=vmin, vmax=vmax)
            elif plot_type == 'contour':
                bmap.contourf(x[line0:line0+lines, :], y[line0:line0+lines, :], swath.squeeze()[line0:line0+lines, :], 100, vmin=vmin, vmax=vmax)
            else:
                raise ValueError('Plot type not supported: %s' % plot_type)
            line0 += lines
            
    draw_coastlines_parallels_meridians(bmap, nmx_lon, nmx_lat, color=BORDER_COLOR)
    if label:
        if font_size:
            plt.title(label, size=font_size)
        else:
            plt.title(label)
    if pngname:
        fig.savefig(pngname, dpi=dpi)
    return fig, bmap
    
