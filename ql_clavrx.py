#!/usr/bin/env python
# encoding: utf-8
"""
ql_clavrx.py
$Id$

Purpose: Create Quicklook PNGs for CLAVR-x Level 2 products

Created by nickb@ssec.wisc.edu Feb 2014.
Copyright (c) 2014 University of Wisconsin SSEC. All rights reserved.
Licensed under GNU GPLv3.
"""

# load libraries 
import matplotlib
matplotlib.use('Agg')
from matplotlib import colors, cm

import numpy as np, glob, os, sys, logging
from netCDF4 import Dataset
from collections import namedtuple
from ql_common_clavrx import *

import string

LOG = logging.getLogger(__name__)

L2Product = namedtuple('L2Product', 'var units longname cmap bg bounds labels')


# Defining some standard cmaps that get used for multiple plots, to keep the table below cleaner
cmap_ug = ['#810541', '#806517', '#C47451', '#F88017', '#FFFF00', '#D4A017', '#00FF00',
           '#347C17', '#54C571', '#99C68E', '#5E5A80', '#6C2DC7', '#8467D7', '#E3E4FA', '#B6B6B4']

# This table holds the basic plot attributes for each of our products

L2_PRODUCTS = {
    "cloud_mask" : L2Product(var="cloud_mask", units="",
                             longname="Cloud Mask",
                             cmap=colors.ListedColormap(['#00FF00', '#00FFFF', '#FF0000', '#FFFFFF']),
                             bg='black',
                             bounds=[0,1,2,3,4],
                             labels=["Clear", "Probably Clear", "Probably Cloudy", "Cloudy"]),

    "cloud_type" : L2Product(var="cloud_type", units="",
                             longname="Cloud Type",
                             cmap=colors.ListedColormap(['#C0C0C0', 'gray', 'green', 'blue', 'cyan', 'pink', 'red', 'yellow', 'orange', 'brown', 'black']),
                             bg='white',
                             bounds=[0,1,2,3,4,5,6,7,8,9,10,11],
                             labels=["Clear", "Prob. Clear", "Fog", "Water", "Supercooled Water", "Mixed", "Opaque Ice", "Cirrus", "Overlapping", "Overshooting", "Unknown"]),

    "cld_temp_acha" : L2Product(var="cld_temp_acha",  units="K",
                                longname='Cloud-top Temperature from AWG Cloud Height Algorithm',
                                cmap=colors.ListedColormap(cmap_ug[::-1]),
                                bg='white',
                                bounds=[295,290,285,280,275,270,265,260,250,240,230,220,210,200,190,180][::-1],
                                labels=None),

    "cld_press_acha" : L2Product(var="cld_press_acha",  units="hPa",
                                 longname='Cloud-top Pressure from AWG Cloud Height Algorithm',
                                 cmap=colors.ListedColormap(cmap_ug[::-1]),
                                 bg='white',
                                 bounds=[1100,900,800,700,650,600,550,500,450,400,350,300,250,200,0][::-1],
                                 labels=None),

    "cld_height_acha" : L2Product(var="cld_height_acha",  units="km",
                                  longname='Cloud-top Height from AWG Cloud Height Algorithm',
                                  cmap=colors.ListedColormap(cmap_ug[::-1]),
                                  bg='white',
                                  bounds=[0,0.5,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,20.0],
                                  labels=None),

    "cld_emiss_acha"    : L2Product(var="cld_emiss_acha",  units="",
                                    longname='Cloud Emissivity from AWG Cloud Height Algorithm',
                                    cmap=cm.jet,
                                    bg='white',
                                    bounds=np.arange(0,1.1,0.1),
                                    labels=None),

#    "cld_opd_acha" : L2Product(var="cld_opd_acha",  units="",
#                               longname='Cloud Optical Depth from AWG Cloud Height Algorithm',
#                               cmap=cm.jet,
#                               bounds=[0,1,2,3,4,5,6,7,8,9,10],
#                               labels=None),
#
#    "cld_reff_acha" : L2Product(var="cld_reff_acha",  units="micron",
#                                longname='Cloud Effective Radius from AWG Cloud Height Algorithm',
#                                cmap=cm.jet,
#                                bounds=[0,4,6,8,10,12,14,16,18,20,24,28,32,36,40],
#                                labels=None),

    "cld_opd_dcomp" : L2Product(var="cld_opd_dcomp",  units="",
                               longname='Cloud Optical Depth from AWG DCOMP Algorithm',
                               cmap=cm.jet,
                               bg='white',
                               bounds=[0,1,2,4,6,10,12,14,18,20,25,30,40,50,60],
                               labels=None),

    "cld_reff_dcomp" : L2Product(var="cld_reff_dcomp",  units="micron",
                                longname='Cloud Effective Radius from AWG DCOMP Algorithm',
                                cmap=cm.jet,
                                bg='white',
                                bounds=[0,4,6,8,10,12,14,16,18,20,24,28,32,36,40],
                                labels=None),

    "cloud_probability" : L2Product(var="cloud_probability",  units="",
                                    longname='Cloud Probability',
                                    cmap=cm.jet,
                                    bg='white',
                                    bounds=np.arange(0,1.1,0.1),
                                    labels=None),

    "cloud_phase" : L2Product(var="cloud_phase",  units="",
                              longname="Cloud Phase",
                              cmap=colors.ListedColormap(['gray', 'blue', 'cyan', 'purple', 'pink', 'black']),
                              bg='white',
                              bounds=[0,1,2,3,4,5,6],
                              labels=['Clear', 'Water', 'Supercooled', 'Mixed', 'Ice', 'Unknown'])
}


DEFAULT_PNG_FMT = 'CLAVRx_%(var)s_%(filedate)s.png'
DEFAULT_LABEL_FMT = 'CSPP CLAVR-x %(longname)s %(date)s %(start_time)s-%(end_time)s'


def clavrx_file_times(clavrx_file):
    stats = {}

    pf = Dataset(clavrx_file, 'r', format="NETCDF3_CLASSIC")
    def fattr(str):
        return getattr(pf, str)

    stats['start_time'] =  datetime.strptime("%s-%s" % (fattr('START_YEAR'), fattr('START_DAY')), "%Y-%j") + timedelta(hours=float(fattr('START_TIME')))
    stats['end_time'] =    datetime.strptime("%s-%s" % (fattr('END_YEAR'), fattr('END_DAY')), "%Y-%j") + timedelta(hours=float(fattr('END_TIME')))

    return stats


MAX_CONTIGUOUS_DELTA = timedelta(seconds=5)
DTZ = timedelta(0)
def are_clavrx_files_contiguous(clavrx_files, groupdicts=None, tolerance=MAX_CONTIGUOUS_DELTA):
    """
    return sequence of booleans saying whether neighboring granules are contiguous or not
    """
    if not groupdicts:
        groupdicts = [clavrx_file_times(f) for f in clavrx_files]
    time_ranges = [(x['start_time'], x['end_time']) for x in groupdicts]
    nfo = [(filename, start, end) for (filename, (start, end)) in zip(clavrx_files, time_ranges)]

    for (na,sa,ea),(nb,sb,eb) in zip(nfo[:-1], nfo[1:]):
        dt = sb - ea
        yield (dt >= DTZ) and (dt <= tolerance)


def clavrx_info(clavrx_files):
    """return a dictionary of information about one or more CLAVR-x files
    takes filenames, but will require opening the files to get this information"""
    nfos = [clavrx_file_times(f) for f in clavrx_files]

    neighbors_are_contig = list(are_clavrx_files_contiguous(clavrx_files, nfos)) + [False]
    swath_breaks_after_filenames = set()
    for flows_into_successor, filename in zip(neighbors_are_contig, clavrx_files):
        if not flows_into_successor:
            swath_breaks_after_filenames.add(filename)
    LOG.debug('swath breaks after %s' % repr(swath_breaks_after_filenames))

    return dict(date=nfos[0]['start_time'].strftime('%Y-%m-%d'),
                start_time=nfos[0]['start_time'].strftime('%H:%M:%S'),
                end_time=nfos[-1]['end_time'].strftime('%H:%M:%S'),
                break_after=swath_breaks_after_filenames)


def clavrx_quicklooks(output_dir, input_paths,
                      break_after_paths=None,
                      png_fmt=DEFAULT_PNG_FMT,
                      label_fmt = DEFAULT_LABEL_FMT,
                      channels = None,
                      dpi=150,
                      swath_lengths=None,
                      products=L2_PRODUCTS.keys(),
                      **kwargs):

    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    #    nfo = info(input_paths)

    # open a directory with a pass of CSPP SaDR files in time order

    h4_pathnames = tuple(input_paths)
    if not break_after_paths:
        break_after_paths = set(h4_pathnames[-1])

    clavrxs = []

    for pathname in h4_pathnames:
        LOG.debug('reading %s' % pathname)
        try:
            pf = Dataset(pathname, 'r', format="NETCDF3_CLASSIC")
        except:
            raise ValueError('Could not open %s - is it a valid HDF4 file?' % pathname)
        clavrxs.append(pf)

    if len(h4_pathnames) == 0:
        LOG.warn("No inputs")
        return None

    infos = clavrx_info(h4_pathnames)

    # load latitude and longitude arrays
    def read_clavrx_variable(f, vname):
        v = f.variables[vname]
        d = v[:]
        return d

    lines_in_swath = 0
    swath_lengths = []
    end_swath_after_these_lines = False
    for f in clavrxs:
        lines_in_swath += len(f.dimensions['scan_lines_along_track_direction'])
        if lines_in_swath and (pathname in break_after_paths):
            swath_lengths.append(lines_in_swath)
            lines_in_swath = 0
    if lines_in_swath:
        swath_lengths.append(lines_in_swath)
    LOG.debug('swath lengths: %r' % swath_lengths)

    lat = np.ma.concatenate([read_clavrx_variable(f, 'latitude')  for f in clavrxs])
    lon = np.ma.concatenate([read_clavrx_variable(f, 'longitude') for f in clavrxs])

    for productname in products:
        if productname not in L2_PRODUCTS.keys():
            LOG.error("%s is not a valid product key for plotting, skipping..." % (productname))
            continue
        product = L2_PRODUCTS[productname]
        LOG.debug("plotting %s" % (product.longname))
        product_data = np.ma.concatenate([read_clavrx_variable(f, product.var) for f in clavrxs])

        pinfo = dict(product._asdict().items() + infos.items())
        pinfo['filedate'] = string.replace(pinfo['date'], '-', '') + "_" + \
                            string.replace(pinfo['start_time'], ':', '') + '-' + \
                            string.replace(pinfo['end_time'], ':', '')
        pngname = os.path.join(output_dir, png_fmt % pinfo)
        label = label_fmt % pinfo

        map_clavrx(pngname, product_data, lon, lat, vmin=min(product.bounds), vmax=max(product.bounds),
                   label = label, dpi=dpi, cmap=product.cmap, bounds=product.bounds, bg=product.bg,
                   units=product.units, labels=product.labels)


"""
A modified version of map_image from ql_common_clavrx
putting it here so all CLAVR-x specific functions stay out of common for now.
"""
def map_clavrx(pngname, swath, lon, lat, size = (DEFAULT_X_SIZE, DEFAULT_Y_SIZE), label=None, dpi=175,
               vmin=None, vmax=None, area_thresh=500000, cmap=None, bg='white', bounds=None, swath_lengths=None, units=None, labels=None):
    """
    map product image to base map
    """
    from matplotlib import colors
    n,m,x = min_median_max(swath)

    if n==x or np.isnan(n) or np.isnan(x):
        LOG.warning('no data in array, min-max check found flat or empty field with minimum %s' % n)
        return None, None

    mpl.rc('font', size=6)

    norm = colors.BoundaryNorm(bounds, cmap.N)
#    cmap.set_over('w',1)
#    cmap.set_under('w',1)
    cmap.set_bad(bg,1)

    fig = plt.figure()

    nmx_lon = min_median_max(lon)
    nmx_lat = min_median_max(lat)
    area_def = optimal_projection(nmx_lon, nmx_lat, size)
    bmap = pr.plot.area_def2basemap(area_def, resolution='i',area_thresh=area_thresh)

#    draw_coastlines_parallels_meridians(bmap, nmx_lon, nmx_lat, bg='black', color=BORDER_COLOR)
    bmap.drawcoastlines(color='grey', linewidth=0.5)
    bmap.drawcountries(color='grey', linewidth=0.5)
    bmap.drawmapboundary(fill_color=bg)
    pstride = np.ceil((nmx_lat.max-nmx_lat.min)/20.0)*5.0
    parallels = np.arange(-90.,90,pstride)
    bmap.drawparallels(parallels, labels=[1,0,0,1], color='grey')
    mstride = np.ceil((nmx_lon.max-nmx_lon.min)/20.0)*5.0
    meridians = np.arange(0.,360.,mstride)
    bmap.drawmeridians(meridians, labels=[1,0,0,1], color='grey')

    if not swath_lengths:
        image = swath2grid(area_def, swath, lon, lat)
        widget = bmap.imshow(image, origin='upper',interpolation='nearest',vmin=vmin,vmax=vmax,cmap=cmap, norm=norm)
    else:
        start=0
        for lines in swath_lengths:
            stop = start + lines
            image = swath2grid(area_def, swath[start:stop,:], lon[start:stop,:], lat[start:stop,:])
            widget = bmap.imshow(image, origin='upper', interpolation='nearest', vmin=vmin, vmax=vmax, cmap=cmap, norm=norm)
            start += lines

    cbar = plt.colorbar(label=units)

    if labels is None:
        labels = bounds
        ticks = bounds
    else:
        # if we have labels, put ticks in the center of each bound range
        ticks = [ ((bounds[n] + bounds[n-1]) / 2.0) for n in range(1, len(bounds)) ]

    cbar.set_ticks(ticks)
    cbar.set_ticklabels(labels)

    if label:
        plt.title(label)

    fig.savefig(pngname, dpi=dpi, format="png")


def main():
    import optparse
    usage = """
%prog [options] ...
This program creates PNG files of instrument quick-look data.

If given a directory instead of filenames, it will find all input files in the directory
and order them by time.

If a series of directories are listed, all swaths (each swath represented as a directory)
will be placed on a single plot.

The output directory will be created if it does not exist.

Example: 
%prog -o /tmp/clavrx-quicklooks /path/to/cspp-output /path/to/cspp-output2 /path/to/cspp-output3



"""
    parser = optparse.OptionParser(usage)
    parser.add_option('-t', '--test', dest="self_test",
                      action="store_true", default=False, help="run self-tests")
    parser.add_option('-v', '--verbose', dest='verbosity', action="count", default=0,
                      help='each occurrence increases verbosity 1 level through ERROR-WARNING-INFO-DEBUG')
    parser.add_option('-D', '--dpi', dest='dpi', type='int', default=175,
                      help='dots per inch of plots to produce')
    parser.add_option('-o', '--output', dest='output',  default='.',
                      help='directory in which to store output')
    parser.add_option('-F', '--format', dest='format', default=DEFAULT_PNG_FMT,
                      help='format string for output filenames')
    parser.add_option('-L', '--label', dest='label', default=DEFAULT_LABEL_FMT,
                      help='format string for labels')
    parser.add_option('-p', '--products', dest='products', default=",".join(L2_PRODUCTS.keys()),
                      help='comma-separated list of products to plot')

    (options, args) = parser.parse_args()

    # FUTURE: validating the format strings is advisable

    # make options a globally accessible structure, e.g. OPTS.
    global OPTS
    OPTS = options

    if options.self_test:
        # FIXME - run any self-tests
        # import doctest
        # doctest.testmod()
        sys.exit(2)

    levels = [logging.ERROR, logging.WARN, logging.INFO, logging.DEBUG]
    logging.basicConfig(level = levels[min(3,options.verbosity)])

    if not args:
        parser.error( 'incorrect arguments, try -h or --help.' )
        return 9

    # convert list of files and/or directories to a flat sequence of pathnames
    pathnames = list(expand_paths(args, 'CLAVR*.hdf'))
    if not pathnames:
        return 1
        # collect information about how we should label these and break them into swaths
    #    nfo = info(*pathnames)

    # load the swath data by actually reading the files, noting number of scanlines in each swath
    clavrx_quicklooks(options.output, pathnames, png_fmt=options.format, label_fmt=options.label,
                      dpi=options.dpi, products=options.products.split(','))

    return 0



if __name__=='__main__':
    sys.exit(main())
    
