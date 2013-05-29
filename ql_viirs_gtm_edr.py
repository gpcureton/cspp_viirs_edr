#!/usr/bin/env python
# encoding: utf-8
"""
ql_viirs_sdr.py
$Id$

Purpose: Generate quicklook images of VIIRS SDR data

Copyright (c) 2012 University of Wisconsin SSEC. All rights reserved.
"""

import matplotlib
matplotlib.use('Agg')

import os, sys, logging
from collections import namedtuple
from glob import glob
import numpy as np
import re
import pyresample as pr
#from datetime import datetime, timedelta

LOG = logging.getLogger(__name__)

import h5py
from ql_common import *

LOG = logging.getLogger(__name__)

h5py.File.var = lambda h5,pth: reduce(lambda x,a: x[a] if a else x, pth.split('/'), h5)

GtmSwath = namedtuple('Swath', 'lat lon data raw paths info')


VIIRS_GTM_EDR_GUIDE = [('VM01O',
                        'GMGTO',
                        'VIIRS-M1ST-EDR',
                        'BrightnessTemperatureOrReflectance'),
                       ('VM02O',
                        'GMGTO',
                        'VIIRS-M2ND-EDR',
                        'BrightnessTemperatureOrReflectance'
                       ),
                       ('VM03O',
                        'GMGTO',
                        'VIIRS-M3RD-EDR',
                        'BrightnessTemperatureOrReflectance'
                       ),
                       ('VM04O',
                        'GMGTO',
                        'VIIRS-M4TH-EDR',
                        'BrightnessTemperatureOrReflectance'
                       ),
                       ('VM05O',
                        'GMGTO',
                        'VIIRS-M5TH-EDR',
                        'BrightnessTemperatureOrReflectance'
                       ),
                       ('VM06O',
                        'GMGTO',
                        'VIIRS-M6TH-EDR',
                        'BrightnessTemperatureOrReflectance'
                       ),
                       ('GMGTO',
                        None,
                        'VIIRS-MOD-GTM-EDR-GEO',
                        None),
                       ('VI1BO',
                        'GIGTO',
                        'VIIRS-I1-IMG-EDR',
                        'Radiance'),
                       ('VI2BO',
                        'GIGTO',
                        'VIIRS-I2-IMG-EDR',
                        'Radiance'
                       ),
                       ('VI3BO',
                        'GIGTO',
                        'VIIRS-I3-IMG-EDR',
                        'Radiance'
                       ),
                       ('VI4BO',
                        'GIGTO',
                        'VIIRS-I4-IMG-EDR',
                        'BrightnessTemperature'
                       ),
                       ('VI5BO',
                        'GIGTO',
                        'VIIRS-I5-IMG-EDR',
                        'BrightnessTemperature'
                       ),
                       ('GIGTO',
                        None,
                        'VIIRS-IMG-GTM-EDR-GEO',
                        None),
                       ('VNCCO',
                        'GNCCO',
                        'VIIRS-NCC-EDR',
                        'Albedo'),
                       ('GNCCO',
                        None,
                        'VIIRS-NCC-EDR-GEO',
                        None)]

# modified regular expression grabs information from EDR files
RE_NPP = re.compile('(?P<kind>[A-Z0-9]+)(?P<band>[0-9]*)_(?P<sat>[A-Za-z0-9]+)_d(?P<date>\d+)'
                    '_t(?P<start_time>\d+)_e(?P<end_time>\d+)_b(?P<orbit>\d+)_c(?P<created_time>\d+)'
                    '_(?P<site>[a-zA-Z0-9]+)_(?P<domain>[a-zA-Z0-9]+)\.h5')

def _variables(h5):
    """
    seek out the data content for a given product file

    :param h5: hdf5 object
    :return: (data-variable-path, factors-variable-path, qf-variable-path-or-None) for this file
    """
    fn = os.path.split(h5.filename)[-1]
    # grab information about science
    nfo, = [x for x in VIIRS_GTM_EDR_GUIDE if fn.startswith(x[0])]
    pfx, gpfx, cn, vn = nfo
    # grab information about geo
    gnfo, = [x for x in VIIRS_GTM_EDR_GUIDE if x[0] == nfo[1]]
    _, _, gcn, _ = gnfo
    fvn = vn.replace('Temperature', '')  # BrightnessTemperature -> BrightnessFactors
    zult = ('/All_Data/%(cn)s_All/%(vn)s' % locals(),
            '/All_Data/%(cn)s_All/%(fvn)sFactors' % locals(),
            '/All_Data/%(gcn)s_All/' % locals(),
            '/All_Data/VIIRS-NCC-EDR_All/QF1_VIIRSNCCEDR' if fn.startswith('VNCCO') else None)
    LOG.debug('%s => %s' % (fn, repr(zult)))
    return zult


# def _startend(date, start_time, end_time, **kwargs):
#     s = datetime.strptime('%sT%s' % (date, start_time[:6]), '%Y%m%dT%H%M%S')
#     e = datetime.strptime('%sT%s' % (date, end_time[:6]), '%Y%m%dT%H%M%S')
#     if (e < s):
#         e += timedelta(hours=24)
#     return s,e


def _info(h5s):
    """
    return dictionary of info given a series of hdf5 objects
    :param h5s: hdf5 object sequence
    :return: dictionary of strings
    """
    s = RE_NPP.match(os.path.split(h5s[0].filename)[-1]).groupdict()
    e = RE_NPP.match(os.path.split(h5s[-1].filename)[-1]).groupdict()
    nfo = dict(s)
    # ss,se = _startend(s)
    # es,ee = _startend(e)
    # nfo['start_time'] = ss
    nfo['end_time'] = e['end_time']
    return nfo


def gtm_swath(*edr_filenames, **kwargs):
    """Load a swath from a series of input files.
    If given a directory name, will load all files in the directory and sort them into lex order.
    returns GtmSwath named_tuple
    """
    # open a directory with a pass of CSPP SDR files in time order
    LOG.debug('loading from %s' % repr(edr_filenames))
    if len(edr_filenames)==1 and os.path.isdir(edr_filenames[0]):
        edr_filenames = glob.glob(os.path.join(edr_filenames[0], 'VNCCO*'))
    edr_filenames = list(sorted(edr_filenames))

    if len(edr_filenames) == 0:
        LOG.warn("No inputs")
        return None

    edrs= [h5py.File(filename, 'r') for filename in edr_filenames]
    data_var, factor_var, geo_var_pfx, qf_var = _variables(edrs[0])
    info = _info(edrs)

    # read all unscaled BTs, and their scaling slope and intercept
    unscaled = [f.var(data_var)[:] for f in edrs]
    scale_factors = [f.var(factor_var)[:] for f in edrs]

    # FUTURE: handle masking off missing values using np.ma.masked_array

    # scale them and concatenate into a contiguous array
    scaled = [piecewise_scaled(t, s, dtype=np.float32) for (t,s) in zip(unscaled,scale_factors)]
    missing = [(x >= 65528) for x in unscaled] # from CDFCB-X Vol3 p12, missing data sentinels
    mask = np.concatenate(missing)
    data = np.concatenate(scaled)
    data[mask] = np.nan
    # FIXME : if QF is provided, use that (VNCCO)

    # load latitude and longitude arrays
    def geo_filename(pn, hp):
        dirname = os.path.split(pn)[0]
        LOG.debug('reading N_GEO_Ref from %s' % pn)
        return os.path.join(dirname, hp.attrs['N_GEO_Ref'][0][0])
    geo_filenames = [geo_filename(pn, hp) for pn,hp in zip(edr_filenames, edrs)]
    geos = [h5py.File(filename, 'r') for filename in list(geo_filenames)]

    lat_var = geo_var_pfx + 'Latitude'
    lon_var = geo_var_pfx + 'Longitude'

    lat = np.concatenate([f.var(lat_var)[:] for f in geos])
    lon = np.concatenate([f.var(lon_var)[:] for f in geos])

    lat[lat <= -999] = np.nan
    lon[lon <= -999] = np.nan
    LOG.debug(str(type(data)))
    LOG.debug(np.nanmax(data.flatten()))

    return GtmSwath(lat=lat, lon=lon, data=data, raw=np.concatenate(unscaled), paths=edr_filenames, info=info)


SCALE_TYPE = {
    'VNCCO': ('black2white', '0', '1.2'),  # albedo
    'VI1BO': ('black2white', '0', '1.2'),  # radiance
    'VI2BO': ('black2white', '0', '1.2'),  # radiance
    'VI3BO': ('black2white', '0', '1.2'),  # radiance
    'VI4BO': ('white2black','180','320'),  # BT
    'VI5BO': ('white2black','180','320'),  # BT
    # FIXME add VM##O
}


def gtm_quicklook(pathnames,
    png_fmt = 'viirs_%(kind)s %(date)s.%(start_time)s-%(end_time)s.png',
    label_fmt = 'Suomi s%(sat)  %(kind)s %(date)s.%(start_time)s-%(end_time)s',scale='COLOR',vmin=None,vmax=None,std=False,nosqrt=False,raw=False):

    LOG.info("Loading swath data")
    viirs = gtm_swath(*pathnames)
    eva = evaluator(**viirs.info)
    png_path = png_fmt % eva
 
    if std == True :
        std_max = (viirs.data).mean(dtype=np.float64) + 3*(viirs.data).std(dtype=np.float64)
        std_min = (viirs.data).mean(dtype=np.float64) -  3*(viirs.data).std(dtype=np.float64)
        vmin=std_min
        vmax=std_max
        nosqrt = True

        LOG.debug("mean: " + str( (viirs.data).mean(dtype=np.float64)) +" max: "+str((viirs.data).max())+" min: "+str((viirs.data).min())+" std "+str((viirs.data).std(dtype=np.float64))+" vmax: "+str(max)+" smax: "+str(std_max)+" smin: "+str(std_min))

    sqrt_enhance = False
    if scale == 'black2white' :
        sqrt_enhance = True
        
    if nosqrt:
        sqrt_enhance = False
 
    label = label_fmt % eva
    label = label.replace("Var", "")

    if (viirs.data).max() == -999. and viirs.data.min() == -999.:
        LOG.warn('No data for "%s"' % (png_path) )
        return 1

    LOG.info('rendering %s with label "%s"' % (png_path, label))

    if raw:
        raw_image(png_path, viirs.data, label = label,
                  scale=scale,sqrt_enhance=sqrt_enhance,
                  vmin=vmin, vmax=vmax)

    map_image(png_path,
              viirs.data, viirs.lon, viirs.lat, label=label,
              scale=scale,sqrt_enhance=sqrt_enhance,
              vmin=vmin, vmax=vmax, size=(600, 600))


def get_names():
    LOG.debug("Get names")
#    newpaths = list(glob(os.path.join(name, reg)))

# 
# pyresample to a map
# 

def viirs_quicklook(pathnames, 
                    png_fmt = 'viirs_%(kind)s%(band)s%(date)s.%(start_time)s-%(end_time)s.png', 
                    label_fmt = 'Suomi NPP %(kind)s%(band)s %(date)s.%(start_time)s-%(end_time)s',std=False,nosqrt=False,raw=False):

#    default_path=None
    if len(pathnames)==1 and os.path.isdir(*pathnames):
        raise NotImplementedError('directory iteration not yet supported')
        # FIXME : go through all the available datasets by iterating across SCALE_TYPE keys
        # default_path = pathnames[0]
        # for reg , scale , vmin,vmax in sorted ( SCALE_TYPE ):
        #    # get the files
        #     newpaths = sorted( list(glob(os.path.join(default_path, reg+"*.h5"))))
        #     if len ( newpaths ) > 0:
        #         gtm_quicklook(newpaths,
        #             png_fmt = png_fmt,
        #             label_fmt = label_fmt,scale=scale,vmin=vmin,vmax=vmax,std=std,nosqrt=nosqrt,raw=raw)
    else :
        file_to_match=pathnames[0]
        mdict = RE_NPP.match(os.path.split(pathnames[0])[-1]).groupdict()
        prefix = mdict['kind']

        scale_to_use, vmin_to_use, vmax_to_use = SCALE_TYPE[prefix]

        pathnames = sorted(pathnames)
        for file_to_check in pathnames:
            if not file_to_check.startswith(prefix):
                LOG.error("All files must be same band"+file_to_check+" "+prefix)
                sys.exit(1)
        gtm_quicklook(pathnames,
            png_fmt = png_fmt, 
            label_fmt = label_fmt,scale=scale_to_use,vmin=vmin_to_use,vmax=vmax_to_use,std=std,nosqrt=nosqrt,raw=raw)



def main():
    import optparse
    usage = """
%prog [options] ...

"""
    parser = optparse.OptionParser(usage)
    parser.add_option('-t', '--test', dest="self_test",
                    action="store_true", default=False, help="run self-tests") 
    parser.add_option('-v', '--verbose', dest='verbosity', action="count", default=0,
                    help='each occurrence increases verbosity 1 level through ERROR-WARNING-INFO-DEBUG')
    parser.add_option('-o', '--output', dest='output', default='.',
                     help='location to store output')
    
    parser.add_option('-s', '--std',
                    action="store_true", default=False, help="scale image to 3*std deviation")
    
    parser.add_option('-n', '--nosqrt',
                    action="store_true", default=False, help="disable sqrt enhancement")
    
    parser.add_option('-r', '--raw',
                    action="store_true", default=False, help="create unmapped quick look")



    # parser.add_option('-I', '--include-path', dest="includes",
    #                 action="append", help="include path to append to GCCXML call")                           
    (options, args) = parser.parse_args()

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
 
    std=False
    if options.std == True :
        std=True
        
    nosqrt=False
    if options.nosqrt == True :
        nosqrt=True
        
    raw=False
    if options.raw == True :
        raw=True

    if not args:
        parser.error( 'incorrect arguments, try -h or --help.' )
        return 9
    
    png_fmt = 'viirs_%(kind)s%(band)s_%(date)s.%(start_time)s-%(end_time)s.png'
    if options.output is not None:
        if os.path.isdir(options.output):
            png_fmt = os.path.join(options.output, png_fmt)
        else:
            png_fmt = options.output
    
    viirs_quicklook(args, png_fmt=png_fmt,std=std,nosqrt=nosqrt,raw=raw)

    return 0



if __name__=='__main__':
    sys.exit(main())
    
