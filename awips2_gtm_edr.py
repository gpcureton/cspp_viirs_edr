#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
awips2_gtm_edr
~~~~~~~~~~~~~~

Transcode HDF5 format Ground-Track Mercator (GTM) Imagery EDR files to NetCDF4 compatible with AWIPS2.

Note on .nc.wmo files:
These are NetCDF4 files which have been gzipped, then prepended with a 21-character header e.g. "TIBP99 KNES 000000\r\r\n"
To unpack them back to NetCDF files at the command line,

dd ibs=1 iseek=21 if=VIIRS_I4_IMG_EDR_TIPB04_KNES_npp_s201307141706188_e201307141707531_c201307311502170.nc.gz.wmo |gunzip >awips2.nc
ncdump -h awips2.nc

To do a quick visualization using ipython, try this:

ipython --pylab
 from netCDF4 import Dataset
 nc = Dataset('awips2.nc', 'r+')
 bt = nc.variables['BrightnessTemperature@VIIRS-I4-IMG-EDR']   # variable name will vary by band; see ncdump output
 bt.missing_value = 65535   # this is a hack to work around multiple-missing-value issue, note it modifies the actual file
 imshow(bt[:,:])

From code written by DJHoese, Apr2013

:copyright: 2013 by University of Wisconsin Regents, see AUTHORS for more details
:license: GPLv3, see LICENSE for more details
:author: RKGarcia

"""
__author__ = 'rayg'
__docformat__ = 'reStructuredText'

import os
import sys
import re
import logging
import datetime

import h5py
from netCDF4 import Dataset
from collections import namedtuple
import numpy as np
from subprocess import PIPE, Popen


LOG = logging.getLogger(__name__)

# This is modified from adl_geo_ref:RE_NPP due to things like VI1BO instead of SVI01; kind and band are combined
RE_NPP_EDR = re.compile('(?P<kindband>[A-Z0-9]+)_(?P<sat>[A-Za-z0-9]+)_d(?P<date>\d+)_t(?P<start_time>\d+)_e(?P<end_time>\d+)_b(?P<orbit>\d+)_c(?P<created_time>\d+)_(?P<site>[a-zA-Z0-9]+)_(?P<domain>[a-zA-Z0-9]+)\.h5')


def h5path(elf, path, groups=None):
    LOG.debug('path search for %s' % repr(path))
    if not path:
        return elf
    if isinstance(path, str):
        return h5path(elf, path.split('/'), groups)
    LOG.debug('compiling %r' % path)
    assert(isinstance(path[0], str))
    rx = re.compile(path[0])
    for k, v in elf.iteritems():
        m = rx.match(k)
        LOG.debug('checking %s' % k)
        if not m:
            continue
        if groups is not None:
            LOG.debug('updating information to include %s' % repr(m.groupdict()))
            groups.update(m.groupdict())
        return h5path(v, path[1:], groups) if len(path) > 1 else v
    LOG.warning('no match for %s' % path[0])
    return None

EDR_PATH = r'All_Data/(?P<collection>VIIRS-.*-EDR)_All/(?P<kind>(?:BrightnessTemperature)|(?:Reflectance)|(?:Albedo))'
# where to find Latitude and Longitude
GEO_PATH = r'All_Data/(?P<geo_collection>VIIRS-.*-EDR-GEO)_All'
TIME_PATH = r'Data_Products/VIIRS-.*-EDR-GEO/VIIRS-.*-EDR-GEO_Aggr'
GRING_PATH = r'Data_Products/VIIRS-.*-EDR-GEO/VIIRS-.*-EDR-GEO_Gran_0'


latlon = namedtuple('g_ring', 'lat lon')


# Time helpers
class UTC(datetime.tzinfo):
    """Time zone class for UTC
    """
    ZERO = datetime.timedelta(0)

    def utcoffset(self, dt):
        return self.ZERO

    def tzname(self, dt):
        return "UTC"

    def dst(self, dt):
        return self.ZERO


def utc_now():
    return datetime.datetime.utcnow().replace(tzinfo=UTC())


def nppdatetime(d, t, e=None):
    """
    given d,t,e strings from filename, return datetime object or objects
    """
    assert(len(d) == 8)
    t = t.split('.')[0]
    assert(len(t) in (6,7))
    start = datetime.datetime.strptime(d + 'T' + t[:6], '%Y%m%dT%H%M%S')
    start = start.replace(tzinfo=UTC())
    if len(t) == 7:
        start = start.replace(microsecond=int(t[6]) * 100000)
    if not e:
        return start
    end = nppdatetime(d,e)
    if end < start:
        end += datetime.timedelta(days=1)
    return start, end


class Granule(object):
    """
    Tool for accessing necessary components of EDR + GEO HDF5 files in accordance to CDFCB
    """
    kind = None
    collection = None
    edr_path = None
    edr = None
    geo_path = None
    geo = None

    def __init__(self, edr_path, geo_path=None):
        self.edr_path = edr_path
        self.edr = h5py.File(edr_path, 'r')
        if geo_path is None:
            geo_ref = self.edr.attrs.get('N_GEO_Ref', None)
            if geo_ref:
                geo_ref = str(geo_ref[0,0])
                LOG.debug("N_GEO_Ref is {0:s}".format(geo_ref))
                _, geo_filename = os.path.split(geo_ref)
                geo_dir, _ = os.path.split(edr_path)
                geo_path = os.path.join(geo_dir, geo_filename)
                LOG.debug("using {0:s} as geo path via N_GEO_Ref".format(geo_path))
            else:
                LOG.info('no N_GEO_Ref; assuming integrated geolocation')
        if geo_path:
            self.geo_path = geo_path
            self.geo = h5py.File(geo_path, 'r')
        else:
            self.geo_path = edr_path
            self.geo = self.edr
        self._info = {}
        h5v = h5path(self.edr, EDR_PATH, self._info)
        self._data = h5v[:]
        self._geo_var = h5path(self.geo, GEO_PATH, self._info)

    @property
    def factors(self):
        """
        Find the RadianceFactors or BrightnessFactors hdf5 variable, depending on self.kind
        :return: numpy array
        """
        # this gets a little tricky
        # search for variable ending with Factors but starting with first 3 letters of kind
        pat = r'All_Data/(?P<collection>VIIRS-.*-EDR)_All/{0:s}.*Factors'.format(self.kind[:3])
        return h5path(self.edr, pat)[:]

    @property
    def data(self):
        return self._data

    @property
    def collection(self):
        return self._info['collection']

    @property
    def geo_collection(self):
        return self._info['geo_collection']

    @property
    def kind(self):
        return self._info['kind']

    @property
    def along_track_pixels(self):
        return self._data.shape[0]

    @property
    def cross_track_pixels(self):
        return self._data.shape[1]

    @property
    def factor_count(self):
        return self.factors.shape[0]

    @property
    def lat_lon(self):
        base = self._geo_var
        return latlon(base['Latitude'][:], base['Longitude'][:])

    @property
    def lat_lon_envelope(self):
        # Make data manipulations
        # FIXME: missing value testing? short granule testing?
        lat_data, lon_data = self.lat_lon
        mid_idx = lat_data.shape[-1] / 2
        lat_envelope = lat_data[:, (0, mid_idx, -1)]
        lon_envelope = lon_data[:, (0, mid_idx, -1)]
        return latlon(lat_envelope, lon_envelope)

    @property
    def start_end(self):
        base = h5path(self.geo, TIME_PATH)
        start_date = base.attrs["AggregateBeginningDate"][0, 0]
        start_time = base.attrs["AggregateBeginningTime"][0, 0]
        end_date = base.attrs["AggregateEndingDate"][0, 0]
        end_time = base.attrs["AggregateEndingTime"][0, 0]
        return nppdatetime(start_date, start_time), nppdatetime(end_date, end_time)

    @property
    def gring_lat_lon(self):
        # FIXME this doesn't deal with aggregates having more than one granule, consistency could be improved on h5 paths
        base = h5path(self.geo, GRING_PATH)
        return latlon(base.attrs["G-Ring_Latitude"][:].astype(np.float64), base.attrs["G-Ring_Longitude"][:].astype(np.float64))


class AWIPS2_NetCDF4(object):
    """
    Tool for creating AWIPS2 NetCDF files
    """
    _nc = None
    row_dim_name, col_dim_name = None, None
    fac_dim_name, env_dim_name = None, None

    def create_dimensions(self, along, cross, factors):
        # Create Dimensions
        _nc = self._nc
        self.row_dim_name = "AlongTrack-%d" % along
        self.fac_dim_name = "Granule-%d" % factors  # FIXME???
        self.col_dim_name = "CrossTrack-%d" % cross
        self.env_dim_name = "CrossTrack-3"   # FIXME??
        LOG.debug('along-track %d' % along)
        LOG.debug('cross-track %d' % cross)
        LOG.debug('factors %d' % factors)
        _nc.createDimension(self.row_dim_name, along)
        _nc.createDimension(self.fac_dim_name, factors)
        _nc.createDimension(self.col_dim_name, cross)
        _nc.createDimension(self.env_dim_name, 3)

    def create_time_attrs(self, sdt, edt):
        # Create Global Attributes
        self._nc.time_coverage_start = sdt.strftime("%Y-%m-%dT%H:%M:%SZ")
        self._nc.time_coverage_end = edt.strftime("%Y-%m-%dT%H:%M:%SZ")
        self._nc.date_created = utc_now().strftime("%Y-%m-%dT%H:%M:%SZ")

    def create_g_ring_attrs(self, g_ring_lat, g_ring_lon):
        g_ring_lat = np.append(g_ring_lat, g_ring_lat[0])
        g_ring_lon = np.append(g_ring_lon, g_ring_lon[0])
        self._nc.setncattr("g_ring_latitude", g_ring_lat)
        self._nc.setncattr("g_ring_longitude", g_ring_lon)

    def create_image_vars(self, var_stem, collection, data, factors):
        # Create and write Variables
        # Image data
        LOG.debug('data shape is {0:s}'.format(repr(data.shape)))
        bt_var = self._nc.createVariable("{0:s}@{1:s}".format(var_stem, collection), 'u2',
                                    dimensions=(self.row_dim_name, self.col_dim_name))
        bt_var[:, :] = data
        bt_var.setncattr("missing_value", "65535 65534 65533 65532 65531 65530 65529 65528")  # FUTURE: fix this, it can break NetCDF readers. Do we really need it?
        # Scaling Factors
        prefix = re.match(r'^([A-Z][a-z]+).*', var_stem).group(1)   # BrightnessTemperature -> Brightness
        LOG.debug('{0:s} is prefix'.format(prefix))

        bt_factors_var = self._nc.createVariable(
            "{0:s}Factors@{1:s}".format(prefix, collection),
            'f4',
            dimensions=(self.fac_dim_name,)
        )
        # Factors should be 2 but this aggregate has more
        bt_factors_var[:] = factors

    def create_geo_vars(self, collection, lat_envelope, lon_envelope):
        # remove -M##, -I# in geo collection
        # alpha = [x for x in collection.split('-') if not x[-1].isdigit()]
        # collection = '-'.join(alpha)

        # create navigation variables
        lat_var = self._nc.createVariable("Latitude@{0:s}".format(collection), 'f4',
                                          dimensions=(self.row_dim_name, self.env_dim_name))
        lat_var[:, :] = lat_envelope
        lat_var.setncattr("missing_value", "-999.9 -999.8 -999.5 -999.4 -999.3")

        lon_var = self._nc.createVariable("Longitude@{0:s}".format(collection), 'f4',
                                          dimensions=(self.row_dim_name, self.env_dim_name))
        lon_var[:, :] = lon_envelope
        lon_var.setncattr("missing_value", "-999.9 -999.8 -999.5 -999.4 -999.3")

    def __init__(self, filename):
        self._nc = Dataset(filename, 'w')

    def close(self):
        self._nc.sync()
        self._nc.close()
        self._nc = None


TITANIUM_LEAD = {
    'VI1BO': ('TIPB01', 'Reflectance', 'VIIRS_I1_IMG_EDR'),
    'VI2BO': ('TIPB02', 'Reflectance', 'VIIRS_I2_IMG_EDR'),
    'VI3BO': ('TIPB03', 'Reflectance', 'VIIRS_I3_IMG_EDR'),
    'VI4BO': ('TIPB04', 'Reflectance', 'VIIRS_I4_IMG_EDR'),
    'VI5BO': ('TIPB05', 'Reflectance', 'VIIRS_I5_IMG_EDR'),
    'VM01O': ('TIPB11', 'BrightnessTemperature', 'VIIRS_M1_EDR'),
    'VM02O': ('TIPB14', 'BrightnessTemperature', 'VIIRS_M4_EDR'),
    'VM03O': ('TIPB19', 'BrightnessTemperature', 'VIIRS_M9_EDR'),
    'VM04O': ('TIPB24', 'BrightnessTemperature', 'VIIRS_M14_EDR'),
    'VM05O': ('TIPB25', 'BrightnessTemperature', 'VIIRS_M15_EDR'),
    'VM06O': ('TIPB26', 'BrightnessTemperature', 'VIIRS_M16_EDR'),
    'VNCCO': ('TIPB27', 'Albedo', 'VIIRS_NCC_EDR'),
}


def _ncdatefmt(dt):
    return dt.strftime('%Y%m%d%H%M%S') + ('%1d' % (dt.microsecond / 100000))


nc_info = namedtuple('nc_info', 'filename_stem tipb_id station ddhhmm variable_stem')


def _nc_info_from_edr_path(edr_path, station=None):
    "convert granule and source information into a netcdf filename and supporting header information"
    dn, fn = os.path.split(edr_path)
    m = RE_NPP_EDR.match(fn)
    if not m:
        raise ValueError('{0:s} is not a valid CDFCB-compliant NPP pathname'.format(gran.edr_path))
    g = m.groupdict()
    sat, d, t, e, c, site, kind_band = map(lambda x: g[x], ('sat', 'date', 'start_time', 'end_time', 'created_time', 'site', 'kindband'))
    tipb_id, nc_var_stem, stem = TITANIUM_LEAD.get(kind_band, (None, None, None))
    if not tipb_id:
        raise ValueError('{0:s} is not a known EDR type'.format(kind_band))
    sdt, edt = nppdatetime(d,t,e)
    ncs = _ncdatefmt(sdt)
    nce = _ncdatefmt(edt)
    creation = c[:15]  # truncate
    ddhhmm = sdt.strftime('%d%H%M')
    return nc_info(
            '{stem:s}_{tipb_id:s}_{station:s}_{sat:s}_s{ncs:s}_e{nce:s}_c{creation:s}'.format(**locals()),
            tipb_id,
            station or site.upper(),
            ddhhmm,
            nc_var_stem)


WMO_SUFFIX = '.gz.wmo'


def wmo_wrap(nc_path, wmo_header="TIPB99 KNES 000000", wmo_path=None):
    """
    create AWIPS2 wrapper, which involves a text header followed by a gzipped netcdf file.
    Yeah, about that...
    """
    wmo_header = wmo_header + '\r\r\n'
    wmo_path = nc_path + WMO_SUFFIX if wmo_path is None else wmo_path
    LOG.debug('writing {0:s} with header {1!r:s}'.format(wmo_path, wmo_header))
    gzfp = open(wmo_path, 'wb')
    gzfp.write(wmo_header)
    gzfp.flush()
    ncfp = open(nc_path, 'rb')
    gz = Popen(['gzip'], stdin=ncfp, stdout=gzfp)
    _,_ = gz.communicate()


def transform(edr_path, output_dir=None, geo_path=None, station=None, wmo_path=None, cleanup=True):
    LOG.info('opening files')
    gran = Granule(edr_path, geo_path)
    ncstem, tipb_id, station, ddhhmm, varname = _nc_info_from_edr_path(edr_path, station=station)
    ncfn = ncstem + '.nc'
    if wmo_path is None:
        wmo_path = ncfn + WMO_SUFFIX
    if output_dir is not None:
        wmo_path = os.path.join(output_dir, os.path.split(wmo_path)[-1])
    start, end = gran.start_end
    LOG.debug('start, end = {0:s}, {1:s}'.format(start, end))

    LOG.info('creating output file {0:s}'.format(ncfn))
    nc = AWIPS2_NetCDF4(ncfn)
    LOG.debug('adding dimensions')
    nc.create_dimensions(gran.along_track_pixels, gran.cross_track_pixels, gran.factor_count)
    LOG.debug('adding time attributes')
    nc.create_time_attrs(start, end)
    LOG.debug('accessing G-Ring navigation')
    gr = gran.gring_lat_lon
    LOG.debug("writing G-Ring attributes")
    nc.create_g_ring_attrs(g_ring_lat=gr.lat, g_ring_lon=gr.lon)
    LOG.debug('transferring image data')
    nc.create_image_vars(varname, gran.collection, gran.data, gran.factors)
    LOG.debug('transferring lat-lon envelope')
    env = gran.lat_lon_envelope
    nc.create_geo_vars(gran.geo_collection, lat_envelope=env.lat, lon_envelope=env.lon)
    LOG.debug('writing NetCDF4 file')
    nc.close()

    LOG.debug('wrapping NetCDF4 as WMO')
    wmo_header = '{0:s} {1:s} {2:s}'.format(tipb_id, station, ddhhmm)
    wmo_wrap(ncfn, wmo_header=wmo_header, wmo_path=wmo_path)
    if cleanup:
        LOG.debug('cleaning out intermediate file {0:s}'.format(ncfn))
        os.unlink(ncfn)
    LOG.info('done!')
    return wmo_path


def main():
    import optparse
    usage = """
%prog [options] VIIRS-imagery-file.h5 {VIIRS-geo-file.h5} {output-nc-file}
"""
    parser = optparse.OptionParser(usage)
    #parser.add_option('-t', '--test', dest="self_test",
    #                action="store_true", default=False, help="run self-tests")
    parser.add_option('-v', '--verbose', dest='verbosity', action="count", default=0,
                    help='each occurrence increases verbosity 1 level through ERROR-WARNING-INFO-DEBUG')
    parser.add_option('-d', '--debug', dest='debug', action="count", default=0,
                      help='enable debug mode where clean-up does not occur (results in .nc file creation)')
    parser.add_option('-s', '--station', dest='station', default='KNES', type='str',
                      help='Station of origin, which is placed in headers')
    parser.add_option('-o', '--output', dest='output', default=None,
                     help='destination directory to store output to')
    # parser.add_option('-I', '--include-path', dest="includes",
    #                 action="append", help="include path to append to GCCXML call")
    (options, args) = parser.parse_args()

    # make options a globally accessible structure, e.g. OPTS.
    global OPTS
    OPTS = options

    #if options.self_test:a
    #    # import doctest
    #    # doctest.testmod()
    #    sys.exit(2)

    levels = [logging.ERROR, logging.WARN, logging.INFO, logging.DEBUG]
    logging.basicConfig(level=levels[min(3,options.verbosity)])

    if not args:
        parser.error( 'incorrect arguments, try -h or --help.' )
        return 9

    for arg in args:
        if os.path.isfile(arg):
            transform(arg, station=options.station, output_dir=options.output, cleanup=not options.debug)
        elif os.path.isdir(arg):
            from glob import glob
            pat = os.path.join(arg, 'V???O*h5')   # VI?BO, VM??O, VNCCO
            for edr_path in glob(pat):
                LOG.info('processing {0:s}'.format(edr_path))
                transform(edr_path, station=options.station, output_dir=options.output, cleanup=not options.debug)
        else:
            LOG.warning('really not sure what to do with {0!r:s} - ignoring'.format(arg))

    return 0


if __name__=='__main__':
    sys.exit(main())



DJHoese = """
# Read the input files
img_file = h5py.File("./VI4BO_npp_d20130116_t0944041_e0945402_b06328_c20130305060723097136_noaa_ops.h5")
geo_file = h5py.File("./GIGTO_npp_d20130116_t0944041_e0945402_b06328_c20130305060441825409_noaa_ops.h5")

img_data = img_file["All_Data"]["VIIRS-I4-IMG-EDR_All"]["BrightnessTemperature"][:,:]
img_factors = img_file["All_Data"]["VIIRS-I4-IMG-EDR_All"]["BrightnessFactors"][:]
lat_data = geo_file["All_Data"]["VIIRS-IMG-GTM-EDR-GEO_All"]["Latitude"][:,:]
lon_data = geo_file["All_Data"]["VIIRS-IMG-GTM-EDR-GEO_All"]["Longitude"][:,:]
start_date = geo_file["Data_Products"]["VIIRS-IMG-GTM-EDR-GEO"]["VIIRS-IMG-GTM-EDR-GEO_Aggr"].attrs["AggregateBeginningDate"][0,0]
start_time = geo_file["Data_Products"]["VIIRS-IMG-GTM-EDR-GEO"]["VIIRS-IMG-GTM-EDR-GEO_Aggr"].attrs["AggregateBeginningTime"][0,0]
end_date = geo_file["Data_Products"]["VIIRS-IMG-GTM-EDR-GEO"]["VIIRS-IMG-GTM-EDR-GEO_Aggr"].attrs["AggregateEndingDate"][0,0]
end_time = geo_file["Data_Products"]["VIIRS-IMG-GTM-EDR-GEO"]["VIIRS-IMG-GTM-EDR-GEO_Aggr"].attrs["AggregateEndingTime"][0,0]
g_ring_lat = geo_file["Data_Products"]["VIIRS-IMG-GTM-EDR-GEO"]["VIIRS-IMG-GTM-EDR-GEO_Gran_0"].attrs["G-Ring_Latitude"][:].astype(numpy.float64)
g_ring_lon = geo_file["Data_Products"]["VIIRS-IMG-GTM-EDR-GEO"]["VIIRS-IMG-GTM-EDR-GEO_Gran_0"].attrs["G-Ring_Longitude"][:].astype(numpy.float64)

# Make data manipulations
mid_idx = lat_data.shape[-1]/2
lat_envelope = lat_data[:,(0,mid_idx,-1)]
lon_envelope = lon_data[:,(0,mid_idx,-1)]

# Write the output file
_nc = Dataset("VIIRS_I4_IMG_EDR_TIPB99_KNES_npp_s201301160944041_e201301160945402_c20130305060441825409.nc", mode="w")
nc_file = Dataset("viirs_img_edr_20121001.nc", mode="w")

# Create Dimensions
row_dim_name = "AlongTrack-%d" % img_data.shape[0]
fac_dim_name = "Granule-%d" % img_factors.shape[0]
col_dim_name = "CrossTrack-%d" % img_data.shape[1]
env_dim_name = "CrossTrack-%d" % 3
_nc.createDimension(row_dim_name, img_data.shape[0])
_nc.createDimension(fac_dim_name, img_factors.shape[0])
_nc.createDimension(col_dim_name, img_data.shape[1])
_nc.createDimension(env_dim_name, 3)

# Create Global Attributes
_nc.time_coverage_start = datetime.datetime.strptime(start_date + start_time.split(".")[0], "%Y%m%d%H%M%S").strftime("%Y-%m-%dT%H:%M:%SZ")
_nc.time_coverage_end   = datetime.datetime.strptime(end_date + end_time.split(".")[0], "%Y%m%d%H%M%S").strftime("%Y-%m-%dT%H:%M:%SZ")
_nc.date_created = utc_now().strftime("%Y-%m-%dT%H:%M:%SZ")
print g_ring_lat.dtype
print g_ring_lat
print g_ring_lon
g_ring_lat = np.append(g_ring_lat, g_ring_lat[0])
g_ring_lon = np.append(g_ring_lon, g_ring_lon[0])
_nc.setncattr("g_ring_latitude", g_ring_lat)
_nc.setncattr("g_ring_longitude", g_ring_lon)

# Create and write Variables
# Image data
bt_var = _nc.createVariable("BrightnessTemperature@VIIRS-I4-IMG-EDR", 'u2', dimensions=(row_dim_name,col_dim_name))
bt_var[:,:] = img_data
bt_var.setncattr("missing_value", "65535 65534 65533 65532 65531 65530 65529 65528")
# XXX: Need missing_valuename?

# Scaling Factors
bt_factors_var = _nc.createVariable("BrightnessFactors@VIIRS-I4-IMG-EDR", 'f4', dimensions=(fac_dim_name,))
# Factors should be 2 but this aggregate has more
bt_factors_var[:] = img_factors

# Navigation
lat_var = _nc.createVariable("Latitude@VIIRS-IMG-GTM-EDR-GEO", 'f4', dimensions=(row_dim_name,env_dim_name))
lat_var[:,:] = lat_envelope
lat_var.setncattr("missing_value", "-999.9 -999.8 -999.5 -999.4 -999.3")
lon_var = _nc.createVariable("Longitude@VIIRS-IMG-GTM-EDR-GEO", 'f4', dimensions=(row_dim_name,env_dim_name))
lon_var[:,:] = lon_envelope
lon_var.setncattr("missing_value", "-999.9 -999.8 -999.5 -999.4 -999.3")

_nc.sync()
_nc.close()

"""
