#!/usr/bin/env python
# encoding: utf-8
"""
viirs_cm_ip_ac_int.py

Purpose: Edit a VIIRS-CM-IP-AC-Int blob using data from its XML PCT

Input:
    * Guide file defining the internal structure of the blob file.
    * Blob file endianness.
    * New blob file copied from existing file with the same structure.
    * Data Quality Monitoring (DQM) XML file giving the blob data values.

Copyright 2013, University of Wisconsin Regents.

Created by Geoff Cureton <geoff.cureton@ssec.wisc.edu> on 2013-04-29.
Copyright (c) 2013-2013 University of Wisconsin Regents. All rights reserved.

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

import os,sys,logging
import xml.etree.ElementTree as ET
from os import path

import adl_blob2 as adl

# every module should have a LOG object
LOG = logging.getLogger(__file__)


def expr_iter(xmldataname):
    '''
    Trawls through the xml file and finds all instances of the "name":"value" 
    attribute pairs in the "coeff" blocks.
    '''
    xml = ET.fromstring(file(xmldataname, 'rt').read())
    (coeffs_node,) = xml.findall('coefficients')
    coeffs = coeffs_node.findall('coeff')
    for coeff in coeffs:
        lhs = coeff.find('name').text
        rhs = coeff.find('value').text
        yield lhs, rhs

# these are used by the expressions in the XML
CONST = dict( true = 1, false = 0 )


class blob_as_mapping(object):
    it = None
    def __init__(self,thing):
        self.it = thing
    def __getitem__(self, name):
        if name in CONST: return CONST[name]
        return getattr(self.it, name)
    def __setitem__(self, name, value):
        return setattr(self.it, name, value)


def process(xmlguidename, xmldataname, blobname, blobEndianness):
    if blobEndianness=='little':
        endian=adl.LITTLE_ENDIAN
    elif blobEndianness=='big':
        endian=adl.BIG_ENDIAN
    else :
        LOG.error('Incorrect endianness {} specified. Should either be "big" or"little"'
                .format(blobEndianness))
        sys.exit(1)

    xmlguidename = path.expanduser(xmlguidename)
    xmldataname = path.expanduser(xmldataname)
    blobname = path.expanduser(blobname)

    blob_obj = adl.map(xmlguidename, blobname, writable=True, endian=endian)

    for name, value in expr_iter(xmldataname):
        LOG.debug('Updating the field {} in file {}'.format(name,blobname))
        LOG.debug('Data XML {} = {}'.format(name,value))

        try:
            LOG.debug("altering field {}".format(name))
            old_value = eval(name, None, blob_obj[0])
            LOG.info('Old {} = {}'.format(name,old_value))
        except IndexError as huhwhut:
            LOG.error('unable to access field "{}", present in data but missing from schema'
                    .format(name))
            continue

        # Make a string containing a python expression
        expr = "{} = {}".format(name,value)

        exec expr in globals(), blob_obj[0]
        
        new_value = eval(name, None, blob_obj[0])
        LOG.info('New {} = {}'.format(name,new_value))

        if old_value != new_value :
            LOG.warn('"{} = {}" changed to {}'.format(name,old_value,new_value))

    blob_obj.sync()


def _argparse():
    '''
    Method to encapsulate the option parsing and various setup tasks.
    '''

    import argparse

    endianChoices = ['little','big']

    defaults = {
                'blobEndianness':'big',
                }

    description = '''Edit a VIIRS-CM-IP-AC-Int blob using data from its XML PCT.'''
    usage = "usage: %prog [mandatory args] [options]"
    version = __version__

    parser = argparse.ArgumentParser(
                                     #usage=usage,
                                     #version=version,
                                     description=description
                                     )

    # Mandatory arguments

    parser.add_argument('-g','--guide_xml',
                      action="store",
                      dest="xmlguidename",
                      type=str,
                      required=True,
                      help='''Guide file defining the internal structure of the
                      blob file.'''
                      )

    parser.add_argument('-d','--data_xml',
                      action="store",
                      dest="xmldataname",
                      type=str,
                      required=True,
                      help='''Data Quality Monitoring (DQM) XML file giving the 
                      blob data values.'''
                      )

    parser.add_argument('-b','--data_blob',
                      action="store",
                      dest="blobname",
                      type=str,
                      required=True,
                      help='''New blob file copied from existing file with the 
                      same structure.'''
                      )

    # Optional arguments 

    parser.add_argument('-e','--endianness',
                      action="store",
                      dest="blobEndianness",
                      type=str,
                      default=defaults['blobEndianness'],
                      choices=endianChoices,
                      help='''The endiannessof the coefficient blob file.\n\n
                              Possible values are...
                              {}. [default: '{}']
                           '''.format(endianChoices.__str__()[1:-1],
                               defaults['blobEndianness'])
                      )

    parser.add_argument('-v', '--verbose',
                      dest='verbosity',
                      action="count",
                      default=0,
                      help='''Each occurrence increases 
                      verbosity 1 level from INFO: -v=DEBUG'''
                      )

    args = parser.parse_args()

    # Set up the logging
    console_logFormat = '{} : {}:{}:{}:{}:  {}'.format(
            '%(asctime)s',
            '(%(levelname)s)',
            '%(filename)s',
            '%(funcName)s',
            '%(lineno)d',
            '  %(message)s',
            )
    dateFormat = '%Y-%m-%d %H:%M:%S'
    levels = [logging.ERROR, logging.WARN, logging.INFO, logging.DEBUG]
    logging.basicConfig(level = levels[args.verbosity], 
            format = console_logFormat, 
            datefmt = dateFormat)


    return args


def main():

    options = _argparse()

    process(options.xmlguidename, 
            options.xmldataname, 
            options.blobname, 
            options.blobEndianness)


if __name__=='__main__':
    sys.exit(main())  
    
