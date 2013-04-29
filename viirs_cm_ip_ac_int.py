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
Licensed under GNU Public License (GPL) v3. See http://www.gnu.org/licenses/gpl-3.0-standalone.html

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
import optparse as optparse
from os import path
import uuid
from datetime import datetime,timedelta


import adl_blob as adl

# every module should have a LOG object
sourcename= file_Id.split(" ")
LOG = logging.getLogger(sourcename[1])
from adl_common import configure_logging


def getURID() :
    '''
    Create a new URID to be used in making the asc filenames
    '''
    
    URID_dict = {}

    URID_timeObj = datetime.utcnow()
    
    creationDateStr = URID_timeObj.strftime("%Y-%m-%d %H:%M:%S.%f")
    creationDate_nousecStr = URID_timeObj.strftime("%Y-%m-%d %H:%M:%S.000000")
    
    tv_sec = int(URID_timeObj.strftime("%s"))
    tv_usec = int(URID_timeObj.strftime("%f"))
    hostId_ = uuid.getnode()
    thisAddress = id(URID_timeObj)
    
    l = tv_sec + tv_usec + hostId_ + thisAddress
    
    URID = '-'.join( ('{0:08x}'.format(tv_sec)[:8],
                      '{0:05x}'.format(tv_usec)[:5],
                      '{0:08x}'.format(hostId_)[:8],
                      '{0:08x}'.format(l)[:8]) )
    
    URID_dict['creationDateStr'] = creationDateStr
    URID_dict['creationDate_nousecStr'] = creationDate_nousecStr
    URID_dict['tv_sec'] = tv_sec
    URID_dict['tv_usec'] = tv_usec
    URID_dict['hostId_'] = hostId_
    URID_dict['thisAddress'] = thisAddress
    URID_dict['URID'] = URID
    
    return URID_dict


def expr_iter(xmldataname):
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
        LOG.error('Incorrect endianness %r specified. Should either be "big" or"little"'%(blobEndianness))
        sys.exit(1)

    xmlguidename = path.expanduser(xmlguidename)
    xmldataname = path.expanduser(xmldataname)
    blobname = path.expanduser(blobname)

    cc = adl.map(xmlguidename, blobname, writable=True, endian=endian)
    ccm = blob_as_mapping(cc)
    for lhs, rhs in expr_iter(xmldataname):
        expr = "%s = %s" % (lhs,rhs)
        old_value = eval(lhs, None, ccm)
        exec expr in globals(), ccm
        changed = eval(lhs, None, ccm) != old_value
        LOG.info(repr(expr) + ' previously %r %s' % (old_value, '*'*32 if changed else ''))
    cc.sync()


def main():

    endianChoices = ['little','big']

    description = '''Edit a VIIRS-CM-IP-AC-Int blob using data from its XML PCT.'''
    usage = "usage: %prog [mandatory args] [options]"
    version = "%prog "+__version__

    parser = optparse.OptionParser(description=description,usage=usage,version=version)

    # Mandatory arguments

    mandatoryGroup = optparse.OptionGroup(parser, "Mandatory Arguments",
                        "At a minimum these arguments must be specified")

    mandatoryGroup.add_option('-g','--guide_xml',
                      action="store",
                      dest="xmlguidename",
                      type="string",
                      help="Guide file defining the internal structure of the blob file.")

    mandatoryGroup.add_option('-d','--data_xml',
                      action="store",
                      dest="xmldataname",
                      type="string",
                      help="Data Quality Monitoring (DQM) XML file giving the blob data values.")

    mandatoryGroup.add_option('-b','--data_blob',
                      action="store",
                      dest="blobname",
                      type="string",
                      help="New blob file copied from existing file with the same structure.")

    parser.add_option_group(mandatoryGroup)


    # Optional arguments 

    optionalGroup = optparse.OptionGroup(parser, "Extra Options",
                         "These options may be used to customize behaviour of this program.")

    optionalGroup.add_option('-e','--sdr_endianness',
                      action="store",
                      dest="blobEndianness",
                      type="choice",
                      default='big',
                      choices=endianChoices,
                      help='''The input VIIRS SDR endianness.\n\n
                              Possible values are...
                              %s. [default: 'little']
                           ''' % (endianChoices.__str__()[1:-1]))

    optionalGroup.add_option('-v', '--verbose',
                      dest='verbosity',
                      action="count",
                      default=0,
                      help='each occurrence increases verbosity 1 level from ERROR: -v=WARNING -vv=INFO -vvv=DEBUG')

    parser.add_option_group(optionalGroup)


    # Parse the arguments from the command line
    (options, args) = parser.parse_args()

    # Check that all of the mandatory options are given. If one or more 
    # are missing, print error message and exit...
    mandatories = ['xmlguidename','xmldataname','blobname']
    mand_errors = ["Missing mandatory argument [-g guide_xml --guide_xml=guide_xml]",
                   "Missing mandatory argument [-d data_xml --data_xml=data_xml]",
                   "Missing mandatory argument [-b data_blob --data_blob=data_blob]"
                   ]
    isMissingMand = False
    for m,m_err in zip(mandatories,mand_errors):
        if not options.__dict__[m]:
            isMissingMand = True
            parser.error(m_err)
    if isMissingMand :
        parser.error("Incomplete mandatory arguments, aborting...")

    # Set up the logging
    levels = [logging.ERROR, logging.WARN, logging.INFO, logging.DEBUG]
    configure_logging(level = levels[min(options.verbosity,3)])

    process(options.xmlguidename, \
            options.xmldataname, \
            options.blobname, \
            options.blobEndianness)


if __name__=='__main__':
    sys.exit(main())  
    #print sys.argv
    #if len(sys.argv)<5:
        #print "Usage: viirs-cm-ip-ac-int.py xml-blob-guide-file VIIRS-CM-IP-AC-Int_blob blobEndianness  xml-pct-file"
    #else:
        #process(*sys.argv[1:])
    
