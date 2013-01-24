#!/usr/bin/env python
# encoding: utf-8
"""
HDF4File.py

Created by Geoff Cureton on 2012-09-25.
Copyright (c) 2012 University of Wisconsin SSEC. All rights reserved.
"""

file_Date = '$Date$'
file_Revision = '$Revision$'
file_Author = '$Author$'
file_HeadURL = '$HeadURL$'
file_Id = '$Id$'

__author__ = 'Geoff Cureton <geoff.cureton@ssec.wisc.edu>'
__version__ = '$Id$'
__docformat__ = 'Epytext'


import string

import pyhdf.HDF as HDF
import pyhdf.V   as V
import pyhdf.VS  as VS
import pyhdf.SD  as SD

# every module should have a LOG object
import logging
#sourcename= file_Id.split(" ")
sourcename=__name__
LOG = logging.getLogger(sourcename[1])

class HDF4File():
    """
    Common interface for (reading) HDF4 files.
    """

    HDFobjTag = {
            'dataset':  HDF.HC.DFTAG_NDG,
            'vdata':    HDF.HC.DFTAG_VH,
            'vgroup':   HDF.HC.DFTAG_VG
            }

    def __init__(self,filename):
        self.filename = filename

    def readData(self, varName):
        inFile = SD.SD(self.filename, SD.SDC.READ)
        return inFile.select(varName).get()
  

    """
    Quick fix for large datasets that we don't want to use 'get' with right away.
    """
    def readLargeData(self, varName):
        inFile = SD.SD(self.filename, SD.SDC.READ)
        return inFile.select(varName)


    def readAttribute(self, varName, attrName):
        inFile = SD.SD(self.filename, SD.SDC.READ)
    
        inAttrs = inFile.select(varName).attributes()
    
        if attrName in inAttrs:
            return inAttrs[attrName]
        else:
            return None
    
  
    # Small wrapper for getSDS
    def getDataset(self, varName):
        return self._getSDS(varName).get()
    
    
    # Module for specifying the whole vgroup path to an SDS
    # in a UNIX-style path name.  Returns the SDS.
    #
    # Once trusted, this should probably be rolled in with readData
    # so that only one needs to be used!
    def _getSDS(self, varPath):
    
        # Open HDF file in readonly mode.    
        hdf = HDF.HDF(self.filename)
      
        # Initialize the SD, V and VS interfaces on the file.
        sd = SD.SD(self.filename)
        vs = hdf.vstart()
        v  = hdf.vgstart()
      
        vars = string.split(varPath, "/")
        LOG.info("vars : %r" %(vars))
      
        if len(vars) > 0:
            LOG.info("Specifying '%s' as root group." % vars[0])
            vgroup = v.attach(vars[0])
            
            for i in range(0, len(vars)):
                varname = vars[i]
                LOG.info("Traversing level %s with varname '%s'"%(str(i),varname))
                
                members = vgroup.tagrefs()
                LOG.info("%s has members: %r" % (varname,members))

                for tag, ref in members:
                    LOG.info("\ttag,ref = (%r,%r)" %(tag,ref))
                
                    # Vdata tag
                    if tag == HDF.HC.DFTAG_VH:
                        LOG.info("\tWe have a vdata tag (DFTAG_VH)")
                        vd = vs.attach(ref)
                        nrecs, intmode, fields, size, name = vd.inquire()
                        LOG.info("\tvdata: %s %s,%s" % (name,tag,ref))
                        LOG.info("\tfields: %s" %(fields))
                        LOG.info("\tnrecs: %d" %(nrecs))
                        vd.detach()
                   
                    # SDS tag
                    elif tag == HDF.HC.DFTAG_NDG:
                        LOG.info("\tWe have a dataset tag (DFTAG_NDG)")
                        sds = sd.select(sd.reftoindex(ref))
                        name, rank, dims, type, nattrs = sds.info()
                        LOG.info("\tWe have dataset %s" %(name))
                        if name == varname:
                            # Check that if we reach our desired name, that it is the 
                            # last one (and hence likely to be a dataset rather than a group).
                            # TODO: Make this more understandable and robust.
                            if  i != (len(vars)-1):
                                LOG.warn("Only at '%s' and reached SDS!  Returning  '%s'\n" %(varname,name))
                            return sds
                        else:
                            sds.endaccess()
                   
                    # VG tag
                    elif tag == HDF.HC.DFTAG_VG:
                        LOG.info("\tWe have a vgroup tag (DFTAG_VG)")
                        vg0 = v.attach(ref)
                        if vg0._name == varname:
                            vgroup.detach()
                            vgroup = vg0
                   
                    # Unhandled tag
                    else:
                        LOG.warn("unhandled (tag,ref): %r,%r" % (tag,ref))
      
        # Terminate V, VS and SD interfaces.
        v.end()
        vs.end()
        sd.end()
      
        # Close HDF file.
        hdf.close()

