#!/usr/bin/env python
# encoding: utf-8
"""
ANC_subclass.py

 * DESCRIPTION:  Subclass of the ANC class, whch be used to implement the individual ancillary
                 classes. 

Created by Geoff Cureton on 2013-02-28.
Copyright (c) 2011 University of Wisconsin SSEC. All rights reserved.
"""

file_Date = '$Date$'
file_Revision = '$Revision$'
file_Author = '$Author$'
file_HeadURL = '$HeadURL$'
file_Id = '$Id$'

__author__ = 'G.P. Cureton <geoff.cureton@ssec.wisc.edu>'
__version__ = '$Id$'
__docformat__ = 'Epytext'

import logging

# every module should have a LOG object
LOG = logging.getLogger('ANC_subclass')

# import the superclass
from ANC import ANC

class ANC_subclass(ANC) :

    def __init__(self):

        self.collectionShortName = 'VIIRS-ANC-My-Anc-Type-Mod-Gran'
        self.xmlName = 'VIIRS_ANC_MY_ANC_TYPE_MOD_GRAN.xml'
        self.blobDatasetName = 'myAncDataType'
        self.dataType = 'float64'

    def granulate(self):
        '''
        Granulate the ancillary dataset.
        '''
        pass
