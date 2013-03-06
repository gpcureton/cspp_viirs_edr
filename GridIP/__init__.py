#!/usr/bin/env python
# encoding: utf-8
"""
This module re-implements the GridIP gridded ingest and granulation
in the Algorithm Development Package (ADL).

Created by Geoff Cureton on 2013-03-05.
Copyright (c) 2013 University of Wisconsin SSEC. All rights reserved.
"""

file_Date = '$Date$'
file_Revision = '$Revision$'
file_Author = '$Author$'
file_HeadURL = '$HeadURL$'
file_Id = '$Id$'

__author__ = 'G.P. Cureton <geoff.cureton@ssec.wisc.edu>'
__version__ = '$Id$'
__docformat__ = 'Epytext'


#import ProCmnPhysConst
from LandWaterMask        import LandWaterMask
from QuarterlySurfaceType import QuarterlySurfaceType
from QstLwm               import QstLwm
from NbarNdvi17Day        import NbarNdvi17Day
from SnowIceCover         import SnowIceCover

classNames = {}
classNames['VIIRS-GridIP-VIIRS-Lwm-Mod-Gran'] = 'LandWaterMask'
classNames['VIIRS-GridIP-VIIRS-Qst-Mod-Gran'] = 'QuarterlySurfaceType'
classNames['VIIRS-GridIP-VIIRS-Qst-Lwm-Mod-Gran'] = 'QstLwm'
classNames['VIIRS-GridIP-VIIRS-Nbar-Ndvi-Mod-Gran'] = 'NbarNdvi17Day'
classNames['VIIRS-GridIP-VIIRS-Snow-Ice-Cover-Mod-Gran'] = 'SnowIceCover'
