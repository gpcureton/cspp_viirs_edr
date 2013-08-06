#!/usr/bin/env python
# encoding: utf-8
"""
This module contains data required for any VIIRS EDR algorithms.

Created by Geoff Cureton on 2013-02-27.
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

import CloudMaskIP
import AerosolOpticalThicknessIP
import SeaSurfaceTemperatureEDR
import SurfaceReflectanceIP

modules = {}
modules['VCM'] = 'CloudMaskIP'
modules['AOT'] = 'AerosolOpticalThicknessIP'
modules['SST'] = 'SeaSurfaceTemperatureEDR'
modules['SRFREF'] = 'SurfaceReflectanceIP'

crossGranules = {}
crossGranules['VCM'] = 1
crossGranules['AOT'] = 1
crossGranules['SST'] = 0
crossGranules['SRFREF'] = 0

prerequisites = {}
prerequisites['VCM'] = []
prerequisites['AOT'] = ['VCM']
prerequisites['SST'] = ['VCM','AOT']
prerequisites['SRFREF'] = ['VCM','AOT']
