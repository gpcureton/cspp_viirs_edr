#!/usr/bin/env python
# encoding: utf-8
"""
This module re-implements the ancillary gridded ingest and granulation
in the Algorithm Development Package (ADL).

Created by Geoff Cureton on 2013-02-25.
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

#import ANC
#import ANC_subclass
#from ANC_subclass import ANC_subclass

#import ProCmnPhysConst
from PrecipWater            import PrecipWater
from SurfTemp               import SurfTemp
from WindSpeed              import WindSpeed
from TerrainGeopotentialHeight import TerrainGeopotentialHeight

classNames = {}
classNames['VIIRS-ANC-Temp-Surf2M-Mod-Gran'] = 'SurfTemp'
classNames['VIIRS-ANC-Preci-Wtr-Mod-Gran'] = 'PrecipWater'
classNames['VIIRS-ANC-Wind-Speed-Mod-Gran'] = 'WindSpeed'
classNames['VIIRS-ANC-Surf-Ht-Mod-Gran'] = 'TerrainGeopotentialHeight'

#import ProCmnPhysConst
#import PrecipWater
#import SurfGeopotentialHeight
#import SurfTemp
#import WindSpeed

#import AotClimatology
#import Bathymetry
#import GeopotentialHeight
#import NhfOzone
#import NhfPresLevelTemp
#import NhfSpecSurfHumidity
#import NhfSurfPres
#import NhfSurfTemp
#import NhfWaterVaporMixRatio
#import NitrateDepletion
#import OpticalDepth
#import Ozone
#import PresLevelTemp
#import SkinTemp
#import SpecSurfHumidity
#import SurfPresCorrection
#import SurfPres
#import TerrainGeopotentialHeight
#import TropoGeopotentialHeight
#import WaterVaporMixRatio
#import WindDirection


