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


#import ProCmnPhysConst
from PrecipWater               import PrecipWater
from SurfTemp                  import SurfTemp
from WindSpeed                 import WindSpeed
from TerrainGeopotentialHeight import TerrainGeopotentialHeight

classNames = {}
classNames['VIIRS-ANC-Temp-Surf2M-Mod-Gran'] = 'SurfTemp'
classNames['VIIRS-ANC-Preci-Wtr-Mod-Gran'] = 'PrecipWater'
classNames['VIIRS-ANC-Wind-Speed-Mod-Gran'] = 'WindSpeed'
classNames['VIIRS-ANC-Surf-Ht-Mod-Gran'] = 'TerrainGeopotentialHeight'

#from AotClimatology            import AotClimatology
#from Bathymetry                import Bathymetry
#from GeopotentialHeight        import GeopotentialHeight
#from NhfOzone                  import NhfOzone
#from NhfPresLevelTemp          import NhfPresLevelTemp
#from NhfSpecSurfHumidity       import NhfSpecSurfHumidity
#from NhfSurfPres               import NhfSurfPres
#from NhfSurfTemp               import NhfSurfTemp
#from NhfWaterVaporMixRatio     import NhfWaterVaporMixRatio
#from NitrateDepletion          import NitrateDepletion
#from OpticalDepth              import OpticalDepth
#from Ozone                     import Ozone
#from PresLevelTemp             import PresLevelTemp
#from SkinTemp                  import SkinTemp
#from SpecSurfHumidity          import SpecSurfHumidity
#from SurfPresCorrection        import SurfPresCorrection
#from SurfPres                  import SurfPres
#from TerrainGeopotentialHeight import TerrainGeopotentialHeight
#from TropoGeopotentialHeight   import TropoGeopotentialHeight
#from WaterVaporMixRatio        import WaterVaporMixRatio
#from WindDirection             import WindDirection


