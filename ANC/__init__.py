#!/usr/bin/env python
# encoding: utf-8
"""
This module re-implements the ANC gridded ingest and granulation
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
from Utils                     import retrieve_NCEP_grib_files
from Utils                     import create_NCEP_grid_blobs
from Utils                     import create_NAAPS_grid_blobs
from PrecipWater               import PrecipWater
from SurfTemp                  import SurfTemp
from WindSpeed                 import WindSpeed
from WindDirection             import WindDirection
from TerrainGeopotentialHeight import TerrainGeopotentialHeight
from SurfPres                  import SurfPres
from Ozone                     import Ozone
from OpticalDepth              import OpticalDepth
from SurfGeopotentialHeight    import SurfGeopotentialHeight
from SpecSurfHumidity          import SpecSurfHumidity
from SkinTemp                  import SkinTemp

classNames = {}
classNames['VIIRS-ANC-Preci-Wtr-Mod-Gran'] = 'PrecipWater'
classNames['VIIRS-ANC-Temp-Surf2M-Mod-Gran'] = 'SurfTemp'
classNames['VIIRS-ANC-Wind-Speed-Mod-Gran'] = 'WindSpeed'
classNames['VIIRS-ANC-Wind-Direction-Mod-Gran'] = 'WindDirection'
classNames['VIIRS-ANC-Surf-Ht-Mod-Gran'] = 'TerrainGeopotentialHeight'
classNames['VIIRS-ANC-Press-Surf-Mod-Gran'] = 'SurfPres'
classNames['VIIRS-ANC-Tot-Col-Mod-Gran'] = 'Ozone'
classNames['VIIRS-ANC-Optical-Depth-Mod-Gran'] = 'OpticalDepth'
classNames['VIIRS-ANC-Geopot-Ht-Lev-Mod-Gran'] = 'SurfGeopotentialHeight'
classNames['VIIRS-ANC-Sp-Humd-Surf-Mod-Gran'] = 'SpecSurfHumidity'
classNames['VIIRS-ANC-Temp-Skin-Mod-Gran'] = 'SkinTemp'

#from AotClimatology            import AotClimatology
#from Bathymetry                import Bathymetry
#from NhfOzone                  import NhfOzone
#from NhfPresLevelTemp          import NhfPresLevelTemp
#from NhfSpecSurfHumidity       import NhfSpecSurfHumidity
#from NhfSurfPres               import NhfSurfPres
#from NhfSurfTemp               import NhfSurfTemp
#from NhfWaterVaporMixRatio     import NhfWaterVaporMixRatio
#from NitrateDepletion          import NitrateDepletion
#from PresLevelTemp             import PresLevelTemp
#from SurfPresCorrection        import SurfPresCorrection
#from TropoGeopotentialHeight   import TropoGeopotentialHeight
#from WaterVaporMixRatio        import WaterVaporMixRatio


