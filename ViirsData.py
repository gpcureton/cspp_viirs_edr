#!/usr/bin/env python
# encoding: utf-8
"""
ViirsData.py

Purpose: Provide required data for VIIRS Cloud and Aerosol products.

Created by Geoff Cureton on 2012-11-13.
Copyright (c) 2012 University of Wisconsin SSEC. All rights reserved.
"""

file_Date = '$Date$'
file_Revision = '$Revision$'
file_Author = '$Author$'
file_HeadURL = '$HeadURL$'
file_Id = '$Id$'

__author__ = 'G.P. Cureton <geoff.cureton@ssec.wisc.edu>'
__version__ = '$Id$'
__docformat__ = 'Epytext'


from matplotlib.colors import Colormap, normalize, LinearSegmentedColormap,ListedColormap
from matplotlib import cm as cm
import numpy as np
import time

"""
#############
ViirsData.py
#############

This is a set of band definitions for the VIIRS instrument, for the band names
['M01','M02',...,'M16']. Using these band names as keys, we define the dictionaries...

ViirsBandCenter:

    The center of each of the VIIRS bands (microns)

ViirsBandWidth:

    The center of each of the VIIRS bands (microns)

ViirsBandType:

    The band type of each of the VIIRS bands, which may take either the value
    'reflective' for bands M1 to M11, or 'emissive' for bands M12 to M16.

ViirsBandModisEquiv:

    The equivalant MODIS bands for each of the VIIRS bands, i.e.: from the 
    ipython shell we specify the VIIRS band 'M15'...

        In [13]: ViirsData.ViirsBandModisEquiv['M15']
        Out[13]: '31'

    ... and get the MODIS band '31'.

Example code to plot solar spectra and VIIRS bands...

import ViirsData as VD
import Science.solar as solar
from matplotlib.pylab import plot,legend,xlabel,ylabel,title,broken_barh,xlim,ylim,show,ioff

ioff()
wavelength,Etr,Direct_Circumsolar = solar.ASTM[:,0],solar.ASTM[:,1],solar.ASTM[:,3]
plot(wavelength,Etr,'r',wavelength,Direct_Circumsolar,'g')
legend(("Extraterrestrial","Sea Level"))
xlabel("$\lambda$"+" (nm)")
ylabel("Spectral Irradiance "+"$\mathrm{(W m^{-2} nm^{-1})}$")
title("Viirs Spectral Bands")
refBands = []
emissBands = []
for bands in VD.ViirsBandCenter.keys():
    bandWidth = VD.ViirsBandWidth[bands]*1000.
    bandStart = VD.ViirsBandCenter[bands]*1000. - bandWidth/2.
    if (bands.rfind('M') != -1 and VD.ViirsBandType[bands] == 'reflective'):
        refBands.append((bandStart,bandWidth))
    else :
        emissBands.append((bandStart,bandWidth))

broken_barh(refBands,(0.,2.5),facecolors='lightblue',edgecolors='lightblue')
#broken_barh(emissBands,(0.,2.5),facecolors='green',edgecolors='green')



xlim(250,2500)
ylim(0,2.5)

##############################
Geoff Cureton, SSEC (Oct 2008)
##############################
    
"""


ViirsBandCenter = {}
ViirsBandCenter['DNB'] = 0.700
ViirsBandCenter['I01'] = 0.640
ViirsBandCenter['I02'] = 0.865
ViirsBandCenter['I03'] = 1.610
ViirsBandCenter['I04'] = 3.740
ViirsBandCenter['I05'] = 11.450
ViirsBandCenter['M01'] = 0.412
ViirsBandCenter['M02'] = 0.445
ViirsBandCenter['M03'] = 0.488
ViirsBandCenter['M04'] = 0.555
ViirsBandCenter['M05'] = 0.672
ViirsBandCenter['M06'] = 0.746
ViirsBandCenter['M07'] = 0.865
ViirsBandCenter['M08'] = 1.240
ViirsBandCenter['M09'] = 1.378
ViirsBandCenter['M10'] = 1.610
ViirsBandCenter['M11'] = 2.250
ViirsBandCenter['M12'] = 3.700
ViirsBandCenter['M13'] = 4.050
ViirsBandCenter['M14'] = 8.550
ViirsBandCenter['M15'] = 10.763
ViirsBandCenter['M16'] = 12.0125

ViirsBandWidth = {}
ViirsBandWidth['DNB'] = 0.400
ViirsBandWidth['I01'] = 0.080
ViirsBandWidth['I02'] = 0.039
ViirsBandWidth['I03'] = 0.060
ViirsBandWidth['I04'] = 0.380
ViirsBandWidth['I05'] = 1.900
ViirsBandWidth['M01'] = 0.020
ViirsBandWidth['M02'] = 0.018
ViirsBandWidth['M03'] = 0.020
ViirsBandWidth['M04'] = 0.020
ViirsBandWidth['M05'] = 0.020
ViirsBandWidth['M06'] = 0.015
ViirsBandWidth['M07'] = 0.039
ViirsBandWidth['M08'] = 0.020
ViirsBandWidth['M09'] = 0.015
ViirsBandWidth['M10'] = 0.060
ViirsBandWidth['M11'] = 0.050
ViirsBandWidth['M12'] = 0.180
ViirsBandWidth['M13'] = 0.155
ViirsBandWidth['M14'] = 0.300
ViirsBandWidth['M15'] = 1.000
ViirsBandWidth['M16'] = 0.950

ViirsBandDynamicRadRange = {}
ViirsBandDynamicRadRange['DNB'] = [  3.0E-5, 2.0E2 ]
ViirsBandDynamicRadRange['I01'] = [  5     , 718   ]
ViirsBandDynamicRadRange['I02'] = [  10.3  , 349   ]
ViirsBandDynamicRadRange['I03'] = [  1.2   , 72.5  ]
ViirsBandDynamicRadRange['I04'] = [  0.0016411044 , 4.6132823   ]
ViirsBandDynamicRadRange['I05'] = [  0.0067865545 , 23.057737   ]
ViirsBandDynamicRadRange['M01'] = [  135   , 615   ]
ViirsBandDynamicRadRange['M02'] = [  127   , 687   ]
ViirsBandDynamicRadRange['M03'] = [  107   , 702   ]
ViirsBandDynamicRadRange['M04'] = [  78    , 667   ]
ViirsBandDynamicRadRange['M05'] = [  59    , 651   ]
ViirsBandDynamicRadRange['M06'] = [  5.3   , 41.0  ]
ViirsBandDynamicRadRange['M07'] = [  29    , 349   ]
ViirsBandDynamicRadRange['M08'] = [  3.5   , 164.9 ]
ViirsBandDynamicRadRange['M09'] = [  0.6   , 77.1  ]
ViirsBandDynamicRadRange['M10'] = [  1.2   , 71.2  ]
ViirsBandDynamicRadRange['M11'] = [  0.12  , 31.8  ]
ViirsBandDynamicRadRange['M12'] = [  1.6135e-07 , 3.50999   ]
ViirsBandDynamicRadRange['M13'] = [  0.00214046 , 404.05   ]
ViirsBandDynamicRadRange['M14'] = [  0.0051186468 , 26.182466   ]
ViirsBandDynamicRadRange['M15'] = [ 0.005033703 , 25.330209   ]
ViirsBandDynamicRadRange['M16'] = [  0.00421426 , 21.742956 ]

ViirsBandDynamicBtempRange = {}
ViirsBandDynamicBtempRange['DNB'] = None
ViirsBandDynamicBtempRange['I01'] = None
ViirsBandDynamicBtempRange['I02'] = None
ViirsBandDynamicBtempRange['I03'] = None
ViirsBandDynamicBtempRange['I04'] = [  210   , 353   ]
ViirsBandDynamicBtempRange['I05'] = [  190   , 340   ]
ViirsBandDynamicBtempRange['M01'] = None
ViirsBandDynamicBtempRange['M02'] = None
ViirsBandDynamicBtempRange['M03'] = None
ViirsBandDynamicBtempRange['M04'] = None
ViirsBandDynamicBtempRange['M05'] = None
ViirsBandDynamicBtempRange['M06'] = None
ViirsBandDynamicBtempRange['M07'] = None
ViirsBandDynamicBtempRange['M08'] = None
ViirsBandDynamicBtempRange['M09'] = None
ViirsBandDynamicBtempRange['M10'] = None
ViirsBandDynamicBtempRange['M11'] = None
ViirsBandDynamicBtempRange['M12'] = [  230   , 353   ]
ViirsBandDynamicBtempRange['M13'] = [  343   , 634   ]
ViirsBandDynamicBtempRange['M14'] = [  190   , 336   ]
ViirsBandDynamicBtempRange['M15'] = [  190   , 343   ]
ViirsBandDynamicBtempRange['M16'] = [  190   , 340   ]

ViirsBandType = {}
ViirsBandType['DNB'] = 'radiance'
ViirsBandType['I01'] = 'reflective'
ViirsBandType['I02'] = 'reflective'
ViirsBandType['I03'] = 'reflective'
ViirsBandType['I04'] = 'emissive'
ViirsBandType['I05'] = 'emissive'    
ViirsBandType['M01'] = 'reflective'
ViirsBandType['M02'] = 'reflective'
ViirsBandType['M03'] = 'reflective'
ViirsBandType['M04'] = 'reflective'
ViirsBandType['M05'] = 'reflective'
ViirsBandType['M06'] = 'reflective'
ViirsBandType['M07'] = 'reflective'
ViirsBandType['M08'] = 'reflective'
ViirsBandType['M09'] = 'reflective'
ViirsBandType['M10'] = 'reflective'
ViirsBandType['M11'] = 'reflective'
ViirsBandType['M12'] = 'emissive'
ViirsBandType['M13'] = 'emissive'
ViirsBandType['M14'] = 'emissive'
ViirsBandType['M15'] = 'emissive'
ViirsBandType['M16'] = 'emissive'

# TODO : Add correct Group names for VIIRS geolocation
ViirsGeoGroup = {}
ViirsGeoGroup['GDNB'] = '/All_Data/VIIRS-DNB-GEO_All'
ViirsGeoGroup['GMOD'] = '/All_Data/VIIRS-MOD-GEO_All'
ViirsGeoGroup['GIMG'] = '/All_Data/VIIRS-IMG-GEO_All'
#ViirsGeoGroup['GMTC'] = '/All_Data/VIIRS-MOD-GEO_All'
#ViirsGeoGroup['GITC'] = '/All_Data/VIIRS-IMG-GEO_All'

ViirsBandGroup = {}
ViirsBandGroup['DNB'] = '/All_Data/VIIRS-DNB-SDR_All'
ViirsBandGroup['I01'] = '/All_Data/VIIRS-I1-SDR_All'
ViirsBandGroup['I02'] = '/All_Data/VIIRS-I2-SDR_All'
ViirsBandGroup['I03'] = '/All_Data/VIIRS-I3-SDR_All'
ViirsBandGroup['I04'] = '/All_Data/VIIRS-I4-SDR_All'
ViirsBandGroup['I05'] = '/All_Data/VIIRS-I5-SDR_All'
ViirsBandGroup['M01'] = '/All_Data/VIIRS-M1-SDR_All'
ViirsBandGroup['M02'] = '/All_Data/VIIRS-M2-SDR_All'
ViirsBandGroup['M03'] = '/All_Data/VIIRS-M3-SDR_All'
ViirsBandGroup['M04'] = '/All_Data/VIIRS-M4-SDR_All'
ViirsBandGroup['M05'] = '/All_Data/VIIRS-M5-SDR_All'
ViirsBandGroup['M06'] = '/All_Data/VIIRS-M6-SDR_All'
ViirsBandGroup['M07'] = '/All_Data/VIIRS-M7-SDR_All'
ViirsBandGroup['M08'] = '/All_Data/VIIRS-M8-SDR_All'
ViirsBandGroup['M09'] = '/All_Data/VIIRS-M9-SDR_All'
ViirsBandGroup['M10'] = '/All_Data/VIIRS-M10-SDR_All'
ViirsBandGroup['M11'] = '/All_Data/VIIRS-M11-SDR_All'
ViirsBandGroup['M12'] = '/All_Data/VIIRS-M12-SDR_All'
ViirsBandGroup['M13'] = '/All_Data/VIIRS-M13-SDR_All'
ViirsBandGroup['M14'] = '/All_Data/VIIRS-M14-SDR_All'
ViirsBandGroup['M15'] = '/All_Data/VIIRS-M15-SDR_All'
ViirsBandGroup['M16'] = '/All_Data/VIIRS-M16-SDR_All'

ViirsBandModisEquiv = {}
ViirsBandModisEquiv['M01'] = '8'     
ViirsBandModisEquiv['M02'] = '9'     
ViirsBandModisEquiv['M03'] = '10'     
ViirsBandModisEquiv['M04'] = '12'     
ViirsBandModisEquiv['M05'] = '13'     
ViirsBandModisEquiv['M06'] = '15'     
ViirsBandModisEquiv['M07'] = '16'     
ViirsBandModisEquiv['M08'] = '5'  
ViirsBandModisEquiv['M09'] = '26'
ViirsBandModisEquiv['M10'] = '6'  
ViirsBandModisEquiv['M11'] = '7'  
ViirsBandModisEquiv['M12'] = '20'
ViirsBandModisEquiv['M13'] = '23'
ViirsBandModisEquiv['M14'] = '29'
ViirsBandModisEquiv['M15'] = '31'
ViirsBandModisEquiv['M16'] = '32'

class ViirsTrimTable:
    """
    This class defines the VIIRS bow-tie deletion pixels, and has methods to
    generate instantiations of the trim table.
    """

    def __init__(self):
        """
            Class init defines the two lists, modTrimTable and imgTrimTable.
            These lists define the start and end columns of the non-trimmed pixels
            for each row in the moderate (16 rows) and imagery (32 rows) scans
            respectively.
        """

        self.modTrimTable = [[1090, 2109, 1008, 2191],
                             [820 , 2379, 640 , 2559],
                             [520 , 2679, 0   , 3199],
                             [130 , 3069, 0   , 3199],
                             [0   , 3199, 0   , 3199],
                             [0   , 3199, 0   , 3199],
                             [0   , 3199, 0   , 3199],
                             [0   , 3199, 0   , 3199],
                             [0   , 3199, 0   , 3199],
                             [0   , 3199, 0   , 3199],
                             [0   , 3199, 0   , 3199],
                             [0   , 3199, 0   , 3199],
                             [130 , 3069, 0   , 3199],
                             [520 , 2679, 0   , 3199],
                             [820 , 2379, 640 , 2559],
                             [1090, 2109, 1008, 2191]]

        
        self.imgTrimTable = [[2180, 4219, 2016, 4383],
                             [2180, 4219, 2016, 4383],
                             [1640, 4759, 1280, 5119],
                             [1640, 4759, 1280, 5119],
                             [1040, 5359, 0   , 6399],
                             [1040, 5359, 0   , 6399],
                             [260 , 6139, 0   , 6399],
                             [260 , 6139, 0   , 6399],
                             [0   , 6399, 0   , 6399],
                             [0   , 6399, 0   , 6399],
                             [0   , 6399, 0   , 6399],
                             [0   , 6399, 0   , 6399],
                             [0   , 6399, 0   , 6399],
                             [0   , 6399, 0   , 6399],
                             [0   , 6399, 0   , 6399],
                             [0   , 6399, 0   , 6399],
                             [0   , 6399, 0   , 6399],
                             [0   , 6399, 0   , 6399],
                             [0   , 6399, 0   , 6399],
                             [0   , 6399, 0   , 6399],
                             [0   , 6399, 0   , 6399],
                             [0   , 6399, 0   , 6399],
                             [0   , 6399, 0   , 6399],
                             [0   , 6399, 0   , 6399],
                             [260 , 6139, 0   , 6399],
                             [260 , 6139, 0   , 6399],
                             [1040, 5359, 0   , 6399],
                             [1040, 5359, 0   , 6399],
                             [1640, 4759, 1280, 5119],
                             [1640, 4759, 1280, 5119],
                             [2180, 4219, 2016, 4383],
                             [2180, 4219, 2016, 4383]]

        ### Some fill values to be used with the pixel trim, depending
        ### on the datatype of the dataset.
        ### (IDPS header: ProCmnDefs.h)
        
        # CDFDB-X_Volume1
        # Table 3.5.5-1, Pixel Level Fill Values

        self.sdrTypeFill = { 
            'NA_FILL' : {
                    ## Algorithm exclusions, N/A
                    ## The pixel was not to be computed because it is not aplicable to situation
                    'float64'   : -999.9e0             ,       # NA_FLOAT64_FILL     
                    'int64'     : -999                 ,       # NA_INT64_FILL       
                    'uint64'    : 18446744073709551615 ,       # NA_UINT64_FILL      
                    'float32'   : -999.9               ,       # NA_FLOAT32_FILL     
                    'int32'     : -999                 ,       # NA_INT32_FILL       
                    'uint32'    : 424967295            ,       # NA_UINT32_FILL      
                    'int16'     : -999                 ,       # NA_INT16_FILL       
                    'uint16'    : 65535                ,       # NA_UINT16_FILL      
                    'int8'      : 127                  ,       # NA_INT8_FILL        
                    'uint8'     : 255                          # NA_UINT8_FILL       
                },
            'MISS_FILL': {
                    ##Missing - C3S provided a fill value or AP missing
                    'float64'   : -999.8e0             ,     # MISS_FLOAT64_FILL   
                    'int64'     : -998                 ,     # MISS_INT64_FILL     
                    'uint64'    : 18446744073709551614 ,     # MISS_UINT64_FILL    
                    'float32'   : -999.8               ,     # MISS_FLOAT32_FILL   
                    'int32'     : -998                 ,     # MISS_INT32_FILL     
                    'uint32'    : 424967294            ,     # MISS_UINT32_FILL    
                    'int16'     : -998                 ,     # MISS_INT16_FILL     
                    'uint16'    : 65534                ,     # MISS_UINT16_FILL    
                    'int8'      : 126                  ,     # MISS_INT8_FILL      
                    'uint8'     : 254                        # MISS_UINT8_FILL     
                },
            'ONBOARD_PT_FILL': {
                    ## Onboard pixel trim 
                    ## The VIIRS pixel was trimmed on the S/C
                    'float64'   : -999.7e0             ,     # ONBOARD_PT_FLOAT64_FILL   
                    'int64'     : -997                 ,     # ONBOARD_PT_INT64_FILL     
                    'uint64'    : 18446744073709551613 ,     # ONBOARD_PT_UINT64_FILL    
                    'float32'   : -999.7               ,     # ONBOARD_PT_FLOAT32_FILL   
                    'int32'     : -997                 ,     # ONBOARD_PT_INT32_FILL     
                    'uint32'    : 424967293            ,     # ONBOARD_PT_UINT32_FILL    
                    'int16'     : -997                 ,     # ONBOARD_PT_INT16_FILL     
                    'uint16'    : 65533                ,     # ONBOARD_PT_UINT16_FILL    
                    'int8'      : 125                  ,     # ONBOARD_PT_INT8_FILL      
                    'uint8'     : 253                        # ONBOARD_PT_UINT8_FILL     
                },
            'ONGROUND_PT_FILL': {
                    ## On-ground pixel trim
                    ## The VIIRS pixel was trimmed during processing.
                    'float64'   : -999.6e0             ,     # ONGROUND_PT_FLOAT64_FILL  
                    'int64'     : -996                 ,     # ONGROUND_PT_INT64_FILL    
                    'uint64'    : 18446744073709551612 ,     # ONGROUND_PT_UINT64_FILL   
                    'float32'   : -999.6               ,     # ONGROUND_PT_FLOAT32_FILL  
                    'int32'     : -996                 ,     # ONGROUND_PT_INT32_FILL    
                    'uint32'    : 424967292            ,     # ONGROUND_PT_UINT32_FILL   
                    'int16'     : -996                 ,     # ONGROUND_PT_INT16_FILL    
                    'uint16'    : 65532                ,     # ONGROUND_PT_UINT16_FILL   
                    'int8'      : 124                  ,     # ONGROUND_PT_INT8_FILL     
                    'uint8'     : 252                        # ONGROUND_PT_UINT8_FILL    
                },
            'ERR_FILL': {
                    ## Cannot calculate 
                    ## Algorithm could not compute the pixel because of a software problem 
                    ## ie. couldn't converge
                    'float64'   : -999.5e0             ,     # ERR_FLOAT64_FILL          
                    'int64'     : -995                 ,     # ERR_INT64_FILL            
                    'uint64'    : 18446744073709551611 ,     # ERR_UINT64_FILL           
                    'float32'   : -999.5               ,     # ERR_FLOAT32_FILL          
                    'int32'     : -995                 ,     # ERR_INT32_FILL            
                    'uint32'    : 424967291            ,     # ERR_UINT32_FILL           
                    'int16'     : -995                 ,     # ERR_INT16_FILL            
                    'uint16'    : 65531                ,     # ERR_UINT16_FILL           
                    'int8'      : 123                  ,     # ERR_INT8_FILL             
                    'uint8'     : 251                        # ERR_UINT8_FILL            
                },
            'VDNE_FLOAT64_FILL': {
                    ## Value Does Not Exist/Missing Scan
                    'float64'   : -999.3e0             ,     # VDNE_FLOAT64_FILL          
                    'int64'     : -993                 ,     # VDNE_INT64_FILL            
                    'uint64'    : 18446744073709551609 ,     # VDNE_UINT64_FILL           
                    'float32'   : -999.3               ,     # VDNE_FLOAT32_FILL          
                    'int32'     : -993                 ,     # VDNE_INT32_FILL            
                    'uint32'    : 424967289            ,     # VDNE_UINT32_FILL           
                    'int16'     : -993                 ,     # VDNE_INT16_FILL            
                    'uint16'    : 65529                ,     # VDNE_UINT16_FILL           
                    'int8'      : 121                  ,     # VDNE_INT8_FILL             
                    'uint8'     : 249                        # VDNE_UINT8_FILL            
                },
            'SOUB_FLOAT64_FILL': {
                    ##Scale Out of Bounds
                    'float64'   : -999.2e0             ,     # SOUB_FLOAT64_FILL          
                    'int64'     : -992                 ,     # SOUB_INT64_FILL            
                    'uint64'    : 18446744073709551608 ,     # SOUB_UINT64_FILL           
                    'float32'   : -999.2               ,     # SOUB_FLOAT32_FILL          
                    'int32'     : -992                 ,     # SOUB_INT32_FILL            
                    'uint32'    : 424967288            ,     # SOUB_UINT32_FILL           
                    'int16'     : -992                 ,     # SOUB_INT16_FILL            
                    'uint16'    : 65528                ,     # SOUB_UINT16_FILL           
                    'int8'      : 120                  ,     # SOUB_INT8_FILL             
                    'uint8'     : 248                        # SOUB_UINT8_FILL
                },
        }

    def createModTrimArray(self,nscans=48,trimType=bool):
        """
            Creates an array with nDetectors*nscans pixel rows, with the trimmed
            pixels set to True.
        """
        nDetectors = 16
        trimScanArray = np.ones((nDetectors,3200),dtype=trimType)
        for row in range(len(self.modTrimTable)):
            colStart = self.modTrimTable[row][0]
            colEnd = self.modTrimTable[row][1] + 1
            trimScanArray[row,colStart:colEnd] = False

        trimArray = np.ones((nDetectors*nscans,3200),dtype=trimType)
        for scan in range(nscans):
            startRow = nDetectors * scan
            endRow = startRow + nDetectors
            trimArray[startRow:endRow,:] = trimScanArray

        return trimArray

    def createImgTrimArray(self,nscans=48,trimType=bool):
        """
            Creates an array with nDetectors*nscans pixel rows, with the trimmed
            pixels set to True.
        """
        nDetectors = 32
        trimScanArray = np.ones((nDetectors,6400),dtype=trimType)
        for row in range(len(self.imgTrimTable)):
            colStart = self.imgTrimTable[row][0]
            colEnd = self.imgTrimTable[row][1] + 1
            trimScanArray[row,colStart:colEnd] = False

        trimArray = np.ones((nDetectors*nscans,6400),dtype=trimType)
        for scan in range(nscans):
            startRow = nDetectors * scan
            endRow = startRow + nDetectors
            trimArray[startRow:endRow,:] = trimScanArray

        return trimArray

    def createOnboardModTrimArray(self,nscans=48,trimType=bool):
        """
            Creates an array with nDetectors*nscans pixel rows, with the trimmed
            pixels set to True.
        """
        nDetectors = 16
        trimScanArray = np.ones((nDetectors,3200),dtype=trimType)
        for row in range(len(self.modTrimTable)):
            colStart = self.modTrimTable[row][2]
            colEnd = self.modTrimTable[row][3] + 1
            trimScanArray[row,colStart:colEnd] = False

        trimArray = np.ones((nDetectors*nscans,3200),dtype=trimType)
        for scan in range(nscans):
            startRow = nDetectors * scan
            endRow = startRow + nDetectors
            trimArray[startRow:endRow,:] = trimScanArray

        return trimArray

    def createOnboardImgTrimArray(self,nscans=48,trimType=bool):
        """
            Creates an array with nDetectors*nscans pixel rows, with the trimmed
            pixels set to True.
        """
        nDetectors = 32
        trimScanArray = np.ones((nDetectors,6400),dtype=trimType)
        for row in range(len(self.imgTrimTable)):
            colStart = self.imgTrimTable[row][2]
            colEnd = self.imgTrimTable[row][3] + 1
            trimScanArray[row,colStart:colEnd] = False

        trimArray = np.ones((nDetectors*nscans,6400),dtype=trimType)
        for scan in range(nscans):
            startRow = nDetectors * scan
            endRow = startRow + nDetectors
            trimArray[startRow:endRow,:] = trimScanArray

        return trimArray


class CloudMaskData:
    """
    This class contains static data for the interpretation of the 
    VIIRS Cloud Mask.
    """

    ViirsCMbitMasks = ((3,12,16,32,192),
                       (7,8,16,32,64,128),
                       (1,2,4,8,16,32,64,128),
                       (3,4,8,16,32,64,128),
                       (256,),
                       (7,8,16,32,64,128))

    ViirsCMbitShift = ((0,2,4,5,6),
                       (0,3,4,5,6,7),
                       (0,1,2,3,4,5,6,7),
                       (0,2,3,4,5,6,7),
                       (0,),
                       (0,3,4,5,6,7))

    ViirsCMbitMaskNames = (('CM Quality','Cloud Mask','Day/Night','Snow/Ice','Sunglint'),
                           ('Land/Water','Shadow','Heavy Aerosol','Fire','Cirrus (RM9)','Cirrus (BTM15-BTM16)'),
                           ('IR Threshold Cloud Test (BTM15)','High Cloud (BTM12-BTM16) Test','IR Temperature Difference  Test (BTM14 - BTM15 BTM15 - BTM16)','Temperature Difference Test (BTM15 - BTM12)','Temperature Difference Test (BTM12 - BTM13)','Visible Reflectance Test (RM5)','Visible Reflectance Test (RM7)','Visible Ratio Test (RM7 / RM5)'),
                           ('Adjacent Pixel Cloud Confident Value','Conifer Boreal Forest','Spatial Uniformity','Dust','Smoke','Dust/Vol. Ash','Spare'),
                           ('Spare',),
                           ('Cloud Phase','Thin Cirrus','Ephemeral Water','Degraded TOC NDVI Flag','Degraded Sun Glint Flag','Degraded Polar Night Flag' ))
    
    ViirsCMvalues = ((  (0,1,2,3),(0,1,2,3),(0,1),(0,1),(0,1,2,3)   ),           # Byte 0
                     (  (0,1,2,3,4,5),(0,1),(0,1),(0,1),(0,1),(0,1)        ),      # Byte 1
                     (  (0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1),(0,1)  ),      # Byte 2
                     (  (0,1,2,3),(0,1),(0,1),(0,1),(0,1),(0,1),(0)      ),      # Byte 3
                     (  (0),                                             ),      # Byte 4
                     #(  (2,3,4,5,6,7),(0,1),(0,1),(0)                    ))      # Byte 5
                     (  (0,1,2,3,4,5,6,7),(0,1),(0,1),(0,1),(0,1),(0,1)   ))      # Byte 5

    ViirsCMfillBoundaries = [
                             [ [-0.5,0.5,1.5,2.5,3.5],[-0.5,0.5,1.5,2.5,3.5],[-0.5,0.5,1.5],[-0.5,0.5,1.5],[-0.5,0.5,1.5,2.5,3.5] ], # Byte 0
                             [ [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5],[-0.5,0.5,1.5],[-0.5,0.5,1.5],[-0.5,0.5,1.5],[-0.5,0.5,1.5],[-0.5,0.5,1.5] ], # Byte 1
                             [ [-0.5,0.5,1.5],[-0.5,0.5,1.5],[-0.5,0.5,1.5],[-0.5,0.5,1.5],[-0.5,0.5,1.5],[-0.5,0.5,1.5],[-0.5,0.5,1.5],[-0.5,0.5,1.5] ], # Byte 2
                             [ [-0.5,0.5,1.5,2.5,3.5],[-0.5,0.5,1.5],[-0.5,0.5,1.5],[-0.5,0.5,1.5],[-0.5,0.5,1.5],[-0.5,0.5,1.5],[0] ], # Byte 3
                             [ [0], ], # Byte 4
                             #[ [2.5,3.5,4.5,5.5,6.5,7.5],[-0.5,0.5,1.5],[-0.5,0.5,1.5],[0] ] # Byte 5
                             [ [-0.5,0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5],[-0.5,0.5,1.5],[-0.5,0.5,1.5],[-0.5,0.5,1.5],[-0.5,0.5,1.5],[-0.5,0.5,1.5] ] # Byte 5
                            ]

    ViirsCMfillColours = [
                            [
                                ['k','r','b','#00ff00'],['#00ff00','#00ffff','#ff0000','w'],['k','yellow'],['k','g'],['#000080','#0012ff','#00a4ff','#b6ff41'] ], # Byte 0
                            [   ['brown','green','cyan','blue','w','yellow'],['k','g'],['k','g'],['k','g'],['k','g'],['k','g']   ], # Byte 1
                            [   ['k','g'],['k','g'],['k','g'],['k','g'],['k','g'],['k','g'],['k','g'],['k','g']  ], # Byte 2
                            [   ['#00ff00','#00ffff','#ff0000','w'],['k','g'],['k','g'],['k','g'],['k','g'],['k','g'],['k','g'] ], # Byte 3
                            [   [0], ], # Byte 4
                            #[   ['#00A5FF','#00A5FF','#00A5FF','#FFB900','#FFB900','#FF3000'],[0],[0],[0]   ]   # Byte 5
                            [   ['#000080','#0011FF','#00A5FF','#00A5FF','#B7FF40','#FFB900','#FF3000','#FF3000'],['k','g'],['k','g'],['k','g'],['k','g'],['k','g']   ]   # Byte 5
                         ]

    # The tick labels of the colourbar categories
    ViirsCMtickNames = [[['poor','low','medium','high'],['confident clear','probably clear','probably cloudy','confident cloudy'],['night','day'],['no snow','snow/ice'],['none','geometry','wind','geometry/wind']],      # Byte 0
                     [  ['Land/Desert','Land','Inland Water','Sea Water','None','Coastal'],['No','Yes'],['No','Yes'],['No','Yes'],['No Cloud','Cloud'],['No Cloud','Cloud']],      # Byte 1
                     [  ['No Cloud','Cloud'],['No Cloud','Cloud'],['No Cloud','Cloud'],['No Cloud','Cloud'],['No Cloud','Cloud'],['No Cloud','Cloud'],['No Cloud','Cloud'],['No Cloud','Cloud']  ],      # Byte 2
                     [  ['confident clear','probably clear','probably cloudy','confident cloudy'],['No','Yes'],['No','Yes'],['No','Yes'],['No','Yes'],['No','Yes'],['No','Yes']],      # Byte 3
                     [  [0],],      # Byte 4
                     #[  ['Part Cld','Water','Mixed','OP Ice','Cirrus','Overlap'],[0],[0],[0]]]      # Byte 5
                     [  ['Not Exe','Clear','Part Cld','Water','Mixed','OP Ice','Cirrus','Overlap'],['No','Yes'],['No','Yes'],['No','Yes'],['No','Yes'],['No','Yes']      ]]      # Byte 5

    # Each of the rows in each of red,green,blue triples defines a boundary. There
    # are boundaries at each end (at 0 and 1 points of colormap), and a boundary
    # between each category. For 4 categories, this gives up 5 boundaries. Each row
    # gives the point along the colorscale of the boundary, and the colour value
    # above and below the boundary, for that Red, Blue or Green tuple.

    # This is a segmented colourmap for the four-category
    # cloudmask product

    cdict_CloudMask = {'red': (
                                (0.0 , 0.0, 0.0),       # 1 max boundary
                                (0.25, 0.0, 0.0),       # 2 category boundary
                                (0.5 , 0.0, 1.0),       # 3 category boundary
                                (0.75, 1.0, 1.0),       # 4 category boundary
                                (1.0 , 1.0, 1.0)        # 5 min boundary
                              ),
                    'green': (
                                (0.0 , 1.0, 1.0),       # 1 max boundary
                                (0.25, 1.0, 1.0),       # 2 category boundary
                                (0.5 , 1.0, 0.0),       # 3 category boundary
                                (0.75, 0.0, 1.0),       # 4 category boundary
                                (1.0 , 1.0, 1.0)        # 5 min boundary
                              ),
                    'blue': (
                                (0.0 , 0.0, 0.0),       # 1 max boundary
                                (0.25, 0.0, 1.0),       # 2 category boundary
                                (0.5 , 1.0, 0.0),       # 3 category boundary
                                (0.75, 0.0, 1.0),       # 4 category boundary
                                (1.0 , 1.0, 1.0)        # 5 min boundary
                            )
                        }

    #cmap_CloudMask = LinearSegmentedColormap('colormap',cdict_CloudMask,1024)

    # This is a segmented colourmap the the two-category
    # cloudmask product
    cdict_binary =  {'red': (   (0.0, 1.0, 1.0),        # min boundary
                                (0.5, 1.0, 1.0),        # category boundary
                                (1.0, 1.0, 0.0)),       # max boundary
                    'green': (  (0.0, 1.0, 1.0),
                                (0.5, 1.0, 0.0),
                                (1.0, 0.0, 1.0)),
                    'blue': (   (0.0, 1.0, 1.0),
                                (0.5, 1.0, 0.0),
                                (1.0, 0.0, 1.0))}

    #cmap_binary = LinearSegmentedColormap('colormap',cdict_binary,1024)


class CloudProdData:
    """
    This class contains static data for the VIIRS Cloud Products.
    """

    class CloudProd:
        """
        This class provides a data structure to hold the plot attributes
        for a particular cloud product.
        """
        def __init__(self,
            SDSname         = None,\
            SDSphaseName    = None,\
            SDStitleString  = None,\
            vmin            = None,\
            vmax            = None,\
            cmap            = None,\
            cbarTickPos     = None,\
            cbarTickNames   = None,\
            cbarTitle       = None,\
            vmin_water      = None,\
            vmax_water      = None,\
            dv_water        = None,\
            cbarTitle_water = None,\
            cmap_water      = None,\
            vmin_ice        = None,\
            vmax_ice        = None,\
            dv_ice          = None,\
            cbarTitle_ice   = None,\
            cmap_ice        = None,\
            logScale        = False
            ):
            self.SDSname         = SDSname        
            self.SDSphaseName    = SDSphaseName        
            self.SDStitleString  = SDStitleString 
            self.vmin            = vmin           
            self.vmax            = vmax           
            self.cmap            = cmap           
            self.cbarTickPos     = cbarTickPos    
            self.cbarTickNames   = cbarTickNames  
            self.cbarTitle       = cbarTitle      
            self.vmin_water      = vmin_water     
            self.vmax_water      = vmax_water     
            self.dv_water        = dv_water       
            self.cbarTitle_water = cbarTitle_water
            self.cmap_water      = cmap_water     
            self.vmin_ice        = vmin_ice       
            self.vmax_ice        = vmax_ice       
            self.dv_ice          = dv_ice         
            self.cbarTitle_ice   = cbarTitle_ice  
            self.cmap_ice        = cmap_ice
            self.logScale        = logScale

            # These are the water/ice colormaps that are closest to what is used
            # for MODIS in LAADS.
            ice_cmdata = {
                'red'  :  (
                             (0.0 , 0.72, 0.72), 
                             (0.67, 0.0 , 0.0 ), 
                             (1.0 , 0.0 , 0.0 )
                          ),
                'green':  (
                             (0.0 , 0.0 , 0.0 ),
                             (0.358, 0.0 , 0.0 ),
                             (1.0 , 1.0 , 1.0 )
                          ),
                'blue' :  (
                             (0.0 , 0.72, 0.72),
                             (0.35, 0.56 , 0.56 ),
                             (0.67, 1.0 , 1.0 ),
                             (1.0 , 0.0 , 0.0 )
                          )
            }
            #cmap_ice = LinearSegmentedColormap('ice_cm',ice_cmdata,256)
            #self.cmap_ice = cmap_ice

            # this colormap goes from yellow to white-ish orange and to red
            water_cmdata = {
                'red'  :  (
                             (0.0 , 1.0 , 1.0 ), 
                             (0.82 , 1.0 , 1.0 ),
                             (1.0  ,0.5 , 0.5 )
                          ),
                'green':  (
                             (0.0 , 1.0 , 1.0 ), 
                             (0.325, 1.0 , 1.0 ),
                             (0.68 , 0.5 , 0.5 ),
                             (0.816 , 0.0 , 0.0 ),
                             (1.0 , 0.0 , 0.0 )
                          ),
                'blue' :  (
                             (0.0 , 0.0 , 0.0 ), 
                             (0.358, 0.74, 0.74),
                             (0.667, 0.0 , 0.0 ),
                             (1.0 , 0.0 , 0.0 )
                          )
            }
            #cmap_water = LinearSegmentedColormap('water_cm',water_cmdata,256)
            #self.cmap_water = cmap_water     

            # this colormap is a concatenation of the previous two
            ice_water_cmdata = {
                'red'  :  (
                             (0.0   , 0.67, 0.67), 
                             (0.193 ,0.0 , 0.0 ), 
                             (0.445 , 0.0 , 0.0 ),
                             (0.536 , 1.0 , 1.0 ), 
                             (0.912 , 1.0 , 1.0 ),
                             (1.0  ,0.697 , 0.697 )
                          ),
                'green':  (
                             (0.0 , 0.0 , 0.0 ),
                             (0.139, 0.0 , 0.0 ),
                             (0.507 , 1.0 , 1.0 ),
                             (0.605, 1.0 , 1.0 ),
                             (1.0 , 0.0 , 0.0 )
                          ),
                'blue' :  (
                             (0.0 , 0.67, 0.67),
                             (0.234, 1.0 , 1.0 ),
                             (0.256, 0.92 , 0.92 ),
                             (0.298 , 0.942 , 0.942 ),
                             (0.489 , 0.0 , 0.0 ), 
                             (1.0 , 0.0 , 0.0 )
                          )
            }
            cmap_ice_water = LinearSegmentedColormap('ice_water_cm',ice_water_cmdata,256)
            self.cmap_ice_water = cmap_ice_water


    CloudProduct = {
        'cphase_cop': CloudProd(
            # mask with 224, shift by 5
            #SDSphaseName = '/All_Data/VIIRS-Cd-Opt-Prop-IP_All/copQF0'  # mini-IDPS 1.5.0.48
            SDSphaseName = '/All_Data/VIIRS-Cd-Opt-Prop-IP_All/QF1_VIIRSCOPIP',
        ),
        'cphase_ctp': CloudProd(
            # mask with 7
            #SDSphaseName = '/All_Data/VIIRS-Cd-Top-Parm-IP_All/ctpQF1'  # mini-IDPS 1.5.0.48
            SDSphaseName = '/All_Data/VIIRS-Cd-Top-Parm-IP_All/QF2_VIIRSCTPIP',
        ),
        'ctp': CloudProd(
            SDSname = '/All_Data/VIIRS-Cd-Top-Parm-IP_All/ctp',
            SDStitleString = "VIIRS Cloud Top\nPressure",
            #SDSphaseName = '/All_Data/VIIRS-Cd-Top-Parm-IP_All/ctpQF1',  # mini-IDPS 1.5.0.48
            SDSphaseName = '/All_Data/VIIRS-Cd-Top-Parm-IP_All/QF2_VIIRSCTPIP',
            vmin_water = 300.,
            vmax_water = 1200.,
            dv_water = 300.0 ,
            cbarTitle_water = "$ctp_{water} (\mathrm{hPa})$",
            #cmap_water = cm.hot_r,
            vmin_ice = 0.,
            vmax_ice = 300.,
            dv_ice = 50.0 ,
            cbarTitle_ice = '$ctp_{ice} (\mathrm{hPa}) $',
            #cmap_ice = cm.winter
        ),
        'ctt': CloudProd(
            SDSname = '/All_Data/VIIRS-Cd-Top-Parm-IP_All/ctt',
            SDStitleString = "VIIRS Cloud Top\nTemperature",
            #SDSphaseName = '/All_Data/VIIRS-Cd-Top-Parm-IP_All/ctpQF1',  # mini-IDPS 1.5.0.48
            SDSphaseName = '/All_Data/VIIRS-Cd-Top-Parm-IP_All/QF2_VIIRSCTPIP',
            vmin_water = 180.,
            vmax_water = 320.,
            dv_water = 20.0,
            cbarTitle_water = r'$ctt_{water} (\mathrm{K})$',
            #cmap_water = cm.hot_r,
            vmin_ice = 180.,
            vmax_ice = 320.,
            dv_ice = 20.0,
            cbarTitle_ice = r'$ctt_{ice} (\mathrm{K})$',
            #cmap_ice = cm.winter
        ),
        'cth': CloudProd(
            SDSname = '/All_Data/VIIRS-Cd-Top-Parm-IP_All/cth',
            SDStitleString = "VIIRS Cloud Top\nHeight",
            #SDSphaseName = '/All_Data/VIIRS-Cd-Top-Parm-IP_All/ctpQF1',  # mini-IDPS 1.5.0.48
            SDSphaseName = '/All_Data/VIIRS-Cd-Top-Parm-IP_All/QF2_VIIRSCTPIP',
            vmin_water = 0.,
            vmax_water = 20.,
            dv_water = 2.0,
            cbarTitle_water = r'$cth_{water} (\mathrm{km})$',
            #cmap_water = cm.hot_r,
            vmin_ice = 0.,
            vmax_ice = 20.,
            dv_ice = 2.0,
            cbarTitle_ice = r'$cth_{ice} (\mathrm{km})$',
            #cmap_ice = cm.winter
        ),
        'cot':  {
            'log': CloudProd(
                SDSname = '/All_Data/VIIRS-Cd-Opt-Prop-IP_All/cot',
                SDStitleString = "VIIRS Cloud\nOptical Thickness",
                #SDSphaseName = '/All_Data/VIIRS-Cd-Opt-Prop-IP_All/copQF0',  # mini-IDPS 1.5.0.48
                SDSphaseName = '/All_Data/VIIRS-Cd-Opt-Prop-IP_All/QF1_VIIRSCOPIP',
                vmin_water = np.log10(0.1),
                vmax_water = np.log10(100.),
                dv_water = 0.5,
                #cbarTitle_water = r"$\log_{10}(\tau_{water})$",
                cbarTitle_water = r"$\tau_{water}$",
                #cmap_water = cm.hot_r,
                vmin_ice = np.log10(0.1),
                vmax_ice = np.log10(100.),
                dv_ice = 0.5,
                cbarTitle_ice = r"$\tau_{ice}$",
                #cmap_ice = cm.winter,
                logScale = True
            ),
            'linear' : CloudProd(
                SDSname = '/All_Data/VIIRS-Cd-Opt-Prop-IP_All/cot',
                SDStitleString = "VIIRS Cloud\nOptical Thickness",
                #SDSphaseName = '/All_Data/VIIRS-Cd-Opt-Prop-IP_All/copQF0',  # mini-IDPS 1.5.0.48
                SDSphaseName = '/All_Data/VIIRS-Cd-Opt-Prop-IP_All/QF1_VIIRSCOPIP',
                vmin_water = 0.,
                vmax_water = 100.,
                dv_water = 20.0,
                cbarTitle_water = r"$\tau_{water}$",
                #cmap_water = cm.hot_r,
                vmin_ice = 0.,
                vmax_ice = 100.,
                dv_ice = 20.0,
                cbarTitle_ice = r'$\tau_{ice}$',
                #cmap_ice = cm.winter
            )
        },
        'eps': {
            'log' : CloudProd(
                SDSname = '/All_Data/VIIRS-Cd-Opt-Prop-IP_All/eps',
                SDStitleString = "VIIRS Cloud Effective\nParticle Radius",
                #SDSphaseName = '/All_Data/VIIRS-Cd-Opt-Prop-IP_All/copQF0',  # mini-IDPS 1.5.0.48
                SDSphaseName = '/All_Data/VIIRS-Cd-Opt-Prop-IP_All/QF1_VIIRSCOPIP',
                vmin_water = np.log10(1),
                vmax_water = np.log10(100.),
                dv_water = 0.5,
                cbarTitle_water = r"$r_{water}$ $(\mu\mathrm{m})$",
                #cmap_water = cm.hot_r,
                vmin_ice = np.log10(1),
                vmax_ice = np.log10(100.),
                dv_ice = 0.5,
                cbarTitle_ice = r'$r_{ice}$ $(\mu\mathrm{m})$',
                #cmap_ice = cm.winter,
                logScale = True
            ),
            'linear' : CloudProd(
                SDSname = '/All_Data/VIIRS-Cd-Opt-Prop-IP_All/eps',
                SDStitleString = "VIIRS Cloud Effective\nParticle Radius",
                #SDSphaseName = '/All_Data/VIIRS-Cd-Opt-Prop-IP_All/copQF0',  # mini-IDPS 1.5.0.48
                SDSphaseName = '/All_Data/VIIRS-Cd-Opt-Prop-IP_All/QF1_VIIRSCOPIP',
                vmin_water = 0.,
                vmax_water = 100.,
                dv_water = 10.0,
                cbarTitle_water = r"$r_{water}$ $(\mu\mathrm{m})$",
                #cmap_water = cm.hot_r,
                vmin_ice = 0.,
                vmax_ice = 60.,
                dv_ice = 20.0,
                cbarTitle_ice = r'$r_{ice}$ $(\mu\mathrm{m})$',
                #cmap_ice = cm.winter
            )
        }
    }


class AerosolProdData:
    """
    This class contains static data for the VIIRS Aerosol Products.
    """

    ViirsAeroQualBitMasks = ((3, 12,  48, 192             ),
                             (3, 12, 112, 128             ),
                             (3, 12, 112                  ),
                             (1,  2,   4,   8, 48, 64, 128),
                             (1,  2,   4,   8))

    ViirsAeroQualBitShift = ((0,  2,   4,   6             ),
                             (0,  2,   4,   7             ),
                             (0,  2,   4                  ),
                             (0,  1,   2,   3,  4,  6,   7),
                             (0,  1,   2,   3))

    ViirsAeroQualBitMaskNames = (('Aerosol Optical Thickness Quality','Angstrom Exponent Quality','Suspended Matter Type Quality','Cloud Mask Quality'),
                                 ('Cloud Detection Result & Confidence Indicator','Adjacent Pixel Cloud Confidence Value','Land/Water Background','Bad SDR'),
                                 ('Day/Night Flag','Interpolation/NAAPS/Climatology Processing','Sun Glint'),
                                 ('Snow / Ice', 'Cirrus', 'Cloud Shadow', 'Fire', 'Bright Pixel', 'Turbid / Shallow Water', 'Volcanic Ash'),
                                 ('Low AOT, SM Typing Excluded', 'Low AOT, SM Detection Excluded', 'AOT or APSP Out of Spec Range', 'Residual Threshold Exceeded')
                                )
    
    ViirsAeroQualvalues = ((  (0,1,2,3), (0,1,2,3), (0,1,2,3), (0,1,2,3)                ),      # Byte 0
                           (  (0,1,2,3), (0,1,2,3), (0,1,2,3,4,5,6), (0,1)              ),      # Byte 1
                           (  (0,1,2,3), (0,1,2,3), (0,1,2,3,4,5,6,7)                   ),      # Byte 2
                           (  (0,1), (0,1), (0,1), (0,1), (0,1,2), (0,1), (0,1)         ),      # Byte 3
                           (  (0,1), (0,1), (0,1), (0,1), (0,1)                         )       # Byte 4
                          )

    ViirsAeroQualFillBoundaries = [
                                   [ [-0.5, 0.5, 1.5, 2.5, 3.5],[-0.5, 0.5, 1.5, 2.5, 3.5],[-0.5, 0.5, 1.5, 2.5, 3.5],[-0.5, 0.5, 1.5, 2.5, 3.5]   ], # Byte 0
                                   [ [-0.5, 0.5, 1.5, 2.5, 3.5],[-0.5, 0.5, 1.5, 2.5, 3.5], [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5], [-0.5, 0.5, 1.5] ], # Byte 1
                                   [ [-0.5, 0.5, 1.5, 2.5, 3.5],[-0.5, 0.5, 1.5, 2.5, 3.5], [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5]  ], # Byte 2
                                   [ [-0.5, 0.5, 1.5], [-0.5, 0.5, 1.5], [-0.5, 0.5, 1.5], [-0.5, 0.5, 1.5], [-0.5, 0.5, 1.5, 2.5], [-0.5, 0.5, 1.5], [-0.5, 0.5, 1.5] ], # Byte 3
                                   [ [-0.5, 0.5, 1.5],[-0.5, 0.5, 1.5],[-0.5, 0.5, 1.5],[-0.5, 0.5, 1.5],[-0.5, 0.5, 1.5]   ] # Byte 4
                                  ]

    ViirsAeroQualFillColours = [
                                   [ ['#00ff00','b','r','k'],['#00ff00','b','r','k'],['#00ff00','b','r','k'],['r','b','#00ff00','w'] ], # Byte 0
                                   [ ['#00ff00','#00ffff','#ff0000','w'],['#00ff00','#00ffff','#ff0000','w'],['brown','green','cyan','blue','w','yellow','r'], ['k','g'] ], # Byte 1
                                   [ ['yellow','orange','brown','k'],['w','#000080','#0012ff','#00a4ff'],['#000080','#0012ff','#00a4ff','#b6ff41','#ffb900','#ff3200','#800000','yellow'] ], # Byte 2
                                   [ ['k','g'], ['k','g'], ['k','g'], ['k','g'], ['brown','green','yellow'], ['k','g'], ['k','g'] ], # Byte 3
                                   [ ['k','g'], ['k','g'], ['k','g'], ['k','g']   ] # Byte 4
                               ]

    # The tick labels of the colourbar categories
    ViirsAeroQualTickNames = [
                              [ ['High Quality','Degraded Quality','Excluded Quality','Not Produced'],['High Quality','Degraded Quality','Excluded Quality','Not Produced'],['High Quality','Degraded Quality','Excluded Quality','Not Produced'],['Poor','Low','Medium','High'] ],      # Byte 0
                              [ ['Confident Clear','Probably Clear','Probably Cloudy','Confident Cloudy'],['Confident Clear','Probably Clear','Probably Cloudy','Confident Cloudy'],['Land & Desert','Land, No Desert','Inland Water','Sea Water','None','Coastal','Ephemeral Water'], ['No','Yes'] ],      # Byte 1
                              [ ['Day','Low Sun, Degraded','Twilight, Excluded','Night'],['None','Interpolation only','Interpolation & Climatology/NAAPS','Climatology/NAAPS'],['None','Geometry Based','Wind Speed Based','Geometry & Wind','Internal','Internal & Geometry','Internal & Wind','All'] ],      # Byte 2
                              [ ['No','Yes'], ['No','Yes'], ['No','Yes'], ['No','Yes'],[ 'Dark','Soil Dominated','Bright'], ['No','Yes'], ['No','Yes'] ],      # Byte 3
                              [ ['No','Yes'], ['No','Yes'], ['No','Yes'], ['No','Yes'] ]       # Byte 4
                             ]

    class AerosolProd:
        """
        This class provides a data structure to hold the plot attributes
        for a particular aerosol product.
        """
        def __init__(self,
            SDSname         = None,\
            SDSfactorsName  = None,\
            gridType        = None,\
            SDStitleString  = None,\
            vmin            = None,\
            vmax            = None,\
            cmap            = None,\
            cbarTickPos     = None,\
            cbarTickNames   = None,\
            cbarTitle       = None,\
            logScale        = False
            ):
            self.SDSname         = SDSname        
            self.SDSfactorsName  = SDSfactorsName        
            self.gridType        = gridType       
            self.SDStitleString  = SDStitleString 
            self.vmin            = vmin           
            self.vmax            = vmax           
            self.cmap            = cmap           
            self.cbarTickPos     = cbarTickPos    
            self.cbarTickNames   = cbarTickNames  
            self.cbarTitle       = cbarTitle      
            self.logScale        = logScale

    AerosolProduct = {
        'aot': AerosolProd(
            SDSname = '/All_Data/VIIRS-Aeros-Opt-Thick-IP_All/faot550',
            gridType = 'fine',
            SDStitleString = "VIIRS Aerosol\nOptical Depth",
            vmin = 0.,
            vmax = 1.,
            cbarTitle = "AOT",
            #cmap = cm.jet,
        ),
        'angexp': AerosolProd(
            SDSname = '/All_Data/VIIRS-Aeros-Opt-Thick-IP_All/angexp',
            gridType = 'fine',
            SDStitleString = "VIIRS Aerosol\nAngstrom Exponent",
            vmin = 0.,
            vmax = 3.,
            cbarTitle = "Angstrom Exponent",
            #cmap = cm.jet,
        ),
        'aot_edr': AerosolProd(
            SDSname = '/All_Data/VIIRS-Aeros-Opt-Thick-IP_All/faot550',
            gridType = 'fine',
            SDStitleString = "VIIRS Aerosol\nOptical Depth",
            vmin = 0.,
            vmax = 1.,
            cbarTitle = "AOT",
            #cmap = cm.jet,
        ),
    }


class SeaSurfaceTempProdData:
    """
    This class contains static data for the VIIRS SST Products.
    """

    ViirsSSTqualBitMasks = ((3, 12,  48, 64 , 128        ),
                            (1, 2, 12, 48, 64, 128       ),
                            (1, 2, 4, 8, 16, 32, 64, 128 ),
                            (1, 2, 4, 8, 16, 32, 64, 128 ))

    ViirsSSTqualBitShift = ((0, 2, 4, 6, 7         ),
                            (0, 1, 2, 4, 6, 7      ),
                            (0, 1, 2, 3, 4, 5, 6, 7),
                            (0, 1, 2, 3, 4, 5, 6, 7))

    ViirsSSTqualBitMaskNames = (('Skin SST Quality','Bulk SST Quality','SST State','Algorithm','Day / Night'),
                                 ('Bad LWIR Pixel','Bad SWIR Pixel','Cloud Confidence','Adjacent Pixel Cloud Confident Value','Thin Cirrus','Sea Ice'),
                                 ('Sun Glint','Exclusion, AOT > 1','Degraded, AOT > 0.6','Exclusion, Not Ocean','Degraded, HCS limit','Degraded, Sensor Zenith Angle > 40','Skin SST Outside Range','Bulk SST Outside Range'),
                                 ('Skin SST Degraded, T > 305 K','Bulk SST Degraded, T > 305 K','Spare','Spare','Spare','Spare','Spare','Spare'),
                                )
    
    ViirsSSTqualvalues = ( (  (0,1,2,3), (0,1,2,3), (0,1,2), (0,1), (0,1)             ),      # Byte 0
                           (  (0,1), (0,1), (0,1,2,3), (0,1,2,3), (0,1), (0,1)        ),      # Byte 1
                           (  (0,1), (0,1), (0,1), (0,1), (0,1), (0,1), (0,1), (0,1)  ),      # Byte 2
                           (  (0,1), (0,1), (0,1), (0,1), (0,1), (0,1), (0,1), (0,1)  )       # Byte 3
                          )

    ViirsSSTqualFillBoundaries = [
                                   [ [-0.5, 0.5, 1.5, 2.5, 3.5],[-0.5, 0.5, 1.5, 2.5, 3.5],[-0.5, 0.5, 1.5, 2.5],[-0.5, 0.5, 1.5],[-0.5, 0.5, 1.5]                          ], # Byte 0
                                   [ [-0.5, 0.5, 1.5],[-0.5, 0.5, 1.5],[-0.5, 0.5, 1.5, 2.5, 3.5],[-0.5, 0.5, 1.5, 2.5, 3.5], [-0.5, 0.5, 1.5],[-0.5, 0.5, 1.5]             ], # Byte 1
                                   [ [-0.5, 0.5, 1.5],[-0.5, 0.5, 1.5],[-0.5, 0.5, 1.5],[-0.5, 0.5, 1.5],[-0.5, 0.5, 1.5],[-0.5, 0.5, 1.5],[-0.5, 0.5, 1.5],[-0.5, 0.5, 1.5]], # Byte 3
                                   [ [-0.5, 0.5, 1.5],[-0.5, 0.5, 1.5],[-0.5, 0.5, 1.5],[-0.5, 0.5, 1.5],[-0.5, 0.5, 1.5],[-0.5, 0.5, 1.5],[-0.5, 0.5, 1.5],[-0.5, 0.5, 1.5]] # Byte 4
                                 ]

    ViirsSSTqualFillColours = [
                                [ ['k', 'r', 'b', '#00ff00'],['k', 'r', 'b', '#00ff00'],['brown','green','yellow'],['k','g'],['k','yellow']       ], # Byte 0
                                [ ['g','k'],['g','k'],['#00ff00','#00ffff','#ff0000','w'],['#00ff00','#00ffff','#ff0000','w'],['k','g'],['k','g'] ], # Byte 2
                                [ ['k','g'],['k','g'],['k','g'],['b','k'],['g','k'],['k','g'],['g','k'],['g','k']                                 ], # Byte 3
                                [ ['k','g'],['k','g'],['k','g'],['k','g'],['k','g'],['k','g'],['k','g'],['k','g']                                 ] # Byte 4
                              ]

    # The tick labels of the colourbar categories
    ViirsSSTqualTickNames = [
                              [ ['Not retrieved','Excluded','Degraded','High Quality'],['Not retrieved','Excluded','Degraded','High Quality'],[' Dry / None','Moist','Average'],['Non-linear Split Window','Triple Window'],['Night','Day'] ],      # Byte 0
                              [ ['Good SDR','Bad SDR'],['Good SDR','Bad SDR'],['Confident Clear','Probably Clear','Probably Cloudy','Confident Cloudy'],['Confident Clear','Probably Clear','Probably Cloudy','Confident Cloudy'],['No Thin Cirrus','Thin Cirrus'],['No Sea Ice','Sea Ice'] ],      # Byte 1
                              [ ['No sun glint','Sun glint'],['No','Yes'],['No','Yes'],['Ocean','Not ocean'],['Within HCS limit','Past HCS limit'],['No','Yes'],['In range','Out of range'],['In range','Out of range'] ],      # Byte 3
                              [ ['Not degraded','Degraded'],['Not degraded','Degraded'],['',''],['',''],['',''],['',''],['',''],['',''] ]       # Byte 4
                            ]
                              


    class SeaSurfaceTempProd:
        """
        This class provides a data structure to hold the plot attributes
        for a particular aerosol product.
        """
        def __init__(self,
            SDSname         = None,\
            SDSfactorsName  = None,\
            gridType        = None,\
            SDStitleString  = None,\
            vmin            = None,\
            vmax            = None,\
            cmap            = None,\
            cbarTickPos     = None,\
            cbarTickNames   = None,\
            cbarTitle       = None,\
            logScale        = False
            ):
            self.SDSname         = SDSname        
            self.SDSfactorsName  = SDSfactorsName        
            self.gridType        = gridType       
            self.SDStitleString  = SDStitleString 
            self.vmin            = vmin           
            self.vmax            = vmax           
            self.cmap            = cmap           
            self.cbarTickPos     = cbarTickPos    
            self.cbarTickNames   = cbarTickNames  
            self.cbarTitle       = cbarTitle      
            self.logScale        = logScale

    SeaSurfaceTempProduct = {
        'bulk_sst': SeaSurfaceTempProd(
            SDSname = '/All_Data/VIIRS-SST-EDR_All/BulkSST',
            SDStitleString = "VIIRS Bulk Sea Surface Temperature",
            vmin = 273.,
            vmax = 300.,
            cbarTitle = "SST ($K$)",
            #cmap = cm.jet,
        ),
        'skin_sst': SeaSurfaceTempProd(
            SDSname = '/All_Data/VIIRS-SST-EDR_All/SkinSST',
            SDStitleString = "VIIRS Skin Sea Surface Temperature",
            vmin = 273.,
            vmax = 300.,
            cbarTitle = "SST ($K$)",
            #cmap = cm.jet,
        ),
    }


def toJulianDate(inTime):
    '''
    Takes time in regular yyyymmdd format and returns time string in Julian yyyyddd format.
    '''
    inTime = str(inTime)
    try :
        return time.strftime("%Y%j",time.strptime(inTime, "%Y%m%d"))
    except :
        print "error: incorrect data format (%s). Should conform to yyyymmdd." % (inTime)
        return 1

def fromJulianDate(inTime):
    '''
    Takes time in Julian yyyyddd format and returns time string in regular yyyymmdd format .
    '''
    inTime = str(inTime)
    try :
        return time.strftime("%Y%m%d",time.strptime(inTime,"%Y%j"))
    except :
        print "error: incorrect data format (%s). Should conform to yyyyddd." % (inTime)
        return 1

