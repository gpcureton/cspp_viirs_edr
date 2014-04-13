#!/usr/bin/env python
# encoding: utf-8
"""
viirs_edr_data.py

Purpose: Provide required data for VIIRS Cloud and Aerosol products.

Created by Geoff Cureton <geoff.cureton@ssec.wisc.edu> on 2012-11-13.
Copyright (c) 2012-2013 University of Wisconsin Regents. All rights reserved.

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

__author__ = 'G.P. Cureton <geoff.cureton@ssec.wisc.edu>'
__version__ = '$Id$'
__docformat__ = 'Epytext'


from matplotlib.colors import Colormap, normalize, LinearSegmentedColormap,ListedColormap
from matplotlib import cm as cm
import numpy as np
import time


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
                           #(  (0,1,2,3), (0,1,2,3), (0,1,2,3,4,5,6), (0,1)              ),      # Byte 1
                           (  (0,1,2,3), (0,1,2,3), (0,1,2,3,4,5), (0,1)              ),      # Byte 1
                           (  (0,1,2,3), (0,1,2,3), (0,1,2,3,4,5,6,7)                   ),      # Byte 2
                           (  (0,1), (0,1), (0,1), (0,1), (0,1,2), (0,1), (0,1)         ),      # Byte 3
                           (  (0,1), (0,1), (0,1), (0,1), (0,1)                         )       # Byte 4
                          )

    ViirsAeroQualFillBoundaries = [
                                   [ [-0.5, 0.5, 1.5, 2.5, 3.5],[-0.5, 0.5, 1.5, 2.5, 3.5],[-0.5, 0.5, 1.5, 2.5, 3.5],[-0.5, 0.5, 1.5, 2.5, 3.5]   ], # Byte 0
                                   #[ [-0.5, 0.5, 1.5, 2.5, 3.5],[-0.5, 0.5, 1.5, 2.5, 3.5], [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5], [-0.5, 0.5, 1.5] ], # Byte 1
                                   [ [-0.5, 0.5, 1.5, 2.5, 3.5],[-0.5, 0.5, 1.5, 2.5, 3.5], [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5], [-0.5, 0.5, 1.5] ], # Byte 1
                                   [ [-0.5, 0.5, 1.5, 2.5, 3.5],[-0.5, 0.5, 1.5, 2.5, 3.5], [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5]  ], # Byte 2
                                   [ [-0.5, 0.5, 1.5], [-0.5, 0.5, 1.5], [-0.5, 0.5, 1.5], [-0.5, 0.5, 1.5], [-0.5, 0.5, 1.5, 2.5], [-0.5, 0.5, 1.5], [-0.5, 0.5, 1.5] ], # Byte 3
                                   [ [-0.5, 0.5, 1.5],[-0.5, 0.5, 1.5],[-0.5, 0.5, 1.5],[-0.5, 0.5, 1.5],[-0.5, 0.5, 1.5]   ] # Byte 4
                                  ]

    ViirsAeroQualFillColours = [
                                   [ ['#00ff00','b','r','k'],['#00ff00','b','r','k'],['#00ff00','b','r','k'],['r','b','#00ff00','w'] ], # Byte 0
                                   #[ ['#00ff00','#00ffff','#ff0000','w'],['#00ff00','#00ffff','#ff0000','w'],['brown','green','cyan','blue','w','yellow','r'], ['k','g'] ], # Byte 1
                                   [ ['#00ff00','#00ffff','#ff0000','w'],['#00ff00','#00ffff','#ff0000','w'],['brown','green','cyan','blue','w','yellow'], ['k','g'] ], # Byte 1
                                   [ ['yellow','orange','brown','k'],['w','#000080','#0012ff','#00a4ff'],['#000080','#0012ff','#00a4ff','#b6ff41','#ffb900','#ff3200','#800000','yellow'] ], # Byte 2
                                   [ ['k','g'], ['k','g'], ['k','g'], ['k','g'], ['brown','green','yellow'], ['k','g'], ['k','g'] ], # Byte 3
                                   [ ['k','g'], ['k','g'], ['k','g'], ['k','g']   ] # Byte 4
                               ]

    # The tick labels of the colourbar categories
    ViirsAeroQualTickNames = [
                              [ ['High Quality','Degraded Quality','Excluded Quality','Not Produced'],['High Quality','Degraded Quality','Excluded Quality','Not Produced'],['High Quality','Degraded Quality','Excluded Quality','Not Produced'],['Poor','Low','Medium','High'] ],      # Byte 0
                              #[ ['Confident Clear','Probably Clear','Probably Cloudy','Confident Cloudy'],['Confident Clear','Probably Clear','Probably Cloudy','Confident Cloudy'],['Land/Desert','Land','Inland Water','Sea Water','None','Coastal','Ephemeral Water'], ['No','Yes'] ],      # Byte 1
                              [ ['Confident Clear','Probably Clear','Probably Cloudy','Confident Cloudy'],['Confident Clear','Probably Clear','Probably Cloudy','Confident Cloudy'],['Land/Desert','Land','Inland Water','Sea Water','None','Coastal'], ['No','Yes'] ],      # Byte 1
                              [ ['Day','Low Sun, Degraded','Twilight, Excluded','Night'],['None','Interpolation only','Interpolation \n& Climatology/NAAPS','Climatology/NAAPS'],['None','Geometry \nBased','Wind Speed \nBased','Geometry \n& Wind','Internal','Internal \n& Geometry','Internal \n& Wind','All'] ],      # Byte 2
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

            sst_cmdata = {
                'blue': [ (0.0, 0.46274509803921571, 0.46274509803921571),
                          (0.0039215686274509803, 0.54509803921568623, 0.54509803921568623),
                          (0.0078431372549019607, 0.6705882352941176, 0.6705882352941176),
                          (0.011764705882352941, 0.792156862745098, 0.792156862745098),
                          (0.015686274509803921, 0.87450980392156863, 0.87450980392156863),
                          (0.019607843137254902, 0.91764705882352937, 0.91764705882352937),
                          (0.023529411764705882, 0.93333333333333335, 0.93333333333333335),
                          (0.027450980392156862, 0.93333333333333335, 0.93333333333333335),
                          (0.031372549019607843, 0.93725490196078431, 0.93725490196078431),
                          (0.035294117647058823, 0.93725490196078431, 0.93725490196078431),
                          (0.039215686274509803, 0.93725490196078431, 0.93725490196078431),
                          (0.043137254901960784, 0.93725490196078431, 0.93725490196078431),
                          (0.047058823529411764, 0.93333333333333335, 0.93333333333333335),
                          (0.050980392156862744, 0.92549019607843142, 0.92549019607843142),
                          (0.054901960784313725, 0.9137254901960784, 0.9137254901960784),
                          (0.058823529411764705, 0.89411764705882357, 0.89411764705882357),
                          (0.062745098039215685, 0.8784313725490196, 0.8784313725490196),
                          (0.066666666666666666, 0.8666666666666667, 0.8666666666666667),
                          (0.070588235294117646, 0.85098039215686272, 0.85098039215686272),
                          (0.074509803921568626, 0.83137254901960789, 0.83137254901960789),
                          (0.078431372549019607, 0.81176470588235294, 0.81176470588235294),
                          (0.082352941176470587, 0.78823529411764703, 0.78823529411764703),
                          (0.086274509803921567, 0.7686274509803922, 0.7686274509803922),
                          (0.090196078431372548, 0.74901960784313726, 0.74901960784313726),
                          (0.094117647058823528, 0.72549019607843135, 0.72549019607843135),
                          (0.098039215686274508, 0.70196078431372544, 0.70196078431372544),
                          (0.10196078431372549, 0.67843137254901964, 0.67843137254901964),
                          (0.10588235294117647, 0.65098039215686276, 0.65098039215686276),
                          (0.10980392156862745, 0.62352941176470589, 0.62352941176470589),
                          (0.11372549019607843, 0.59999999999999998, 0.59999999999999998),
                          (0.11764705882352941, 0.57647058823529407, 0.57647058823529407),
                          (0.12156862745098039, 0.54509803921568623, 0.54509803921568623),
                          (0.12549019607843137, 0.50980392156862742, 0.50980392156862742),
                          (0.12941176470588237, 0.47450980392156861, 0.47450980392156861),
                          (0.13333333333333333, 0.44705882352941179, 0.44705882352941179),
                          (0.13725490196078433, 0.42745098039215684, 0.42745098039215684),
                          (0.14117647058823529, 0.41568627450980394, 0.41568627450980394),
                          (0.14509803921568629, 0.40784313725490196, 0.40784313725490196),
                          (0.14901960784313725, 0.40784313725490196, 0.40784313725490196),
                          (0.15294117647058825, 0.40784313725490196, 0.40784313725490196),
                          (0.15686274509803921, 0.41568627450980394, 0.41568627450980394),
                          (0.16078431372549021, 0.41960784313725491, 0.41960784313725491),
                          (0.16470588235294117, 0.42745098039215684, 0.42745098039215684),
                          (0.16862745098039217, 0.4392156862745098, 0.4392156862745098),
                          (0.17254901960784313, 0.44705882352941179, 0.44705882352941179),
                          (0.17647058823529413, 0.45882352941176469, 0.45882352941176469),
                          (0.1803921568627451, 0.47843137254901963, 0.47843137254901963),
                          (0.18431372549019609, 0.50980392156862742, 0.50980392156862742),
                          (0.18823529411764706, 0.54509803921568623, 0.54509803921568623),
                          (0.19215686274509805, 0.57647058823529407, 0.57647058823529407),
                          (0.19607843137254902, 0.59607843137254901, 0.59607843137254901),
                          (0.20000000000000001, 0.60784313725490191, 0.60784313725490191),
                          (0.20392156862745098, 0.61568627450980395, 0.61568627450980395),
                          (0.20784313725490197, 0.61960784313725492, 0.61960784313725492),
                          (0.21176470588235294, 0.62352941176470589, 0.62352941176470589),
                          (0.21568627450980393, 0.62352941176470589, 0.62352941176470589),
                          (0.2196078431372549, 0.62745098039215685, 0.62745098039215685),
                          (0.22352941176470589, 0.63529411764705879, 0.63529411764705879),
                          (0.22745098039215686, 0.63921568627450975, 0.63921568627450975),
                          (0.23137254901960785, 0.6470588235294118, 0.6470588235294118),
                          (0.23529411764705882, 0.65098039215686276, 0.65098039215686276),
                          (0.23921568627450981, 0.65490196078431373, 0.65490196078431373),
                          (0.24313725490196078, 0.66274509803921566, 0.66274509803921566),
                          (0.24705882352941178, 0.67450980392156867, 0.67450980392156867),
                          (0.25098039215686274, 0.69411764705882351, 0.69411764705882351),
                          (0.25490196078431371, 0.71372549019607845, 0.71372549019607845),
                          (0.25882352941176473, 0.73333333333333328, 0.73333333333333328),
                          (0.2627450980392157, 0.75294117647058822, 0.75294117647058822),
                          (0.26666666666666666, 0.77254901960784317, 0.77254901960784317),
                          (0.27058823529411763, 0.78431372549019607, 0.78431372549019607),
                          (0.27450980392156865, 0.80000000000000004, 0.80000000000000004),
                          (0.27843137254901962, 0.81960784313725488, 0.81960784313725488),
                          (0.28235294117647058, 0.83529411764705885, 0.83529411764705885),
                          (0.28627450980392155, 0.84705882352941175, 0.84705882352941175),
                          (0.29019607843137257, 0.85882352941176465, 0.85882352941176465),
                          (0.29411764705882354, 0.87058823529411766, 0.87058823529411766),
                          (0.29803921568627451, 0.88235294117647056, 0.88235294117647056),
                          (0.30196078431372547, 0.8901960784313725, 0.8901960784313725),
                          (0.30588235294117649, 0.90196078431372551, 0.90196078431372551),
                          (0.30980392156862746, 0.9137254901960784, 0.9137254901960784),
                          (0.31372549019607843, 0.92156862745098034, 0.92156862745098034),
                          (0.31764705882352939, 0.92549019607843142, 0.92549019607843142),
                          (0.32156862745098042, 0.91764705882352937, 0.91764705882352937),
                          (0.32549019607843138, 0.89803921568627454, 0.89803921568627454),
                          (0.32941176470588235, 0.87450980392156863, 0.87450980392156863),
                          (0.33333333333333331, 0.84705882352941175, 0.84705882352941175),
                          (0.33725490196078434, 0.81960784313725488, 0.81960784313725488),
                          (0.3411764705882353, 0.792156862745098, 0.792156862745098),
                          (0.34509803921568627, 0.75294117647058822, 0.75294117647058822),
                          (0.34901960784313724, 0.69803921568627447, 0.69803921568627447),
                          (0.35294117647058826, 0.63529411764705879, 0.63529411764705879),
                          (0.35686274509803922, 0.57647058823529407, 0.57647058823529407),
                          (0.36078431372549019, 0.52156862745098043, 0.52156862745098043),
                          (0.36470588235294116, 0.4823529411764706, 0.4823529411764706),
                          (0.36862745098039218, 0.45490196078431372, 0.45490196078431372),
                          (0.37254901960784315, 0.44705882352941179, 0.44705882352941179),
                          (0.37647058823529411, 0.44705882352941179, 0.44705882352941179),
                          (0.38039215686274508, 0.45098039215686275, 0.45098039215686275),
                          (0.3843137254901961, 0.45882352941176469, 0.45882352941176469),
                          (0.38823529411764707, 0.46274509803921571, 0.46274509803921571),
                          (0.39215686274509803, 0.46274509803921571, 0.46274509803921571),
                          (0.396078431372549, 0.45882352941176469, 0.45882352941176469),
                          (0.40000000000000002, 0.45098039215686275, 0.45098039215686275),
                          (0.40392156862745099, 0.4392156862745098, 0.4392156862745098),
                          (0.40784313725490196, 0.41960784313725491, 0.41960784313725491),
                          (0.41176470588235292, 0.38823529411764707, 0.38823529411764707),
                          (0.41568627450980394, 0.35686274509803922, 0.35686274509803922),
                          (0.41960784313725491, 0.32941176470588235, 0.32941176470588235),
                          (0.42352941176470588, 0.31764705882352939, 0.31764705882352939),
                          (0.42745098039215684, 0.31372549019607843, 0.31372549019607843),
                          (0.43137254901960786, 0.31372549019607843, 0.31372549019607843),
                          (0.43529411764705883, 0.31372549019607843, 0.31372549019607843),
                          (0.4392156862745098, 0.31372549019607843, 0.31372549019607843),
                          (0.44313725490196076, 0.30980392156862746, 0.30980392156862746),
                          (0.44705882352941179, 0.30980392156862746, 0.30980392156862746),
                          (0.45098039215686275, 0.30980392156862746, 0.30980392156862746),
                          (0.45490196078431372, 0.30588235294117649, 0.30588235294117649),
                          (0.45882352941176469, 0.29411764705882354, 0.29411764705882354),
                          (0.46274509803921571, 0.27450980392156865, 0.27450980392156865),
                          (0.46666666666666667, 0.24705882352941178, 0.24705882352941178),
                          (0.47058823529411764, 0.2196078431372549, 0.2196078431372549),
                          (0.47450980392156861, 0.18823529411764706, 0.18823529411764706),
                          (0.47843137254901963, 0.14509803921568629, 0.14509803921568629),
                          (0.4823529411764706, 0.094117647058823528, 0.094117647058823528),
                          (0.48627450980392156, 0.054901960784313725, 0.054901960784313725),
                          (0.49019607843137253, 0.035294117647058823, 0.035294117647058823),
                          (0.49411764705882355, 0.031372549019607843, 0.031372549019607843),
                          (0.49803921568627452, 0.031372549019607843, 0.031372549019607843),
                          (0.50196078431372548, 0.031372549019607843, 0.031372549019607843),
                          (0.50588235294117645, 0.031372549019607843, 0.031372549019607843),
                          (0.50980392156862742, 0.031372549019607843, 0.031372549019607843),
                          (0.51372549019607838, 0.031372549019607843, 0.031372549019607843),
                          (0.51764705882352946, 0.031372549019607843, 0.031372549019607843),
                          (0.52156862745098043, 0.031372549019607843, 0.031372549019607843),
                          (0.52549019607843139, 0.031372549019607843, 0.031372549019607843),
                          (0.52941176470588236, 0.031372549019607843, 0.031372549019607843),
                          (0.53333333333333333, 0.031372549019607843, 0.031372549019607843),
                          (0.53725490196078429, 0.031372549019607843, 0.031372549019607843),
                          (0.54117647058823526, 0.031372549019607843, 0.031372549019607843),
                          (0.54509803921568623, 0.031372549019607843, 0.031372549019607843),
                          (0.5490196078431373, 0.031372549019607843, 0.031372549019607843),
                          (0.55294117647058827, 0.031372549019607843, 0.031372549019607843),
                          (0.55686274509803924, 0.031372549019607843, 0.031372549019607843),
                          (0.5607843137254902, 0.031372549019607843, 0.031372549019607843),
                          (0.56470588235294117, 0.031372549019607843, 0.031372549019607843),
                          (0.56862745098039214, 0.031372549019607843, 0.031372549019607843),
                          (0.5725490196078431, 0.031372549019607843, 0.031372549019607843),
                          (0.57647058823529407, 0.031372549019607843, 0.031372549019607843),
                          (0.58039215686274515, 0.031372549019607843, 0.031372549019607843),
                          (0.58431372549019611, 0.031372549019607843, 0.031372549019607843),
                          (0.58823529411764708, 0.031372549019607843, 0.031372549019607843),
                          (0.59215686274509804, 0.031372549019607843, 0.031372549019607843),
                          (0.59607843137254901, 0.031372549019607843, 0.031372549019607843),
                          (0.59999999999999998, 0.031372549019607843, 0.031372549019607843),
                          (0.60392156862745094, 0.031372549019607843, 0.031372549019607843),
                          (0.60784313725490191, 0.031372549019607843, 0.031372549019607843),
                          (0.61176470588235299, 0.031372549019607843, 0.031372549019607843),
                          (0.61568627450980395, 0.031372549019607843, 0.031372549019607843),
                          (0.61960784313725492, 0.031372549019607843, 0.031372549019607843),
                          (0.62352941176470589, 0.031372549019607843, 0.031372549019607843),
                          (0.62745098039215685, 0.031372549019607843, 0.031372549019607843),
                          (0.63137254901960782, 0.031372549019607843, 0.031372549019607843),
                          (0.63529411764705879, 0.031372549019607843, 0.031372549019607843),
                          (0.63921568627450975, 0.031372549019607843, 0.031372549019607843),
                          (0.64313725490196083, 0.031372549019607843, 0.031372549019607843),
                          (0.6470588235294118, 0.031372549019607843, 0.031372549019607843),
                          (0.65098039215686276, 0.031372549019607843, 0.031372549019607843),
                          (0.65490196078431373, 0.031372549019607843, 0.031372549019607843),
                          (0.6588235294117647, 0.031372549019607843, 0.031372549019607843),
                          (0.66274509803921566, 0.031372549019607843, 0.031372549019607843),
                          (0.66666666666666663, 0.031372549019607843, 0.031372549019607843),
                          (0.6705882352941176, 0.031372549019607843, 0.031372549019607843),
                          (0.67450980392156867, 0.031372549019607843, 0.031372549019607843),
                          (0.67843137254901964, 0.031372549019607843, 0.031372549019607843),
                          (0.68235294117647061, 0.031372549019607843, 0.031372549019607843),
                          (0.68627450980392157, 0.031372549019607843, 0.031372549019607843),
                          (0.69019607843137254, 0.031372549019607843, 0.031372549019607843),
                          (0.69411764705882351, 0.031372549019607843, 0.031372549019607843),
                          (0.69803921568627447, 0.031372549019607843, 0.031372549019607843),
                          (0.70196078431372544, 0.031372549019607843, 0.031372549019607843),
                          (0.70588235294117652, 0.031372549019607843, 0.031372549019607843),
                          (0.70980392156862748, 0.031372549019607843, 0.031372549019607843),
                          (0.71372549019607845, 0.031372549019607843, 0.031372549019607843),
                          (0.71764705882352942, 0.031372549019607843, 0.031372549019607843),
                          (0.72156862745098038, 0.031372549019607843, 0.031372549019607843),
                          (0.72549019607843135, 0.031372549019607843, 0.031372549019607843),
                          (0.72941176470588232, 0.031372549019607843, 0.031372549019607843),
                          (0.73333333333333328, 0.031372549019607843, 0.031372549019607843),
                          (0.73725490196078436, 0.031372549019607843, 0.031372549019607843),
                          (0.74117647058823533, 0.031372549019607843, 0.031372549019607843),
                          (0.74509803921568629, 0.031372549019607843, 0.031372549019607843),
                          (0.74901960784313726, 0.031372549019607843, 0.031372549019607843),
                          (0.75294117647058822, 0.031372549019607843, 0.031372549019607843),
                          (0.75686274509803919, 0.031372549019607843, 0.031372549019607843),
                          (0.76078431372549016, 0.031372549019607843, 0.031372549019607843),
                          (0.76470588235294112, 0.031372549019607843, 0.031372549019607843),
                          (0.7686274509803922, 0.031372549019607843, 0.031372549019607843),
                          (0.77254901960784317, 0.031372549019607843, 0.031372549019607843),
                          (0.77647058823529413, 0.031372549019607843, 0.031372549019607843),
                          (0.7803921568627451, 0.031372549019607843, 0.031372549019607843),
                          (0.78431372549019607, 0.031372549019607843, 0.031372549019607843),
                          (0.78823529411764703, 0.031372549019607843, 0.031372549019607843),
                          (0.792156862745098, 0.031372549019607843, 0.031372549019607843),
                          (0.79607843137254897, 0.050980392156862744, 0.050980392156862744),
                          (0.80000000000000004, 0.074509803921568626, 0.074509803921568626),
                          (0.80392156862745101, 0.086274509803921567, 0.086274509803921567),
                          (0.80784313725490198, 0.10196078431372549, 0.10196078431372549),
                          (0.81176470588235294, 0.11764705882352941, 0.11764705882352941),
                          (0.81568627450980391, 0.12941176470588237, 0.12941176470588237),
                          (0.81960784313725488, 0.14901960784313725, 0.14901960784313725),
                          (0.82352941176470584, 0.16078431372549021, 0.16078431372549021),
                          (0.82745098039215681, 0.17647058823529413, 0.17647058823529413),
                          (0.83137254901960789, 0.19215686274509805, 0.19215686274509805),
                          (0.83529411764705885, 0.20784313725490197, 0.20784313725490197),
                          (0.83921568627450982, 0.2196078431372549, 0.2196078431372549),
                          (0.84313725490196079, 0.23529411764705882, 0.23529411764705882),
                          (0.84705882352941175, 0.24705882352941178, 0.24705882352941178),
                          (0.85098039215686272, 0.2627450980392157, 0.2627450980392157),
                          (0.85490196078431369, 0.27450980392156865, 0.27450980392156865),
                          (0.85882352941176465, 0.29019607843137257, 0.29019607843137257),
                          (0.86274509803921573, 0.30588235294117649, 0.30588235294117649),
                          (0.8666666666666667, 0.32156862745098042, 0.32156862745098042),
                          (0.87058823529411766, 0.33333333333333331, 0.33333333333333331),
                          (0.87450980392156863, 0.35294117647058826, 0.35294117647058826),
                          (0.8784313725490196, 0.36470588235294116, 0.36470588235294116),
                          (0.88235294117647056, 0.38039215686274508, 0.38039215686274508),
                          (0.88627450980392153, 0.39215686274509803, 0.39215686274509803),
                          (0.8901960784313725, 0.40784313725490196, 0.40784313725490196),
                          (0.89411764705882357, 0.42352941176470588, 0.42352941176470588),
                          (0.89803921568627454, 0.43529411764705883, 0.43529411764705883),
                          (0.90196078431372551, 0.45098039215686275, 0.45098039215686275),
                          (0.90588235294117647, 0.46666666666666667, 0.46666666666666667),
                          (0.90980392156862744, 0.4823529411764706, 0.4823529411764706),
                          (0.9137254901960784, 0.49411764705882355, 0.49411764705882355),
                          (0.91764705882352937, 0.50980392156862742, 0.50980392156862742),
                          (0.92156862745098034, 0.52549019607843139, 0.52549019607843139),
                          (0.92549019607843142, 0.53725490196078429, 0.53725490196078429),
                          (0.92941176470588238, 0.55294117647058827, 0.55294117647058827),
                          (0.93333333333333335, 0.56862745098039214, 0.56862745098039214),
                          (0.93725490196078431, 0.58431372549019611, 0.58431372549019611),
                          (0.94117647058823528, 0.59999999999999998, 0.59999999999999998),
                          (0.94509803921568625, 0.61176470588235299, 0.61176470588235299),
                          (0.94901960784313721, 0.62352941176470589, 0.62352941176470589),
                          (0.95294117647058818, 0.63921568627450975, 0.63921568627450975),
                          (0.95686274509803926, 0.65490196078431373, 0.65490196078431373),
                          (0.96078431372549022, 0.6705882352941176, 0.6705882352941176),
                          (0.96470588235294119, 0.68235294117647061, 0.68235294117647061),
                          (0.96862745098039216, 0.69803921568627447, 0.69803921568627447),
                          (0.97254901960784312, 0.71372549019607845, 0.71372549019607845),
                          (0.97647058823529409, 0.72549019607843135, 0.72549019607843135),
                          (0.98039215686274506, 0.74117647058823533, 0.74117647058823533),
                          (0.98431372549019602, 0.75686274509803919, 0.75686274509803919),
                          (0.9882352941176471, 0.77254901960784317, 0.77254901960784317),
                          (0.99215686274509807, 0.78823529411764703, 0.78823529411764703),
                          (0.99607843137254903, 0.80000000000000004, 0.80000000000000004),
                          (1.0, 0.81176470588235294, 0.81176470588235294)],

               'green': [ (0.0, 0.039215686274509803, 0.039215686274509803),
                          (0.0039215686274509803, 0.035294117647058823, 0.035294117647058823),
                          (0.0078431372549019607, 0.027450980392156862, 0.027450980392156862),
                          (0.011764705882352941, 0.027450980392156862, 0.027450980392156862),
                          (0.015686274509803921, 0.027450980392156862, 0.027450980392156862),
                          (0.019607843137254902, 0.027450980392156862, 0.027450980392156862),
                          (0.023529411764705882, 0.027450980392156862, 0.027450980392156862),
                          (0.027450980392156862, 0.027450980392156862, 0.027450980392156862),
                          (0.031372549019607843, 0.027450980392156862, 0.027450980392156862),
                          (0.035294117647058823, 0.027450980392156862, 0.027450980392156862),
                          (0.039215686274509803, 0.027450980392156862, 0.027450980392156862),
                          (0.043137254901960784, 0.027450980392156862, 0.027450980392156862),
                          (0.047058823529411764, 0.031372549019607843, 0.031372549019607843),
                          (0.050980392156862744, 0.031372549019607843, 0.031372549019607843),
                          (0.054901960784313725, 0.031372549019607843, 0.031372549019607843),
                          (0.058823529411764705, 0.031372549019607843, 0.031372549019607843),
                          (0.062745098039215685, 0.031372549019607843, 0.031372549019607843),
                          (0.066666666666666666, 0.031372549019607843, 0.031372549019607843),
                          (0.070588235294117646, 0.031372549019607843, 0.031372549019607843),
                          (0.074509803921568626, 0.031372549019607843, 0.031372549019607843),
                          (0.078431372549019607, 0.031372549019607843, 0.031372549019607843),
                          (0.082352941176470587, 0.031372549019607843, 0.031372549019607843),
                          (0.086274509803921567, 0.031372549019607843, 0.031372549019607843),
                          (0.090196078431372548, 0.031372549019607843, 0.031372549019607843),
                          (0.094117647058823528, 0.031372549019607843, 0.031372549019607843),
                          (0.098039215686274508, 0.031372549019607843, 0.031372549019607843),
                          (0.10196078431372549, 0.031372549019607843, 0.031372549019607843),
                          (0.10588235294117647, 0.031372549019607843, 0.031372549019607843),
                          (0.10980392156862745, 0.031372549019607843, 0.031372549019607843),
                          (0.11372549019607843, 0.031372549019607843, 0.031372549019607843),
                          (0.11764705882352941, 0.043137254901960784, 0.043137254901960784),
                          (0.12156862745098039, 0.062745098039215685, 0.062745098039215685),
                          (0.12549019607843137, 0.090196078431372548, 0.090196078431372548),
                          (0.12941176470588237, 0.11764705882352941, 0.11764705882352941),
                          (0.13333333333333333, 0.14901960784313725, 0.14901960784313725),
                          (0.13725490196078433, 0.1803921568627451, 0.1803921568627451),
                          (0.14117647058823529, 0.21176470588235294, 0.21176470588235294),
                          (0.14509803921568629, 0.24313725490196078, 0.24313725490196078),
                          (0.14901960784313725, 0.27450980392156865, 0.27450980392156865),
                          (0.15294117647058825, 0.30588235294117649, 0.30588235294117649),
                          (0.15686274509803921, 0.33725490196078434, 0.33725490196078434),
                          (0.16078431372549021, 0.36470588235294116, 0.36470588235294116),
                          (0.16470588235294117, 0.3843137254901961, 0.3843137254901961),
                          (0.16862745098039217, 0.39215686274509803, 0.39215686274509803),
                          (0.17254901960784313, 0.38823529411764707, 0.38823529411764707),
                          (0.17647058823529413, 0.3843137254901961, 0.3843137254901961),
                          (0.1803921568627451, 0.38039215686274508, 0.38039215686274508),
                          (0.18431372549019609, 0.3843137254901961, 0.3843137254901961),
                          (0.18823529411764706, 0.396078431372549, 0.396078431372549),
                          (0.19215686274509805, 0.40784313725490196, 0.40784313725490196),
                          (0.19607843137254902, 0.42352941176470588, 0.42352941176470588),
                          (0.20000000000000001, 0.44313725490196076, 0.44313725490196076),
                          (0.20392156862745098, 0.46274509803921571, 0.46274509803921571),
                          (0.20784313725490197, 0.4823529411764706, 0.4823529411764706),
                          (0.21176470588235294, 0.50196078431372548, 0.50196078431372548),
                          (0.21568627450980393, 0.52156862745098043, 0.52156862745098043),
                          (0.2196078431372549, 0.53333333333333333, 0.53333333333333333),
                          (0.22352941176470589, 0.5490196078431373, 0.5490196078431373),
                          (0.22745098039215686, 0.56862745098039214, 0.56862745098039214),
                          (0.23137254901960785, 0.58823529411764708, 0.58823529411764708),
                          (0.23529411764705882, 0.61176470588235299, 0.61176470588235299),
                          (0.23921568627450981, 0.63137254901960782, 0.63137254901960782),
                          (0.24313725490196078, 0.65098039215686276, 0.65098039215686276),
                          (0.24705882352941178, 0.67450980392156867, 0.67450980392156867),
                          (0.25098039215686274, 0.69411764705882351, 0.69411764705882351),
                          (0.25490196078431371, 0.71372549019607845, 0.71372549019607845),
                          (0.25882352941176473, 0.73333333333333328, 0.73333333333333328),
                          (0.2627450980392157, 0.75294117647058822, 0.75294117647058822),
                          (0.26666666666666666, 0.77254901960784317, 0.77254901960784317),
                          (0.27058823529411763, 0.78431372549019607, 0.78431372549019607),
                          (0.27450980392156865, 0.80000000000000004, 0.80000000000000004),
                          (0.27843137254901962, 0.81960784313725488, 0.81960784313725488),
                          (0.28235294117647058, 0.83529411764705885, 0.83529411764705885),
                          (0.28627450980392155, 0.84705882352941175, 0.84705882352941175),
                          (0.29019607843137257, 0.85882352941176465, 0.85882352941176465),
                          (0.29411764705882354, 0.87058823529411766, 0.87058823529411766),
                          (0.29803921568627451, 0.88235294117647056, 0.88235294117647056),
                          (0.30196078431372547, 0.8901960784313725, 0.8901960784313725),
                          (0.30588235294117649, 0.90196078431372551, 0.90196078431372551),
                          (0.30980392156862746, 0.9137254901960784, 0.9137254901960784),
                          (0.31372549019607843, 0.92156862745098034, 0.92156862745098034),
                          (0.31764705882352939, 0.92941176470588238, 0.92941176470588238),
                          (0.32156862745098042, 0.93333333333333335, 0.93333333333333335),
                          (0.32549019607843138, 0.93333333333333335, 0.93333333333333335),
                          (0.32941176470588235, 0.93333333333333335, 0.93333333333333335),
                          (0.33333333333333331, 0.92941176470588238, 0.92941176470588238),
                          (0.33725490196078434, 0.92156862745098034, 0.92156862745098034),
                          (0.3411764705882353, 0.91764705882352937, 0.91764705882352937),
                          (0.34509803921568627, 0.90980392156862744, 0.90980392156862744),
                          (0.34901960784313724, 0.90588235294117647, 0.90588235294117647),
                          (0.35294117647058826, 0.90196078431372551, 0.90196078431372551),
                          (0.35686274509803922, 0.89803921568627454, 0.89803921568627454),
                          (0.36078431372549019, 0.8901960784313725, 0.8901960784313725),
                          (0.36470588235294116, 0.88627450980392153, 0.88627450980392153),
                          (0.36862745098039218, 0.8784313725490196, 0.8784313725490196),
                          (0.37254901960784315, 0.8666666666666667, 0.8666666666666667),
                          (0.37647058823529411, 0.85882352941176465, 0.85882352941176465),
                          (0.38039215686274508, 0.84705882352941175, 0.84705882352941175),
                          (0.3843137254901961, 0.83529411764705885, 0.83529411764705885),
                          (0.38823529411764707, 0.81960784313725488, 0.81960784313725488),
                          (0.39215686274509803, 0.80000000000000004, 0.80000000000000004),
                          (0.396078431372549, 0.7803921568627451, 0.7803921568627451),
                          (0.40000000000000002, 0.75686274509803919, 0.75686274509803919),
                          (0.40392156862745099, 0.73725490196078436, 0.73725490196078436),
                          (0.40784313725490196, 0.72156862745098038, 0.72156862745098038),
                          (0.41176470588235292, 0.70980392156862748, 0.70980392156862748),
                          (0.41568627450980394, 0.69803921568627447, 0.69803921568627447),
                          (0.41960784313725491, 0.68627450980392157, 0.68627450980392157),
                          (0.42352941176470588, 0.67450980392156867, 0.67450980392156867),
                          (0.42745098039215684, 0.66274509803921566, 0.66274509803921566),
                          (0.43137254901960786, 0.64313725490196083, 0.64313725490196083),
                          (0.43529411764705883, 0.62352941176470589, 0.62352941176470589),
                          (0.4392156862745098, 0.59999999999999998, 0.59999999999999998),
                          (0.44313725490196076, 0.58039215686274515, 0.58039215686274515),
                          (0.44705882352941179, 0.5607843137254902, 0.5607843137254902),
                          (0.45098039215686275, 0.54117647058823526, 0.54117647058823526),
                          (0.45490196078431372, 0.52549019607843139, 0.52549019607843139),
                          (0.45882352941176469, 0.51764705882352946, 0.51764705882352946),
                          (0.46274509803921571, 0.52549019607843139, 0.52549019607843139),
                          (0.46666666666666667, 0.53725490196078429, 0.53725490196078429),
                          (0.47058823529411764, 0.5490196078431373, 0.5490196078431373),
                          (0.47450980392156861, 0.5490196078431373, 0.5490196078431373),
                          (0.47843137254901963, 0.54509803921568623, 0.54509803921568623),
                          (0.4823529411764706, 0.54117647058823526, 0.54117647058823526),
                          (0.48627450980392156, 0.54509803921568623, 0.54509803921568623),
                          (0.49019607843137253, 0.55294117647058827, 0.55294117647058827),
                          (0.49411764705882355, 0.5725490196078431, 0.5725490196078431),
                          (0.49803921568627452, 0.59215686274509804, 0.59215686274509804),
                          (0.50196078431372548, 0.61176470588235299, 0.61176470588235299),
                          (0.50588235294117645, 0.63137254901960782, 0.63137254901960782),
                          (0.50980392156862742, 0.65098039215686276, 0.65098039215686276),
                          (0.51372549019607838, 0.67450980392156867, 0.67450980392156867),
                          (0.51764705882352946, 0.69411764705882351, 0.69411764705882351),
                          (0.52156862745098043, 0.71372549019607845, 0.71372549019607845),
                          (0.52549019607843139, 0.73333333333333328, 0.73333333333333328),
                          (0.52941176470588236, 0.75294117647058822, 0.75294117647058822),
                          (0.53333333333333333, 0.7686274509803922, 0.7686274509803922),
                          (0.53725490196078429, 0.7803921568627451, 0.7803921568627451),
                          (0.54117647058823526, 0.78823529411764703, 0.78823529411764703),
                          (0.54509803921568623, 0.80000000000000004, 0.80000000000000004),
                          (0.5490196078431373, 0.80784313725490198, 0.80784313725490198),
                          (0.55294117647058827, 0.81960784313725488, 0.81960784313725488),
                          (0.55686274509803924, 0.82745098039215681, 0.82745098039215681),
                          (0.5607843137254902, 0.83921568627450982, 0.83921568627450982),
                          (0.56470588235294117, 0.85098039215686272, 0.85098039215686272),
                          (0.56862745098039214, 0.85882352941176465, 0.85882352941176465),
                          (0.5725490196078431, 0.87058823529411766, 0.87058823529411766),
                          (0.57647058823529407, 0.88235294117647056, 0.88235294117647056),
                          (0.58039215686274515, 0.8901960784313725, 0.8901960784313725),
                          (0.58431372549019611, 0.89803921568627454, 0.89803921568627454),
                          (0.58823529411764708, 0.90196078431372551, 0.90196078431372551),
                          (0.59215686274509804, 0.90196078431372551, 0.90196078431372551),
                          (0.59607843137254901, 0.90196078431372551, 0.90196078431372551),
                          (0.59999999999999998, 0.89803921568627454, 0.89803921568627454),
                          (0.60392156862745094, 0.8901960784313725, 0.8901960784313725),
                          (0.60784313725490191, 0.8784313725490196, 0.8784313725490196),
                          (0.61176470588235299, 0.8666666666666667, 0.8666666666666667),
                          (0.61568627450980395, 0.85098039215686272, 0.85098039215686272),
                          (0.61960784313725492, 0.82745098039215681, 0.82745098039215681),
                          (0.62352941176470589, 0.80392156862745101, 0.80392156862745101),
                          (0.62745098039215685, 0.77647058823529413, 0.77647058823529413),
                          (0.63137254901960782, 0.74901960784313726, 0.74901960784313726),
                          (0.63529411764705879, 0.72549019607843135, 0.72549019607843135),
                          (0.63921568627450975, 0.70588235294117652, 0.70588235294117652),
                          (0.64313725490196083, 0.68235294117647061, 0.68235294117647061),
                          (0.6470588235294118, 0.65098039215686276, 0.65098039215686276),
                          (0.65098039215686276, 0.62352941176470589, 0.62352941176470589),
                          (0.65490196078431373, 0.59215686274509804, 0.59215686274509804),
                          (0.6588235294117647, 0.5607843137254902, 0.5607843137254902),
                          (0.66274509803921566, 0.52941176470588236, 0.52941176470588236),
                          (0.66666666666666663, 0.49803921568627452, 0.49803921568627452),
                          (0.6705882352941176, 0.46666666666666667, 0.46666666666666667),
                          (0.67450980392156867, 0.43529411764705883, 0.43529411764705883),
                          (0.67843137254901964, 0.40392156862745099, 0.40392156862745099),
                          (0.68235294117647061, 0.37254901960784315, 0.37254901960784315),
                          (0.68627450980392157, 0.34901960784313724, 0.34901960784313724),
                          (0.69019607843137254, 0.32941176470588235, 0.32941176470588235),
                          (0.69411764705882351, 0.30588235294117649, 0.30588235294117649),
                          (0.69803921568627447, 0.27450980392156865, 0.27450980392156865),
                          (0.70196078431372544, 0.24705882352941178, 0.24705882352941178),
                          (0.70588235294117652, 0.21568627450980393, 0.21568627450980393),
                          (0.70980392156862748, 0.1803921568627451, 0.1803921568627451),
                          (0.71372549019607845, 0.14509803921568629, 0.14509803921568629),
                          (0.71764705882352942, 0.10980392156862745, 0.10980392156862745),
                          (0.72156862745098038, 0.074509803921568626, 0.074509803921568626),
                          (0.72549019607843135, 0.047058823529411764, 0.047058823529411764),
                          (0.72941176470588232, 0.035294117647058823, 0.035294117647058823),
                          (0.73333333333333328, 0.031372549019607843, 0.031372549019607843),
                          (0.73725490196078436, 0.031372549019607843, 0.031372549019607843),
                          (0.74117647058823533, 0.031372549019607843, 0.031372549019607843),
                          (0.74509803921568629, 0.031372549019607843, 0.031372549019607843),
                          (0.74901960784313726, 0.031372549019607843, 0.031372549019607843),
                          (0.75294117647058822, 0.031372549019607843, 0.031372549019607843),
                          (0.75686274509803919, 0.031372549019607843, 0.031372549019607843),
                          (0.76078431372549016, 0.031372549019607843, 0.031372549019607843),
                          (0.76470588235294112, 0.031372549019607843, 0.031372549019607843),
                          (0.7686274509803922, 0.031372549019607843, 0.031372549019607843),
                          (0.77254901960784317, 0.031372549019607843, 0.031372549019607843),
                          (0.77647058823529413, 0.031372549019607843, 0.031372549019607843),
                          (0.7803921568627451, 0.031372549019607843, 0.031372549019607843),
                          (0.78431372549019607, 0.031372549019607843, 0.031372549019607843),
                          (0.78823529411764703, 0.031372549019607843, 0.031372549019607843),
                          (0.792156862745098, 0.031372549019607843, 0.031372549019607843),
                          (0.79607843137254897, 0.050980392156862744, 0.050980392156862744),
                          (0.80000000000000004, 0.074509803921568626, 0.074509803921568626),
                          (0.80392156862745101, 0.086274509803921567, 0.086274509803921567),
                          (0.80784313725490198, 0.10196078431372549, 0.10196078431372549),
                          (0.81176470588235294, 0.11764705882352941, 0.11764705882352941),
                          (0.81568627450980391, 0.13333333333333333, 0.13333333333333333),
                          (0.81960784313725488, 0.14901960784313725, 0.14901960784313725),
                          (0.82352941176470584, 0.16078431372549021, 0.16078431372549021),
                          (0.82745098039215681, 0.17647058823529413, 0.17647058823529413),
                          (0.83137254901960789, 0.19215686274509805, 0.19215686274509805),
                          (0.83529411764705885, 0.20392156862745098, 0.20392156862745098),
                          (0.83921568627450982, 0.2196078431372549, 0.2196078431372549),
                          (0.84313725490196079, 0.23529411764705882, 0.23529411764705882),
                          (0.84705882352941175, 0.25098039215686274, 0.25098039215686274),
                          (0.85098039215686272, 0.2627450980392157, 0.2627450980392157),
                          (0.85490196078431369, 0.27843137254901962, 0.27843137254901962),
                          (0.85882352941176465, 0.29411764705882354, 0.29411764705882354),
                          (0.86274509803921573, 0.30980392156862746, 0.30980392156862746),
                          (0.8666666666666667, 0.32156862745098042, 0.32156862745098042),
                          (0.87058823529411766, 0.33333333333333331, 0.33333333333333331),
                          (0.87450980392156863, 0.34901960784313724, 0.34901960784313724),
                          (0.8784313725490196, 0.36470588235294116, 0.36470588235294116),
                          (0.88235294117647056, 0.38039215686274508, 0.38039215686274508),
                          (0.88627450980392153, 0.39215686274509803, 0.39215686274509803),
                          (0.8901960784313725, 0.40784313725490196, 0.40784313725490196),
                          (0.89411764705882357, 0.41960784313725491, 0.41960784313725491),
                          (0.89803921568627454, 0.4392156862745098, 0.4392156862745098),
                          (0.90196078431372551, 0.45098039215686275, 0.45098039215686275),
                          (0.90588235294117647, 0.46666666666666667, 0.46666666666666667),
                          (0.90980392156862744, 0.4823529411764706, 0.4823529411764706),
                          (0.9137254901960784, 0.49411764705882355, 0.49411764705882355),
                          (0.91764705882352937, 0.50980392156862742, 0.50980392156862742),
                          (0.92156862745098034, 0.52549019607843139, 0.52549019607843139),
                          (0.92549019607843142, 0.54117647058823526, 0.54117647058823526),
                          (0.92941176470588238, 0.55294117647058827, 0.55294117647058827),
                          (0.93333333333333335, 0.56862745098039214, 0.56862745098039214),
                          (0.93725490196078431, 0.58431372549019611, 0.58431372549019611),
                          (0.94117647058823528, 0.59999999999999998, 0.59999999999999998),
                          (0.94509803921568625, 0.61176470588235299, 0.61176470588235299),
                          (0.94901960784313721, 0.62352941176470589, 0.62352941176470589),
                          (0.95294117647058818, 0.63921568627450975, 0.63921568627450975),
                          (0.95686274509803926, 0.65490196078431373, 0.65490196078431373),
                          (0.96078431372549022, 0.66666666666666663, 0.66666666666666663),
                          (0.96470588235294119, 0.68235294117647061, 0.68235294117647061),
                          (0.96862745098039216, 0.69803921568627447, 0.69803921568627447),
                          (0.97254901960784312, 0.71372549019607845, 0.71372549019607845),
                          (0.97647058823529409, 0.72549019607843135, 0.72549019607843135),
                          (0.98039215686274506, 0.74117647058823533, 0.74117647058823533),
                          (0.98431372549019602, 0.75294117647058822, 0.75294117647058822),
                          (0.9882352941176471, 0.77254901960784317, 0.77254901960784317),
                          (0.99215686274509807, 0.78431372549019607, 0.78431372549019607),
                          (0.99607843137254903, 0.80000000000000004, 0.80000000000000004),
                          (1.0, 0.81176470588235294, 0.81176470588235294)],

                 'red': [ (0.0, 0.35686274509803922, 0.35686274509803922),
                          (0.0039215686274509803, 0.38823529411764707, 0.38823529411764707),
                          (0.0078431372549019607, 0.4392156862745098, 0.4392156862745098),
                          (0.011764705882352941, 0.48627450980392156, 0.48627450980392156),
                          (0.015686274509803921, 0.50980392156862742, 0.50980392156862742),
                          (0.019607843137254902, 0.50196078431372548, 0.50196078431372548),
                          (0.023529411764705882, 0.47058823529411764, 0.47058823529411764),
                          (0.027450980392156862, 0.43137254901960786, 0.43137254901960786),
                          (0.031372549019607843, 0.39215686274509803, 0.39215686274509803),
                          (0.035294117647058823, 0.34509803921568627, 0.34509803921568627),
                          (0.039215686274509803, 0.28627450980392155, 0.28627450980392155),
                          (0.043137254901960784, 0.21176470588235294, 0.21176470588235294),
                          (0.047058823529411764, 0.13333333333333333, 0.13333333333333333),
                          (0.050980392156862744, 0.070588235294117646, 0.070588235294117646),
                          (0.054901960784313725, 0.043137254901960784, 0.043137254901960784),
                          (0.058823529411764705, 0.031372549019607843, 0.031372549019607843),
                          (0.062745098039215685, 0.031372549019607843, 0.031372549019607843),
                          (0.066666666666666666, 0.031372549019607843, 0.031372549019607843),
                          (0.070588235294117646, 0.031372549019607843, 0.031372549019607843),
                          (0.074509803921568626, 0.031372549019607843, 0.031372549019607843),
                          (0.078431372549019607, 0.031372549019607843, 0.031372549019607843),
                          (0.082352941176470587, 0.031372549019607843, 0.031372549019607843),
                          (0.086274509803921567, 0.031372549019607843, 0.031372549019607843),
                          (0.090196078431372548, 0.031372549019607843, 0.031372549019607843),
                          (0.094117647058823528, 0.031372549019607843, 0.031372549019607843),
                          (0.098039215686274508, 0.031372549019607843, 0.031372549019607843),
                          (0.10196078431372549, 0.031372549019607843, 0.031372549019607843),
                          (0.10588235294117647, 0.031372549019607843, 0.031372549019607843),
                          (0.10980392156862745, 0.031372549019607843, 0.031372549019607843),
                          (0.11372549019607843, 0.031372549019607843, 0.031372549019607843),
                          (0.11764705882352941, 0.031372549019607843, 0.031372549019607843),
                          (0.12156862745098039, 0.031372549019607843, 0.031372549019607843),
                          (0.12549019607843137, 0.031372549019607843, 0.031372549019607843),
                          (0.12941176470588237, 0.031372549019607843, 0.031372549019607843),
                          (0.13333333333333333, 0.031372549019607843, 0.031372549019607843),
                          (0.13725490196078433, 0.031372549019607843, 0.031372549019607843),
                          (0.14117647058823529, 0.031372549019607843, 0.031372549019607843),
                          (0.14509803921568629, 0.031372549019607843, 0.031372549019607843),
                          (0.14901960784313725, 0.031372549019607843, 0.031372549019607843),
                          (0.15294117647058825, 0.031372549019607843, 0.031372549019607843),
                          (0.15686274509803921, 0.031372549019607843, 0.031372549019607843),
                          (0.16078431372549021, 0.031372549019607843, 0.031372549019607843),
                          (0.16470588235294117, 0.031372549019607843, 0.031372549019607843),
                          (0.16862745098039217, 0.031372549019607843, 0.031372549019607843),
                          (0.17254901960784313, 0.031372549019607843, 0.031372549019607843),
                          (0.17647058823529413, 0.031372549019607843, 0.031372549019607843),
                          (0.1803921568627451, 0.031372549019607843, 0.031372549019607843),
                          (0.18431372549019609, 0.031372549019607843, 0.031372549019607843),
                          (0.18823529411764706, 0.031372549019607843, 0.031372549019607843),
                          (0.19215686274509805, 0.031372549019607843, 0.031372549019607843),
                          (0.19607843137254902, 0.031372549019607843, 0.031372549019607843),
                          (0.20000000000000001, 0.031372549019607843, 0.031372549019607843),
                          (0.20392156862745098, 0.031372549019607843, 0.031372549019607843),
                          (0.20784313725490197, 0.031372549019607843, 0.031372549019607843),
                          (0.21176470588235294, 0.031372549019607843, 0.031372549019607843),
                          (0.21568627450980393, 0.031372549019607843, 0.031372549019607843),
                          (0.2196078431372549, 0.031372549019607843, 0.031372549019607843),
                          (0.22352941176470589, 0.031372549019607843, 0.031372549019607843),
                          (0.22745098039215686, 0.031372549019607843, 0.031372549019607843),
                          (0.23137254901960785, 0.031372549019607843, 0.031372549019607843),
                          (0.23529411764705882, 0.031372549019607843, 0.031372549019607843),
                          (0.23921568627450981, 0.031372549019607843, 0.031372549019607843),
                          (0.24313725490196078, 0.031372549019607843, 0.031372549019607843),
                          (0.24705882352941178, 0.031372549019607843, 0.031372549019607843),
                          (0.25098039215686274, 0.031372549019607843, 0.031372549019607843),
                          (0.25490196078431371, 0.031372549019607843, 0.031372549019607843),
                          (0.25882352941176473, 0.031372549019607843, 0.031372549019607843),
                          (0.2627450980392157, 0.031372549019607843, 0.031372549019607843),
                          (0.26666666666666666, 0.031372549019607843, 0.031372549019607843),
                          (0.27058823529411763, 0.031372549019607843, 0.031372549019607843),
                          (0.27450980392156865, 0.031372549019607843, 0.031372549019607843),
                          (0.27843137254901962, 0.031372549019607843, 0.031372549019607843),
                          (0.28235294117647058, 0.031372549019607843, 0.031372549019607843),
                          (0.28627450980392155, 0.031372549019607843, 0.031372549019607843),
                          (0.29019607843137257, 0.031372549019607843, 0.031372549019607843),
                          (0.29411764705882354, 0.031372549019607843, 0.031372549019607843),
                          (0.29803921568627451, 0.031372549019607843, 0.031372549019607843),
                          (0.30196078431372547, 0.031372549019607843, 0.031372549019607843),
                          (0.30588235294117649, 0.031372549019607843, 0.031372549019607843),
                          (0.30980392156862746, 0.031372549019607843, 0.031372549019607843),
                          (0.31372549019607843, 0.031372549019607843, 0.031372549019607843),
                          (0.31764705882352939, 0.031372549019607843, 0.031372549019607843),
                          (0.32156862745098042, 0.031372549019607843, 0.031372549019607843),
                          (0.32549019607843138, 0.031372549019607843, 0.031372549019607843),
                          (0.32941176470588235, 0.031372549019607843, 0.031372549019607843),
                          (0.33333333333333331, 0.031372549019607843, 0.031372549019607843),
                          (0.33725490196078434, 0.031372549019607843, 0.031372549019607843),
                          (0.3411764705882353, 0.031372549019607843, 0.031372549019607843),
                          (0.34509803921568627, 0.031372549019607843, 0.031372549019607843),
                          (0.34901960784313724, 0.031372549019607843, 0.031372549019607843),
                          (0.35294117647058826, 0.031372549019607843, 0.031372549019607843),
                          (0.35686274509803922, 0.031372549019607843, 0.031372549019607843),
                          (0.36078431372549019, 0.031372549019607843, 0.031372549019607843),
                          (0.36470588235294116, 0.031372549019607843, 0.031372549019607843),
                          (0.36862745098039218, 0.031372549019607843, 0.031372549019607843),
                          (0.37254901960784315, 0.031372549019607843, 0.031372549019607843),
                          (0.37647058823529411, 0.031372549019607843, 0.031372549019607843),
                          (0.38039215686274508, 0.031372549019607843, 0.031372549019607843),
                          (0.3843137254901961, 0.031372549019607843, 0.031372549019607843),
                          (0.38823529411764707, 0.031372549019607843, 0.031372549019607843),
                          (0.39215686274509803, 0.031372549019607843, 0.031372549019607843),
                          (0.396078431372549, 0.031372549019607843, 0.031372549019607843),
                          (0.40000000000000002, 0.031372549019607843, 0.031372549019607843),
                          (0.40392156862745099, 0.031372549019607843, 0.031372549019607843),
                          (0.40784313725490196, 0.031372549019607843, 0.031372549019607843),
                          (0.41176470588235292, 0.031372549019607843, 0.031372549019607843),
                          (0.41568627450980394, 0.031372549019607843, 0.031372549019607843),
                          (0.41960784313725491, 0.031372549019607843, 0.031372549019607843),
                          (0.42352941176470588, 0.031372549019607843, 0.031372549019607843),
                          (0.42745098039215684, 0.031372549019607843, 0.031372549019607843),
                          (0.43137254901960786, 0.031372549019607843, 0.031372549019607843),
                          (0.43529411764705883, 0.031372549019607843, 0.031372549019607843),
                          (0.4392156862745098, 0.031372549019607843, 0.031372549019607843),
                          (0.44313725490196076, 0.031372549019607843, 0.031372549019607843),
                          (0.44705882352941179, 0.031372549019607843, 0.031372549019607843),
                          (0.45098039215686275, 0.031372549019607843, 0.031372549019607843),
                          (0.45490196078431372, 0.031372549019607843, 0.031372549019607843),
                          (0.45882352941176469, 0.031372549019607843, 0.031372549019607843),
                          (0.46274509803921571, 0.031372549019607843, 0.031372549019607843),
                          (0.46666666666666667, 0.031372549019607843, 0.031372549019607843),
                          (0.47058823529411764, 0.031372549019607843, 0.031372549019607843),
                          (0.47450980392156861, 0.031372549019607843, 0.031372549019607843),
                          (0.47843137254901963, 0.031372549019607843, 0.031372549019607843),
                          (0.4823529411764706, 0.035294117647058823, 0.035294117647058823),
                          (0.48627450980392156, 0.054901960784313725, 0.054901960784313725),
                          (0.49019607843137253, 0.094117647058823528, 0.094117647058823528),
                          (0.49411764705882355, 0.14509803921568629, 0.14509803921568629),
                          (0.49803921568627452, 0.18431372549019609, 0.18431372549019609),
                          (0.50196078431372548, 0.20392156862745098, 0.20392156862745098),
                          (0.50588235294117645, 0.21568627450980393, 0.21568627450980393),
                          (0.50980392156862742, 0.22352941176470589, 0.22352941176470589),
                          (0.51372549019607838, 0.23529411764705882, 0.23529411764705882),
                          (0.51764705882352946, 0.25490196078431371, 0.25490196078431371),
                          (0.52156862745098043, 0.27450980392156865, 0.27450980392156865),
                          (0.52549019607843139, 0.29411764705882354, 0.29411764705882354),
                          (0.52941176470588236, 0.31372549019607843, 0.31372549019607843),
                          (0.53333333333333333, 0.33333333333333331, 0.33333333333333331),
                          (0.53725490196078429, 0.34901960784313724, 0.34901960784313724),
                          (0.54117647058823526, 0.37254901960784315, 0.37254901960784315),
                          (0.54509803921568623, 0.40392156862745099, 0.40392156862745099),
                          (0.5490196078431373, 0.43137254901960786, 0.43137254901960786),
                          (0.55294117647058827, 0.46274509803921571, 0.46274509803921571),
                          (0.55686274509803924, 0.49411764705882355, 0.49411764705882355),
                          (0.5607843137254902, 0.52549019607843139, 0.52549019607843139),
                          (0.56470588235294117, 0.55686274509803924, 0.55686274509803924),
                          (0.56862745098039214, 0.59215686274509804, 0.59215686274509804),
                          (0.5725490196078431, 0.62745098039215685, 0.62745098039215685),
                          (0.57647058823529407, 0.66274509803921566, 0.66274509803921566),
                          (0.58039215686274515, 0.70196078431372544, 0.70196078431372544),
                          (0.58431372549019611, 0.73333333333333328, 0.73333333333333328),
                          (0.58823529411764708, 0.76470588235294112, 0.76470588235294112),
                          (0.59215686274509804, 0.80392156862745101, 0.80392156862745101),
                          (0.59607843137254901, 0.84313725490196079, 0.84313725490196079),
                          (0.59999999999999998, 0.87450980392156863, 0.87450980392156863),
                          (0.60392156862745094, 0.88627450980392153, 0.88627450980392153),
                          (0.60784313725490191, 0.88235294117647056, 0.88235294117647056),
                          (0.61176470588235299, 0.8784313725490196, 0.8784313725490196),
                          (0.61568627450980395, 0.8784313725490196, 0.8784313725490196),
                          (0.61960784313725492, 0.8784313725490196, 0.8784313725490196),
                          (0.62352941176470589, 0.8784313725490196, 0.8784313725490196),
                          (0.62745098039215685, 0.8784313725490196, 0.8784313725490196),
                          (0.63137254901960782, 0.8784313725490196, 0.8784313725490196),
                          (0.63529411764705879, 0.8784313725490196, 0.8784313725490196),
                          (0.63921568627450975, 0.8784313725490196, 0.8784313725490196),
                          (0.64313725490196083, 0.8784313725490196, 0.8784313725490196),
                          (0.6470588235294118, 0.8784313725490196, 0.8784313725490196),
                          (0.65098039215686276, 0.8784313725490196, 0.8784313725490196),
                          (0.65490196078431373, 0.8784313725490196, 0.8784313725490196),
                          (0.6588235294117647, 0.8784313725490196, 0.8784313725490196),
                          (0.66274509803921566, 0.8784313725490196, 0.8784313725490196),
                          (0.66666666666666663, 0.8784313725490196, 0.8784313725490196),
                          (0.6705882352941176, 0.8784313725490196, 0.8784313725490196),
                          (0.67450980392156867, 0.8784313725490196, 0.8784313725490196),
                          (0.67843137254901964, 0.8784313725490196, 0.8784313725490196),
                          (0.68235294117647061, 0.8784313725490196, 0.8784313725490196),
                          (0.68627450980392157, 0.8784313725490196, 0.8784313725490196),
                          (0.69019607843137254, 0.8784313725490196, 0.8784313725490196),
                          (0.69411764705882351, 0.8784313725490196, 0.8784313725490196),
                          (0.69803921568627447, 0.87450980392156863, 0.87450980392156863),
                          (0.70196078431372544, 0.87450980392156863, 0.87450980392156863),
                          (0.70588235294117652, 0.87450980392156863, 0.87450980392156863),
                          (0.70980392156862748, 0.87058823529411766, 0.87058823529411766),
                          (0.71372549019607845, 0.8666666666666667, 0.8666666666666667),
                          (0.71764705882352942, 0.85882352941176465, 0.85882352941176465),
                          (0.72156862745098038, 0.84705882352941175, 0.84705882352941175),
                          (0.72549019607843135, 0.83529411764705885, 0.83529411764705885),
                          (0.72941176470588232, 0.81960784313725488, 0.81960784313725488),
                          (0.73333333333333328, 0.80000000000000004, 0.80000000000000004),
                          (0.73725490196078436, 0.78431372549019607, 0.78431372549019607),
                          (0.74117647058823533, 0.7686274509803922, 0.7686274509803922),
                          (0.74509803921568629, 0.74509803921568629, 0.74509803921568629),
                          (0.74901960784313726, 0.71372549019607845, 0.71372549019607845),
                          (0.75294117647058822, 0.68627450980392157, 0.68627450980392157),
                          (0.75686274509803919, 0.65490196078431373, 0.65490196078431373),
                          (0.76078431372549016, 0.62352941176470589, 0.62352941176470589),
                          (0.76470588235294112, 0.59215686274509804, 0.59215686274509804),
                          (0.7686274509803922, 0.5607843137254902, 0.5607843137254902),
                          (0.77254901960784317, 0.52941176470588236, 0.52941176470588236),
                          (0.77647058823529413, 0.49803921568627452, 0.49803921568627452),
                          (0.7803921568627451, 0.47058823529411764, 0.47058823529411764),
                          (0.78431372549019607, 0.45098039215686275, 0.45098039215686275),
                          (0.78823529411764703, 0.43137254901960786, 0.43137254901960786),
                          (0.792156862745098, 0.40784313725490196, 0.40784313725490196),
                          (0.79607843137254897, 0.42352941176470588, 0.42352941176470588),
                          (0.80000000000000004, 0.4392156862745098, 0.4392156862745098),
                          (0.80392156862745101, 0.44705882352941179, 0.44705882352941179),
                          (0.80784313725490198, 0.45098039215686275, 0.45098039215686275),
                          (0.81176470588235294, 0.46274509803921571, 0.46274509803921571),
                          (0.81568627450980391, 0.46666666666666667, 0.46666666666666667),
                          (0.81960784313725488, 0.47450980392156861, 0.47450980392156861),
                          (0.82352941176470584, 0.4823529411764706, 0.4823529411764706),
                          (0.82745098039215681, 0.49019607843137253, 0.49019607843137253),
                          (0.83137254901960789, 0.49803921568627452, 0.49803921568627452),
                          (0.83529411764705885, 0.50196078431372548, 0.50196078431372548),
                          (0.83921568627450982, 0.50980392156862742, 0.50980392156862742),
                          (0.84313725490196079, 0.51764705882352946, 0.51764705882352946),
                          (0.84705882352941175, 0.52549019607843139, 0.52549019607843139),
                          (0.85098039215686272, 0.53333333333333333, 0.53333333333333333),
                          (0.85490196078431369, 0.54117647058823526, 0.54117647058823526),
                          (0.85882352941176465, 0.5490196078431373, 0.5490196078431373),
                          (0.86274509803921573, 0.55686274509803924, 0.55686274509803924),
                          (0.8666666666666667, 0.5607843137254902, 0.5607843137254902),
                          (0.87058823529411766, 0.56862745098039214, 0.56862745098039214),
                          (0.87450980392156863, 0.57647058823529407, 0.57647058823529407),
                          (0.8784313725490196, 0.58431372549019611, 0.58431372549019611),
                          (0.88235294117647056, 0.58823529411764708, 0.58823529411764708),
                          (0.88627450980392153, 0.59607843137254901, 0.59607843137254901),
                          (0.8901960784313725, 0.60392156862745094, 0.60392156862745094),
                          (0.89411764705882357, 0.61176470588235299, 0.61176470588235299),
                          (0.89803921568627454, 0.61960784313725492, 0.61960784313725492),
                          (0.90196078431372551, 0.62745098039215685, 0.62745098039215685),
                          (0.90588235294117647, 0.63529411764705879, 0.63529411764705879),
                          (0.90980392156862744, 0.64313725490196083, 0.64313725490196083),
                          (0.9137254901960784, 0.6470588235294118, 0.6470588235294118),
                          (0.91764705882352937, 0.6588235294117647, 0.6588235294117647),
                          (0.92156862745098034, 0.66274509803921566, 0.66274509803921566),
                          (0.92549019607843142, 0.66666666666666663, 0.66666666666666663),
                          (0.92941176470588238, 0.67843137254901964, 0.67843137254901964),
                          (0.93333333333333335, 0.68627450980392157, 0.68627450980392157),
                          (0.93725490196078431, 0.69019607843137254, 0.69019607843137254),
                          (0.94117647058823528, 0.69803921568627447, 0.69803921568627447),
                          (0.94509803921568625, 0.70980392156862748, 0.70980392156862748),
                          (0.94901960784313721, 0.71372549019607845, 0.71372549019607845),
                          (0.95294117647058818, 0.71764705882352942, 0.71764705882352942),
                          (0.95686274509803926, 0.72549019607843135, 0.72549019607843135),
                          (0.96078431372549022, 0.73725490196078436, 0.73725490196078436),
                          (0.96470588235294119, 0.74509803921568629, 0.74509803921568629),
                          (0.96862745098039216, 0.74901960784313726, 0.74901960784313726),
                          (0.97254901960784312, 0.75686274509803919, 0.75686274509803919),
                          (0.97647058823529409, 0.76470588235294112, 0.76470588235294112),
                          (0.98039215686274506, 0.77254901960784317, 0.77254901960784317),
                          (0.98431372549019602, 0.7803921568627451, 0.7803921568627451),
                          (0.9882352941176471, 0.78431372549019607, 0.78431372549019607),
                          (0.99215686274509807, 0.792156862745098, 0.792156862745098),
                          (0.99607843137254903, 0.80000000000000004, 0.80000000000000004),
                          (1.0, 0.80784313725490198, 0.80784313725490198)]
                }

            self.cmap = LinearSegmentedColormap('sst', sst_cmdata)

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



class VegetationIndexProdData:
    """
    This class contains static data for the VIIRS VI Products.
    """

    #ViirsVIqualBitMasks = ((3, 12,  48, 64 , 128        ),
                            #(1, 2, 12, 48, 64, 128       ),
                            #(1, 2, 4, 8, 16, 32, 64, 128 ),
                            #(1, 2, 4, 8, 16, 32, 64, 128 ))

    #ViirsVIqualBitShift = ((0, 2, 4, 6, 7         ),
                            #(0, 1, 2, 4, 6, 7      ),
                            #(0, 1, 2, 3, 4, 5, 6, 7),
                            #(0, 1, 2, 3, 4, 5, 6, 7))

    #ViirsVIqualBitMaskNames = (('Skin SST Quality','Bulk SST Quality','SST State','Algorithm','Day / Night'),
                                 #('Bad LWIR Pixel','Bad SWIR Pixel','Cloud Confidence','Adjacent Pixel Cloud Confident Value','Thin Cirrus','Sea Ice'),
                                 #('Sun Glint','Exclusion, AOT > 1','Degraded, AOT > 0.6','Exclusion, Not Ocean','Degraded, HCS limit','Degraded, Sensor Zenith Angle > 40','Skin SST Outside Range','Bulk SST Outside Range'),
                                 #('Skin SST Degraded, T > 305 K','Bulk SST Degraded, T > 305 K','Spare','Spare','Spare','Spare','Spare','Spare'),
                                #)
    
    #ViirsVIqualvalues = ( (  (0,1,2,3), (0,1,2,3), (0,1,2), (0,1), (0,1)             ),      # Byte 0
                           #(  (0,1), (0,1), (0,1,2,3), (0,1,2,3), (0,1), (0,1)        ),      # Byte 1
                           #(  (0,1), (0,1), (0,1), (0,1), (0,1), (0,1), (0,1), (0,1)  ),      # Byte 2
                           #(  (0,1), (0,1), (0,1), (0,1), (0,1), (0,1), (0,1), (0,1)  )       # Byte 3
                          #)

    #ViirsVIqualFillBoundaries = [
                                   #[ [-0.5, 0.5, 1.5, 2.5, 3.5],[-0.5, 0.5, 1.5, 2.5, 3.5],[-0.5, 0.5, 1.5, 2.5],[-0.5, 0.5, 1.5],[-0.5, 0.5, 1.5]                          ], # Byte 0
                                   #[ [-0.5, 0.5, 1.5],[-0.5, 0.5, 1.5],[-0.5, 0.5, 1.5, 2.5, 3.5],[-0.5, 0.5, 1.5, 2.5, 3.5], [-0.5, 0.5, 1.5],[-0.5, 0.5, 1.5]             ], # Byte 1
                                   #[ [-0.5, 0.5, 1.5],[-0.5, 0.5, 1.5],[-0.5, 0.5, 1.5],[-0.5, 0.5, 1.5],[-0.5, 0.5, 1.5],[-0.5, 0.5, 1.5],[-0.5, 0.5, 1.5],[-0.5, 0.5, 1.5]], # Byte 3
                                   #[ [-0.5, 0.5, 1.5],[-0.5, 0.5, 1.5],[-0.5, 0.5, 1.5],[-0.5, 0.5, 1.5],[-0.5, 0.5, 1.5],[-0.5, 0.5, 1.5],[-0.5, 0.5, 1.5],[-0.5, 0.5, 1.5]] # Byte 4
                                 #]

    #ViirsVIqualFillColours = [
                                #[ ['k', 'r', 'b', '#00ff00'],['k', 'r', 'b', '#00ff00'],['brown','green','yellow'],['k','g'],['k','yellow']       ], # Byte 0
                                #[ ['g','k'],['g','k'],['#00ff00','#00ffff','#ff0000','w'],['#00ff00','#00ffff','#ff0000','w'],['k','g'],['k','g'] ], # Byte 2
                                #[ ['k','g'],['k','g'],['k','g'],['b','k'],['g','k'],['k','g'],['g','k'],['g','k']                                 ], # Byte 3
                                #[ ['k','g'],['k','g'],['k','g'],['k','g'],['k','g'],['k','g'],['k','g'],['k','g']                                 ] # Byte 4
                              #]

    ## The tick labels of the colourbar categories
    #ViirsVIqualTickNames = [
                              #[ ['Not retrieved','Excluded','Degraded','High Quality'],['Not retrieved','Excluded','Degraded','High Quality'],[' Dry / None','Moist','Average'],['Non-linear Split Window','Triple Window'],['Night','Day'] ],      # Byte 0
                              #[ ['Good SDR','Bad SDR'],['Good SDR','Bad SDR'],['Confident Clear','Probably Clear','Probably Cloudy','Confident Cloudy'],['Confident Clear','Probably Clear','Probably Cloudy','Confident Cloudy'],['No Thin Cirrus','Thin Cirrus'],['No Sea Ice','Sea Ice'] ],      # Byte 1
                              #[ ['No sun glint','Sun glint'],['No','Yes'],['No','Yes'],['Ocean','Not ocean'],['Within HCS limit','Past HCS limit'],['No','Yes'],['In range','Out of range'],['In range','Out of range'] ],      # Byte 3
                              #[ ['Not degraded','Degraded'],['Not degraded','Degraded'],['',''],['',''],['',''],['',''],['',''],['',''] ]       # Byte 4
                            #]
                              


    class VegetationIndexProd:
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

            ndvi_cmdata = {
                    'blue': [ (0.0, 0.92156862745098034, 0.92156862745098034),
                              (0.0039215686274509803, 0.92156862745098034, 0.92156862745098034),
                              (0.0078431372549019607, 0.92156862745098034, 0.92156862745098034),
                              (0.011764705882352941, 0.92156862745098034, 0.92156862745098034),
                              (0.015686274509803921, 0.92156862745098034, 0.92156862745098034),
                              (0.019607843137254902, 0.92156862745098034, 0.92156862745098034),
                              (0.023529411764705882, 0.92156862745098034, 0.92156862745098034),
                              (0.027450980392156862, 0.92156862745098034, 0.92156862745098034),
                              (0.031372549019607843, 0.92156862745098034, 0.92156862745098034),
                              (0.035294117647058823, 0.92156862745098034, 0.92156862745098034),
                              (0.039215686274509803, 0.92156862745098034, 0.92156862745098034),
                              (0.043137254901960784, 0.92156862745098034, 0.92156862745098034),
                              (0.047058823529411764, 0.92156862745098034, 0.92156862745098034),
                              (0.050980392156862744, 0.92156862745098034, 0.92156862745098034),
                              (0.054901960784313725, 0.92156862745098034, 0.92156862745098034),
                              (0.058823529411764705, 0.92156862745098034, 0.92156862745098034),
                              (0.062745098039215685, 0.12549019607843137, 0.12549019607843137),
                              (0.066666666666666666, 0.13333333333333333, 0.13333333333333333),
                              (0.070588235294117646, 0.14117647058823529, 0.14117647058823529),
                              (0.074509803921568626, 0.14901960784313725, 0.14901960784313725),
                              (0.078431372549019607, 0.15294117647058825, 0.15294117647058825),
                              (0.082352941176470587, 0.16078431372549021, 0.16078431372549021),
                              (0.086274509803921567, 0.16862745098039217, 0.16862745098039217),
                              (0.090196078431372548, 0.17647058823529413, 0.17647058823529413),
                              (0.094117647058823528, 0.18431372549019609, 0.18431372549019609),
                              (0.098039215686274508, 0.19215686274509805, 0.19215686274509805),
                              (0.10196078431372549, 0.19607843137254902, 0.19607843137254902),
                              (0.10588235294117647, 0.20392156862745098, 0.20392156862745098),
                              (0.10980392156862745, 0.21176470588235294, 0.21176470588235294),
                              (0.11372549019607843, 0.2196078431372549, 0.2196078431372549),
                              (0.11764705882352941, 0.22745098039215686, 0.22745098039215686),
                              (0.12156862745098039, 0.23137254901960785, 0.23137254901960785),
                              (0.12549019607843137, 0.23921568627450981, 0.23921568627450981),
                              (0.12941176470588237, 0.24705882352941178, 0.24705882352941178),
                              (0.13333333333333333, 0.25098039215686274, 0.25098039215686274),
                              (0.13725490196078433, 0.25882352941176473, 0.25882352941176473),
                              (0.14117647058823529, 0.2627450980392157, 0.2627450980392157),
                              (0.14509803921568629, 0.27058823529411763, 0.27058823529411763),
                              (0.14901960784313725, 0.27843137254901962, 0.27843137254901962),
                              (0.15294117647058825, 0.28235294117647058, 0.28235294117647058),
                              (0.15686274509803921, 0.29019607843137257, 0.29019607843137257),
                              (0.16078431372549021, 0.29411764705882354, 0.29411764705882354),
                              (0.16470588235294117, 0.29803921568627451, 0.29803921568627451),
                              (0.16862745098039217, 0.30588235294117649, 0.30588235294117649),
                              (0.17254901960784313, 0.30980392156862746, 0.30980392156862746),
                              (0.17647058823529413, 0.31764705882352939, 0.31764705882352939),
                              (0.1803921568627451, 0.32156862745098042, 0.32156862745098042),
                              (0.18431372549019609, 0.32549019607843138, 0.32549019607843138),
                              (0.18823529411764706, 0.33333333333333331, 0.33333333333333331),
                              (0.19215686274509805, 0.33725490196078434, 0.33725490196078434),
                              (0.19607843137254902, 0.3411764705882353, 0.3411764705882353),
                              (0.20000000000000001, 0.34509803921568627, 0.34509803921568627),
                              (0.20392156862745098, 0.34901960784313724, 0.34901960784313724),
                              (0.20784313725490197, 0.35294117647058826, 0.35294117647058826),
                              (0.21176470588235294, 0.35686274509803922, 0.35686274509803922),
                              (0.21568627450980393, 0.36078431372549019, 0.36078431372549019),
                              (0.2196078431372549, 0.36470588235294116, 0.36470588235294116),
                              (0.22352941176470589, 0.36862745098039218, 0.36862745098039218),
                              (0.22745098039215686, 0.37254901960784315, 0.37254901960784315),
                              (0.23137254901960785, 0.37647058823529411, 0.37647058823529411),
                              (0.23529411764705882, 0.38039215686274508, 0.38039215686274508),
                              (0.23921568627450981, 0.3843137254901961, 0.3843137254901961),
                              (0.24313725490196078, 0.38823529411764707, 0.38823529411764707),
                              (0.24705882352941178, 0.38823529411764707, 0.38823529411764707),
                              (0.25098039215686274, 0.39215686274509803, 0.39215686274509803),
                              (0.25490196078431371, 0.396078431372549, 0.396078431372549),
                              (0.25882352941176473, 0.396078431372549, 0.396078431372549),
                              (0.2627450980392157, 0.40000000000000002, 0.40000000000000002),
                              (0.26666666666666666, 0.40000000000000002, 0.40000000000000002),
                              (0.27058823529411763, 0.40392156862745099, 0.40392156862745099),
                              (0.27450980392156865, 0.40392156862745099, 0.40392156862745099),
                              (0.27843137254901962, 0.40392156862745099, 0.40392156862745099),
                              (0.28235294117647058, 0.40784313725490196, 0.40784313725490196),
                              (0.28627450980392155, 0.40784313725490196, 0.40784313725490196),
                              (0.29019607843137257, 0.40784313725490196, 0.40784313725490196),
                              (0.29411764705882354, 0.17254901960784313, 0.17254901960784313),
                              (0.29803921568627451, 0.17254901960784313, 0.17254901960784313),
                              (0.30196078431372547, 0.17254901960784313, 0.17254901960784313),
                              (0.30588235294117649, 0.17254901960784313, 0.17254901960784313),
                              (0.30980392156862746, 0.17254901960784313, 0.17254901960784313),
                              (0.31372549019607843, 0.17254901960784313, 0.17254901960784313),
                              (0.31764705882352939, 0.17254901960784313, 0.17254901960784313),
                              (0.32156862745098042, 0.17254901960784313, 0.17254901960784313),
                              (0.32549019607843138, 0.17254901960784313, 0.17254901960784313),
                              (0.32941176470588235, 0.17254901960784313, 0.17254901960784313),
                              (0.33333333333333331, 0.17254901960784313, 0.17254901960784313),
                              (0.33725490196078434, 0.17254901960784313, 0.17254901960784313),
                              (0.3411764705882353, 0.17254901960784313, 0.17254901960784313),
                              (0.34509803921568627, 0.17254901960784313, 0.17254901960784313),
                              (0.34901960784313724, 0.17254901960784313, 0.17254901960784313),
                              (0.35294117647058826, 0.074509803921568626, 0.074509803921568626),
                              (0.35686274509803922, 0.074509803921568626, 0.074509803921568626),
                              (0.36078431372549019, 0.074509803921568626, 0.074509803921568626),
                              (0.36470588235294116, 0.074509803921568626, 0.074509803921568626),
                              (0.36862745098039218, 0.074509803921568626, 0.074509803921568626),
                              (0.37254901960784315, 0.074509803921568626, 0.074509803921568626),
                              (0.37647058823529411, 0.074509803921568626, 0.074509803921568626),
                              (0.38039215686274508, 0.074509803921568626, 0.074509803921568626),
                              (0.3843137254901961, 0.074509803921568626, 0.074509803921568626),
                              (0.38823529411764707, 0.074509803921568626, 0.074509803921568626),
                              (0.39215686274509803, 0.074509803921568626, 0.074509803921568626),
                              (0.396078431372549, 0.074509803921568626, 0.074509803921568626),
                              (0.40000000000000002, 0.074509803921568626, 0.074509803921568626),
                              (0.40392156862745099, 0.074509803921568626, 0.074509803921568626),
                              (0.40784313725490196, 0.074509803921568626, 0.074509803921568626),
                              (0.41176470588235292, 0.0, 0.0),
                              (0.41568627450980394, 0.0, 0.0),
                              (0.41960784313725491, 0.0, 0.0),
                              (0.42352941176470588, 0.0, 0.0),
                              (0.42745098039215684, 0.0, 0.0),
                              (0.43137254901960786, 0.0, 0.0),
                              (0.43529411764705883, 0.0, 0.0),
                              (0.4392156862745098, 0.0, 0.0),
                              (0.44313725490196076, 0.0, 0.0),
                              (0.44705882352941179, 0.0, 0.0),
                              (0.45098039215686275, 0.0, 0.0),
                              (0.45490196078431372, 0.0, 0.0),
                              (0.45882352941176469, 0.0, 0.0),
                              (0.46274509803921571, 0.0, 0.0),
                              (0.46666666666666667, 0.0, 0.0),
                              (0.47058823529411764, 0.0, 0.0),
                              (0.47450980392156861, 0.0, 0.0),
                              (0.47843137254901963, 0.0, 0.0),
                              (0.4823529411764706, 0.0, 0.0),
                              (0.48627450980392156, 0.0, 0.0),
                              (0.49019607843137253, 0.0, 0.0),
                              (0.49411764705882355, 0.0, 0.0),
                              (0.49803921568627452, 0.0, 0.0),
                              (0.50196078431372548, 0.0, 0.0),
                              (0.50588235294117645, 0.0, 0.0),
                              (0.50980392156862742, 0.0, 0.0),
                              (0.51372549019607838, 0.0, 0.0),
                              (0.51764705882352946, 0.0, 0.0),
                              (0.52156862745098043, 0.0, 0.0),
                              (0.52549019607843139, 0.0, 0.0),
                              (0.52941176470588236, 0.0, 0.0),
                              (0.53333333333333333, 0.0, 0.0),
                              (0.53725490196078429, 0.0, 0.0),
                              (0.54117647058823526, 0.0, 0.0),
                              (0.54509803921568623, 0.0, 0.0),
                              (0.5490196078431373, 0.0, 0.0),
                              (0.55294117647058827, 0.0, 0.0),
                              (0.55686274509803924, 0.0, 0.0),
                              (0.5607843137254902, 0.0, 0.0),
                              (0.56470588235294117, 0.0, 0.0),
                              (0.56862745098039214, 0.0, 0.0),
                              (0.5725490196078431, 0.0, 0.0),
                              (0.57647058823529407, 0.0, 0.0),
                              (0.58039215686274515, 0.0, 0.0),
                              (0.58431372549019611, 0.0, 0.0),
                              (0.58823529411764708, 0.0039215686274509803, 0.0039215686274509803),
                              (0.59215686274509804, 0.0039215686274509803, 0.0039215686274509803),
                              (0.59607843137254901, 0.0039215686274509803, 0.0039215686274509803),
                              (0.59999999999999998, 0.0039215686274509803, 0.0039215686274509803),
                              (0.60392156862745094, 0.0039215686274509803, 0.0039215686274509803),
                              (0.60784313725490191, 0.0039215686274509803, 0.0039215686274509803),
                              (0.61176470588235299, 0.0039215686274509803, 0.0039215686274509803),
                              (0.61568627450980395, 0.0039215686274509803, 0.0039215686274509803),
                              (0.61960784313725492, 0.0039215686274509803, 0.0039215686274509803),
                              (0.62352941176470589, 0.0039215686274509803, 0.0039215686274509803),
                              (0.62745098039215685, 0.0039215686274509803, 0.0039215686274509803),
                              (0.63137254901960782, 0.0039215686274509803, 0.0039215686274509803),
                              (0.63529411764705879, 0.0039215686274509803, 0.0039215686274509803),
                              (0.63921568627450975, 0.0039215686274509803, 0.0039215686274509803),
                              (0.64313725490196083, 0.0039215686274509803, 0.0039215686274509803),
                              (0.6470588235294118, 0.0039215686274509803, 0.0039215686274509803),
                              (0.65098039215686276, 0.0039215686274509803, 0.0039215686274509803),
                              (0.65490196078431373, 0.0039215686274509803, 0.0039215686274509803),
                              (0.6588235294117647, 0.0039215686274509803, 0.0039215686274509803),
                              (0.66274509803921566, 0.0039215686274509803, 0.0039215686274509803),
                              (0.66666666666666663, 0.0039215686274509803, 0.0039215686274509803),
                              (0.6705882352941176, 0.0039215686274509803, 0.0039215686274509803),
                              (0.67450980392156867, 0.0039215686274509803, 0.0039215686274509803),
                              (0.67843137254901964, 0.0039215686274509803, 0.0039215686274509803),
                              (0.68235294117647061, 0.0039215686274509803, 0.0039215686274509803),
                              (0.68627450980392157, 0.0039215686274509803, 0.0039215686274509803),
                              (0.69019607843137254, 0.0039215686274509803, 0.0039215686274509803),
                              (0.69411764705882351, 0.0039215686274509803, 0.0039215686274509803),
                              (0.69803921568627447, 0.0039215686274509803, 0.0039215686274509803),
                              (0.70196078431372544, 0.0039215686274509803, 0.0039215686274509803),
                              (0.70588235294117652, 0.0039215686274509803, 0.0039215686274509803),
                              (0.70980392156862748, 0.0039215686274509803, 0.0039215686274509803),
                              (0.71372549019607845, 0.0039215686274509803, 0.0039215686274509803),
                              (0.71764705882352942, 0.0039215686274509803, 0.0039215686274509803),
                              (0.72156862745098038, 0.0039215686274509803, 0.0039215686274509803),
                              (0.72549019607843135, 0.0039215686274509803, 0.0039215686274509803),
                              (0.72941176470588232, 0.0039215686274509803, 0.0039215686274509803),
                              (0.73333333333333328, 0.0039215686274509803, 0.0039215686274509803),
                              (0.73725490196078436, 0.0039215686274509803, 0.0039215686274509803),
                              (0.74117647058823533, 0.0039215686274509803, 0.0039215686274509803),
                              (0.74509803921568629, 0.0039215686274509803, 0.0039215686274509803),
                              (0.74901960784313726, 0.0039215686274509803, 0.0039215686274509803),
                              (0.75294117647058822, 0.0039215686274509803, 0.0039215686274509803),
                              (0.75686274509803919, 0.0039215686274509803, 0.0039215686274509803),
                              (0.76078431372549016, 0.0039215686274509803, 0.0039215686274509803),
                              (0.76470588235294112, 0.0, 0.0),
                              (0.7686274509803922, 0.0, 0.0),
                              (0.77254901960784317, 0.0, 0.0),
                              (0.77647058823529413, 0.0, 0.0),
                              (0.7803921568627451, 0.0, 0.0),
                              (0.78431372549019607, 0.0, 0.0),
                              (0.78823529411764703, 0.0, 0.0),
                              (0.792156862745098, 0.0, 0.0),
                              (0.79607843137254897, 0.0, 0.0),
                              (0.80000000000000004, 0.0, 0.0),
                              (0.80392156862745101, 0.0, 0.0),
                              (0.80784313725490198, 0.0, 0.0),
                              (0.81176470588235294, 0.0, 0.0),
                              (0.81568627450980391, 0.0, 0.0),
                              (0.81960784313725488, 0.0, 0.0),
                              (0.82352941176470584, 0.0039215686274509803, 0.0039215686274509803),
                              (0.82745098039215681, 0.0039215686274509803, 0.0039215686274509803),
                              (0.83137254901960789, 0.0039215686274509803, 0.0039215686274509803),
                              (0.83529411764705885, 0.0039215686274509803, 0.0039215686274509803),
                              (0.83921568627450982, 0.0039215686274509803, 0.0039215686274509803),
                              (0.84313725490196079, 0.0039215686274509803, 0.0039215686274509803),
                              (0.84705882352941175, 0.0039215686274509803, 0.0039215686274509803),
                              (0.85098039215686272, 0.0039215686274509803, 0.0039215686274509803),
                              (0.85490196078431369, 0.0039215686274509803, 0.0039215686274509803),
                              (0.85882352941176465, 0.0039215686274509803, 0.0039215686274509803),
                              (0.86274509803921573, 0.0039215686274509803, 0.0039215686274509803),
                              (0.8666666666666667, 0.0039215686274509803, 0.0039215686274509803),
                              (0.87058823529411766, 0.0039215686274509803, 0.0039215686274509803),
                              (0.87450980392156863, 0.0039215686274509803, 0.0039215686274509803),
                              (0.8784313725490196, 0.0039215686274509803, 0.0039215686274509803),
                              (0.88235294117647056, 0.0039215686274509803, 0.0039215686274509803),
                              (0.88627450980392153, 0.0039215686274509803, 0.0039215686274509803),
                              (0.8901960784313725, 0.0039215686274509803, 0.0039215686274509803),
                              (0.89411764705882357, 0.0039215686274509803, 0.0039215686274509803),
                              (0.89803921568627454, 0.0039215686274509803, 0.0039215686274509803),
                              (0.90196078431372551, 0.0039215686274509803, 0.0039215686274509803),
                              (0.90588235294117647, 0.0039215686274509803, 0.0039215686274509803),
                              (0.90980392156862744, 0.0039215686274509803, 0.0039215686274509803),
                              (0.9137254901960784, 0.0039215686274509803, 0.0039215686274509803),
                              (0.91764705882352937, 0.0039215686274509803, 0.0039215686274509803),
                              (0.92156862745098034, 0.0039215686274509803, 0.0039215686274509803),
                              (0.92549019607843142, 0.0039215686274509803, 0.0039215686274509803),
                              (0.92941176470588238, 0.0039215686274509803, 0.0039215686274509803),
                              (0.93333333333333335, 0.0039215686274509803, 0.0039215686274509803),
                              (0.93725490196078431, 0.0039215686274509803, 0.0039215686274509803),
                              (0.94117647058823528, 0.0039215686274509803, 0.0039215686274509803),
                              (0.94509803921568625, 0.0039215686274509803, 0.0039215686274509803),
                              (0.94901960784313721, 0.0039215686274509803, 0.0039215686274509803),
                              (0.95294117647058818, 0.0039215686274509803, 0.0039215686274509803),
                              (0.95686274509803926, 0.0039215686274509803, 0.0039215686274509803),
                              (0.96078431372549022, 0.0039215686274509803, 0.0039215686274509803),
                              (0.96470588235294119, 0.0039215686274509803, 0.0039215686274509803),
                              (0.96862745098039216, 0.0039215686274509803, 0.0039215686274509803),
                              (0.97254901960784312, 0.0039215686274509803, 0.0039215686274509803),
                              (0.97647058823529409, 0.0039215686274509803, 0.0039215686274509803),
                              (0.98039215686274506, 0.0039215686274509803, 0.0039215686274509803),
                              (0.98431372549019602, 0.0039215686274509803, 0.0039215686274509803),
                              (0.9882352941176471, 0.0039215686274509803, 0.0039215686274509803),
                              (0.99215686274509807, 0.0039215686274509803, 0.0039215686274509803),
                              (0.99607843137254903, 0.0039215686274509803, 0.0039215686274509803),
                              (1.0, 0.0039215686274509803, 0.0039215686274509803)],

                    'green': [(0.0, 0.92156862745098034, 0.92156862745098034),
                              (0.0039215686274509803, 0.92156862745098034, 0.92156862745098034),
                              (0.0078431372549019607, 0.92156862745098034, 0.92156862745098034),
                              (0.011764705882352941, 0.92156862745098034, 0.92156862745098034),
                              (0.015686274509803921, 0.92156862745098034, 0.92156862745098034),
                              (0.019607843137254902, 0.92156862745098034, 0.92156862745098034),
                              (0.023529411764705882, 0.92156862745098034, 0.92156862745098034),
                              (0.027450980392156862, 0.92156862745098034, 0.92156862745098034),
                              (0.031372549019607843, 0.92156862745098034, 0.92156862745098034),
                              (0.035294117647058823, 0.92156862745098034, 0.92156862745098034),
                              (0.039215686274509803, 0.92156862745098034, 0.92156862745098034),
                              (0.043137254901960784, 0.92156862745098034, 0.92156862745098034),
                              (0.047058823529411764, 0.92156862745098034, 0.92156862745098034),
                              (0.050980392156862744, 0.92156862745098034, 0.92156862745098034),
                              (0.054901960784313725, 0.92156862745098034, 0.92156862745098034),
                              (0.058823529411764705, 0.92156862745098034, 0.92156862745098034),
                              (0.062745098039215685, 0.42352941176470588, 0.42352941176470588),
                              (0.066666666666666666, 0.43137254901960786, 0.43137254901960786),
                              (0.070588235294117646, 0.44313725490196076, 0.44313725490196076),
                              (0.074509803921568626, 0.45098039215686275, 0.45098039215686275),
                              (0.078431372549019607, 0.46274509803921571, 0.46274509803921571),
                              (0.082352941176470587, 0.47058823529411764, 0.47058823529411764),
                              (0.086274509803921567, 0.4823529411764706, 0.4823529411764706),
                              (0.090196078431372548, 0.49019607843137253, 0.49019607843137253),
                              (0.094117647058823528, 0.50196078431372548, 0.50196078431372548),
                              (0.098039215686274508, 0.50980392156862742, 0.50980392156862742),
                              (0.10196078431372549, 0.52156862745098043, 0.52156862745098043),
                              (0.10588235294117647, 0.52941176470588236, 0.52941176470588236),
                              (0.10980392156862745, 0.53725490196078429, 0.53725490196078429),
                              (0.11372549019607843, 0.5490196078431373, 0.5490196078431373),
                              (0.11764705882352941, 0.55686274509803924, 0.55686274509803924),
                              (0.12156862745098039, 0.56470588235294117, 0.56470588235294117),
                              (0.12549019607843137, 0.57647058823529407, 0.57647058823529407),
                              (0.12941176470588237, 0.58431372549019611, 0.58431372549019611),
                              (0.13333333333333333, 0.59215686274509804, 0.59215686274509804),
                              (0.13725490196078433, 0.60392156862745094, 0.60392156862745094),
                              (0.14117647058823529, 0.61176470588235299, 0.61176470588235299),
                              (0.14509803921568629, 0.61960784313725492, 0.61960784313725492),
                              (0.14901960784313725, 0.62745098039215685, 0.62745098039215685),
                              (0.15294117647058825, 0.63529411764705879, 0.63529411764705879),
                              (0.15686274509803921, 0.6470588235294118, 0.6470588235294118),
                              (0.16078431372549021, 0.65490196078431373, 0.65490196078431373),
                              (0.16470588235294117, 0.66274509803921566, 0.66274509803921566),
                              (0.16862745098039217, 0.6705882352941176, 0.6705882352941176),
                              (0.17254901960784313, 0.67843137254901964, 0.67843137254901964),
                              (0.17647058823529413, 0.68627450980392157, 0.68627450980392157),
                              (0.1803921568627451, 0.69411764705882351, 0.69411764705882351),
                              (0.18431372549019609, 0.70196078431372544, 0.70196078431372544),
                              (0.18823529411764706, 0.70980392156862748, 0.70980392156862748),
                              (0.19215686274509805, 0.71764705882352942, 0.71764705882352942),
                              (0.19607843137254902, 0.72549019607843135, 0.72549019607843135),
                              (0.20000000000000001, 0.73333333333333328, 0.73333333333333328),
                              (0.20392156862745098, 0.74117647058823533, 0.74117647058823533),
                              (0.20784313725490197, 0.74901960784313726, 0.74901960784313726),
                              (0.21176470588235294, 0.75294117647058822, 0.75294117647058822),
                              (0.21568627450980393, 0.76078431372549016, 0.76078431372549016),
                              (0.2196078431372549, 0.7686274509803922, 0.7686274509803922),
                              (0.22352941176470589, 0.77647058823529413, 0.77647058823529413),
                              (0.22745098039215686, 0.7803921568627451, 0.7803921568627451),
                              (0.23137254901960785, 0.78823529411764703, 0.78823529411764703),
                              (0.23529411764705882, 0.79607843137254897, 0.79607843137254897),
                              (0.23921568627450981, 0.80000000000000004, 0.80000000000000004),
                              (0.24313725490196078, 0.80784313725490198, 0.80784313725490198),
                              (0.24705882352941178, 0.81568627450980391, 0.81568627450980391),
                              (0.25098039215686274, 0.81960784313725488, 0.81960784313725488),
                              (0.25490196078431371, 0.82352941176470584, 0.82352941176470584),
                              (0.25882352941176473, 0.83137254901960789, 0.83137254901960789),
                              (0.2627450980392157, 0.83529411764705885, 0.83529411764705885),
                              (0.26666666666666666, 0.84313725490196079, 0.84313725490196079),
                              (0.27058823529411763, 0.84705882352941175, 0.84705882352941175),
                              (0.27450980392156865, 0.85098039215686272, 0.85098039215686272),
                              (0.27843137254901962, 0.85882352941176465, 0.85882352941176465),
                              (0.28235294117647058, 0.86274509803921573, 0.86274509803921573),
                              (0.28627450980392155, 0.8666666666666667, 0.8666666666666667),
                              (0.29019607843137257, 0.87058823529411766, 0.87058823529411766),
                              (0.29411764705882354, 0.61176470588235299, 0.61176470588235299),
                              (0.29803921568627451, 0.61176470588235299, 0.61176470588235299),
                              (0.30196078431372547, 0.61176470588235299, 0.61176470588235299),
                              (0.30588235294117649, 0.61176470588235299, 0.61176470588235299),
                              (0.30980392156862746, 0.61176470588235299, 0.61176470588235299),
                              (0.31372549019607843, 0.61176470588235299, 0.61176470588235299),
                              (0.31764705882352939, 0.61176470588235299, 0.61176470588235299),
                              (0.32156862745098042, 0.61176470588235299, 0.61176470588235299),
                              (0.32549019607843138, 0.61176470588235299, 0.61176470588235299),
                              (0.32941176470588235, 0.61176470588235299, 0.61176470588235299),
                              (0.33333333333333331, 0.61176470588235299, 0.61176470588235299),
                              (0.33725490196078434, 0.61176470588235299, 0.61176470588235299),
                              (0.3411764705882353, 0.61176470588235299, 0.61176470588235299),
                              (0.34509803921568627, 0.61176470588235299, 0.61176470588235299),
                              (0.34901960784313724, 0.61176470588235299, 0.61176470588235299),
                              (0.35294117647058826, 0.71372549019607845, 0.71372549019607845),
                              (0.35686274509803922, 0.71372549019607845, 0.71372549019607845),
                              (0.36078431372549019, 0.71372549019607845, 0.71372549019607845),
                              (0.36470588235294116, 0.71372549019607845, 0.71372549019607845),
                              (0.36862745098039218, 0.71372549019607845, 0.71372549019607845),
                              (0.37254901960784315, 0.71372549019607845, 0.71372549019607845),
                              (0.37647058823529411, 0.71372549019607845, 0.71372549019607845),
                              (0.38039215686274508, 0.71372549019607845, 0.71372549019607845),
                              (0.3843137254901961, 0.71372549019607845, 0.71372549019607845),
                              (0.38823529411764707, 0.71372549019607845, 0.71372549019607845),
                              (0.39215686274509803, 0.71372549019607845, 0.71372549019607845),
                              (0.396078431372549, 0.71372549019607845, 0.71372549019607845),
                              (0.40000000000000002, 0.71372549019607845, 0.71372549019607845),
                              (0.40392156862745099, 0.71372549019607845, 0.71372549019607845),
                              (0.40784313725490196, 0.71372549019607845, 0.71372549019607845),
                              (0.41176470588235292, 0.66666666666666663, 0.66666666666666663),
                              (0.41568627450980394, 0.66666666666666663, 0.66666666666666663),
                              (0.41960784313725491, 0.66666666666666663, 0.66666666666666663),
                              (0.42352941176470588, 0.66666666666666663, 0.66666666666666663),
                              (0.42745098039215684, 0.66666666666666663, 0.66666666666666663),
                              (0.43137254901960786, 0.66666666666666663, 0.66666666666666663),
                              (0.43529411764705883, 0.66666666666666663, 0.66666666666666663),
                              (0.4392156862745098, 0.66666666666666663, 0.66666666666666663),
                              (0.44313725490196076, 0.66666666666666663, 0.66666666666666663),
                              (0.44705882352941179, 0.66666666666666663, 0.66666666666666663),
                              (0.45098039215686275, 0.66666666666666663, 0.66666666666666663),
                              (0.45490196078431372, 0.66666666666666663, 0.66666666666666663),
                              (0.45882352941176469, 0.66666666666666663, 0.66666666666666663),
                              (0.46274509803921571, 0.66666666666666663, 0.66666666666666663),
                              (0.46666666666666667, 0.66666666666666663, 0.66666666666666663),
                              (0.47058823529411764, 0.63137254901960782, 0.63137254901960782),
                              (0.47450980392156861, 0.63137254901960782, 0.63137254901960782),
                              (0.47843137254901963, 0.63137254901960782, 0.63137254901960782),
                              (0.4823529411764706, 0.63137254901960782, 0.63137254901960782),
                              (0.48627450980392156, 0.63137254901960782, 0.63137254901960782),
                              (0.49019607843137253, 0.63137254901960782, 0.63137254901960782),
                              (0.49411764705882355, 0.63137254901960782, 0.63137254901960782),
                              (0.49803921568627452, 0.63137254901960782, 0.63137254901960782),
                              (0.50196078431372548, 0.63137254901960782, 0.63137254901960782),
                              (0.50588235294117645, 0.63137254901960782, 0.63137254901960782),
                              (0.50980392156862742, 0.63137254901960782, 0.63137254901960782),
                              (0.51372549019607838, 0.63137254901960782, 0.63137254901960782),
                              (0.51764705882352946, 0.63137254901960782, 0.63137254901960782),
                              (0.52156862745098043, 0.63137254901960782, 0.63137254901960782),
                              (0.52549019607843139, 0.63137254901960782, 0.63137254901960782),
                              (0.52941176470588236, 0.58039215686274515, 0.58039215686274515),
                              (0.53333333333333333, 0.58039215686274515, 0.58039215686274515),
                              (0.53725490196078429, 0.58039215686274515, 0.58039215686274515),
                              (0.54117647058823526, 0.58039215686274515, 0.58039215686274515),
                              (0.54509803921568623, 0.58039215686274515, 0.58039215686274515),
                              (0.5490196078431373, 0.58039215686274515, 0.58039215686274515),
                              (0.55294117647058827, 0.58039215686274515, 0.58039215686274515),
                              (0.55686274509803924, 0.58039215686274515, 0.58039215686274515),
                              (0.5607843137254902, 0.58039215686274515, 0.58039215686274515),
                              (0.56470588235294117, 0.58039215686274515, 0.58039215686274515),
                              (0.56862745098039214, 0.58039215686274515, 0.58039215686274515),
                              (0.5725490196078431, 0.58039215686274515, 0.58039215686274515),
                              (0.57647058823529407, 0.58039215686274515, 0.58039215686274515),
                              (0.58039215686274515, 0.58039215686274515, 0.58039215686274515),
                              (0.58431372549019611, 0.58039215686274515, 0.58039215686274515),
                              (0.58823529411764708, 0.52549019607843139, 0.52549019607843139),
                              (0.59215686274509804, 0.52549019607843139, 0.52549019607843139),
                              (0.59607843137254901, 0.52549019607843139, 0.52549019607843139),
                              (0.59999999999999998, 0.52549019607843139, 0.52549019607843139),
                              (0.60392156862745094, 0.52549019607843139, 0.52549019607843139),
                              (0.60784313725490191, 0.52549019607843139, 0.52549019607843139),
                              (0.61176470588235299, 0.52549019607843139, 0.52549019607843139),
                              (0.61568627450980395, 0.52549019607843139, 0.52549019607843139),
                              (0.61960784313725492, 0.52549019607843139, 0.52549019607843139),
                              (0.62352941176470589, 0.52549019607843139, 0.52549019607843139),
                              (0.62745098039215685, 0.52549019607843139, 0.52549019607843139),
                              (0.63137254901960782, 0.52549019607843139, 0.52549019607843139),
                              (0.63529411764705879, 0.52549019607843139, 0.52549019607843139),
                              (0.63921568627450975, 0.52549019607843139, 0.52549019607843139),
                              (0.64313725490196083, 0.52549019607843139, 0.52549019607843139),
                              (0.6470588235294118, 0.45098039215686275, 0.45098039215686275),
                              (0.65098039215686276, 0.45098039215686275, 0.45098039215686275),
                              (0.65490196078431373, 0.45098039215686275, 0.45098039215686275),
                              (0.6588235294117647, 0.45098039215686275, 0.45098039215686275),
                              (0.66274509803921566, 0.45098039215686275, 0.45098039215686275),
                              (0.66666666666666663, 0.45098039215686275, 0.45098039215686275),
                              (0.6705882352941176, 0.45098039215686275, 0.45098039215686275),
                              (0.67450980392156867, 0.45098039215686275, 0.45098039215686275),
                              (0.67843137254901964, 0.45098039215686275, 0.45098039215686275),
                              (0.68235294117647061, 0.45098039215686275, 0.45098039215686275),
                              (0.68627450980392157, 0.45098039215686275, 0.45098039215686275),
                              (0.69019607843137254, 0.45098039215686275, 0.45098039215686275),
                              (0.69411764705882351, 0.45098039215686275, 0.45098039215686275),
                              (0.69803921568627447, 0.45098039215686275, 0.45098039215686275),
                              (0.70196078431372544, 0.45098039215686275, 0.45098039215686275),
                              (0.70588235294117652, 0.37254901960784315, 0.37254901960784315),
                              (0.70980392156862748, 0.37254901960784315, 0.37254901960784315),
                              (0.71372549019607845, 0.37254901960784315, 0.37254901960784315),
                              (0.71764705882352942, 0.37254901960784315, 0.37254901960784315),
                              (0.72156862745098038, 0.37254901960784315, 0.37254901960784315),
                              (0.72549019607843135, 0.37254901960784315, 0.37254901960784315),
                              (0.72941176470588232, 0.37254901960784315, 0.37254901960784315),
                              (0.73333333333333328, 0.37254901960784315, 0.37254901960784315),
                              (0.73725490196078436, 0.37254901960784315, 0.37254901960784315),
                              (0.74117647058823533, 0.37254901960784315, 0.37254901960784315),
                              (0.74509803921568629, 0.37254901960784315, 0.37254901960784315),
                              (0.74901960784313726, 0.37254901960784315, 0.37254901960784315),
                              (0.75294117647058822, 0.37254901960784315, 0.37254901960784315),
                              (0.75686274509803919, 0.37254901960784315, 0.37254901960784315),
                              (0.76078431372549016, 0.37254901960784315, 0.37254901960784315),
                              (0.76470588235294112, 0.28235294117647058, 0.28235294117647058),
                              (0.7686274509803922, 0.28235294117647058, 0.28235294117647058),
                              (0.77254901960784317, 0.28235294117647058, 0.28235294117647058),
                              (0.77647058823529413, 0.28235294117647058, 0.28235294117647058),
                              (0.7803921568627451, 0.28235294117647058, 0.28235294117647058),
                              (0.78431372549019607, 0.28235294117647058, 0.28235294117647058),
                              (0.78823529411764703, 0.28235294117647058, 0.28235294117647058),
                              (0.792156862745098, 0.28235294117647058, 0.28235294117647058),
                              (0.79607843137254897, 0.28235294117647058, 0.28235294117647058),
                              (0.80000000000000004, 0.28235294117647058, 0.28235294117647058),
                              (0.80392156862745101, 0.28235294117647058, 0.28235294117647058),
                              (0.80784313725490198, 0.28235294117647058, 0.28235294117647058),
                              (0.81176470588235294, 0.28235294117647058, 0.28235294117647058),
                              (0.81568627450980391, 0.28235294117647058, 0.28235294117647058),
                              (0.81960784313725488, 0.28235294117647058, 0.28235294117647058),
                              (0.82352941176470584, 0.21568627450980393, 0.21568627450980393),
                              (0.82745098039215681, 0.21568627450980393, 0.21568627450980393),
                              (0.83137254901960789, 0.21568627450980393, 0.21568627450980393),
                              (0.83529411764705885, 0.21568627450980393, 0.21568627450980393),
                              (0.83921568627450982, 0.21568627450980393, 0.21568627450980393),
                              (0.84313725490196079, 0.21568627450980393, 0.21568627450980393),
                              (0.84705882352941175, 0.21568627450980393, 0.21568627450980393),
                              (0.85098039215686272, 0.21568627450980393, 0.21568627450980393),
                              (0.85490196078431369, 0.21568627450980393, 0.21568627450980393),
                              (0.85882352941176465, 0.21568627450980393, 0.21568627450980393),
                              (0.86274509803921573, 0.21568627450980393, 0.21568627450980393),
                              (0.8666666666666667, 0.21568627450980393, 0.21568627450980393),
                              (0.87058823529411766, 0.21568627450980393, 0.21568627450980393),
                              (0.87450980392156863, 0.21568627450980393, 0.21568627450980393),
                              (0.8784313725490196, 0.21568627450980393, 0.21568627450980393),
                              (0.88235294117647056, 0.16078431372549021, 0.16078431372549021),
                              (0.88627450980392153, 0.16078431372549021, 0.16078431372549021),
                              (0.8901960784313725, 0.16078431372549021, 0.16078431372549021),
                              (0.89411764705882357, 0.16078431372549021, 0.16078431372549021),
                              (0.89803921568627454, 0.16078431372549021, 0.16078431372549021),
                              (0.90196078431372551, 0.16078431372549021, 0.16078431372549021),
                              (0.90588235294117647, 0.16078431372549021, 0.16078431372549021),
                              (0.90980392156862744, 0.16078431372549021, 0.16078431372549021),
                              (0.9137254901960784, 0.16078431372549021, 0.16078431372549021),
                              (0.91764705882352937, 0.16078431372549021, 0.16078431372549021),
                              (0.92156862745098034, 0.16078431372549021, 0.16078431372549021),
                              (0.92549019607843142, 0.16078431372549021, 0.16078431372549021),
                              (0.92941176470588238, 0.16078431372549021, 0.16078431372549021),
                              (0.93333333333333335, 0.16078431372549021, 0.16078431372549021),
                              (0.93725490196078431, 0.16078431372549021, 0.16078431372549021),
                              (0.94117647058823528, 0.074509803921568626, 0.074509803921568626),
                              (0.94509803921568625, 0.074509803921568626, 0.074509803921568626),
                              (0.94901960784313721, 0.074509803921568626, 0.074509803921568626),
                              (0.95294117647058818, 0.074509803921568626, 0.074509803921568626),
                              (0.95686274509803926, 0.074509803921568626, 0.074509803921568626),
                              (0.96078431372549022, 0.074509803921568626, 0.074509803921568626),
                              (0.96470588235294119, 0.074509803921568626, 0.074509803921568626),
                              (0.96862745098039216, 0.074509803921568626, 0.074509803921568626),
                              (0.97254901960784312, 0.074509803921568626, 0.074509803921568626),
                              (0.97647058823529409, 0.074509803921568626, 0.074509803921568626),
                              (0.98039215686274506, 0.074509803921568626, 0.074509803921568626),
                              (0.98431372549019602, 0.074509803921568626, 0.074509803921568626),
                              (0.9882352941176471, 0.074509803921568626, 0.074509803921568626),
                              (0.99215686274509807, 0.074509803921568626, 0.074509803921568626),
                              (0.99607843137254903, 0.074509803921568626, 0.074509803921568626),
                              (1.0, 0.074509803921568626, 0.074509803921568626)],

                    'red': [  (0.0, 0.92156862745098034, 0.92156862745098034),
                              (0.0039215686274509803, 0.92156862745098034, 0.92156862745098034),
                              (0.0078431372549019607, 0.92156862745098034, 0.92156862745098034),
                              (0.011764705882352941, 0.92156862745098034, 0.92156862745098034),
                              (0.015686274509803921, 0.92156862745098034, 0.92156862745098034),
                              (0.019607843137254902, 0.92156862745098034, 0.92156862745098034),
                              (0.023529411764705882, 0.92156862745098034, 0.92156862745098034),
                              (0.027450980392156862, 0.92156862745098034, 0.92156862745098034),
                              (0.031372549019607843, 0.92156862745098034, 0.92156862745098034),
                              (0.035294117647058823, 0.92156862745098034, 0.92156862745098034),
                              (0.039215686274509803, 0.92156862745098034, 0.92156862745098034),
                              (0.043137254901960784, 0.92156862745098034, 0.92156862745098034),
                              (0.047058823529411764, 0.92156862745098034, 0.92156862745098034),
                              (0.050980392156862744, 0.92156862745098034, 0.92156862745098034),
                              (0.054901960784313725, 0.92156862745098034, 0.92156862745098034),
                              (0.058823529411764705, 0.92156862745098034, 0.92156862745098034),
                              (0.062745098039215685, 0.792156862745098, 0.792156862745098),
                              (0.066666666666666666, 0.79607843137254897, 0.79607843137254897),
                              (0.070588235294117646, 0.80000000000000004, 0.80000000000000004),
                              (0.074509803921568626, 0.80784313725490198, 0.80784313725490198),
                              (0.078431372549019607, 0.81176470588235294, 0.81176470588235294),
                              (0.082352941176470587, 0.81960784313725488, 0.81960784313725488),
                              (0.086274509803921567, 0.82352941176470584, 0.82352941176470584),
                              (0.090196078431372548, 0.82745098039215681, 0.82745098039215681),
                              (0.094117647058823528, 0.83529411764705885, 0.83529411764705885),
                              (0.098039215686274508, 0.83921568627450982, 0.83921568627450982),
                              (0.10196078431372549, 0.84313725490196079, 0.84313725490196079),
                              (0.10588235294117647, 0.85098039215686272, 0.85098039215686272),
                              (0.10980392156862745, 0.85490196078431369, 0.85490196078431369),
                              (0.11372549019607843, 0.85882352941176465, 0.85882352941176465),
                              (0.11764705882352941, 0.8666666666666667, 0.8666666666666667),
                              (0.12156862745098039, 0.87058823529411766, 0.87058823529411766),
                              (0.12549019607843137, 0.87450980392156863, 0.87450980392156863),
                              (0.12941176470588237, 0.8784313725490196, 0.8784313725490196),
                              (0.13333333333333333, 0.88627450980392153, 0.88627450980392153),
                              (0.13725490196078433, 0.8901960784313725, 0.8901960784313725),
                              (0.14117647058823529, 0.89411764705882357, 0.89411764705882357),
                              (0.14509803921568629, 0.89803921568627454, 0.89803921568627454),
                              (0.14901960784313725, 0.90196078431372551, 0.90196078431372551),
                              (0.15294117647058825, 0.90588235294117647, 0.90588235294117647),
                              (0.15686274509803921, 0.9137254901960784, 0.9137254901960784),
                              (0.16078431372549021, 0.91764705882352937, 0.91764705882352937),
                              (0.16470588235294117, 0.92156862745098034, 0.92156862745098034),
                              (0.16862745098039217, 0.92549019607843142, 0.92549019607843142),
                              (0.17254901960784313, 0.92941176470588238, 0.92941176470588238),
                              (0.17647058823529413, 0.93333333333333335, 0.93333333333333335),
                              (0.1803921568627451, 0.93725490196078431, 0.93725490196078431),
                              (0.18431372549019609, 0.94117647058823528, 0.94117647058823528),
                              (0.18823529411764706, 0.94509803921568625, 0.94509803921568625),
                              (0.19215686274509805, 0.94901960784313721, 0.94901960784313721),
                              (0.19607843137254902, 0.94901960784313721, 0.94901960784313721),
                              (0.20000000000000001, 0.95294117647058818, 0.95294117647058818),
                              (0.20392156862745098, 0.95686274509803926, 0.95686274509803926),
                              (0.20784313725490197, 0.96078431372549022, 0.96078431372549022),
                              (0.21176470588235294, 0.96470588235294119, 0.96470588235294119),
                              (0.21568627450980393, 0.96470588235294119, 0.96470588235294119),
                              (0.2196078431372549, 0.96862745098039216, 0.96862745098039216),
                              (0.22352941176470589, 0.97254901960784312, 0.97254901960784312),
                              (0.22745098039215686, 0.97254901960784312, 0.97254901960784312),
                              (0.23137254901960785, 0.97647058823529409, 0.97647058823529409),
                              (0.23529411764705882, 0.98039215686274506, 0.98039215686274506),
                              (0.23921568627450981, 0.98039215686274506, 0.98039215686274506),
                              (0.24313725490196078, 0.98431372549019602, 0.98431372549019602),
                              (0.24705882352941178, 0.98431372549019602, 0.98431372549019602),
                              (0.25098039215686274, 0.9882352941176471, 0.9882352941176471),
                              (0.25490196078431371, 0.9882352941176471, 0.9882352941176471),
                              (0.25882352941176473, 0.99215686274509807, 0.99215686274509807),
                              (0.2627450980392157, 0.99215686274509807, 0.99215686274509807),
                              (0.26666666666666666, 0.99215686274509807, 0.99215686274509807),
                              (0.27058823529411763, 0.99607843137254903, 0.99607843137254903),
                              (0.27450980392156865, 0.99607843137254903, 0.99607843137254903),
                              (0.27843137254901962, 0.99607843137254903, 0.99607843137254903),
                              (0.28235294117647058, 1.0, 1.0),
                              (0.28627450980392155, 1.0, 1.0),
                              (0.29019607843137257, 1.0, 1.0),
                              (0.29411764705882354, 0.49411764705882355, 0.49411764705882355),
                              (0.29803921568627451, 0.49411764705882355, 0.49411764705882355),
                              (0.30196078431372547, 0.49411764705882355, 0.49411764705882355),
                              (0.30588235294117649, 0.49411764705882355, 0.49411764705882355),
                              (0.30980392156862746, 0.49411764705882355, 0.49411764705882355),
                              (0.31372549019607843, 0.49411764705882355, 0.49411764705882355),
                              (0.31764705882352939, 0.49411764705882355, 0.49411764705882355),
                              (0.32156862745098042, 0.49411764705882355, 0.49411764705882355),
                              (0.32549019607843138, 0.49411764705882355, 0.49411764705882355),
                              (0.32941176470588235, 0.49411764705882355, 0.49411764705882355),
                              (0.33333333333333331, 0.49411764705882355, 0.49411764705882355),
                              (0.33725490196078434, 0.49411764705882355, 0.49411764705882355),
                              (0.3411764705882353, 0.49411764705882355, 0.49411764705882355),
                              (0.34509803921568627, 0.49411764705882355, 0.49411764705882355),
                              (0.34901960784313724, 0.49411764705882355, 0.49411764705882355),
                              (0.35294117647058826, 0.58823529411764708, 0.58823529411764708),
                              (0.35686274509803922, 0.58823529411764708, 0.58823529411764708),
                              (0.36078431372549019, 0.58823529411764708, 0.58823529411764708),
                              (0.36470588235294116, 0.58823529411764708, 0.58823529411764708),
                              (0.36862745098039218, 0.58823529411764708, 0.58823529411764708),
                              (0.37254901960784315, 0.58823529411764708, 0.58823529411764708),
                              (0.37647058823529411, 0.58823529411764708, 0.58823529411764708),
                              (0.38039215686274508, 0.58823529411764708, 0.58823529411764708),
                              (0.3843137254901961, 0.58823529411764708, 0.58823529411764708),
                              (0.38823529411764707, 0.58823529411764708, 0.58823529411764708),
                              (0.39215686274509803, 0.58823529411764708, 0.58823529411764708),
                              (0.396078431372549, 0.58823529411764708, 0.58823529411764708),
                              (0.40000000000000002, 0.58823529411764708, 0.58823529411764708),
                              (0.40392156862745099, 0.58823529411764708, 0.58823529411764708),
                              (0.40784313725490196, 0.58823529411764708, 0.58823529411764708),
                              (0.41176470588235292, 0.45882352941176469, 0.45882352941176469),
                              (0.41568627450980394, 0.45882352941176469, 0.45882352941176469),
                              (0.41960784313725491, 0.45882352941176469, 0.45882352941176469),
                              (0.42352941176470588, 0.45882352941176469, 0.45882352941176469),
                              (0.42745098039215684, 0.45882352941176469, 0.45882352941176469),
                              (0.43137254901960786, 0.45882352941176469, 0.45882352941176469),
                              (0.43529411764705883, 0.45882352941176469, 0.45882352941176469),
                              (0.4392156862745098, 0.45882352941176469, 0.45882352941176469),
                              (0.44313725490196076, 0.45882352941176469, 0.45882352941176469),
                              (0.44705882352941179, 0.45882352941176469, 0.45882352941176469),
                              (0.45098039215686275, 0.45882352941176469, 0.45882352941176469),
                              (0.45490196078431372, 0.45882352941176469, 0.45882352941176469),
                              (0.45882352941176469, 0.45882352941176469, 0.45882352941176469),
                              (0.46274509803921571, 0.45882352941176469, 0.45882352941176469),
                              (0.46666666666666667, 0.45882352941176469, 0.45882352941176469),
                              (0.47058823529411764, 0.40392156862745099, 0.40392156862745099),
                              (0.47450980392156861, 0.40392156862745099, 0.40392156862745099),
                              (0.47843137254901963, 0.40392156862745099, 0.40392156862745099),
                              (0.4823529411764706, 0.40392156862745099, 0.40392156862745099),
                              (0.48627450980392156, 0.40392156862745099, 0.40392156862745099),
                              (0.49019607843137253, 0.40392156862745099, 0.40392156862745099),
                              (0.49411764705882355, 0.40392156862745099, 0.40392156862745099),
                              (0.49803921568627452, 0.40392156862745099, 0.40392156862745099),
                              (0.50196078431372548, 0.40392156862745099, 0.40392156862745099),
                              (0.50588235294117645, 0.40392156862745099, 0.40392156862745099),
                              (0.50980392156862742, 0.40392156862745099, 0.40392156862745099),
                              (0.51372549019607838, 0.40392156862745099, 0.40392156862745099),
                              (0.51764705882352946, 0.40392156862745099, 0.40392156862745099),
                              (0.52156862745098043, 0.40392156862745099, 0.40392156862745099),
                              (0.52549019607843139, 0.40392156862745099, 0.40392156862745099),
                              (0.52941176470588236, 0.32156862745098042, 0.32156862745098042),
                              (0.53333333333333333, 0.32156862745098042, 0.32156862745098042),
                              (0.53725490196078429, 0.32156862745098042, 0.32156862745098042),
                              (0.54117647058823526, 0.32156862745098042, 0.32156862745098042),
                              (0.54509803921568623, 0.32156862745098042, 0.32156862745098042),
                              (0.5490196078431373, 0.32156862745098042, 0.32156862745098042),
                              (0.55294117647058827, 0.32156862745098042, 0.32156862745098042),
                              (0.55686274509803924, 0.32156862745098042, 0.32156862745098042),
                              (0.5607843137254902, 0.32156862745098042, 0.32156862745098042),
                              (0.56470588235294117, 0.32156862745098042, 0.32156862745098042),
                              (0.56862745098039214, 0.32156862745098042, 0.32156862745098042),
                              (0.5725490196078431, 0.32156862745098042, 0.32156862745098042),
                              (0.57647058823529407, 0.32156862745098042, 0.32156862745098042),
                              (0.58039215686274515, 0.32156862745098042, 0.32156862745098042),
                              (0.58431372549019611, 0.32156862745098042, 0.32156862745098042),
                              (0.58823529411764708, 0.23921568627450981, 0.23921568627450981),
                              (0.59215686274509804, 0.23921568627450981, 0.23921568627450981),
                              (0.59607843137254901, 0.23921568627450981, 0.23921568627450981),
                              (0.59999999999999998, 0.23921568627450981, 0.23921568627450981),
                              (0.60392156862745094, 0.23921568627450981, 0.23921568627450981),
                              (0.60784313725490191, 0.23921568627450981, 0.23921568627450981),
                              (0.61176470588235299, 0.23921568627450981, 0.23921568627450981),
                              (0.61568627450980395, 0.23921568627450981, 0.23921568627450981),
                              (0.61960784313725492, 0.23921568627450981, 0.23921568627450981),
                              (0.62352941176470589, 0.23921568627450981, 0.23921568627450981),
                              (0.62745098039215685, 0.23921568627450981, 0.23921568627450981),
                              (0.63137254901960782, 0.23921568627450981, 0.23921568627450981),
                              (0.63529411764705879, 0.23921568627450981, 0.23921568627450981),
                              (0.63921568627450975, 0.23921568627450981, 0.23921568627450981),
                              (0.64313725490196083, 0.23921568627450981, 0.23921568627450981),
                              (0.6470588235294118, 0.10980392156862745, 0.10980392156862745),
                              (0.65098039215686276, 0.10980392156862745, 0.10980392156862745),
                              (0.65490196078431373, 0.10980392156862745, 0.10980392156862745),
                              (0.6588235294117647, 0.10980392156862745, 0.10980392156862745),
                              (0.66274509803921566, 0.10980392156862745, 0.10980392156862745),
                              (0.66666666666666663, 0.10980392156862745, 0.10980392156862745),
                              (0.6705882352941176, 0.10980392156862745, 0.10980392156862745),
                              (0.67450980392156867, 0.10980392156862745, 0.10980392156862745),
                              (0.67843137254901964, 0.10980392156862745, 0.10980392156862745),
                              (0.68235294117647061, 0.10980392156862745, 0.10980392156862745),
                              (0.68627450980392157, 0.10980392156862745, 0.10980392156862745),
                              (0.69019607843137254, 0.10980392156862745, 0.10980392156862745),
                              (0.69411764705882351, 0.10980392156862745, 0.10980392156862745),
                              (0.69803921568627447, 0.10980392156862745, 0.10980392156862745),
                              (0.70196078431372544, 0.10980392156862745, 0.10980392156862745),
                              (0.70588235294117652, 0.0, 0.0),
                              (0.70980392156862748, 0.0, 0.0),
                              (0.71372549019607845, 0.0, 0.0),
                              (0.71764705882352942, 0.0, 0.0),
                              (0.72156862745098038, 0.0, 0.0),
                              (0.72549019607843135, 0.0, 0.0),
                              (0.72941176470588232, 0.0, 0.0),
                              (0.73333333333333328, 0.0, 0.0),
                              (0.73725490196078436, 0.0, 0.0),
                              (0.74117647058823533, 0.0, 0.0),
                              (0.74509803921568629, 0.0, 0.0),
                              (0.74901960784313726, 0.0, 0.0),
                              (0.75294117647058822, 0.0, 0.0),
                              (0.75686274509803919, 0.0, 0.0),
                              (0.76078431372549016, 0.0, 0.0),
                              (0.76470588235294112, 0.0, 0.0),
                              (0.7686274509803922, 0.0, 0.0),
                              (0.77254901960784317, 0.0, 0.0),
                              (0.77647058823529413, 0.0, 0.0),
                              (0.7803921568627451, 0.0, 0.0),
                              (0.78431372549019607, 0.0, 0.0),
                              (0.78823529411764703, 0.0, 0.0),
                              (0.792156862745098, 0.0, 0.0),
                              (0.79607843137254897, 0.0, 0.0),
                              (0.80000000000000004, 0.0, 0.0),
                              (0.80392156862745101, 0.0, 0.0),
                              (0.80784313725490198, 0.0, 0.0),
                              (0.81176470588235294, 0.0, 0.0),
                              (0.81568627450980391, 0.0, 0.0),
                              (0.81960784313725488, 0.0, 0.0),
                              (0.82352941176470584, 0.0078431372549019607, 0.0078431372549019607),
                              (0.82745098039215681, 0.0078431372549019607, 0.0078431372549019607),
                              (0.83137254901960789, 0.0078431372549019607, 0.0078431372549019607),
                              (0.83529411764705885, 0.0078431372549019607, 0.0078431372549019607),
                              (0.83921568627450982, 0.0078431372549019607, 0.0078431372549019607),
                              (0.84313725490196079, 0.0078431372549019607, 0.0078431372549019607),
                              (0.84705882352941175, 0.0078431372549019607, 0.0078431372549019607),
                              (0.85098039215686272, 0.0078431372549019607, 0.0078431372549019607),
                              (0.85490196078431369, 0.0078431372549019607, 0.0078431372549019607),
                              (0.85882352941176465, 0.0078431372549019607, 0.0078431372549019607),
                              (0.86274509803921573, 0.0078431372549019607, 0.0078431372549019607),
                              (0.8666666666666667, 0.0078431372549019607, 0.0078431372549019607),
                              (0.87058823529411766, 0.0078431372549019607, 0.0078431372549019607),
                              (0.87450980392156863, 0.0078431372549019607, 0.0078431372549019607),
                              (0.8784313725490196, 0.0078431372549019607, 0.0078431372549019607),
                              (0.88235294117647056, 0.0, 0.0),
                              (0.88627450980392153, 0.0, 0.0),
                              (0.8901960784313725, 0.0, 0.0),
                              (0.89411764705882357, 0.0, 0.0),
                              (0.89803921568627454, 0.0, 0.0),
                              (0.90196078431372551, 0.0, 0.0),
                              (0.90588235294117647, 0.0, 0.0),
                              (0.90980392156862744, 0.0, 0.0),
                              (0.9137254901960784, 0.0, 0.0),
                              (0.91764705882352937, 0.0, 0.0),
                              (0.92156862745098034, 0.0, 0.0),
                              (0.92549019607843142, 0.0, 0.0),
                              (0.92941176470588238, 0.0, 0.0),
                              (0.93333333333333335, 0.0, 0.0),
                              (0.93725490196078431, 0.0, 0.0),
                              (0.94117647058823528, 0.0039215686274509803, 0.0039215686274509803),
                              (0.94509803921568625, 0.0039215686274509803, 0.0039215686274509803),
                              (0.94901960784313721, 0.0039215686274509803, 0.0039215686274509803),
                              (0.95294117647058818, 0.0039215686274509803, 0.0039215686274509803),
                              (0.95686274509803926, 0.0039215686274509803, 0.0039215686274509803),
                              (0.96078431372549022, 0.0039215686274509803, 0.0039215686274509803),
                              (0.96470588235294119, 0.0039215686274509803, 0.0039215686274509803),
                              (0.96862745098039216, 0.0039215686274509803, 0.0039215686274509803),
                              (0.97254901960784312, 0.0039215686274509803, 0.0039215686274509803),
                              (0.97647058823529409, 0.0039215686274509803, 0.0039215686274509803),
                              (0.98039215686274506, 0.0039215686274509803, 0.0039215686274509803),
                              (0.98431372549019602, 0.0039215686274509803, 0.0039215686274509803),
                              (0.9882352941176471, 0.0039215686274509803, 0.0039215686274509803),
                              (0.99215686274509807, 0.0039215686274509803, 0.0039215686274509803),
                              (0.99607843137254903, 0.0039215686274509803, 0.0039215686274509803),
                              (1.0, 0.0039215686274509803, 0.0039215686274509803)]
                              }

            self.cmap = LinearSegmentedColormap('ndvi', ndvi_cmdata)


    #SeaSurfaceTempProduct = {
        #'bulk_sst': SeaSurfaceTempProd(
            #SDSname = '/All_Data/VIIRS-VI-EDR_All/BulkSST',
            #SDStitleString = "VIIRS Bulk Sea Surface Temperature",
            #vmin = 273.,
            #vmax = 300.,
            #cbarTitle = "SST ($K$)",
            ##cmap = cm.jet,
        #),
        #'skin_sst': SeaSurfaceTempProd(
            #SDSname = '/All_Data/VIIRS-VI-EDR_All/SkinSST',
            #SDStitleString = "VIIRS Skin Sea Surface Temperature",
            #vmin = 273.,
            #vmax = 300.,
            #cbarTitle = "SST ($K$)",
            ##cmap = cm.jet,
        #),
    #}

