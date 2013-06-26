#!/usr/bin/env python
# encoding: utf-8
"""
viirs_edr_data.py

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

