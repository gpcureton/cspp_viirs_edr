import numpy as np
from os import path,readlink
import string
from numpy import ma
from glob import glob
import re

import matplotlib
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap
from matplotlib.figure import Figure

matplotlib.use('Agg')
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

# This must come *after* the backend is specified.
import matplotlib.pyplot as ppl
from mpl_toolkits.basemap import Basemap,addcyclic,shiftgrid

import adl_blob
from ViirsData import ViirsTrimTable

# every module should have a LOG object
# e.g. LOG.warning('my dog has fleas')
import logging
logging.basicConfig() 

#dpi=200

### Moderate and Imager resolution trim table arrays. These are 
### bool arrays, and the trim pixels are set to True.
trimObj = ViirsTrimTable()
modTrimMask = trimObj.createModTrimArray(nscans=48,trimType=bool)

#xmlDir = '/home/geoffc/SSEC/sshfs/rh5b/RH5B/CSPP/EDR/common/ADL/xml/VIIRS'
#oldBlobPath = '/home/geoffc/SSEC/sshfs/rh5b/WORK/sample_data/viirs/edr/input/VIIRS_DB_ADL4.1/blobs/single/ANC'
#newBlobPath = '/home/geoffc/SSEC/sshfs/rh5b/WORK/Work_Area/VIIRS_DB_ADL4.1'


def get_blob_asc_dict(xmlDir,blobDir,collShortNames):
    blob_asc_Files = {}
    
    blobDir = path.expanduser(blobDir)

    if path.isdir(blobDir):
        pass
    else :
        blobDir = path.dirname(blobDir)

    for shortName in collShortNames :
        #print '%s --> ' % (shortName)
        fileGlob = path.join(path.expanduser(blobDir),'*.%s'%(shortName))

        blobFiles = glob(fileGlob)
        if blobFiles != []:
            blobFiles.sort()
            blobDict = {}
            print '%s --> ' % (shortName)
            for files in blobFiles :
                blobFile = path.basename(files)
                ascFile = string.replace(blobFile,shortName,'asc')
                ascFullPath = path.join(blobDir,'%s'%(ascFile))
                ascFileObj = open(ascFullPath,"r")
                for line in ascFileObj:
                    if re.search("N_Granule_ID", line):
                        granID = string.split(line,'"')[3]
                ascFileObj.close()
                blobDict[granID] = [blobFile,ascFile]
            blob_asc_Files[shortName] = blobDict

    return blob_asc_Files

def NISE_globalPlot(NISE_fileName,pngDir=None):

    from HDF4File import HDF4File

    if pngDir is None :
        pngDir = path.abspath(path.curdir)

    latStart,latEnd = 2*3000,2*4500
    lonStart,lonEnd = 2*3000,2*5000
    print "NISE_latMinIdx = %d" % (latEnd)
    print "NISE_latMaxIdx = %d" % (latStart)
    print "NISE_lonMinIdx = %d" % (lonStart)
    print "NISE_lonMaxIdx = %d" % (lonEnd)

    try :
        fileObj = HDF4File(NISE_fileName)
    except Exception, err :
        print "%s" % (err)
        print "Problem opening NISE file (%s), aborting." % (NISE_fileName)
        sys.exit(1)

    try :

        northDsetName = "Northern Hemisphere/Data Fields/Extent"
        southDsetName = "Southern Hemisphere/Data Fields/Extent"

        print "Retrieving NISE HDF4 path '%s'" % (northDsetName)
        nHemi = fileObj.get_dataset(northDsetName)

        print "Retrieving NISE HDF4 path '%s'" % (southDsetName)
        sHemi = fileObj.get_dataset(southDsetName)

    except Exception, err :

        print "EXCEPTION: %s" % (err)
        sys.exit(1)

    for data,title in zip([nHemi,sHemi],['Northern','Sourthern']):

        plotTitle =  "NISE %s Hemi : %s" %(title,path.basename(NISE_fileName))
        cbTitle   =  "Snow / Ice"
        #vmin,vmax =  0,1
        vmin,vmax =  None,None

        # Create figure with default size, and create canvas to draw on
        scale=1.5
        fig = Figure(figsize=(scale*8,scale*8))
        canvas = FigureCanvas(fig)

        # Create main axes instance, leaving room for colorbar at bottom,
        # and also get the Bbox of the axes instance
        ax_rect = [0.05, 0.18, 0.9, 0.75  ] # [left,bottom,width,height]
        ax = fig.add_axes(ax_rect)

        # Granule axis title
        ax_title = ppl.setp(ax,title=plotTitle)
        ppl.setp(ax_title,fontsize=12)
        ppl.setp(ax_title,family="sans-serif")

        # Plot the data
        im = ax.imshow(data,axes=ax,interpolation='nearest',vmin=vmin,vmax=vmax)
        #im = ax.imshow(data,axes=ax,interpolation='nearest')
        
        # add a colorbar axis
        cax_rect = [0.05 , 0.05, 0.9 , 0.06 ] # [left,bottom,width,height]
        cax = fig.add_axes(cax_rect,frameon=False) # setup colorbar axes

        # Plot the colorbar.
        cb = fig.colorbar(im, cax=cax, orientation='horizontal')
        ppl.setp(cax.get_xticklabels(),fontsize=9)

        # Colourbar title
        cax_title = ppl.setp(cax,title=cbTitle)
        ppl.setp(cax_title,fontsize=9)

        # Redraw the figure
        canvas.draw()

        # save image 
        pngFile = "%s/%s_%s.png" % (pngDir,path.basename(NISE_fileName),title)
        print "Writing file to ",pngFile
        canvas.print_figure(pngFile,dpi=200)




class NCEPclass():

    def __init__(self,xmlFile,blobFile):
        xmlFile = path.expanduser(xmlFile)
        blobFile = path.expanduser(blobFile)

        # A default 0.5 degree grid...
        degInc = 0.5
        grids = np.mgrid[-90.:90.+degInc:degInc,-180.:180.:degInc]
        self.gridLat,self.gridLon = grids[0],grids[1]

    pngDir = path.abspath(path.curdir)

    ancDataSetTitle = {}
    ancDataSetTitle['geopotentialHeightLayers'] = 'geopotentialHeightLayers'
    ancDataSetTitle['temperatureLayers'] = 'temperatureLayers'
    ancDataSetTitle['waterVaporMixingRatioLayers'] = 'waterVaporMixingRatioLayers'
    ancDataSetTitle['pressureReducedToMSL'] = 'pressureReducedToMSL'
    ancDataSetTitle['uComponentOfWind'] = 'uComponentOfWind'
    ancDataSetTitle['vComponentOfWind'] = 'vComponentOfWind'
    ancDataSetTitle['surfacePressure'] = 'surfacePressure'
    ancDataSetTitle['skinTemperature'] = 'skinTemperature'
    ancDataSetTitle['surfaceTemperature'] = 'surfaceTemperature'
    ancDataSetTitle['totalPrecipitableWater'] = 'totalPrecipitableWater'
    ancDataSetTitle['surfaceGeopotentialHeight'] = 'surfaceGeopotentialHeight'
    ancDataSetTitle['surfaceSpecificHumidity'] = 'surfaceSpecificHumidity'
    ancDataSetTitle['tropopauseGeopotentialHeight'] = 'tropopauseGeopotentialHeight'
    ancDataSetTitle['totalColumnOzone'] = 'totalColumnOzone'

    NCEP_LAYER_LEVELS = {
                          '10mb' : 0,    '20mb' : 1,   '30mb' : 2,   '50mb' : 3,   '70mb' : 4,  '100mb' : 5,
                         '150mb' : 6,   '200mb' : 7,  '250mb' : 8,  '300mb' : 9,  '350mb' : 10, '400mb' : 11, 
                         '450mb' : 12,  '500mb' : 13, '550mb' : 14, '600mb' : 15, '650mb' : 16, '700mb' : 17, 
                         '750mb' : 18,  '800mb' : 19, '850mb' : 20, '900mb' : 21, '925mb' : 22, '950mb' : 23, 
                         '975mb' : 24, '1000mb' : 25
                        }

    NCEP_LAYER_VALUES = np.array([10.0, 20.0, 30.0, 50.0, 70.0, 100.0, 150.0, 200.0, 
                                  250.0, 300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 
                                  600.0, 650.0, 700.0, 750.0, 800.0, 850.0, 900.0, 
                                  925.0, 950.0, 975.0, 1000.0 ]);


    @staticmethod
    def plotAncData(gridLat,gridLon,gridData,plotData) : #ancType,plotTitle,cbTitle) :

        ancType   =  plotData['ancType']
        plotTitle =  plotData['plotTitle']
        cbTitle   =  plotData['cbTitle']
        vmin,vmax =  plotData['plotLims'][0], plotData['plotLims'][1]

        # Create figure with default size, and create canvas to draw on
        scale=1.5
        fig = Figure(figsize=(scale*8,scale*5))
        canvas = FigureCanvas(fig)

        # Create main axes instance, leaving room for colorbar at bottom,
        # and also get the Bbox of the axes instance
        ax_rect = [0.05, 0.18, 0.9, 0.75  ] # [left,bottom,width,height]
        ax = fig.add_axes(ax_rect)

        # Granule axis title
        ax_title = ppl.setp(ax,title=plotTitle)
        ppl.setp(ax_title,fontsize=12)
        ppl.setp(ax_title,family="sans-serif")

        # Create the basemap object
        m = Basemap(projection='cyl',lon_0=0.,ax=ax)
        x,y = m(gridLon,gridLat)

        # Plot the data
        im = m.imshow(gridData.astype('float32'),axes=ax,interpolation='nearest',vmin=vmin,vmax=vmax)
        
        m.drawmapboundary(ax=ax,linewidth=0.01,fill_color='grey')
        m.drawcoastlines(ax=ax,linewidth=0.5)

        # draw parallels
        delat = 30.
        circles = np.arange(-90.,90.+delat,delat)
        m.drawparallels(circles,ax=ax,labelstyle="+/-",labels=[1,0,0,0])

        # draw meridians
        delon = 60.
        meridians = np.arange(-180,180,delon)
        m.drawmeridians(meridians,ax=ax,labelstyle="+/-",labels=[0,0,0,1])

        # add a colorbar axis
        cax_rect = [0.05 , 0.05, 0.9 , 0.06 ] # [left,bottom,width,height]
        cax = fig.add_axes(cax_rect,frameon=False) # setup colorbar axes

        # Plot the colorbar.
        cb = fig.colorbar(im, cax=cax, orientation='horizontal')
        ppl.setp(cax.get_xticklabels(),fontsize=9)

        # Colourbar title
        cax_title = ppl.setp(cax,title=cbTitle)
        ppl.setp(cax_title,fontsize=9)

        # Redraw the figure
        canvas.draw()

        # save image 
        levelStr = "_"+string.replace(string.split(plotTitle,'@')[1],')','') if "@" in plotTitle else ""
        levelStr = string.replace(levelStr," ","")
        print "LevelStr = ",levelStr
        pngFile = "%s/NCEP-ANC-Int_%s%s.png" % (NCEPclass.pngDir,NCEPclass.ancDataSetTitle[ancType],levelStr)
        print "Writing file to ",pngFile
        canvas.print_figure(pngFile,dpi=dpi)

        del(m)
        return 0


class ANCclass():

    def __init__(self,xmlDir,blobPath):

        self.xmlDir = xmlDir
        self.blobPath = path.expanduser(blobPath)

        self.collShortNames = [
                               'VIIRS-ANC-Preci-Wtr-Mod-Gran',
                               'VIIRS-ANC-Temp-Surf2M-Mod-Gran',
                               'VIIRS-ANC-Wind-Speed-Mod-Gran',
                               'VIIRS-ANC-Wind-Direction-Mod-Gran',
                               'VIIRS-ANC-Surf-Ht-Mod-Gran',
                               'VIIRS-ANC-Press-Surf-Mod-Gran',
                               'VIIRS-ANC-Tot-Col-Mod-Gran',
                               'VIIRS-ANC-Optical-Depth-Mod-Gran',
                               'VIIRS-ANC-Geopot-Ht-Lev-Mod-Gran',
                               'VIIRS-ANC-Sp-Humd-Surf-Mod-Gran',
                               'VIIRS-ANC-Temp-Skin-Mod-Gran'
                              ]


        self.xmlName = {}
        self.xmlName['VIIRS-ANC-Preci-Wtr-Mod-Gran'] = 'VIIRS_ANC_PRECI_WTR_MOD_GRAN.xml'
        self.xmlName['VIIRS-ANC-Temp-Surf2M-Mod-Gran'] = 'VIIRS_ANC_TEMP_SURF2M_MOD_GRAN.xml'
        self.xmlName['VIIRS-ANC-Wind-Speed-Mod-Gran'] = 'VIIRS_ANC_WIND_SPEED_MOD_GRAN.xml'
        self.xmlName['VIIRS-ANC-Wind-Direction-Mod-Gran'] = 'VIIRS_ANC_WIND_DIRECTION_MOD_GRAN.xml'
        self.xmlName['VIIRS-ANC-Surf-Ht-Mod-Gran'] = 'VIIRS_ANC_SURF_HT_MOD_GRAN.xml'
        self.xmlName['VIIRS-ANC-Press-Surf-Mod-Gran'] = 'VIIRS_ANC_PRESS_SURF_MOD_GRAN.xml'
        self.xmlName['VIIRS-ANC-Tot-Col-Mod-Gran'] = 'VIIRS_ANC_TOT_COL_MOD_GRAN.xml'
        self.xmlName['VIIRS-ANC-Optical-Depth-Mod-Gran'] = 'VIIRS_ANC_OPTICAL_DEPTH_MOD_GRAN.xml'
        self.xmlName['VIIRS-ANC-Geopot-Ht-Lev-Mod-Gran'] = 'VIIRS_ANC_GEOPOT_HT_LEV_MOD_GRAN.xml'
        self.xmlName['VIIRS-ANC-Sp-Humd-Surf-Mod-Gran'] = 'VIIRS_ANC_SP_HUMD_SURF_MOD_GRAN.xml'
        self.xmlName['VIIRS-ANC-Temp-Skin-Mod-Gran'] = 'VIIRS_ANC_TEMP_SKIN_MOD_GRAN.xml'

        self.plotDescr = {}
        self.plotDescr['VIIRS-ANC-Preci-Wtr-Mod-Gran']   = r'Total Precipitable Water (cm, $\times$10 kg m$^{-2}$)'
        #self.plotDescr['VIIRS-ANC-Preci-Wtr-Mod-Gran']   = r'Total Precipitable Water (kg m$^{-2}$)'
        self.plotDescr['VIIRS-ANC-Temp-Surf2M-Mod-Gran'] = r'Air Temperature @ 2m (K)'
        self.plotDescr['VIIRS-ANC-Wind-Speed-Mod-Gran']  = r'Surface Wind velocity (ms$^{-1}$)'
        self.plotDescr['VIIRS-ANC-Wind-Direction-Mod-Gran'] = r'Surface Wind Direction (degrees)'
        self.plotDescr['VIIRS-ANC-Surf-Ht-Mod-Gran']     = r'GMTCO Surface Height (m)'
        self.plotDescr['VIIRS-ANC-Press-Surf-Mod-Gran']  = r'Surface Pressure (mb, $\times$10$^{2}$ Pa)'
        self.plotDescr['VIIRS-ANC-Tot-Col-Mod-Gran']     = r'Total Column Ozone (Atm cm, $\times$10$^{3}$DU)'
        self.plotDescr['VIIRS-ANC-Optical-Depth-Mod-Gran'] = 'NAAPS Total Column Aerosol Optical Depth'
        self.plotDescr['VIIRS-ANC-Geopot-Ht-Lev-Mod-Gran'] = 'Surface Geopotential Height (gpm)'
        self.plotDescr['VIIRS-ANC-Sp-Humd-Surf-Mod-Gran'] = 'Surface Specific Humidity'
        self.plotDescr['VIIRS-ANC-Temp-Skin-Mod-Gran'] = 'Sea Surface Skin Temperature (K)'

        self.plotLims = {}
        self.plotLims['VIIRS-ANC-Preci-Wtr-Mod-Gran']    = [None,None] #[0.,0.7]
        #self.plotLims['VIIRS-ANC-Preci-Wtr-Mod-Gran']    = [0.,130.]
        self.plotLims['VIIRS-ANC-Temp-Surf2M-Mod-Gran']  = [253.,293.]
        self.plotLims['VIIRS-ANC-Wind-Speed-Mod-Gran']   = [0.,20.]
        self.plotLims['VIIRS-ANC-Wind-Direction-Mod-Gran'] =[0.,360.]
        self.plotLims['VIIRS-ANC-Surf-Ht-Mod-Gran']      = [None, None]
        #self.plotLims['VIIRS-ANC-Press-Surf-Mod-Gran']  = [300.,1080.]
        self.plotLims['VIIRS-ANC-Press-Surf-Mod-Gran']  = [800.,1013.]
        #self.plotLims['VIIRS-ANC-Press-Surf-Mod-Gran']  = [None, None]
        self.plotLims['VIIRS-ANC-Tot-Col-Mod-Gran']     = [None, None]
        self.plotLims['VIIRS-ANC-Optical-Depth-Mod-Gran'] = [None, None]
        self.plotLims['VIIRS-ANC-Geopot-Ht-Lev-Mod-Gran'] = [None, None]
        self.plotLims['VIIRS-ANC-Sp-Humd-Surf-Mod-Gran'] = [None, None]
        self.plotLims['VIIRS-ANC-Temp-Skin-Mod-Gran'] = [270.,320.]

        self.dataName = {}
        self.dataName['VIIRS-ANC-Preci-Wtr-Mod-Gran'] = 'data'
        self.dataName['VIIRS-ANC-Temp-Surf2M-Mod-Gran'] = 'data'
        self.dataName['VIIRS-ANC-Wind-Speed-Mod-Gran'] = 'data'
        self.dataName['VIIRS-ANC-Wind-Direction-Mod-Gran'] = 'data'
        self.dataName['VIIRS-ANC-Surf-Ht-Mod-Gran'] = 'data'
        self.dataName['VIIRS-ANC-Press-Surf-Mod-Gran'] = 'data'
        self.dataName['VIIRS-ANC-Tot-Col-Mod-Gran'] = 'data'
        self.dataName['VIIRS-ANC-Optical-Depth-Mod-Gran'] = 'faot550'
        self.dataName['VIIRS-ANC-Geopot-Ht-Lev-Mod-Gran'] = 'data'
        self.dataName['VIIRS-ANC-Sp-Humd-Surf-Mod-Gran'] = 'data'
        self.dataName['VIIRS-ANC-Temp-Skin-Mod-Gran'] = 'data'

        self.blob_dict = get_blob_asc_dict(self.xmlDir,self.blobPath,self.collShortNames)


    def plot_ANC_granules(self,pngDir=None,endian=adl_blob.LITTLE_ENDIAN):

        if pngDir is None :
            pngDir = path.abspath(path.curdir)

        xmlDir = self.xmlDir
        blobPath = self.blobPath

        xmlName = self.xmlName
        plotDescr = self.plotDescr
        plotLims = self.plotLims

        blobDict = self.blob_dict
        collShortNames = blobDict.keys()

        for shortName in collShortNames :

            dataName = self.dataName[shortName]

            for granID in  blobDict[shortName].keys() :
                print '%s --> %s ' % (shortName, granID)
                xmlFile = path.join(xmlDir,xmlName[shortName])
                blobFile = blobDict[shortName][granID][0]
                blobFile = path.join(blobPath,'%s'%(blobFile))
                blobObj = adl_blob.map(xmlFile,blobFile,endian=endian)
                print "dataName = %s" % (dataName)
                blobArrObj = blobObj.as_arrays()
                data = getattr(blobArrObj,dataName)
                print "%s is of kind %r" % (shortName,data.dtype.kind)
                if (data.dtype.kind == 'O') :
                    data = data.astype(np.float)

                vmin,vmax = plotLims[shortName]
                plotTitle = '%s : %s' % (shortName,granID)
                cbTitle = plotDescr[shortName]

                # Create figure with default size, and create canvas to draw on
                scale=1.5
                fig = Figure(figsize=(scale*8,scale*3))
                canvas = FigureCanvas(fig)

                # Create main axes instance, leaving room for colorbar at bottom,
                # and also get the Bbox of the axes instance
                ax_rect = [0.05, 0.18, 0.9, 0.75  ] # [left,bottom,width,height]
                ax = fig.add_axes(ax_rect)

                # Granule axis title
                ax_title = ppl.setp(ax,title=plotTitle)
                ppl.setp(ax_title,fontsize=12)
                ppl.setp(ax_title,family="sans-serif")

                # Plot the data
                im = ax.imshow(ma.masked_less(data,-800.),axes=ax,interpolation='nearest',vmin=vmin,vmax=vmax)
                #im = ax.imshow(ma.masked_less(data,-800.),axes=ax,interpolation='nearest')
                
                ppl.setp(ax.get_xticklabels(), visible=False)
                ppl.setp(ax.get_yticklabels(), visible=False)
                ppl.setp(ax.get_xticklines(),visible=False)
                ppl.setp(ax.get_yticklines(),visible=False)

                # add a colorbar axis
                cax_rect = [0.05 , 0.05, 0.9 , 0.10 ] # [left,bottom,width,height]
                cax = fig.add_axes(cax_rect,frameon=False) # setup colorbar axes

                # Plot the colorbar.
                cb = fig.colorbar(im, cax=cax, orientation='horizontal')
                ppl.setp(cax.get_xticklabels(),fontsize=9)

                # Colourbar title
                cax_title = ppl.setp(cax,title=cbTitle)
                ppl.setp(cax_title,fontsize=10)

                # Redraw the figure
                canvas.draw()

                # Save the figure to a png file...
                pngFile = path.join(pngDir,'%s_%s.png' % (shortName,granID))
                canvas.print_figure(pngFile,dpi=100)

                ppl.close('all')


class GridIPclass():

    def __init__(self,xmlDir,blobPath):

        self.xmlDir = xmlDir
        self.blobPath = path.expanduser(blobPath)

        self.collShortNames = [
                               'VIIRS-GridIP-VIIRS-Qst-Mod-Gran',
                               'VIIRS-GridIP-VIIRS-Lwm-Mod-Gran',
                               'VIIRS-GridIP-VIIRS-Qst-Lwm-Mod-Gran',
                               'VIIRS-GridIP-VIIRS-Nbar-Ndvi-Mod-Gran',
                               'VIIRS-GridIP-VIIRS-Snow-Ice-Cover-Mod-Gran',
                               'VIIRS-I-Conc-IP'
                              ]

        self.xmlName = {}
        self.xmlName['VIIRS-GridIP-VIIRS-Qst-Mod-Gran'] = 'VIIRS_GRIDIP_VIIRS_QST_MOD_GRAN.xml'
        self.xmlName['VIIRS-GridIP-VIIRS-Lwm-Mod-Gran'] = 'VIIRS_GRIDIP_VIIRS_LWM_MOD_GRAN.xml'
        self.xmlName['VIIRS-GridIP-VIIRS-Qst-Lwm-Mod-Gran'] = 'VIIRS_GRIDIP_VIIRS_QST_LWM_MOD_GRAN.xml'
        self.xmlName['VIIRS-GridIP-VIIRS-Nbar-Ndvi-Mod-Gran'] = 'VIIRS_GRIDIP_VIIRS_NBAR_NDVI_MOD_GRAN.xml'
        self.xmlName['VIIRS-GridIP-VIIRS-Snow-Ice-Cover-Mod-Gran'] = 'VIIRS_GRIDIP_VIIRS_SNOW_ICE_COVER_MOD_GRAN.xml'
        self.xmlName['VIIRS-I-Conc-IP'] = 'VIIRS_I_CONC_IP.xml'

        self.plotDescr = {}
        self.plotDescr['VIIRS-GridIP-VIIRS-Qst-Mod-Gran']   =  r'Quarterly Surface Type'
        self.plotDescr['VIIRS-GridIP-VIIRS-Lwm-Mod-Gran'] = 'Land Sea Mask'
        self.plotDescr['VIIRS-GridIP-VIIRS-Qst-Lwm-Mod-Gran']   =  r'Quarterly Surface Type / Land Water Mask'
        self.plotDescr['VIIRS-GridIP-VIIRS-Nbar-Ndvi-Mod-Gran'] = r'Normalized Difference Vegetation Index'
        self.plotDescr['VIIRS-GridIP-VIIRS-Snow-Ice-Cover-Mod-Gran'] = r'Snow Ice Cover'
        self.plotDescr['VIIRS-I-Conc-IP'] = r'Ice Concentration'

        self.plotLims = {}
        self.plotLims['VIIRS-GridIP-VIIRS-Qst-Mod-Gran']   =  [1, 17]
        self.plotLims['VIIRS-GridIP-VIIRS-Lwm-Mod-Gran'] = [0, 7]
        self.plotLims['VIIRS-GridIP-VIIRS-Qst-Lwm-Mod-Gran']   =  [1, 20]
        self.plotLims['VIIRS-GridIP-VIIRS-Nbar-Ndvi-Mod-Gran'] = [None, None]
        #self.plotLims['VIIRS-GridIP-VIIRS-Snow-Ice-Cover-Mod-Gran'] = [None, 200]
        self.plotLims['VIIRS-GridIP-VIIRS-Snow-Ice-Cover-Mod-Gran'] = [0, 1]
        #self.plotLims['VIIRS-GridIP-VIIRS-Snow-Ice-Cover-Mod-Gran'] = [None,None]
        self.plotLims['VIIRS-I-Conc-IP'] = [0., 1.]

        self.dataName = {}
        self.dataName['VIIRS-GridIP-VIIRS-Qst-Mod-Gran']   =  'igbp'
        self.dataName['VIIRS-GridIP-VIIRS-Lwm-Mod-Gran'] = 'landWaterMask'
        self.dataName['VIIRS-GridIP-VIIRS-Qst-Lwm-Mod-Gran']   =  'qstlwm'
        self.dataName['VIIRS-GridIP-VIIRS-Nbar-Ndvi-Mod-Gran'] = 'nbarNdvi'
        self.dataName['VIIRS-GridIP-VIIRS-Snow-Ice-Cover-Mod-Gran'] = 'snowIceCover'
        self.dataName['VIIRS-I-Conc-IP'] = 'iceFraction'

        self.blob_dict = get_blob_asc_dict(self.xmlDir,self.blobPath,self.collShortNames)


    def plot_GridIP_granules(self,pngDir=None,endian=adl_blob.LITTLE_ENDIAN):

        if pngDir is None :
            pngDir = path.abspath(path.curdir)

        xmlDir = self.xmlDir
        blobPath = self.blobPath

        xmlName = self.xmlName
        plotDescr = self.plotDescr
        plotLims = self.plotLims

        blobDict = self.blob_dict
        collShortNames = blobDict.keys()

        for shortName in collShortNames :

            dataName = self.dataName[shortName]

            for granID in  blobDict[shortName].keys() :
                print '%s --> %s ' % (shortName, granID)
                xmlFile = path.join(xmlDir,xmlName[shortName])
                blobFile = blobDict[shortName][granID][0]
                blobFile = path.join(blobPath,'%s'%(blobFile))
                blobObj = adl_blob.map(xmlFile,blobFile,endian=endian)
                print "dataName = %s" % (dataName)
                blobArrObj = blobObj.as_arrays()
                data = getattr(blobArrObj,dataName)
                if (data.dtype.kind == 'O') :
                    data = data.astype(np.float)

                pixelTrimValue = trimObj.sdrTypeFill['ONBOARD_PT_FILL'][data.dtype.name]
                print "pixelTrimValue is %r" % (pixelTrimValue)

                vmin,vmax = plotLims[shortName]
                plotTitle = '%s : %s' % (shortName,granID)
                cbTitle = plotDescr[shortName]

                # Create figure with default size, and create canvas to draw on
                scale=1.5
                fig = Figure(figsize=(scale*8,scale*3))
                canvas = FigureCanvas(fig)

                # Create main axes instance, leaving room for colorbar at bottom,
                # and also get the Bbox of the axes instance
                ax_rect = [0.05, 0.18, 0.9, 0.75  ] # [left,bottom,width,height]
                ax = fig.add_axes(ax_rect)

                # Granule axis title
                ax_title = ppl.setp(ax,title=plotTitle)
                ppl.setp(ax_title,fontsize=12)
                ppl.setp(ax_title,family="sans-serif")

                # Plot the data
                print "%s is of kind %r" % (shortName,data.dtype.kind)
                #print data
                if (data.dtype.kind =='i' or data.dtype.kind =='u'):
                    #data = ma.masked_equal(data,pixelTrimValue)
                    data = ma.masked_greater(data,200)
                else:
                    data = ma.masked_less(data,-800.)

                #im = ax.imshow(ma.masked_equal(data,pixelTrimValue),axes=ax,interpolation='nearest',vmin=vmin,vmax=vmax)
                im = ax.imshow(data,axes=ax,interpolation='nearest',vmin=vmin,vmax=vmax)
                #im = ax.imshow(ma.masked_less(data,-800.),axes=ax,interpolation='nearest')
                
                ppl.setp(ax.get_xticklabels(), visible=False)
                ppl.setp(ax.get_yticklabels(), visible=False)
                ppl.setp(ax.get_xticklines(),visible=False)
                ppl.setp(ax.get_yticklines(),visible=False)

                # add a colorbar axis
                cax_rect = [0.05 , 0.05, 0.9 , 0.10 ] # [left,bottom,width,height]
                cax = fig.add_axes(cax_rect,frameon=False) # setup colorbar axes

                # Plot the colorbar.
                cb = fig.colorbar(im, cax=cax, orientation='horizontal')
                ppl.setp(cax.get_xticklabels(),fontsize=9)

                # Colourbar title
                cax_title = ppl.setp(cax,title=cbTitle)
                ppl.setp(cax_title,fontsize=10)

                # Redraw the figure
                canvas.draw()

                # Save the figure to a png file...
                pngFile = path.join(pngDir,'%s_%s.png' % (shortName,granID))
                canvas.print_figure(pngFile,dpi=100)

                ppl.close('all')


