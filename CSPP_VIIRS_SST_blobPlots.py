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
from VIIRS import ViirsData as ViirsData

# every module should have a LOG object
# e.g. LOG.warning('my dog has fleas')
import logging
logging.basicConfig() 

dpi=200

### Moderate and Imager resolution trim table arrays. These are 
### bool arrays, and the trim pixels are set to True.
trimObj = ViirsData.ViirsTrimTable()
modTrimMask = trimObj.createModTrimArray(nscans=48,trimType=bool)

xmlDir = '/home/geoffc/SSEC/sshfs/rh5b/RH5B/CSPP/EDR/common/ADL/xml/VIIRS'
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


class SSTclass():

    def __init__(self,xmlDir,blobPath):

        self.xmlDir = xmlDir
        self.blobPath = blobPath

        self.collShortNames = [
                               'VIIRS-SST-EDR',
                              ]

        self.xmlName = {}
        self.xmlName['VIIRS-SST-EDR'] = 'VIIRS_SST_EDR.xml'

        self.plotDescr = {}
        self.plotDescr['VIIRS-SST-EDR'] = ['Skin Sea Surface Temperature (K)','Bulk Sea Surface Temperature (K)']

        self.plotLims = {}
        #self.plotLims['VIIRS-SST-EDR'] = [250., 290.]
        self.plotLims['VIIRS-SST-EDR'] = [None,None]

        self.dataName = {}
        self.dataName['VIIRS-SST-EDR'] = ['skinSST','bulkSST']
        self.dataScales = {}
        self.dataScales['VIIRS-SST-EDR'] = ['Skin_Scale','Bulk_Scale']
        self.dataOffsets = {}
        self.dataOffsets['VIIRS-SST-EDR'] = ['Skin_Offset','Bulk_Offset']

        dataFlags =['sstFlags0','sstFlags1','sstFlags2','sstFlags3']

        self.blob_dict = get_blob_asc_dict(self.xmlDir,self.blobPath,self.collShortNames)


    def plot_SST_granules(self,pngDir=None,endian=adl_blob.LITTLE_ENDIAN):

        if pngDir is None :
            pngDir = path.abspath(path.curdir)

        xmlDir = self.xmlDir
        blobPath = self.blobPath

        xmlName = self.xmlName
        #plotDescr = self.plotDescr
        plotLims = self.plotLims

        blobDict = self.blob_dict
        collShortNames = blobDict.keys()

        for shortName in collShortNames :

            dataName = self.dataName[shortName]

            for granID in  blobDict[shortName].keys() :
                print '%s --> %s ' % (shortName, granID)
                xmlFile = path.join(xmlDir,xmlName[shortName])
                print '%s' % (xmlFile)
                blobFile = blobDict[shortName][granID][0]
                blobFile = path.join(blobPath,'%s'%(blobFile))
                print '%s' % (blobFile)
                blobObj = adl_blob.map(xmlFile,blobFile,endian=endian)
                blobArrObj = blobObj.as_arrays()

                for dataName,dataScales,dataOffsets,plotDescr in zip(self.dataName[shortName],self.dataScales[shortName],self.dataOffsets[shortName],self.plotDescr[shortName]):
                    print dataName,dataScales,dataOffsets

                    data = getattr(blobArrObj,dataName)
                    if (data.dtype.kind == 'O') :
                        data = data.astype(np.float)
                    scale = getattr(blobArrObj,dataScales)
                    offset = getattr(blobArrObj,dataOffsets)

                    print "Scale, Offset = %f %f" % (scale,offset)

                    data = data * scale + offset

                    pixelTrimValue = trimObj.sdrTypeFill['ONGROUND_PT_FILL'][data.dtype.name]
                    print "pixelTrimValue is %r" % (pixelTrimValue)

                    # Apply the moderate pixel trim, so that we can properly mask them out at plot time.
                    data = ma.array(data,mask=modTrimMask,fill_value=trimObj.sdrTypeFill['ONBOARD_PT_FILL'][data.dtype.name])
                    data = data.filled()

                    vmin,vmax = plotLims[shortName]
                    plotTitle = '%s : %s' % (shortName,granID)
                    cbTitle = plotDescr

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
                    pngFile = path.join(pngDir,'%s_%s_%s.png' % (shortName,dataName,granID))
                    canvas.print_figure(pngFile,dpi=100)

                    ppl.close('all')

#def plot_SST_granules_HDF5:
#from os import path
#import numpy as np
#from numpy import ma
#from matplotlib import pyplot as ppl
#from glob import glob
#import tables as pytables
#import ViirsData

#trimObj = ViirsData.ViirsTrimTable()

#trimFillValue = trimObj.sdrTypeFill['ONBOARD_PT_FILL']['float32']
#modTrimMask = trimObj.createModTrimArray(nscans=4*48,trimType=bool)

#naFillValue_uint16 = trimObj.sdrTypeFill['NA_FILL']['uint16']
#naFillValue_float32 = trimObj.sdrTypeFill['NA_FILL']['float32']

#SSTfiles = glob('/home/geoffc/SSEC/sshfs/jpss-cloud/data/rayg/viirs_gtm_test/VSSTO*.h5')
#SSTfiles.sort()

#for files in SSTfiles :
    #print "Opening %s ..." % (files)
    #SSTobj = pytables.openFile(files)

    #skinSST = SSTobj.getNode('/All_Data/VIIRS-SST-EDR_All/SkinSST')[:,:]
    #skinSSTfactors = SSTobj.getNode('/All_Data/VIIRS-SST-EDR_All/SkinSSTFactors')[:]

    #bulkSST = SSTobj.getNode('/All_Data/VIIRS-SST-EDR_All/BulkSST')[:,:]
    #bulkSSTfactors = SSTobj.getNode('/All_Data/VIIRS-SST-EDR_All/BulkSSTFactors')[:]

    #SSTobj.close()

    #naFillMask_skin = ma.masked_equal(skinSST,naFillValue_uint16).mask

    #skinSST = skinSST * skinSSTfactors[0] + skinSSTfactors[1]

    #skinSST = ma.array(skinSST,mask=modTrimMask,fill_value=trimFillValue)
    #skinSST = skinSST.filled()
    #skinSST = ma.array(skinSST,mask=naFillMask_skin,fill_value=naFillValue_float32)
    #skinSST = skinSST.filled()

    #f = ppl.figure();
    #ppl.imshow(ma.masked_less(skinSST,-800.),interpolation='nearest');
    #ppl.colorbar(orientation='horizontal');
    #fileName = path.basename(files)
    #f.savefig('/home/geoffc/SSEC/ADL/CSPP_granTest/PEATE/SST/'+fileName+'_skin.png',dpi=200)
    #ppl.close('all')

    #naFillMask_bulk = ma.masked_equal(bulkSST,naFillValue_uint16).mask

    #bulkSST = bulkSST * bulkSSTfactors[0] + bulkSSTfactors[1]

    #bulkSST = ma.array(bulkSST,mask=modTrimMask,fill_value=trimFillValue)
    #bulkSST = bulkSST.filled()
    #bulkSST = ma.array(bulkSST,mask=naFillMask_bulk,fill_value=naFillValue_float32)
    #bulkSST = bulkSST.filled()

    #f = ppl.figure();
    #ppl.imshow(ma.masked_less(bulkSST,-800.),interpolation='nearest');
    #ppl.colorbar(orientation='horizontal');
    #fileName = path.basename(files)
    #f.savefig('/home/geoffc/SSEC/ADL/CSPP_granTest/PEATE/SST/'+fileName+'_bulk.png',dpi=200)
    #ppl.close('all')

