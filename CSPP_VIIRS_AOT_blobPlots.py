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


class AOTclass():

    def __init__(self,xmlDir,blobPath,endian=adl_blob.LITTLE_ENDIAN):

        self.xmlDir = xmlDir
        self.blobPath = blobPath

        self.collShortNames = [
                               'VIIRS-Aeros-Opt-Thick-IP'
                              ]

        self.xmlName = {}
        self.xmlName['VIIRS-Aeros-Opt-Thick-IP'] = 'VIIRS_AEROS_OPT_THICK_IP.xml'

        self.plotDescr = {}
        self.plotDescr['VIIRS-Aeros-Opt-Thick-IP'] = r'Aerosol Optical Thickness'

        self.plotLims = {}
        self.plotLims['VIIRS-Aeros-Opt-Thick-IP'] = [-0.05,0.8]

        self.dataName = {}
        self.dataName['VIIRS-Aeros-Opt-Thick-IP'] = 'faot550'

        self.blob_dict = get_blob_asc_dict(self.xmlDir,self.blobPath,self.collShortNames)

        for granID in  self.blob_dict['VIIRS-Aeros-Opt-Thick-IP'].keys() :
            xmlFile = path.join(xmlDir,self.xmlName['VIIRS-Aeros-Opt-Thick-IP'])
            blobFile = self.blob_dict['VIIRS-Aeros-Opt-Thick-IP'][granID][0]
            blobFile = path.join(blobPath,'%s'%(blobFile))
            blobObj = adl_blob.map(xmlFile,blobFile,endian=endian)
            self.blob_dict['VIIRS-Aeros-Opt-Thick-IP'][granID].append(blobObj)


    def plot_AOT_granules(self,pngDir=None,annotation='',dpi=300):

        if pngDir is None :
            pngDir = path.abspath(path.curdir)

        plotDescr = self.plotDescr
        plotLims = self.plotLims

        blobDict = self.blob_dict
        collShortNames = blobDict.keys()

        for shortName in collShortNames :

            dataName = self.dataName[shortName]

            for granID in  blobDict[shortName].keys() :
                print '%s --> %s ' % (shortName, granID)
                blobObj = blobDict[shortName][granID][2]
                blobArrObj = blobObj.as_arrays()
                data = getattr(blobArrObj,dataName)

                pixelTrimValue = trimObj.sdrTypeFill['ONBOARD_PT_FILL'][data.dtype.name]
                print "pixelTrimValue is %r" % (pixelTrimValue)

                # Apply the moderate pixel trim, so that we can properly mask them out at plot time.
                data = ma.array(data,mask=modTrimMask,fill_value=trimObj.sdrTypeFill['ONBOARD_PT_FILL'][data.dtype.name])
                data = data.filled()

                vmin,vmax = plotLims[shortName]
                plotTitle = '%s : %s %s' % (shortName,granID,annotation)
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

                im = ax.imshow(ma.masked_less(data,0),interpolation='nearest',vmin=vmin,vmax=vmax)
                #im = ax.imshow(ma.masked_equal(data,pixelTrimValue),axes=ax,interpolation='nearest',vmin=vmin,vmax=vmax)
                #im = ax.imshow(data,axes=ax,interpolation='nearest',vmin=vmin,vmax=vmax)
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
                ppl.setp(cax.get_xticklines(),visible=True)

                # Colourbar title
                cax_title = ppl.setp(cax,title=cbTitle)
                ppl.setp(cax_title,fontsize=10)

                # Turn off the tickmarks on the colourbar
                #ppl.setp(cb.ax.get_xticklines(),visible=False)
                #ppl.setp(cb.ax.get_xticklabels(),fontsize=9)

                # Redraw the figure
                canvas.draw()

                # Save the figure to a png file...
                pngFile = path.join(pngDir,'%s_%s.png' % (shortName,granID))
                canvas.print_figure(pngFile,dpi=dpi)

                ppl.close('all')


    def plot_AOT_tests(self,pngDir=None,annotation='',dpi=300):

        if pngDir is None :
            pngDir = path.abspath(path.curdir)

        plotDescr = self.plotDescr
        plotLims = self.plotLims

        blobDict = self.blob_dict
        collShortNames = blobDict.keys()

        print 'collShortNames = %r' % (collShortNames)

        CMD = ViirsData.AerosolProdData

        for shortName in collShortNames :

            print 'shortName = %s' % (shortName)

            dataName = self.dataName[shortName]

            for granID in  blobDict[shortName].keys() :
                print '%s --> %s ' % (shortName, granID)
                blobObj = blobDict[shortName][granID][2]
                blobArrObj = blobObj.as_arrays()

                for byte in [0,1,2,3,4] :

                    print ""

                    numPlots = len(CMD.ViirsAeroQualBitMasks[byte])
                    numRows = np.ceil(float(numPlots)/2.)

                    figWidth = 15. # inches
                    figHeight = 3. * numPlots/2. # inches
                    fig = Figure(figsize=((figWidth,figHeight)))
                    canvas = FigureCanvas(fig)

                    Ax = []
                    Im = []
                    Txt = []
                    Cb = []

                    plotIdx = 1

                    for dSet in range(np.shape(CMD.ViirsAeroQualBitMasks[byte])[0]):

                        print 'dSet = ',dSet

                        dSetName = 'qf'+str(byte+1)

                        plotTitle = '%s : %s %s (Byte %d)' % (shortName,granID,annotation,byte)
                        fig.text(0.5, 0.95, plotTitle, fontsize=16, color='black', ha='center', va='bottom', alpha=1.0)

                        if (CMD.ViirsAeroQualBitMaskNames[byte][dSet] == 'Spare') :
                            print "Skipping dataset with byte = %d, dSet = %d" % (byte, dSet)
                        else :
                            print "byte = %d, dSet = %d" % (byte, dSet)

                            byteMask = CMD.ViirsAeroQualBitMasks[byte][dSet]
                            byteShift = CMD.ViirsAeroQualBitShift[byte][dSet]

                            print "byteMask = %d, byteShift = %d, dSetName = %s" % (byteMask, byteShift, dSetName)

                            data = np.bitwise_and(getattr(blobArrObj,dSetName)[:,:],byteMask) >> byteShift
                            vmin,vmax = CMD.ViirsAeroQualvalues[byte][dSet][0], CMD.ViirsAeroQualvalues[byte][dSet][-1]
                            print "vmin = %d, vmax = %d" % (vmin, vmax)

                            titleStr = CMD.ViirsAeroQualBitMaskNames[byte][dSet]
                            print "titleStr = %s" % (titleStr)

                            Ax.append(fig.add_subplot(numRows,2,plotIdx))
                            print "data.dtype.__str__() = %s" % (data.dtype.__str__())
                            Im.append(Ax[dSet].imshow(data.astype('int')[::-1,::-1], vmin=vmin, vmax=vmax, interpolation='nearest'))
                            Txt.append(Ax[dSet].set_title(titleStr))

                            ppl.setp(Txt[dSet],fontsize=10)
                            ppl.setp(Ax[dSet].get_xticklines(), visible=False)
                            ppl.setp(Ax[dSet].get_yticklines(), visible=False)
                            ppl.setp(Ax[dSet].get_xticklabels(), visible=False)
                            ppl.setp(Ax[dSet].get_yticklabels(), visible=False)

                            Cb.append(fig.colorbar(Im[dSet], orientation='horizonal', pad=0.05))

                            plotIdx += 1

                    pngFile = path.join(pngDir,'%s_%s_%s.png' % (shortName,granID,dSetName))
                    print "Writing to %s..." % (pngFile)

                    canvas.draw()
                    canvas.print_figure(pngFile,dpi=dpi)

                    ppl.close('all')


