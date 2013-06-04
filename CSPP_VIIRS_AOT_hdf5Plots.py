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

import ViirsData as ViirsData

import tables as pytables

# every module should have a LOG object
# e.g. LOG.warning('my dog has fleas')
import logging
logging.basicConfig() 

dpi=200

### Moderate and Imager resolution trim table arrays. These are 
### bool arrays, and the trim pixels are set to True.
trimObj = ViirsData.ViirsTrimTable()
modTrimMask = trimObj.createModTrimArray(nscans=48,trimType=bool)

#xmlDir = '/home/geoffc/SSEC/sshfs/rh5b/RH5B/CSPP/EDR/common/ADL/xml/VIIRS'
#oldBlobPath = '/home/geoffc/SSEC/sshfs/rh5b/WORK/sample_data/viirs/edr/input/VIIRS_DB_ADL4.1/blobs/single/ANC'
#newBlobPath = '/home/geoffc/SSEC/sshfs/rh5b/WORK/Work_Area/VIIRS_DB_ADL4.1'


def get_hdf5_dict(hdf5Dir,filePrefix):
    hdf5_asc_Files = {}
    
    hdf5Dir = path.abspath(path.expanduser(hdf5Dir))
    print "hdf5Dir = %s" % (hdf5Dir)

    if path.isdir(hdf5Dir):
        pass
    else :
        hdf5Dir = path.dirname(hdf5Dir)
        print "hdf5Dir = %s" % (hdf5Dir)

    #for prefix in filePrefix :
    #print '%s --> ' % (shortName)
    print "prefix = %s" % (filePrefix)
    fileGlob = path.join(path.expanduser(hdf5Dir),'%s_*.h5'%(filePrefix))
    print "fileGlob = %s" % (fileGlob)

    hdf5Files = glob(fileGlob)
    if hdf5Files != []:
        hdf5Files.sort()
        hdf5Dict = {}
        for files in hdf5Files :
            fileObj = pytables.openFile(files)
            VIIRS_Aeros_Opt_Thick_IP_Gran_0 = fileObj.getNode('/Data_Products/VIIRS-Aeros-Opt-Thick-IP/VIIRS-Aeros-Opt-Thick-IP_Gran_0')
            granID =  getattr(VIIRS_Aeros_Opt_Thick_IP_Gran_0.attrs,'N_Granule_ID')[0][0]
            print 'N_Granule_ID = %s' % (granID)
            dayNightFlag =  getattr(VIIRS_Aeros_Opt_Thick_IP_Gran_0.attrs,'N_Day_Night_Flag')[0][0]
            print 'N_Day_Night_Flag = %s' % (dayNightFlag)
            shortName = fileObj.getNodeAttr('/Data_Products/VIIRS-Aeros-Opt-Thick-IP','N_Collection_Short_Name')[0][0]
            print 'N_Collection_Short_Name = %s' % (shortName)
            hdf5File = path.basename(files)
            hdf5Dict[granID] = [hdf5File,fileObj]
            #fileObj.close()
        hdf5_asc_Files[shortName] = hdf5Dict

    return hdf5_asc_Files


class AOTclass():

    def __init__(self,hdf5Dir):

        self.hdf5Dir = hdf5Dir

        self.collShortNames = [
                               'VIIRS-Aeros-Opt-Thick-IP'
                              ]

        self.plotDescr = {}
        self.plotDescr['VIIRS-Aeros-Opt-Thick-IP'] = r'Aerosol Optical Thickness'

        self.plotLims = {}
        self.plotLims['VIIRS-Aeros-Opt-Thick-IP'] = [-0.05,0.8]

        self.dataName = {}
        self.dataName['VIIRS-Aeros-Opt-Thick-IP'] = '/All_Data/VIIRS-Aeros-Opt-Thick-IP_All/faot550'

        self.hdf5_dict = get_hdf5_dict(hdf5Dir,'IVAOT')


    def plot_AOT_granules(self,pngDir=None,annotation='',dpi=300):

        if pngDir is None :
            pngDir = path.abspath(path.curdir)

        plotDescr = self.plotDescr
        plotLims = self.plotLims

        hdf5_dict = self.hdf5_dict
        collShortNames = hdf5_dict.keys()

        print 'collShortNames = %r' % (collShortNames)

        for shortName in collShortNames :

            print 'shortName = %s' % (shortName)

            dataName = self.dataName[shortName]

            for granID in  hdf5_dict[shortName].keys() :
            #for granID in  [hdf5_dict[shortName].keys()[0]] :

                print '%s --> %s ' % (shortName, granID)

                hdf5Obj = hdf5_dict[shortName][granID][1]
                
                VIIRS_Aeros_Opt_Thick_IP_Gran_0 = hdf5Obj.getNode('/Data_Products/VIIRS-Aeros-Opt-Thick-IP/VIIRS-Aeros-Opt-Thick-IP_Gran_0')
                dayNightFlag =  getattr(VIIRS_Aeros_Opt_Thick_IP_Gran_0.attrs,'N_Day_Night_Flag')[0][0]
                print 'N_Day_Night_Flag = %s' % (dayNightFlag)
                orient = -1 if dayNightFlag == 'Day' else 1

                data = hdf5Obj.getNode(dataName)[:,:]
                NAAPSflag = hdf5Obj.getNode('/All_Data/VIIRS-Aeros-Opt-Thick-IP_All/QF3')
                NAAPSflag = np.bitwise_and(NAAPSflag,12) >> 2
                NAAPSflagMask = ma.masked_greater(NAAPSflag,0).mask

                hdf5Obj.close()

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
                if (data.dtype.kind =='i' or data.dtype.kind =='u'):
                    #data = ma.masked_equal(data,pixelTrimValue)
                    #data = ma.masked_greater(data,200)
                    fill_mask = ma.masked_greater(data,200).mask
                else:
                    #data = ma.masked_less(data,-800.)
                    fill_mask = ma.masked_less(data,-800.).mask

                # Construct the total mask

                totalMask = NAAPSflagMask + fill_mask

                # Mask the aerosol so we only have the retrievals
                data = ma.masked_array(data,mask=totalMask)
                
                im = ax.imshow(data[::orient,::orient],interpolation='nearest',vmin=vmin,vmax=vmax)
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

        hdf5_dict = self.hdf5_dict
        collShortNames = hdf5_dict.keys()

        print 'collShortNames = %r' % (collShortNames)

        CMD = ViirsData.AerosolProdData

        for shortName in collShortNames :

            print 'shortName = %s' % (shortName)

            dataName = self.dataName[shortName]

            for granID in  hdf5_dict[shortName].keys() :
            #for granID in  [hdf5_dict[shortName].keys()[0]] :

                print '%s --> %s ' % (shortName, granID)
                hdf5Obj = hdf5_dict[shortName][granID][1]

                VIIRS_Aeros_Opt_Thick_IP_Gran_0 = hdf5Obj.getNode('/Data_Products/VIIRS-Aeros-Opt-Thick-IP/VIIRS-Aeros-Opt-Thick-IP_Gran_0')
                dayNightFlag =  getattr(VIIRS_Aeros_Opt_Thick_IP_Gran_0.attrs,'N_Day_Night_Flag')[0][0]
                print 'N_Day_Night_Flag = %s' % (dayNightFlag)
                orient = -1 if dayNightFlag == 'Day' else 1

                for byte in [0,1,2,3,4] :
                #for byte in [0] :

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

                        print '\ndSet = %d' %(dSet)

                        dSetName = '/All_Data/VIIRS-Aeros-Opt-Thick-IP_All/QF%s'%(str(byte+1))
                        byteData = hdf5Obj.getNode(dSetName)[:,:]

                        plotTitle = '%s : %s %s (Byte %d)' % (shortName,granID,annotation,byte)
                        fig.text(0.5, 0.95, plotTitle, fontsize=16, color='black', ha='center', va='bottom', alpha=1.0)

                        if (CMD.ViirsAeroQualBitMaskNames[byte][dSet] == 'Spare') :
                            print "Skipping dataset with byte = %d, dSet = %d" % (byte, dSet)
                        else :
                            print "byte = %d, dSet = %d" % (byte, dSet)

                            byteMask = CMD.ViirsAeroQualBitMasks[byte][dSet]
                            byteShift = CMD.ViirsAeroQualBitShift[byte][dSet]

                            print "byteMask = %d, byteShift = %d, dSetName = %s" % (byteMask, byteShift, dSetName)

                            data = np.bitwise_and(byteData,byteMask) >> byteShift
                            vmin,vmax = CMD.ViirsAeroQualvalues[byte][dSet][0], CMD.ViirsAeroQualvalues[byte][dSet][-1]
                            print "vmin = %d, vmax = %d" % (vmin, vmax)

                            cmap = ListedColormap(CMD.ViirsAeroQualFillColours[byte][dSet])

                            numCats = np.array(CMD.ViirsAeroQualFillColours[byte][dSet]).size
                            numBounds = numCats + 1

                            tickPos = np.arange(float(numBounds))/float(numCats)
                            tickPos = tickPos[0 :-1] + tickPos[1]/2.

                            print "numCats = ",numCats
                            print "tickPos = ",tickPos

                            titleStr = CMD.ViirsAeroQualBitMaskNames[byte][dSet]
                            print "titleStr = %s" % (titleStr)

                            Ax.append(fig.add_subplot(numRows,2,plotIdx))
                            print "data.dtype.__str__() = %s" % (data.dtype.__str__())
                            Im.append(Ax[dSet].imshow(data.astype('int')[::orient,::orient], vmin=vmin, vmax=vmax, interpolation='nearest',cmap=cmap))
                            Txt.append(Ax[dSet].set_title(titleStr))

                            ppl.setp(Txt[dSet],fontsize=10)
                            ppl.setp(Ax[dSet].get_xticklines(), visible=False)
                            ppl.setp(Ax[dSet].get_yticklines(), visible=False)
                            ppl.setp(Ax[dSet].get_xticklabels(), visible=False)
                            ppl.setp(Ax[dSet].get_yticklabels(), visible=False)

                            Cb.append(fig.colorbar(Im[dSet], orientation='horizonal', pad=0.05))

                            print "Cb byte = %d, dSet = %d" % (byte, dSet)
                            print "CMD.ViirsAeroQualTickNames[%d][%d] = %s" % \
                                    (byte,dSet,CMD.ViirsAeroQualTickNames[byte][dSet])
                            print "CMD.ViirsAeroQualFillColours[%d][%d] = %s" % \
                                    (byte,dSet,CMD.ViirsAeroQualFillColours[byte][dSet])

                            Cb[dSet].set_ticks(vmax*tickPos)
                            ppl.setp(Cb[dSet].ax,xticklabels=CMD.ViirsAeroQualTickNames[byte][dSet])
                            ppl.setp(Cb[dSet].ax.get_xticklabels(),fontsize=6)
                            ppl.setp(Cb[dSet].ax.get_xticklines(),visible=False)

                            plotIdx += 1

                    pngFile = path.join(pngDir,'%s_%s_QF%s.png' % (shortName,granID,str(byte+1)))
                    print "Writing to %s..." % (pngFile)

                    canvas.draw()
                    canvas.print_figure(pngFile,dpi=dpi)

                    ppl.close('all')

                hdf5Obj.close()



