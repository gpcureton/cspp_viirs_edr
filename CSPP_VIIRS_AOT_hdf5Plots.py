#!/usr/bin/env python
# encoding: utf-8
"""
CSPP_VIIRS_AOT_hdf5Plots.py

Purpose: Create swath projection quicklook PNGs from the VIIRS AOT IP HDF5 files.
         Images can be created for the IP product or the associated quality flags.

Minimum commandline...

#export CSPP_HOME=/path/to/CSPP
#source $CSPP_HOME/cspp_edr_env.sh
#source $CSPP_HOME/common/cspp_common.sh

export CSPP_EDR_HOME=$(readlink -f /path/to/EDR)
source $CSPP_EDR_HOME/cspp_edr_env.sh

python CSPP_VIIRS_AOT_hdf5Plots.py -g '/path/to/files/GMTCO*.h5' -i '/path/to/files/IICMO*.h5' -p VCM

  or

python CSPP_VIIRS_AOT_hdf5Plots.py --geo_file=/path/to/files/GMTCO*.h5 --ip_file=/path/to/files/IICMO*.h5 -p VCM


Created by Geoff Cureton on 2013-06-04.
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

#############

import string, sys
from glob import glob
from os import path,uname
from time import time
import traceback

import numpy as np
from  numpy import ma as ma
import scipy as scipy

import matplotlib
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap
from matplotlib.figure import Figure

matplotlib.use('Agg')
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

# This must come *after* the backend is specified.
import matplotlib.pyplot as ppl

import optparse as optparse

import ViirsData

import tables as pytables
from tables import exceptions as pyEx


# every module should have a LOG object
# e.g. LOG.warning('my dog has fleas')
import logging
logging.basicConfig() 

dpi=200

### Moderate and Imager resolution trim table arrays. These are 
### bool arrays, and the trim pixels are set to True.
trimObj = ViirsData.ViirsTrimTable()
modTrimMask = trimObj.createModTrimArray(nscans=48,trimType=bool)


def get_hdf5_dict(hdf5Path,filePrefix):
    hdf5_asc_Files = {}
    
    hdf5Path = path.abspath(path.expanduser(hdf5Path))
    print "hdf5Path = %s" % (hdf5Path)

    hdf5Dir = path.dirname(hdf5Path)
    hdf5Glob = path.basename(hdf5Path)

    print "hdf5Dir = %s" % (hdf5Dir)
    print "hdf5Glob = %s" % (hdf5Glob)

    if (hdf5Glob == '' or hdf5Glob == '*'):
        print "prefix = %s" % (filePrefix)
        hdf5Glob = path.join(hdf5Dir,'%s_*.h5'%(filePrefix))
    else :
        hdf5Glob = path.join(hdf5Dir,'%s'%(hdf5Glob))

    print "Final hdf5Glob = %s" % (hdf5Glob)

    hdf5Files = glob(hdf5Glob)
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
        self.plotDescr['VIIRS-Aeros-Opt-Thick-IP'] = ['AOT @ 550 nm', \
                                                      'Angstrom Exponent', \
                                                      'Slant AOT @ 550 nm']

        self.plotLims = {}
        #self.plotLims['VIIRS-Aeros-Opt-Thick-IP'] = [-0.05,0.8]
        self.plotLims['VIIRS-Aeros-Opt-Thick-IP'] = [None,None]

        self.dataName = {}
        self.dataName['VIIRS-Aeros-Opt-Thick-IP'] = ['/All_Data/VIIRS-Aeros-Opt-Thick-IP_All/faot550', \
                                                     '/All_Data/VIIRS-Aeros-Opt-Thick-IP_All/angexp', \
                                                     '/All_Data/VIIRS-Aeros-Opt-Thick-IP_All/aotSlant550']

        self.hdf5_dict = get_hdf5_dict(hdf5Dir,'IVAOT')


    def plot_AOT_granules(self,plotProd='IP',vmin=None,vmax=None,pngDir=None,pngPrefix=None,annotation='',dpi=300):

        if pngDir is None :
            pngDir = path.abspath(path.curdir)

        plotDescr = self.plotDescr
        plotLims = self.plotLims

        hdf5_dict = self.hdf5_dict
        collShortNames = hdf5_dict.keys()

        print 'collShortNames = %r' % (collShortNames)

        for shortName in collShortNames :

            print 'shortName = %s' % (shortName)

            if (plotProd == 'IP'):

                dataNames = self.dataName[shortName]
                plotDescrs = plotDescr[shortName]
                prodNames = ['faot550','angexp','aotSlant550']

            elif (plotProd == 'faot550'):

                dataNames = [self.dataName[shortName][0]]
                plotDescrs = [plotDescr[shortName][0]]
                prodNames = ['faot550']

            elif (plotProd == 'angexp'):

                dataNames = [self.dataName[shortName][1]]
                plotDescrs = [plotDescr[shortName][1]]
                prodNames = ['angexp']

            elif (plotProd == 'aotSlant550'):

                dataNames = [self.dataName[shortName][2]]
                plotDescrs = [plotDescr[shortName][2]]
                prodNames = ['aotSlant550']


            for granID in  hdf5_dict[shortName].keys() :

                print '%s --> %s ' % (shortName, granID)

                hdf5Obj = hdf5_dict[shortName][granID][1]
                
                VIIRS_Aeros_Opt_Thick_IP_Gran_0 = hdf5Obj.getNode('/Data_Products/VIIRS-Aeros-Opt-Thick-IP/VIIRS-Aeros-Opt-Thick-IP_Gran_0')
                dayNightFlag =  getattr(VIIRS_Aeros_Opt_Thick_IP_Gran_0.attrs,'N_Day_Night_Flag')[0][0]
                print 'N_Day_Night_Flag = %s' % (dayNightFlag)
                orient = -1 if dayNightFlag == 'Day' else 1

                for dataName,plotDescr,prodName in zip(dataNames,plotDescrs,prodNames):

                    data = hdf5Obj.getNode(dataName)[:,:]
                    NAAPSflag = hdf5Obj.getNode('/All_Data/VIIRS-Aeros-Opt-Thick-IP_All/QF3')
                    NAAPSflag = np.bitwise_and(NAAPSflag,12) >> 2
                    NAAPSflagMask = ma.masked_greater(NAAPSflag,0).mask

                    pixelTrimValue = trimObj.sdrTypeFill['ONBOARD_PT_FILL'][data.dtype.name]
                    print "pixelTrimValue is %r" % (pixelTrimValue)

                    # Apply the moderate pixel trim, so that we can properly mask them out at plot time.
                    data = ma.array(data,mask=modTrimMask,fill_value=trimObj.sdrTypeFill['ONBOARD_PT_FILL'][data.dtype.name])
                    data = data.filled()

                    #vmin,vmax = plotLims[shortName]
                    plotTitle = '%s : %s %s' % (shortName,granID,annotation)
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
                    print "%s is of kind %r" % (shortName,data.dtype.kind)
                    if (data.dtype.kind =='i' or data.dtype.kind =='u'):
                        fill_mask = ma.masked_greater(data,200).mask
                    else:
                        fill_mask = ma.masked_less(data,-800.).mask

                    # Construct the total mask

                    totalMask = NAAPSflagMask + fill_mask

                    # Mask the aerosol so we only have the retrievals
                    data = ma.masked_array(data,mask=totalMask)
                    
                    im = ax.imshow(data[::orient,::orient],interpolation='nearest',vmin=vmin,vmax=vmax)
                    
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
                    pngFile = path.join(pngDir,'%s%s_%s_%s.png' % (pngPrefix,shortName,granID,prodName))
                    canvas.print_figure(pngFile,dpi=dpi)

                    ppl.close('all')

                hdf5Obj.close()


    def plot_AOT_tests(self,plotProd='QF',pngDir=None,pngPrefix=None,annotation='',dpi=300):

        if pngDir is None :
            pngDir = path.abspath(path.curdir)

        hdf5_dict = self.hdf5_dict
        collShortNames = hdf5_dict.keys()

        if (plotProd == 'QF'):
            byteList = [0,1,2,3,4]
        else :
            byteList = [int(plotProd.strip('QF'))]

        print 'collShortNames = %r' % (collShortNames)

        CMD = ViirsData.AerosolProdData

        for shortName in collShortNames :

            print 'shortName = %s' % (shortName)

            dataName = self.dataName[shortName]

            for granID in  hdf5_dict[shortName].keys() :

                print '%s --> %s ' % (shortName, granID)
                hdf5Obj = hdf5_dict[shortName][granID][1]

                VIIRS_Aeros_Opt_Thick_IP_Gran_0 = hdf5Obj.getNode('/Data_Products/VIIRS-Aeros-Opt-Thick-IP/VIIRS-Aeros-Opt-Thick-IP_Gran_0')
                dayNightFlag =  getattr(VIIRS_Aeros_Opt_Thick_IP_Gran_0.attrs,'N_Day_Night_Flag')[0][0]
                print 'N_Day_Night_Flag = %s' % (dayNightFlag)
                orient = -1 if dayNightFlag == 'Day' else 1

                for byte in byteList :

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

                    pngFile = path.join(pngDir,'%s%s_%s_QF%s.png' % (pngPrefix,shortName,granID,str(byte+1)))
                    print "Writing to %s..." % (pngFile)

                    canvas.draw()
                    canvas.print_figure(pngFile,dpi=dpi)

                    ppl.close('all')

                hdf5Obj.close()


###################################################
#                  Main Function                  #
###################################################

def main():

    prodChoices=['IP','QF','faot550','angexp','aotslant550','QF1','QF2','QF3','QF4','QF5']

    description = \
    '''
    This is a brief description of %prog
    '''
    usage = "usage: %prog [mandatory args] [options]"
    version = version="%prog"
    parser = optparse.OptionParser(description=description,usage=usage,version=version)

    # Mandatory arguments
    mandatoryGroup = optparse.OptionGroup(parser, "Mandatory Arguments",
                        "At a minimum these arguments must be specified")

    mandatoryGroup.add_option('--dir',
                      action="store",
                      dest="hdf5Dir" ,
                      type="string",
                      help="The directory containing the VSSTO HDF5 files.")

    parser.add_option_group(mandatoryGroup)

    # Optional arguments
    optionalGroup = optparse.OptionGroup(parser, "Extra Options",
                        "These options may be used to customize plot characteristics.")

    optionalGroup.add_option('-r','--svn_revision',
                      action="store",
                      dest="svnRevision",
                      default=string.split(__version__," ")[2],
                      type="string",
                      help="The Subversion revision number/tag of this script")
    optionalGroup.add_option('--plotMin',
                      action="store",
                      type="float",
                      dest="plotMin",
                      help="Minimum value to plot.")
    optionalGroup.add_option('--plotMax',
                      action="store",
                      type="float",
                      dest="plotMax",
                      help="Maximum value to plot.")
    optionalGroup.add_option('-d','--dpi',
                      action="store",
                      dest="dpi",
                      default='200.',
                      type="float",
                      help="The resolution in dots per inch of the output png file. [default: %default]")
    optionalGroup.add_option('-a','--map_annotation',
                      action="store",
                      dest="mapAnn",
                      #default='',
                      type="string",
                      help="The map legend describing the dataset being shown. [default: IPPROD]")
    optionalGroup.add_option('-p','--product',
                      action="store",
                      dest="plotProduct",
                      type="choice",
                      choices=prodChoices,
                      help='''The VIIRS AOT IP or QF datasets to plot.\n\n
                           Possible values are...
                           %s
                           ''' % (prodChoices.__str__()[1:-1]))
    optionalGroup.add_option('--png_dir',
                      action="store",
                      dest="pngDir" ,
                      type="string",
                      help="The directory where png files will be written.")
    optionalGroup.add_option('-o','--output_file_prefix',
                      action="store",
                      dest="outputFilePrefix",
                      default="",
                      type="string",
                      help="""String to prefix to the automatically generated png names, which are of
the form <N_Collection_Short_Name>_<N_Granule_ID>_<dset>.png. [default: %default]""")

    parser.add_option_group(optionalGroup)

    # Parse the arguments from the command line
    (options, args) = parser.parse_args()


    # Check that all of the mandatory options are given. If one or more 
    # are missing, print error message and exit...
    mandatories = ['hdf5Dir']
    mand_errors = ["Missing mandatory argument [-d DIR | --dir=DIR]"
                  ]
    isMissingMand = False
    for m,m_err in zip(mandatories,mand_errors):
        if not options.__dict__[m]:
            isMissingMand = True
            print m_err
    if isMissingMand :
        parser.error("Incomplete mandatory arguments, aborting...")

    vmin = options.plotMin
    vmax = options.plotMax

    hdf5Path = path.abspath(path.expanduser(options.hdf5Dir))

    print "hdf5Path = %s" % (hdf5Path)
    pngDir = '.' if (options.pngDir is None) else options.pngDir
    pngDir = path.abspath(path.expanduser(pngDir))
    print "pngDir = %s" % (pngDir)

    pngPrefix = options.outputFilePrefix
    dpi = options.dpi
    plotProduct = options.plotProduct

    plotEDR = False
    plotQF = False

    if (plotProduct is None):
        plotEDR = True
        plotQF = True
    else :
        if ('IP' in plotProduct) \
           or ('faot550' in plotProduct) \
           or ('angexp' in plotProduct) \
           or ('aotslant550' in plotProduct) :
               plotEDR = True

        if ('QF' in plotProduct) :
            plotQF = True

    if plotEDR :
        try :
            AOTobj = AOTclass(hdf5Path)
            AOTobj.plot_AOT_granules(plotProd=plotProduct,vmin=vmin,vmax=vmax,pngDir=pngDir,pngPrefix=pngPrefix,dpi=dpi)
            pytables.file.close_open_files()
        except Exception, err:
            print "%s" % (str(err))
            pytables.file.close_open_files()

    if plotQF :
        try :
            AOTobj = AOTclass(hdf5Path)
            AOTobj.plot_AOT_tests(plotProd=plotProduct,pngDir=pngDir,pngPrefix=pngPrefix,dpi=dpi)
        except Exception, err:
            print "%s" % (str(err))
            pytables.file.close_open_files()



    print "Exiting..."
    sys.exit(0)


if __name__ == '__main__':
    main()

