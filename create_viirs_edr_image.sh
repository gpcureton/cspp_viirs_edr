#!/bin/bash
# $Id$
# Create quicklook PNGs for VIIRS EDR products.
#
# Environment settings:
# CSPP_EDR_HOME : the location of the CSPP_RT directory
#
# Copyright 2011-2013, University of Wisconsin Regents.
# Licensed under the GNU GPLv3.

if [ -z "$CSPP_EDR_HOME" ]; then
    echo "CSPP_EDR_HOME is not set, but is required for this script to operate."
    exit 9
fi

. ${CSPP_EDR_HOME}/cspp_edr_runtime.sh

#
# Gather the various command line options...
#

GEO_FILES_OPT=
IP_FILES_OPT=
PROD_OPT=
PLOT_MIN_OPT=
PLOT_MAX_OPT=
DPI_OPT=
SCALE_OPT=
LAT_0_OPT=
LON_0_OPT=
STRIDE_OPT=
SCATTER_PLOT_OPT=
POINTSIZE_OPT=
MAP_RES_OPT=
MAP_ANNOTATION_OPT=
OUTPUT_FILE_OPT=

#echo $@

OPTS=`getopt -o "g:i:p:d:s:S:P:m:a:o:vh" -l "geo_files:,ip_files:,product:,plotMin:,plotMax:,dpi:,scale:,lat_0:,lon_0:,stride:,scatter_plot,pointSize:,map_res:,map_annotation:,output_file:,verbose,help" -- "$@"`

# From viirs_edr.sh
#OPTS=`getopt -o "i:w:p:dvhaz" -l "input_files:,alg:,work_directory:,processors:,anc_endianness:,sdr_endianness:,skip_sdr_unpack,skip_aux_linking,skip_ancillary,skip_algorithm,zip,aggregate,no_dummy_granules,debug,no_chain,verbose,help" -- "$@"`

# If returncode from getopt != 0, exit with error.
if [ $? != 0 ]
then
    echo "There was an error with the command line parameters to viirs_edr.sh, aborting..."
    exit 1
fi

# A little magic
eval set -- "$OPTS"

# Now go through all the options
haveFlag=0
helpFlag=0
usageFlag=0

while true ;
do
    case "$1" in

        ### Mandatory

        -g|--geo_files)
            GEO_FILES_OPT="--geo_file=$2"
            #echo "Setting INPUT_FILES_OPT"
            haveFlag=1
            shift 2;;

        -i|--ip_files)
            IP_FILES_OPT="--ip_file=$2"
            #echo "Setting INPUT_FILES_OPT"
            haveFlag=1
            shift 2;;

        -p|--product)
            PROD_OPT="--product=$2"
            #echo "Setting ALG_OPT"
            haveFlag=1
            shift 2;;

        ### Optional

        --plotMin)
            PLOT_MIN_OPT="--plotMin=$2"
            #echo "Setting WORK_DIR_OPT"
            haveFlag=1
            shift 2;;

        --plotMax)
            PLOT_MAX_OPT="--plotMax=$2"
            #echo "Setting WORK_DIR_OPT"
            haveFlag=1
            shift 2;;

        -d|--dpi)
            DPI_OPT="--dpi=$2"
            #echo "Setting WORK_DIR_OPT"
            haveFlag=1
            shift 2;;

        -s|--scale)
            SCALE_OPT="--scale=$2"
            #echo "Setting WORK_DIR_OPT"
            haveFlag=1
            shift 2;;

        --lat_0)
            LAT_0_OPT="--lat_0=$2"
            #echo "Setting WORK_DIR_OPT"
            haveFlag=1
            shift 2;;

        --lon_0)
            LON_0_OPT="--lon_0=$2"
            #echo "Setting WORK_DIR_OPT"
            haveFlag=1
            shift 2;;

        -S|--stride)
            STRIDE_OPT="--stride=$2"
            #echo "Setting WORK_DIR_OPT"
            haveFlag=1
            shift 2;;

        --scatter_plot)
            SCATTER_PLOT_OPT="--scatter_plot"
            #echo "Setting SKIP_SDR_UNPACK_OPT"
            haveFlag=1
            shift ;;

        -P|--pointSize)
            POINTSIZE_OPT="--pointSize=$2"
            #echo "Setting WORK_DIR_OPT"
            haveFlag=1
            shift 2;;

        -m|--map_res)
            MAP_RES_OPT="--map_res=$2"
            #echo "Setting WORK_DIR_OPT"
            haveFlag=1
            shift 2;;

        -a|--map_annotation)
            MAP_ANNOTATION_OPT="--map_annotation=$2"
            #echo "Setting WORK_DIR_OPT"
            haveFlag=1
            shift 2;;

        -o|--output_file)
            OUTPUT_FILE_OPT="--output_file=$2"
            #echo "Setting WORK_DIR_OPT"
            haveFlag=1
            shift 2;;

        -v|--verbose)
            VERBOSITY_OPT="-"$(echo $VERBOSITY_OPT | sed s#-##)"v"
            #echo "Setting VERBOSITY_OPT"
            haveFlag=1
            shift ;;

        -h|--help)
            if [[ $haveFlag -eq 0 ]];
            then
                helpFlag=1
            fi
            shift;
            break ;;

        --)
            if [[ $haveFlag -eq 0 ]];
            then
                usageFlag=1
            fi
            shift;
            break;;

    esac
done

if [[ $helpFlag -eq 1 ]];
then
    $PY $CSPP_EDR_HOME/viirs/ql_viirs_edr.py -h
    exit 0
fi
if [[ $usageFlag -eq 1 ]];
then
    $PY $CSPP_EDR_HOME/viirs/ql_viirs_edr.py -h
    exit 0
fi

#echo "GEO_FILES_OPT      = "$GEO_FILES_OPT
#echo "IP_FILES_OPT       = "$IP_FILES_OPT
#echo "PROD_OPT           = "$PROD_OPT
#echo "PLOT_MIN_OPT       = "$PLOT_MIN_OPT
#echo "PLOT_MAX_OPT       = "$PLOT_MAX_OPT
#echo "DPI_OPT            = "$DPI_OPT
#echo "SCALE_OPT          = "$SCALE_OPT
#echo "LAT_0_OPT          = "$LAT_0_OPT
#echo "LON_0_OPT          = "$LON_0_OPT
#echo "STRIDE_OPT         = "$STRIDE_OPT
#echo "SCATTER_PLOT_OPT   = "$SCATTER_PLOT_OPT
#echo "POINTSIZE_OPT      = "$POINTSIZE_OPT
#echo "MAP_RES_OPT        = "$MAP_RES_OPT
#echo "MAP_ANNOTATION_OPT = "$MAP_ANNOTATION_OPT
#echo "OUTPUT_FILE_OPT    = "$OUTPUT_FILE_OPT
#echo "VERBOSITY_OPT        = "$VERBOSITY_OPT



GDB=''
#GDB='gdb --args'
#$GDB $PY $CSPP_RT_HOME/viirs/edr/adl_viirs_edr.py \


#echo "$PY $CSPP_EDR_HOME/viirs/edr/ql_viirs_edr.py \
    #$GEO_FILES_OPT \
    #$IP_FILES_OPT \
    #$PROD_OPT \
    #$PLOT_MIN_OPT \
    #$PLOT_MAX_OPT \
    #$DPI_OPT \
    #$SCALE_OPT \
    #$LAT_0_OPT \
    #$LON_0_OPT \
    #$STRIDE_OPT \
    #$SCATTER_PLOT_OPT \
    #$POINTSIZE_OPT \
    #$MAP_RES_OPT \
    #$MAP_ANNOTATION_OPT \
    #$OUTPUT_FILE_OPT \
    #$VERBOSITY_OPT

#"

#exit 1

$PY $CSPP_EDR_HOME/viirs/ql_viirs_edr.py \
    $GEO_FILES_OPT \
    $IP_FILES_OPT \
    $PROD_OPT \
    $PLOT_MIN_OPT \
    $PLOT_MAX_OPT \
    $DPI_OPT \
    $SCALE_OPT \
    $LAT_0_OPT \
    $LON_0_OPT \
    $STRIDE_OPT \
    $SCATTER_PLOT_OPT \
    $POINTSIZE_OPT \
    $MAP_RES_OPT \
    $MAP_ANNOTATION_OPT \
    $OUTPUT_FILE_OPT \
    $VERBOSITY_OPT


##############################
#         Packaging          #
##############################

#bash $CSPP_RT_HOME/../CSPP_RT_repo/trunk/scripts/edr/CSPP_RT_ViirsEdrMasks_Package.sh  $CSPP_RT_HOME/viirs/edr/viirs_edr.sh ../../sample_data/viirs/edr/input/VIIRS_OPS_unpackTest/HDF5/
