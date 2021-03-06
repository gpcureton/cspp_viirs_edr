#!/bin/bash
# $Id: ql_viirs_edr.sh 1512 2013-07-03 15:27:11Z geoffc $
# Wrapper script for VIIRS EDR quicklooks python script.
#
# Environment settings:
# CSPP_EDR_HOME : the location of the CSPP EDR main directory
#
# This script is dependent on cspp_edr_env.sh having been sourced

if [ -z "$CSPP_EDR_HOME" ]; then
    echo "CSPP_EDR_HOME is not set, but is required for this script to operate."
    exit 9
fi

. ${CSPP_EDR_HOME}/cspp_edr_runtime.sh

# Check number of arguments
if [ "$#" -ne 3 ]; then
  echo "Usage: ql_viirs_edr.sh <PROD> <GMTCO or GITCO file path> <PROD file path>"
  echo "where"
  echo " PROD is one of the following: " 
  echo "    VCM (cloud mask) uses GMTCO* geolocation files" 
  echo "    AOT (Aerosol Optical Thickness) uses GMTCO* geolocation files"
  echo "    SST (Sea Surface Temperature) uses GMTCO* geolocation files"
  echo "    NDVI (Normalized Difference Vegetation Index) uses GITCO* geolocation files"
  echo "    EVI (Enhanced Vegetation Index) uses GITCO* geolocation files"
  echo "    LST (Land Surface Temperature ) uses GMTCO* geolocation files"
  echo " GMTCO/GITCO file path is the full path to the Geolocation file directory"
  echo " PROD file path is the path to the EDR file directory"
  exit 1
fi

PROD=$1
if [[ "$PROD" != VCM && "$PROD" != AOT && "$PROD" != SST && "$PROD" != NDVI && "$PROD" != EVI && "$PROD" != LST ]] ; then
   echo "Input product is not valid :" $PROD
   exit 1
fi

GEO_PATH=$2
if [ ! -d $GEO_PATH ] ; then
  echo "Path to geolocation GMTCO files does not exist :" $GEO_PATH
  exit 1
fi

EDR_PATH=$3
if [ ! -d $EDR_PATH ] ; then
  echo "Path to VIIRS EDR files does not exist :" $EDR_PATH
  exit 1
fi

#Run python command
if [[ "$PROD" == "VCM" ]] ; then
  IPFILE=IICMO
  OUTFILENAME=VIIRS_Cloud_Mask.png
fi

if [[ "$PROD" == "AOT" ]] ; then
  IPFILE=IVAOT
  OUTFILENAME=VIIRS_Aerosol_Optical_Thickness.png
fi

if [[ "$PROD" == "SST" ]] ; then
  PROD=SST_EDR
  IPFILE=VSSTO
  OUTFILENAME=VIIRS_Sea_Surface_Temperature.png
fi

if [[ "$PROD" == "NDVI" ]] ; then
  IPFILE=VIVIO
  OUTFILENAME=VIIRS_Normalized_Difference_Vegetation_Index.png
fi

if [[ "$PROD" == "EVI" ]] ; then
  IPFILE=VIVIO
  OUTFILENAME=VIIRS_Enhanced_Vegetation_Index.png
fi

if [[ "$PROD" == "LST" ]] ; then
  IPFILE=VLSTO
  OUTFILENAME=VIIRS_Land_Surface_Temperature.png
fi

if [[ "$PROD" != "NDVI" && "$PROD" != "EVI" ]] ; then

   $PY $CSPP_EDR_HOME/viirs/ql_viirs_edr.py --geo_file="${GEO_PATH}/GMTCO*" --ip_file="${EDR_PATH}/${IPFILE}*.h5" -p ${PROD} -d 300 --stride=5 -m 'l' -o ${OUTFILENAME}

else

   $PY $CSPP_EDR_HOME/viirs/ql_viirs_edr.py --geo_file="${GEO_PATH}/GITCO*" --ip_file="${EDR_PATH}/${IPFILE}*.h5" -p ${PROD} -d 300 --stride=5 -m 'l' -o ${OUTFILENAME}

fi

if [ $? -eq  0 ] ;  then
   echo "VIIRS EDR quick look script successfully finished"
   exit 0
else 
   echo "VIIRS EDR quick look script failed"
   exit 1
fi

