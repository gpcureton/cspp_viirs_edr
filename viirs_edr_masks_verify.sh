#!/bin/bash
#
# This script runs 'hdiff' on portions of self-test output for the ADL 3.1 based VIIRS EDR. 
#
# Copyright 2012, University of Wisconsin Regents.
# Licensed under GNU GPL v2.

test -d "$CSPP_HOME" || echo "ERROR: CSPP_HOME should be set"

echo "Checking cloud Mask output..."
echo "The arguments are "$@

EXPECTED_ARGS=2

if [ $# -ne $EXPECTED_ARGS ]
then
    echo "Usage: `basename $0` dir1 dir2"
    exit -1
fi

sampleDir=$1
workDir=$2

if [ ! -d "$sampleDir" ]; then
    # Will enter here if $DIRECTORY exists, even if it contains spaces
    echo "Directory "$sampleDir" does not exist, aborting"
    exit -1
fi

if [ ! -d "$workDir" ]; then
    # Will enter here if $DIRECTORY exists, even if it contains spaces
    echo "Directory "$workDir" does not exist, aborting"
    exit -1
fi

echo "sampleDir = "$sampleDir
echo "workDir = "$workDir

for fileName in $workDir/IICMO*.h5;
do
    fileName=$(echo $fileName | sed -e 's/\/\//\//')
    stem=$(basename $(echo $fileName | cut -d 'c' -f 1))
    #echo "stem = "$stem
    #echo "workDir = "$workDir
    newFile=$(find $sampleDir -name $stem*h5 -print)
    echo $'\n''####################################'
    echo '####################################'
    echo $'\n''Comparing...'$'\n\t'$fileName$'\n''with...'$'\n\t'$newFile$'\n'
    echo '####################################'$'\n'
    
    for dataSet in \
            /All_Data/VIIRS-CM-IP_All/GranuleAllOcean \
            /All_Data/VIIRS-CM-IP_All/GranuleNoOcean \
            /All_Data/VIIRS-CM-IP_All/QF1_VIIRSCMIP \
            /All_Data/VIIRS-CM-IP_All/QF2_VIIRSCMIP \
            /All_Data/VIIRS-CM-IP_All/QF3_VIIRSCMIP \
            /All_Data/VIIRS-CM-IP_All/QF4_VIIRSCMIP \
            /All_Data/VIIRS-CM-IP_All/QF5_VIIRSCMIP \
            /All_Data/VIIRS-CM-IP_All/QF6_VIIRSCMIP \
            /All_Data/VIIRS-CM-IP_All/ScanAllOcean \
            /All_Data/VIIRS-CM-IP_All/ScanNoOcean \
            /Data_Products/VIIRS-CM-IP/VIIRS-CM-IP_Aggr \
            /Data_Products/VIIRS-CM-IP/VIIRS-CM-IP_Gran_0

    do
        #$CSPP_EDR_HOME/common/ShellB3/bin/h5diff  -v $newFile $fileName $dataSet \
        SSEC/sshfs/rh5b/RH5B/CSPP/common/local/bin/h5diff  -v $newFile $fileName $dataSet \
        || echo ERROR: $dataSet mismatch in $fileName
        echo ""
    done

    echo $'\n''>>>>>> END FILE >>>>>>'$'\n'

done

echo "Done."
