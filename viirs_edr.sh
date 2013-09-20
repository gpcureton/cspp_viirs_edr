#!/bin/bash
# $Id$
# Wrapper environment script for VIIRS EDR components from ADL 3.1
#
# Environment settings:
# CSPP_RT_HOME : the location of the CSPP_RT directory
#
# Copyright 2011-2012, University of Wisconsin Regents.
# Licensed under the GNU GPLv3.

# This script is dependent on CSPP_RT_env.sh having been sourced
if [ -z "$CSPP_EDR_HOME" ]; then
    echo "CSPP_EDR_HOME is not set, but is required for CSPP_EDR to operate."
    exit 9
fi

export CSPP_RT_HOME=${CSPP_EDR_HOME}
export CSPP_RT_ANC_CACHE_DIR=${CSPP_EDR_ANC_CACHE_DIR}
# static ancillary data including default algorithm settings
export CSPP_RT_ANC_HOME=${CSPP_EDR_ANC_HOME}
# default location of static ancillary tiles, which we use in-place rather than linking into workspace
export CSPP_RT_ANC_TILE_PATH=${CSPP_EDR_ANC_TILE_PATH}
# load in commonly-used routines and derived environment
# sources CSPP_RT_env.sh if needed
#source $CSPP_RT_HOME/common/cspp_common.sh
. ${CSPP_EDR_HOME}/common/cspp_common.sh


# set up CSPP_RT_ANC_PATH to find VIIRS default configuration
# ancillary tile directory is directly referenced in input XML files
# dynamic ancillary cache is handled by adl_anc_retrieval.py (this last script)
# is currently used for SDR controllers for the moment. VIIRS EDR controllers
# are supplied ancillary data by running bash ancillary scripts directly.
if [ -z "$CSPP_RT_ANC_PATH" ]; then
    export CSPP_RT_ANC_PATH=$CSPP_RT_HOME/anc/static #:$CSPP_RT_ANC_HOME/ADL/data/repositories/cache
else
    echo "INFO: CSPP_RT_ANC_PATH changed by user to $CSPP_RT_ANC_PATH"
fi

# test that we are reasonably sure we have what we need installed
test -f "$CSPP_RT_HOME/viirs/adl_viirs_edr.py" \
    || oops "$CSPP_RT_HOME/viirs/adl_viirs_edr.py not found"
test -x "$CSPP_RT_HOME/common/ADL/bin/ProEdrViirsMasksController.exe" \
    || oops "$CSPP_RT_HOME/common/ADL/bin/ProEdrViirsMasksController.exe not found"
test -x "$CSPP_RT_HOME/common/ADL/bin/ProEdrViirsAerosolController.exe" \
    || oops "$CSPP_RT_HOME/common/ADL/bin/ProEdrViirsAerosolController.exe not found"
test -x "$CSPP_RT_HOME/common/ADL/bin/ProEdrViirsSstController.exe" \
    || oops "$CSPP_RT_HOME/common/ADL/bin/ProEdrViirsSstController.exe not found"
test -x "$CSPP_RT_HOME/common/ADL/bin/ProEdrViirsSurfReflectController.exe" \
    || oops "$CSPP_RT_HOME/common/ADL/bin/ProEdrViirsSurfReflectController.exe not found"
test -x "$CSPP_RT_HOME/common/ADL/bin/ProEdrViirsVI.exe" \
    || oops "$CSPP_RT_HOME/common/ADL/bin/ProEdrViirsVI.exe not found"
test -x "$PY" \
    || oops "Common CSPP_RT python interpreter $PY not found"
test -w "$CSPP_RT_ANC_CACHE_DIR" \
    || warn "CSPP_RT_ANC_CACHE_DIR is not writable"

help_usage() {
  cat <<EOF

Run the ADL VIIRS EDR.

Usage: 
    export CSPP_EDR_HOME=/path/to/CSPP_EDR/dir
    source \$CSPP_EDR_HOME/cspp_edr_env.sh
    \$CSPP_EDR_HOME/viirs/viirs_edr.sh [mandatory args] [options]


Options:

  --version             Show program's version number and exit

  -h, --help            Show this help message and exit

  Mandatory Arguments:
    At a minimum these arguments must be specified

    -i INPUTFILES, --input_files=INPUTFILES
                        The fully qualified path to the input files. 
                        May be a directory or a file glob.

    --alg=ALG           The VIIRS algorithm to run.
                        Possible values are...
                        'VCM','AOT','SST','SRFREF','VI','ATMOS','LAND','OCEAN','ALL'
    
  Extra Options:
    These options may be used to customize behaviour of this program.

    -w WORK_DIR, --work_directory=WORK_DIR
                        The directory which all activity will occur in,
                        defaults to the current directory.

    --skip_sdr_unpack   Skip the unpacking of the VIIRS SDR HDF5 files.

    --skip_aux_linking  Skip the the linking to the auxillary files.
    
    --skip_ancillary    Skip the retrieval and granulation of ancillary data.

    --skip_algorithm    Skip running the VIIRS Masks algorithm.
    
    --debug             Enable debug mode on ADL and avoid cleaning workspace

    --no_chain          Do not run prerequisite algorithms.
    
    --sdr_endianness=SDR_ENDIANNESS
                        The input VIIRS SDR endianness.
                        Possible values are...
                        'little', 'big'. [default: 'little']

    --anc_endianness=ANC_ENDIANNESS
                        The input VIIRS ancillary endianness.
                        Possible values are...
                        'little', 'big'. [default: 'little']
    
    -v, --verbose       each occurrence increases verbosity 1 level from
                        ERROR: -v=WARNING -vv=INFO -vvv=DEBUG"
EOF

}

usage() {
  cat <<EOF

Run the ADL VIIRS EDR.

Usage: 
    export CSPP_EDR_HOME=/path/to/CSPP_EDR/dir
    source \$CSPP_EDR_HOME/cspp_edr_env.sh
    \$CSPP_EDR_HOME/viirs/viirs_edr.sh [mandatory args] [options]

  -h, --help            Show the mandatory args and options and exit.

EOF

}

#
# Gather the various command line options...
#

INPUT_FILES_OPT=
ALG_OPT=
WORK_DIR_OPT=
SKIP_SDR_UNPACK_OPT=
SKIP_ANCILLARY_OPT=
SKIP_ALGORITHM_OPT=
ANC_ENDIANNESS_OPT=
SDR_ENDIANNESS_OPT=
DEBUG_OPT=
CHAIN_OPT=
VERBOSITY_OPT=

#echo $@

OPTS=`getopt -o "i:w:dvh" -l "input_files:,alg:,work_directory:,anc_endianness:,sdr_endianness:,skip_sdr_unpack,skip_aux_linking,skip_ancillary,skip_algorithm,debug,no_chain,verbose,help" -- "$@"`

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
        -i|--input_files)
            INPUT_FILES_OPT="--input_files=$2"
            #echo "Setting INPUT_FILES_OPT"
            haveFlag=1
            shift 2;;

        --alg)
            ALG_OPT="--alg=$2"
            #echo "Setting ALG_OPT"
            haveFlag=1
            shift 2;;

        -w|--work_directory)
            WORK_DIR_OPT="--work_directory=$2"
            #echo "Setting WORK_DIR_OPT"
            haveFlag=1
            shift 2;;

        --sdr_endianness)
            SDR_ENDIANNESS_OPT="--sdr_endianness=$2"
            #echo "Setting SDR_ENDIANNESS_OPT"
            haveFlag=1
            shift 2;;

        --anc_endianness)
            ANC_ENDIANNESS_OPT="--anc_endianness=$2"
            #echo "Setting ANC_ENDIANNESS_OPT"
            haveFlag=1
            shift 2;;

        --skip_sdr_unpack)
            SKIP_SDR_UNPACK_OPT="--skip_sdr_unpack"
            #echo "Setting SKIP_SDR_UNPACK_OPT"
            haveFlag=1
            shift ;;

        --skip_aux_linking)
            SKIP_AUX_LINKING_OPT="--skip_aux_linking"
            #echo "Setting SKIP_AUX_LINKING_OPT"
            haveFlag=1
            shift ;;

        --skip_ancillary)
            SKIP_ANCILLARY_OPT="--skip_ancillary"
            #echo "Setting SKIP_ANCILLARY_OPT"
            haveFlag=1
            shift ;;

        --skip_algorithm)
            SKIP_ALGORITHM_OPT="--skip_algorithm"
            #echo "Setting SKIP_ALGORITHM_OPT"
            haveFlag=1
            shift ;;

        -d|--debug)
            DEBUG_OPT="--debug"
            #echo "Setting DEBUG_OPT"
            haveFlag=1
            shift ;;

        --no_chain)
            CHAIN_OPT="--no_chain"
            #echo "Setting DEBUG_OPT"
            haveFlag=1
            shift ;;

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
    help_usage
    exit 0
fi
if [[ $usageFlag -eq 1 ]];
then
    usage
    exit 0
fi

#echo "INPUT_FILES_OPT      = "$INPUT_FILES_OPT
#echo "ALG_OPT              = "$ALG_OPT
#echo "WORK_DIR_OPT         = "$WORK_DIR_OPT
#echo "SKIP_SDR_UNPACK_OPT  = "$SKIP_SDR_UNPACK_OPT
#echo "SKIP_AUX_LINKING_OPT = "$SKIP_AUX_LINKING_OPT
#echo "SKIP_ANCILLARY_OPT   = "$SKIP_ANCILLARY_OPT
#echo "SKIP_ALGORITHM_OPT   = "$SKIP_ALGORITHM_OPT
#echo "SDR_ENDIANNESS_OPT   = "$SDR_ENDIANNESS_OPT
#echo "ANC_ENDIANNESS_OPT   = "$ANC_ENDIANNESS_OPT
#echo "DEBUG_OPT            = "$DEBUG_OPT
#echo "CHAIN_OPT            = "$CHAIN_OPT
#echo "VERBOSITY_OPT        = "$VERBOSITY_OPT


GDB=''
#GDB='gdb --args'
#$GDB $PY $CSPP_RT_HOME/viirs/edr/adl_viirs_edr.py \


#echo "$PY $CSPP_RT_HOME/viirs/edr/adl_viirs_edr.py \
    #$INPUT_FILES_OPT \
    #$ALG_OPT \
    #$WORK_DIR_OPT \
    #$SKIP_SDR_UNPACK_OPT \
    #$SKIP_AUX_LINKING_OPT \
    #$SKIP_ANCILLARY_OPT \
    #$SKIP_ALGORITHM_OPT \
    #$SDR_ENDIANNESS_OPT \
    #$ANC_ENDIANNESS_OPT \
    #$DEBUG_OPT \
    #$VERBOSITY_OPT
#"

#exit 1

$PY $CSPP_RT_HOME/viirs/adl_viirs_edr.py \
    $INPUT_FILES_OPT \
    $ALG_OPT \
    $WORK_DIR_OPT \
    $SKIP_SDR_UNPACK_OPT \
    $SKIP_AUX_LINKING_OPT \
    $SKIP_ANCILLARY_OPT \
    $SKIP_ALGORITHM_OPT \
    $SDR_ENDIANNESS_OPT \
    $ANC_ENDIANNESS_OPT \
    $DEBUG_OPT \
    $CHAIN_OPT \
    $VERBOSITY_OPT

##############################
#         Packaging          #
##############################

#bash $CSPP_RT_HOME/../CSPP_RT_repo/trunk/scripts/edr/CSPP_RT_ViirsEdrMasks_Package.sh  $CSPP_RT_HOME/viirs/edr/viirs_edr.sh ../../sample_data/viirs/edr/input/VIIRS_OPS_unpackTest/HDF5/
