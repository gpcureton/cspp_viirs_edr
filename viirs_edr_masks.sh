#!/bin/env bash
# $Id$
# Wrapper environment script for VIIRS EDR components from ADL 3.1
#
# Environment settings:
# CSPP_HOME : the location of the CSPP directory
#
# Copyright 2011-2012, University of Wisconsin Regents.
# Licensed under the GNU GPLv3.

# This script is dependent on cspp_env.sh having been sourced
if [ -z "$CSPP_HOME" ]; then
    echo "CSPP_HOME is not set, but is required for CSPP to operate."
    exit 9
fi

# load in commonly-used routines and derived environment
# sources cspp_env.sh if needed
source $CSPP_HOME/common/cspp_common.sh


# set up CSPP_ANC_PATH to find VIIRS default configuration
# ancillary tile directory is directly referenced in input XML files
# dynamic ancillary cache is handled by adl_anc_retrieval.py (this last script)
# is currently used for SDR controllers for the moment. VIIRS EDR controllers
# are supplied ancillary data by running bash ancillary scripts directly.
if [ -z "$CSPP_ANC_PATH" ]; then
    export CSPP_ANC_PATH=$CSPP_HOME/viirs/edr:$CSPP_ANC_HOME/ADL/data/repositories/cache
else
    echo "INFO: CSPP_ANC_PATH changed by user to $CSPP_ANC_PATH"
fi

# test that we are reasonably sure we have what we need installed
test -f "$CSPP_HOME/viirs/edr/adl_viirs_edr_masks.py" \
    || oops "$CSPP_HOME/viirs/edr/adl_viirs_edr_masks.py not found"
test -x "$CSPP_HOME/ADL/bin/ProEdrViirsMasksController.exe" \
    || oops "$CSPP_HOME/ADL/bin/ProEdrViirsMasksController.exe not found"
test -x "$PY" \
    || oops "Common CSPP python interpreter $PY not found"
test -w "$CSPP_ANC_CACHE_DIR" \
    || warn "CSPP_ANC_CACHE_DIR is not writable"

help_usage() {
  cat <<EOF

Run the ADL VIIRS EDR.

Usage: 
    export CSPP_HOME=/path/to/CSPP/dir
    \$CSPP_HOME/viirs/edr/viirs_edr_masks.sh [mandatory args] [options]


Options:

  --version             Show program's version number and exit

  -h, --help            Show this help message and exit

  Mandatory Arguments:
    At a minimum these arguments must be specified

    -i INPUTFILES, --input_files=INPUTFILES
                        The base directory where input are placed
    
  Extra Options:
    These options may be used to customize behaviour of this program.

    -w WORK_DIR, --workDirectory=WORK_DIR
                        The directory which all activity will occur in,
                        defaults to the current directory.

    --skip_sdr_unpack   Skip the unpacking of the VIIRS SDR HDF5 files.

    --skip_aux_linking  Skip the the linking to the auxillary files.
    
    --skip_ancillary    Skip the retrieval and granulation of ancillary data.

    --skip_algorithm    Skip running the VIIRS Masks algorithm.
    
    --debug             Enable debug mode on ADL and avoid cleaning workspace
    
    --sdr_endianness=SDR_ENDIANNESS
                        The input VIIRS SDR endianness.
                        Possible values are...
                        'little', 'big'

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
    export CSPP_HOME=/path/to/CSPP/dir
    \$CSPP_HOME/viirs/edr/viirs_edr_masks.sh [mandatory args] [options]

  -h, --help            Show the mandatory args and options and exit.

EOF

}

#
# Gather the various command line options...
#

INPUT_FILES_OPT=
WORK_DIR_OPT=
SKIP_SDR_UNPACK_OPT=
SKIP_ANCILLARY_OPT=
SKIP_ALGORITHM_OPT=
ANC_ENDIANNESS_OPT=
SDR_ENDIANNESS_OPT=
DEBUG_OPT=
VERBOSITY_OPT=

#echo $@

OPTS=`getopt -o "i:w:dvh" -l "input_files:,work_directory:,anc_endianness:,sdr_endianness:,skip_sdr_unpack,skip_aux_linking,skip_ancillary,skip_algorithm,debug,verbose,help" -- "$@"`

# If returncode from getopt != 0, exit with error.
if [ $? != 0 ]
then
    echo "There was an error with getopt, aborting.."
    exit 1
fi

# A little magic
eval set -- "$OPTS"

# Now go through all the options
while true ;
do
    #echo "Checking opts"
    case "$1" in
        -i|--input_files)
            INPUT_FILES_OPT="--input_files=$2"
            shift 2;;

        -w|--work_directory)
            WORK_DIR_OPT="--work_directory=$2"
            shift 2;;

        --sdr_endianness)
            SDR_ENDIANNESS_OPT="--sdr_endianness=$2"
            shift 2;;

        --anc_endianness)
            ANC_ENDIANNESS_OPT="--anc_endianness=$2"
            shift 2;;

        --skip_sdr_unpack)
            SKIP_SDR_UNPACK_OPT="--skip_sdr_unpack"
            shift ;;

        --skip_aux_linking)
            SKIP_AUX_LINKING_OPT="--skip_aux_linking"
            shift ;;

        --skip_ancillary)
            SKIP_ANCILLARY_OPT="--skip_ancillary"
            shift ;;

        --skip_algorithm)
            SKIP_ALGORITHM_OPT="--skip_algorithm"
            shift ;;

        -d|--debug)
            DEBUG_OPT="--debug"
            shift ;;

        -v|--verbose)
            VERBOSITY_OPT="-"$(echo $VERBOSITY_OPT | sed s#-##)"v"
            shift ;;

        -h|--help)
            help_usage
            shift;
            break ;;

        --) 
            usage
            shift;
            break;;
    esac
done


echo "INPUT_FILES_OPT      = "$INPUT_FILES_OPT
echo "WORK_DIR_OPT         = "$WORK_DIR_OPT
echo "SKIP_SDR_UNPACK_OPT  = "$SKIP_SDR_UNPACK_OPT
echo "SKIP_AUX_LINKING_OPT = "$SKIP_AUX_LINKING_OPT
echo "SKIP_ANCILLARY_OPT   = "$SKIP_ANCILLARY_OPT
echo "SKIP_ALGORITHM_OPT   = "$SKIP_ALGORITHM_OPT
echo "SDR_ENDIANNESS_OPT   = "$SDR_ENDIANNESS_OPT
echo "ANC_ENDIANNESS_OPT   = "$ANC_ENDIANNESS_OPT
echo "DEBUG_OPT            = "$DEBUG_OPT
echo "VERBOSITY_OPT        = "$VERBOSITY_OPT


GDB=''
#GDB='gdb --args'
#$GDB $PY $CSPP_HOME/viirs/edr/adl_viirs_edr_masks.py \


#echo "$PY $CSPP_HOME/viirs/edr/adl_viirs_edr_masks.py \
    #$INPUT_FILES_OPT \
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

$PY $CSPP_HOME/viirs/edr/adl_viirs_edr_masks.py \
    $INPUT_FILES_OPT \
    $WORK_DIR_OPT \
    $SKIP_SDR_UNPACK_OPT \
    $SKIP_AUX_LINKING_OPT \
    $SKIP_ANCILLARY_OPT \
    $SKIP_ALGORITHM_OPT \
    $SDR_ENDIANNESS_OPT \
    $ANC_ENDIANNESS_OPT \
    $DEBUG_OPT \
    $VERBOSITY_OPT

##############################
#       No Algorithm         #
##############################

# Run everything except the algorithm...
#$PY $CSPP_HOME/viirs/edr/adl_viirs_edr_masks.py --inputDirectory=$1 --sdr_endianness=little -w ./  --skip_algorithm -vvv --debug

# Skip the SDR unpacking...
#$PY $CSPP_HOME/viirs/edr/adl_viirs_edr_masks.py --inputDirectory=$1 --sdr_endianness=little -w ./  --skip_sdr_unpack --skip_algorithm -vvv --debug
#$GDB $PY $CSPP_HOME/viirs/edr/adl_viirs_edr_masks.py --inputDirectory=$1 --sdr_endianness=little -w ./  --skip_sdr_unpack --skip_algorithm -vvv --debug

# Skip the ancillary generation...
#$PY $CSPP_HOME/viirs/edr/adl_viirs_edr_masks.py --inputDirectory=$1 --sdr_endianness=little -w ./ --skip_ancillary --skip_algorithm -vvv --debug

# Skip the sdr unpacking AND ancillary generation...
#$PY $CSPP_HOME/viirs/edr/adl_viirs_edr_masks.py --inputDirectory=$1 --sdr_endianness=little -w ./ --skip_sdr_unpack --skip_ancillary --skip_algorithm -vvv --debug

##############################
#       With Algorithm       #
##############################

# Skip the SDR unpacking...
#$GDB $PY $CSPP_HOME/viirs/edr/adl_viirs_edr_masks.py --inputDirectory=$1 --sdr_endianness=little -w ./  --skip_sdr_unpack #-vvv --debug

# Skip the ancillary generation...
#$PY $CSPP_HOME/viirs/edr/adl_viirs_edr_masks.py --inputDirectory=$1 --sdr_endianness=little -w ./  --skip_ancillary -vvv --debug

# Skip the sdr unpacking AND ancillary generation...
#$GDB $PY $CSPP_HOME/viirs/edr/adl_viirs_edr_masks.py --inputDirectory=$1 --sdr_endianness=little -w ./  --skip_sdr_unpack --skip_ancillary -vvv  --debug

# Run the whole thing...
#$PY $CSPP_HOME/viirs/edr/adl_viirs_edr_masks.py --inputDirectory=$1 --sdr_endianness=little -w ./ # -vvv  --debug

##############################
#         Packaging          #
##############################

#bash $CSPP_HOME/../CSPP_repo/trunk/scripts/edr/CSPP_ViirsEdrMasks_Package.sh  $CSPP_HOME/viirs/edr/viirs_edr_masks.sh ../../sample_data/viirs/edr/input/VIIRS_OPS_unpackTest/HDF5/
