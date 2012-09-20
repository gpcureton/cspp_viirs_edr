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


#
# help text
#

if [ -z "$1" ]; then
  cat <<EOF
Usage:

  export CSPP_HOME=/path/to/CSPP3.1
  source $CSPP_HOME/cspp_env.sh

  mkdir work
  cd work
  viirs_edr_masks.sh /path/to/testdata

  Only input is the directory containing the VIIRS SDR HDF5 files.

  Work-directory will be created if it doesn't exist; defaults to '.'
  Intermediate and output files will be written to the work directory, which defaults to '.'.
  You can override work directory with '-W new_work_dir'.

  HDF5 EDR output files will be written to the work directory.

  Work directory should be cleaned between uses.
Example:
  \$CSPP_HOME/viirs/edr/viirs_edr_masks.sh /path/to/testdata
EOF
  exit 3
fi

# set LD_LIBRARY_PATH as late as possible
# no longer needed as of CSPP3.1beta2
# export LD_LIBRARY_PATH=$CSPP_HOME/common/local/lib:$CSPP_HOME/common/local/lib64:$ADL_HOME/lib:$CSPP_HOME/common/ShellB3/lib

#$PY $CSPP_HOME/viirs/edr/viirsEdrMasks.py -h || oops "VIIRS EDR Masks Controller did not complete without errors."

GDB=''
#GDB='gdb --args'

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
#$GDB $PY $CSPP_HOME/viirs/edr/adl_viirs_edr_masks.py --inputDirectory=$1 --sdr_endianness=little -w ./  --skip_sdr_unpack -vvv --debug

# Skip the ancillary generation...
#$PY $CSPP_HOME/viirs/edr/adl_viirs_edr_masks.py --inputDirectory=$1 --sdr_endianness=little -w ./  --skip_ancillary -vvv --debug

# Skip the sdr unpacking AND ancillary generation...
$GDB $PY $CSPP_HOME/viirs/edr/adl_viirs_edr_masks.py --inputDirectory=$1 --sdr_endianness=little -w ./  --skip_sdr_unpack --skip_ancillary -vvv  --debug

# Run the whole thing...
#$PY $CSPP_HOME/viirs/edr/adl_viirs_edr_masks.py --inputDirectory=$1 --sdr_endianness=little -w ./  -vvv  #--debug


#$PY $CSPP_HOME/viirs/edr/NCEPtoBlob.py -x $CSPP_HOME/ADL/xml/ANC/NCEP_ANC_Int.xml -g $CSPP_HOME/cache/2012_02_17_048/gdas1.pgrb00.1p0deg.20120217_18_000.grib2 -e little -o $CSPP_HOME/cache/2012_02_17_048/gdas1.pgrb00.1p0deg.20120217_18_000.grib2_blob.le
#$PY $CSPP_HOME/viirs/edr/NCEPtoBlob.py -x $CSPP_HOME/ADL/xml/ANC/NCEP_ANC_Int.xml -g $CSPP_HOME/cache/2012_09_05_249/gdas1.pgrb00.1p0deg.20120905_00_000.grib2 -e little -o $CSPP_HOME/cache/2012_09_05_249/gdas1.pgrb00.1p0deg.20120905_00_000.grib2_blob.le
