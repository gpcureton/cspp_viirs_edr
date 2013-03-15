#!/bin/sh
# Wrapper script for VIIRS SDR components from ADL 3.1
#
# Copyright 2011, University of Wisconsin Regents.
# Licensed under the GNU GPLv3.

if [ -z "$CSPP_SDR_HOME" ]; then
  echo "CSPP_SDR_HOME must be set to the path where the CSPP software was installed."
  echo "export CSPP_SDR_HOME=/home/me/SDR_x_x"
  exit 1
fi


usage() {
  cat <<EOF
Usage:
  export CSPP_SDR_HOME=/path/to/here
  \$CSPP_SDR_HOME/viirs/sdr/viirs_sdr.sh [OPTIONS] viirs-rdr-1.h5 [viirs-rdr-2.h5 [viirs-rdr-3.h5 [...]]]

   -h                      Print this message
   -d                      Optional. If specified, intermediate files will always be 
                          retained
   -z                     Optional.  Compress with h5repack zip
   -a                     Optional.  Aggregate using nagg tool
   -l                     Optional.  Inhibit automatic download of ancillary files

  
  One or more VIIRS RDR HDF5 files containing 1 or more  granules is needed. 
  Wildcards are allowed. Files can be specified in any order.

Example:
  \$CSPP_SDR_HOME/viirs/sdr/viirs_sdr.sh  /data/rdr/RVIIRS-RNSCA_npp_d20030125_t06*h5 
EOF

}

. ${CSPP_SDR_HOME}/common/cspp_common.sh

if [ -z "$1" ]; then
    usage
  exit 3
fi

$PY $CSPP_RT_HOME/viirs/adl_viirs_sdr.py -vv "$@"
