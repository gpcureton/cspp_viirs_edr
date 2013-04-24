#!/bin/sh
# Wrapper script for VIIRS SDR components from ADL 3.1
#
# Copyright 2013, University of Wisconsin Regents.
# Licensed under the GNU GPLv3.

if [ -z "$CSPP_EDR_HOME" ]; then
  echo "CSPP_EDR_HOME must be set to the path where the CSPP software was installed."
  echo "export CSPP_EDR_HOME=/home/me/EDR_x_x"
  exit 1
fi


usage() {
  cat <<EOF
Usage:
  export CSPP_EDR_HOME=/path/to/here
  \$CSPP_EDR_HOME/viirs/sdr/viirs_gtm_edr.sh [OPTIONS] viirs-rdr-1.h5 [viirs-rdr-2.h5 [viirs-rdr-3.h5 [...]]]

   -h                      Print this message
   -d                      Optional. If specified, intermediate files will always be 
                          retained
   -z                     Optional.  Compress with h5repack
   -a                     Optional.  Aggregate using nagg tool
   -l                     Optional.  Inhibit automatic download of ancillary files

  
  One or more VIIRS SDR HDF5 files containing 1 or more  granules is needed.
  Wildcards are allowed. Files can be specified in any order.

Example:
  \$CSPP_EDR_HOME/viirs/sdr/viirs_gtm_edr.sh /data/sdr/SV*h5
EOF

}

. ${CSPP_EDR_HOME}/common/cspp_common.sh

if [ -z "$1" ]; then
    usage
  exit 3
fi

$PY $CSPP_RT_HOME/viirs/adl_viirs_gtm_edr.py -vv "$@"
