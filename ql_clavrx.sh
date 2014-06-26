#!/bin/bash
# Wrapper environment script for CLAVR-x Quicklook
#
# Environment settings:
# CSPP_CLAVRX_HOME : the location of the CSPP-CLAVRx directory
#
# Copyright 2011-2012, University of Wisconsin Regents.
# Licensed under the GNU GPLv3.

if [ -z "$CSPP_CLAVRX_HOME" ]; then
  echo "CSPP_CLAVRX_HOME must be set to the path where the CSPP software was installed."
  echo "export CSPP_CLAVRX_HOME=/home/me/CSPP_CLAVRx_v1.0"
  exit 1
fi


#
# help text
#

if [ -z "$1" ]; then
  cat <<EOF
Usage:

    source $CSPP_CLAVRX_HOME/env/set_clavrx_env.sh
    ql_clavrx.sh /path/to/input
  
  Input can be a directory or an individual file. By default, all available
  quicklooks will be created and placed in the current directory.

  You can specify an output directory with the -o flag:

    ql_clavrx.sh -o /path/to/output /path/to/input 
 
  You can create quicklooks for a subset of the available products with the -p flag:
 
    ql_clavrx.sh -p cloud_mask,cloud_type -o /path/to/output /path/to/input
    
    Valid product names are:

      cld_emiss_acha
      cld_height_acha
      cld_press_acha
      cld_opd_dcomp
      cld_reff_dcomp
      cld_temp_acha
      cloud_mask
      cloud_phase
      cloud_probability
      cloud_type

EOF
  exit 3
fi

export LD_LIBRARY_PATH=${CSPP_CLAVRX_HOME}/ShellB3/lib:$LD_LIBRARY_PATH
$PY $CSPP_CLAVRX_HOME/quicklooks/ql_clavrx.py -vv "$@" || echo "Quicklook generator had errors"

