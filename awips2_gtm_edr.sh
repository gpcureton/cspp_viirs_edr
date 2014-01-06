#!/bin/bash
# Wrapper script for VIIRS EDR components from ADL 3.1
#
# Copyright 2013, University of Wisconsin Regents.
# Licensed under the GNU GPLv3.

if [ -z "$CSPP_GTM_HOME" ]; then
    echo "CSPP_GTM_HOME is not set, but is required for this script to operate."
    exit 9
fi


. ${CSPP_GTM_HOME}/cspp_gtm_runtime.sh

if [ -z "$1" ]; then
  $PY $CSPP_GTM_HOME/viirs/awips2_gtm_edr.py -h
  exit 3
fi

$PY $CSPP_GTM_HOME/viirs/awips2_gtm_edr.py -v "$@"
