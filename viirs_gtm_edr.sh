#!/bin/sh
# Wrapper script for VIIRS EDR components from ADL 3.1
#
# Copyright 2013, University of Wisconsin Regents.
# Licensed under the GNU GPLv3.

if [ -z "$CSPP_EDR_HOME" ]; then
  echo "CSPP_EDR_HOME must be set to the path where the CSPP software was installed."
  echo "export CSPP_EDR_HOME=/home/me/EDR_x_x"
  exit 1
fi

. ${CSPP_EDR_HOME}/common/cspp_common.sh

if [ -z "$1" ]; then
  $PY $CSPP_RT_HOME/viirs/adl_viirs_gtm_edr.py -h
  exit 3
fi

$PY $CSPP_RT_HOME/viirs/adl_viirs_gtm_edr.py -vv "$@"
