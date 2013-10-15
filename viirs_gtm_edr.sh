#!/bin/bash
# Wrapper script for VIIRS EDR components from ADL 3.1
#
# Copyright 2013, University of Wisconsin Regents.
# Licensed under the GNU GPLv3.

if [ -z "$CSPP_GTM_HOME" ]; then
  echo "CSPP_GTM_HOME must be set to the path where the CSPP software was installed."
  echo "export CSPP_GTM_HOME=/home/me/GTM_x_x"
  exit 1
fi

. ${CSPP_GTM_HOME}/common/cspp_common.sh

if [ -z "$1" ]; then
  $PY $CSPP_GTM_HOME/viirs/adl_viirs_gtm_edr.py -h
  exit 3
fi

$PY $CSPP_GTM_HOME/viirs/adl_viirs_gtm_edr.py -vv "$@"
