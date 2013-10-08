#!/bin/sh
# $Id$
# Environment script for CSPP / ADL 3.1.


test -n "$CSPP_EDR_HOME" || echo "CSPP_EDR_HOME is not set. Please set this environment variable to the install location of CSPP software packages. (When installed, \$CSPP_EDR_HOME/ADL is a directory.)"

test -d "$CSPP_EDR_HOME/common/ADL" || echo "CSPP_EDR_HOME does not appear to be set properly. See cspp_edr_env.sh"

# revision string for this CSPP release, which we set if we have reasonable expectation that environment is correct
test -d "$CSPP_EDR_HOME/common/ADL" && export CSPP_REV="20120215"


#
# derived CSPP default locations (site installs may move these under some circumstances)
#

export CSPP_RT_HOME=$CSPP_EDR_HOME

# read-write directory into which new ancillary data can be downloaded
export CSPP_EDR_ANC_CACHE_DIR=${CSPP_EDR_HOME}/anc/cache

# read-write directory for initial EDR luts and download EDR luts
export CSPP_EDR_LUTS=${CSPP_EDR_ANC_CACHE_DIR}/luts

# static ancillary data including default algorithm settings
export CSPP_EDR_ANC_HOME=${CSPP_EDR_HOME}/anc/static

# default location of static ancillary tiles, which we use in-place rather than linking into workspace
export CSPP_EDR_ANC_TILE_PATH=${CSPP_EDR_ANC_HOME}

export LD_LIBRARY_PATH=${CSPP_SDR_HOME}/common/local/lib64

#
# user path environment settings, making it easy to invoke wrapper scripts
#

#export DCONFIG=${CSPP_EDR_HOME}/common/cspp_cfg/cfg
#export DCONFIG=${CSPP_EDR_HOME}/common/ADL/cfg
unset DCONFIG

export PATH=${CSPP_EDR_HOME}/common:$PATH
export PATH=${CSPP_EDR_HOME}/common/ShellB3/bin:$PATH
export PATH=${CSPP_EDR_HOME}/viirs:$PATH


