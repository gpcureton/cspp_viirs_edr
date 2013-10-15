#!/bin/bash
# $Id$
# Environment script for CSPP / ADL 3.1.


test -n "$CSPP_GTM_HOME" || echo "CSPP_GTM_HOME is not set. Please set this environment variable to the install location of CSPP software packages. (When installed, \$CSPP_GTM_HOME/ADL is a directory.)"

test -d "$CSPP_GTM_HOME/common/ADL" || echo "CSPP_GTM_HOME does not appear to be set properly. See cspp_edr_env.sh"

# revision string for this CSPP release, which we set if we have reasonable expectation that environment is correct
test -d "$CSPP_GTM_HOME/common/ADL" && export CSPP_REV="20120215"


#
# derived CSPP default locations (site installs may move these under some circumstances)
#


# read-write directory into which new ancillary data can be downloaded
#export CSPP_RT_ANC_CACHE_DIR=${CSPP_RT_HOME}/anc/cache


# static ancillary data including default algorithm settings
export CSPP_RT_ANC_HOME=${CSPP_RT_HOME}/anc/static
export CSPP_RT_ANC_TILE_PATH=$CSPP_RT_ANC_HOME
# read-write directory for initial EDR luts and download EDR luts
#export CSPP_RT_LUTS=${CSPP_RT_ANC_HOME}




# default location of static ancillary tiles, which we use in-place rather than linking into workspace
#export CSPP_RT_ANC_TILE_PATH=${CSPP_RT_ANC_HOME}
#export CSPP_RT_ANC_PATH=${CSPP_RT_ANC_HOME}

#
# user path environment settings, making it easy to invoke wrapper scripts
#

export PATH=${CSPP_GTM_HOME}/common:$PATH
export PATH=${CSPP_GTM_HOME}/common/ShellB3/bin:$PATH
export PATH=${CSPP_GTM_HOME}/viirs:$PATH


