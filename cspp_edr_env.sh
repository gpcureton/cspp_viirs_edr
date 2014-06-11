#!/bin/bash
# $Id$
# Environment script for CSPP / ADL 3.1.


test -n "$CSPP_EDR_HOME" || echo "CSPP_EDR_HOME is not set. Please set this environment variable to the install location of CSPP software packages. (When installed, \$CSPP_EDR_HOME/ADL is a directory.)"

test -d "$CSPP_EDR_HOME/common/ADL" || echo "CSPP_EDR_HOME does not appear to be set properly. See cspp_edr_env.sh"

# revision string for this CSPP release, which we set if we have reasonable expectation that environment is correct
test -d "$CSPP_EDR_HOME/common/ADL" && export CSPP_REV="20120215"


#
# user path environment settings, making it easy to invoke wrapper scripts
#

#export DCONFIG=${CSPP_EDR_HOME}/common/cspp_cfg/cfg

export PATH=${CSPP_EDR_HOME}/common:$PATH
export PATH=${CSPP_EDR_HOME}/common/ShellB3/bin:$PATH
export PATH=${CSPP_EDR_HOME}/viirs:$PATH


