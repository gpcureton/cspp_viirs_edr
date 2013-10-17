#!/bin/bash
# $Id$
# Environment script for CSPP / ADL 3.1.


test -n "$CSPP_GTM_HOME" || echo "CSPP_GTM_HOME is not set. Please set this environment variable to the install location of CSPP software packages. (When installed, \$CSPP_GTM_HOME/ADL is a directory.)"

test -d "$CSPP_GTM_HOME/common/ADL" || echo "CSPP_GTM_HOME does not appear to be set properly. See cspp_edr_env.sh"

# revision string for this CSPP release, which we set if we have reasonable expectation that environment is correct
test -d "$CSPP_GTM_HOME/common/ADL" && export CSPP_REV="20120215"

#
# user path environment settings, making it easy to invoke wrapper scripts
#

export PATH=${CSPP_GTM_HOME}/common:$PATH
export PATH=${CSPP_GTM_HOME}/common/ShellB3/bin:$PATH
export PATH=${CSPP_GTM_HOME}/viirs:$PATH


