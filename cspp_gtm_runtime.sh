test -n "$CSPP_GTM_HOME" || echo "CSPP_GTM_HOME is not set. Please set this environment variable to the install location of CSPP software packages. (When installed, \$CSPP_GTM_HOME/ADL is a directory.)"



# the adl-common.py module will assign defaults
# these variables should only be set for custom installations
unset CSPP_RT_ANC_CACHE_DIR
unset CSPP_RT_ANC_PATH
unset CSPP_RT_ANC_HOME
unset CSPP_RT_ANC_TILE_PATH
unset CSPP_RT_ANC_HOME
unset CSPP_RT_HOME
unset DCONFIG
unset JPSS_REMOTE_ANC_DIR
export DPE_VER=CSPP_GTM_2_0
#
# derived CSPP default locations (site installs may move these under some circumstances)
#
#

#
# scripting environment settings
#

# python interpreter including numpy, h5py, pytables, scipy; used by CSPP scripts
export PY=${CSPP_GTM_HOME}/common/ShellB3/bin/python

# common modules location used by CSPP scripts
export PYTHONPATH=$CSPP_GTM_HOME/common

#environment cleanups
unset LD_PRELOAD

test -x "$PY" || echo "Python interpreter not available; please source cspp_env.sh"

# Linux execution configuration
export OSTYPE=`uname`

# make the stack size unlimited
ulimit -s unlimited

# Make the core file size unlimited, so that if the algorithm does have a
# segmentation fault, it'll write a core file that can be investigated.
ulimit -c unlimited

# Make the data size unlimited
ulimit -d unlimited

# insurance

export LD_LIBRARY_PATH=${CSPP_GTM_HOME}/common/ADL/lib:${CSPP_GTM_HOME}/common/local/lib64:${CSPP_GTM_HOME}/common/local/lib


