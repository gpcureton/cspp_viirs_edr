#!/bin/bash
# $Id$
# Contains some package specific setup.
#
# Copyright 2011-2013, University of Wisconsin Regents.
# Licensed under the GNU GPLv3.

test -n "$CSPP_EDR_HOME" || echo "CSPP_EDR_HOME is not set. Please set this environment variable to the install location of CSPP software packages. (When installed, \$CSPP_EDR_HOME/ADL is a directory.)"


# the adl-common.py module will assign defaults
# these variables should only be set for custom installations

if [ ! -z "${CSPP_RT_ANC_CACHE_DIR}" ];
then
    echo "Warning overridden default: CSPP_RT_ANC_CACHE_DIR="${CSPP_RT_ANC_CACHE_DIR}
fi

if [ ! -z "${CSPP_RT_ANC_PATH}" ];
then
    echo "Warning overridden default: CSPP_RT_ANC_PATH="${CSPP_RT_ANC_PATH}
fi

if [ ! -z "${CSPP_RT_ANC_HOME}" ];
then
    echo "Warning overridden default: CSPP_RT_ANC_HOME="${CSPP_RT_ANC_HOME}

fi

if [ ! -z "${CSPP_RT_ANC_TILE_PATH}" ];
then
    echo "Warning overridden default: CSPP_RT_ANC_TILE_PATH="${CSPP_RT_ANC_TILE_PATH}

fi

if [ ! -z "${CSPP_RT_ANC_HOME}" ];
then
    echo "Warning overridden default: CSPP_RT_ANC_HOME="${CSPP_RT_ANC_HOME}

fi

if [ ! -z "${CSPP_RT_HOME}" ];
then
    echo "Warning overridden default: CSPP_RT_HOME="${CSPP_RT_HOME}

fi

export DCONFIG=${CSPP_EDR_HOME}/common/cspp_cfg/cfg


if [ ! -z "${DCONFIG}" ];
then
      echo "Warning overridden default: DCONFIG="${DCONFIG}

fi

if [ ! -z "${JPSS_REMOTE_ANC_DIR}" ];
then
     echo "Warning overridden default: JPSS_REMOTE_ANC_DIR="${JPSS_REMOTE_ANC_DIR}
fi


export DPE_VER=CSPP_EDR_2_0

#
# derived CSPP default locations (site installs may move these under some circumstances)
#
#

#
# scripting environment settings
#

# python interpreter including numpy, h5py, pytables, scipy; used by CSPP scripts
export PY=${CSPP_EDR_HOME}/common/ShellB3/bin/python

# common modules location used by CSPP scripts
export PYTHONPATH=$CSPP_EDR_HOME/common:${CSPP_EDR_HOME}/viirs

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

export LD_LIBRARY_PATH=${CSPP_EDR_HOME}/common/ADL/lib:${CSPP_EDR_HOME}/common/local/lib64:${CSPP_EDR_HOME}/common/local/lib


