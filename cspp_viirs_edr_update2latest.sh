#!/bin/bash
# $Id$
#
# Script to download the latest version of the CSPP VIIRS EDR software
# and LUTS. Checks if the installed version is the latest version.
# 
#
# Environment settings:
# CSPP_EDR_HOME : the location of the CSPP EDR installation
#
# Written by Geoff Cureton UW/SSEC, January 2014
#
# Copyright 2011-2014, University of Wisconsin Regents.
# Licensed under the GNU GPLv3.


if [ -z "${CSPP_EDR_HOME}" ] ; then
    echo "Please set CSPP_EDR_HOME "
    exit 1
fi

currentDir=$(pwd)
echo "The current directory is ${currentDir}"
CSPP_ROOT=$(dirname ${CSPP_EDR_HOME})
echo "The CSPP root directory is ${CSPP_ROOT}"
CSPP_EDR_LUTS=${CSPP_EDR_HOME}/anc/cache/luts/viirs
echo "The CSPP VIIRS LUTS are in ${CSPP_EDR_LUTS}"
echo ""

#export CSPP_EDR_REMOTE_DIR=http://www.ssec.wisc.edu/~geoffc/CSPP_VIIRS_EDR_packages
export CSPP_EDR_REMOTE_DIR=ftp://ftp.ssec.wisc.edu:/pub/geoffc/CSPP_VIIRS_EDR_packages

################################
#  Updating the cache package  #
################################

checksum_file="CSPP_EDR_CACHE_LATEST_Checklist.md5"
echo "Downloading latest CSPP EDR VIIRS cache checksum ${CSPP_EDR_REMOTE_DIR}/${checksum_file} ..."

lftpget -c "${CSPP_EDR_REMOTE_DIR}/${checksum_file}"
if [ $? -ne 0 ] ; then
    echo ""
    echo "Download of latest CSPP VIIRS CACHE package was not successful."
    exit 1
else
    echo "Download of ${checksum_file} successful"
fi

# Check that the checksum file is actually a file that exists
if [[ ! -e ${checksum_file} && ! -f ${checksum_file} ]]; then
    echo "the file ${checksum_file} does not exist, aborting..."
    exit 1
fi
if [[ ! -f ${checksum_file} ]]; then
    echo "the file ${checksum_file} is not a file, aborting..."
    exit 1
fi

#
# Check that the installed cache package matches the checksum...
#

cache_tarfile="CSPP_EDR_CACHE_LATEST.tar.gz"

cd ${CSPP_ROOT}
md5sum -c ${currentDir}/${checksum_file} --status
md5_success=$?
cd ${currentDir}

if [ ${md5_success} -ne 0 ] ; then
    echo ""
    echo "Current CSPP EDR VIIRS cache installation does not match latest package, reinstalling."

    echo "Downloading latest CSPP EDR VIIRS cache tarball..."
    lftpget -c "${CSPP_EDR_REMOTE_DIR}/${cache_tarfile}"
    if [ $? -ne 0 ] ; then
        echo ""
        echo "Download of latest CSPP EDR VIIRS cache package was not successful."
        exit 1
    else
        echo ""
        echo "Download successful"
    fi

    echo "Removing old CSPP EDR VIIRS cache installation..."
    rm -rvf ${CSPP_EDR_HOME}/anc/cache/luts/viirs
    rm -rvf ${CSPP_EDR_HOME}/anc/cache/NAAPS-ANC-Int

    echo "Installing latest CSPP EDR VIIRS cache package..."
    tar xzvf ${cache_tarfile} -C ${CSPP_ROOT}

    if [ $? -ne 0 ] ; then
        echo ""
        echo "Installation of latest CSPP EDR VIIRS cache package was not successful."
        exit 1
    else
        echo ""
        echo "Installation of latest CSPP EDR VIIRS cache package was successful."
    fi

    # Check that the newly installed cache matches the checksum
    cd ${CSPP_ROOT}
    md5sum -c ${currentDir}/${checksum_file} --status
    md5_success=$?
    cd ${currentDir}

    if [ ${md5_success} -ne 0 ] ; then
        echo ""
        echo "Latest CSPP EDR VIIRS cache installation does not match associated checksum, aborting..."
        exit 1
    else
        echo ""
        echo "Latest CSPP EDR VIIRS cache installation matches associated checksum."
        echo ""
    fi

else
    echo ""
    echo "Current CSPP EDR VIIRS cache installation matches latest package."
    echo ""
fi

###############################
#  Updating the code package  #
###############################

checksum_file="CSPP_EDR_LATEST_Checklist.md5"
echo "Downloading latest CSPP EDR VIIRS code checksum ${CSPP_EDR_REMOTE_DIR}/${checksum_file} ..."

lftpget -c "${CSPP_EDR_REMOTE_DIR}/${checksum_file}"
if [ $? -ne 0 ] ; then
    echo ""
    echo "Download of latest CSPP VIIRS code package was not successful."
    exit 1
else
    echo "Download of ${checksum_file} successful"
fi

# Check that the checksum file is actually a file that exists
if [[ ! -e ${checksum_file} && ! -f ${checksum_file} ]]; then
    echo "the file ${checksum_file} does not exist, aborting..."
    exit 1
fi
if [[ ! -f ${checksum_file} ]]; then
    echo "the file ${checksum_file} is not a file, aborting..."
    exit 1
fi

#
# Check that the installed code package matches the checksum...
#

code_tarfile="CSPP_EDR_LATEST.tar.gz"

cd ${CSPP_ROOT}
md5sum -c ${currentDir}/${checksum_file} --status
md5_success=$?
cd ${currentDir}

if [ ${md5_success} -ne 0 ] ; then
    echo ""
    echo "Current CSPP EDR VIIRS code installation does not match latest package, reinstalling."

    echo "Downloading latest CSPP EDR VIIRS code tarball..."
    lftpget -c "${CSPP_EDR_REMOTE_DIR}/${code_tarfile}"
    if [ $? -ne 0 ] ; then
        echo ""
        echo "Download of latest CSPP EDR VIIRS code package was not successful."
        exit 1
    else
        echo ""
        echo "Download successful"
    fi

    echo "Removing old CSPP EDR VIIRS code installation..."
    rm -rvf $CSPP_EDR_HOME/common
    rm -rvf $CSPP_EDR_HOME/viirs
    rm -rvf $CSPP_EDR_HOME/cspp_edr_*.sh

    echo "Installing latest CSPP EDR VIIRS code package..."
    tar xzvf ${code_tarfile} -C ${CSPP_ROOT}

    if [ $? -ne 0 ] ; then
        echo ""
        echo "Installation of latest CSPP EDR VIIRS code package was not successful."
        exit 1
    else
        echo ""
        echo "Installation of latest CSPP EDR VIIRS code package was successful."
    fi

    # Check that the newly installed code matches the checksum
    cd ${CSPP_ROOT}
    md5sum -c ${currentDir}/${checksum_file} --status
    md5_success=$?
    cd ${currentDir}

    if [ ${md5_success} -ne 0 ] ; then
        echo ""
        echo "Latest CSPP EDR VIIRS code installation does not match associated checksum, aborting..."
    else
        echo ""
        echo "Latest CSPP EDR VIIRS code installation matches associated checksum."
        exit 1
    fi

else
    echo ""
    echo "Current CSPP EDR VIIRS code installation matches latest package."
    echo ""
fi

exit 0
