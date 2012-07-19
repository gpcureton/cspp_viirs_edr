#######################
# University of Wisconsin Madison,  Space Science Engineering Center (SSEC)
#file_Date = '$Date: 2011-09-20 15:55:43 -0500 (Tue, 20 Sep 2011) $'
#file_Revision = '$Revision: 160 $'
#file_Author = '$Author: scottm $'
#file_HeadURL = '$HeadURL: https://svn.ssec.wisc.edu/repos/jpss_adl/trunk/add-ons/scripts/adlController.sh $'
#file_Id = '$Id: adlController.sh 160 2011-09-20 20:55:43Z scottm $'
#######################

trace_cmd()
{
    echo "trace"

    #ln -s /mnt/hgfs/CSPP_ANC_HOME/ADL ${CSPP_HOME}/static
    #ln -s  ln -s /mnt/hgfs/CSPP_ANC_CACHE/* ${CSPP_HOME}/cache

    echo "rm -rf ${CSPP_HOME}/cache/*"
    rm -rf ${CSPP_HOME}/cache/*
    echo "strace -o ${TRACE} -f bash $CMD $PARMS"
    strace -o ${TRACE} -f bash $CMD $PARMS

}

get_opens_and_xml_and_exe()
{
    echo "opens"
    rm *.txt

    # Get all the file opens. This catches a lot of junk we don't want, so the output from this 
    # will be further filtered by filter_assorted().
    cat ${TRACE} | grep open \
        | python -c 'import re,sys; txt=sys.stdin.read(); print "\n".join(re.findall( r"""open\("(.*?)",""", txt))' \
        | sort | uniq > $ALLOPENS

    # Get the executables
    cat ${TRACE} | grep "execve(" \
        | python -c 'import re,sys; txt=sys.stdin.read(); print "\n".join(re.findall( r"""execve\("(.*?)",""", txt))' \
        | sort | uniq | grep -v ldd  | grep -v ^/bin/  | grep -v ^/usr/bin/  > exe.txt

    # Get the xml files
    cat $ALLOPENS | grep \\.xml > xml.txt
}

filter_assorted()
{
    # remove system stuff 
    cat $ALLOPENS \
        | grep -v "/dev" \
        | grep -v "/etc" \
        | grep -v "/proc" \
        | grep -v "selinux" \
        | grep -v "/tmp" \
        | grep -v "/usr" \
        | grep -v "linked_data" \
        | grep -v ${WORK_DIR} \
        | grep -v "${CSPP_HOME}/static" \
        | grep -v \.pyc \
        | grep -v "/common/ShellB3"  > filterOutput.txt
    
}

combine_lib_and_xml_and_exe()
{
    # combine opened, executed and miscellaneous lists.

    cat exe.txt  > fileList.txt
    cat xml.txt  >> fileList.txt
    cat filterOutput.txt >> fileList.txt

}

check_if_files_exist()
{

    # check that files realy exist. strace records all opens,  this includes search path attempts
    fileList=$(cat fileList.txt)
    for f in $fileList ;
    do
        if [ -f $f ] ;
        then
            if [ `dirname $f` != "." ] ;
            then
            echo "Good: " $f  
            echo $f >> filesExist.txt
            fi
        fi
    
    done
}

add_python()
{
    # Add the ShellB3 python distro...
    echo "${CSPP_HOME}/common/ShellB3/include/" >> filesExist.txt
    echo "${CSPP_HOME}/common/ShellB3/bin/" >> filesExist.txt
    echo "${CSPP_HOME}/common/ShellB3/sbin/" >> filesExist.txt
    echo "${CSPP_HOME}/common/ShellB3/share/" >> filesExist.txt
    echo "${CSPP_HOME}/common/ShellB3/lib/" >> filesExist.txt
}

# Add various things that might not have been captured by the strace...
add_misc_files()
{
    echo "${CSPP_HOME}/cspp_env.sh" >> filesExist.txt
}

prepare_txt()
{
    # Create a list of the static data directories

    leave=$(basename ${CSPP_HOME})

    cat filesExist.txt \
        | sed s#${CSPP_HOME}#${leave}# \
        | sed s#$leave/ADL/bin/../../common#$leave/common# \
        | sed s#$leave/ADL/bin/../lib#$leave/ADL/lib# \
        | sed s#$leave/ADL/tools/bin/../../../common#$leave/common# \
        | sed s#$leave/ADL/tools/bin/../../lib#$leave/ADL/lib# \
        | sed s#${REPO_HOME}/contrib/scripts#$leave/common# \
        | grep  -v -x -e "/lib64/.*" \
        | grep  -x -v -e "/lib/.*" \
        | grep -v -e "cache" \
        | grep -v -e "GM[OT][DC]O.*\.h5" \
        | grep -v -e "SV[IM][0-9][0-9].*\.h5" \
        | grep -v -e "edr_viirs_masks_NPP.*\.xml" \
        | sort | uniq > $FILELIST

}


create_viirs_edr_masks_tarball()
{

    # Create a tarball containing the parts of ADL required to run the 
    # ProEdrViirsMasks constroller, as well as the supporting python and 
    # bash scripts.
    echo 
    if [ ! -e package ] ;
    then
        mkdir package
    fi

    mv -f *.txt package

    cd ${CSPP_HOME}
    cd ..
    tar -cvz --dereference -T ${WORK_DIR}/package/$FILELIST -f ${CNAME}.tar.gz
    
    cp ${CNAME}.tar.gz $CSPP_PACKAGES
    cp -f ${WORK_DIR}/package/$FILELIST $CSPP_PACKAGES/${FILELIST%.txt}.lst
    chmod gu+rw $CSPP_PACKAGES/${FILELIST%.txt}.lst
}

create_static_data_list()
{
    # Create a list of the static data directories

    leave=$(basename ${CSPP_HOME})

    echo "${CSPP_HOME}/static/asc_templates/"$'\n'\
         "${CSPP_HOME}/static/IGBP/"$'\n'\
         "${CSPP_HOME}/static/LSM/"$'\n'\
         "${CSPP_HOME}/static/NAAPS-ANC-Int/"$'\n'\
         "${CSPP_HOME}/static/NDVI/"$'\n'\
         "${CSPP_HOME}/static/ViirsEdrMasks_Aux/"$'\n'\
         "${CSPP_HOME}/static/VIIRS-MOD-GRC/" \
         | sort | uniq | sed s/^\ // | sed s#${CSPP_HOME}#${leave}# >> $STATIC_FILELIST
}

create_viirs_edr_masks_static_tarball()
{

    # Create a tarball containing the static data required to run the 
    # ProEdrViirsMasks constroller, and the asc and blob file templates used
    # generate custom versions of these files with the correct metadata.
    echo 
    if [ ! -e package ] ;
    then
        mkdir package
    fi

    mv -f $STATIC_FILELIST package

    cd ${CSPP_HOME}
    cd ..
    tar -cvz --dereference -T ${WORK_DIR}/package/packing_list.txt -f ${CNAME}.tar.gz
    
    cp ${CNAME}.tar.gz $CSPP_PACKAGES
    cp -f ${WORK_DIR}/package/packing_list.txt $CSPP_PACKAGES/viirs_edr_masks_packing_list.lst
    chmod gu+rw $CSPP_PACKAGES/viirs_edr_masks_packing_list.lst
}


echo "Run command from CSPP work dir"
if [ -z "$1" ]; then
  cat <<EOF
Usage:
  export DPE_DOMAIN=mydomain
  export DPE_SITE_ID=myorg
  export CSPP_HOME=/path/to/CSPP
#  export CSPP_ANC_PATH=/path/to/ADLData/ADL/data/tiles/Terrain-Eco-ANC-Tile/withMetadata:/path/to/ADLData/ADL/data/repositories/cache
  \$CSPP_HOME/atms/sdr/atms_sdr.sh -W work-directory atms-rdr-1.h5 atms-rdr-2.h5 atms-rdr-3.h5 ...
  
  This command runs the CSPP command using strace.  It then uses the strace ouput to direct tar in creating a relocatible tar file.
  The resulting tar file is copied to jpss-cloud.
  This has been tested with the VIIRS SDR.  It may work unchanged for other packages.
  

Example:
  CSPP_ViirsEdrMasks_Package.sh viirs_edr_masks.sh \$HOME/sample_data/viirs/edr/input/VIIRS_OPS_unpackTest

EOF
  exit 3
fi

main ()
{
    echo "Starting"
    #trace_cmd 
    get_opens_and_xml_and_exe
    filter_assorted
    combine_lib_and_xml_and_exe
    check_if_files_exist
    add_python
    add_misc_files
    prepare_txt
    create_static_data_list
    #create_viirs_edr_masks_tarball
    #create_viirs_edr_masks_static_tarball
}

CMD=$1
PARMS=$2

export WORK_DIR=$(pwd)
echo "WORK_DIR = "$WORK_DIR

export CNAME=$(basename ${CMD%.sh})
echo "CNAME = "$CNAME

export REPO_HOME='/mnt/data/geoffc/CSPP_BETA2_2-15-2012/CSPP_repo'
echo "REPO_HOME = "$REPO_HOME

export TRACE=${WORK_DIR}/${CNAME}.strace.log
echo "TRACE = "$TRACE

echo "Strace command..."$'\n\t'$CMD $PARMS $TRACE

export ALLOPENS=opens.txt
export FILELIST=$CNAME-packing_list.txt
export STATIC_FILELIST=$CNAME-static_packing_list.txt

#export CSPP_PACKAGES="/mnt/hgfs/LATESTEST/packing_lists"
export CSPP_PACKAGES="/mnt/data/geoffc/CSPP_BETA2_2-15-2012/packing_lists"


main 
