#!/bin/bash

#
# Move ADL blob and asc files which match the supplied pattern, to 
# the specified directory.
# i.e.:
#
# move_adl_blobs "./*.VIIRS-M[2]-*" ~/blobs/SDR/
#
function move_adl_blobs() {
    if [ $# != 2 ] ;
    then 
        return 1 ;
    else
        blobFiles=$1
        destDir=$2
        curDir=$(pwd)
        echo "Current directory is: "$curDir
        echo "Data directory is: "$dataDir$'\n'

        for files in $(ls $blobFiles --color=none);
        do
            dataDir=$(dirname $files)
            blobFile=$(basename $files)
            urid=$(echo $blobFile | gawk -F\. '{print $1}')
            ascFile=$urid".asc"

            echo $urid"  -->  "$(ls $dataDir/$urid.*)
            mv -vf $(ls $dataDir/$urid.*) $destDir
            
        done

		return 0 ;
    fi
}

#
# Remove ADL blob and asc files which match the supplied pattern.
# i.e.:
#
# remove_adl_blobs "./*.VIIRS-M[2]-*"
#
function remove_adl_blobs() {
    if [ $# != 1 ] ;
    then 
        return 1 ;
    else
        blobFiles=$1
        #destDir=$2
        curDir=$(pwd)
        echo "Current directory is: "$curDir
        echo "Data directory is: "$dataDir$'\n'

        for files in $(ls $blobFiles --color=none);
        do
            dataDir=$(dirname $files)
            blobFile=$(basename $files)
            urid=$(echo $blobFile | gawk -F\. '{print $1}')
            ascFile=$urid".asc"

            echo $urid"  -->  "$(ls $dataDir/$urid.*)
            rm -vf $(ls $dataDir/$urid.*)

        done

		return 0 ;
    fi
}

#
# Creates soft links from ADL blob and asc files which match the supplied pattern, to 
# the current directory.
# i.e.:
#
# link_adl_blobs "~/blobs/SDR/*.VIIRS-M[2]-*"
#
function link_adl_blobs() {
    if [ $# != 1 ] ;
    then 
        echo "Wrong number of arguments" ;
        echo "Number of arguments: "$#", should be 1"
        return 1 ;
    else

        blobFiles=$1
        curDir=$(pwd)
        echo "Current directory is: "$curDir
        echo "Data directory is: "$dataDir$'\n'

        for files in $(ls $blobFiles --color=none);
        do
            dataDir=$(dirname $files)
            blobFile=$(basename $files)
            urid=$(echo $blobFile | gawk -F\. '{print $1}')
            ascFile=$urid".asc"
            ln -s $files $curDir/$blobFile;
            ln -s $dataDir/$ascFile $curDir"/"$ascFile ;
        done

		return 0 ;
    fi
}

#
# Prints out the lines in ADL asc files which match the supplied pattern.
# i.e.:
#
# probe_adl_asc "~/blobs/SDR/*.VIIRS-M[2]-*" N_Granule*
#
function probe_adl_asc() {
    if [ $# != 2 ] ;
    then 
        return 1 ;
    else
        blobFiles=$1
        metadataGlob=$2
        curDir=$(pwd)
        echo "Current directory is: "$curDir
        echo "Data directory is: "$dataDir$'\n'

        for files in $(ls $blobFiles --color=none);
        do
            dataDir=$(dirname $files)
            blobFile=$(basename $files)
            urid=$(echo $blobFile | gawk -F\. '{print $1}')
            ascFile=$urid".asc"

            echo $'\n'$urid"  -->  "$ascFile"  :::  "$blobFile$'\n'
            
            grepVar="";
            grep -E "$metadataGlob" $dataDir/$ascFile  | while read line ; 
            do
                echo $line
            done

        done

		return 0 ;
    fi
}

function probe_adl_asc_alt() {
    if [ $# != 2 ] ;
    then 
        return 1 ;
    else
        blobFiles=$1
        metadataGlob=$2
        curDir=$(pwd)
        echo "Current directory is: "$curDir
        echo "Data directory is: "$dataDir$'\n'

        for files in $(ls $blobFiles --color=none);
        do
            dataDir=$(dirname $files)
            blobFile=$(basename $files)
            urid=$(echo $blobFile | gawk -F\. '{print $1}')
            ascFile=$urid".asc"

            echo $'\n'$urid"  -->  "$ascFile"  :::  "$blobFile$'\n'
            
            grepVar="";
            grep -E "$metadataGlob" $dataDir/$ascFile  | while read line ; 
            do
                grepVar=$grepVar" "$(echo $line)
                echo $grepVar;
            done

        done

		return 0 ;
    fi
}

#
# Prints out the granule ID and various times for a collection of HDF5 files matching the supplied pattern.
# i.e.:
#
# probe_npp_h5_times "~/HDF5/SDR/GMTCO*"
#
function probe_npp_h5_times() {
    if [ $# != 2 ] ;
    then 
        return 1 ;
    else
        h5Files=$1
        shortName=$2
        curDir=$(pwd)
        echo "Current directory is: "$curDir
        echo "Data directory is: "$dataDir$'\n'

        echo "  N_Granule_ID  (Beginning_Date Beginning_Time) (Ending_Date Ending_Time) , <<Agg>> (AggBeginningDate AggBeginningTime)  (AggEndingDate AggEndingTime)    N_Beginning_Time_IET  N_Ending_Time_IET"
        for files in $(ls $h5Files --color=none);
        do
            begDate=$(h5dump -a /Data_Products/"$shortName"/"$shortName"_Gran_0/Beginning_Date $files | grep "(0,0)" | gawk -F\" '{print $2}');
            begTime=$(h5dump -a /Data_Products/"$shortName"/"$shortName"_Gran_0/Beginning_Time $files | grep "(0,0)" | gawk -F\" '{print $2}');
            endDate=$(h5dump -a /Data_Products/"$shortName"/"$shortName"_Gran_0/Ending_Date $files | grep "(0,0)" | gawk -F\" '{print $2}');
            endTime=$(h5dump -a /Data_Products/"$shortName"/"$shortName"_Gran_0/Ending_Time $files | grep "(0,0)" | gawk -F\" '{print $2}');
            aggbegDate=$(h5dump -a /Data_Products/"$shortName"/"$shortName"_Aggr/AggregateBeginningDate $files | grep "(0,0)" | gawk -F\" '{print $2}');
            aggbegTime=$(h5dump -a /Data_Products/"$shortName"/"$shortName"_Aggr/AggregateBeginningTime $files | grep "(0,0)" | gawk -F\" '{print $2}');
            aggendDate=$(h5dump -a /Data_Products/"$shortName"/"$shortName"_Aggr/AggregateEndingDate $files | grep "(0,0)" | gawk -F\" '{print $2}');
            aggendTime=$(h5dump -a /Data_Products/"$shortName"/"$shortName"_Aggr/AggregateEndingTime $files | grep "(0,0)" | gawk -F\" '{print $2}');
            startIETTime=$(h5dump -a /Data_Products/"$shortName"/"$shortName"_Gran_0/N_Beginning_Time_IET $files | grep "(0,0)" | gawk -F\: '{print $2}');
            endIETTime=$(h5dump -a /Data_Products/"$shortName"/"$shortName"_Gran_0/N_Ending_Time_IET $files | grep "(0,0)" | gawk -F\: '{print $2}');
            granID=$(h5dump -a /Data_Products/"$shortName"/"$shortName"_Gran_0/N_Granule_ID $files | grep "(0,0)" | gawk -F\" '{print $2}');
            echo $granID"     "$begDate"    "$begTime"   "$endDate" "$endTime" ,    Agg:     "$aggbegDate"       "$aggbegTime"        "$aggendDate"   "$aggendTime"      "$startIETTime"   "$endIETTime;
        done

		return 0 ;
    fi
}

#
# Prints out the information from /proc/meminfo
#
#function probe_meminfo() {}

function list_CM_luts() {
    if [ $# != 1 ] ;
    then 
        return 1 ;
    else
        AuxPath=$1
    fi

    for names in VIIRS-CM-IP-AC-Int VIIRS-AF-EDR-AC-Int VIIRS-AF-EDR-DQTT-Int;
    do
        newName=${names%-Int};
        #ls $AuxPath/*.*$names*; 
        ls $AuxPath/*.*$newName*; 
    done
}

function link_CM_luts() {
    if [ $# != 1 ] ;
    then 
        return 1 ;
    else
        AuxPath=$1
    fi

    for names in VIIRS-CM-IP-AC-Int VIIRS-AF-EDR-AC-Int VIIRS-AF-EDR-DQTT-Int;
    do
        newName=${names%-Int};
        ln -f -s $(ls $AuxPath/*.*$newName*) template.$names
    done
}

function list_AOT_luts() {
    if [ $# != 1 ] ;
    then 
        return 1 ;
    else
        AuxPath=$1
    fi

    for names in   VIIRS-Aeros-EDR-AC-Int \
                   VIIRS-Aeros-EDR-DQTT-Int \
                   AOT-ANC \
                   VIIRS-AOT-LUT \
                   VIIRS-AOT-Sunglint-LUT \
                   VIIRS-SusMat-EDR-DQTT-Int ;
    do
        newName=${names%-Int};
        ls $AuxPath/*.*$newName*; 
    done
}

function link_AOT_luts() {
    if [ $# != 1 ] ;
    then 
        return 1 ;
    else
        AuxPath=$1
    fi

    for names in   VIIRS-Aeros-EDR-AC-Int \
                   VIIRS-Aeros-EDR-DQTT-Int \
                   AOT-ANC \
                   VIIRS-AOT-LUT \
                   VIIRS-AOT-Sunglint-LUT \
                   VIIRS-SusMat-EDR-DQTT-Int ;
    do
        newName=${names%-Int};
        ln -f -s $(ls $AuxPath/*.*$newName*) template.$names
    done
}

function list_SST_luts() {
    if [ $# != 1 ] ;
    then 
        return 1 ;
    else
        AuxPath=$1
    fi

    for names in VIIRS-SST-EDR-AC-Int VIIRS-SST-EDR-DQTT-Int VIIRS-SST-Coef-LUT;
    do
        newName=${names%-Int};
        ls $AuxPath/*.*$newName*; 
    done
}

function link_SST_luts() {
    if [ $# != 1 ] ;
    then 
        return 1 ;
    else
        AuxPath=$1
    fi

    for names in VIIRS-SST-EDR-AC-Int VIIRS-SST-EDR-DQTT-Int VIIRS-SST-Coef-LUT;
    do
        newName=${names%-Int};
        ln -f -s $(ls $AuxPath/*.*$newName*) template.$names
    done
}

function list_SR_luts() {
    if [ $# != 1 ] ;
    then 
        return 1 ;
    else
        AuxPath=$1
    fi

    for names in VIIRS-SR-IP-AC-Int \
                 VIIRS-SR-AOTValues-LUT \
                 VIIRS-SR-SolZenAngles-LUT \
                 VIIRS-SR-SatZenAngles-LUT \
                 VIIRS-SR-IncScatAngles-LUT \
                 VIIRS-SR-ScatAngDims-LUT \
                 VIIRS-SR-DownTrans-LUT \
                 VIIRS-SR-SphAlb-LUT \
                 VIIRS-SR-AtmReflect-LUT;
    do
        newName=${names%-Int};
        ls $AuxPath/*.*$newName*; 
    done
}

function link_SR_luts() {
    if [ $# != 1 ] ;
    then 
        return 1 ;
    else
        AuxPath=$1
    fi

    for names in VIIRS-SR-IP-AC-Int \
                 VIIRS-SR-AOTValues-LUT \
                 VIIRS-SR-SolZenAngles-LUT \
                 VIIRS-SR-SatZenAngles-LUT \
                 VIIRS-SR-IncScatAngles-LUT \
                 VIIRS-SR-ScatAngDims-LUT \
                 VIIRS-SR-DownTrans-LUT \
                 VIIRS-SR-SphAlb-LUT \
                 VIIRS-SR-AtmReflect-LUT; 
    do
        newName=${names%-Int};
        ln -f -s $(ls $AuxPath/*.*$newName*) template.$names
    done
}

function list_VI_luts() {
    if [ $# != 1 ] ;
    then 
        return 1 ;
    else
        AuxPath=$1
    fi

    for names in VIIRS-VI-EDR-AC-Int VIIRS-VI-EDR-DQTT-Int;
    do
        newName=${names%-Int};
        ls $AuxPath/*.*$newName*; 
    done
}

function link_VI_luts() {
    if [ $# != 1 ] ;
    then 
        return 1 ;
    else
        AuxPath=$1
    fi

    for names in VIIRS-VI-EDR-AC-Int VIIRS-VI-EDR-DQTT-Int;
    do
        newName=${names%-Int};
        ln -f -s $(ls $AuxPath/*.*$newName*) template.$names
    done
}

function list_ST_luts() {
    if [ $# != 1 ] ;
    then 
        return 1 ;
    else
        AuxPath=$1
    fi

    for names in VIIRS-ST-EDR-AC-Int VIIRS-ST-EDR-DQTT-Int;
    do
        newName=${names%-Int};
        ls $AuxPath/*.*$newName*; 
    done
}

function link_ST_luts() {
    if [ $# != 1 ] ;
    then 
        return 1 ;
    else
        AuxPath=$1
    fi

    for names in VIIRS-ST-EDR-AC-Int VIIRS-ST-EDR-DQTT-Int;
    do
        newName=${names%-Int};
        ln -f -s $(ls $AuxPath/*.*$newName*) template.$names
    done
}

function list_LST_luts() {
    if [ $# != 1 ] ;
    then 
        return 1 ;
    else
        AuxPath=$1
    fi

    for names in VIIRS-LST-EDR-AC-Int VIIRS-LST-EDR-DQTT-Int VIIRS-LST-Coef-LUT;
    do
        newName=${names%-Int};
        ls $AuxPath/*.*$newName*; 
    done
}

function link_LST_luts() {
    if [ $# != 1 ] ;
    then 
        return 1 ;
    else
        AuxPath=$1
    fi

    for names in VIIRS-LST-EDR-AC-Int VIIRS-LST-EDR-DQTT-Int VIIRS-LST-Coef-LUT;
    do
        newName=${names%-Int};
        ln -f -s $(ls $AuxPath/*.*$newName*) template.$names
    done
}

