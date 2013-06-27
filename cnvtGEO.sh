#!/bin/tcsh -f

# Setup command arguements
set inFile = $1
set replaceFile = "" 
if($#argv == 2) then
   set replaceFile = "Yes" 
endif

if(($#argv == 0) || (! -e "${inFile}")) then
   echo ""
   echo "Converts GEO products from Mx6.2 format to Mx6.3 and vice versa."
   echo "Converting to 6.3 inserts a zero filled ScanQuality2 Field"
   echo "Converting to 6.2 removes ScanQuality2 Field"
   echo ""
   echo "Usage: $0 fileName [replaceFile]"
   echo "fileName: Must be a fully qualified path and filename to be converted"
   echo "replaceFile: Any non blank value will cause the input file to be replaced "
   echo "             by the converted file."
   echo ""
   exit
endif

# Setup variables
set inFileSize = `ls -l ${inFile} | awk '{print $5}'`
set Mx63Suffix = _6_3
set Mx62Suffix = _6_2
set suffix = "unknown"

#Check the in file size against known sizes to determine direction and cut into parts
if("${inFileSize}" == "127970224") then # 6.2,VIIRS_DNB_GEO
   set suffix = ${Mx63Suffix}
   split -b 124849072 ${inFile} ${Mx63Suffix}
else if("${inFileSize}" == "324406184") then # 6.2,VIIRS_IMG_GEO, TC, RGEO,
   set suffix = ${Mx63Suffix}
   split -b 314575736 ${inFile} ${Mx63Suffix}
else if("${inFileSize}" == "81103784") then # 6.2,VIIRS_MOD_GEO, TC, RGEO,
   set suffix = ${Mx63Suffix}
   split -b 78646184 ${inFile} ${Mx63Suffix}
else if("${inFileSize}" == "127331240") then # 6.2,VIIRS_MOD_UNAGG, RGEO
   set suffix = ${Mx63Suffix}
   split -b 123472808 ${inFile} ${Mx63Suffix}
else if("${inFileSize}" == "127970272") then # 6.3,VIIRS_DNB_GEO
   set suffix = ${Mx62Suffix}
   split -b 124849120 ${inFile} ${Mx62Suffix}
   split -b 124849072 ${Mx62Suffix}aa 1${Mx62Suffix}
else if("${inFileSize}" == "324406232") then # 6.3,VIIRS_IMG_GEO...
   set suffix = ${Mx62Suffix}
   split -b 314575832 ${inFile} ${Mx62Suffix}
   split -b 314575784 ${Mx62Suffix}aa 1${Mx62Suffix}
else if("${inFileSize}" == "81103832") then # 6.3,VIIRS_MOD_GEO...
   set suffix = ${Mx62Suffix}
   split -b  78646232 ${inFile} ${Mx62Suffix}
   split -b  78646184 ${Mx62Suffix}aa 1${Mx62Suffix}
else if("${inFileSize}" == "127331288") then # 6.3,VIIRS_MOD_UNAGG_GEO...
   set suffix = ${Mx62Suffix}
   split -b  123472856 ${inFile} ${Mx62Suffix}
   split -b  123472808 ${Mx62Suffix}aa 1${Mx62Suffix}
endif


# Join parts back into a file and replace file if need be
set outFile = "${inFile}${suffix}"
if ("${suffix}" == "${Mx63Suffix}") then
   #Create zero filled binary.
   dd if=/dev/zero of=output.bin bs=48 count=1 

   cat ${suffix}aa   > ${outFile}
   cat output.bin    >> ${outFile} 
   cat ${suffix}ab   >> ${outFile} 
   if("${replaceFile}" != "") then
      mv ${outFile} ${inFile}
   endif
else if ("${suffix}" == "${Mx62Suffix}") then
   cat 1${suffix}aa   > ${outFile} 
   cat ${suffix}ab   >> ${outFile}
   if("${replaceFile}" != "") then
      mv ${outFile} ${inFile}
   endif
else 
   echo "Skipping unknown conversion"
   exit
endif

