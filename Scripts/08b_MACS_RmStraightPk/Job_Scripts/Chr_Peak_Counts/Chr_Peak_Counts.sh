#!/bin/bash
set -eu
#################################################################################
# Andy Rampersaud, 01.11.2016
#This script would be used to call the Chr_Peak_Counts.R
#################################################################################
#--------------------------------------------------------------------------------
#This script would be placed in the same dir containing the following:
#Chr_Peak_Counts.R
#*.bed file(s)
#--------------------------------------------------------------------------------
#Way to run script:
#Usage: ./Chr_Peak_Counts.sh
#Example: 
#./Chr_Peak_Counts.sh
#################################################################################
INPUT_DIR=$(pwd)
cd $INPUT_DIR
file_list=*.bed
for file in ${file_list}
do
echo ${file}
Rscript Chr_Peak_Counts.R ${file}
done
#--------------------------------------------------------------------------------
OUTPUT_FILE_PDF=${INPUT_DIR}/'Chr_Peak_Counts.pdf'
######################
if [ -f ${OUTPUT_FILE_PDF} ]
then 
rm ${OUTPUT_FILE_PDF}
else
touch ${OUTPUT_FILE_PDF}
fi
######################
echo "-----------------------"
echo 'Create montage of PNG files'
montage -geometry 500x500 *.png ${OUTPUT_FILE_PDF}
echo "-----------------------"
echo "Removing individual *.png files"
rm *.png
echo "-----------------------"
echo 'Check out '${OUTPUT_FILE_PDF}'!'
echo "-----------------------"
#################################################################################
