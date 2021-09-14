#!/bin/bash
set -eu
##################################################################################
# Andy Rampersaud, 07.19.17
#This script would be used to summarize FASTQC jobs
#Way to run script:
#Usage: 
#./FASTQC_Summary.sh
##################################################################################
#---------------------------------------------------------------------------------
#Source the setup file to initialize variables
source ../00_Setup_Pipeline/01_Pipeline_Setup.sh
#---------------------------------------------------------------------------------
#Check that each variable prints a value to the terminal:
echo "-----------------------"
echo "Start of variable list:"
echo "-----------------------"
echo "Dataset_DIR:"
echo ${Dataset_DIR}
echo "Sample_Labels_DIR:"
echo ${Sample_Labels_DIR}
echo "BU_User:"
echo ${BU_User}
echo "VM_DIR_FASTQC:"
echo ${VM_DIR_FASTQC}
echo "TIME_LIMIT:"
echo ${TIME_LIMIT}
echo "-----------------------"
echo "End of variable list"
echo "-----------------------"
#---------------------------------------------------------------------------------
#Retrieve the job name for this step:
Current_DIR=$(pwd)
#Extract the folder name:
DIR_name=`basename ${Current_DIR}`
#---------------------------------------------------------------------------------
#Do not hard-code /Scripts/
#Instead use: ${SCRIPT_DIR}/Job_Summary
OUTPUT_DIR=${SCRIPT_DIR}/Job_Summary
#---------------------------------------------------------------------------------
SCRIPT_DIR=$(pwd)
cp $Sample_Labels_DIR/Sample_Labels.txt $SCRIPT_DIR
################################################
#The text file is formatted like the following:
#----------------------------------------------
#Sample_DIR	Sample_ID	Description
#Sample_Waxman-TP17	G83_M1	Male 8wk-pool 1
#Sample_Waxman-TP18	G83_M2	Male 8wk-pool 2
#Sample_Waxman-TP19	G83_M3	Female 8wk-pool 1
#Sample_Waxman-TP20	G83_M4	Female 8wk-pool 2
#----------------------------------------------
#The 1st column: The Sample_DIR name
#The 2nd column: Waxman Lab Sample_ID 
#The 3rd column: Sample's description 
################################################
#Text file has a header line to ignore:
tail -n +2 Sample_Labels.txt > Sample_Labels.temp
#Use a while loop to run jobs
while IFS=$'\t' read -r -a myArray
do
#---------------------------
##Check that text file is read in properly:
#echo 'Sample_DIR:'
Sample_DIR=${myArray[0]}
#echo 'Sample_ID:'
Sample_ID=${myArray[1]}
#echo $Sample_ID
#echo 'Description:'
Description=${myArray[2]}
#echo $Description
#---------------------------
echo
echo $Sample_ID
echo
#Copy html file to VM_DIR_FASTQC
#Note: We are using the waxman-server mount point on the SCC
set +eu
(set -x; cp ${OUTPUT_DIR}/$Sample_ID'_FASTQC'/*'_fastqc.html' ${VM_DIR_FASTQC})
set -eu
done < Sample_Labels.temp
#Remove the temp file:
rm Sample_Labels.temp
echo '#--------------------------------------------------------------------------'
echo 'Check the URL for the waxmanlabvm HTML files!'
echo 'Load the following link:'
echo 'Note: use full paths (instead of using tilde notation)'
echo 'http://waxmanlabvm.bu.edu/waxmanlab/FASTQC/'
echo '#--------------------------------------------------------------------------'
##################################################################################
