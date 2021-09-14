#!/bin/bash
set -eu
##################################################################################
#Andy Rampersaud, 11.24.15
#This script would be used to rename folders
#Way to run script:
#Usage: ./Rename_Folders.sh
#Example: 
#./Rename_Folders.sh 
##################################################################################
#Remove *.o files from previous jobs
#Need the 2>/dev/null to supress the "No such file or directory"
count=`ls -1 *.o* 2>/dev/null  | wc -l`
if [ $count != 0 ]
then 
rm *.o*
rm *.e*
fi 
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
echo "-----------------------"
echo "End of variable list"
echo "-----------------------"
#---------------------------------------------------------------------------------
#A text file (Sample_Labels.txt) is needed to run this script
SCRIPT_DIR=$(pwd)
cp $Sample_Labels_DIR/Sample_Labels.txt $SCRIPT_DIR
#---------------------------------------------------------------------------------
#This will guarantee that only tabs exists between columns:
awk '{print $1"\t"$2"\t"$3}' Sample_Labels.txt > Sample_Labels.temp1
mv Sample_Labels.temp1 Sample_Labels.txt
#---------------------------------------------------------------------------------
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
echo "-----------------------"
echo "Start renaming folders:"
echo "-----------------------"
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
#Rename folders:
mv ${Dataset_DIR}/${Sample_DIR} ${Dataset_DIR}/${Sample_ID}
#---------------------------
done < Sample_Labels.temp
#Remove the temp file:
rm Sample_Labels.temp
echo "-----------------------"
echo "End renaming folders"
echo "-----------------------"
##################################################################################
