#!/bin/bash
set -eu
##################################################################################
# Andy Rampersaud, 03.21.17
#This script would be used to summarize BAM file read counts (after RemoveStraightPeak job)
#Way to run script:
#Usage: 
#./RmStraightPk_Summary.sh.sh
################################################
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
echo "SCRIPT_DIR"
echo ${SCRIPT_DIR}
echo "-----------------------"
echo "End of variable list"
echo "-----------------------"
#---------------------------------------------------------------------------------
SCRIPT_DIR=$(pwd)
DATA_DIR=${Dataset_DIR}
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
###############################
if [ ! -d ${OUTPUT_DIR} ]; then
mkdir -p ${OUTPUT_DIR}
#Don't remove dir (storing files from jobs)
fi
###############################
#OUTPUT_DIR already contains the summary from each sample's job
OUTPUT_FILE=$OUTPUT_DIR/StrgtPksRead_count.txt
OUTPUT_FILE_2=$OUTPUT_DIR/StrgtPksRead_confirm_count.txt
######################
if [ -f $OUTPUT_FILE ]
then 
rm $OUTPUT_FILE
else
touch $OUTPUT_FILE
fi
######################
if [ -f $OUTPUT_FILE_2 ]
then 
rm $OUTPUT_FILE_2
else
touch $OUTPUT_FILE_2
fi
######################
cd $OUTPUT_DIR
################################################
#-----------------------------------------------
#A text file (Sample_Labels.txt) is needed to run this script
SCRIPT_DIR=$(pwd)
cp $Sample_Labels_DIR/Sample_Labels.txt $SCRIPT_DIR
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
#Print header to OUTPUT_FILE:
echo 'SAMPLE_ID' $'\t''MAPPED_READ_COUNT(R1 + R2)' $'\t''MAPPED_READ_COUNT_(After_StrgtPks_removal)(R1 + R2)' $'\t''StrgtPks_READ_RATIO' $'\t''NUMBER_StrgtPks_READS' >> $OUTPUT_FILE
#Print header to OUTPUT_FILE_2:
echo 'SAMPLE_ID' $'\t''Description' $'\t''BAM_FILE_MAPPED_READ_COUNT_(After_StrgtPks_removal)(R1 + R2)' >> $OUTPUT_FILE_2
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
echo ${Sample_ID}
cat ${Sample_ID}'_StrgtPksRead_count.txt' >> $OUTPUT_FILE
#Get BAM file read count
Total_Read_Count=$(head -1 $DATA_DIR/${Sample_ID}/fastq/bowtie2/StrgtPksRm/'flagstat_'${Sample_ID}'_StrgtPksRm.txt' | awk '{print $1}')
echo 'Printing to output file'
echo ${Sample_ID} $'\t'$Description $'\t'$Total_Read_Count >> $OUTPUT_FILE_2
done < Sample_Labels.temp
#Remove the temp file:
rm Sample_Labels.temp
rm Sample_Labels.txt
################################################
echo '#-------------------------------------------------------------------------'
echo 'Check out '$OUTPUT_FILE
echo 'Check out '$OUTPUT_FILE_2
echo '#-------------------------------------------------------------------------'
##################################################################################
