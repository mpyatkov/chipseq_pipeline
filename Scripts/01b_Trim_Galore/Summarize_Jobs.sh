#!/bin/bash
set -eu
##################################################################################
# Andy Rampersaud, 03.21.17
#This script would be used to summarize Trim_Galore statistics 
#Way to run script:
#Usage: 
#./Summarize_Jobs.sh
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
echo "SCRIPT_DIR:"
echo ${SCRIPT_DIR}
echo "TIME_LIMIT:"
echo ${TIME_LIMIT}
echo "-----------------------"
echo "End of variable list"
echo "-----------------------"
#---------------------------------------------------------------------------------
#Minor issue with jobs using multiple cores
#Empty job output log files (*.pe* and *.po*) are created
#Remove them if they exist:
echo 'Removing  *.pe* files...'
#Need the 2>/dev/null to supress the "No such file or directory"
count=`ls -1 *.pe* 2>/dev/null  | wc -l`
if [ $count != 0 ]
then 
rm *.pe*
fi 
echo 'Removing  *.po* files...'
#Need the 2>/dev/null to supress the "No such file or directory"
count=`ls -1 *.po* 2>/dev/null  | wc -l`
if [ $count != 0 ]
then 
rm *.po*
fi 
#---------------------------------------------------------------------------------
##################################################################################
#---------------------------------------------------------------------------------
#Retrieve the job name for this step:
#Extract the folder name:
DIR_name=`basename ${SCRIPT_DIR}`
#---------------------------------------------------------------------------------
#Do not hard-code /Scripts/
#Instead use: ${SCRIPT_DIR}/Job_Summary
OUTPUT_DIR=${SCRIPT_DIR}/Job_Summary
#---------------------------------------------------------------------------------
###############################
if [ ! -d ${OUTPUT_DIR} ]; then
mkdir -p ${OUTPUT_DIR}
else
#Remove dir:
rm -r ${OUTPUT_DIR}
#Make new dir:
mkdir -p ${OUTPUT_DIR}
fi
###############################
OUTPUT_FILE=$OUTPUT_DIR/Trim_Galore_Stats.txt
######################
if [ -f ${OUTPUT_FILE} ]
then 
rm ${OUTPUT_FILE}
else
touch ${OUTPUT_FILE}
fi
######################
cd ${SCRIPT_DIR}
################################################
#A text file (Sample_Labels.txt) is needed to run this script
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
#Print header to output file:
echo Sample_ID $'\t'Description $'\t'FASTQ_Name $'\t'Total_Reads_Processed $'\t'Reads_With_Adapters $'\t'Reads_Passing_Filters $'\t'Total_bp_processed $'\t'Quality_Trimmed $'\t'Total_Written_to_Output $'\t'Total_Sequences_Processed >> ${OUTPUT_FILE}
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
echo $Sample_ID
#Need to cd to sample specific CollectInsertSizeMetrics folder
cd ${Dataset_DIR}/$Sample_ID/fastq/trim_galore_output
#---------------------------------------------------------------------------------
#Use a list of files to deal with either 1 (SE) or 2 (PE) report files
Report_List=*'_trimming_report.txt'
for Report_File in ${Report_List}
do
echo ${Report_File}
#Extract the FASTQ filename:
FASTQ_Name=${Report_File%\_trimming_report.txt}
#echo ${FASTQ_Name}
#Extract statistics from the Report_File
Total_Reads_Processed=$(grep 'Total reads processed:' ${Report_File} | awk '{print $4}')
Reads_With_Adapters=$(grep 'Reads with adapters:' ${Report_File} | awk '{print $4" "$5}')
Reads_Passing_Filters=$(grep 'Reads written (passing filters)' ${Report_File} | awk '{print $5" "$6}')
Total_bp_processed=$(grep 'Total basepairs processed' ${Report_File} | awk '{print $4" "$5}')
Quality_Trimmed=$(grep 'Quality-trimmed' ${Report_File} | awk '{print $2" "$3" "$4}')
Total_Written_to_Output=$(grep 'Total written (filtered)' ${Report_File}  | awk '{print $4" "$5" "$6}')
Total_Sequences_Processed=$(grep 'sequences processed in total' ${Report_File} | awk '{print $1}')
#Print to output file
echo "${Sample_ID}\t${Description}t${FASTQ_Name}\t${Total_Reads_Processed}\t${Reads_With_Adapters}\t${Reads_Passing_Filters}\t${Total_bp_processed}\t${Quality_Trimmed}\t${Total_Written_to_Output}\t${Total_Sequences_Processed}" >> ${OUTPUT_FILE}
#End of loop over Report_List
done
#End of loop over Sample_Labels.temp
done < Sample_Labels.temp
#---------------------------------------------------------------------------------
cd ${SCRIPT_DIR}
#Remove the temp file:
rm Sample_Labels.temp
echo '#-------------------------------------------------------------------------'
echo 'Check out '${OUTPUT_FILE}
echo '#-------------------------------------------------------------------------'
##################################################################################
