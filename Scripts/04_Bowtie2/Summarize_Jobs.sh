#!/bin/bash
set -eu
##################################################################################
# Andy Rampersaud, 03.21.17
# Dana Lau, 10.16.14 
# Nicholas J. Lodato, 11.03.15
# Andy Rampersaud, 11.03.15
#This script would be used to summarize Bowtie2 mapping statistics 
#(need the G*_M*.o files from the Bowtie2 job)
#Way to run script:
#Usage: 
#cd scripts directory
#./Bowtie2_Summary.sh
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
echo "Bowtie2Index_DIR:"
echo ${Bowtie2Index_DIR}
echo "-----------------------"
echo "End of variable list"
echo "-----------------------"
#---------------------------------------------------------------------------------
INPUT_DIR=$(pwd)
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
######################
if [ ! -d $OUTPUT_DIR ]
then 
mkdir $OUTPUT_DIR
fi
######################
#OUTPUT_DIR already contains the summary from each sample's job
OUTPUT_FILE=$OUTPUT_DIR/Bowtie2_Stats.txt
######################
if [ -f $OUTPUT_FILE ]
then 
rm $OUTPUT_FILE
else
touch $OUTPUT_FILE
fi
######################
cd $INPUT_DIR
#-----------------------------------------------
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
#-----------------------------------------------
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
echo SAMPLE_ID $'\t'Description $'\t'TOTAL_MAP_RATE $'\t'TOTAL_PE_READ_COUNT $'\t'TOTAL_PAIRS_MADE $'\t'ALIGN_CONCORDANTLY_ONCE_PAIR_COUNT $'\t'PAIRS_ALIGN_CONCORDANTLY_ONE_READ_RATIO $'\t'ALIGN_CONCORDANTLY_ZERO_PAIR_COUNT $'\t'PAIRS_ALIGN_CONCORDANTLY_ZERO_READ_RATIO $'\t'ALIGN_CONCORDANTLY_MORE_ONE_READ_COUNT $'\t'PAIRS_ALIGN_CONCORDANTLY_MORE_ONE_READ_RATIO $'\t'ALIGN_DISCORDANTLY_ONCE_READ_COUNT $'\t'MATES_NOT_DISCORDANT_OR_CONCORDANT $'\t'SINGLETS_NOT_ALIGNED $'\t'SINGLETS_ALIGNED_ONCE $'\t'SINGLETS_ALIGNED_MULTIPLE  >> $OUTPUT_FILE
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
#Get counts from G*_M*.e file
#Use grep extract lines matching a pattern:
TOTAL_PE_READ_COUNT=$(grep 'reads; of these:' *${Sample_ID}.o*  | awk '{print $1}')
TOTAL_PAIRS_MADE=$(grep 'were paired; of these:' *${Sample_ID}.o*  | awk '{print $1}')
ALIGN_CONCORDANTLY_ZERO_PAIR_COUNT=$(grep 'pairs aligned concordantly 0 times; of these:' *${Sample_ID}.o* | awk '{print $1}')
ALIGN_CONCORDANTLY_ONCE_PAIR_COUNT=$(grep 'aligned concordantly exactly 1 time' *${Sample_ID}.o* | awk '{print $1}')
ALIGN_CONCORDANTLY_MORE_ONE_READ_COUNT=$(grep 'aligned concordantly >1 times' *${Sample_ID}.o* | awk '{print $1}')
ALIGN_DISCORDANTLY_ONCE_READ_COUNT=$(grep 'aligned discordantly 1 time' *${Sample_ID}.o* | awk '{print $1}')
MATES_NOT_DISCORDANT_OR_CONCORDANT=$(grep 'mates make up the pairs; of these:' *${Sample_ID}.o* | awk '{print $1}')
SINGLETS_ALIGNED_ONCE=$(grep 'aligned exactly 1 time' *${Sample_ID}.o* | awk '{print $1}')
SINGLETS_ALIGNED_MULTIPLE=$(grep 'aligned >1 times' *${Sample_ID}.o* | awk '{print $1}')
TOTAL_MAP_RATE=$(grep 'overall alignment rate' *${Sample_ID}.o* | awk '{print $1}')
#Calculate percentages
PAIRS_ALIGN_CONCORDANTLY_ZERO_READ_RATIO=$(echo "scale=4;$ALIGN_CONCORDANTLY_ZERO_PAIR_COUNT/$TOTAL_PAIRS_MADE" | bc)
PAIRS_ALIGN_CONCORDANTLY_ONE_READ_RATIO=$(echo "scale=4;$ALIGN_CONCORDANTLY_ONCE_PAIR_COUNT/$TOTAL_PAIRS_MADE" | bc)
PAIRS_ALIGN_CONCORDANTLY_MORE_ONE_READ_RATIO=$(echo "scale=4;$ALIGN_CONCORDANTLY_MORE_ONE_READ_COUNT/$TOTAL_PAIRS_MADE" | bc)
SINGLETS_NOT_ALIGNED=$(echo "scale=4;$MATES_NOT_DISCORDANT_OR_CONCORDANT-$SINGLETS_ALIGNED_ONCE-$SINGLETS_ALIGNED_MULTIPLE" | bc)
#Print to output file
echo ${Sample_ID} $'\t'$Description $'\t'$TOTAL_MAP_RATE $'\t'$TOTAL_PE_READ_COUNT $'\t'$TOTAL_PAIRS_MADE $'\t'$ALIGN_CONCORDANTLY_ONCE_PAIR_COUNT $'\t'$PAIRS_ALIGN_CONCORDANTLY_ONE_READ_RATIO $'\t'$ALIGN_CONCORDANTLY_ZERO_PAIR_COUNT $'\t'$PAIRS_ALIGN_CONCORDANTLY_ZERO_READ_RATIO $'\t'$ALIGN_CONCORDANTLY_MORE_ONE_READ_COUNT $'\t'$PAIRS_ALIGN_CONCORDANTLY_MORE_ONE_READ_RATIO $'\t'$ALIGN_DISCORDANTLY_ONCE_READ_COUNT $'\t'$MATES_NOT_DISCORDANT_OR_CONCORDANT $'\t'$SINGLETS_NOT_ALIGNED $'\t'$SINGLETS_ALIGNED_ONCE $'\t'$SINGLETS_ALIGNED_MULTIPLE  >> $OUTPUT_FILE
#---------------------------
done < Sample_Labels.temp
#Remove the temp file:
rm Sample_Labels.temp
rm Sample_Labels.txt
echo '#-------------------------------------------------------------------------'
echo 'Check out '$OUTPUT_FILE
echo '#-------------------------------------------------------------------------'
##################################################################################
