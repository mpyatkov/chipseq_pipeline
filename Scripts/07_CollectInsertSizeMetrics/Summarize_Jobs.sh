#!/bin/bash
set -eu
##################################################################################
# Andy Rampersaud, 03.21.17
#This script would be used to summarize CollectInsertSizeMetrics statistics 
#Way to run script:
#Usage: 
#./CollectInsertSizeMetrics_Summary.sh
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
##################################################################################
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
OUTPUT_FILE=$OUTPUT_DIR/CollectInsertSizeMetrics_Stats.txt
OUTPUT_FILE_PDF=$OUTPUT_DIR/CollectInsertSizeMetrics_Plots.pdf
######################
if [ -f $OUTPUT_FILE ]
then 
rm $OUTPUT_FILE
else
touch $OUTPUT_FILE
fi
######################
######################
if [ -f $OUTPUT_FILE_PDF ]
then 
rm $OUTPUT_FILE_PDF
else
touch $OUTPUT_FILE_PDF
fi
######################
cd $INPUT_DIR
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
echo Sample_ID $'\t'Description $'\t'MEDIAN_INSERT_SIZE $'\t'MEDIAN_ABSOLUTE_DEVIATION $'\t'MIN_INSERT_SIZE $'\t'MAX_INSERT_SIZE $'\t'MEAN_INSERT_SIZE $'\t'STANDARD_DEVIATION $'\t'READ_PAIRS $'\t'PAIR_ORIENTATION >> $OUTPUT_FILE
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
#Need to cd to sample specific CollectInsertSizeMetrics folder
cd ${Dataset_DIR}/$Sample_ID/fastq/bowtie2/CollectInsertSizeMetrics
#---------------------------------------------------------------------------------
#Copy hist file to OUTPUT_DIR
cp ${Sample_ID}'_hist.pdf' ${OUTPUT_DIR}
#---------------------------------------------------------------------------------
#Extract data from temp1.txt
#Extract lines 7 and 8 from file
sed -n 7,8p $Sample_ID'_metrics' > temp1.txt
MEDIAN_INSERT_SIZE=$(awk '{print $1}' temp1.txt | sed -n 2p)
#echo $MEDIAN_INSERT_SIZE
MEDIAN_ABSOLUTE_DEVIATION=$(awk '{print $3}' temp1.txt | sed -n 2p)
#echo $MEDIAN_ABSOLUTE_DEVIATION
MIN_INSERT_SIZE=$(awk '{print $4}' temp1.txt | sed -n 2p)
#echo $MIN_INSERT_SIZE
MAX_INSERT_SIZE=$(awk '{print $5}' temp1.txt | sed -n 2p)
#echo $MAX_INSERT_SIZE
MEAN_INSERT_SIZE=$(awk '{print $6}' temp1.txt | sed -n 2p)
#echo $MEAN_INSERT_SIZE
STANDARD_DEVIATION=$(awk '{print $7}' temp1.txt | sed -n 2p)
#echo $STANDARD_DEVIATION
READ_PAIRS=$(awk '{print $8}' temp1.txt | sed -n 2p)
#echo $READ_PAIRS
PAIR_ORIENTATION=$(awk '{print $9}' temp1.txt | sed -n 2p)
#echo $PAIR_ORIENTATION
#Remove temp1.txt
rm temp1.txt
#Print to output file
echo $Sample_ID $'\t'$Description $'\t'$MEDIAN_INSERT_SIZE $'\t'$MEDIAN_ABSOLUTE_DEVIATION $'\t'$MIN_INSERT_SIZE $'\t'$MAX_INSERT_SIZE $'\t'$MEAN_INSERT_SIZE $'\t'$STANDARD_DEVIATION $'\t'$READ_PAIRS $'\t'$PAIR_ORIENTATION >> $OUTPUT_FILE
done < Sample_Labels.temp
#---------------------------------------------------------------------------------
#Concatenate the PDF files:
cd ${OUTPUT_DIR}
pdftk *'_hist.pdf' cat output ${OUTPUT_FILE_PDF}
#Remove sample files:
rm *'_hist.pdf'
#---------------------------------------------------------------------------------
cd $INPUT_DIR
#Remove the temp file:
rm Sample_Labels.temp
echo '#-------------------------------------------------------------------------'
echo 'Check out '$OUTPUT_FILE
echo '#-------------------------------------------------------------------------'
##################################################################################
