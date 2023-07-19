#!/bin/bash
set -eu

##################################################################################
# Andy Rampersaud, 03.29.17
##This script would be used to run macs2 in parallel for different BAM files
#Way to run script:
#Usage: ./macs2_Parallel.sh
#Example: 
#./macs2_Parallel.sh 

#Remove *.o* files from previous jobs
#Need the 2>/dev/null to supress the "No such file or directory"
count=`ls -1 *.o* 2>/dev/null  | wc -l`
if [ $count != 0 ]
then 
    rm *.o*
fi 

#Source the setup file to initialize variables
source ../00_Setup_Pipeline/01_Pipeline_Setup.sh

#Create a job name that's a function of the folder name:
#Extract the folder name:
DIR_name=`basename ${SCRIPT_DIR}`
Step_Num=$(echo ${DIR_name} | cut -d'_' -f 1)
Job_Name=$(echo 'Step_'${Step_Num}'_'${Dataset_Label})

#Check that each variable prints a value to the terminal:

echo "Dataset_DIR: ${Dataset_DIR}"
echo "Sample_Labels_DIR: ${Sample_Labels_DIR}"
echo "DHS_DATA: ${DHS_DATA}"
echo "SCRIPT_DIR: ${SCRIPT_DIR}"
echo "TIME_LIMIT: ${TIME_LIMIT}"

#A text file (Sample_Labels.txt) is needed to run this script
SCRIPT_DIR=$(pwd)
cp $Sample_Labels_DIR/Sample_Labels.txt $SCRIPT_DIR

#This will guarantee that only tabs exists between columns:
awk '{print $1"\t"$2"\t"$3}' Sample_Labels.txt > Sample_Labels.temp1
mv Sample_Labels.temp1 Sample_Labels.txt

echo "Start of qsub commands:"
#Text file has a header line to ignore:
tail -n +2 Sample_Labels.txt > Sample_Labels.temp
#Use a while loop to run jobs
while IFS=$'\t' read -r -a myArray
do
    Sample_DIR=${myArray[0]}
    Sample_ID=${myArray[1]}
    Description=${myArray[2]}
    (set -x; qsub -N ${Job_Name}'_'${Sample_ID} -P wax-es -l h_rt=${TIME_LIMIT} -q !linga macs2.qsub ${Sample_ID} ${Dataset_DIR} ${Sample_Labels_DIR} ${DHS_DATA} ${SCRIPT_DIR})
done < Sample_Labels.temp

rm Sample_Labels.temp
echo "End of qsub commands"

