#!/bin/bash
set -eu

##################################################################################
#Andy Rampersaud, 03.29.17
#This script would be used to run FASTQC.pbs in parallel for different FASTQ files
#Way to run script:
#Usage: ./FASTQC.sh
#Example: 
#./FASTQC.sh 
##################################################################################
#Remove *.o files from previous jobs
#Need the 2>/dev/null to supress the "No such file or directory"
count=`ls -1 *.o* 2>/dev/null  | wc -l`
if [ $count != 0 ]
then 
rm *.o*
fi 
##################################################################################
#---------------------------------------------------------------------------------
#Source the setup file to initialize variables
source ../00_Setup_Pipeline/01_Pipeline_Setup.sh
#---------------------------------------------------------------------------------
#Create a job name that's a function of the folder name:
#Extract the folder name:
DIR_name=`basename ${SCRIPT_DIR}`
#Split by "_" grab the 1st part (use cut command)
Step_Num=$(echo ${DIR_name} | cut -d'_' -f 1)
#Create the Job_Name:
Job_Name=$(echo 'Step_'${Step_Num}'_'${Dataset_Label})
#Use ${Job_Name} for qsub -N parameter
#--------------------------------------------------------------------------------
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
#Output dir to store text files:
#Parallel jobs will save output to the same folder
#Need to provide the full path (going to be using that path in the scratch dir)
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
echo "Start of qsub commands:"
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
#Now use arguments in the PBS script call:
	qsub -N ${Job_Name}'_'${Sample_ID} -P wax-es -l h_rt=${TIME_LIMIT} FASTQC.qsub ${Sample_ID} ${Dataset_DIR} ${Sample_Labels_DIR} ${BU_User} ${VM_DIR_FASTQC} ${OUTPUT_DIR}
#Printing the qsub command which submits the job
	echo "qsub -N ${Job_Name}'_'${Sample_ID} -P wax-es -l h_rt=${TIME_LIMIT} FASTQC.qsub ${Sample_ID} ${Dataset_DIR} ${Sample_Labels_DIR} ${BU_User} ${VM_DIR_FASTQC} ${OUTPUT_DIR}"
#---------------------------
done < Sample_Labels.temp
#Remove the temp file:
rm Sample_Labels.temp
echo "-----------------------"
echo "End of qsub commands"
echo "-----------------------"
##################################################################################
