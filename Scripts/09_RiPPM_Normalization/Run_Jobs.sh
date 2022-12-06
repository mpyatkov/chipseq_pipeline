#!/bin/bash
set -eu
##################################################################################
# Andy Rampersaud, 03.29.17
##This script would be used to run Peak_Union_Count in parallel for different BAM files
#Way to run script:
#Usage: ./Run_Jobs.sh
#Example: 
#./Run_Jobs.sh 
##################################################################################
#Remove *.o* files from previous jobs
#Need the 2>/dev/null to supress the "No such file or directory"

count=`ls -1 *.o* 2>/dev/null  | wc -l`
if [ $count != 0 ]
then 
    rm -rf *.o*
fi 
##################################################################################
#---------------------------------------------------------------------------------
#Source the setup file to initialize variables
source ../00_Setup_Pipeline/01_Pipeline_Setup.sh
#--------------------------------------------------------------------------------
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
echo "DHS_DATA:"
echo ${DHS_DATA}
echo "SCRIPT_DIR:"
echo ${SCRIPT_DIR}
echo "TIME_LIMIT:"
echo ${TIME_LIMIT}
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
    qsub -N ${Job_Name}'_'${Sample_ID} -P wax-es -l h_rt=${TIME_LIMIT} Peak_Union_Count.qsub ${Sample_ID} ${Dataset_DIR} ${Sample_Labels_DIR} ${DHS_DATA} ${SCRIPT_DIR}
    #Printing the qsub command which submits the job
    echo "qsub -N ${Job_Name}'_'${Sample_ID} -P wax-es -l h_rt=${TIME_LIMIT} -q !linga Peak_Union_Count.qsub ${Sample_ID} ${Dataset_DIR} ${Sample_Labels_DIR} ${DHS_DATA} ${SCRIPT_DIR}"
    #---------------------------
done < Sample_Labels.temp
#Remove the temp file:
rm Sample_Labels.temp
echo "-----------------------"
echo "End of qsub commands"
echo "-----------------------"
##################################################################################
