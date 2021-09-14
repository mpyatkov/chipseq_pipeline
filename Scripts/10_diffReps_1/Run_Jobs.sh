#!/bin/bash
set -eu
##################################################################################
#Andy Rampersaud, 09.26.18
#This script would be used to run diffReps.qsub 
#Way to run script:
#Usage: ./diffReps.sh
#Example: 
#./diffReps.sh 
##################################################################################
#Remove *.o* files from previous jobs
#Need the 2>/dev/null to supress the "No such file or directory"
count=`ls -1 *.o* 2>/dev/null  | wc -l`
if [ $count != 0 ]
then 
rm *.o*
#Remove previous diffReps output folders
rm -r Output_*
#Remove files from using multiple cores
rm *.po*
fi 
##################################################################################
#---------------------------------------------------------------------------------
#Source the setup file to initialize variables
source ../00_Setup_Pipeline/01_Pipeline_Setup.sh
#Source job-specific variables:
source setup_diffReps.sh
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
echo "Dataset_Label:"
echo ${Dataset_Label}
echo "Control_Samples_NAME:"
echo ${Control_Samples_NAME}
echo "Treatment_Samples_NAME:"
echo ${Treatment_Samples_NAME}
echo "COMPAR_NUM:"
echo ${COMPAR_NUM}
echo "WINDOW_SIZE:"
echo ${WINDOW_SIZE}
echo "FRAG_SIZE:"
echo ${FRAG_SIZE}
echo "Input_Sites_RiPPM:"
echo ${Input_Sites_RiPPM}
echo "SCRIPT_DIR:"
echo ${SCRIPT_DIR}
echo "TIME_LIMIT:"
echo ${TIME_LIMIT}
echo "-----------------------"
echo "End of variable list"
echo "-----------------------"
#---------------------------------------------------------------------------------
echo "-----------------------"
echo "Start of qsub commands:"
echo "-----------------------"
##################################################################################
	#Now use arguments in the qsub script call:
	qsub -N ${Job_Name}'_diffReps_'${COMPAR_NUM} -P wax-es -q !linga diffReps.qsub ${Dataset_DIR} ${Dataset_Label} ${Control_Samples_NAME} ${Treatment_Samples_NAME} ${COMPAR_NUM} ${WINDOW_SIZE} ${FRAG_SIZE} ${Input_Sites_RiPPM} ${SCRIPT_DIR}
	#Printing the qsub command which submits the job
	echo "qsub -N ${Job_Name}'_diffReps_'${COMPAR_NUM} -P wax-es -q !linga diffReps.qsub ${Dataset_DIR} ${Dataset_Label} ${Control_Samples_NAME} ${Treatment_Samples_NAME} ${COMPAR_NUM} ${WINDOW_SIZE} ${FRAG_SIZE} ${Input_Sites_RiPPM} ${SCRIPT_DIR}"
##################################################################################
