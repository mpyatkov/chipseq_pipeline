##################################################################################
#Andy Rampersaud, 08.04.17
#This 03_Run_Pipeline.sh script is used to run all steps in the pipeline
#The purpose of this script is to stream-line the pipeline so that 
#	1. A single command executes all steps sequentially
#	2. Minimize user intervention by automatically moving to the next step
#Assumptions for this script to work: 
#	1. This script is located in the "00_Setup_Pipeline" folder
#	2. The 00_Setup_Pipeline folder is in a "Scripts"
#	3. The Scripts folder contains all of the pipeline steps to run
#Way to run script:
#Usage: ./03_Run_Pipeline.sh <Start_Step>
#	<Start_Step> = Needs to be either "Full_Pipeline" or a specific pipeline step 
#Example: 
#User wants to run the full pipeline:
#./03_Run_Pipeline.sh Full_Pipeline
#User wants to run a subset of the pipeline starting at a specific step (04_TopHat_Paired_End):
#./03_Run_Pipeline.sh 04_TopHat_Paired_End
##################################################################################
#---------------------------------------------------------------------------------
#Process system argument:
set -e

if [ $# -ne 1 ]
then
    echo "Usage: ./03_Run_Pipeline.sh <Start_Step>"
    echo "<Start_Step> = Needs to be either \"Full_Pipeline\" or a specific pipeline step"
    echo "See 03_Run_Pipeline.sh for details."	
    exit 
fi
Start_Step=$1
#---------------------------------------------------------------------------------
#Use to calculate job time:
#Start_Time in seconds
Start_Time=$(date +"%s")
#---------------------------------------------------------------------------------
#Source the setup file to initialize variables
source ./01_Pipeline_Setup.sh

# if we have only comments lines in diffreps config
# then just skip the creation of steps 10
numlines_diffreps_config=`cat diffreps_config.txt | sed -e '/^#/d' | wc -l`
if [ "${numlines_diffreps_config}" -ne "0" ]; then
  ./diffreps_parser.sh
fi

if ( [ ! -d ${VM_DIR_FASTQC} ] ); then
    echo "Creating: ${VM_DIR_FASTQC}"
    mkdir -p ${VM_DIR_FASTQC}
fi

if ( [ ! -d ${VM_DIR_UCSC} ] ); then
    echo "Creating: ${VM_DIR_UCSC}"
    mkdir -p ${VM_DIR_UCSC}
fi

#Save 00_Setup_Pipeline location:
Setup_Pipeline_DIR=$(pwd)
#---------------------------------------------------------------------------------
#The Scripts folder is 1 level up from 00_Setup_Pipeline:
#--------------------------------------------------------------------------------
cd ..
Scripts_DIR=$(pwd)
#--------------------------------------------------------------------------------
##################################################################################
#Start If statement for ${Start_Step}
if [ ${Start_Step} == "Full_Pipeline" ]
then
    #For the 01_Rename_Folders - I just want to run this step and move to the rest of the pipeline
    Step_DIR=$(ls -l | awk 'FNR>1 {print $9}' | sed -n '/Rename_Folders/p')
    echo '#--------------------------------------------------------------------------'
    echo 'Running: '${Step_DIR}'...' 
    echo '#--------------------------------------------------------------------------'
    #Run the step:
    # cd ${Step_DIR}
    # ./Rename_Folders.sh
    # cd ..
    # echo '#--------------------------------------------------------------------------'
    # echo 'Done: '${Step_DIR}
    # echo '#--------------------------------------------------------------------------'
    # #--------------------------------------------------------------------------------
    #Full list (omit the 1st blank line):
    #ls -l | awk 'FNR>1 {print $9}'
    #http://www.tutorialspoint.com/unix/unix-regular-expressions.htm
    #http://www.grymoire.com/Unix/Sed.html
    #Need the -n option:
    #The "-n" option will not print anything unless an explicit request to print is found. 
    #The pattern is a 2-digit number and then a bunch of characters
    #Add sed command for regular expression (2-digit number at the beginning):
    Pipeline_Steps=$(ls -l | awk 'FNR>1 {print $9}' | sed -n '/[0-9]/p' | sed '/00_Setup_Pipeline/d' | sed '/TEMPLATE_10_diffreps/d' | sed '/Rename_Folders/d' | sed '/Generate_Tracks/d')
    #For loop over steps in the pipeline
    for step in ${Pipeline_Steps}
    do
	echo '#--------------------------------------------------------------------------'
	echo 'Running: '${step}'...' 
	echo '#--------------------------------------------------------------------------'
	cd ${step}
	#Run_Jobs:
	./Run_Jobs.sh
	#--------------------------------------------------------------------------------
	#Need a way to periodically check that all jobs have completed
	#Count the number of jobs submitted by the user
	#Once the count goes to zero then summarize this job and move to the next step
	#Use the "${BU_User}" variable from the 01_Pipeline_Setup.sh
	#Need to omit the 1st 2 lines of qstat command:
	#Job_Count=$(qstat -u ${BU_User} | awk 'FNR>2 {print $0}' | wc -l)
	#echo ${Job_Count} 
	#--------------------------------------------------------------------------------
	#Retrieve the job name for this step:
	Current_DIR=$(pwd)
	#Extract the folder name:
	DIR_name=`basename ${Current_DIR}`
	#Split by "_" grab the 1st part (use cut command)
	Step_Num=$(echo ${DIR_name} | cut -d'_' -f 1)
	#Create the Job_Name:
	Job_Name=$(echo 'Step_'${Step_Num}'_'${Dataset_Label})
	#Use ${Job_Name} for qsub -N parameter
	#--------------------------------------------------------------------------------
	#Better to check: qstat -r -u ${BU_User}
	#The (-r) option gives the Full jobname
	#Extract lines with "Full jobname:"
	#Extract lines with step-specific name
	#Count the number of lines (number of running jobs)
	Job_Count=$(qstat -r -u ${BU_User} | grep 'Full jobname:' | grep ${Job_Name} | wc -l)
	#echo ${Job_Count}
	#--------------------------------------------------------------------------------
	echo ${BU_User}" running "${Job_Count}" jobs"
	echo "Note: the user may have other jobs or pipeline instances running."
	echo "Only jobs in this pipeline instance will be monitored or counted."
	echo "Periodically checking until jobs complete (please wait)..."
	#Use a while loop to check ${Job_Count}
	while [ ${Job_Count} -ne 0 ]
	do
	    #Wait 01 minute then check ${Job_Count} again
	    sleep 1m
	    Job_Count=$(qstat -r -u ${BU_User} | grep 'Full jobname:' | grep ${Job_Name} | wc -l)
	done
	#--------------------------------------------------------------------------------
	#Summarize_Jobs:
	./Summarize_Jobs.sh
	echo '#--------------------------------------------------------------------------'
	echo 'Done: '${step}
	echo '#--------------------------------------------------------------------------'
	cd ..
    done
    #--------------------------------------------------------------------------------
    #For the 12_Generate_Tracks - I just want to run this step 
    Step_DIR=$(ls -l | awk 'FNR>1 {print $9}' | sed -n '/Generate_Tracks/p')
    echo '#--------------------------------------------------------------------------'
    echo 'Running: '${Step_DIR}'...' 
    echo '#--------------------------------------------------------------------------'
    #Run the step:
    cd ${Step_DIR}
    ./Generate_Tracks.sh
    cd ..
    echo '#--------------------------------------------------------------------------'
    echo 'Done: '${Step_DIR}
    echo '#--------------------------------------------------------------------------'
    #End of If statement for ${Start_Step}
fi
##################################################################################
#Start If statement for ${Start_Step}
if [ ${Start_Step} != "Full_Pipeline" ]
then
    #At the very least the user will be skipping 01_Rename_Folders (already omitted)
    #Only print steps after (but including) ${Start_Step}
    #http://stackoverflow.com/questions/7103531/how-to-get-the-part-of-file-after-the-line-that-matches-grep-expression-first
    #Need: sed -n -e '/TERMINATE/,$p'
    #If using variable within sed command need to enclose the variable in single quotes
    Pipeline_Steps=$(ls -l | awk 'FNR>1 {print $9}' | sed -n '/[0-9]/p' | sed '/00_Setup_Pipeline/d' | sed '/TEMPLATE_10_diffreps/d' | sed '/Rename_Folders/d' | sed '/Generate_Tracks/d' | sed -n -e '/'${Start_Step}'/,$p')
    #Now the pipeline will start at ${Start_Step} and continue to the end
    #For loop over steps in the pipeline
    for step in ${Pipeline_Steps}
    do
	echo '#--------------------------------------------------------------------------'
	echo 'Running: '${step}'...' 
	echo '#--------------------------------------------------------------------------'
	cd ${step}
	#Run_Jobs:
	./Run_Jobs.sh
	#--------------------------------------------------------------------------------
	#Need a way to periodically check that all jobs have completed
	#Count the number of jobs submitted by the user
	#Once the count goes to zero then summarize this job and move to the next step
	#Use the "${BU_User}" variable from the 01_Pipeline_Setup.sh
	#Need to omit the 1st 2 lines of qstat command:
	#Job_Count=$(qstat -u ${BU_User} | awk 'FNR>2 {print $0}' | wc -l)
	#echo ${Job_Count} 
	#--------------------------------------------------------------------------------
	#Retrieve the job name for this step:
	Current_DIR=$(pwd)
	#Extract the folder name:
	DIR_name=`basename ${Current_DIR}`
	#Split by "_" grab the 1st part (use cut command)
	Step_Num=$(echo ${DIR_name} | cut -d'_' -f 1)
	#Create the Job_Name:
	Job_Name=$(echo 'Step_'${Step_Num}'_'${Dataset_Label})
	#Use ${Job_Name} for qsub -N parameter
	#--------------------------------------------------------------------------------
	#Better to check: qstat -r -u ${BU_User}
	#The (-r) option gives the Full jobname
	#Extract lines with "Full jobname:"
	#Extract lines with step-specific name
	#Count the number of lines (number of running jobs)
	Job_Count=$(qstat -r -u ${BU_User} | grep 'Full jobname:' | grep ${Job_Name} | wc -l)
	#echo ${Job_Count}
	#--------------------------------------------------------------------------------
	echo ${BU_User}" running "${Job_Count}" jobs"
	echo "Note: the user may have other jobs or pipeline instances running."
	echo "Only jobs in this pipeline instance will be monitored or counted."
	echo "Periodically checking until jobs complete (please wait)..."
	#Use a while loop to check ${Job_Count}
	while [ ${Job_Count} -ne 0 ]
	do
	    #Wait 01 minute then check ${Job_Count} again
	    sleep 1m
	    Job_Count=$(qstat -r -u ${BU_User} | grep 'Full jobname:' | grep ${Job_Name} | wc -l)
	done
	#--------------------------------------------------------------------------------
	#Summarize_Jobs:
	./Summarize_Jobs.sh
	echo '#--------------------------------------------------------------------------'
	echo 'Done: '${step}
	echo '#--------------------------------------------------------------------------'
	cd ..
    done
    #--------------------------------------------------------------------------------
    #For the 12_Generate_Tracks - I just want to run this step 
    Step_DIR=$(ls -l | awk 'FNR>1 {print $9}' | sed -n '/Generate_Tracks/p')
    echo '#--------------------------------------------------------------------------'
    echo 'Running: '${Step_DIR}'...' 
    echo '#--------------------------------------------------------------------------'
    #Run the step:
    cd ${Step_DIR}
    ./Generate_Tracks.sh
    cd ..
    echo '#--------------------------------------------------------------------------'
    echo 'Done: '${Step_DIR}
    echo '#--------------------------------------------------------------------------'
    #End of If statement for ${Start_Step}
fi
##################################################################################
echo '#--------------------------------------------------------------------------'
echo 'Printing job duration for all steps...'
OUTPUT_FILE=${Setup_Pipeline_DIR}/Pipeline_Runtime.txt
######################
if [ -f ${OUTPUT_FILE} ]
then 
    rm ${OUTPUT_FILE}
else
    touch ${OUTPUT_FILE}
fi
######################
#Print header to output file:
echo '#--------------------------------------------------------------' >> ${OUTPUT_FILE}
echo "Run time for each submitted job:" >> ${OUTPUT_FILE}
echo "Job run times that deviate from the average should be inspected for possible errors (check the job log files)" >> ${OUTPUT_FILE}
echo '#--------------------------------------------------------------' >> ${OUTPUT_FILE}
cd ${Scripts_DIR}
Pipeline_Steps=$(ls -l | awk 'FNR>1 {print $9}' | sed -n '/[0-9]/p' | sed '/00_Setup_Pipeline/d' | sed '/TEMPLATE_10_diffreps/d' | sed '/Rename_Folders/d' | sed '/Generate_Tracks/d')
#For loop over steps in the pipeline
for step in ${Pipeline_Steps}
do
    #Print job runtimes to file
    cd ${step}
    echo ${step} >> ${OUTPUT_FILE}
    set +eu
    grep 'elapsed' *.o* >> ${OUTPUT_FILE}
    set -eu
    cd ..
done
#Also want to print the time to run this script:
echo '#--------------------------------------------------------------' >> ${OUTPUT_FILE}
echo '03_Run_Pipeline.sh run time:' >> ${OUTPUT_FILE}
echo "Check out: "${OUTPUT_FILE}
echo '#--------------------------------------------------------------------------'
##################################################################################
echo 'Pipeline is done, check out your results!'
echo "=========================================================="
echo "Finished on : $(date)"
#Use to calculate job time:
#End_Time in seconds
End_Time=$(date +"%s")
diff=$(($End_Time-$Start_Time))
echo "$(($diff / 3600)) hours, $((($diff / 60) % 60)) minutes and $(($diff % 60)) seconds elapsed."
echo "$(($diff / 3600)) hours, $((($diff / 60) % 60)) minutes and $(($diff % 60)) seconds elapsed." >> ${OUTPUT_FILE}
echo '#--------------------------------------------------------------' >> ${OUTPUT_FILE}
echo "=========================================================="
##################################################################################
