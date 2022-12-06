set -eu
##################################################################################
#Andy Rampersaud, 12.13.18
#This 02_Review_Pipeline_Parameters.sh script is used to print all parameters used in the pipeline
#The purpose of this script is to organize pipeline parameters:
#	1. Create a single text file containing all pipeline parameters
#	2. Give the user a chance to review these parameters before running the pipeline
#	3. Provide notation about how the pipeline was run
#Assumptions for this script to work: 
#	1. This script is located in the "00_Setup_Pipeline" folder
#	2. The 00_Setup_Pipeline folder is in a "Scripts" folder
#	3. The Scripts folder contains all of the pipeline steps to run
#Output:
#	1. Text file: Pipeline_Parameters.txt
#Way to run the script:
#Usage: ./02_Review_Pipeline_Parameters.sh > <Data_Label>_Pipeline_Parameters.txt
#Example command: 
#./02_Review_Pipeline_Parameters.sh > G113_Pipeline_Parameters.txt
##################################################################################
#---------------------------------------------------------------------------------
#Print header to output file:

echo '###################################################################' 
echo "List of Pipeline Parameters" 
echo "Please check that the following parameters are correct before running the pipeline." 
echo "02_Review_Pipeline_Parameters.sh script run date:"
echo $(date)
Pipeline_Version='v2.5.0'
echo "Pipeline version:"
echo ${Pipeline_Version}
echo 'Refer to Pipeline_Version_History.txt for details.'
echo '###################################################################' 

#---------------------------------------------------------------------------------
#Source the setup file to initialize variables
source ./01_Pipeline_Setup.sh
#---------------------------------------------------------------------------------
#Print variables to output file:

echo '#--------------------------------------------------------------' 
echo '"Global Variables"  = variables used by multiple steps' 
echo '#--------------------------------------------------------------' 
echo 'Dataset_DIR='
echo ${Dataset_DIR}
echo 'BU_User='
echo ${BU_User}
echo 'VM_DIR_FASTQC='
echo ${VM_DIR_FASTQC}
echo '#--------------------------------------------------------------'
echo 'Check if VM_DIR_FASTQC exists:'

set +eu
if ( [ -d ${VM_DIR_FASTQC} ] )
then
    echo 'VM_DIR_FASTQC exists!'
else
    echo 'Need to create VM_DIR_FASTQC!'
    mkdir -p ${VM_DIR_FASTQC}
    echo 'Just created the remote directory.'
fi
set -eu

echo '#--------------------------------------------------------------'
echo 'Dataset_Label='
echo ${Dataset_Label}
echo 'GTF_Files_DIR='
echo ${GTF_Files_DIR}
echo '#--------------------------------------------------------------' 
echo '"Step-specific Variables"  = variables used by a particular step'
echo '#--------------------------------------------------------------'
echo '#---------------------------'
echo 'Job: TopHat_Paired_End'
echo 'Bowtie2Index_DIR='
echo ${Bowtie2Index_DIR}
echo '#---------------------------'
echo 'Job: MACS2'
echo 'DHS_DATA='
echo ${DHS_DATA}
echo '#---------------------------'
echo 'Job: MEME_ChIP'
echo 'Input_FASTA_DIR='
echo ${Input_FASTA_DIR}
echo 'SEQ_LEN='
echo ${SEQ_LEN}
echo 'MEME_MOD='
echo ${MEME_MOD}
echo '#---------------------------'
echo 'Job: bamCorrelate'
echo 'zMin='
echo ${zMin}
echo 'zMax'
echo ${zMax}
echo '#---------------------------'
echo '#--------------------------------------------------------------'
echo '"Additional Variables" = variables not user-specified'
echo '#--------------------------------------------------------------'
echo '#---------------------------'
echo 'Sample_Labels_DIR'
echo ${Sample_Labels_DIR}
echo 'SCRIPT_DIR'
echo ${SCRIPT_DIR}
echo '#---------------------------'

echo 'Job: UCSC_BigWig'
echo 'VM_DIR_UCSC'
echo ${VM_DIR_UCSC}
echo '#--------------------------------------------------------------'
echo 'Check if VM_DIR_UCSC exists:'
set +eu
if ( [ -d ${VM_DIR_UCSC} ] )
then
    echo 'VM_DIR_UCSC_FASTQC exists!'
else
    echo 'Need to create VM_DIR_UCSC!'
    mkdir -p ${VM_DIR_UCSC}
    echo 'Just created the remote directory.'
fi
set -eu
echo '#--------------------------------------------------------------'
echo '#---------------------------'
echo 'TIME_LIMIT'
echo ${TIME_LIMIT}
echo '#--------------------------------------------------------------'
echo 'diffReps-specific variables'
echo '#--------------------------------------------------------------'
#---------------------------------------------------------------------------------
#The Scripts folder is 1 level up from 00_Setup_Pipeline:
#--------------------------------------------------------------------------------
# cd ..
# Scripts_DIR=$(pwd)
# #--------------------------------------------------------------------------------
# #Need a list of steps (same way in 03_Run_Pipeline.sh):
# Pipeline_Steps=$(ls -l | awk 'FNR>1 {print $9}' | sed -n '/[0-9]/p' | grep 'diffReps')
# #For loop over steps in the pipeline
# for step in ${Pipeline_Steps}
# do
#     echo '#---------------------------'
#     echo 'Job: '${step} 
#     echo '#---------------------------'
#     cd ${step}
#     #Source the setup file to initialize variables
#     source ./setup_diffReps.sh
#     echo 'Control_Samples_NAME'
#     echo ${Control_Samples_NAME}
#     echo 'Treatment_Samples_NAME'
#     echo ${Treatment_Samples_NAME}
#     echo 'COMPAR_NUM'
#     echo ${COMPAR_NUM}
#     echo 'WINDOW_SIZE'
#     echo ${WINDOW_SIZE}
#     echo 'FRAG_SIZE'
#     echo ${FRAG_SIZE}
#     cd ${Scripts_DIR}
# done
echo '#--------------------------------------------------------------'
echo 'Generate_Tracks:'
echo 'No parameters to print.'
echo 'Check that the Sample_Labels.txt has a color column.'
echo '#--------------------------------------------------------------'
echo '###################################################################' 
