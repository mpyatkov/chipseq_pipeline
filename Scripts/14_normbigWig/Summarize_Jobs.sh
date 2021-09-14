#!/bin/bash
set -eu
##################################################################################
# Andy Rampersaud, 07.19.17
#This script would be used to summarize normbigWig
#(need the *.o* files from the normbigWig job)
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
echo "BU_User:"
echo ${BU_User}
echo "VM_DIR_UCSC:"
echo ${VM_DIR_UCSC}
echo "TIME_LIMIT:"
echo ${TIME_LIMIT}
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
cd $INPUT_DIR
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
#---------------------------------------------------------------------------------
echo 'Copying over BAM, bigWig, and bigBed files to waxmanlabvm...'
#BAM file of reads:
cp ${Dataset_DIR}/${Sample_ID}/fastq/bowtie2/$Sample_ID'_sorted_mapped.bam' ${VM_DIR_UCSC}
#BAM file index:
cp ${Dataset_DIR}/${Sample_ID}/fastq/bowtie2/$Sample_ID'_sorted_mapped.bam.bai' ${VM_DIR_UCSC}
#---------------------------------------------------------------------------------
#RiPPM BigWig file of reads:
cp ${Dataset_DIR}/${Sample_ID}/fastq/bowtie2/normbigWig/*.bw ${VM_DIR_UCSC}
#---------------------------------------------------------------------------------
#BigBed file of peaks:
#macs2:
cp ${Dataset_DIR}/${Sample_ID}/fastq/bowtie2/$Sample_ID'_MACS2_output'/*.bb ${VM_DIR_UCSC}
#SICER:
cp ${Dataset_DIR}/${Sample_ID}/fastq/bowtie2/$Sample_ID'_SICER_output'/*.bb ${VM_DIR_UCSC}
#---------------------------------------------------------------------------------
done < Sample_Labels.temp
#Remove the temp file:
rm Sample_Labels.temp
rm Sample_Labels.txt
##################################################################################
echo '#--------------------------------------------------------------------------'
echo 'Check out '${INPUT_DIR}'/Input for the Norm_Factors.txt'
echo 'UCSC files for data visualization should now be copied to the VM.'
echo '#--------------------------------------------------------------------------'
##################################################################################
