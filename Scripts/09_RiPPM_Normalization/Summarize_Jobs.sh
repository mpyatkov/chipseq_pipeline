#!/bin/bash
##################################################################################
# Andy Rampersaud, 12.13.18
#This script would be used to summarize Peak_Union_Count jobs
#Way to run script:
#Usage: 
#./Summarize_Jobs.sh
##################################################################################
#---------------------------------------------------------------------------------
set -eu
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
DATA_DIR=${Dataset_DIR}
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
OUTPUT_FILE=${OUTPUT_DIR}/Peak_Union_Count_Stats.txt
#Need to create output files:
######################
if [ -f $OUTPUT_FILE ]
then 
rm $OUTPUT_FILE
else
touch $OUTPUT_FILE
fi
######################
#---------------------------------------------------------------------------------
cd $DATA_DIR
################################################
#-----------------------------------------------
#A text file (Sample_Labels.txt) is needed to run this script
Current_DIR=$(pwd)
cp $Sample_Labels_DIR/Sample_Labels.txt ${Current_DIR}
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
echo SAMPLE_ID $'\t'Description $'\t'FRAGMENT_COUNT $'\t'Peak_Union_Count $'\t'FRAGMENT_IN_PEAK_COUNT $'\t'FRAGMENT_IN_PEAK_RATIO >> $OUTPUT_FILE
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
###############################################
#Peak_Union_Count Summary
###############################################
#FRAGMENT_COUNT will be controlled by the value of ${Input_Sites_RiPPM}
#Check the value of ${Input_Sites_RiPPM}
#---------------------------------------------------------------------------------
if [ "${Input_Sites_RiPPM}" == "MACS2" ] ;
then
cd ${Sample_ID}
cd fastq
cd bowtie2
cd ${Sample_ID}'_MACS2_output'
echo 'Getting read count...'
#For single-end data search for: "# total tags in treatment:"
#For paired-end data search for: "total fragments in treatment:"
FRAGMENT_COUNT=$(grep 'total fragments in treatment:' ${Sample_ID}'_MACS2_peaks.xls' | awk '{print $NF}')
fi
#---------------------------------------------------------------------------------
if [ "${Input_Sites_RiPPM}" == "SICER" ] ;
then
cd ${SCRIPT_DIR}
#Need a general way to obtain the job name for the "SICER" step
#List all step folders and find the "SICER" step
SICER_Job_Name=$(ls -al ../ | awk '{print $9}' | grep 'SICER')
cd ../${SICER_Job_Name}
#Search for echo statement and extract the line after the match:
FRAGMENT_COUNT=$(grep -A 1 'Line count of fragments BED file:' *${Sample_ID}'.o'* | awk 'FNR==2{print $0}')
fi
#---------------------------------------------------------------------------------
#Check if the ${Input_Sites_RiPPM} was set:
#Need the "and" conditional (&&):
#I want to check if the variable does not match both values
if [ "${Input_Sites_RiPPM}" != "MACS2" ] && [ "${Input_Sites_RiPPM}" != "SICER" ]
then
echo 'WARNING: ${Input_Sites_RiPPM} was not set to either MACS2 or SICER.'
echo 'The FRAGMENT_COUNT variable will be blank.'
echo 'Specify either MACS2 or SICER for the FRAGMENT_COUNT.'
fi
#---------------------------------------------------------------------------------
echo 'Getting peak count...'
#cd ${Dataset_DIR}/${Sample_ID}/fastq/bowtie2/Peak_Union_Count/Peaks_Called/Peak_Union
#Now saving the Peak_Union.bed in the ${SCRIPT_DIR} (rather than sample-specific folders)
cd ${SCRIPT_DIR}/Peak_Union
Peak_count=$(wc -l < Peak_Union.bed)
cd ${Dataset_DIR}/${Sample_ID}/fastq/bowtie2/Peak_Union_Count
echo 'Getting read in peak count...'
FRAGMENT_IN_PEAK_COUNT=$(awk '{n+=$5;} ; END {print n;}' ${Sample_ID}'_read_Peak_Union.out1')
echo 'Calculating read in peak ratio...'
FRAGMENT_IN_PEAK_RATIO=$(echo "scale=4;$FRAGMENT_IN_PEAK_COUNT/$FRAGMENT_COUNT" | bc)
#---------------------------------------------------------------------------------
###############################################
echo 'Printing to output file'
echo
echo ${Sample_ID} $'\t'$Description $'\t'$FRAGMENT_COUNT $'\t'$Peak_count $'\t'$FRAGMENT_IN_PEAK_COUNT $'\t'$FRAGMENT_IN_PEAK_RATIO >> $OUTPUT_FILE
cd $DATA_DIR
done < Sample_Labels.temp
#Remove the temp file:
rm Sample_Labels.temp
rm Sample_Labels.txt
#---------------------------------------------------------------------------------
##################################################################################
#---------------------------------------------------------------------------------
#Copy all the sample-specific *.out1 to a single dir
#Run scripts to generate the dataset summary 
#---------------------------------------------------------------------------------
OUTPUT_DIR_2=${OUTPUT_DIR}/Peak_Union_Count
######################
if [ ! -d ${OUTPUT_DIR_2} ]
then 
mkdir ${OUTPUT_DIR_2}
else
#Remove dir:
rm -r ${OUTPUT_DIR_2}
#Make new dir:
mkdir -p ${OUTPUT_DIR_2}
fi
######################
cd ${OUTPUT_DIR_2}
################################################
#-----------------------------------------------
#A text file (Sample_Labels.txt) is needed to run this script
Current_DIR=$(pwd)
cp $Sample_Labels_DIR/Sample_Labels.txt ${Current_DIR}
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
#Copy sample-specific *.out1 file to ${OUTPUT_DIR_2}:
cp ${DATA_DIR}/${Sample_ID}/fastq/bowtie2/Peak_Union_Count/${Sample_ID}'_read_Peak_Union.out1' ${OUTPUT_DIR_2}
done < Sample_Labels.temp
#Remove the temp file:
rm Sample_Labels.temp
rm Sample_Labels.txt
#Copy Chr_Peak_Counts scripts to this dir:
#cp -r ${SCRIPT_DIR}/Scripts/Chr_Peak_Counts/* ${OUTPUT_DIR_2}
#cd ${OUTPUT_DIR_2}
#echo
#echo 'Running Chr_Peak_Counts.* scripts'
#echo
#./Chr_Peak_Counts.sh
#Remove script files:
#rm -r ${OUTPUT_DIR_2}/genomeIndex
#rm ${OUTPUT_DIR_2}/Chr_Peak_Counts.R
#rm ${OUTPUT_DIR_2}/Chr_Peak_Counts.sh
#---------------------------
##################################################################################
echo '#-------------------------------------------------------------------------'
echo 'Compiling R Markdown code (loading modules and rendering Rmd script)'
echo 'Note: this may take some time (please allow command to complete)...'
echo '#-------------------------------------------------------------------------'
#Run R_Markdown code to generate report
#Copy R_Markdown files to run in Summary folder:
cp -r ${SCRIPT_DIR}/Job_Scripts/R_Markdown/* ${OUTPUT_DIR}
#Compile the R_Mardown code:
#---------------------------------------------------------------------------------
#Load modules:
#Katia (Wed, May 30, 2018):
#load gcc module first:
#Some libraries were recompiled and some packages in R now require gcc to be loaded.
module load gcc/7.4.0
module load R
module load texlive
module load pandoc/2.2.1
cd ${OUTPUT_DIR}
Rscript -e 'rmarkdown::render("Peak_Union_Count_Stats.Rmd")'
#Clean up the Summary folder:
rm *.tex
rm *.Rmd
echo '#-------------------------------------------------------------------------'
echo 'Copying Norm_Factors.txt to the normbigWig step Input folder'
echo 'Note: this Norm_Factors.txt file will be used for bigWig file creation.'
echo '#-------------------------------------------------------------------------'
#Need a general way to obtain the job name for the "normbigWig" step
#List all step folders and find the "normbigWig" step
normbigWig_Job_Name=$(ls -al ../../ | awk '{print $9}' | grep 'normbigWig')
cp Norm_Factors.txt ../../${normbigWig_Job_Name}/Input
echo '#-------------------------------------------------------------------------'
echo 'Copying Norm_Factors.txt to the RiPPM_Boxplots step Input folder'
echo 'Note: this Norm_Factors.txt file will be used for plotting.'
echo '#-------------------------------------------------------------------------'
RiPPM_Boxplots_Job_Name=$(ls -al ../../ | awk '{print $9}' | grep 'RiPPM_Boxplots')
cp Norm_Factors.txt ../../${RiPPM_Boxplots_Job_Name}/Input/RiPPM
#Go back to job dir:
cd ${SCRIPT_DIR}
##################################################################################
echo '#-------------------------------------------------------------------------'
echo 'Check out '${OUTPUT_DIR}
echo '#-------------------------------------------------------------------------'
##################################################################################
