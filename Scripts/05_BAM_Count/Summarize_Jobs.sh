#!/bin/bash
set -eu
##################################################################################
# Andy Rampersaud, 03.21.17
#This script would be used to summarize BAM_Count jobs
#Way to run script:
#Usage: 
#./BAM_Count_Summary.sh
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
echo "-----------------------"
echo "End of variable list"
echo "-----------------------"
#---------------------------------------------------------------------------------
Current_DIR=$(pwd)
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
    #Don't remove dir (storing files from jobs)
fi
###############################
##################################################################################
#---------------------------------------------------------------------------------
OUTPUT_FILE=$OUTPUT_DIR/BAM_Count_Stats.txt
######################
if [ -f $OUTPUT_FILE ]
then 
    rm $OUTPUT_FILE
else
    touch $OUTPUT_FILE
fi
######################
cd $DATA_DIR
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
#Print header to output file:
echo 'SAMPLE_ID' $'\t''Description' $'\t''TOTAL_READ_COUNT(R1 + R2)' $'\t''MAPPED_READ_COUNT(properly paired uniquely mapped)(R1 + R2)' $'\t''MAPPED_READ_RATIO'  >> $OUTPUT_FILE
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
    cd ${Sample_ID}/
    cd fastq/
    cd bowtie2/
    #Get total read count:
    Total_Read_Count=$(head -1 flagstat_${Sample_ID}.txt | awk '{print $1}')
    #Need to modify idxstats.out file to exclude chr*_random and chrM (and * (has a zero read count))
    grep -v '_random' ${Sample_ID}'_idxstats.out' | grep -v 'chrM' | grep -v '*' > ${Sample_ID}'_temp_idxstats.out'
    #The output is TAB delimited with each line consisting of reference sequence name, sequence length, # mapped reads and # unmapped reads. 
    #Need to get sum of 3rd column
    echo 'Getting read count...'
    Mapped_Read_Count=$(awk '{n+=$3;} ; END {print n;}' ${Sample_ID}'_temp_idxstats.out')
    Mapped_Read_Ratio=$(echo "scale=4;($Mapped_Read_Count/$Total_Read_Count)" | bc)
    echo 'Printing to output file'
    echo ${Sample_ID} $'\t'$Description $'\t'$Total_Read_Count $'\t'$Mapped_Read_Count $'\t'$Mapped_Read_Ratio >> $OUTPUT_FILE
    #Remove the *_temp_idxstats.out file
    rm ${Sample_ID}'_temp_idxstats.out'
    cd ..
    cd ..
    cd ..
done < Sample_Labels.temp
#Remove the temp file:
rm Sample_Labels.temp
rm Sample_Labels.txt
#---------------------------------------------------------------------------------
##################################################################################
echo 'Summarize the Read_Strand_Count'
OUTPUT_DIR_2=${OUTPUT_DIR}/Read_Strand_Count
#OUTPUT_DIR_2 already contains the summary from each sample's job
OUTPUT_FILE_2=$OUTPUT_DIR_2/Read_strand_count_Stats.txt
######################
if [ -f $OUTPUT_FILE_2 ]
then 
    rm $OUTPUT_FILE_2
else
    touch $OUTPUT_FILE_2
fi
######################
cd $OUTPUT_DIR_2
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
#Print header to output file:
echo 'SAMPLE_ID' $'\t''Positive_Strand_READ_COUNT' $'\t''Negative_Strand_READ_COUNT' $'\t''Total_Read_Count' $'\t''Count_Difference(Postive-Negative)' $'\t''Count_Difference_Ratio' >> $OUTPUT_FILE_2
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
    cat ${Sample_ID}'_read_strand_count.txt' >> $OUTPUT_FILE_2
done < Sample_Labels.temp
#Remove the temp file:
rm Sample_Labels.temp
rm Sample_Labels.txt
#---------------------------------------------------------------------------------
##################################################################################
echo 'Generate the *_Chr_Read_Counts.pdf'
OUTPUT_DIR_3=${OUTPUT_DIR}/Chr_Read_Counts
######################
if [ ! -d ${OUTPUT_DIR_3} ]
then 
    mkdir ${OUTPUT_DIR_3}
else
    #Remove dir:
    rm -r ${OUTPUT_DIR_3}
    #Make new dir:
    mkdir -p ${OUTPUT_DIR_3}
fi
######################
cd ${OUTPUT_DIR_3}
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
    #Copy each sample's *_idxstats.out file to ${OUTPUT_DIR_3}:
    cp ${DATA_DIR}/${Sample_ID}/fastq/bowtie2/*_idxstats.out ${OUTPUT_DIR_3}
done < Sample_Labels.temp
#Remove the temp file:
rm Sample_Labels.temp
rm Sample_Labels.txt
#Copy Chr_Read_Counts scripts to this dir:
cp -r ${SCRIPT_DIR}/Job_Scripts/Chr_Read_Counts/* ${OUTPUT_DIR_3}
cd ${OUTPUT_DIR_3}
echo
echo 'Running Chr_Read_Counts.* scripts'
echo
module load R/3.6.2
./Chr_Read_Counts.sh
#Remove script files:
rm -r ${OUTPUT_DIR_3}/genomeIndex
rm ${OUTPUT_DIR_3}/Chr_Read_Counts.R
rm ${OUTPUT_DIR_3}/Chr_Read_Counts.sh
#---------------------------
##################################################################################
#---------------------------------------------------------------------------------
echo 'Summarize the Fragment_Count'
OUTPUT_DIR_4=${OUTPUT_DIR}/Fragment_Count
#OUTPUT_DIR_4 already contains the summary from each sample's job
OUTPUT_FILE_4=$OUTPUT_DIR_4/Fragment_count_Stats.txt
######################
if [ -f $OUTPUT_FILE_4 ]
then 
    rm $OUTPUT_FILE_4
else
    touch $OUTPUT_FILE_4
fi
######################
cd $OUTPUT_DIR_4
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
#Print header to output file:
echo 'SAMPLE_ID' $'\t''Fragment Count (properly paired - BEDPE file line count)' >> $OUTPUT_FILE_4
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
    cat ${Sample_ID}'_fragment_count.txt' >> $OUTPUT_FILE_4
done < Sample_Labels.temp
#Remove the temp file:
rm Sample_Labels.temp
rm Sample_Labels.txt
#---------------------------------------------------------------------------------
##################################################################################
echo '#--------------------------------------------------------------------------'
echo 'Check out '${OUTPUT_DIR}
echo '#--------------------------------------------------------------------------'
##################################################################################
