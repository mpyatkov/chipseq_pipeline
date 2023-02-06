#!/bin/bash -l
set -eu
##################################################################################
# Andy Rampersaud, 07.19.17
#This script would be used to summarize SICER jobs
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
echo "TIME_LIMIT:"
echo ${TIME_LIMIT}
echo "-----------------------"
echo "End of variable list"
echo "-----------------------"
#---------------------------------------------------------------------------------
echo
echo 'Loading required modules...'
echo
#Make sure the shebang line = #!/bin/bash -l
set -eu
#Need the -l option to load modules
#Search for latest program installed:
#module avail -t 2>&1 | grep -i bedtools
module load bedtools/2.27.1
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
OUTPUT_FILE=$OUTPUT_DIR/SICER_Stats.txt
OUTPUT_FILE_2=$OUTPUT_DIR/SICER_Peak_Width.txt
#Need to create output files:
######################
if [ -f $OUTPUT_FILE ]
then 
    rm $OUTPUT_FILE
else
    touch $OUTPUT_FILE
fi
######################
######################
if [ -f $OUTPUT_FILE_2 ]
then 
    rm $OUTPUT_FILE_2
else
    touch $OUTPUT_FILE_2
fi
######################
#---------------------------------------------------------------------------------
Overlap_Summary=$OUTPUT_DIR/ENCODE_Blacklist_Overlap_Summary.txt
######################
if [ -f $Overlap_Summary ]
then 
    rm $Overlap_Summary
else
    touch $Overlap_Summary
fi
######################
#---------------------------------------------------------------------------------
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
echo SAMPLE_ID $'\t'Description $'\t'FRAGMENT_COUNT $'\t'SICER_PEAK_COUNT $'\t'FRAGMENT_IN_PEAK_COUNT $'\t'FRAGMENT_IN_PEAK_RATIO $'\t'Ratio_Genome_Covered >> $OUTPUT_FILE
#Data header:
echo 'Peak Width (bp) Distribution' >> $OUTPUT_FILE_2
#Column headers:
echo Peak Set $'\t'Min. $'\t'1st Qu. $'\t'Median $'\t'Mean $'\t'3rd Qu. $'\t'Max. >> $OUTPUT_FILE_2
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
    #SICER_Stats
    ###############################################
    echo 'Getting read count...'
    #For single-end data search for: "# total tags in treatment:"
    #For paired-end data search for: "total fragments in treatment:"
    #---------------------------------------------------------------------------------
    #I need a way to count the total number of reads
    #Save the line count of the fragments BED file
    #Search the job log file
    #---------------------------------------------------------------------------------
    cd ${SCRIPT_DIR}
    #Search for echo statement and extract the line after the match:
    FRAGMENT_COUNT=$(grep -A 1 'Line count of fragments BED file:' *${Sample_ID}'.o'* | awk 'FNR==2{print $0}')
    #---------------------------------------------------------------------------------
    cd $DATA_DIR
    cd ${Sample_ID}
    cd fastq
    cd bowtie2
    #---------------------------
    cd ${Sample_ID}'_SICER_output'
    #---------------------------
    echo 'Calculating peak count...'
    #General way to refer to the main output file
    Peaks_Called_File=$(ls *.scoreisland)
    Peak_count=$(wc -l < ${Peaks_Called_File})
    echo 'Calculating read in peak count...'
    FRAGMENT_IN_PEAK_COUNT=$(awk '{n+=$5;} ; END {print n;}' ${Sample_ID}'_read_SICER.out1')
    echo 'Calculating read in peak ratio...'
    FRAGMENT_IN_PEAK_RATIO=$(echo "scale=4;$FRAGMENT_IN_PEAK_COUNT/$FRAGMENT_COUNT" | bc)
    echo 'Printing peak width distribution...'
    #Only print the 2nd line of the text file:
    awk 'FNR==2 {print $0}' ${Sample_ID}*'_Width_Stats.txt' >> $OUTPUT_FILE_2
    #---------------------------------------------------------------------------------
    echo 'Calculating the ratio of genome covered...'
    peak_width_sum=$(awk '{peak_width_sum+=$3 - $2;} ; END {print peak_width_sum;}' ${Peaks_Called_File})
    #For macs2 (mm:	1.87e9):
    #https://github.com/taoliu/MACS
    #Definition:
    #It's the mappable genome size or effective genome size which is defined as the genome size which can be sequenced. Because of the repetitive features on the chromsomes, the actual mappable genome size will be smaller than the original size, about 90% or 70% of the genome size. The default hs — 2.7e9 is recommended for UCSC human hg18 assembly. Here are all precompiled parameters for effective genome size
    #https://www.biostars.org/p/19380/
    #macs2 README:
    #http://liulab.dfci.harvard.edu/MACS/00README.html
    Effective_Genome_Size=1870000000
    Ratio_Genome_Covered=$(echo "scale=4;${peak_width_sum} /${Effective_Genome_Size}" | bc)
    #---------------------------------------------------------------------------------
    echo 'Summarizing overlap comparisons:'
    #---------------------------------------------------------------------------------
    folder_list=*_Output
    for folder in $folder_list
    do
	echo $folder
	sed -n '1,5p' $folder/Overlap_Summary.txt >> $Overlap_Summary
	echo >> $Overlap_Summary
    done
    #---------------------------------------------------------------------------------
    ###############################################
    echo 'Printing to output file'
    echo
    echo ${Sample_ID} $'\t'$Description $'\t'$FRAGMENT_COUNT $'\t'$Peak_count $'\t'$FRAGMENT_IN_PEAK_COUNT $'\t'$FRAGMENT_IN_PEAK_RATIO $'\t' ${Ratio_Genome_Covered} >> $OUTPUT_FILE
    cd $DATA_DIR
done < Sample_Labels.temp
#Remove the temp file:
cd ${SCRIPT_DIR}
rm Sample_Labels.temp
#---------------------------------------------------------------------------------
##################################################################################
#---------------------------------------------------------------------------------
#Copy all the sample specific *.scoreisland to a single dir
#Run Chr_Peak_Counts.* scripts to generate the dataset summary (Chr_Peak_Counts.pdf)
#---------------------------------------------------------------------------------
#Generate the *_Chr_Peak_Counts.pdf
OUTPUT_DIR_2=${OUTPUT_DIR}/Chr_Peak_Counts
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
    #Copy each sample's *.scoreisland file to ${OUTPUT_DIR_2}:
    cd ${DATA_DIR}/${Sample_ID}/fastq/bowtie2/${Sample_ID}'_SICER_output'
    #General way to refer to the main output file
    Peaks_Called_File=$(ls *.scoreisland)
    #Note: the Chr_Peak_Counts.sh looks for a .bed extension
    cp ${Peaks_Called_File} ${OUTPUT_DIR_2}/${Peaks_Called_File}'.bed'
    cd ${OUTPUT_DIR_2}
    #---------------------------
done < Sample_Labels.temp
cd ${OUTPUT_DIR_2}
#Remove the temp file:
rm Sample_Labels.temp
rm Sample_Labels.txt
#Copy Chr_Peak_Counts scripts to this dir:
cp -r ${SCRIPT_DIR}/Job_Scripts/Chr_Peak_Counts/* ${OUTPUT_DIR_2}
cd ${OUTPUT_DIR_2}
echo
echo 'Running Chr_Peak_Counts.* scripts'
echo
module load R/3.6.2
./Chr_Peak_Counts.sh
#Remove script files:
rm -r ${OUTPUT_DIR_2}/genomeIndex
rm ${OUTPUT_DIR_2}/Chr_Peak_Counts.R
rm ${OUTPUT_DIR_2}/Chr_Peak_Counts.sh
#---------------------------
##################################################################################
cd ${OUTPUT_DIR}
echo '#---------------------------------------------------------------------------'
echo
echo 'Create and populate Peaks_Called folder'
echo
#---------------------------------------------------------------------------------
#Create Peaks_Called_DIR:
Peaks_Called_DIR=${OUTPUT_DIR}/Peaks_Called
#Using an if statement to make the output folder if it does not exists
if [ ! -d ${Peaks_Called_DIR} ]; then
    mkdir ${Peaks_Called_DIR}; 
else
    rm ${Peaks_Called_DIR}/*
fi
#---------------------------------------------------------------------------------
echo
echo 'Copy sample-specific peak BED files'
echo
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
    #Copy sample-specific sample peak files 
    cd ${DATA_DIR}/${Sample_ID}/fastq/bowtie2/${Sample_ID}'_SICER_output'
    #General way to refer to the main output file
    Peaks_Called_File=$(ls *.scoreisland)
    #Note: the Chr_Peak_Counts.sh looks for a .bed extension
    cp ${Peaks_Called_File}  ${Peaks_Called_DIR}/${Peaks_Called_File}'.bed'
    cd ${OUTPUT_DIR}
    #---------------------------
done < Sample_Labels.temp
#Remove the temp file:
rm Sample_Labels.temp
rm Sample_Labels.txt
###############################################
#---------------------------------------------------------------------------------
#Generate a Peak_Union BED file from sample specific BED files:
#---------------------------------------------
Union_DIR=${Peaks_Called_DIR}/Peak_Union
if [ ! -d ${Union_DIR} ]
then 
    mkdir ${Union_DIR}
    touch ${Union_DIR}/temp.txt
else
    touch ${Union_DIR}/temp.txt
fi
#---------------------------------------------
echo 'Generate Peak_Union sites...'
cd ${Peaks_Called_DIR}
echo '#---------------------------------------------'
echo 'Within the Peaks_Called directory'
echo 'List of BED files and line counts:'
wc -l *.bed
echo '#---------------------------------------------'
file_list=*.bed
echo 'Concatenating files ...'
for file in $file_list
do
    #echo $file
    cat $file ${Union_DIR}/temp.txt >> ${Union_DIR}/temp2.txt
    mv ${Union_DIR}/temp2.txt ${Union_DIR}/temp.txt
done
#Confirmed this worked:
#echo 'Line count for temp.txt:'
#wc -l ${Union_DIR}/temp.txt | awk '{print $1}'
#echo 'This is the sum of all the peak counts'
#cd $Peaks_Called_DIR
#wc -l *.bed | grep total
cd ${Union_DIR}
echo 'Sorting ...'
sort -k1,1 -k2,2n temp.txt > temp2.txt
echo 'Running bedtools merge...'
bedtools merge -i temp2.txt > temp3.txt
mv temp3.txt Peak_Union.bed
rm temp*.txt
#Add a peak ID for easier processing (sorting, table merge, etc...)
#http://askubuntu.com/questions/231995/how-to-separate-fields-with-space-or-tab-in-awk
awk -v OFS='\t' '{print $1, $2, $3, "Union_Peak_"NR}' Peak_Union.bed > temp.bed
mv temp.bed Peak_Union.bed
echo 'Created Peak_Union.bed, the line count:'
wc -l Peak_Union.bed | awk '{print $1}'
echo '#-------------------------------------------------------------------------'
echo 'Copying Peak_Union.bed to the RiPPM_Normalization step Input folder'
echo 'Note: this Peak_Union.bed file will be used for RiPPM_Normalization.'
echo '#-------------------------------------------------------------------------'
cd ${OUTPUT_DIR}
##################################################################################
#This copy command will be controlled by the value of ${Input_Sites_RiPPM}
#Check the value of ${Input_Sites_RiPPM}
#---------------------------------------------------------------------------------
if [ "${Input_Sites_RiPPM}" == "MACS2" ] ;
then
    echo 'WARNING: ${Input_Sites_RiPPM} was set to MACS2.'
    echo 'Therefore, the Peak_Union.bed generated by this job (SICER) will not be used.'
    echo 'Please run the MACS2 job with the correspnding Summarize_Jobs.sh script to generate the appropriate Peak_Union.bed file.'
fi
#---------------------------------------------------------------------------------
if [ "${Input_Sites_RiPPM}" == "SICER" ] ;
then
    RiPPM_Normalization_Job_Name=09_RiPPM_Normalization
    cp ${OUTPUT_DIR}/Peaks_Called/Peak_Union/Peak_Union.bed ../../${RiPPM_Normalization_Job_Name}/Peak_Union
fi
#---------------------------------------------------------------------------------
#Check if the ${Input_Sites_RiPPM} was set:
#Need the "and" conditional (&&):
#I want to check if the variable does not match both values
if [ "${Input_Sites_RiPPM}" != "MACS2" ] && [ "${Input_Sites_RiPPM}" != "SICER" ]
then
    echo 'WARNING: ${Input_Sites_RiPPM} was not set to either MACS2 or SICER.'
    echo 'Zero input sites will be used for RiPPM normalization.'
    echo 'Specify either MACS2 or SICER for the Input_Sites_RiPPM variable.'
fi
#---------------------------------------------------------------------------------
##################################################################################
echo '#-------------------------------------------------------------------------'
echo 'Check out '$OUTPUT_DIR
#echo 'Check out '$OUTPUT_FILE
#echo 'Check out '$OUTPUT_FILE_2
#echo 'Check out '${OUTPUT_DIR_2}
#echo 'Check out '${Overlap_Summary}
echo '#-------------------------------------------------------------------------'
##################################################################################
