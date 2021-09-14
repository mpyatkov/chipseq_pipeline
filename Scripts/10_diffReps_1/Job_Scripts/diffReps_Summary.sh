#! /bin/bash

###################################################################################
# Andy Rampersaud, 10.26.2017
#This script would be used to summarize diffReps output
#In this context I use differential to mean "condition-specific"
#In other words, differential sites = condition-specific sites
#Way to run script:
#----------------------------------------------------------------------------------
#Usage: 
#Rscript diffReps_Summary.R <Output File Name> <Sample1 Name> <Sample2 Name> <Peak Caller/Differential program used> <output folder name> <log2FC_cutoff>
#Note:
#If you want 2-fold delta-sites, use <log2FC_cutoff> = 1
#If you want 1.5-fold delta-sites, use <log2FC_cutoff> = 0.5849
#----------------------------------------------------------------------------------
#Options explained in the diffReps_Summary.R
#----------------------------------------------------------------------------------
#Example command for this script: 
#./diffReps_Summary.sh diffReps_1.annotated STAT5_Low STAT5_High diffReps_1 1
#----------------------------------------------------------------------------------
#Notes:
#The Peaks_Called folder must be already created and populated with all of the MACS2 BED files of peaks (can be separate BED file for each sample)
#----------------------------------------------------------------------------------
#Output folders:
#Unfiltered_Summary: 	parsed diffReps output without peak filtering
#Peaks_Filtered_Summary:parsed diffReps output with peak filtering
###################################################################################
#----------------------------------------------------------------------------------
#Check to make sure that all arguments are present:
if [ $# -ne 5 ]
then
	echo "Need 5 arguments:"
	echo "Usage: `basename $0` <Output File Name> <Sample1 Name> <Sample2 Name> <Peak Caller/Differential program used> <output folder name> <log2FC_cutoff>"
	exit 
fi
#Initialize Current_DIR
Current_DIR=$(pwd)
#----------------------------------------------------------------------------------
#Need a diffReps_overlap_peaks_summary.txt
Overlap_Summary=${Current_DIR}/diffReps_overlap_peaks_summary.txt
######################
if [ -f ${Overlap_Summary} ]
then 
rm ${Overlap_Summary}
else
touch ${Overlap_Summary}
fi
######################
#----------------------------------------------------------------------------------
#Need an Overlap_Output_DIR folder:
Overlap_Output_DIR=$Current_DIR/Overlap_Output
if [ ! -d $Overlap_Output_DIR ]
then 
mkdir $Overlap_Output_DIR
fi
#---------------------------------------------------------------------------
#Check that the Peaks_Called folder is present:
Peaks_Called_DIR=$Current_DIR/Peaks_Called
if [ ! -d $Peaks_Called_DIR ]
then
echo '#############################################'
echo 'Warning: Peaks_Called folder must be created!'
echo 'This folder is required to generate the Peaks_Filtered_Summary.'
echo '#############################################'
#Exit script if folder not present
exit
fi
#----------------------------------------------------------------------------------
##Need to get the name of the diffReps output file:
#diffReps_name=${1%\.annotated};
#echo 'This is the file prefix:'
#echo $diffReps_name
echo '#---------------------------------------------------------------------------'
#Generate a Peak_Union BED file from sample specific BED files:
#---------------------------------------------
Union_DIR=$Peaks_Called_DIR/Peak_Union
if [ ! -d $Union_DIR ]
then 
mkdir $Union_DIR
touch $Union_DIR/temp.txt
else
touch $Union_DIR/temp.txt
fi
#---------------------------------------------
#---------------------------------------------------------------------------
#Check that the Input/BED_Regions folder is present:
BED_Regions_DIR=$Current_DIR/Input/BED_Regions
if [ ! -d ${BED_Regions_DIR} ]
then
echo '#############################################'
echo 'Warning: Input/BED_Regions folder is missing!'
echo 'This folder is useful for including a more expansive called-peak list.'
echo 'It is not required for running diffReps but maybe useful for maximizing delta-site overlap with called-peaks.'
echo '#############################################'
fi
#If the folder is present, indicate whether BED files exist
if [ -d ${BED_Regions_DIR} ]
then
cd ${BED_Regions_DIR}
#Want to suppress the error message (when there are zero BED files):
#http://www.unix.com/unix-for-dummies-questions-and-answers/17868-suppressing-error-message-using-ls-command.html

BED_count=$(ls 2>/dev/null -al *.bed  | wc -l)

echo 'Regarding creation of the Peak_Union BED file:'
echo 'Number of additional BED files provided: '${BED_count}
cd ${Current_DIR}
fi
#If there are zero BED files:
#https://ryanstutorials.net/bash-scripting-tutorial/bash-if-statements.php
#Need an "OR" conditional for the if statement
#-z will check if the BED_count variable is blank due to the missing Input/BED_Regions
if [[ ${BED_count} == 0 ]] || [[ -z ${BED_count} ]]
then
echo "No additional BED files will be used to create the Peak_Union BED file."
fi
#If there are BED files present:
#Need an "OR" conditional for the if statement
#-n will check that the length of the variable is greater than zero
if [[ ${BED_count} != 0 ]] && [[ -n ${BED_count} ]]
then
echo "Additional BED files will be used to create the Peak_Union BED file."
#Copy extra user provided BED files from Input/BED_Regions into Peaks_Called_DIR
cp ${BED_Regions_DIR}/*.bed ${Peaks_Called_DIR}
fi
#---------------------------------------------------------------------------
echo 'Generate Peak_Union sites...'
cd $Peaks_Called_DIR
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
cat $file $Union_DIR/temp.txt >> $Union_DIR/temp2.txt
mv $Union_DIR/temp2.txt $Union_DIR/temp.txt
done
#Confirmed this worked:
#echo 'Line count for temp.txt:'
#wc -l $Union_DIR/temp.txt | awk '{print $1}'
#echo 'This is the sum of all the peak counts'
#cd $Peaks_Called_DIR
#wc -l *.bed | grep total
cd $Union_DIR
echo 'Sorting ...'
sort -k1,1 -k2,2n temp.txt > temp2.txt
echo 'Running bedtools merge...'
bedtools merge -i temp2.txt > temp3.txt
mv temp3.txt Peak_Union.bed
rm temp*.txt
echo 'Created Peak_Union.bed, the line count:'
wc -l Peak_Union.bed | awk '{print $1}'
echo '#---------------------------------------------------------------------------'
echo '############################################################################'
echo 'Running Rscript to get Unfiltered_Summary...'
echo '############################################################################'
echo 'Running Rscript ...'
#Use the system arguments
cd $Current_DIR
#----------------------------------------------------------------------------------
Rscript diffReps_Summary.R $1 $2 $3 $4 Unfiltered_Summary $5
#----------------------------------------------------------------------------------
echo '#---------------------------------------------------------------------------'
echo 'Generate BG sites (Peak_Union sites that do not overlap condition-specific sites) ...'
#---------------------------------------------
BG_DIR=$Current_DIR/Unfiltered_Summary/BG_Sites
if [ ! -d $BG_DIR ]
then 
mkdir $BG_DIR
touch $BG_DIR/temp.txt
else
touch $BG_DIR/temp.txt
fi
#---------------------------------------------
diffReps_BED_Files_DIR=$Current_DIR/Unfiltered_Summary/BED_Files
#---------------------------------------------
#The 1st temp file should be the Union peak file (just want to start with a full peak list)
cp $Union_DIR/Peak_Union.bed $BG_DIR
mv $BG_DIR/Peak_Union.bed $BG_DIR/Peak_Union_temp.bed
cd $diffReps_BED_Files_DIR
file_list=*.bed
for file in $file_list
do
echo $file
#---------------------------------------------
#---------------------------------------------
#intersectBed
#-v	Only report those entries in A that have _no overlaps_ with B.
#		- Similar to "grep -v" (an homage).
#---------------------------------------------
intersectBed -v -a $BG_DIR/Peak_Union_temp.bed -b $file >> $BG_DIR/temp.txt
#Take the result of the intersectBed command and use for the next command:
mv $BG_DIR/temp.txt $BG_DIR/Peak_Union_temp.bed
#Sort the Peak_Union_temp.bed before next intersectBed command
sort -k1,1 -k2,2n $BG_DIR/Peak_Union_temp.bed > $BG_DIR/temp2.bed
mv $BG_DIR/temp2.bed $BG_DIR/Peak_Union_temp.bed
done
#---------------------------------------------
#Rename file:
mv $BG_DIR/Peak_Union_temp.bed $BG_DIR/'BG_Sites_'$4'.bed'
#Check that the intersection worked:
echo 'These are the BG_Sites:'
BG_Sites_Count=$(wc -l < $BG_DIR/'BG_Sites_'$4'.bed')
echo $BG_Sites_Count
#The 'BG_Sites_'$4'.bed' should then be used for subsequent enrichment calculations
###################################################################################
#Print statements to both the terminal and Overlap_Summary file (${Overlap_Summary})
echo '#---------------------------------------------------------------------------'
echo '#---------------------------------------------------------------------------' >> ${Overlap_Summary}
#Now I need to overlap the diffReps output with the Peak_Union.bed
cd $Current_DIR
echo 'Pre-process the diffReps output file ...'
echo 'Pre-process the diffReps output file ...' >> ${Overlap_Summary}
#Before running intersectBed I need to sort input files
#Get rid of the header for the *.annotated file
awk 'NR > 1' $Current_DIR/$1 > $Current_DIR/$1'.temp1'
#Print the header for later:
awk 'NR == 1' $Current_DIR/$1 > $Current_DIR/$1'.header'
#Sort file:
sort -k1,1 -k2,2n $Current_DIR/$1'.temp1' > $Current_DIR/$1'.temp2'
echo 'These are the total number of diffReps sites:'
echo 'These are the total number of diffReps sites:' >> ${Overlap_Summary}
diff_total=$(wc -l < $Current_DIR/$1'.temp2')
#http://unix.stackexchange.com/questions/113795/add-thousands-separator-in-a-number
diff_total_Formatted=$(echo ${diff_total} | sed ':a;s/\B[0-9]\{3\}\>/,&/;ta')
echo $diff_total_Formatted
echo $diff_total_Formatted >> ${Overlap_Summary}
echo 'Running bedtools intersect ...'
echo 'Running bedtools intersect ...' >> ${Overlap_Summary}
#Run intersectBed (want the full row from the diffReps output file)
#-wa	Write the original entry in A for each overlap.
#-u	Write the original A entry _once_ if _any_ overlaps found in B.
#		- In other words, just report the fact >=1 hit was found.
#		- Overlaps restricted by -f and -r.
bedtools intersect -wa -u -a  $Current_DIR/$1'.temp2' -b $Union_DIR/Peak_Union.bed > ${Overlap_Output_DIR}/$1'.peaks'
echo 'These are the diffReps sites that overlap peaks:'
echo 'These are the diffReps sites that overlap peaks:' >> ${Overlap_Summary}
diff_peak_overlap=$(wc -l < ${Overlap_Output_DIR}/$1'.peaks')
#http://unix.stackexchange.com/questions/113795/add-thousands-separator-in-a-number
diff_peak_overlap_Formatted=$(echo ${diff_peak_overlap} | sed ':a;s/\B[0-9]\{3\}\>/,&/;ta')
echo $diff_peak_overlap_Formatted
echo $diff_peak_overlap_Formatted >> ${Overlap_Summary}
diff_peak_overlap_Ratio=$(echo "scale=2;$diff_peak_overlap/$diff_total" | bc)
#echo $diff_peak_overlap_Ratio
diff_peak_overlap_Percent=$(echo "scale=2;$diff_peak_overlap_Ratio*100" | bc)
echo $diff_peak_overlap_Percent'% of the total'
echo $diff_peak_overlap_Percent'% of the total' >> ${Overlap_Summary}
#Want to know how many diffReps sites don't overlap peaks:
#-v	Only report those entries in A that have _no overlaps_ with B.
#		- Similar to "grep -v" (an homage).
bedtools intersect -v -a $Current_DIR/$1'.temp2' -b $Union_DIR/Peak_Union.bed > ${Overlap_Output_DIR}/$1'.nopeaks'
echo 'These are the diffReps sites that DO NOT overlap peaks:'
echo 'These are the diffReps sites that DO NOT overlap peaks:' >> ${Overlap_Summary}
diff_peak_NO_overlap=$(wc -l < ${Overlap_Output_DIR}/$1'.nopeaks')
#http://unix.stackexchange.com/questions/113795/add-thousands-separator-in-a-number
diff_peak_NO_overlap_Formatted=$(echo ${diff_peak_NO_overlap} | sed ':a;s/\B[0-9]\{3\}\>/,&/;ta')
echo $diff_peak_NO_overlap_Formatted
echo $diff_peak_NO_overlap_Formatted >> ${Overlap_Summary}
diff_peak_NO_overlap_Ratio=$(echo "scale=2;$diff_peak_NO_overlap/$diff_total" | bc)
#echo $diff_peak_NO_overlap_Ratio
diff_peak_NO_overlap_Percent=$(echo "scale=2;$diff_peak_NO_overlap_Ratio*100" | bc)
echo $diff_peak_NO_overlap_Percent'% of the total'
echo $diff_peak_NO_overlap_Percent'% of the total' >> ${Overlap_Summary}
#Add header to *.peaks' file
cat $Current_DIR/$1'.header' ${Overlap_Output_DIR}/$1'.peaks' > ${Overlap_Output_DIR}/$1'.temp3'
mv ${Overlap_Output_DIR}/$1'.temp3' ${Overlap_Output_DIR}/$1'.peaks'
#Might as well add header to *'.nopeaks' file
cat $Current_DIR/$1'.header' ${Overlap_Output_DIR}/$1'.nopeaks' > ${Overlap_Output_DIR}/$1'.temp4'
mv ${Overlap_Output_DIR}/$1'.temp4' ${Overlap_Output_DIR}/$1'.nopeaks'
#Remove temp files:
rm $Current_DIR/$1'.temp'*
rm $Current_DIR/$1'.header'
#----------------------------------------------------------------------------------
#I should also parse the *.hotspot file
#I only use the *.hotspot file for generating a UCSC Browser BED file
#I'll hold off on parsing the *.hotspot file
#----------------------------------------------------------------------------------
echo '#---------------------------------------------------------------------------'
echo '#---------------------------------------------------------------------------' >> ${Overlap_Summary}
###################################################################################
#Call the Rscript to parse the *'.peaks' file
echo '############################################################################'
echo 'Running Rscript to get Peaks_Filtered_Summary...'
echo '############################################################################'
echo 'Running Rscript ...'
#Use the system arguments
cd $Current_DIR
#Copy the $1'.peaks' file for the R script:
cp ${Overlap_Output_DIR}/$1'.peaks' $Current_DIR
#----------------------------------------------------------------------------------
Rscript diffReps_Summary.R $1'.peaks' $2 $3 $4 Peaks_Filtered_Summary $5
#----------------------------------------------------------------------------------
rm $Current_DIR/$1'.peaks'
echo '#---------------------------------------------------------------------------'
echo 'Generate BG sites (Peak_Union sites that do not overlap condition-specific sites) ...'
#---------------------------------------------
BG_DIR=$Current_DIR/Peaks_Filtered_Summary/BG_Sites
if [ ! -d $BG_DIR ]
then 
mkdir $BG_DIR
touch $BG_DIR/temp.txt
else
touch $BG_DIR/temp.txt
fi
#---------------------------------------------
diffReps_BED_Files_DIR=$Current_DIR/Peaks_Filtered_Summary/BED_Files
#---------------------------------------------
#The 1st temp file should be the Union peak file (just want to start with a full peak list)
cp $Union_DIR/Peak_Union.bed $BG_DIR
mv $BG_DIR/Peak_Union.bed $BG_DIR/Peak_Union_temp.bed
cd $diffReps_BED_Files_DIR
file_list=*.bed
for file in $file_list
do
echo $file
#---------------------------------------------
#---------------------------------------------
#intersectBed
#-v	Only report those entries in A that have _no overlaps_ with B.
#		- Similar to "grep -v" (an homage).
#---------------------------------------------
intersectBed -v -a $BG_DIR/Peak_Union_temp.bed -b $file >> $BG_DIR/temp.txt
#Take the result of the intersectBed command and use for the next command:
mv $BG_DIR/temp.txt $BG_DIR/Peak_Union_temp.bed
#Sort the Peak_Union_temp.bed before next intersectBed command
sort -k1,1 -k2,2n $BG_DIR/Peak_Union_temp.bed > $BG_DIR/temp2.bed
mv $BG_DIR/temp2.bed $BG_DIR/Peak_Union_temp.bed
done
#---------------------------------------------
#Rename file:
mv $BG_DIR/Peak_Union_temp.bed $BG_DIR/'BG_Sites_'$4'.bed'
#Check that the intersection worked:
echo 'These are the BG_Sites:'
BG_Sites_Count=$(wc -l < $BG_DIR/'BG_Sites_'$4'.bed')
echo $BG_Sites_Count
#The 'BG_Sites_'$4'.bed' should then be used for subsequent enrichment calculations
echo '#---------------------------------------------------------------------------'
echo 'Rename Peaks_Filtered_Summary BED files...'
cd ${diffReps_BED_Files_DIR}
file_list=*.bed
for file in $file_list
do
#echo $file
file_name=${file%\.bed}
#echo $file_name
#Rename the file:
mv ${file} ${file_name}'.peak_filtered.bed'
done
#---------------------------------------------
diffReps_PCA_BED_Files_DIR=$Current_DIR/Peaks_Filtered_Summary/PCA_BED_Files
cd ${diffReps_PCA_BED_Files_DIR}
file_list=*.bed
for file in $file_list
do
#echo $file
file_name=${file%\.bed}
#echo $file_name
#Rename the file:
mv ${file} ${file_name}'.peak_filtered.bed'
done
echo '#---------------------------------------------------------------------------'
echo 'Done!'
echo '############################################################################'
###################################################################################
