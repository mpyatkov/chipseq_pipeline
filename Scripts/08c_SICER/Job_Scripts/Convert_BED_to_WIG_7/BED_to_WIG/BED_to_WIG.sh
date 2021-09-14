#! /bin/bash
set -eu
##########################################################################################
#The following code is used to convert BED files to BigBed files for visualization on the UCSC Browser

#The current procedure for conversion involves the following:
# a. Pre-process the MACS BED files
# b. bed to bigBed using UCSC - bedToBigBed

#Modifiy permissions of this bash script:
#chmod a+x BED_to_WIG.sh
#Command to run this script:
#time ./BED_to_WIG.sh
#This script needs to be executed in the parent folder containing the following folders:
#genomeIndex
#Input_files
#Input_files_2
#Perl_scripts
#UCSC_Tools

# mkdir Input_files
#Transfer BED files you want to convert to the Input_files folder

#Decompress the input BED files:
echo 'Start decompressing BED files'
cd Input_files
#Check if any compressed files exist:
count=`ls -1 *.gz 2>/dev/null  | wc -l`
if [ $count != 0 ]
then 
gunzip *.gz
fi 
cd ..
echo 'Done decompressing BED files'
#Now I have the following in the BED_to_WIG folder:
#genomeIndex  Input_files  UCSC_Tools
#I still need some temp folders to hold the intermediary files

mkdir temp_bb

#Generate list of input files, perform preprocessing for each file
echo 'Started Pre-processing'
cd Input_files
Input_list=*.bed
for file in $Input_list
do
file_name=${file%\.bed};
#Only want the first 4 columns of MACS BED file
#MACS format:
# chr1	3406014	3406713	MACS_peak_1	83.92
# chr1	3466424	3468614	MACS_peak_2	355.79
# chr1	3504190	3505963	MACS_peak_3	2755.71
# chr1	3541097	3542564	MACS_peak_4	144.58
#UCSC only accept socres 0 and 1000 but MACS scores can be > 1000
awk '{print $1"\t"$2"\t"$3"\t"$4}' $file > $file_name.temp1
#Remove track lines
grep -v "track" $file_name.temp1 > $file_name.temp2
#MACS2 peaks could contain chr13_random entries
grep -v "random" $file_name.temp2 > $file_name.temp3
grep -v "chrM" $file_name.temp3 > $file_name.temp4
#Need to sort BED file
sort -k1,1 -k2,2n $file_name.temp4 > $file_name.temp5
mv $file_name.temp5 $file_name.bed
rm $file_name.temp*
done
#Need to remove peaks outside of chromosome size limits
cp ../Perl_scripts/Remove_Outbound_Peaks.pl .
#The Remove_Outbound_Peaks.pl operates on all BED files in the Input_files dir
#The Remove_Outbound_Peaks.pl creates *_Inbound_peaks.bed in the Input_files_2 dir     
perl Remove_Outbound_Peaks.pl
echo
echo 'Line count of peak BED files (before Remove_Outbound_Peaks.pl):'
echo
wc -l *.bed
echo
cd ..
echo 'Finished Pre-processing'

#Now I  need to run bedToBigBed and send the output to temp_bb folder

echo 'Started bedToBigBed'
cd Input_files_2
echo
echo 'Line count of peak BED files (after Remove_Outbound_Peaks.pl):'
echo
wc -l *.bed
echo
Input_list=*_Inbound_peaks.bed
for file in $Input_list
do
file_name=${file%\_Inbound_peaks.bed};
../UCSC_Tools/bedToBigBed $file ../genomeIndex/mm9.chrom.sizes ../temp_bb/$file_name'.bb'
done
cd ..
echo 'Finished bedToBigBed'
#Above loop works

echo 'Finished BED_to_WIG conversion!'
##########################################################################################
