#! /bin/bash
set -eu
#Andy Rampersaud
#modified by Amy Li
#08.19.13
#This script is paired with the Perl script (RemoveStraightPeak.pl)
#The Perl script was originally provided from Yijing
#I modified the script so that it works for files from my pipeline
#This script is supposed to pre-process BED files for input into the RemoveStraightPeak.pl script
#Ideally the input BED file would be produced from a similar command:
#bamToBed -i 'G92_M'$M_number'_sorted_mapped.bam' > 'G92_M'$M_number'.bed'
#Where the sorting and mapped read extraction happens with the BAM file
#Example command to run this script:
#time ./remove_straight_peaks.sh
Input_list=*.bed
#For all BED files in the current dir:
for file in $Input_list
do
file_name=${file%\.bed};
echo
echo 'Processing file: '$file
echo
#Input BED file of peaks should be a BED file with the standard 9 columns
#The 9th column indicates the peak color
Pos_Read_Color="0,0,128"
Neg_Read_Color="255,0,0"
awk '{if ($6 == "+") print $1"\t"$2"\t"$3"\t"$4"\t"1000"\t"$6"\t"0"\t"0"\t""'$Pos_Read_Color'"}' $file > temp1.bed
awk '{if ($6 == "-") print $1"\t"$2"\t"$3"\t"$4"\t"1000"\t"$6"\t"0"\t"0"\t""'$Neg_Read_Color'"}' $file > temp2.bed
cat temp1.bed temp2.bed > temp3.bed
#Filter our track, chrRandom, chrM lines
#grep -v "track" temp3.bed | grep -v "random" | grep -v "chrM" > temp4.bed
grep -v "track" temp3.bed > temp4.bed
grep -v "random" temp4.bed > temp5.bed 
grep -v "chrM" temp5.bed > temp6.bed
#Sort BED file
sort -k1,1 -k2,2n temp6.bed > temp7.bed
mv temp7.bed $file_name'_sort_tab.bed'
## KK added a perl line to remove the TAB in the eol
sed 's/\t$//' $file_name'_sort_tab.bed' > $file_name'_sort.bed'
rm temp*.bed
rm $file_name'_sort_tab.bed'
#Need to get the read length from the input BED file to pipe into Perl script
Read_End=$(awk '{print $3}' $file_name'_sort.bed' | head -1)
Read_Start=$(awk '{print $2}' $file_name'_sort.bed' | head -1)
Read_Length=$(echo "$Read_End - $Read_Start" | bc)
echo
echo 'Running RemoveStraightPeak.py '$Read_Length
echo
#Perl script command:
time python remove_straight_peaks.py $Read_Length
echo
echo 'RemoveStraightPeak.py Complete'
echo
#End of for loop
done
