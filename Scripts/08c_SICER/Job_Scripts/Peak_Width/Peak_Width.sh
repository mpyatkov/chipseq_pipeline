#! /bin/bash

set -eu
module load R/3.6.2
#This script is used to summarize the peak width distribution for a single BED file of peaks
#Sample command: 
#./Peak_Width.sh ./G90_M1/G90_M1_MACS2_peaks.bed
#OR
#./Peak_Width.sh /home/aramp10/Desktop/Peak_Width_Test/G90_M1/G90_M1_MACS2_peaks.bed
if [ $# -ne 1 ]
then
  echo "Usage: `basename $0` peak1.bed"
  exit 
fi
# echo $1
peak1_base=`basename $1`
#echo $peak1_base
peak1_name=${peak1_base%\.bed};
echo $peak1_name
echo "Step1: clean input"
#sed 's/\s$//g' $1 | awk 'BEGIN {OFS="\t"}
awk 'BEGIN {OFS="\t"} 
     {if ($1~/chr/ && $1 !="chrM" && $1 !~/random/ && $3>$2 && $2>0 && $3>0)
          print $1,$2,$3,$4,$5>"peak1.bed";
      else 
          print $0 > "peak1_dump.bed"}' $1
echo "Step2: Running Peak_Width_Stats.R script to get peak width statistics ..."
Rscript Peak_Width_Stats.R $peak1_name
rm peak1*
echo "Step3: Done!"
echo
