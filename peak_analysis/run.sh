#!/bin/bash

set -eu
## Find overlapping peaks 

module load bedtools
module load R

# echo "Clean all bed files inside current directory"
# rm -rf *.bed
rm -rf combined.* *.xlsx

# echo "Copying files"
find ../G* -name "*narrowPeak*bed" | grep -v Scripts | grep -v blacklist | grep -v StrgtPksRm | xargs -n1 -I{} cp {} ./

# echo "Copying xls"
find ../G* -name "*xls" | grep -v Scripts | grep -v blacklist | grep -v StrgtPksRm | xargs -n1 -I{} cp {} ./

echo "Create combined.bed"
## combine all narrowPeaks inside current directory in one file
cat *.bed > combined.bed

echo "Sort combined.bed to combined.sorted.bed"
sort -k1,1 -k2,2n combined.bed > combined.sorted.bed

echo "Merge combined.sorted.bed to combined.sorted.merged.bed"
bedtools merge -i combined.sorted.bed > combined.sorted.merged.bed

echo "Make intersection of all bed files with combined.sorted.merged.bed"
input_files=`find . -maxdepth 1 -name "*narrowPeak.bed" | xargs -n1 basename |sort | paste -s -d " "`

## only M numbers. ex. M9,M10,... (convinient when need short names)
# input_names=`find . -name "*narrowPeak.bed" | xargs -n1 basename |sort | grep -E -o 'M[0-9]+' | paste -s -d " "`

## GM numbers. Ex. G196_M9, G196_M10 (more convinient if Sample_Labels.txt will be parsed)
input_names=`find . -maxdepth 1 -name "*narrowPeak.bed" | xargs -n1 basename |sort | grep -E -o 'G[0-9]+_?M[0-9]+M?[0-9]+?' | paste -s -d " "`

bedtools intersect -wa -wb\
    -a combined.sorted.merged.bed \
    -b ${input_files} \
    -sorted \
    -names ${input_names} > combined.output

## just aggregate and create xls file
## hardcoded combined.output as input
echo "Run process_output.R"
Rscript ./process_output_compact.R

echo "Run process_output_detailed.R"
Rscript ./process_output_detailed.R

echo "Run process_output_extradetailed.R"
Rscript ./process_output_extradetailed.R
## almost work (does not remove duplicates)
# cat combined.output | bedtools groupby -i combined.output -g 1-3 -c 4 -o collapse,count > combined.output1
