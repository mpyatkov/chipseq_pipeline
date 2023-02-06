#!/bin/bash

## Find overlapping peaks 

module load bedtools
module load R/4.1.2

find ../G* -name "*narrowPeak*bed" | grep -v Scripts | grep -v blacklist | grep -v StrgtPksRm | xargs -n1 -I{} cp {} ./

## combine all narrowPeaks inside current directory in one file
cat *.bed > combined.bed

sort -k1,1 -k2,2n combined.bed > combined.sorted.bed
bedtools merge -i combined.sorted.bed > combined.sorted.merged.bed

input_files=`find . -name "*narrowPeak.bed" | xargs -n1 basename |sort | paste -s -d " "`
input_names=`find . -name "*narrowPeak.bed" | xargs -n1 basename |sort | grep -E -o 'M[0-9]+' | paste -s -d " "`
bedtools intersect -wa -wb\
    -a combined.sorted.merged.bed \
    -b ${input_files} \
    -sorted \
    -names ${input_names} > combined.output

## just aggregate and create xls file
## hardcoded combined.output as input
Rscript ./process_output.R

## almost work (does not remove duplicates)
# cat combined.output | bedtools groupby -i combined.output -g 1-3 -c 4 -o collapse,count > combined.output1
