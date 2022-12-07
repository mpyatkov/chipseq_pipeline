#!/bin/bash

module load bedtools
module load R

input=combined.sorted.merged.bed

bedtools intersect -wa -wb\
    -a ${input} \
    -b G196_M9_MACS2_peaks.narrowPeak.bed G196_M10_MACS2_peaks.narrowPeak.bed G196_M11_MACS2_peaks.narrowPeak.bed G196_M12_MACS2_peaks.narrowPeak.bed G196_M13_MACS2_peaks.narrowPeak.bed G196_M14_MACS2_peaks.narrowPeak.bed \
    -sorted \
    -names M9 M10 M11 M12 M13 M14 > combined.output

## just aggregate and create xls file
## hardcoded combined.output as input
Rscript ./process_output.R

## almost work (does not remove duplicates)
# cat combined.output | bedtools groupby -i combined.output -g 1-3 -c 4 -o collapse,count > combined.output1

