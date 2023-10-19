#!/bin/bash


source ../00_Setup_Pipeline/01_Pipeline_Setup.sh

set -eu

rm -rf Job_Summary && mkdir -p Job_Summary

DIR_name=`basename $PWD`
Step_Num=$(echo ${DIR_name} | cut -d'_' -f 1)
Job_Name="Step_${Step_Num}_${Dataset_Label}"

(set -x; qsub -N "${Job_Name}" -P wax-es fastq_reads.qsub ${Dataset_DIR})
