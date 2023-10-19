#!/bin/bash

set -eu

source ../00_Setup_Pipeline/01_Pipeline_Setup.sh

module load multiqc 
multiqc ${Dataset_DIR} -o "Job_Summary/${Dataset_Label}_MULTIQC"

## copy to server
(set -x; cp -a "Job_Summary/${Dataset_Label}_MULTIQC" ${VM_DIR_FASTQC})

