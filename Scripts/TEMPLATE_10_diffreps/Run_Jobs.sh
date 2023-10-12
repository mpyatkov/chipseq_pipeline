#!/bin/bash
set -eu

rm -rf *.o* *.po*

source ../00_Setup_Pipeline/01_Pipeline_Setup.sh

Control_Samples_NAME=TEMPLATE_CONTROL_NAME
Treatment_Samples_NAME=TEMPLATE_TREATMENT_NAME
COMPAR_NUM=TEMPLATE_COMPARISON_NUMBER
WINDOW_SIZE=TEMPLATE_DIFFREPS_WINDOW
NORMALIZATION=TEMPLATE_NORMALIZATION
FRAG_SIZE=0

CONTROL_SAMPLES=TEMPLATE_CONTROL_SAMPLES
TREATMENT_SAMPLES=TEMPLATE_TREATMENT_SAMPLES

## create control and treatment samples files
cat ../00_Setup_Pipeline/Sample_Labels.txt | grep -E "Sample_DIR|${CONTROL_SAMPLES}" > Control_Samples.txt
cat ../00_Setup_Pipeline/Sample_Labels.txt | grep -E "Sample_DIR|${TREATMENT_SAMPLES}" > Treatment_Samples.txt

## if normalization is MACS2 than norm.txt should be present in this
## directory

if [ "$NORMALIZATION" = "RIPPM" ]; then
    
    #Extract normalization factors and represent them in the 1/factor
    # form, because diffreps works only with such representation
    
    control_norm_line=$(cat ../09_RiPPM_Normalization/Job_Summary/Norm_Factors.txt | \
	grep -v FRAG | \
	cut -f1,4 | \
	awk '{printf "%s %.2f\n", $1, 1/$2}' | \
	grep -E "${CONTROL_SAMPLES}" | \
	cut -d " " -f2 | \
	paste -s -d " ")

    treatment_norm_line=$(cat ../09_RiPPM_Normalization/Job_Summary/Norm_Factors.txt | \
	grep -v FRAG | \
	cut -f1,4 | \
	awk '{printf "%s %.2f\n", $1, 1/$2}' | \
	grep -E "${TREATMENT_SAMPLES}" | \
	cut -d " " -f2 | \
	paste -s -d " ")

    echo "treatment ${treatment_norm_line}" > norm.txt
    echo "control ${control_norm_line}" >> norm.txt
fi

DIR_name=`basename ${SCRIPT_DIR}`
Step_Num=$(echo ${DIR_name} | cut -d'_' -f 1)
Job_Name="Step_${Step_Num}_${Dataset_Label}"

(set -x; qsub -N ${Job_Name}'_diffReps_'${COMPAR_NUM} -P wax-es -q !linga diffReps.qsub ${Dataset_DIR} ${Dataset_Label} ${Control_Samples_NAME} ${Treatment_Samples_NAME} ${COMPAR_NUM} ${WINDOW_SIZE} ${FRAG_SIZE} ${Input_Sites_RiPPM} ${SCRIPT_DIR})
