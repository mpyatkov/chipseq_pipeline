#!/bin/bash

set -eu

config="./diffreps_config.txt"

while IFS=',' read -r num treatment_name control_name treatment_samples control_samples norm diffreps_window; do
    echo $num
    echo $treatment_name
    echo $control_name
    echo $treatment_samples
    echo $control_samples
    echo $norm
    echo $diffreps_window
    
    template_name="10_diffreps_${num}_${treatment_name}_vs_${control_name}_${norm}"
    rm -rf "../${template_name}"
    
    cp -a ../TEMPLATE_10_diffreps "../${template_name}"

    pushd "../${template_name}"

    sed -i "s/TEMPLATE_CONTROL_NAME/${control_name}/g" Summarize_Jobs.sh
    sed -i "s/TEMPLATE_TREATMENT_NAME/${treatment_name}/g" Summarize_Jobs.sh
    sed -i "s/TEMPLATE_NORMALIZATION/\"${norm}\"/g" Summarize_Jobs.sh
    sed -i "s/TEMPLATE_CONTROL_SAMPLES/\"${control_samples}\"/g" Summarize_Jobs.sh
    sed -i "s/TEMPLATE_TREATMENT_SAMPLES/\"${treatment_samples}\"/g" Summarize_Jobs.sh

    sed -i "s/TEMPLATE_CONTROL_NAME/${control_name}/g" Run_Jobs.sh
    sed -i "s/TEMPLATE_TREATMENT_NAME/${treatment_name}/g" Run_Jobs.sh
    sed -i "s/TEMPLATE_NORMALIZATION/\"${norm}\"/g" Run_Jobs.sh
    sed -i "s/TEMPLATE_COMPARISON_NUMBER/${num}/g" Run_Jobs.sh
    sed -i "s/TEMPLATE_DIFFREPS_WINDOW/${diffreps_window}/g" Run_Jobs.sh
    
    #ctrl_samples=`echo "${control_samples}" | sed -e 's/,/|/g'`
    #treatment_samples=`echo "${treatment_samples}" | sed -e 's/,/|/g'`
    sed -i "s/TEMPLATE_CONTROL_SAMPLES/\"${control_samples}\"/g" Run_Jobs.sh
    sed -i "s/TEMPLATE_TREATMENT_SAMPLES/\"${treatment_samples}\"/g" Run_Jobs.sh
    
    popd
    echo "---"
done < <(cat $config | tr -d " " | sed -e '/^#/d')


    # for control in $(echo $control_samples | tr "," "\n"); do
    #     echo "control: $control"
    # done

    # for treatment in $(echo $treatment_samples | tr "," "\n"); do
    #     echo "treatment: $treatment"
    # done
    


# copy directory with different name
# 10_diffreps_Treatment_vs_Control_DIFFREPS_1
# 10_diffreps_Treatment_vs_Control_MACS2_2

# go inside
# replace TEMPLATE variables to the real one in Run_Jobs.sh
# replace Control_Samples.txt
# replace Treatment_Samples.txt
# create norm file if normalization is MACS2
# go back
