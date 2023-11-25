#!/bin/bash

module load bedtools/2.27.1

set -eu

source ../00_Setup_Pipeline/01_Pipeline_Setup.sh

Control_Samples_NAME=TEMPLATE_CONTROL_NAME
Treatment_Samples_NAME=TEMPLATE_TREATMENT_NAME
NORMALIZATION=TEMPLATE_NORMALIZATION
CONTROL_SAMPLES=TEMPLATE_CONTROL_SAMPLES
TREATMENT_SAMPLES=TEMPLATE_TREATMENT_SAMPLES
WINDOW_SIZE=TEMPLATE_DIFFREPS_WINDOW

# Control_Samples_NAME=G215K27ac_Female
# Treatment_Samples_NAME=G215K27ac_Male
# NORMALIZATION=DIFFREPS
# CONTROL_SAMPLES="G215M1|G215M2"
# TREATMENT_SAMPLES="G215M3|G215M4"

# TREATMENT_SAMPLES=$(cat Treatment_Samples.txt | grep -v Sample_DIR | cut -f 2 | paste -s -d ",")
# CONTROL_SAMPLES=$(cat Control_Samples.txt | grep -v Sample_DIR | cut -f 2 | paste -s -d ",")
# COMPAR_NUM=3

OUTPUT_DIR="${SCRIPT_DIR}/Output_diffReps"

#Make the rest of the script easier to read: save long output file name as another variable
OUTPUT_NAME="diffReps_${Treatment_Samples_NAME}.vs.${Control_Samples_NAME}"
Count1_Before=$(wc -l < $OUTPUT_DIR/${OUTPUT_NAME})
Count2_Before=$(wc -l < $OUTPUT_DIR/${OUTPUT_NAME}'.annotated')
Count3_Before=$(wc -l < $OUTPUT_DIR/${OUTPUT_NAME}'.hotspot')

function filter_blacklisted_inplace() {
  local bl=$1
  local fn=$2
  local comments=$3

  if [[ $comments == 1 ]]; then
    grep '^#' ${fn} > comments.txt
  fi

  grep -v '^#' ${fn} > coord_with_colnames.txt
  head -1 coord_with_colnames.txt > col_names.txt
  awk 'FNR > 1' coord_with_colnames.txt > coord_without_col_names.txt
  intersectBed -v -a coord_without_col_names.txt -b ${bl} > filtered.txt

  if [[ $comments == 1 ]]; then
      cat comments.txt col_names.txt filtered.txt > "${fn}.filtered"
  else
      cat col_names.txt filtered.txt > "${fn}.filtered"
  fi

  mv "${fn}.filtered" ${fn}

  rm coord*.txt filtered.txt
}

#ENCODE has a list of blacklisted sites that should be removed:
#https://sites.google.com/site/anshulkundaje/projects/blacklists

filter_blacklisted_inplace ./Job_Scripts/ENCODE_Blacklist/mm9-blacklist.bed.gz "${OUTPUT_DIR}/${OUTPUT_NAME}.annotated" 0
filter_blacklisted_inplace ./Job_Scripts/ENCODE_Blacklist/mm9-blacklist.bed.gz "${OUTPUT_DIR}/${OUTPUT_NAME}" 1
filter_blacklisted_inplace ./Job_Scripts/ENCODE_Blacklist/mm9-blacklist.bed.gz "${OUTPUT_DIR}/${OUTPUT_NAME}.hotspot" 1

echo 'AFTER ENCODE_Blacklist filter:'

Count1_Filename=${OUTPUT_NAME}
Count1_After=$(wc -l < $OUTPUT_DIR/${OUTPUT_NAME})
Count1_Diff=$(echo "$Count1_Before - $Count1_After" | bc)
echo "Output file line count: ${Count1_After}"
echo "Number of regions lost: ${Count1_Diff}"

Count2_Filename=${OUTPUT_NAME}'.annotated'
Count2_After=$(wc -l < $OUTPUT_DIR/${OUTPUT_NAME}'.annotated')
Count2_Diff=$(echo "$Count2_Before - $Count2_After" | bc)
echo "Output.annotated file line count: ${Count2_After}"
echo "Number of regions lost: ${Count2_Diff}"

Count3_Filename=${OUTPUT_NAME}'.hotspot'
Count3_After=$(wc -l < $OUTPUT_DIR/${OUTPUT_NAME}'.hotspot')
Count3_Diff=$(echo "$Count3_Before - $Count3_After" | bc)
echo "Output.hotspot file line count: ${Count3_After}"
echo "Number of regions lost: ${Count3_Diff}"

BLACKLIST_OVERLAP_SUMMARY=$OUTPUT_DIR/ENCODE_Blacklist_BLACKLIST_OVERLAP_SUMMARY.txt
rm -rf ${BLACKLIST_OVERLAP_SUMMARY}

echo 'Creating '${BLACKLIST_OVERLAP_SUMMARY}

echo 'Filtering ENCODE blacklisted sites' > ${BLACKLIST_OVERLAP_SUMMARY}
echo Filename $'\t'Before_Filter $'\t'After_Filter $'\t'Sites_Lost >> ${BLACKLIST_OVERLAP_SUMMARY}
echo 'Printing to output file'
echo -e "${Count1_Filename}\t${Count1_Before}\t${Count1_After}\t${Count1_Diff}" >> ${BLACKLIST_OVERLAP_SUMMARY}
echo -e "${Count2_Filename}\t${Count2_Before}\t${Count2_After}\t${Count2_Diff}" >> ${BLACKLIST_OVERLAP_SUMMARY}
echo -e "${Count3_Filename}\t${Count3_Before}\t${Count3_After}\t${Count3_Diff}" >> ${BLACKLIST_OVERLAP_SUMMARY}

echo 'Copy control samples peak BED files'
function copy_sample() {
    
    local samples=("$@");

    for Sample_ID in ${samples[@]}; do
        if [ "${Input_Sites_RiPPM}" == "MACS2" ] ;
        then
            cp "${Dataset_DIR}/${Sample_ID}/fastq/bowtie2/${Sample_ID}_MACS2_output/${Sample_ID}_narrow_MACS2_peaks.xls" ./
        fi

        if [ "${Input_Sites_RiPPM}" == "SICER" ] ;
        then
            # echo 'Note: ${Input_Sites_RiPPM} was set to SICER.'
            # echo 'Therefore, the corresponding SICER BED file will be copied to the Peaks_Called folder.'
            # find ${Dataset_DIR}/${Sample_ID}/fastq/bowtie2/${Sample_ID}_SICER_output -name "*.scoreiland" | xargs -I {} cp {}  ${PEAKS_CALLED_DIR}/{}.bed

            pushd ${Dataset_DIR}/${Sample_ID}/fastq/bowtie2/${Sample_ID}'_SICER_output'
            Peaks_Called_File=$(ls *.scoreisland)
            cp ${Dataset_DIR}/${Sample_ID}/fastq/bowtie2/${Sample_ID}'_SICER_output'/${Peaks_Called_File} ${PEAKS_CALLED_DIR}/${Peaks_Called_File}'.bed'
            popd
        fi
    done
} 

## copy control and treatment samples
control_treatment_samples=($(cat Control_Samples.txt Treatment_Samples.txt | grep -v Sample_DIR | cut -f1 | paste -s -d " "))
mkdir -p ${OUTPUT_DIR}/XLSfiles/
pushd ${OUTPUT_DIR}/XLSfiles/
copy_sample ${control_treatment_samples[@]}
sed -i -E 's/-log10\(/minus_log10_/g;s/value\)/value/g' *.xls
popd

## copy Sample_Labels.txt
cp ../00_Setup_Pipeline/Sample_Labels.txt ${OUTPUT_DIR}
cp ./Job_Scripts/diffreps_genomicRanges.R ${OUTPUT_DIR}

module load R

pushd ${OUTPUT_DIR}
(set -x; Rscript ./diffreps_genomicRanges.R \
                 --annotated_path ${OUTPUT_NAME}.annotated \
                 --hotspot_path ${OUTPUT_NAME}.hotspot \
                 --macs2_xls_dir_path ./XLSfiles/ \
                 --sample_labels_path ./Sample_Labels.txt \
                 --min_avg_count 20 \
                 --log2fc_cutoff 1 \
                 --control_name ${Control_Samples_NAME} \
                 --treatment_name ${Treatment_Samples_NAME} \
                 --peak_caller ${Input_Sites_RiPPM} \
                 --histone_mark ${Dataset_Label} \
                 --normalization_caller "${NORMALIZATION}_${WINDOW_SIZE}" \
                 --treatment_samples ${TREATMENT_SAMPLES} \
                 --control_samples ${CONTROL_SAMPLES})

mkdir -p plots && mv *.pdf ./plots
mkdir -p fullreport && mv *.xlsx ./fullreport
mkdir -p ucsc_tracks && mv *.bed ./ucsc_tracks

rm -rf XLSfiles *.R *.sh Sample_Labels.txt
popd




