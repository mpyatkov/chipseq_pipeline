#!/bin/bash

set -eu
module load bedops
module load R
module load poppler

source ../00_Setup_Pipeline/01_Pipeline_Setup.sh

rm -rf Job_Summary && mkdir -p Job_Summary

##### Aggregate pdf files
#####
## extract only group names without DIFFREPS and RIPPM suffixes and prefixes
## run only if directory diffreps directories exist

if [[ -n `find ../ -maxdepth 1 -type d | grep -v TEMPLATE|grep 10` ]];
then
    echo "Aggregate pdf files"

    groups=(`find ../ -maxdepth 1 -type d | grep -v TEMPLATE| grep 10 | sed -E 's/10_diffreps_[1-9]+_//g;s/_RIPPM_w[0-9]+//g;s/_DIFFREPS_w[0-9]+//g' | xargs basename -a | sort | uniq| paste -s -d " "`)

    for group in "${groups[@]}"; do
        echo $group
        pdfunite `find ../ -name "*.pdf" | grep -v TEMPLATE | grep ${group} | grep -E "FDR"   | paste -s -d " "` "${group}_FDR_Barcharts_MACS2.pdf"
        pdfunite `find ../ -name "*.pdf" | grep -v TEMPLATE | grep ${group} | grep -E "Hist"  | paste -s -d " "` "${group}_Hist_Barcharts_MACS2.pdf"
        pdfunite `find ../ -name "*.pdf" | grep -v TEMPLATE | grep ${group} | grep -E "\/Bar" | paste -s -d " "` "${group}_Barcharts_MACS2.pdf"

        mkdir -p "Job_Summary/${group}_PDFs" && mv *.pdf "Job_Summary/${group}_PDFs"
    done


    echo "Create diffReps track lines"
    echo
    #### convert bed to bb and create tracks
    function diffreps_bed_to_bb2(){
        local bed_path=$1
        local dataset_label=$2
        local bed=`basename $1`
        
        bb_name=`echo $bed | sed -e 's/UCSC_FILTERED_track_//g;s/\.bed//g'`

        ## swap colors + ignore chrM and header
        # cat $bed | sed -E 's/0,0,255/tmp/g;s/255,0,0/0,0,255/g;s/tmp/255,0,0/g' | grep -vE "track|chrM" > tmp1

        # ignore chrM and header
        cat ${bed_path} | grep -vE "track|chrM" > tmp1
        sort-bed tmp1 > tmp2 && rm tmp1
        ./Job_Scripts/bedToBigBed -allow1bpOverlap tmp2 ./Job_Scripts/mm9.chrom.sizes "${bb_name}.bb" && rm tmp2

        track_line="track type=bigBed name=${bb_name} description=${bb_name} visibility=dense itemRgb=on bigDataUrl=http://waxmanlabvm.bu.edu/mpyatkov/${dataset_label}/${bb_name}.bb"
        echo "${track_line}" > "${bb_name}.bb.ucsc"
    }

    for f in `find ../ -name  "*.bed" | grep 10_ | grep _FILTERED`; do
        echo "Processing: $f"
        diffreps_bed_to_bb2 $f ${Dataset_Label}
    done

    cp *.bb ${VM_DIR_UCSC}
    cat *.ucsc >  ./Job_Summary/all_diffreps_ucsc_tracks.txt
    rm -rf *.bb *.ucsc

    echo "Create aggregated xls files"
    echo
    #### create aggregated xls
    Rscript ./Job_Scripts/diffreps_output_parser.R --path ../
    mv *.xlsx ./Job_Summary
fi


