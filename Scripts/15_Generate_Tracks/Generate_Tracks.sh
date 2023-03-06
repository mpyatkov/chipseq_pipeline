#!/bin/bash
set -eu
##################################################################################
#Andy Rampersaud, 07.19.17
#This script would be used to run Generate_Tracks for different samples
#Way to run script:
#Usage: ./Generate_Tracks.sh
#Example: 
#./Generate_Tracks.sh 
##################################################################################
#---------------------------------------------------------------------------------
#Source the setup file to initialize variables
source ../00_Setup_Pipeline/01_Pipeline_Setup.sh
#---------------------------------------------------------------------------------
#Check that each variable prints a value to the terminal:
echo "-----------------------"
echo "Start of variable list:"
echo "-----------------------"
echo "SCRIPT_DIR:"
echo ${SCRIPT_DIR}
echo "Dataset_Label:"
echo ${Dataset_Label}
echo "BU_User:"
echo ${BU_User}
echo "VM_DIR_UCSC:"
echo ${VM_DIR_UCSC}
echo "-----------------------"
echo "End of variable list"
echo "-----------------------"
#---------------------------------------------------------------------------------
#A text file (Sample_Labels.txt) is needed to run this script
SCRIPT_DIR=$(pwd)
cp $Sample_Labels_DIR/Sample_Labels.txt $SCRIPT_DIR
#---------------------------------------------------------------------------------
#Output dir to store text files:
OUTPUT_DIR=${SCRIPT_DIR}/UCSC_Track_Lines
###############################
if [ ! -d ${OUTPUT_DIR} ]; then
mkdir -p ${OUTPUT_DIR}
else
rm ${OUTPUT_DIR}/*.txt
fi
###############################
#---------------------------------------------------------------------------------
cd ${SCRIPT_DIR}
OUTPUT_FILE_PileUp=${OUTPUT_DIR}/${Dataset_Label}'_Tracks_PileUp.txt'
#---------------------------
if [ -f ${OUTPUT_FILE_PileUp} ]; then
rm ${OUTPUT_FILE_PileUp}; 
fi
#---------------------------
echo
echo 'Start *_Tracks_PileUp.txt'
echo
#---------------------------------------------------------------------------------
#Add TAD regions track line (generated by Bryan):
echo "track type=bigBed name="TADS_ABv4" description="Liver_TADs" visibility=pack itemRgb="On" bigDataUrl=http://waxmanlabvm.bu.edu/aramp10/Lab_Files/TADS_ABv4.bb" >> ${OUTPUT_FILE_PileUp}
#---------------------------------------------------------------------------------
#More than likely for a typical UCSC Session you'll want the chromatin state maps with your screen shots:
echo 'Start chromatin state tracks'
#Add the chromatin state tracks
while IFS=$'\t' read -r -a myArray
do
#---------------------------
#Check that text file is read in properly:
#echo 'Name:'
#echo ${myArray[0]}
#echo 'Description:
#echo ${myArray[1]}
#---------------------------
#visibility options for bb files:
#hide, dense, squish, pack, full
#"pack" works for showing the emission label
#"dense" works for screenshots
#Choose one:
#BB_Visual=hide
BB_Visual=dense
#BB_Visual=squish
#BB_Visual=pack
#BB_Visual=full
#---------------------------
echo 'track type=bigBed name="'${myArray[0]}'" description="'${myArray[1]}'" visibility='${BB_Visual}' itemRgb="On" bigDataUrl=http://waxmanlabvm.bu.edu/aramp10/Lab_Files/'${myArray[0]}'.bb' 
done < BigBed_Name_Color.txt >> ${OUTPUT_FILE_PileUp}
echo 'End chromatin state tracks'
echo 'Add BAM and BigWig file tracks'
################################################
#The text file is formatted like the following:
#----------------------------------------------
#Sample_DIR	Sample_ID	Description	Color
#Sample_Waxman-TP17	G83_M1	Male 8wk-pool 1	0,0,255
#Sample_Waxman-TP18	G83_M2	Male 8wk-pool 2	0,0,255
#Sample_Waxman-TP19	G83_M3	Female 8wk-pool 1	255,0,0
#Sample_Waxman-TP20	G83_M4	Female 8wk-pool 2	255,0,0
#----------------------------------------------
#The 1st column: The Sample_DIR name
#The 2nd column: Waxman Lab Sample_ID 
#The 3rd column: Sample's description 
#The 4th column: BigWig Color
################################################
#Text file has a header line to ignore:
tail -n +2 Sample_Labels.txt > Sample_Labels.temp
#Use a while loop to run jobs
while IFS=$'\t' read -r -a myArray
do
#---------------------------
##Check that text file is read in properly:
#echo 'Sample_DIR:'
Sample_DIR=${myArray[0]}
#echo 'Sample_ID:'
Sample_ID=${myArray[1]}
#echo $Sample_ID
#echo 'Description:'
Description=${myArray[2]}
#echo $Description
Color=${myArray[3]}
#echo $Color
#---------------------------
#visibility options for bw files:
#full, dense, hide
#We typically want to see the read pile ups more than the wiggle tracks:
#Choose one:
BW_Visual=hide
#BW_Visual=dense
#BW_Visual=full
#---------------------------
#Make bigWig tracks (reads):
#---------------------------------------------------------------------------------
#For un-stranded data we get a single *.bw file
#---------------------------------------------------------------------------------
echo 'track type=bigWig name="'${Sample_ID}'_RiPPM_norm" description="'${Description}'" db=mm9 visibility='${BW_Visual}' autoScale=on viewLimits=0.0:100.0 color="'${Color}'" yLineOnOff=off windowingFunction=mean smoothingWindow=3 maxHeightPixels=100:64:8 bigDataUrl=http://waxmanlabvm.bu.edu/'${BU_User}'/'${Dataset_Label}'/'${Sample_ID}'_RiPPM_norm.bw'
#---------------------------------------------------------------------------------
#Make bigBed tracks (peaks):
#macs2:
echo "track type=bigBed name=\"${Sample_ID}_MACS2_narrowPeak_peaks\" description=\"${Sample_ID}_MACS2_narrowPeak_peaks\" visibility=squish color=0,0,0 bigDataUrl=http://waxmanlabvm.bu.edu/${BU_User}/${Dataset_Label}/${Sample_ID}_MACS2_peaks.narrowPeak.bb"
#---------------------------------------------------------------------------------
#SICER:
if [ -d "${Dataset_DIR}/${Sample_ID}/fastq/bowtie2/${Sample_ID}_SICER_output" ]
then
    cd ${Dataset_DIR}/${Sample_ID}/fastq/bowtie2/${Sample_ID}'_SICER_output'
    #General way to refer to the main output file
    Peaks_Called_File=$(ls *.bb)
    echo 'track type=bigBed name="'${Sample_ID}'_SICER" description="'${Sample_ID}'_SICER" visibility=hide color=0,0,0 bigDataUrl=http://waxmanlabvm.bu.edu/'${BU_User}'/'${Dataset_Label}'/'${Peaks_Called_File}
    cd ${SCRIPT_DIR}
fi

#---------------------------------------------------------------------------------
#Make BAM tracks (reads):
#---------------------------
#visibility=<display_mode>
#0 - hide, 1 - dense, 2 - full, 3 - pack, and 4 - squish
#squish looks neater (features on the same horizontal line)
#Read pile ups look better with: pack
#ChIPSeq: #Read pile ups look better with: squish 
#Choose one:
#BAM_Visual=hide
#BAM_Visual=dense
#BAM_Visual=full
#BAM_Visual=pack
BAM_Visual=squish
#---------------------------
echo 'track type=bam name="'${Sample_ID}'_'${Description}'" description="'${Sample_ID}'_'${Description}'" bamColorMode=strand db=mm9 visibility='${BAM_Visual}' bigDataUrl=http://waxmanlabvm.bu.edu/'${BU_User}'/'${Dataset_Label}'/'${Sample_ID}'_sorted_mapped.bam'
done < Sample_Labels.temp >> ${OUTPUT_FILE_PileUp}
echo
echo 'End *_Tracks_PileUp.txt'
echo
##################################################################################
OUTPUT_FILE_Wiggle=${OUTPUT_DIR}/${Dataset_Label}'_Tracks_Wiggle.txt'
#---------------------------
if [ -f ${OUTPUT_FILE_Wiggle} ]; then
rm ${OUTPUT_FILE_Wiggle}; 
fi
#---------------------------
echo
echo 'Start *_Tracks_Wiggle.txt'
echo
#---------------------------------------------------------------------------------
#Add TAD regions track line (generated by Bryan):
echo "track type=bigBed name="TADS_ABv4" description="Liver_TADs" visibility=pack itemRgb="On" bigDataUrl=http://waxmanlabvm.bu.edu/aramp10/Lab_Files/TADS_ABv4.bb" >> ${OUTPUT_FILE_Wiggle}
#---------------------------------------------------------------------------------
#More than likely for a typical UCSC Session you'll want the chromatin state maps with your screen shots:
echo 'Start chromatin state tracks'
#Add the chromatin state tracks
while IFS=$'\t' read -r -a myArray
do
#---------------------------
#Check that text file is read in properly:
#echo 'Name:'
#echo ${myArray[0]}
#echo 'Description:
#echo ${myArray[1]}
#---------------------------
#visibility options for bb files:
#hide, dense, squish, pack, full
#"pack" works for showing the emission label
#"dense" works for screenshots
#Choose one:
#BB_Visual=hide
BB_Visual=dense
#BB_Visual=squish
#BB_Visual=pack
#BB_Visual=full
#---------------------------
echo 'track type=bigBed name="'${myArray[0]}'" description="'${myArray[1]}'" visibility='${BB_Visual}' itemRgb="On" bigDataUrl=http://waxmanlabvm.bu.edu/aramp10/Lab_Files/'${myArray[0]}'.bb' 
done < BigBed_Name_Color.txt >> ${OUTPUT_FILE_Wiggle}
echo 'End chromatin state tracks'
echo 'Add BAM and BigWig file tracks'
################################################
#The text file is formatted like the following:
#----------------------------------------------
#Sample_DIR	Sample_ID	Description	Color
#Sample_Waxman-TP17	G83_M1	Male 8wk-pool 1	0,0,255
#Sample_Waxman-TP18	G83_M2	Male 8wk-pool 2	0,0,255
#Sample_Waxman-TP19	G83_M3	Female 8wk-pool 1	255,0,0
#Sample_Waxman-TP20	G83_M4	Female 8wk-pool 2	255,0,0
#----------------------------------------------
#The 1st column: The Sample_DIR name
#The 2nd column: Waxman Lab Sample_ID 
#The 3rd column: Sample's description 
#The 4th column: BigWig Color
################################################
#Text file has a header line to ignore:
tail -n +2 Sample_Labels.txt > Sample_Labels.temp
#Use a while loop to run jobs
while IFS=$'\t' read -r -a myArray
do
#---------------------------
##Check that text file is read in properly:
#echo 'Sample_DIR:'
Sample_DIR=${myArray[0]}
#echo 'Sample_ID:'
Sample_ID=${myArray[1]}
#echo $Sample_ID
#echo 'Description:'
Description=${myArray[2]}
#echo $Description
Color=${myArray[3]}
#echo $Color
#---------------------------
#visibility options for bw files:
#full, dense, hide
#We typically want to see the read pile ups more than the wiggle tracks:
#Choose one:
#BW_Visual=hide
#BW_Visual=dense
BW_Visual=full
#---------------------------
#Make bigWig tracks (reads):
#---------------------------------------------------------------------------------
#For un-stranded data we get a single *.bw file
#---------------------------------------------------------------------------------
echo 'track type=bigWig name="'${Sample_ID}'_RiPPM_norm" description="'${Description}'" db=mm9 visibility='${BW_Visual}' autoScale=off viewLimits=0.0:100.0 color="'${Color}'" yLineOnOff=off windowingFunction=mean smoothingWindow=3 maxHeightPixels=100:64:8 bigDataUrl=http://waxmanlabvm.bu.edu/'${BU_User}'/'${Dataset_Label}'/'${Sample_ID}'_RiPPM_norm.bw'
#---------------------------------------------------------------------------------
#Make bigBed tracks (peaks):
#macs2:
echo 'track type=bigBed name="'${Sample_ID}'_MACS2_narrowPeak_peaks" description="'${Sample_ID}'_MACS2_narrowPeak_peaks" visibility=squish color=0,0,0 bigDataUrl=http://waxmanlabvm.bu.edu/'${BU_User}'/'${Dataset_Label}'/'${Sample_ID}'_MACS2_peaks.narrowPeak.bb'
#---------------------------------------------------------------------------------
#SICER:

if [ -d "${Dataset_DIR}/${Sample_ID}/fastq/bowtie2/${Sample_ID}_SICER_output" ]
then

    cd ${Dataset_DIR}/${Sample_ID}/fastq/bowtie2/${Sample_ID}'_SICER_output'
    #General way to refer to the main output file
    Peaks_Called_File=$(ls *.bb)
    echo 'track type=bigBed name="'${Sample_ID}'_SICER" description="'${Sample_ID}'_SICER" visibility=hide color=0,0,0 bigDataUrl=http://waxmanlabvm.bu.edu/'${BU_User}'/'${Dataset_Label}'/'${Peaks_Called_File}
    cd ${SCRIPT_DIR}
fi

#---------------------------------------------------------------------------------
#Make BAM tracks (reads):
#---------------------------
#visibility=<display_mode>
#0 - hide, 1 - dense, 2 - full, 3 - pack, and 4 - squish
#squish looks neater (features on the same horizontal line)
#Read pile ups look better with: pack
#Choose one:
BAM_Visual=hide
#BAM_Visual=dense
#BAM_Visual=full
#BAM_Visual=pack
#BAM_Visual=squish
#---------------------------
echo 'track type=bam name="'${Sample_ID}'_'${Description}'" description="'${Sample_ID}'_'${Description}'" bamColorMode=strand db=mm9 visibility='${BAM_Visual}' bigDataUrl=http://waxmanlabvm.bu.edu/'${BU_User}'/'${Dataset_Label}'/'${Sample_ID}'_sorted_mapped.bam'
done < Sample_Labels.temp >> ${OUTPUT_FILE_Wiggle}
echo
echo 'End *_Tracks_Wiggle.txt'
echo
##################################################################################
OUTPUT_FILE_Wiggle_autoON=${OUTPUT_DIR}/${Dataset_Label}'_Tracks_Wiggle_autoON.txt'
#---------------------------
if [ -f ${OUTPUT_FILE_Wiggle_autoON} ]; then
rm ${OUTPUT_FILE_Wiggle_autoON}; 
fi
#---------------------------
echo
echo 'Start *_Tracks_Wiggle_autoON.txt'
echo
#---------------------------------------------------------------------------------
#Add TAD regions track line (generated by Bryan):
echo "track type=bigBed name="TADS_ABv4" description="Liver_TADs" visibility=pack itemRgb="On" bigDataUrl=http://waxmanlabvm.bu.edu/aramp10/Lab_Files/TADS_ABv4.bb" >> ${OUTPUT_FILE_Wiggle_autoON}
#---------------------------------------------------------------------------------
#More than likely for a typical UCSC Session you'll want the chromatin state maps with your screen shots:
echo 'Start chromatin state tracks'
#Add the chromatin state tracks
while IFS=$'\t' read -r -a myArray
do
#---------------------------
#Check that text file is read in properly:
#echo 'Name:'
#echo ${myArray[0]}
#echo 'Description:
#echo ${myArray[1]}
#---------------------------
#visibility options for bb files:
#hide, dense, squish, pack, full
#"pack" works for showing the emission label
#"dense" works for screenshots
#Choose one:
#BB_Visual=hide
BB_Visual=dense
#BB_Visual=squish
#BB_Visual=pack
#BB_Visual=full
#---------------------------
echo 'track type=bigBed name="'${myArray[0]}'" description="'${myArray[1]}'" visibility='${BB_Visual}' itemRgb="On" bigDataUrl=http://waxmanlabvm.bu.edu/aramp10/Lab_Files/'${myArray[0]}'.bb' 
done < BigBed_Name_Color.txt >> ${OUTPUT_FILE_Wiggle_autoON}
echo 'End chromatin state tracks'
echo 'Add BAM and BigWig file tracks'
################################################
#The text file is formatted like the following:
#----------------------------------------------
#Sample_DIR	Sample_ID	Description	Color
#Sample_Waxman-TP17	G83_M1	Male 8wk-pool 1	0,0,255
#Sample_Waxman-TP18	G83_M2	Male 8wk-pool 2	0,0,255
#Sample_Waxman-TP19	G83_M3	Female 8wk-pool 1	255,0,0
#Sample_Waxman-TP20	G83_M4	Female 8wk-pool 2	255,0,0
#----------------------------------------------
#The 1st column: The Sample_DIR name
#The 2nd column: Waxman Lab Sample_ID 
#The 3rd column: Sample's description 
#The 4th column: BigWig Color
################################################
#Text file has a header line to ignore:
tail -n +2 Sample_Labels.txt > Sample_Labels.temp
#Use a while loop to run jobs
while IFS=$'\t' read -r -a myArray
do
#---------------------------
##Check that text file is read in properly:
#echo 'Sample_DIR:'
Sample_DIR=${myArray[0]}
#echo 'Sample_ID:'
Sample_ID=${myArray[1]}
#echo $Sample_ID
#echo 'Description:'
Description=${myArray[2]}
#echo $Description
Color=${myArray[3]}
#echo $Color
#---------------------------
#visibility options for bw files:
#full, dense, hide
#We typically want to see the read pile ups more than the wiggle tracks:
#Choose one:
#BW_Visual=hide
#BW_Visual=dense
BW_Visual=full
#---------------------------
#Make bigWig tracks (reads):
#---------------------------------------------------------------------------------
#For un-stranded data we get a single *.bw file
#---------------------------------------------------------------------------------
echo 'track type=bigWig name="'${Sample_ID}'_RiPPM_norm" description="'${Description}'" db=mm9 visibility='${BW_Visual}' autoScale=on viewLimits=0.0:100.0 color="'${Color}'" yLineOnOff=off windowingFunction=mean smoothingWindow=3 maxHeightPixels=100:64:8 bigDataUrl=http://waxmanlabvm.bu.edu/'${BU_User}'/'${Dataset_Label}'/'${Sample_ID}'_RiPPM_norm.bw'
#---------------------------------------------------------------------------------
#Make bigBed tracks (peaks):
#macs2:
echo 'track type=bigBed name="'${Sample_ID}'_MACS2_narrowPeak_peaks" description="'${Sample_ID}'_MACS2_narrowPeak_peaks" visibility=squish color=0,0,0 bigDataUrl=http://waxmanlabvm.bu.edu/'${BU_User}'/'${Dataset_Label}'/'${Sample_ID}'_MACS2_peaks.narrowPeak.bb'
#---------------------------------------------------------------------------------
#SICER:
if [ -d "${Dataset_DIR}/${Sample_ID}/fastq/bowtie2/${Sample_ID}_SICER_output" ]
then
    cd ${Dataset_DIR}/${Sample_ID}/fastq/bowtie2/${Sample_ID}'_SICER_output'
    #General way to refer to the main output file
    Peaks_Called_File=$(ls *.bb)
    echo 'track type=bigBed name="'${Sample_ID}'_SICER" description="'${Sample_ID}'_SICER" visibility=hide color=0,0,0 bigDataUrl=http://waxmanlabvm.bu.edu/'${BU_User}'/'${Dataset_Label}'/'${Peaks_Called_File}
    cd ${SCRIPT_DIR}
fi

#---------------------------------------------------------------------------------
#Make BAM tracks (reads):
#---------------------------
#visibility=<display_mode>
#0 - hide, 1 - dense, 2 - full, 3 - pack, and 4 - squish
#squish looks neater (features on the same horizontal line)
#Read pile ups look better with: pack
#Choose one:
BAM_Visual=hide
#BAM_Visual=dense
#BAM_Visual=full
#BAM_Visual=pack
#BAM_Visual=squish
#---------------------------
echo 'track type=bam name="'${Sample_ID}'_'${Description}'" description="'${Sample_ID}'_'${Description}'" bamColorMode=strand db=mm9 visibility='${BAM_Visual}' bigDataUrl=http://waxmanlabvm.bu.edu/'${BU_User}'/'${Dataset_Label}'/'${Sample_ID}'_sorted_mapped.bam'
done < Sample_Labels.temp >> ${OUTPUT_FILE_Wiggle_autoON}
echo
echo 'End *_Tracks_Wiggle_autoON.txt'
echo
##################################################################################
#Remove the temp file:
rm *.temp
echo '#---------------------------------------------------------------------------'
echo 'Check out '${OUTPUT_FILE_PileUp}
echo 'Check out '${OUTPUT_FILE_Wiggle}
echo 'Check out '${OUTPUT_FILE_Wiggle_autoON}
echo '#---------------------------------------------------------------------------'
##################################################################################

cp -a UCSC_Track_Lines ${VM_DIR_UCSC}/
