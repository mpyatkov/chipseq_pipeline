#! /bin/bash
set -eu
##################################################################################
# Andy Rampersaud, 11.20.2017
#The following code is used to convert BAM files to BigWig files for visualization on the UCSC Browser
#------------------------------------------------------------
#The current procedure for conversion involves the following:
# a. Convert bam to bedGraph to bigWig
# b. bam to bedGraph using BedTools - genomeCoverageBed
# c. bedGraph to bigWig using UCSC - bedGraphToBigWig program
#------------------------------------------------------------
#Command to run this script:
#time ./BAM_to_WIG.sh
#This script needs to be executed in the parent folder containing the following folders:
#-----------
#genomeIndex
#Input_files
#Perl_Scripts
#UCSC_Tools
#-----------
#genomeIndex: contains mm9.chrom.sizes file
#Input_files: contain the BAM files you want to convert
#UCSC_Tools: contain the UCSC executable programs (bedGraphToBigWig and bigWigToWig)
##################################################################################
#Need to run this shell script with a system argument
if [ $# -ne 1 ]
then
    echo "Usage: `basename $0` <sample_norm_factor>"
    exit 
fi
sample_norm_factor=$1
echo
echo 'Value of sample_norm_factor:'
echo $sample_norm_factor
echo
#Decompress the read files:
echo 'Start decompressing BAM files'
cd Input_files
#Check if any compressed files exist:
count=`ls -1 *.gz 2>/dev/null  | wc -l`
if [ $count != 0 ]
then 
    gunzip *.gz
fi 
cd ..
echo 'Done decompressing BAM files'
#--------------------------------------------------------------------------------
echo 'Start pre-process input BAM file'
cd Input_files
Input_list=*.bam
for file in $Input_list
do
    Sample_ID=${file%\_sorted_mapped.bam};
    echo 'Start read count'
    #Get per chromosome count (need to eliminate chrM/chrRandom)
    #Need the BAM file index before doing idxstats command:
    #Generate BAM index
    samtools index ${file}
    samtools idxstats ${file} > ${Sample_ID}'_idxstats.out'
    #Need to modify idxstats.out file to exclude chr*_random and chrM (and * (has a zero read count))
    sed -i '/random/d;/chrM/d;/*/d' ${Sample_ID}'_idxstats.out'
    #--------------------------------------------------------------------------------
done
echo 'End pre-process input BAM file'
#--------------------------------------------------------------------------------
cd ..
#Now I have the following in the BAM_to_WIG folder:
#genomeIndex  Input_files  UCSC_Tools
#I still need some temp folders to hold the intermediary files
mkdir temp_bedGraph
mkdir temp_bedGraph_sorted
mkdir temp_bw
#Generate list of input files, for each file execute genomeCoverageBed and put output file to temp_bedGraph folder
echo 'Started genomeCoverageBed'
cd Input_files
Input_list=*.bam
for file in $Input_list
do
    Sample_ID=${file%\_sorted_mapped.bam};
    #http://bedtools.readthedocs.io/en/latest/content/tools/genomecov.html
    #Use the (-i) option for BED/GFF/VCF files
    #Use the (-ibam) option for BAM files
    genomeCoverageBed -bg -ibam $file -g ../genomeIndex/mm9.chrom.sizes > ../temp_bedGraph/${Sample_ID}'.bedGraph'
done
cd ..
echo 'Finished genomeCoverageBed'
#The above works for generating bedGraph files for each input file
#Now I need to sort the bedGraph files and put into temp_bedGraph_sorted folder:
echo 'Started sorting bedGraph files'
cd temp_bedGraph
Input_list=*.bedGraph
for file in $Input_list
do
    Sample_ID=${file%\.bedGraph};
    #Filter command (want to use the -i option):
    #--------------------------------------------------------------------------------
    #-i[SUFFIX], --in-place[=SUFFIX]
    #edit files in place (makes backup if extension supplied). The default operation mode is to break symbolic and hard links. This can be changed with --follow-symlinks and --copy.
    #--------------------------------------------------------------------------------
    #The (-i) option will avoid creating a temporary file (takes up space)
    sed -i '/track/d;/random/d;/chrM/d' ${Sample_ID}'.bedGraph'
    #--------------------------------------------------------------------------------
    sort -k1,1 -k2,2n $file > ../temp_bedGraph_sorted/${Sample_ID}'_sorted.bedGraph'
done
#--------------------------------------------------------------------------------
cd ..
echo 'Finished sorting bedGraph files'
#Above loop works for sorting all input files
#Need to call normalize_read_counts.pl for BedGraph normalization
echo 'Started bedGraph normalization'
cd temp_bedGraph_sorted
########################################################
#---------------------------------------------------------------------------------
#If you are using the RiPPM_Factor to normalize raw counts (use division):
#	1. Peak_Union_Count_Stats.Rmd: Normalized counts = (raw count) / (RiPPM_Factor)
#	2. BAM_to_WIG.sh: $4/sample_norm_factor
#If you are using the Norm_Factor (need to multiply rather than divide):
#	1. Peak_Union_Count_Stats.Rmd: Normalized counts = (raw count) * (Norm_Factor)
#	2. BAM_to_WIG.sh: $4 * sample_norm_factor
#---------------------------------------------------------------------------------
echo
#echo 'Using awk to divide the 4th column of the bedgraph file by the norm_factor'
echo 'Using awk to multiply the 4th column of the bedgraph file by the norm_factor'
echo
########################################################
#Directory to store normalized file:
mkdir norm_bedGraph
#Awk arithmetic:
#https://www.gnu.org/software/gawk/manual/html_node/Arithmetic-Ops.html
#Need the (-v) for each variable and enclose in double quotes
#time awk -v OFS='\t' -v sample_norm_factor="${sample_norm_factor}" '{print $1, $2, $3, $4/sample_norm_factor}' ${Sample_ID}'_sorted.bedGraph' > norm_bedGraph/${Sample_ID}'_sorted.bedGraph'
time awk -v OFS='\t' -v sample_norm_factor="${sample_norm_factor}" '{print $1, $2, $3, $4 * sample_norm_factor}' ${Sample_ID}'_sorted.bedGraph' > norm_bedGraph/${Sample_ID}'_sorted.bedGraph'
########################################################
echo
echo 'Confirm that the awk command is normalizing read counts:'
echo 'Sample of bedGraph_sorted file (non-normalized input to the awk command)'
head -10 *.bedGraph
echo 'Sample of the norm_bedGraph (normalized output of the awk command)'
head -10 norm_bedGraph/*.bedGraph
cd ..
echo 'Finished bedGraph normalization'
#Now I  need to run bedGraphToBigWig and send the output to temp_bw folder
echo 'Started bedGraphToBigWig'
cd temp_bedGraph_sorted/norm_bedGraph
Input_list=*_sorted.bedGraph
for file in $Input_list
do
    Sample_ID=${file%\_sorted.bedGraph};
    #Add "RiPPM" label to indicate normalization method
    ../../UCSC_Tools/bedGraphToBigWig $file ../../genomeIndex/mm9.chrom.sizes ../../temp_bw/${Sample_ID}'_RiPPM_norm.bw'
done
cd ../../
echo 'Finished bedGraphToBigWig'
#Above loop works
#echo 'Start Compress bigWig file'
#cd temp_bw
#gzip *.bw
#cd ..
#echo 'End Compress bigWig file'
echo 'Removing temporary files and folders...'
rm -r temp_bedGraph
rm -r temp_bedGraph_sorted
echo 'Finished BAM_to_WIG conversion!'
##################################################################################
