##################################################################################
#Andy Rampersaud, 07.19.17
#This 01_Pipeline_Setup.sh file is used to initialize all variables used in the pipeline
#The purpose of this file is to streamline the pipeline so that 
#	1. Variables are initialized from a single file
#	2. Minimize redundancy between pipeline steps
#This file is organized by 
#	1. "Global Variables"  		= variables used by multiple steps
#	2. "Step-specific Variables"  	= variables used by a particular step
##################################################################################
#---------------------------------------------------------------------------------
#"Global Variables"  = variables used by multiple steps
#---------------------------------------------------------------------------------
##changed the directory path, BU_User and Dataset_Label
Dataset_DIR=/projectnb/wax-es/test_DNAChIP
BU_User="mpyatkov"
Dataset_Label="G182"
GTF_Files_DIR=/projectnb/wax-es/aramp10/GTF_Files
#---------------------------------------------------------------------------------
##################################################################################
#---------------------------------------------------------------------------------
#"Step-specific Variables"  = variables used by a particular step
#---------------------------------------------------------------------------------
#Trim_Galore
#Running the Trim_Galore step is optional
#This job is useful for identifying issues with adaptor sequence
#The read length is used in the following way:
#Discard reads that became shorter than length INT because of either quality or adapter trimming.
#Input the read length (bp)
#Should have this from the Read_Length job
READ_LEN=150
#---------------------------------------------------------------------------------
#Bowtie2
Bowtie2Index_DIR=/projectnb/wax-es/aramp10/Bowtie2
#---------------------------------------------------------------------------------
#MACS2
#Is the pipeline being run on DHS data? ("YES" or "NO")
DHS_DATA="YES"
#---------------------------------------------------------------------------------
#RiPPM_Normalization
#Which type of sites should be used as input for the RiPPM normalization?
#MACS2 is used for narrow genomic features:
#	For example: DHS, TFBS, chromatin marks: H3K4me1, H3K4me3, H3K27ac, etc...
#SICER is used for broad genomic features:
#	For example: chromatin marks: H3K27me3, H3K36me3, etc...
#Choose one: ("MACS2" or "SICER")
Input_Sites_RiPPM="MACS2"
#---------------------------------------------------------------------------------
#MEME_ChIP
#This needs to be the dir location that contains the mouse_genome.fa.gz file
Input_FASTA_DIR=/projectnb/wax-es/aramp10/Motif_Files
#Sequence length
#Using a value of 500 will create 500 bp wide sites from the midpoint of each BED record 
SEQ_LEN=500
#---------------------------------------------------------------------------------
#meme-chip option: 
#MEME Specific Options:
#Options explained here:
#http://bioweb2.pasteur.fr/docs/modules/meme/4.10.1/
#Contributing Site Distribution
#This option is used to describe the distribution of motif sites. 
#-meme-mod [oops|zoops|anr]: sites used in a single sequence
#--------------------------------------------------------------------------------
#Choose one:
#0="oops"
#1="zoops"
#2="anr"
MEME_MOD="2"
#--------------------------------------------------------------------------------
#diffReps_1
#Note: 
#	Modify setup_diffReps.sh (located in job-specific folder) as needed
#	Modify *_Samples.txt files as needed
#---------------------------------------------------------------------------------
#	For users with multiple differential gene expression (DE) comparisons:
#	Each comparison is run in a separate folder; for instance (regardless of the number of replicates per condition):
#	Condition_A vs. Condition_B => calculated in 08_diffReps_1
#	Condition_C vs. Condition_D => calculated in 08_diffReps_2
#	Condition_E vs. Condition_F => calculated in 08_diffReps_3
#	etc ...
#	The user needs to copy a 08_diffReps_* folder and use it as a template for any additional DE comparisons 
#	Then the pipeline or individual step can be run as usual
#---------------------------------------------------------------------------------
#bamCorrelate
zMin=0.6
zMax=1.0
#---------------------------------------------------------------------------------
#Generate_Tracks
#Notes: 
#	Please refer to README_Generate_Tracks.sh (located in job-specific folder)
#	Modify Sample_Labels_Color.txt (located in job-specific folder) as needed
#---------------------------------------------------------------------------------
##################################################################################
#DO NOT EDIT CODE BELOW THIS LINE
##################################################################################
#---------------------------------------------------------------------------------
#Need if statements to process variables with <0, 1, or 2>
#For multiple if statements - use elif
#Need to re-assign these variables using the user defined value and convert them to valid arguments for its respective program
#---------------------------------------------------------------------------------
##MEME_ChIP: MEME_MOD
if [ ${MEME_MOD} -eq 0 ]
then
MEME_MOD="oops"
elif [ ${MEME_MOD} -eq 1 ]
then
MEME_MOD="zoops"
elif [ ${MEME_MOD} -eq 2 ]
then
MEME_MOD="anr"
fi
#echo ${MEME_MOD}
#---------------------------------------------------------------------------------
##################################################################################
#---------------------------------------------------------------------------------
#Location of Sample_Labels.txt (text file indicating the samples to process):
#Trying to save the location of this bash script but it will be sourced from another location
#http://stackoverflow.com/questions/59895/can-a-bash-script-tell-which-directory-it-is-stored-in
Sample_Labels_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
#---------------------------------------------------------------------------------
#Need to get the current dir
#This variable is used when creating the job name for qsub commands
SCRIPT_DIR=$(pwd)
#---------------------------------------------------------------------------------
#Webserver location(s) for hosting files
#Note: We are using the waxman-server mount point on the SCC
VM_DIR_FASTQC=/net/waxman-server/mnt/data/waxmanlabvm_home/${BU_User}/FASTQC/${Dataset_Label}
VM_DIR_UCSC=/net/waxman-server/mnt/data/waxmanlabvm_home/${BU_User}/${Dataset_Label}
#---------------------------------------------------------------------------------
#Time hour limit
#On SCC a 12-hour runtime limit is enforced on all jobs, unless specified explicitly. 
#A runtime limit can be specified in the format "hh:mm:ss"
TIME_LIMIT="24:00:00"
#---------------------------------------------------------------------------------
##################################################################################
