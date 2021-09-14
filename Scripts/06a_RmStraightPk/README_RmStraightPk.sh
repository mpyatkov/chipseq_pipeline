##################################################################################
#Andy Rampersaud, 11.03.15
#This README contains information about running the RmStraightPk job
##################################################################################
#Goal: 		use RmStraightPk to count/filter reads in BAM files
#Input:		BAM file (output from Bowtie2)
#Output:	StrgtPksRm folder (see list below) 
#NOTE: This pipeline allows for processing of samples in parallel. See the setup_RmStraightPk.sh for details.  
#---------------------------------------------------------------------------------
#The order of script calls:
#1. RmStraightPk.sh (sources the setup_RmStraightPk.sh)
#2. RmStraightPk.qsub is called by RmStraightPk.sh
#3. Wait until all jobs have completed running
#4. RmStraightPk_Summary.sh then summarizes the jobs for the data set
#---------------------------------------------------------------------------------
#The template folder "RmStraightPk" contains the required scripts for running this RmStraightPk job
#Location of RmStraightPk:
#/restricted/projectnb/waxmanlab/routines/
#1) Copy this folder to your corresponding "Scripts" folder in your data set dir
#2) Update the setup_RmStraightPk.sh as needed
#3) Navigate (cd) to the directory where you have all of your files from step (1) and change permissions by running the following on the terminal:
chmod 700 *.sh
chmod 700 *.qsub
#4) Run the RmStraightPk.sh:
./RmStraightPk.sh
#5) As mentioned above, wait until all jobs have completed running. Then run RmStraightPk_Summary.sh:
./RmStraightPk_Summary.sh
#This should create a text file summarizing the RmStraightPk job
#Note:
#If you want to know each job's duration run the following command in the job folder:
grep 'elapsed' *.o*
#---------------------------------------------------------------------------------
#Job output (StrgtPksRm folder):
#1. G*_M*_StrgtPksRm.bam:
#	Output BAM file after straight peak removal
#2. G*_M*_StrgtPksRm.bam.bai:
#	BAM file index (typically used for visualizing the reads)
#3. flagstat_G*_M*_StrgtPksRm.txt
#	Statistics for the BAM file 
#	This file is not very detailed most likely due to a BED -> BAM conversion as part of the RmStraightPk job (information about flags are excluded) 
#---------------------------------------------------------------------------------
##################################################################################
#RmStraightPk_Summary.sh:
#This script will generate 2 summary tables (showing all samples):
#1. StrgtPksRead_count.txt:
#	Indicates the presence of straight read/peaks
#2. StrgtPksRead_confirm_count.txt:
#	Confirms the total read count 
#---------------------------------------------------------------------------------
##################################################################################
