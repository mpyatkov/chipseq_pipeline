##################################################################################
#Andy Rampersaud, 11.03.15
#This README contains information about running the BAM_Count job
##################################################################################
#Goal: 		use BAM_Count to count/filter reads in BAM files
#Input:		BAM file (output from Bowtie2)
#Output:	multiple files (see list below) 
#NOTE: This pipeline allows for processing of samples in parallel. See the setup_BAM_Count.sh for details.  
#---------------------------------------------------------------------------------
#The order of script calls:
#1. BAM_Count.sh (sources the setup_BAM_Count.sh)
#2. BAM_Count.qsub is called by BAM_Count.sh
#3. Wait until all jobs have completed running
#4. BAM_Count_Summary.sh then summarizes the jobs for the data set
#---------------------------------------------------------------------------------
#The template folder "BAM_Count" contains the required scripts for running this BAM_Count job
#Location of BAM_Count:
#/restricted/projectnb/waxmanlab/routines/
#1) Copy this folder to your corresponding "Scripts" folder in your data set dir
#2) Update the setup_BAM_Count.sh as needed
#3) Navigate (cd) to the directory where you have all of your files from step (1) and change permissions by running the following on the terminal:
chmod 700 *.sh
chmod 700 *.qsub
#4) Run the BAM_Count.sh:
./BAM_Count.sh
#5) As mentioned above, wait until all jobs have completed running. Then run BAM_Count_Summary.sh:
./BAM_Count_Summary.sh
#This should create a text file summarizing the BAM_Count job
#Note:
#If you want to know each job's duration run the following command in the job folder:
grep 'elapsed' *.o*
#---------------------------------------------------------------------------------
#Job output:
#1. flagstat_G*_M*.txt:
#	Statistics for the BAM file from Bowtie2 mapping (contains mapped and un-mapped reads)
#2. G*_M*_sorted_mapped.bam:
#	BAM file filtered for mapped reads only (for paired-end data these are the "Read-Pairs Aligned Concordantly Exactly 1 Time")
#3. flagstat_G*_M*_sorted_mapped.txt
#	Statistics for the sorted_mapped BAM file 
#4. G*_M*_idxstats.out:
#	The idxstats is used to get the Mapped_Read_Count in (BAM_Count_Summary.sh)
#---------------------------------------------------------------------------------
##################################################################################
#BAM_Count_Summary.sh:
#This script will generate 2 summary tables (showing all samples):
#1. BAM_Count_Stats.txt:
#	Indicates the total read count, mapped reads for downstream analysis, and percentage 
#2. Read_strand_count_Stats.txt
#	Indicates any bias, if any, between positive or negative strand reads 
#---------------------------------------------------------------------------------
#Note about Bowtie2 output (log file of stats) vs. BAM_Count_Stats.txt:
#The filtering command should work for obtaining the Bowtie's reported reads that "aligned concordantly exactly 1 time"
#But the numbers reported in the Bowtie output don't exactly match the numbers from this filtering
#I spent some time trying to de-bug this issue but the read counts are similar enough between the Bowtie2 output and filtering by this script
#It's possible the Bowtie2 calculation for "aligned concordantly exactly 1 time" is different from what we expect
#This BAM_Count job returns slightly fewer reads than what is reported in the Bowtie2 output statistics
#---------------------------------------------------------------------------------
##################################################################################
