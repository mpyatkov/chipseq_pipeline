##################################################################################
#Andy Rampersaud, 11.11.15
#This README contains information about running the diffReps job
##################################################################################
#Goal: 		use the diffReps program to identify differential sites based on read counts
#Input:		BAM file(s) of reads
#Output:	Output_diffReps folder (see list below)
#NOTE: This pipeline allows for processing of samples in parallel. See below for details.  
#---------------------------------------------------------------------------------
#The order of script calls:
#1. diffReps.sh (sources the setup_diffReps.sh)
#2. diffReps.qsub is called by diffReps.sh
#3. Wait until all jobs have completed running
#4. diffReps_Summary.sh then summarizes the job output
#---------------------------------------------------------------------------------
#The template folder "diffReps" contains the required scripts for running this diffReps job
#Location of diffReps:
#/restricted/projectnb/waxmanlab/routines/
#---------------------------------------------------------------------------------
#1) Copy this folder to your corresponding "Scripts" folder in your data set dir
#2) Update the setup_diffReps.sh as needed
#3) Update the Condition_1.txt as needed
#4) Update the Condition_2.txt as needed
#5) Navigate (cd) to the directory where you have all of your files from step (1) and change permissions by running the following on the terminal:
chmod 700 *.sh
chmod 700 *.qsub
#6) Run the diffReps.sh:
./diffReps.sh
#7) As mentioned above, wait until all jobs have completed running. Then run diffReps_Summary.sh:
./diffReps_Summary.sh
#Note:
#If you want to know each job's duration run the following command in the job folder:
grep 'elapsed' *.o*
#---------------------------------------------------------------------------------
#Replicates per condition
#To handle multiple replicates per condition there are corresponding "Condition_1.txt" and "Condition_2.txt" text files
#The purpose of these text files is to indicate which samples belong to each condition
################################################
#The text file is formatted like the following:
#----------------------------------------------
#Sample_DIR	Sample_ID	Description
#Sample_Waxman-TP17	G83_M1	Male 8wk-pool 1
#Sample_Waxman-TP18	G83_M2	Male 8wk-pool 2
#Sample_Waxman-TP19	G83_M3	Female 8wk-pool 1
#Sample_Waxman-TP20	G83_M4	Female 8wk-pool 2	
#----------------------------------------------
#The 1st column: The Sample_DIR name
#The 2nd column: Waxman Lab Sample_ID 
#The 3rd column: Sample's description 
################################################
#As mentioned in the setup_diffReps.sh:
#Assumptions for this job to run correctly:
#1. You have already run the 04_Bowtie2 job
#2. Your data is organized in the following way:
#You have a data set dir such as:
#/projectnb/wax-es/aramp10/G83_Samples
#Within this dir you have sample specific folders such as:
#G83_M1
#G83_M2
#G83_M3
#G83_M4
#Within each sample specific folder you have a *_R1_*.fastq.gz file and *_R2_*.fastq.gz such as:
#Waxman-TP17_CGATGT_L007_R1_001.fastq.gz
#Waxman-TP17_CGATGT_L007_R2_001.fastq.gz
#Within each sample specific folder you have a "fastq/bowtie2" folder 
#---------------------------------------------------------------------------------
##################################################################################
#Note about n=1 comparisons (only 1 sample per condition)
#If the user tries to give diffReps conditions with only 1 sample it returns the following:
#---------------------------------------------------------------------------------
#To use Negative Binomial, you must have at least two replicates per condition.
#Use G-test(preferred) or Chisquare test instead. Exit.
#---------------------------------------------------------------------------------
#The diffReps.qsub script automatically counts the number of replicates per condition
#If only 1 sample is present in both or either condition, diffReps will automatically be run with the suggested G-test statistics
##################################################################################
#Job output:
#A "Output_diffReps_*" folder will be created within the script folder
#The "Output_diffReps_*" folder will contain:
#https://github.com/shenlab-sinai/diffreps
#(1) diffReps_1:
#	a) Main output file: https://github.com/shenlab-sinai/diffreps/wiki/diffRepsOutput
#(2) diffReps_1.annotated:
#	a) File with annotation of the differential sites based on their locations to the nearest genes
#	b) Differential sites will be assigned to categories: https://github.com/shenlab-sinai/diffreps
#(3) diffReps_1.hotspot:
#	a) Column descriptions: https://github.com/shenlab-sinai/diffreps/wiki/diffRepsOutput
#	b) Algorithm description: https://github.com/shenlab-sinai/diffreps 
#---------------------------------------------------------------------------------
##################################################################################
