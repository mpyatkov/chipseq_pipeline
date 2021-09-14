##################################################################################
#Andy Rampersaud, 01.23.17
#This README contains information about running the Trim_Galore job
##################################################################################
#Goal: 		use Trim_Galore to detect presence of adaptor sequence in FASTQ files
#Input:		*.fastq.gz for read data
#Output:	*_trimming_report.txt file
#NOTE: This pipeline allows for processing of samples in parallel. See the README_Paired_End_Reads.txt for details.  
#---------------------------------------------------------------------------------
#The order of script calls:
#1. Run_Jobs.sh (sources the 01_Pipeline_Setup.sh)
#2. Trim_Galore.qsub is called by Run_Jobs.sh
#3. Wait until all jobs have completed running
#4. Summarize_Jobs.sh then summarizes the jobs for the data set
#---------------------------------------------------------------------------------
#The template folder "Trim_Galore" contains the required scripts for running this Trim_Galore job
#Location of Trim_Galore:
#/restricted/projectnb/waxmanlab/routines
#1) Copy this folder to your corresponding "Scripts" folder in your data set dir
#2) Update the 01_Pipeline_Setup.sh as needed
#3) Navigate (cd) to the directory where you have all of your files from step (1) and change permissions by running the following on the terminal:
chmod 700 *.sh
chmod 700 *.qsub
#4) Run the Run_Jobs.sh:
./Run_Jobs.sh
#5) As mentioned above, wait until all jobs have completed running. Then run Summarize_Jobs.sh:
./Summarize_Jobs.sh
#This should create a text file summarizing the Trim_Galore job
#Note:
#If you want to know each job's duration run the following command in the job folder:
grep 'elapsed' *.o*
#---------------------------------------------------------------------------------
#Job output:
#A "trim_galore_output" folder will be created within each sample specific folder
#---------------------------------------------------------------------------------
##################################################################################
#Summarize_Jobs.sh:
#The output of this script is a Trim_Galore_Stats.txt file that summarizes the Trim_Galore stats for each sample
##################################################################################
