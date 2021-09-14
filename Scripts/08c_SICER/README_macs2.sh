##################################################################################
#Andy Rampersaud, 11.09.15
#This README contains information about running the macs2 job
##################################################################################
#Goal: 		use macs2 for peak discovery in BAM files
#Input:		BAM file (output from Bowtie2)
#Output:	macs2 folder (see list below) 
#NOTE: This pipeline allows for processing of samples in parallel. See the setup_macs2.sh for details.  
#---------------------------------------------------------------------------------
#The order of script calls:
#1. macs2.sh (sources the setup_macs2.sh)
#2. macs2.qsub is called by macs2.sh
#3. Wait until all jobs have completed running
#4. macs2_Summary.sh then summarizes the jobs for the data set
#---------------------------------------------------------------------------------
#The template folder "macs2" contains the required scripts for running this macs2 job
#Location of macs2:
#/restricted/projectnb/waxmanlab/routines/
#1) Copy this folder to your corresponding "Scripts" folder in your data set dir
#2) Update the setup_macs2.sh as needed
#3) Navigate (cd) to the directory where you have all of your files from step (1) and change permissions by running the following on the terminal:
chmod 700 *.sh
chmod 700 *.qsub
#4) Run the macs2.sh:
./macs2.sh
#5) As mentioned above, wait until all jobs have completed running. Then run macs2_Summary.sh:
./macs2_Summary.sh
#This should create a text file summarizing the macs2 job
#Note:
#If you want to know each job's duration run the following command in the job folder:
grep 'elapsed' *.o*
#---------------------------------------------------------------------------------
#Job output (macs2 folder):

#---------------------------------------------------------------------------------
##################################################################################
#macs2_Summary.sh:
#This script will generate # summary tables (showing all samples):

#---------------------------------------------------------------------------------
##################################################################################
