##################################################################################
#Andy Rampersaud, 11.09.15
#This README contains information about running the normbigWig job
##################################################################################
#Goal: 		use normbigWig for peak discovery in BAM files
#Input:		BAM file (output from Bowtie2)
#Output:	normbigWig folder (see list below) 
#NOTE: This pipeline allows for processing of samples in parallel. See the setup_normbigWig.sh for details.  
#---------------------------------------------------------------------------------
#The order of script calls:
#1. normbigWig.sh (sources the setup_normbigWig.sh)
#2. normbigWig.qsub is called by normbigWig.sh
#3. Wait until all jobs have completed running
#4. normbigWig_Summary.sh then summarizes the jobs for the data set
#---------------------------------------------------------------------------------
#The template folder "normbigWig" contains the required scripts for running this normbigWig job
#Location of normbigWig:
#/restricted/projectnb/waxmanlab/routines/
#1) Copy this folder to your corresponding "Scripts" folder in your data set dir
#2) Update the setup_normbigWig.sh as needed
#3) Navigate (cd) to the directory where you have all of your files from step (1) and change permissions by running the following on the terminal:
chmod 700 *.sh
chmod 700 *.qsub
#4) Run the normbigWig.sh:
./normbigWig.sh
#5) As mentioned above, wait until all jobs have completed running. Then run normbigWig_Summary.sh:
./normbigWig_Summary.sh
#This should create a text file summarizing the normbigWig job
#Note:
#If you want to know each job's duration run the following command in the job folder:
grep 'elapsed' *.o*
#---------------------------------------------------------------------------------
#Job output (normbigWig folder):

#---------------------------------------------------------------------------------
##################################################################################
#normbigWig_Summary.sh:
#This script will generate # summary tables (showing all samples):

#---------------------------------------------------------------------------------
##################################################################################
