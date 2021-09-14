##################################################################################
#Andy Rampersaud, 11.02.15
#This README contains information about running the FASTQC job
##################################################################################
#Goal: 		QC analysis
#Input:		*.fastq.gz for read1 data and *.fastq.gz for read2 data
#Output:	*_fastqc.html (contains various QC metrics)
#NOTE: This pipeline allows for processing of samples in parallel. See the setup_FASTQC.sh for details.  
#---------------------------------------------------------------------------------
#The order of script calls:
#1. FASTQC.sh (sources the setup_FASTQC.sh)
#2. FASTQC.qsub is called by FASTQC.sh
#3. Wait until all jobs have completed running
#4. FASTQC_Summary.sh then summarizes the jobs for the data set
#---------------------------------------------------------------------------------
#The template folder "FASTQC" contains the required scripts for running this FASTQC job
#Location of FASTQC:
#/restricted/projectnb/waxmanlab/routines/RNASeq_Scripts
#1) Copy this folder to your corresponding "Scripts" folder in your data set dir
#2) Update the setup_FASTQC.sh as needed
#3) Navigate (cd) to the directory where you have all of your files from step (1) and change permissions by running the following on the terminal:
chmod 700 *.sh
chmod 700 *.qsub
#4) Run the FASTQC.sh:
./FASTQC.sh
#5) As mentioned above, wait until all jobs have completed running. Then run FASTQC_Summary.sh:
./FASTQC_Summary.sh
#This should create a text file summarizing the FASTQC job
#Note:
#If you want to know each job's duration run the following command in the job folder:
grep 'elapsed' *.o*
#---------------------------------------------------------------------------------
#The FASTQC_Summary.sh actually copies the *_fastqc.html to the waxmanlabvm where QC analysis can be accessed by all lab members
#---------------------------------------------------------------------------------
##################################################################################
