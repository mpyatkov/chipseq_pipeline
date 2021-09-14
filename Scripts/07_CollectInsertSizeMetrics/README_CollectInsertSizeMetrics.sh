##################################################################################
#Andy Rampersaud, 11.13.15
#This README contains information about running the CollectInsertSizeMetrics job
##################################################################################
#Goal: 		Use CollectInsertSizeMetrics to collect metrics about the fragment length distribution
#Input:		*_sorted_mapped.bam
#Output:	*_metrics folder
#NOTE: This pipeline allows for processing of samples in parallel. See the setup_CollectInsertSizeMetrics.sh for details.  
#---------------------------------------------------------------------------------
#The order of script calls:
#1. CollectInsertSizeMetrics.sh (sources the setup_CollectInsertSizeMetrics.sh)
#2. CollectInsertSizeMetrics.qsub is called by CollectInsertSizeMetrics.sh
#3. Wait until all jobs have completed running
#4. CollectInsertSizeMetrics_Summary.sh then summarizes the jobs for the data set
#---------------------------------------------------------------------------------
#The template folder "CollectInsertSizeMetrics" contains the required scripts for running this CollectInsertSizeMetrics job
#Location of CollectInsertSizeMetrics:
#/restricted/projectnb/waxmanlab/routines
#---------------------------------------------------------------------------------
#1) Copy this folder to your corresponding "Scripts" folder in your data set dir
#2) Update the setup_CollectInsertSizeMetrics.sh as needed
#3) Navigate (cd) to the directory where you have all of your files from step (1) and change permissions by running the following on the terminal:
chmod 700 *.sh
chmod 700 *.qsub
#4) Run the CollectInsertSizeMetrics.sh:
./CollectInsertSizeMetrics.sh
#5) As mentioned above, wait until all jobs have completed running. Then run CollectInsertSizeMetrics_Summary.sh:
./CollectInsertSizeMetrics_Summary.sh
#This should create a text file summarizing the CollectInsertSizeMetrics job
#Note:
#If you want to know each job's duration run the following command in the job folder:
grep 'elapsed' *.o*
#---------------------------------------------------------------------------------
#Job output:
#A *_metrics folder will be saved in each sample-specific folder
#An overall summary text file (*CollectInsertSizeMetrics_Stats.txt) will be created
#---------------------------------------------------------------------------------
##################################################################################
