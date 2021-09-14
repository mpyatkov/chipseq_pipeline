##################################################################################
#Andy Rampersaud, 11.02.15
#This README contains information about running the Bowtie2 job
##################################################################################
#Goal: 		use Bowtie2 to map paired end reads
#Input:		*.fastq.gz for read1 data and *.fastq.gz for read2 data
#Output:	G*_M*_sorted.bam file and mapping statistics
#NOTE: This pipeline allows for processing of samples in parallel. See the setup_Bowtie2.sh for details.  
#---------------------------------------------------------------------------------
#The order of script calls:
#1. Bowtie2.sh (sources the setup_Bowtie2.sh)
#2. Bowtie2.qsub is called by Bowtie2.sh
#3. Wait until all jobs have completed running
#4. Bowtie2_Summary.sh then summarizes the jobs for the data set
#---------------------------------------------------------------------------------
#The template folder "Bowtie2" contains the required scripts for running this Bowtie2 job
#Location of Bowtie2:
#/restricted/projectnb/waxmanlab/routines
#1) Copy this folder to your corresponding "Scripts" folder in your data set dir
#2) Update the setup_Bowtie2.sh as needed
#3) Navigate (cd) to the directory where you have all of your files from step (1) and change permissions by running the following on the terminal:
chmod 700 *.sh
chmod 700 *.qsub
#4) Run the Bowtie2.sh:
./Bowtie2.sh
#5) As mentioned above, wait until all jobs have completed running. Then run Bowtie2_Summary.sh:
./Bowtie2_Summary.sh
#This should create a text file summarizing the Bowtie2 job
#Note:
#If you want to know each job's duration run the following command in the job folder:
grep 'elapsed' *.o*
#---------------------------------------------------------------------------------
#Job output:
#A "bowtie2" folder will be created within each sample specific folder
#Output file descriptions:
#1. G*_M*_sorted.bam:
#	Sorted alignments directly from Bowtie2 (unfiltered file)
#---------------------------------------------------------------------------------
##################################################################################
#Bowtie2_Summary.sh:
#The output of this script is a Bowtie2_Stats.txt file that summarizes the Bowtie2_Stats for each sample
#Example output:
#Bowtie2 output log file (G1*_M*.e*)
#---------------------------------------------------------------------------------
#250000 reads; of these:
#  250000 (100.00%) were paired; of these:
#    18892 (7.56%) aligned concordantly 0 times
#    175074 (70.03%) aligned concordantly exactly 1 time
#    56034 (22.41%) aligned concordantly >1 times
#    ----
#    18892 pairs aligned concordantly 0 times; of these:
#      2587 (13.69%) aligned discordantly 1 time
#    ----
#    16305 pairs aligned 0 times concordantly or discordantly; of these:
#      32610 mates make up the pairs; of these:
#        24186 (74.17%) aligned 0 times
#        4078 (12.51%) aligned exactly 1 time
#        4346 (13.33%) aligned >1 times
#95.16% overall alignment rate
#---------------------------------------------------------------------------------
#Definitions:
#http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
#Concordant pairs match pair expectations, discordant pairs don't
#A pair that aligns with the expected relative mate orientation and with the expected range of distances between mates is said to align "concordantly". If both mates have unique alignments, but the alignments do not match paired-end expectations (i.e. the mates aren't in the expected relative orientation, or aren't within the expected distance range, or both), the pair is said to align "discordantly". Discordant alignments may be of particular interest, for instance, when seeking structural variants.
#---------------------------------------------------------------------------------
#How to interpret the output:
#There are 3 sections:
#Section 1:
#	Statistics for pairs
#	"aligned concordantly 0 times" = (the mates aren't in the expected relative orientation, or aren't within the expected distance range, or both)
#	18,892 pairs aligned discordantly
#	175,074 pairs aligned concordantly once
#	56,034 pairs aligned concordantly more than once
#Section 2:
#	Statistics for "aligned concordantly 0 times"
#	2,587 pairs aligned discordantly once
#Section 3:
#	Statistics for mates/singlets that align either "0 times concordantly or discordantly"
#	24,186 mates don't align to the genome
#	4,078 mates align uniquely
#	4,346 mates align to multiple locations
#Calculation for overall alignment rate:
#(overall alignment rate) = (number of aligned reads or mates)/(total number of reads or mates)
#(number of aligned reads or mates) = (175074 * 2) + (56034 * 2) + (2587 * 2) + (4078) + (4346)
#(number of aligned reads or mates) = 475814
#(total number of reads or mates) = (250000 * 2) = 500000
#(overall alignment rate) = (475814)/(500000) = 0.9516
#Which is the percentage reported in the Bowtie2 output log file (95.16%)
#---------------------------------------------------------------------------------
##################################################################################
