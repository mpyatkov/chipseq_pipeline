#################################################################################
Andy Rampersaud, 03.21.17
This README contains information about running scripts/jobs for Paired_End_Reads
#################################################################################
The goal of this pipeline is to analyze ChIP/DNase-Seq data starting from the raw FASTQ files from the sequencing facility down to differential region analysis and data visualization
---------------------------------------------------------------------------------
There are multiple steps involved in the pipeline, the following provides an overview of the order/sequence of steps and approximate time per job (sample)
More details for each job can be found within README files within each folder
---------------------------------------------------------------------------------
File organization:
	Files must be organized in the following way:
	You have a data set dir such as:
	wax-es/aramp10/G83_Samples
	Within this dir you have sample specific folders such as:
	G83_M1
	G83_M2
	G83_M3
	G83_M4
	Within each sample specific folder you have a "fastq" folder with:
	File: *.fastq.gz such as:
	G83_M1.fastq.gz
	FASTQ file names are:
	1. User specified as long as they are organized correctly 
	2. Have a ".fastq.gz" file name extension
---------------------------------------------------------------------------------
STEP_ONE:
	Samples to process:
	The user needs to fill out a Sample_Labels.txt file in the 00_Setup_Pipeline directory
	To facilitate processing of samples in parallel we can use a text file that lists the samples to analyze
	Note: this text file is still valid even if there is only one sample to process
	###############################################
	The text file is formatted like the following:
	----------------------------------------------
	Sample_DIR	Sample_ID	Description	Color
	Sample_Waxman-TP17	G83_M1	Male 8wk-pool 1	0,0,255
	Sample_Waxman-TP18	G83_M2	Male 8wk-pool 2	0,0,255
	Sample_Waxman-TP19	G83_M3	Female 8wk-pool 1	255,0,0
	Sample_Waxman-TP20	G83_M4	Female 8wk-pool 2	255,0,0
	----------------------------------------------
	The 1st column: The Sample_DIR name
	The 2nd column: Waxman Lab Sample_ID 
	The 3rd column: Sample's description 
	The 4th column: Sample's wiggle track color (for use in the UCSC Browser)
	###############################################
	See current Sample_Labels.txt as an example
---------------------------------------------------------------------------------
STEP_TWO:
00_Setup_Pipeline
	To run this pipeline: update/fill out the Update 01_Pipeline_Setup.sh
	Variables in the Update 01_Pipeline_Setup.sh will be used in subsequent steps
	Instructions for 03_Run_Pipeline.sh script:
	#Way to run script:
	#Usage: ./03_Run_Pipeline.sh <Start_Step>
	#	<Start_Step> = Needs to be either "Full_Pipeline" or a specific pipeline step 
	#Example: 
	#User wants to run the full pipeline:
	./03_Run_Pipeline.sh Full_Pipeline
	Time: (G113 full dataset) about 6 hours
	#User wants to run a subset of the pipeline starting at a specific step
	#User may want to jump straight to mapping reads:
	./03_Run_Pipeline.sh 04_TopHat_Paired_End
---------------------------------------------------------------------------------
Users have the option to run the full pipeline as described above or run individual steps detailed below:
---------------------------------------------------------------------------------
01_Rename_Folders
	Time: quick step to allow parallel processing of samples
02_FASTQC
	Input: *.fastq.gz files
	Output: *_fastqc.html
	Time per sample: about 10-15 minutes
	*_fastqc.html (contains various QC metrics) can be accessed by all lab members
03_Read_Length
	Input: *.fastq.gz files
	Output: Read_Length_Summary table
	Time per sample: about 10 minutes
	confirm our expected read length 
04_Bowtie2
	Input: *.fastq.gz files
	Output: *_alignments.bam file
	Time per sample: about 1.5 - 2.5 hours
	use Bowtie2 to map paired end reads
05_BAM_Count
	Input: *_alignments.bam file
	Output: *_sorted_mapped.bam
	Time per sample: about 1.5 hours
	Check statistics on the BAM file, extract uniquely mapped proper pairs, check (+/-) strand bias
06_RmStraightPk
	Input: *_sorted_mapped.bam
	Output: *_StrgtPksRm.bam
	Time per sample: about 30 minutes
	Check for bad samples (presence of straight peak reads)
07_CollectInsertSizeMetrics
	Input: *_sorted_mapped.bam
	Output: CollectInsertSizeMetrics_Stats.txt
	Time per sample: about 4 minutes
	Confirm insert size length for paired-end data
08_MACS2
	Input: *_sorted_mapped.bam
	Output: *_MACS2_output folder
	Time per sample: about 15 minutes
	Use macs2 to discover peaks 
09_diffReps_*
	Input: *_sorted_mapped.bam
	Output: Output_diffReps_* folder
	Time per sample: about 30 minutes
	Use diffReps to identify differential sites (delta-TFBS, delta-DHS)
10_MEME_ChIP
	Input: collection of BED files
	Output: Folder containing meme-chip output files
	Time per sample: minutes (about 2 minutes) to hours (about 2.5 hours) 
	The time per sample depends on the number of BED records
	use meme-chip command on BED file(s) regions for motif analysis
11_bamCorrelate
	Input:	BAM file(s) of reads
	Output: Output_bamCorrelate folder (see list below) 
	Time per job: about 20 minutes for full genome scan, about 1 minute for BED file regions
	use the bamCorrelate program to evaluate correlation between samples or replicates
12_normbigWig
	Input: *_sorted_mapped.bam
	Output: G*_M*_norm.bw file
	Time per sample: about 30 - 45 minutes
	Convert from BAM format into wiggle format for UCSC Browser visualization
	Copies required files to the waxmanlabvm for UCSC visualization
13_Generate_Tracks
	Time: quick step 
	Create track lines to provide the UCSC Browser to view data
---------------------------------------------------------------------------------
Pipeline Notes:
---------------------------------------------------------------------------------
I'm running the pipeline but I'm prompted for my password for scp commands?
The scp commands are used for copying files from SCC to waxmanlabvm.
Your password is needed because the user is connecting between servers
If you don't want to type your password for scp commands; try the following:
Use your own <BU_ID>
#--------------------------------------------------
#Following command(s) from SCC:
exec ssh-agent bash
ssh-keygen -t rsa
#Don't overwrite:
#Overwrite (y/n)? n
ssh <BU_ID>@waxmanlabvm.bu.edu mkdir -p .ssh
#Password prompt
cat .ssh/id_rsa.pub | ssh <BU_ID>@waxmanlabvm.bu.edu 'cat >> .ssh/authorized_keys'
#Password prompt
ssh-add ~/.ssh/id_rsa
#Identity added: /usr3/graduate/<BU_ID>/.ssh/id_rsa (/usr3/graduate/<BU_ID>/.ssh/id_rsa)
ssh <BU_ID>@waxmanlabvm.bu.edu
#Works!
#Can now login into waxmanlabvm w/o password prompt
#--------------------------------------------------
I'm prompted for my password when previously I did not have to type my password
I need to restore this feature
I found a shortcut for restoring this feature:
	1. Log into the SCC
	2. Obtain the SCC_ID (a long string): cat .ssh/id_rsa.pub
	3. On the VM: copy and paste this SCC_ID into the .ssh/authorized_keys file
	4. Test the ssh command on the SCC -> VM to confirm password-less login
---------------------------------------------------------------------------------
I'm running the pipeline but I can't keep my terminal/computer open for the full duration of the pipeline?
The pipeline is run with the 03_Run_Pipeline.sh script - ideally the user's terminal can stay open and running until the end of the pipeline's analysis
If you cannot keep your terminal/computer open you can use tmux - terminal multiplexer
On the SCC:
module load tmux
#Start a tmux session:
tmux attach <session_name>
#End a tmux session:
tmux detach
#These tmux sessions act like any other terminal except the terminal can be closed and the pipeline still runs within the tmux terminal
#See the following page for more information:
#http://code.tutsplus.com/tutorials/intro-to-tmux--net-33889
---------------------------------------------------------------------------------
The "Sample_Labels" directory is explained in the job-specific README files
It would be a good idea to test that your text file in your "Sample_Labels" directory is in the correct format
To avoid issues with DOS/MAC characters, run the following:
dos2unix - DOS/MAC to UNIX text file format converter
Example command:
dos2unix Sample_Labels.txt
This should correctly format the Sample_Labels.txt file
Also make sure you only have tab (not spaces) characters between columns
---------------------------------------------------------------------------------
How are counts of reads or fragments handled for this paired-end pipeline?
For any step that collects counts of reads or fragments, I include some description in job-specific summary files (column headers)
For instance, the following jobs have modified column headers:
	*_Read_Length
	*_Bowtie2
	*_BAM_Count
	*_RmStraightPk
	*_MACS2
If there's ambiguity about reads or fragments being counted, refer to the column header for that job summary.
---------------------------------------------------------------------------------
Concluding comments:
I've coded each job/script so that they are relatively easy to run (running a single command will submit multiple jobs) 
If you have questions, requests, or bugs to report, please email aramp10@bu.edu
---------------------------------------------------------------------------------
#################################################################################
