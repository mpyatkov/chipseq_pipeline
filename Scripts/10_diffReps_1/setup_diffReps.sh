#!/bin/bash
set -eu
##################################################################################
#Andy Rampersaud, 04.10.17
##################################################################################
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
#Fill in the following information:
##################################################################################
#Information about the conditions being compared
#The differential expression will output the fold change in the following
#format: Treatment_Samples/Control_Samples
#Note:
#"Up" events will be when Treatment>Control
#"Down" events will be when Treatment<Control
#Avoid spaces and parentheses for these variable names: only use underscores for separating words
Control_Samples_NAME="CTCF_Control"
Treatment_Samples_NAME="CTCF_TCPOBOP_1mg"
##################################################################################
#Comparison Number
#If you are doing multiple comparisons within the same dataset, it's helpful to number each comparison
#For example:
#-------------------------------------------------
#Comparison: CAR_ChIP_Vehicles_VS_CAR_ChIP_TCPOBOP
#Output Folders: Output_diffReps_1
#-------------------------------------------------
#Comparison: CEBPa_ChIP_Vehicles_VS_CEBPa_ChIP_TCPOBOP
#Output Folders: Output_diffReps_2
#-------------------------------------------------
COMPAR_NUM=1
##################################################################################
#Specify the step_size (bp) for the diffReps comparison
#As described by the diffReps program:
#---------------------------------------------------------------------------------
#- Genomic region parameters:
#    --mode(peak)     Scanning mode: a selection implies a different window size.
#                     Set window and step size manually to override.
#                     (p)eak      (=1000)  Histone mark peak (Default).
#                     (n)ucleosome(=200)   Single nucleosome (+DNAlinker).
#                     (b)lock     (=10000) Large chromatin modification block.
#    --window(1000)   Window size (default=Histone mark peak size).
#    --step(1/10 win) Window moving step size.
#    --gap(0)         Gap allowed between two consecutive windows.
#---------------------------------------------------------------------------------
#I've been using the following values for delta-DHS analysis:
WINDOW_SIZE=200
#The STEP_SIZE being (1/10) of the WINDOW_SIZE is appropriate
##################################################################################
#Specify the fragment size for the library(s):
#--frag(100)      ChIP-seq library fragment size. Use to shift read positions.
#---------------------------------------------------------------------------------
#For ChIPSeq data:
#Best to use the average fragment size for the samples used in the comparison 
#Example data set:
#Ill forward you the bioanalyzer tracings but here are the numbers I pulled off them when we originally got them:
#G113_M5- 308
#G113_M6- 298
#G113_M7- 313
#G113_M8- 299
#and yeah, there are 60bp adapters ligated to both sides so subtract 120 from the BA results to get the mean insert size
#Sample_ID	#bioanalyzer	#insert size
#G113_M5	308		188
#G113_M6	298		178
#G113_M7	313		193
#G113_M8	299		179
#These numbers can be compared to the output from 06b_CollectInsertSizeMetrics
#---------------------------------------------------------------------------------
#For delta-DHS analysis:
#Regarding diffReps options:
#https://groups.google.com/forum/?hl=en#!searchin/diffreps-discuss/fragment$20size/diffreps-discuss/YdM_uY7klZY/B7bDORP-0jQJ
#p.s.: if you don't want to do any shift, you can set fragment size to 0. 
#For our DHS data, we don't want any read shifting
#---------------------------------------------------------------------------------
#Choose one:
#For ChIPSeq
#FRAG_SIZE=180
#For DHS data:
FRAG_SIZE=0
##################################################################################
#DO NOT EDIT CODE BELOW THIS LINE
##################################################################################
#---------------------------------------------------------------------------------
#Reformat variable names to be compatible for this job
#Replace spaces and parentheses with underscores:
#Control_Samples_NAME
Control_Samples_NAME=$(echo ${Control_Samples_NAME} | sed -e 's/ /_/g' | sed -e 's/(/_/g' | sed -e 's/)/_/g')
#Treatment_Samples_NAME
Treatment_Samples_NAME=$(echo ${Treatment_Samples_NAME} | sed -e 's/ /_/g' | sed -e 's/(/_/g' | sed -e 's/)/_/g')
#---------------------------------------------------------------------------------
##################################################################################
