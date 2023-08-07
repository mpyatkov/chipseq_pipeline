#!/bin/bash
set -eu
##################################################################################
#Andy Rampersaud, 04.10.17
##################################################################################

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

Control_Samples_NAME="Female_Control"
Treatment_Samples_NAME="Male_Treatment"

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
WINDOW_SIZE=200
#The STEP_SIZE being (1/10) of the WINDOW_SIZE is appropriate

##################################################################################
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

#For delta-DHS analysis:
#Regarding diffReps options:
#https://groups.google.com/forum/?hl=en#!searchin/diffreps-discuss/fragment$20size/diffreps-discuss/YdM_uY7klZY/B7bDORP-0jQJ
#p.s.: if you don't want to do any shift, you can set fragment size to 0. 
#For our DHS data, we don't want any read shifting

#Choose one:
#For ChIPSeq
#FRAG_SIZE=180
#For DHS data:
FRAG_SIZE=0

##################################################################################
#DO NOT EDIT CODE BELOW THIS LINE
##################################################################################
#Reformat variable names to be compatible for this job
#Replace spaces and parentheses with underscores:
#Control_Samples_NAME
Control_Samples_NAME=$(echo ${Control_Samples_NAME} | sed -e 's/ /_/g' | sed -e 's/(/_/g' | sed -e 's/)/_/g')
#Treatment_Samples_NAME
Treatment_Samples_NAME=$(echo ${Treatment_Samples_NAME} | sed -e 's/ /_/g' | sed -e 's/(/_/g' | sed -e 's/)/_/g')

