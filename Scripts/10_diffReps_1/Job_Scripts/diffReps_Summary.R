####################################################################################
#Andy Rampersaud
#10.26.2017
#This R script is used to parse the diffReps output files
#In this context I use differential to mean "condition-specific"
#In other words, differential sites = condition-specific sites
#This output is parsed to determine the number of differential peaks
#It also generates BED files of differential peaks
#This script would be placed in the same directory as the output files
#Output file description:
#https://code.google.com/p/diffreps/wiki/diffRepsOutput
#----------------------------------------------------------------------------------
#Notes:
#<Output File Name>	refers to the name of the file with *.annotated extension (26 columns)
#<Sample1 Name>		Control group name (from the diffReps run)
#<Sample2 Name> 	Treatment group name (from the diffReps run)
#<Peak Caller/Differential program used>	can be MACS2/MAnorm/diffReps
#<output folder name>	will be the name of the output folder for the parsed diffReps files
#<log2FC_cutoff>	Cutoff for the |log2FC| for delta-sites  
#			If you want 2-fold delta-sites, use <log2FC_cutoff> = 1
#			Note the FC_cutoff_label = "2-fold"
#			If you want 1.5-fold delta-sites, use <log2FC_cutoff> = 0.5849
#			Note the FC_cutoff_label = "1.5-fold"
#----------------------------------------------------------------------------------
#Usage: 
#Rscript diffReps_Summary.R <Output File Name> <Sample1 Name> <Sample2 Name> <Peak Caller/Differential program used> <output folder name> <log2FC_cutoff>
#----------------------------------------------------------------------------------
#Example command to run script:
#Rscript diffReps_Summary.R diffReps_1.annotated STAT5_Low STAT5_High diffReps_1 Peaks_Filtered_Summary 1
#----------------------------------------------------------------------------------
#Example command to run script:
#Actually call the shell script: diffReps_Summary.sh
#cd /media/Internal_Data_01/Documents/Waxman_Lab/Notes/Scripts/diffReps_Summary/Scripts/
##Initialize parameters:
#Control_Samples_NAME=Male_3hr_Control_K4me3
#Treatment_Samples_NAME=Male_3hr_TCPOBOP_K4me3
#COMPAR_NUM=1
#./diffReps_Summary.sh 'diffReps_'${Control_Samples_NAME}'.vs.'${Treatment_Samples_NAME}.annotated ${Control_Samples_NAME} ${Treatment_Samples_NAME} diffReps_${COMPAR_NUM}
#----------------------------------------------------------------------------------
#Example command to run script:
#Actually call the shell script: diffReps_Summary.sh
#cd /media/Internal_Data_01/Documents/Waxman_Lab/Notes/Scripts/diffReps_Summary/Scripts/
##Initialize parameters:
#Control_Samples_NAME=Male_3hr_Control_K4me3
#Treatment_Samples_NAME=Male_3hr_TCPOBOP_K4me3_exclude_M11
#COMPAR_NUM=2
#./diffReps_Summary.sh 'diffReps_'${Control_Samples_NAME}'.vs.'${Treatment_Samples_NAME}.annotated ${Control_Samples_NAME} ${Treatment_Samples_NAME} diffReps_${COMPAR_NUM}
####################################################################################
#---------------------------------------------------------------------------
#Load library:
library(ggplot2)
#---------------------------------------------------------------------------
#Pass arguments from sh into R
args <- commandArgs(trailingOnly = TRUE)
Output_File <- args[1]
S1_diff <- args[2]
S2_diff <- args[3]
Peak_Caller <- args[4]
output_dir <- args[5]
log2FC_cutoff <- as.numeric(args[6])
#Need a label to indicate the FC cutoff being used:
FC_cutoff_label <- round(2^log2FC_cutoff, digits=2)
min_avg_count <- 20
#---------------------------------------------------------------------------
print("Print arguments:")
print("-----------------")
print("Output_File:")
paste(Output_File,sep="")
print("S1_diff:")
paste(S1_diff,sep="")
print("S2_diff:")
paste(S2_diff,sep="")
print("Peak_Caller:")
paste(Peak_Caller,sep="")
print("output_dir:")
paste(output_dir,sep="")
print("min_avg_count:")
paste(min_avg_count,sep="")
print("-----------------")
#---------------------------------------------------------------------------
#Use a pattern match in an "if" statement 
#to remind user to use the  *.annotated file
#Spaces/brackets matter for "if" statement syntax
#Need invisible() to avoid printing NULL
invisible(if ( grepl("annotated", Output_File) ) {
#Do nothing
} else {
print("WARNING: Output_File is not the *.annotated file!")
}#End of else statement
)#End of invisible
#---------------------------------------------------------------------------
#Need to set the dir so that this script can work in any folder:
dir <- getwd()
setwd(dir)
#---------------------------------------------------------------------------
#Check to make sure *.annotated file exists
#Need invisible() to avoid printing NULL
invisible(if ( length(list.files(pattern = "annotated")) >= 1)  {
#Do nothing
} else {
print("diffReps job most likely failed.")
print("Quitting R now.")
quit()
}#End of else statement
)#End of invisible
#---------------------------------------------------------------------------
diffReps_output <- read.table(Output_File, as.is = TRUE, header = TRUE, sep = "\t", fill = TRUE)
#colnames(diffReps_output)
#Gives the correct 26 columns:
#[1] "Chrom"         "Start"         "End"           "Length"       
# [5] "Treatment.cnt" "Control.cnt"   "Treatment.avg" "Control.avg"  
# [9] "Treatment.enr" "Control.enr"   "Event"         "log2FC"       
#[13] "pval"          "padj"          "winSta"        "winEnd"       
#[17] "winFC"         "winP"          "winQ"          "GName"        
#[21] "TName"         "Strand"        "TSS"           "TES"          
#[25] "Feature"       "D2TSS"     
#dim(diffReps_output)
#head(diffReps_output)
#factor(diffReps_output$Event)[0]
#Always be:
#Levels: Down Up
############################################################
#Also read in *.hotspot data
#Need to get Output_File name without the ".annotation"
#Output_File_Prefix <- gsub(Output_File,pattern=".annotated",replacement="")
#Added functionality to only interogate differential sites overlapping called peaks:
#Output_File_Prefix <- gsub(Output_File,pattern=".annotated.peaks",replacement="")
#Need perl=TRUE 
#".*" stands for any number of any character
#".annotated.*": any number of any character after ".annotated"
#http://stackoverflow.com/questions/11776287/remove-pattern-from-string-with-gsub
Output_File_Prefix <- gsub(Output_File,pattern=".annotated.*",replacement="", perl=TRUE)
diffReps_hotspot_output <- read.table(paste(Output_File_Prefix,".hotspot",sep=""), as.is = TRUE, header = TRUE, sep = "\t", fill = TRUE)
#colnames(diffReps_hotspot_output)
#Gives the correct 11 columns:
#[1] "Chrom"    "Start"    "End"      "Length"   "diff"   "pval"    
# [7] "padj"     "nsite"    "Sites"    "ntype"    "MarkType"
############################################################
#Put Summary files in a Summary folder
output_dir <- paste(dir,"/",output_dir,sep="")
#Create dir:
system(paste("if [ ! -d ",output_dir," ]; then mkdir ",output_dir,"; else rm -r ",output_dir,"/*","; fi",sep=""))
#Set working dir to be output_dir
setwd(output_dir)
#Organize output files within Summary folder
#-------------------------------------------
#Need UCSC_Files folder
UCSC_Files_dir <- paste(output_dir,"/","UCSC_Files", sep="")
#Create dir:
system("if [ ! -d UCSC_Files ]; then mkdir UCSC_Files; fi")
#-------------------------------------------
#Need BED_Files folder
BED_Files_dir <- paste(output_dir,"/","BED_Files", sep="")
#Create dir:
system("if [ ! -d BED_Files ]; then mkdir BED_Files; fi")
#-------------------------------------------
#Need Annotated_Files folder
Annotated_Files_dir <- paste(output_dir,"/","Annotated_Files", sep="")
#Create dir:
system("if [ ! -d Annotated_Files ]; then mkdir Annotated_Files; fi")
#-------------------------------------------
#Need PCA_BED_Files folder
PCA_BED_Files_dir <- paste(output_dir,"/","PCA_BED_Files", sep="")
#Create dir:
system("if [ ! -d PCA_BED_Files ]; then mkdir PCA_BED_Files; fi")
#-------------------------------------------
############################################################
#Hotspot regions are a mix of induced and repressed sites
#Just make a UCSC BED file of these regions
#Sort data_frames
diffReps_hotspot_output <- diffReps_hotspot_output[order(diffReps_hotspot_output$"Chrom",diffReps_hotspot_output$"Start"),]
diffReps_hotspot_BED <- cbind(diffReps_hotspot_output$"Chrom",diffReps_hotspot_output$"Start",diffReps_hotspot_output$"End","diffReps_hotspot","1000",".","0","0","0,0,0")
#Make BED files:
write.table(diffReps_hotspot_BED, file = paste(Peak_Caller,"_","Hotspots_UCSC",".bed", sep =""), quote=FALSE, sep = "\t",row.names = FALSE,col.names = FALSE)
#Remember: Still need a track line to view BED file in UCSC Browser
#Use system to execute Linux commands within R
system(paste("echo track name=",Peak_Caller,"_","Hotspots visibility=4 itemRgb=On > Header.txt",sep =""))
system(paste("cat Header.txt ",Peak_Caller,"_Hotspots_UCSC.bed > temp1.bed",sep =""))
system(paste("mv temp1.bed ",Peak_Caller,"_Hotspots_UCSC.bed",sep =""))
system('rm Header.txt')
#Move *_UCSC.bed files to UCSC_Files_dir 
system(paste("mv *_UCSC.bed ",UCSC_Files_dir,sep=""))
############################################################
print("The number of condition-specific sites between samples:")
#S1_diff_data will be sites which are "down"
#Other words: up in Control data, down in treatment data
#Also want |Fold Change| > 2
#Note that delta DHS with log2FC of -Inf will be included
#--------------------------------------
#With spotty data, it's possible to get many marginally differential sites
#Better to implement a minimum count
#Also want a min Avg. count cutoff (20)
#--------------------------------------
S1_diff_data <- diffReps_output[diffReps_output$Event =="Down" & diffReps_output$log2FC < -(log2FC_cutoff) & diffReps_output$Control.avg > min_avg_count & diffReps_output$padj < 0.05,]
paste(S1_diff,":",sep="")
dim(S1_diff_data)
#S2_diff_data will be sites which are "up"
#Other words: down in Control data, up in treatment data
#Also want |Fold Change| > 2
#Note that delta DHS with log2FC of Inf will be included
#--------------------------------------
#With spotty data, it's possible to get many marginally differential sites
#Better to implement a minimum count
#Also want a min Avg. count cutoff (20)
#--------------------------------------
#Also need to use a FDR cutoff (0.05)
#This incorporates information about sample variability
#Whether the differnce observed is significant
#--------------------------------------
S2_diff_data <- diffReps_output[diffReps_output$Event =="Up" & diffReps_output$log2FC > (log2FC_cutoff) & diffReps_output$Treatment.avg > min_avg_count & diffReps_output$padj < 0.05,]
paste(S2_diff,":",sep="")
dim(S2_diff_data)
print("-----------------")
#---------------------------------------------------------------------------
#Need to check that differential peaks exist (non-zero count)
#If there are zero delta sites, this R script needs to exit
if(dim(S1_diff_data)[1] == 0 | dim(S2_diff_data)[1] == 0){
	print("There were zero delta sites found in 1 or both conditions.")
	print("Quitting R now.")
	quit()
}
#End of if statement
#---------------------------------------------------------------------------
#Need to create PCA BED files (top 200 and top 600 diffReps sites based on log2FC)
#Sort diffReps_output by log2FC (descending)
#--------------------------------------
#Also need to use a FDR cutoff (0.05)
#This incorporates information about sample variability
#Whether the differnce observed is significant
#--------------------------------------
diffReps_output_sorted <- diffReps_output[order(-diffReps_output$"log2FC") & diffReps_output$padj < 0.05,]
#----------------------------------------
#Get the top 200:
Top_200 <- diffReps_output_sorted[1:200,]
#Sort this data_frame by chr and start:
Top_200_sorted <- Top_200[order(Top_200$"Chrom",Top_200$"Start"),]
#Make BED file of chr, start, end:
Top_200_BED <- cbind(Top_200_sorted$"Chrom",Top_200_sorted$"Start",Top_200_sorted$"End")
#----------------------------------------
#Get the top 600:
Top_600 <- diffReps_output_sorted[1:600,]
#Sort this data_frame by chr and start:
Top_600_sorted <- Top_600[order(Top_600$"Chrom",Top_600$"Start"),]
#Make BED file of chr, start, end:
Top_600_BED <- cbind(Top_600_sorted$"Chrom",Top_600_sorted$"Start",Top_600_sorted$"End")
#----------------------------------------
write.table(Top_200_BED, file = paste(S2_diff,"_","diffReps_Top_200",".bed", sep =""), quote=FALSE, sep = "\t",row.names = FALSE,col.names = FALSE)
write.table(Top_600_BED, file = paste(S2_diff,"_","diffReps_Top_600",".bed", sep =""), quote=FALSE, sep = "\t",row.names = FALSE,col.names = FALSE)
#Move *.bed files to PCA_BED_Files_dir 
system(paste("mv *_Top_*.bed ",PCA_BED_Files_dir,sep=""))
#---------------------------------------------------------------------------
#Additional average count filter to visualize strong differential sites:
#Arbitrary cutoffs of 100 and 200 returned handful of sites
#--------------------------------------
#Also need to use a FDR cutoff (0.05)
#This incorporates information about sample variability
#Whether the differnce observed is significant
#--------------------------------------
S1_diff_screenshot <- diffReps_output[diffReps_output$Event =="Down" & diffReps_output$log2FC < -(log2FC_cutoff) & diffReps_output$Control.avg > 100 & diffReps_output$padj < 0.05,]
#Sort by log2FC (most different at the top)
#When sorting the down sites, sorting negative FC values, want to sort ascending
S1_diff_screenshot <- S1_diff_screenshot[order(S1_diff_screenshot$"log2FC"),]
#--------------------------------------
S2_diff_screenshot <- diffReps_output[diffReps_output$Event =="Up" & diffReps_output$log2FC > (log2FC_cutoff) & diffReps_output$Treatment.avg > 200 & diffReps_output$padj < 0.05,]
#Sort by log2FC (most different at the top)
S2_diff_screenshot <- S2_diff_screenshot[order(-S2_diff_screenshot$"log2FC"),]
#--------------------------------------
write.table(S1_diff_screenshot, file = paste(S1_diff,"_",Peak_Caller,"_screenshot.txt", sep =""), quote=FALSE, sep = "\t",row.names = FALSE,col.names = TRUE)
write.table(S2_diff_screenshot, file = paste(S2_diff,"_",Peak_Caller,"_screenshot.txt", sep =""), quote=FALSE, sep = "\t",row.names = FALSE,col.names = TRUE)
#Move *_screenshot.txt files to UCSC_Files_dir 
system(paste("mv *_screenshot.txt ",UCSC_Files_dir,sep=""))
########################################################################################
#I want to see the effect of FDR cutoff on the number of delta sites discovered
#Implement filters on the original output file:
#dim(diffReps_output)
#[1] 8065   26
#For this 1st filter I need an "OR" statement (not "AND")
Filter_01 <- diffReps_output[diffReps_output$Control.avg > min_avg_count | diffReps_output$Treatment.avg > min_avg_count,]
#dim(Filter_01)
#[1] 6100   26
#R boolean for "or" is "|"
Filter_02 <- Filter_01[Filter_01$log2FC < -(log2FC_cutoff) | Filter_01$log2FC > (log2FC_cutoff),]
#dim(Filter_02)
#[1] 3082   26
#Now want to subset Filter_02 by different FDR values
#-----------------------------------------------------
#Even though there are a handful of FDR cutoffs, it would be best to use an R function:
#---------------------------------------------------------------------------
#R function to parse data frame for a single FDR cutoff
FDR_filter <- function(FDR_cutoff){
FDR_cutoff <- as.numeric(FDR_cutoff)
#print(FDR_cutoff)
Filtered_Data <- Filter_02[Filter_02$padj <= FDR_cutoff ,]
Filtered_Table <-  as.data.frame(table(Filtered_Data$"Event"))
#Data frames with 0 data rows will return an error
#Need to check that data exists:
#http://stackoverflow.com/questions/28556658/why-does-an-empty-dataframe-fail-an-is-null-test
#Need an if statement:
if (is.data.frame(Filtered_Table) && nrow(Filtered_Table)==0) {
#If the table is empty (i.e. zero delta sites)
print(paste("WARNING: zero delta sites for the FDR cutoff: ", FDR_cutoff, sep=""))
#Then manually create the table
Var1 <- c("Down", "Up")
Freq <- c(0, 0)
#Avoid converting data type: use cbind.data.frame (not cbind)
#http://stackoverflow.com/questions/19535996/avoid-rbind-cbind-conversion-from-numeric-to-factor
Filtered_Table <- cbind.data.frame(Var1, Freq)
}
#End of if statement
Filtered_Table$"FDR_Label" <- paste("FDR<=", FDR_cutoff, sep="")
#Function return:
output_list <- list(Filtered_Data, Filtered_Table)
return(output_list)
}#end of FDR_filter function
#---------------------------------------------------------------------------
#Creating list objects:
#Access the Filtered_Data: 
#FDR_lt_0.05[[1]]
#Access the Filtered_Table:
#FDR_lt_0.05[[2]]
#---------------------------------------------------------------------------
#Call the function:
FDR_lt_0.05 <- FDR_filter(0.05)
FDR_lt_0.01 <- FDR_filter(0.01)
FDR_lt_0.005 <- FDR_filter(0.005)
FDR_lt_0.001 <- FDR_filter(0.001)
#---------------------------------------------------------------------------
#Now I can concatenate these tables:
FDR_Counts <- rbind(FDR_lt_0.05[[2]], FDR_lt_0.01[[2]], FDR_lt_0.005[[2]], FDR_lt_0.001[[2]])
#Make a grouped bar plot:
#http://stackoverflow.com/questions/17721126/simplest-way-to-do-grouped-barplot
FDR_plot <- ggplot(FDR_Counts, aes(factor(FDR_Label), Freq, fill=Var1)) +
geom_bar(stat="identity", position = "dodge") + 
scale_fill_brewer(palette = "Set1") +
ggtitle(paste("Effect of FDR cutoff on number condition-specific sites")) + 
ylab("Count of Condition-specific Regions") + 
xlab("FDR Threshold") +
guides(fill=guide_legend(title="Event")) +
geom_text(aes(label = Freq, ymax=0), size = 3, position=position_dodge(width=0.9), vjust=-0.25)
#http://stackoverflow.com/questions/13576139/how-to-display-value-in-a-stacked-bar-chart-by-using-geom-text
#Save the plot:
ggsave(FDR_plot, file=paste(Peak_Caller,"_FDR_Barchart.pdf",sep=""),width = 7, height = 7,units = "in")
########################################################################################
#Make a histogram of diffReps differential sites
#Use ggplot2 for histogram:
#http://stackoverflow.com/questions/7027448/change-histogram-bar-colours-greater-than-a-certain-value
#qplot(log2FC,data=diffReps_output,geom="histogram",fill=color)
#The qplot function is supposed make the same graphs as ggplot, but with a simpler syntax. However, in practice, it's often easier to just use ggplot because the options for qplot can be more confusing to use.
#http://www.cookbook-r.com/Graphs/Plotting_distributions_%28ggplot2%29/
#ggplot(diffReps_output, aes(x=log2FC)) + geom_histogram(fill=diffReps_output$color)
#Add "event" (DHS_opening, DHS_closing, Less_2fold) to diffReps_output data_frame
#-----------------------------------------------------
#To be more general, DHS_event should be Site_Category
#"DHS_closing" should be S1_diff_sites
#"DHS_opening" should be S2_diff_sites
#-----------------------------------------------------
#Need a modified histogram to show marginally differential sites
#The diffReps_output data_frame needs to be split into the following:
#1. S2_real_delta
#2. S2_marginal_delta
#3. Less_2fold
#4. S1_marginal_delta
#5. S1_real_delta
S2_real_delta <- diffReps_output[diffReps_output$Event =="Down" & diffReps_output$log2FC < -(log2FC_cutoff) & diffReps_output$Control.avg > min_avg_count,]
#---------------------------------------------------------------------------
#Regarding the colors:
#Some of the Site_Category groups may not exist (i.e. zero rows)
#I should add a color column to the diffReps_output_hist data frame
#Make use of the if statement below:
#---------------------------------------------------------------------------
S2_real_delta_label <- paste("1.",S1_diff,"_","sites", sep ="")
#Check to make sure rows exist in S2_real_delta
#Need invisible() to avoid printing NULL
invisible(if ( dim(S2_real_delta)[1] == 0)  {
print("Warning: There are zero rows in S2_real_delta!")
} else {
#add the Site_Category labels
S2_real_delta$Site_Category <- S2_real_delta_label
}#End of else statement
)#End of invisible
#---------------------------------------------------------------------------
S2_marginal_delta <- diffReps_output[diffReps_output$Event =="Down" & diffReps_output$log2FC < -(log2FC_cutoff) & diffReps_output$Control.avg <= min_avg_count,]
#---------------------------------------------------------------------------
S2_marginal_delta_label <- paste("2.","Marginal","_","sites", sep ="")
#Check to make sure rows exist in S2_marginal_delta
#Need invisible() to avoid printing NULL
invisible(if ( dim(S2_marginal_delta)[1] == 0)  {
print("Warning: There are zero rows in S2_marginal_delta!")
} else {
#add the Site_Category labels
S2_marginal_delta$Site_Category <- S2_marginal_delta_label
}#End of else statement
)#End of invisible
#---------------------------------------------------------------------------
Less_2fold <- diffReps_output[diffReps_output$log2FC >= -(log2FC_cutoff) & diffReps_output$log2FC <= (log2FC_cutoff),]
#---------------------------------------------------------------------------
Less_2fold_label <- paste0("3.Less_",FC_cutoff_label,"-fold")
#Check to make sure rows exist in Less_2fold
#Need invisible() to avoid printing NULL
invisible(if ( dim(Less_2fold)[1] == 0)  {
print("Warning: There are zero rows in Less_FC_cutoff!")
} else {
#add the Site_Category labels
Less_2fold$Site_Category <- Less_2fold_label
}#End of else statement
)#End of invisible
#---------------------------------------------------------------------------
S1_marginal_delta <- diffReps_output[diffReps_output$Event =="Up" & diffReps_output$log2FC > (log2FC_cutoff) & diffReps_output$Treatment.avg <= min_avg_count,]
#---------------------------------------------------------------------------
S1_marginal_delta_label <- paste("4.","Marginal","_","sites", sep ="")
#Check to make sure rows exist in S1_marginal_delta
#Need invisible() to avoid printing NULL
invisible(if ( dim(S1_marginal_delta)[1] == 0)  {
print("Warning: There are zero rows in S1_marginal_delta!")
} else {
#add the Site_Category labels
S1_marginal_delta$Site_Category <- S1_marginal_delta_label
}#End of else statement
)#End of invisible
#---------------------------------------------------------------------------
S1_real_delta <- diffReps_output[diffReps_output$Event =="Up" & diffReps_output$log2FC > (log2FC_cutoff) & diffReps_output$Treatment.avg > min_avg_count,]
#---------------------------------------------------------------------------
S1_real_delta_label <- paste("5.",S2_diff,"_","sites", sep ="")
#Check to make sure rows exist in S1_real_delta
#Need invisible() to avoid printing NULL
invisible(if ( dim(S1_real_delta)[1] == 0)  {
print("Warning: There are zero rows in S1_real_delta!")
} else {
#add the Site_Category labels
S1_real_delta$Site_Category <- S1_real_delta_label
}#End of else statement
)#End of invisible
#---------------------------------------------------------------------------
print("-----------------")
#Confirmed that sum of the rows of the above equals the original number of rows:
#dim(S2_real_delta)[1] + dim(S2_marginal_delta)[1] + dim(Less_2fold)[1] + dim(S1_marginal_delta)[1] + dim(S1_real_delta)[1]
#equals: dim(diffReps_output)[1]
#-----------------------------------------------------
#Now concatenate data_frames for histogram creation
#Need to combine by rows (use rbind)
diffReps_output_hist <- rbind(S2_real_delta, S2_marginal_delta, Less_2fold, S1_marginal_delta, S1_real_delta)
#dim(diffReps_output_hist)
#equals: dim(diffReps_output)[1]
#-----------------------------------------------------
#diffReps_output$Site_Category[diffReps_output$log2FC < -(log2FC_cutoff) ] <- paste("1.",S1_diff,"_","sites", sep ="")
#diffReps_output$Site_Category[diffReps_output$log2FC > (log2FC_cutoff) ] <- paste("3.",S2_diff,"_","sites", sep ="")
#diffReps_output$Site_Category[diffReps_output$log2FC >= -(log2FC_cutoff) & diffReps_output$log2FC <= (log2FC_cutoff)] <- "2.Less_2fold"
#-----------------------------------------------------
#Use different colors from ggplot default
#http://www.cookbook-r.com/Graphs/Colors_%28ggplot2%29/
# If the title is long, it can be split into multiple lines with \n
#http://www.cookbook-r.com/Graphs/Titles_%28ggplot2%29/
#Set y-axis label
#http://www.cookbook-r.com/Graphs/Axes_%28ggplot2%29/
#Various legend options:
#http://docs.ggplot2.org/0.9.2.1/theme.html
#control the legend font size using:
#+ theme(legend.text=element_text(size=X))
#Default size looks like size=10
#---------------------------------------------------------------------------
#Manually setting group colors for ggplot2
#http://stackoverflow.com/questions/17180115/manually-setting-group-colors-for-ggplot2
#Need a "named vector" to do this:
#Initialize list of colors:
hist.colors <- data.frame("blue","gray","black","gray","red")
#Assign names using the above labels:
names(hist.colors) <- c(S2_real_delta_label, S2_marginal_delta_label, Less_2fold_label, S1_marginal_delta_label, S1_real_delta_label)
#transpose data frame:
hist.colors <- t(hist.colors)
hist.colors <- as.data.frame(hist.colors)
#Create the named vector:
#http://stackoverflow.com/questions/19265172/converting-two-columns-of-a-data-frame-to-a-named-vector
#Note the order of the vector (first colors, then group name)
hist.colors.named.vector <- setNames(as.character(hist.colors$V1), as.character(rownames(hist.colors)))
#Plot command:
hist_plot <- ggplot(diffReps_output_hist, aes(x=log2FC,  fill=Site_Category)) + geom_histogram(binwidth=.1) + scale_fill_manual(values=hist.colors.named.vector) + ggtitle(paste("Fold Change for diffReps condition-specific sites","\n"," (",dim(diffReps_output_hist)[1],")"," Total Sites",sep="")) + ylab("Count of Condition-specific Regions") + xlab("log2(Fold Change)") + theme(legend.background = element_rect(colour = "black"), legend.text=element_text(size=10))
#---------------------------------------------------------------------------
#hist_plot <- ggplot(diffReps_output_hist, aes(x=log2FC,  fill=Site_Category)) + geom_histogram(binwidth=.1) + scale_fill_manual(values=c("blue","gray","black","gray","red")) + ggtitle(paste("Fold Change for diffReps condition-specific sites","\n"," (",dim(diffReps_output_hist)[1],")"," Total Sites",sep="")) + ylab("Count of Condition-specific Regions") + xlab("log2(Fold Change)") + theme(legend.background = element_rect(colour = "black"), legend.text=element_text(size=10))
#Save images to PDF file for best quality
#http://docs.ggplot2.org/current/ggsave.html
#Want to avoid creating empty Rplots.pdf files (need dimensions):
#http://stackoverflow.com/questions/17348359/how-to-stop-r-from-creating-empty-rplots-pdf-file-when-using-ggsave-and-rscript
ggsave(hist_plot, file=paste(Peak_Caller,"_Histogram.pdf",sep=""),width = 7, height = 7,units = "in")
########################################################################################
#############################################################
##Single pie chart:
#############################################################
##Want to create a pie chart summarizing the "Feature" column
##Need to initialize a device to save the plot to a png file:
#png(filename=paste(Peak_Caller,"_Feature_Distribution.png",sep=""))
##Useful to know if the majority of differential sites are proximal to gene/gene regions
##Convert column to table
#mytable <- table(diffReps_output$Feature)
##Get percentages:
#pct <- round(prop.table(mytable) * 100)
##lbls will be the labels for the pie chart
#lbls <- paste(names(mytable), "\n", mytable, sep="")
#lbls <- paste(lbls," (",pct,sep="") # add percents to labels 
#lbls <- paste(lbls,"%",")",sep="") # add % sign to labels 
#pie(mytable, labels = lbls, main = paste("Feature Distribution of Differential Sites"," (",dim(diffReps_output)[1],")",sep=""))
##Save plot to the png file
#dev.off()
############################################################
############################################################
#Want to make separate pie charts for the induced and repressed sites on the same image
#Use par to get 2 charts in 2 rows and 1 column
#par(mfrow=c(2,1))
#par does not work well with pie charts (labels and legends overlap)
#Just make separate png files
########################################################################################
##Make a make_piechart function
#make_piechart <- function(data_frame,Peak_name){
##Need to initialize a device to save the plot to a png file:
#png(filename=paste(Peak_name,"_Feature_Dist.png",sep=""))
##Convert column to table
#mytable <- table(data_frame$Feature)
##Get percentages:
#pct <- round(prop.table(mytable) * 100)
##lbls will be the labels for the pie chart
#lbls <- paste(names(mytable), "\n", mytable, sep="")
#lbls <- paste(lbls," (",pct,sep="") # add percents to labels 
#lbls <- paste(lbls,"%",")",sep="") # add % sign to labels 
##Assign colors for legend:
#colors <- c("dodgerblue","firebrick","forestgreen","darkorchid4","dodgerblue4","darkorange","firebrick","gold","limegreen")
#suppressWarnings(pie(mytable, labels = lbls, col=colors, cex =0.8, fill=colors,main = paste("Feature Distribution of ",Peak_name," Sites"," (",dim(data_frame)[1],")",sep="")))
##Legend on the left:
#legend(-1.2, 1, names(mytable), cex=0.8, fill=colors)
##Save plot to the png file
#dev.off()
#}#end of make_piechart function
########################################################################################
##Call the function:
#make_piechart(S1_diff_data, paste(S1_diff,sep =""))
#make_piechart(S2_diff_data, paste(S2_diff,sep =""))
#######################################################################################
#Make a make_stacked_bar_plot function
make_stacked_bar_plot <- function(data_frame,Peak_name){
#Need to initialize a device to save the plot to a png file:
png(filename=paste(Peak_name,"_Feature_Dist.png",sep=""))
#Convert column to table
mytable <- table(data_frame$Feature)
#Get percentages:
pct <- round(prop.table(mytable) * 100)
###################################################
#Sort the pct to get all the high percentages first
#pct <- sort(pct, decreasing=TRUE)
#Sorting the contingency table will make it hard to compare different barplots
#Still get the issue of the legend covering the last bar
###################################################
#--------------------------------------------------
#Want to add color to the data_frame based on the genomic feature
#Only 8 genomic features, separate command for each feature
#https://code.google.com/p/diffreps/
#-----------------------------------
#ProximalPromoter 	+/- 250bp of TSS
#Promoter1k 	+/- 1kbp of TSS
#Promoter3k 	+/- 3kbp of TSS
#Genebody 	Anywhere between a gene's promoter and up to 1kbp downstream of the TES.
#Genedeserts 	Genomic regions that are depleted with genes and are at least 1Mbp long.
#Pericentromere 	Between the boundary of a centromere and the closest gene minus 10kbp of that gene's regulatory region.
#Subtelomere 	Similary defined as pericentromere.
#OtherIntergenic 	Any region that does not belong to the above categories
#-----------------------------------
#Add a color to pct data_frame
#Need to convert to data_frame
pct <- data.frame(pct)
#Add column names
colnames(pct) <- c("Feature","Percentage")
pct$color[pct$Feature == "ProximalPromoter"] <- "dodgerblue"
pct$color[pct$Feature == "Promoter1k"] <- "firebrick"
pct$color[pct$Feature == "Promoter3k"] <- "forestgreen"
pct$color[pct$Feature == "Genebody"] <- "darkorchid4"
pct$color[pct$Feature == "Genedesert"] <- "dodgerblue4"
pct$color[pct$Feature == "Pericentromere"] <- "darkorange"
pct$color[pct$Feature == "Subtelomere"] <- "midnightblue"
pct$color[pct$Feature == "OtherIntergenic"] <- "gold"
#colors <- c("dodgerblue","firebrick","forestgreen","darkorchid4","dodgerblue4","darkorange","midnightblue","gold","limegreen")
#--------------------------------------------------
#Want the legend outside of plot:
# Add extra space to right of plot area; change clipping to figure
#mar parameters change margin space
#The vector is ordered, the first value corresponding to the bottom. The entire array is c(bottom, left, top, right)
par(mar=c(5.1, 4.1, 4.1, 10.1), xpd=TRUE)
#Option xaxt="n" will remove x-axis labels
suppressWarnings(barplot(pct$Percentage, col=pct$color,pch=1,xaxt="n",ylab="Percentage of Sites (%)",main = paste("Feature Distribution of ","\n",Peak_name," Sites"," (",dim(data_frame)[1],")",sep="")))
#inset parameters move legend
legend("topright", levels(pct$Feature), inset=c(-0.5,0), fill=pct$color);
#x <- suppressWarnings(barplot(pct, col=colors,xaxt="n",ylab="Percentage of Sites (%)"))
#Since I have a legend the text labels are redundant:
##Put x-axis labels back on
#labs <- paste(names(pct))
##Playing around with the options, these work:
#text(cex=0.9, x=x-.70, y=-4.5, labs, xpd=TRUE, srt=45)
#Save plot to the png file
dev.off()
}#end of make_stacked_bar_plot function
#######################################################################################
##Call the function:
make_stacked_bar_plot(S1_diff_data, paste(S1_diff,sep =""))
make_stacked_bar_plot(S2_diff_data, paste(S2_diff,sep =""))
#---------------------------------------------------------------------------
#Combine the *.png files:
#-append stacks vertically
#+append stack horizontally
system(paste("convert +append *_Feature_Dist.png +append ",Peak_Caller,"_Feature_Distribution.png",sep=""))
system("rm *_Feature_Dist.png")
#---------------------------------------------------------------------------
#Sort data_frames
S1_diff_data <- S1_diff_data[order(S1_diff_data$"Chrom",S1_diff_data$"Start"),]
S2_diff_data <- S2_diff_data[order(S2_diff_data$"Chrom",S2_diff_data$"Start"),]
#Make text files of delta peaks (with gene annotation information)
write.table(S1_diff_data, file = paste(S1_diff,"_",Peak_Caller,"_annotated.txt", sep =""), quote=FALSE, sep = "\t",row.names = FALSE,col.names = TRUE)
write.table(S2_diff_data, file = paste(S2_diff,"_",Peak_Caller,"_annotated.txt", sep =""), quote=FALSE, sep = "\t",row.names = FALSE,col.names = TRUE)
#Move *_annotated.txt files to Annotated_Files_dir 
system(paste("mv *_annotated.txt ",Annotated_Files_dir,sep=""))
#---------------------------------------------------------------------------
#Add color to BED files to get a single BED file of all peaks with 3 colors
#S1_diff color = Blue (0,0,255)
#S2_diff color = Red 	(255,0,0)
#BG color = Green (0,128,0)
#Need BED files in BED9 format:
#http://genome.ucsc.edu/FAQ/FAQformat.html#format1
#chrom
#chromStart
#chromEnd
#name 
#score
#strand
#thickStart
#thickEnd
#itemRgb
#---------------------------------------------------------------------------
S1_diff_data_BED <- cbind(S1_diff_data$"Chrom",S1_diff_data$"Start",S1_diff_data$"End",S1_diff,"1000",".","0","0","0,0,255")
S2_diff_data_BED <- cbind(S2_diff_data$"Chrom",S2_diff_data$"Start",S2_diff_data$"End",S2_diff,"1000",".","0","0","255,0,0")
#---------------------------------------------------------------------------
#Make BED files of delta regions, run Peak_Width.sh, summarize output 
write.table(S1_diff_data_BED, file = paste(S1_diff,"_",Peak_Caller,".bed", sep =""), quote=FALSE, sep = "\t",row.names = FALSE,col.names = FALSE)
write.table(S2_diff_data_BED, file = paste(S2_diff,"_",Peak_Caller,".bed", sep =""), quote=FALSE, sep = "\t",row.names = FALSE,col.names = FALSE)
#Move *.bed files to BED_Files_dir 
system(paste("mv *.bed ",BED_Files_dir,sep=""))
#---------------------------------------------------------------------------
#Calculate peak widths:
#Make an R function to calculate the peak width and print to output file
peak_width <- function(data_frame,Peak_name){
DHS_width <- data_frame[,3] - data_frame[,2]
data_frame[,dim(data_frame)[2]+1] <- DHS_width
colnames(data_frame)[dim(data_frame)[2]] <- "Peak_Width"
Summary_table <- rbind(summary(data_frame$"Peak_Width"))
#Use round() to format the numbers (no need for decimals)
Summary_table_out <- rbind(round(Summary_table[1:6]))
peakSet <- c(paste(Peak_name," (",dim(data_frame)[1],")",sep =""))
Summary_table_out2 <- cbind("Peak Width Distribution (bp)"=peakSet, Summary_table_out)
colnames(Summary_table_out2)[2:ncol(Summary_table_out2)] <- colnames(Summary_table)
write.table(Summary_table_out2, file = (paste(Peak_name,"_Width_Stats.txt",sep="")), quote=FALSE, sep = "\t",row.names = FALSE,col.names = TRUE)
#print("Check out Stats.txt!")
return(data_frame)
}#end of peak_width function
#Call the function:
S1_diff_Call <- peak_width(S1_diff_data, paste(S1_diff,"_",Peak_Caller,sep =""))
S2_diff_Call <- peak_width(S2_diff_data, paste(S2_diff,"_",Peak_Caller,sep =""))
#---------------------------------------------------------------------------
#Combine data_frames, sort, output to single BED file
#Need list of data_frames
Name_List <- c('S1_diff_data_BED','S2_diff_data_BED')
BED_List <- lapply(Name_List, get)
Full_BED <- do.call(rbind, BED_List)
colnames(Full_BED) <- c("chr","start","end","name","score","strand","thickStart","thickEnd","itemRgb")
#Sort Full_BED
#Need to convert to data_frame to be able to sort
Full_BED=as.data.frame(Full_BED)
Full_BED <- Full_BED[order(Full_BED$"chr",Full_BED$"start"),]
write.table(Full_BED, file = paste(Peak_Caller,"_Peaks_UCSC.bed", sep =""), quote=FALSE, sep = "\t",row.names = FALSE,col.names = FALSE)
#Remember: Still need a track line to view BED file in UCSC Browser
#Use system to execute Linux commands within R
system(paste("echo track name=",Peak_Caller,"_Peaks visibility=4 itemRgb=On > Header.txt",sep =""))
system(paste("cat Header.txt ",Peak_Caller,"_Peaks_UCSC.bed > temp1.bed",sep =""))
system(paste("mv temp1.bed ",Peak_Caller,"_Peaks_UCSC.bed",sep =""))
system('rm Header.txt')
#Move *_UCSC.bed files to UCSC_Files_dir 
system(paste("mv *_UCSC.bed ",UCSC_Files_dir,sep=""))
#---------------------------------------------------------------------------       
S1_diff_peaks <- dim(S1_diff_data)[1]
S2_diff_peaks <- dim(S2_diff_data)[1]
#---------------------------------------------------------------------------
output_table <- matrix(c(S1_diff_peaks, S2_diff_peaks), ncol=1,byrow=TRUE)
colnames(output_table) <- c(paste("Number of Condition-specific Sites","\n","\t","(at least ",FC_cutoff_label,"-fold different)","\n","\t","(min_avg_count > ",min_avg_count,")","\n","\t","(FDR < 0.05)" ,sep=''))
rownames(output_table) <- c(paste(S1_diff,"_",Peak_Caller,sep=""),paste(S2_diff,"_",Peak_Caller,sep=""))
output_table <- as.table(output_table)
write.table(output_table, file = paste(Peak_Caller,"_peak_count.txt",sep=""), quote=FALSE, sep = "\t",col.names=NA)
#Append all the *_Width_Stats.txt files to the diffReps_peak_count.txt
system(paste("cat ",Peak_Caller,"_peak_count.txt ",S1_diff,"_",Peak_Caller,"_Width_Stats.txt ",S2_diff,"_",Peak_Caller,"_Width_Stats.txt > temp1.txt",sep =""))
system(paste("mv temp1.txt ",Peak_Caller,"_peak_count.txt",sep=""))
system("rm *_Width_Stats.txt")
print("Check out the Summary folder!")
#####################################################################################
