####################################################################################
#Andy Rampersaud
#01.08.2016
#This R script is used to generate a histogram of read counts across chromosomes for a single sample
#Input: *_idxstats.out file
#Output: histogram plot
#----------------------------------------------------------------------------------
#Usage: 
#Rscript Chr_Read_Counts.R <idxstats File Name> 
#Example command to run script:
#Rscript Chr_Read_Counts.R G113_M5_idxstats.out
#----------------------------------------------------------------------------------
#Notes:
#Sample from G113_M5_idxstats.out:
#----------------------------------------------------------------------------------
#chr10	129993255	16913	0
#chr11	121843856	18977	0
#chr12	121257530	14791	0
#chr13	120284312	15015	0
#chr14	125194864	14636	0
#----------------------------------------------------------------------------------
#http://www.htslib.org/doc/samtools.html
#idxstats
#    samtools idxstats in.sam|in.bam|in.cram
#    Retrieve and print stats in the index file corresponding to the input file. Before calling idxstats, the input BAM file must be indexed by samtools index.
#    The output is TAB-delimited with each line consisting of reference sequence name, sequence length, # mapped reads and # unmapped reads. It is written to stdout. 
####################################################################################
#---------------------------------------------------------------------------
#Load library:
library(ggplot2)
#---------------------------------------------------------------------------
#Pass arguments from sh into R
args <- commandArgs(trailingOnly = TRUE)
idxstats_File <- args[1]
#---------------------------------------------------------------------------
print("Print arguments:")
print("-----------------")
print("idxstats_File:")
paste(idxstats_File,sep="")
print("-----------------")
#---------------------------------------------------------------------------
#Use a pattern match in an "if" statement 
#to remind user to use the  *.annotated file
#Spaces/brackets matter for "if" statement syntax
#Need invisible() to avoid printing NULL
invisible(if ( grepl("idxstats", idxstats_File) ) {
#Do nothing
} else {
print("WARNING: The input file is not the *.idxstats file!")
}#End of else statement
)#End of invisible
#---------------------------------------------------------------------------
#Need to set the dir so that this script can work in any folder:
dir <- getwd()
setwd(dir)
#---------------------------------------------------------------------------
#Read in idxstats:
idxstats_data <- read.table(idxstats_File, as.is = TRUE, header = FALSE, sep = "\t")
colnames(idxstats_data) <- c("Seq_Name","Seq_Len", "#_Mapped_Reads", "#_Unmapped_Reads")
#View(idxstats_data)
Sample_Name <- gsub(idxstats_File,pattern="_idxstats.out",replacement="", perl=TRUE)
#---------------------------------------------------------------------------
#Read in mm9.chrom.sizes file:
chrom.sizes <- read.table(paste(dir,"/","genomeIndex","/","mm9.chrom.sizes", sep=""), as.is = TRUE, header = FALSE, sep = "\t")
colnames(chrom.sizes) <- c("Chr","Chr_Len")
#View(chrom.sizes)
#---------------------------------------------------------------------------
#https://www.biostars.org/p/17224/
#Using the factor command to redefine the levels of the chr column in the data frame
#Using the "levels" argument in the factor() command
#Once the levels are in a certain sequence, the order() command recapitulates that same sequence
#---------------------------------------------------------------------------
#Need to get chr in the correct order:
chrOrder<-c(paste("chr",1:19,sep=""),"chrX","chrY")
#---------------------------------------------------------------------------
#Filter and sort the data frame:
idxstats_data$"Seq_Name" <- factor(idxstats_data$"Seq_Name", levels=chrOrder)
idxstats_data <- idxstats_data[order(idxstats_data$"Seq_Name"),]
#Filter the chr list:
idxstats_data <- idxstats_data[idxstats_data$"Seq_Name" %in% chrOrder,]
#View(idxstats_data)
#---------------------------------------------------------------------------
#Filter and sort the data frame:
chrom.sizes$"Chr" <- factor(chrom.sizes$"Chr", levels=chrOrder)
chrom.sizes <- chrom.sizes[order(chrom.sizes$"Chr"),]
#Filter the chr list:
chrom.sizes <- chrom.sizes[chrom.sizes$"Chr" %in% chrOrder,]
#View(chrom.sizes)
#---------------------------------------------------------------------------
#Merge data frames:
merge1 <- merge(idxstats_data, chrom.sizes, by.x = "Seq_Name", by.y = "Chr")
#View(merge1)
merge1$"norm_Mapped_Reads" <- (merge1$"#_Mapped_Reads")/(merge1$"Chr_Len")
#View(merge1)
#Express as a percentage:
Total_norm_count <- sum(merge1$"norm_Mapped_Reads")
merge1$"norm_Mapped_Reads_Percent" <- ((merge1$"norm_Mapped_Reads")/(Total_norm_count)) * 100
#View(merge1)
#Check that the sum of percentages = 100
#sum(merge1$"norm_Mapped_Reads_Percent")
#[1] 100
#---------------------------------------------------------------------------
#Format the sum of reads for the plot title:
total_reads <- sum(merge1$"#_Mapped_Reads")
formated_total_reads <- formatC(total_reads, format="d", big.mark=',')
#---------------------------------------------------------------------------
#Now that I have the percentages, I need a bar plot
#Need to rotate x-axis labels:
##http://stackoverflow.com/questions/1330989/rotating-and-spacing-axis-labels-in-ggplot2
#Set y-axis limit
#http://stackoverflow.com/questions/3606697/how-to-set-limits-for-axes-in-ggplot2-r-plots
#---------------------------------------------------------------------------
bar_plot <- ggplot(data=merge1, aes(x=Seq_Name, y=norm_Mapped_Reads_Percent, fill=Seq_Name)) +
#--------------------------------------------------------------------------- 
#Want all the bars to be the same color: #DD8888 (salmon color)
#geom_bar(colour="black", fill="#DD8888", width=.8, stat="identity") +
#---------------------------------------------------------------------------
geom_bar(colour="black", width=.8, stat="identity") +  
guides(fill=FALSE) +
xlab("Chromosome") + ylab("norm_Mapped_Read_Count (%)") +
ggtitle(paste("Mapped Read Count per Chromosome for", "\n", Sample_Name, "\n", "Total Read Count: ", formated_total_reads,sep="")) +
#Minor fix for x-axis labels and tick marks alignment: vjust = .5
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) + 
#ylim(0, 10) 
#http://www.sthda.com/english/wiki/ggplot2-axis-ticks-a-guide-to-customize-tick-marks-and-labels
scale_y_continuous(limits=c(0, 10),breaks=seq(0,10,1))
#---------------------------------------------------------------------------
#Save images to PDF file for best quality
#http://docs.ggplot2.org/current/ggsave.html
#Want to avoid creating empty Rplots.pdf files (need dimensions):
#http://stackoverflow.com/questions/17348359/how-to-stop-r-from-creating-empty-rplots-pdf-file-when-using-ggsave-and-rscript
#Pdf files produce the best quality for single images
#But if you want a montage of images: need to save as png files
#---------------------------------------------------------------------------
ggsave(bar_plot, file=paste(Sample_Name,"_bar_plot.png",sep=""),width = 7, height = 7,units = "in")
#---------------------------------------------------------------------------
print(paste("Check out: ", Sample_Name,"_bar_plot.png",sep=""))
####################################################################################
