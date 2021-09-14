####################################################################################
#Andy Rampersaud
#01.11.2016
#This R script is used to generate a histogram of peak counts across chromosomes for a single sample
#Input: *.bed file
#Output: histogram plot
#----------------------------------------------------------------------------------
#Usage: 
#Rscript Chr_Peak_Counts.R <BED File Name> 
#Example command to run script:
#Rscript Chr_Peak_Counts.R G113_M5_MACS2_peaks.narrowPeak.bed
#----------------------------------------------------------------------------------
#Notes:
#Sample from G113_M5_MACS2_peaks.narrowPeak.bed:
#----------------------------------------------------------------------------------
#chr1	197067792	197068019	G113_M5_MACS2_output/G113_M5_MACS2_peak_1
#chr10	21862288	21862571	G113_M5_MACS2_output/G113_M5_MACS2_peak_2
#chr11	3032016	3032208	G113_M5_MACS2_output/G113_M5_MACS2_peak_3
#chr11	53953599	53953797	G113_M5_MACS2_output/G113_M5_MACS2_peak_4
#chr11	108872999	108873407	G113_M5_MACS2_output/G113_M5_MACS2_peak_5
#----------------------------------------------------------------------------------
#These are MACS2 BED files using the UCSC standard BED format:
#https://genome.ucsc.edu/FAQ/FAQformat.html#format1
#The first three required BED fields are:
#    chrom - The name of the chromosome (e.g. chr3, chrY, chr2_random) or scaffold (e.g. scaffold10671).
#    chromStart - The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.
#    chromEnd - The ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display of the feature. For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99. 
#The 9 additional optional BED fields are:
#    name - Defines the name of the BED line. This label is displayed to the left of the BED line in the Genome Browser window when the track is open to full display mode or directly to the left of the item in pack mode. 
#----------------------------------------------------------------------------------
####################################################################################
#---------------------------------------------------------------------------
#Load library:
library(ggplot2)
#---------------------------------------------------------------------------
#Pass arguments from sh into R
args <- commandArgs(trailingOnly = TRUE)
BED_File <- args[1]
#---------------------------------------------------------------------------
print("Print arguments:")
print("-----------------")
print("BED_File:")
paste(BED_File,sep="")
print("-----------------")
#---------------------------------------------------------------------------
#Use a pattern match in an "if" statement 
#to remind user to use the  *.annotated file
#Spaces/brackets matter for "if" statement syntax
#Need invisible() to avoid printing NULL
invisible(if ( grepl("bed", BED_File) ) {
#Do nothing
} else {
print("WARNING: The input file is not the *.bed file!")
}#End of else statement
)#End of invisible
#---------------------------------------------------------------------------
#Need to set the dir so that this script can work in any folder:
dir <- getwd()
setwd(dir)
#---------------------------------------------------------------------------
#Read in BED file:
BED_data <- read.table(BED_File, as.is = TRUE, header = FALSE, sep = "\t")
colnames(BED_data) <- c("chrom","chromStart", "chromEnd", "name")
#View(BED_data)
Sample_Name <- gsub(BED_File,pattern=".bed",replacement="", perl=TRUE)
#---------------------------------------------------------------------------
#Read in mm9.chrom.sizes file:
chrom.sizes <- read.table(paste(dir,"/","genomeIndex","/","mm9.chrom.sizes", sep=""), as.is = TRUE, header = FALSE, sep = "\t")
colnames(chrom.sizes) <- c("Chr","Chr_Len")
#View(chrom.sizes)
#---------------------------------------------------------------------------
#Need to get a peak count per chromosome
#Use the table() command to a count table
peak_count_table <- as.data.frame(table(BED_data$"chrom"))
colnames(peak_count_table) <- c("chrom", "peak_count")
#View(peak_count_table)
#Confirm that the sum(peak_count_table$"peak_count") equals the total number of peaks
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
peak_count_table$"chrom" <- factor(peak_count_table$"chrom", levels=chrOrder)
peak_count_table <- peak_count_table[order(peak_count_table$"chrom"),]
#Filter the chr list:
peak_count_table <- peak_count_table[peak_count_table$"chrom" %in% chrOrder,]
#View(peak_count_table)
#---------------------------------------------------------------------------
#Filter and sort the data frame:
chrom.sizes$"Chr" <- factor(chrom.sizes$"Chr", levels=chrOrder)
chrom.sizes <- chrom.sizes[order(chrom.sizes$"Chr"),]
#Filter the chr list:
chrom.sizes <- chrom.sizes[chrom.sizes$"Chr" %in% chrOrder,]
#View(chrom.sizes)
#---------------------------------------------------------------------------
#Merge data frames:
merge1 <- merge(peak_count_table, chrom.sizes, by.x = "chrom", by.y = "Chr")
#View(merge1)
merge1$"norm_Peaks" <- (merge1$"peak_count")/(merge1$"Chr_Len")
#View(merge1)
#Express as a percentage:
Total_norm_count <- sum(merge1$"norm_Peaks")
merge1$"norm_Peaks_Percent" <- ((merge1$"norm_Peaks")/(Total_norm_count)) * 100
#View(merge1)
#Check that the sum of percentages = 100
#sum(merge1$"norm_Peaks_Percent")
#[1] 100
#---------------------------------------------------------------------------
#Need an adaptive ylim_plot:
#For samples where few peaks are called sometimes the chr percentage is high
#Want ylim = 100
#For samples where many peaks are called, each chr is contributing about 5% (expected)
#Want ylim = 10
#Need a check on the max_percent:
max_percent <- max(merge1$"norm_Peaks_Percent")
if (max_percent >= 10) {
#print("max_percent is greater than or equal to 10")
ylim_plot <- 100
breaks_value <- 10
} else {
#print("max_percent is NOT greater than or equal to 10")
ylim_plot <- 10
breaks_value <- 1
}
#End of if/else statement
#---------------------------------------------------------------------------
#Format the Total Peak Count for the plot title:
total_peaks <- dim(BED_data)[1]
formated_total_peaks <- formatC(total_peaks, format="d", big.mark=',')
#---------------------------------------------------------------------------
#Now that I have the percentages, I need a bar plot
#Need to rotate x-axis labels:
##http://stackoverflow.com/questions/1330989/rotating-and-spacing-axis-labels-in-ggplot2
#Set y-axis limit
#http://stackoverflow.com/questions/3606697/how-to-set-limits-for-axes-in-ggplot2-r-plots
#---------------------------------------------------------------------------
bar_plot <- ggplot(data=merge1, aes(x=chrom, y=norm_Peaks_Percent, fill=chrom)) +
#--------------------------------------------------------------------------- 
#Want all the bars to be the same color: #DD8888 (salmon color)
#geom_bar(colour="black", fill="#DD8888", width=.8, stat="identity") +
#---------------------------------------------------------------------------
geom_bar(colour="black", width=.8, stat="identity") +
#Setting luminance and saturation (chromaticity)
#http://www.cookbook-r.com/Graphs/Colors_%28ggplot2%29/
scale_fill_hue(l=40) +
guides(fill=FALSE) +
xlab("Chromosome") + ylab("norm_Peak_Count (%)") +
ggtitle(paste("Peak Count per Chromosome for", "\n", Sample_Name, "\n", "Total Peak Count: ", formated_total_peaks, sep="")) +
#Minor fix for x-axis labels and tick marks alignment: vjust = .5
theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) + 
#---------------------------------------------------------------------------
#ylim(0, 10) 
#http://www.sthda.com/english/wiki/ggplot2-axis-ticks-a-guide-to-customize-tick-marks-and-labels
#Adaptive y limits:
scale_y_continuous(limits=c(0, ylim_plot),breaks=seq(0,ylim_plot,breaks_value))
#For typical sample:
#scale_y_continuous(limits=c(0, 10),breaks=seq(0,10,1))
#Want to see full range:
#scale_y_continuous(limits=c(0, 100),breaks=seq(0,100,10))
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
