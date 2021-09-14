#This R script is used to calculate peak width statistics
#This script should be placed in the same folder as the MACS Excel file

#Ouptut: Summary table peak width statistics (peak_name_Stats.txt)
#Sample command:
#Rscript Peak_Width_Stats.R G74C_All_MACS_peaks

#Pass arguments from sh into R
args <- commandArgs(trailingOnly = TRUE)

Peak_name <- args[1]

#Need to set the dir so that this script can work in any folder:
dir <- getwd()
wdir <- setwd(dir)

Peak_Data <- read.table(paste(dir, "/", "peak1", ".bed", sep =""), sep = "\t", header = FALSE)
#Since no control was used for MACS there is no FDR calculation
#colnames(Peak_Data) <- c("chr","start","end","length","abs_summit","pileup","MACS_score","fold_enrichment","q-value")
#print("The dim of Peak_Data:")
#print(dim(Peak_Data))

#Calculate peak widths:
#Make an R function to calculate the width of each peak and append to dataframe
peak_width <- function(dataframe){
DHS_width <- dataframe[,3] - dataframe[,2]
dataframe[,dim(dataframe)[2]+1] <- DHS_width
colnames(dataframe)[dim(dataframe)[2]] <- "Peak_Width"
return(dataframe)
}#end of peak_width function
#Call the function:
Peak_Data <- peak_width(Peak_Data)

#print("The dim of Peak_Data (after peak_width function):")
#print(dim(Peak_Data))

#Make Summary Table:
Summary_table <- rbind(summary(Peak_Data$"Peak_Width"))
#> Summary_table
     # Min. 1st Qu. Median Mean 3rd Qu.  Max.
# [1,]  315     970   1167 1614    1561 52890

#Need to parse Summary_table so that output text file is formatted correctly

Summary_table_out <- rbind(Summary_table[1:6])
# > Summary_table_out
     # [,1] [,2] [,3] [,4] [,5]  [,6]
# [1,]  315  970 1167 1614 1561 52890

peakSet <- c(paste(Peak_name," (",dim(Peak_Data)[1],")",sep =""))
# peakSet
# [1] "TCPOBOP_Induced_1_Union (928)"

Summary_table_out2 <- cbind("Peak Set"=peakSet, Summary_table_out)
# Summary_table_out2
     # Peak Set                                                                
# [1,] "TCPOBOP_Induced_1_Union (928)" "315" "970" "1167" "1614" "1561" "52890"

colnames(Summary_table_out2)[2:ncol(Summary_table_out2)] <- colnames(Summary_table)
# > Summary_table_out2
     # Peak Set                        Min.  1st Qu. Median Mean   3rd Qu. Max.   
# [1,] "TCPOBOP_Induced_1_Union (928)" "315" "970"   "1167" "1614" "1561"  "52890"

#Need to set row.names = FALSE to avoid the "[1,]" in the second row

write.table(Summary_table_out2, file = (paste(dir,"/",Peak_name,"_Width_Stats.txt",sep="")), quote=FALSE, sep = "\t",row.names = FALSE,col.names = TRUE)
print("Check out Stats.txt!")
