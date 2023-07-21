library(dplyr)
library(readr)

args = commandArgs(trailingOnly=TRUE)

stats_path <- args[1]
norm_factors <- read_tsv(stats_path, col_names = T) %>% 
  select(-FRAGMENT_COUNT, -Peak_Union_Count, -FRAGMENT_IN_PEAK_RATIO) %>% 
  mutate(Norm_Factor = round(min(FRAGMENT_IN_PEAK_COUNT)/FRAGMENT_IN_PEAK_COUNT,3)) %>% 
  write_tsv("Norm_Factors.txt", col_names = T)

## for debug
# stats_path <- "/projectnb/wax-dk/max/G205_SG_CHIRP/Scripts/09_RiPPM_Normalization/Job_Summary/Peak_Union_Count_Stats.txt"
## expected output
# SAMPLE_ID       Description     FRAGMENT_IN_PEAK_COUNT  Norm_Factor
# G205M1  LFAR1_ChIRP11_A1_EVEN   703715  0.881
# G205M2  LFAR1_ChIRP11_A1_ODD    619739  1
# G205M3  LFAR1_ChIRP11_B1_EVEN   1248363 0.496
# G205M4  LFAR1_ChIRP11_B1_ODD    1816785 0.341
# G205M5  Lnc14873_ChIRP11_A1_EVEN        1455399 0.425
# G205M6  Lnc14873_ChIRP11_A1_ODD         696552  0.89
# G205M7  Lnc14873_ChIRP11_B1_EVEN        1618677 0.383
# G205M8  Lnc14873_ChIRP11_B1_ODD         2382971 0.261  

