## Create xls with overlapping peaks stats

library(tidyverse)
library(writexl)

bed <- read_tsv("~/tmp/mnt/max/G196/peak_analysis/combined.output", col_names = F) %>% 
    select(chr = 1, start = 2, end = 3, who = 4) %>%
    distinct() %>% 
    group_by(chr, start, end) %>% 
    summarise(who = paste0(who, collapse = ","), n = n()) %>% 
    select(chr, start, end, samples_intersects = who, how_many_samples_intersects = n) %>% 
    arrange(desc(how_many_samples_intersects))

#write_xlsx(bed,'~/tmp/mnt/max/G196/peak_analysis/peaks_intersections.xlsx')
write_xlsx(bed,'peaks_intersections.xlsx')
