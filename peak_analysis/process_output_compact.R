## Create xls with overlapping peaks stats

args = commandArgs(trailingOnly=TRUE)
output_prefix <- args[1]

## install writexl because it does not present in default module repo
if (!require(writexl)){
    ## create lib directory for local user if it is not created
    if (!file.exists(Sys.getenv("R_LIBS_USER"))) {
        dir.create(path = Sys.getenv("R_LIBS_USER"), showWarnings = FALSE, recursive = TRUE)
    }
    ## install package
    install.packages("writexl", lib = Sys.getenv("R_LIBS_USER"), repos = "https://cran.rstudio.com/")
}

require(writexl)
library(tidyverse)

bed <- read_tsv("combined.output", col_names = F) %>% 
    select(chr = 1, start = 2, end = 3, who = 4) %>%
    distinct() %>% 
    group_by(chr, start, end) %>% 
    summarise(who = paste0(who, collapse = ","), n = n()) %>% 
    select(chr, start, end, samples_intersects = who, how_many_samples_intersects = n) %>% 
    arrange(desc(how_many_samples_intersects))

write_xlsx(bed, str_glue("{output_prefix}_peaks_intersections_compact.xlsx"))
