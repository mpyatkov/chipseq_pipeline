## Create xls with overlapping peaks stats

## install writexl because it does not present in default module repo
if (!require(writexl)){
  ## create lib directory for local user if it is not created
  if (!file.exists(Sys.getenv("R_LIBS_USER"))) {
    dir.create(path = Sys.getenv("R_LIBS_USER"), showWarnings = FALSE, recursive = TRUE)
  }
  ## install package
  install.packages("writexl", lib = Sys.getenv("R_LIBS_USER"), repos = "https://cran.rstudio.com/")
}

## install writexl because it does not present in default module repo
if (!require(argparser)){
  ## create lib directory for local user if it is not created
  if (!file.exists(Sys.getenv("R_LIBS_USER"))) {
    dir.create(path = Sys.getenv("R_LIBS_USER"), showWarnings = FALSE, recursive = TRUE)
  }
  ## install package
  install.packages("argparser", lib = Sys.getenv("R_LIBS_USER"), repos = "https://cran.rstudio.com/")
}

require(writexl)
require(argparser)
library(tidyverse)

ParseArguments <- function() {
  p <- arg_parser('Peaks processing')
  p <- add_argument(p,'--combined_bed', default="combined.output", help="intersections with all samples")
  p <- add_argument(p,'--sample_labels', default="../Scripts/00_Setup_Pipeline/Sample_Labels.txt", help="to extract addition information we can use tab separated Sample_Labels.txt")
  p <- add_argument(p,'--output_prefix', default="output", help="output prefix to distinguish narrow/broad peaks")
  return(parse_args(p))
}

argv <- ParseArguments()

sample_labels <- read_tsv(argv$sample_labels, col_names=T) %>% 
  select(who = Sample_ID, desc = Description)

bed <- read_tsv("combined.output", col_names = F) %>% 
  select(chr = 1, start = 2, end = 3, who = 4) %>%
  distinct() %>% 
  group_by(chr, start, end) %>% 
  mutate(how_many_samples_intersect = n()) %>%
  #filter(how_many_samples_intersect > 1) %>% 
  ungroup() 
  
bed <- left_join(bed, sample_labels) %>% 
  mutate(who = str_glue("{who}\n{desc}")) %>% 
  select(-desc)

bed <- bed %>%  
  mutate(fake_value = 1) %>% 
  pivot_wider(names_from = who, values_from = fake_value) %>% 
  mutate(ucsc_line = str_glue("{chr}:{start}-{end}")) %>% 
  relocate(ucsc_line,.before = chr) %>% 
  mutate_all(~replace(., is.na(.), 0)) %>% 
  arrange(desc(how_many_samples_intersect))

## write_xlsx(bed,'peaks_intersections_detailed.xlsx')
write_xlsx(bed, str_glue("{argv$output_prefix}_peaks_intersections_detailed.xlsx"))

## create xlsx workbook
# library(openxlsx)
# wb <- createWorkbook()
# addWorksheet(wb, sheetName = "summary")
# ## colNames false for header because we do not need this information inside xlsx
# writeData(wb, sheet = "summary", binded_samples_header, startRow = 1, startCol = 1, colNames = FALSE)
# writeData(wb, sheet = "summary", binded_samples, startRow = 3, startCol = 1)
# 
# ## apply bold style for rpkm/tpm columns
# rpkm_tpm_cols <- which(str_detect(colnames(binded_samples), "rpkm|tpm") == TRUE)
# boldStyle <- createStyle(textDecoration = "Bold") ## create bold style
# addStyle(wb, "summary", boldStyle, 3, rpkm_tpm_cols)
# 
# ## apply number style for column 'total'
# comma_style_format <- createStyle(numFmt = "#,##0") # create thousands format
# addStyle(wb, "summary", comma_style_format, rows = 1, cols = 2:ncol(binded_samples))
# 
# saveWorkbook(wb, output_xls, overwrite = T)


