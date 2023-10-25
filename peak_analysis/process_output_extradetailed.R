## Create xls with overlapping peaks stats

remotes::install_cran("writexl")
remotes::install_cran("argparser")

library(argparser)
library(writexl)

# ## install writexl because it does not present in default module repo
# if (!require(writexl)){
#   ## create lib directory for local user if it is not created
#   if (!file.exists(Sys.getenv("R_LIBS_USER"))) {
#     dir.create(path = Sys.getenv("R_LIBS_USER"), showWarnings = FALSE, recursive = TRUE)
#   }
#   ## install package
#   install.packages("writexl", lib = Sys.getenv("R_LIBS_USER"), repos = "https://cran.rstudio.com/")
# }
# 
# ## install writexl because it does not present in default module repo
# if (!require(argparser)){
#   ## create lib directory for local user if it is not created
#   if (!file.exists(Sys.getenv("R_LIBS_USER"))) {
#     dir.create(path = Sys.getenv("R_LIBS_USER"), showWarnings = FALSE, recursive = TRUE)
#   }
#   ## install package
#   install.packages("argparser", lib = Sys.getenv("R_LIBS_USER"), repos = "https://cran.rstudio.com/")
# }

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
  select(who = Sample_ID, desc = Description) %>% 
  mutate(ix = str_pad(row_number(), width = 2, pad = "0"),
         desc = str_glue("{ix}_{who}\n{desc}")) %>% 
  select(-ix)
  
xls <- map_dfr(list.files(pattern = "*peaks.xls"), function (xls_fname){
  # fname <- str_extract(xls_fname,"^G\\d+\\w+\\d+(?=_MACS)")
  tmp <- read_tsv(xls_fname, comment = "#", col_names = T) %>% 
    select(peak_id = name, length, pileup, fold_enrichment, log10_qvalue = `minus_log10_qvalue`) #-log10(qvalue) -- was originaly in MACS2 file
})

extra_details <- left_join(read_tsv("combined.output", col_names = F), xls, by = c("X8" = "peak_id")) %>% 
  left_join(., sample_labels, by=c("X4" = "who")) %>% 
  mutate(X4 = desc) %>% select(-desc) %>% 
  # group_by(X1, X2, X3, X4) %>% 
  # mutate(max_param = max(fold_enrichment), n = n()) %>%
  # arrange(desc(max_param), .by_group = T) %>% 
  # ungroup() %>% 
  group_by(X1, X2, X3, X4) %>% 
  # mutate(rowid = row_number()) %>% 
  ## because sometimes enrichment can be the same for two peaks
  ## we need to extract only first
  filter(fold_enrichment == max(fold_enrichment) & pileup == max(pileup) & length == max(length)) %>% 
  # select(-max_param, -rowid) %>% 
  ungroup() %>% 
  select(-X5,-X6,-X7, -X8) %>%
  distinct() %>% 
  # group_by(X1, X2, X3) %>% mutate(how_many_samples_intersect = n()) %>% ungroup() %>% 
  pivot_longer(c("length", "pileup", "fold_enrichment", "log10_qvalue"), values_to = "param") %>% 
  # arrange(X4) %>% 
  pivot_wider(names_from = c("X4","name"), values_from = param, names_sort = T)

##### PARTIALY WORKING


bed <- read_tsv("combined.output", col_names = F) %>% 
  select(chr = 1, start = 2, end = 3, who = 4) %>%
  distinct() %>% 
  group_by(chr, start, end) %>% 
  mutate(how_many_samples_intersect = n()) %>%
  #filter(how_many_samples_intersect > 1) %>% 
  ungroup() 

bed <- left_join(bed, sample_labels) %>% 
  #mutate(desc = str_glue("{who}\n{desc}")) %>% 
  select(-who)

bed <- bed %>%  
  mutate(fake_value = 1) %>% 
  pivot_wider(names_from = desc, values_from = fake_value, names_sort = T, values_fill = 0) %>% 
  mutate(ucsc_line = str_glue("{chr}:{start}-{end}")) %>% 
  relocate(ucsc_line,.before = chr) %>% 
  #mutate_all(~replace(., is.na(.), 0)) %>% 
  arrange(desc(how_many_samples_intersect))

# extra_details <- left_join(bed,xls, by = c("X8" = "peak_id")) %>% 
#   group_by(X1, X2, X3,X4) %>% 
#   mutate(max_param = max(fold_enrichment)) %>% filter(fold_enrichment == max_param) %>% select(-max_param) %>% ungroup() %>% 
#   select(-X5,-X6,-X7, -X8) %>%
#   group_by(X1, X2, X3) %>% mutate(how_many_samples_intersect = n()) %>% ungroup() %>% 
#   ungroup %>% 
#   pivot_longer(c("length", "pileup", "fold_enrichment", "log10_qvalue"), values_to = "param") %>% 
#   pivot_wider(names_from = c("X4","name"), values_from = param)

result <- left_join(bed, extra_details, by=c("chr"= "X1","start" = "X2","end" = "X3"))

## write_xlsx(result,'peaks_intersections_extradetailed.xlsx')
write_xlsx(result, str_glue("{argv$output_prefix}_peaks_intersections_extradetailed.xlsx"))



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


