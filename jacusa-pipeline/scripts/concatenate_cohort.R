library(tidyverse)
library(optparse)

option_list <- list(make_option(c('--inDir'), help = 'The path to the Jacusa output directory', default = '.') )

option.parser <- OptionParser(option_list = option_list)
opt <- parse_args(option.parser)

projectDir <- opt$inDir
inDir <- paste0(projectDir, "filter/")

chroms <- paste0('chr', c(1:22,"X","Y") ) 
#BI_MyND_results/filter/chr13_coverage.tsv.gz, BI_MyND_results/filter/chr13_ratio.tsv.gz
rat_files <- paste0(inDir, chroms, "_coverage.tsv.gz")
cov_files <- paste0(inDir, chroms, "_ratio.tsv.gz")

all_cov <- map_df( rat_files, read_tsv) 
all_rat <- map_df( cov_files,read_tsv)

cov_out <- paste0(projectDir, "all_sites_coverage.tsv.gz")
rat_out <- paste0(projectDir, "all_sites_ratio.tsv.gz")
annovar_out <- paste0(projectDir, "avinput.txt")
message( " * ", nrow(all_rat), " sites in total concatenated" )

message(" * writing out coverage")
write_tsv(all_cov, file = cov_out)
message(" * writing out ratios")
write_tsv(all_rat, file = rat_out)

#formatting for annovar
annovar <- tibble(ESid = all_rat$ESid) %>%
    tidyr::separate(col = ESid, into = c("Chr", "Stop", "Ref", "Alt"), sep = ":") %>%
    mutate(Start = Stop) %>%
    select(Chr, Start, Stop, Ref, Alt)

message(" * writing out in VCF format for Annovar")
print(head(annovar) )

write_tsv(annovar, file = annovar_out, col_names = FALSE)

