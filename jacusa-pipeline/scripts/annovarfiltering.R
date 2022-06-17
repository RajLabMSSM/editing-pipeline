#!/usr/bin/env Rscript

library(optparse)
library(dplyr)

option_list <- list(make_option(c('--inAnno'), help = 'Name of the annovar annotated file', default = ''),
                    make_option(c('--inRat'), help = 'Name of the ratio matrix file before SNP filtering', default = ''),
                    make_option(c('--inCov'), help = 'Name of the coverage matrix file before SNP filtering', default = ''),
                    make_option(c('--outRat'), help = 'Name of the ratio matrix final filtered output file', default = ''),
                    make_option(c('--outCov'), help = 'Name of the coverage matrix final filtered output file', default = ''),
                    make_option(c('--outAnno'), help = 'Name of the output file with annotations for editing sites', default = ''))

option.parser <- OptionParser(option_list = option_list)
opt <- parse_args(option.parser)

anno_file <- opt$inAnno
covMatRaw <- opt$inCov
ratMatRaw <- opt$inRat
annotationsOut <- opt$outAnno
covMatOutFilt <- opt$outCov
ratMatOutFilt <- opt$outRat

library(tidyverse)
anno_df <- read_tsv(anno_file)
message( " * read in ", nrow(anno_df), " sites")  
anno_df <- filter(anno_df, dbsnp153CommonSNV == "." & (AF == "." | AF < 0.05) ) 
message( " * keeping ", nrow(anno_df), " sites not overlapping common SNPs")

anno_df <- anno_df %>% mutate(ESid = paste0(Chr, ":", Start, ":", Ref, ":", Alt)) %>% 
    select(ESid, everything() )
message( " * filtering coverage and ratio matrices")

coverage_df <- read_tsv(covMatRaw) %>% column_to_rownames("ESid")
ratio_df <- read_tsv(ratMatRaw) %>% column_to_rownames("ESid")

coverage_df <- coverage_df[ anno_df$ESid,] %>% rownames_to_column("ESid")
ratio_df <- ratio_df[ anno_df$ESid,] %>% rownames_to_column("ESid")

write_tsv(coverage_df, file = covMatOutFilt)
write_tsv(ratio_df, file = ratMatOutFilt)
write_tsv(anno_df, file = annotationsOut)

