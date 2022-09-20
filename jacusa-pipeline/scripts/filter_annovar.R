#!/usr/bin/env Rscript
# filter ANNOVAR output
library(optparse)
library(dplyr)
library(tidyverse)

option_list <- list(make_option(c('--inAnno'), help = 'Name of the annovar annotated file', default = ''),
                    make_option(c('--inRat'), help = 'Name of the ratio matrix file before SNP filtering', default = ''),
                    make_option(c('--inCov'), help = 'Name of the coverage matrix file before SNP filtering', default = ''),
                    make_option(c('--gencode'), help = 'path to GENCODE GTF', default = ''),
                    make_option(c('--outAnno'), help = 'Name of the output file with annotations for editing sites', default = ''),
                    make_option(c('--outBed'), help = 'Name of the output file with coordinates for editing sites', default = '')
)

option.parser <- OptionParser(option_list = option_list)
opt <- parse_args(option.parser)

anno_file <- opt$inAnno
covMatRaw <- opt$inCov
ratMatRaw <- opt$inRat
gencode <- opt$gencode
annotationsOut <- opt$outAnno
outBed <- opt$outBed

anno_df <- read_tsv(anno_file)
message( " * read in ", nrow(anno_df), " sites")  
anno_df <- filter(anno_df, dbsnp153CommonSNV == "." & AF == "." ) 
message( " * keeping ", nrow(anno_df), " sites not overlapping common and rare SNPs")


# recreate ESid
anno_df <- anno_df %>% mutate(ESid = paste0(Chr, ":", Start, ":", Ref, ":", Alt)) %>% 
    select(ESid, everything() )

# add in gene ID and gene strand
message( " * adding in GENCODE info")

gencode_meta <- read_tsv(gencode)
# remove duplicated gene names
gencode_meta <- gencode_meta[ !duplicated(gencode_meta$gene_name),]

anno_df <- left_join(anno_df, gencode_meta, by = c("Gene.refGene" = "gene_name" ) )

anno_df <- filter(anno_df, !is.na(strand) )

message(" * keeping ", nrow(anno_df), " sites mapped to genes")

message( " * filtering coverage and ratio matrices")

# update ESid to reflect orientations
anno_df <- anno_df %>% mutate(ESid = paste0(Chr, ":", Start, ":", Ref, ":", Alt)) 

# create BED file for pileup
anno_bed <- anno_df %>% select(Chr, Start, End, ESid) %>%
    mutate(Start = Start - 1)

write_tsv(anno_bed, file = outBed, col_names = FALSE)
write_tsv(anno_df, file = annotationsOut)

