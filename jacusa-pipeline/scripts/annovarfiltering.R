#!/usr/bin/env Rscript

library(optparse)
library(dplyr)

option_list <- list(make_option(c('--inAnno'), help = 'Name of the annovar annotated file', default = ''),
                    make_option(c('--inRat'), help = 'Name of the ratio matrix file before SNP filtering', default = ''),
                    make_option(c('--inCov'), help = 'Name of the coverage matrix file before SNP filtering', default = ''),
                    make_option(c('--outRat'), help = 'Name of the ratio matrix final filtered output file', default = ''),
                    make_option(c('--outCov'), help = 'Name of the coverage matrix final filtered output file', default = ''),
                    make_option(c('--gencode'), help = 'path to GENCODE GTF', default = ''),
                    make_option(c('--outAnno'), help = 'Name of the output file with annotations for editing sites', default = ''),
                    make_option(c('--outBed'), help = 'Name of the output file with coordinates for editing sites', default = '')
)

option.parser <- OptionParser(option_list = option_list)
opt <- parse_args(option.parser)

anno_file <- opt$inAnno
covMatRaw <- opt$inCov
ratMatRaw <- opt$inRat
annotationsOut <- opt$outAnno
covMatOutFilt <- opt$outCov
ratMatOutFilt <- opt$outRat
gencode <- opt$gencode
outBed <- opt$outBed
library(tidyverse)

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

anno_df <- left_join(anno_df, gencode_meta, by = c("Gene.refGene" = "gene_name" ) )

anno_df <- filter(anno_df, !is.na(strand) )

message(" * keeping ", nrow(anno_df), " sites mapped to genes")


# flip allele orientation if gene is on negative strand
revcomp <- function(x){
  data.frame(y = c("A","C","G","T"),row.names =  c("T","G","C","A"))[x,]
}
message(" * flipping orientations" )

anno_df <- mutate(anno_df, 
                  Ref = ifelse(!is.na(strand) & strand == "-", revcomp(Ref), Ref ),
                  Alt = ifelse(!is.na(strand) & strand == "-", revcomp(Alt), Alt )
                  )
message( " * filtering coverage and ratio matrices")

coverage_df <- read_tsv(covMatRaw) %>% column_to_rownames("ESid")
ratio_df <- read_tsv(ratMatRaw) %>% column_to_rownames("ESid")

coverage_df <- coverage_df[ anno_df$ESid,] %>% rownames_to_column("ESid")
ratio_df <- ratio_df[ anno_df$ESid,] %>% rownames_to_column("ESid")

# update ESid to reflect orientations
anno_df <- anno_df %>% mutate(ESid = paste0(Chr, ":", Start, ":", Ref, ":", Alt)) 

coverage_df$ESid <- anno_df$ESid
ratio_df$ESid <- anno_df$ESid

# create BED file for pileup
anno_bed <- anno_df %>% select(Chr, Start, End, ESid) %>%
    mutate(Start = Start - 1)

write_tsv(anno_bed, file = outBed, col_names = FALSE)

write_tsv(coverage_df, file = covMatOutFilt)
write_tsv(ratio_df, file = ratMatOutFilt)
write_tsv(anno_df, file = annotationsOut)

