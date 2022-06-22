#!/usr/bin/env Rscript
# merge sites together
# Winston Cuddleston, Jack Humphrey, Hyomin Seo
# 2022

library(optparse)
library(purrr)
library(dplyr)
library(tidyverse)

option_list <- list(make_option(c('--inDir'), help = 'The path to the Jacusa output directory', default = '.'),
                    make_option(c('--anno'), help = 'The site annotation following annovar', default = ''))

option.parser <- OptionParser(option_list = option_list)
opt <- parse_args(option.parser)

inDir <- opt$inDir
anno_file <- opt$anno

anno_out <- paste0(inDir,"all_sites_pileup_annotation.tsv.gz")
cov_out <- paste0(inDir, "all_sites_pileup_coverage.tsv.gz")
rat_out <- paste0(inDir, "all_sites_pileup_editing.tsv.gz")

# read in annotation file to get known ESids
# currently some ESids repeated on multiple lines due to GENCODE join
anno_df <- read_tsv(anno_file)
stopifnot( nrow(anno_df) == length(unique(anno_df$ESid) ) )

# sites are ordered by appearance in annotation
all_sites <- anno_df$ESid

files <- list.files(path = inDir, pattern = "*pileup_parsed.txt$",full.names = TRUE, recursive = TRUE)

message(" * found ", length(files) , " jacusa pileup files")

message(" * reading files")
sample_ids <- gsub("_pileup_parsed.txt", "", basename(files) )

# read in all files together as list
all_data <-  map(files,~{
        read_tsv(.x, col_types = "ccnnnnnn")
    })

names(all_data) <- sample_ids

message(" * merging files!")

message(" * merging ", length(all_sites), " unique sites")

# coverage matrix - total site coverage in each sample 
joint_sites <- map(all_data, ~{
    left_join(data.frame(ESid = all_sites), distinct(.x), by = "ESid")
    })

## now just extract columns and bind - no fancy join needed
coverage_df <- map2(joint_sites, sample_ids, ~{
    d <- select(.x, total_cov)
    names(d)[1] <- .y
    return(d)
}) %>% reduce(cbind)

row.names(coverage_df) <- all_sites

ratio_df <- map2(joint_sites, sample_ids, ~{
        d <- select(.x, edit_rate)
        names(d)[1] <- .y
        return(d)
}) %>% reduce(cbind)
row.names(ratio_df) <- all_sites

ratio_df <- rownames_to_column(ratio_df, "ESid")
coverage_df <- rownames_to_column(coverage_df, "ESid")

# flip orientation of bases if annotated to negative strand gene
revcomp <- function(x){
  data.frame(y = c("A","C","G","T"),row.names =  c("T","G","C","A"))[x,]
}
message(" * flipping orientations" )

anno_df <- mutate(anno_df,
                  Ref = ifelse(!is.na(strand) & strand == "-", revcomp(Ref), Ref ),
                  Alt = ifelse(!is.na(strand) & strand == "-", revcomp(Alt), Alt )
                  )
# update ESids on annotation, coverage and ratio matrices
anno_df <- anno_df %>% mutate(ESid = paste0(Chr, ":", Start, ":", Ref, ":", Alt))

coverage_df$ESid <- anno_df$ESid
ratio_df$ESid <- anno_df$ESid

# write out
write_tsv(coverage_df, file = cov_out)
write_tsv(ratio_df, file = rat_out)
write_tsv(anno_df, file = anno_out)

