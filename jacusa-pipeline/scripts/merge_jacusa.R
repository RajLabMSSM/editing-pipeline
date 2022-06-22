#!/usr/bin/env Rscript
# merge sites together
# Winston Cuddleston, Jack Humphrey, Hyomin Seo
# 2022

library(optparse)
library(purrr)
library(dplyr)
library(tidyverse)

option_list <- list(make_option(c('--inDir'), help = 'The path to the Jacusa output directory', default = '.') )

option.parser <- OptionParser(option_list = option_list)
opt <- parse_args(option.parser)

inDir <- opt$inDir
coverage_df_out <- opt$cov
ratio_df_out <- opt$rat
annovar_out <- opt$av
perc_samples <- opt$percSamples
min_edrate <- opt$minER
chromosome <- opt$chr
#aggregating across samples

temp_cov <- paste0(inDir, "/all_sites_coverage.tsv.gz")
temp_rat <- paste0(inDir, "/all_sites_ratio.tsv.gz")

files <- list.files(path = inDir, pattern = "*.filt$",full.names = TRUE, recursive = TRUE)


message(" * found ", length(files) , " jacusa files")

message(" * reading files")
sample_ids <- gsub(".filt", "", basename(files) )

chroms <- paste0("chr", c(1:22, "X","Y","M") )

# read in all files together as list
all_data <-  map(files,~{
        read_tsv(.x, col_types = "ccnnnnnn")
    })

# iterate through chromosomes
for(chromosome in chroms){
    print(chromosome)
    
    data <- map(all_data,~{
        filter(.x, chr == chromosome)
    })

    names(data) <- sample_ids

    message(" * merging files!")


# get list of  all sites but remove singletons!
    all_sites <- map(data, "ESid") %>% unlist() 
    all_sites <- all_sites[duplicated(all_sites)]
    all_sites <- unique(all_sites)

    message(" * ", length(all_sites), " unique sites found")

# coverage matrix - total site coverage in each sample 
    joint_sites <- map2(data, sample_ids, ~{
#    print(.y)
        left_join(data.frame(ESid = all_sites), .x, by = "ESid")
    })

    message(" * filling matrices ")

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
    
    write_tsv(coverage_df, file = temp_cov, col_names = chromosome == "chr1", append = chromosome != "chr1")
    write_tsv(ratio_df, file = temp_rat, col_names = chromosome == "chr1", append = chromosome != "chr1")
    gc()
}


