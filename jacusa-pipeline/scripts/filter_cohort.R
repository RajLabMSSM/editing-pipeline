#!/usr/bin/env Rscript
# merge sites together
# Winston Cuddleston, Jack Humphrey, Hyomin Seo
# 2022

library(optparse)
library(purrr)
library(dplyr)
library(tidyverse)
library(vroom)

option_list <- list(make_option(c('--inDir'), help = 'The path to the Jacusa output directory', default = '.'),
                    make_option(c('--chr'), help = 'the chromosome to work on', default = 'chr21'),
                    make_option(c('--minSamples'), default = 60, help = '# of samples editing site is required to validate across'),
                    make_option(c('--minER'), default = 0.1, help = 'minimum editing ratio'),
                    make_option(c('--group'), help = 'metadata file', default = ''))

option.parser <- OptionParser(option_list = option_list)
opt <- parse_args(option.parser)

inDir <- opt$inDir
minSampleN <- opt$minSamples
min_edrate <- opt$minER
chrom <- opt$chr
group <- opt$group

cov_out <- paste0(inDir, "/filter/", chrom, "_coverage.tsv.gz")
rat_out <- paste0(inDir, "/filter/", chrom, "_ratio.tsv.gz")

#cohort level-filtering: within a defined group, editing sites must validate across 60 samples and have mean editing efficiency of at least 10% with default settings
message(" * reading in coverage file")
cov_in <- read_tsv(paste0(inDir, "/merge/", chrom, "_coverage.tsv.gz"), col_names = T)
coverage_df <- as.data.frame(cov_in)
rownames(coverage_df) <- cov_in$ESid
coverage_df <- select(coverage_df, -ESid)
coverage_DF <- data.frame(t(coverage_df))

group_in <- read_tsv(paste0(inDir, group), col_names = T)
group_df <- as.data.frame(group_in)
rownames(group_df) <- group_in$sample
group_df <- select(group_df, -c(sample, bam_path, library))

cov_DF <- merge(coverage_DF, group_df, by = 0)
grouped_covDFs <- cov_DF %>% group_by(group) %>% group_split()

coverageDF <- map(grouped_covDFs, ~{
  xDF <- .x %>% select(-Row.names, -group)
  xT <- data.frame(t(xDF))
  colnames(xT) <- .x$Row.names
  xFilt <- xT[which(rowSums(!is.na(xT)) >= minSampleN),] #where the --minSamples flag comes in
  xFilt <- rownames_to_column(xFilt, "ESid")
  return(xFilt)}) %>% reduce(full_join, by = "ESid")

message( " * ", nrow(coverageDF), " sites are found in at least ", minSampleN, "of each group" )

covThresholdSites <- gsub("\\.",":",coverageDF$ESid)

message(" * reading in ratio file")
rat_in <- read_tsv(paste0(inDir, "/merge/", chrom, "_ratio.tsv.gz"), col_names = T)
rat_in_cleaned <- subset(rat_in, (rat_in$ESid %in% covThresholdSites))
ratio_df <- as.data.frame(rat_in_cleaned)
rownames(ratio_df) <- rat_in_cleaned$ESid
ratio_df <- select(ratio_df, -ESid)
ratio_DF <- data.frame(t(ratio_df))

rat_DF <- merge(ratio_DF, group_df, by = 0)
grouped_ratDFs <- rat_DF %>% group_by(group) %>% group_split()

ratioDF <- map(grouped_ratDFs, ~{
  xDF <- .x %>% select(-Row.names, -group)
  xT <- data.frame(t(xDF))
  colnames(xT) <- .x$Row.names
  xFilt <- xT[which(rowMeans(xT, na.rm = TRUE) >= min_edrate),] #this is where the --minER flag comes in
  xFilt <- rownames_to_column(xFilt, "ESid")
  return(xFilt)}) %>% reduce(full_join, by = "ESid")

message( " * ", nrow(ratioDF), " sites have at least ", min_edrate * 100, "% mean editing rate within each group" )

clean_sites <- intersect(ratioDF$ESid, coverageDF$ESid)

message( " * keeping ", length(clean_sites), " editing sites total" )

clean_sites <- clean_sites[order(clean_sites)]

coverageDF <- subset(coverageDF, (coverageDF$ESid %in% clean_sites))
ratioDF <- subset(ratioDF, (ratioDF$ESid %in% clean_sites))

coverageDF$ESid <- gsub("\\.",":",coverageDF$ESid)
ratioDF$ESid <- gsub("\\.",":",ratioDF$ESid)

write_tsv(coverageDF, file = paste0(inDir, "/filter/", chrom, "_coverage.tsv.gz"))
write_tsv(ratioDF, file = paste0(inDir, "/filter/", chrom, "_ratio.tsv.gz"))

