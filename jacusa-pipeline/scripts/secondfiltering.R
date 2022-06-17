#!/usr/bin/env Rscript

library(optparse)
library(purrr)
library(dplyr)
library(tidyverse)

option_list <- list(make_option(c('--inDir'), help = 'The path to the Jacusa output directory', default = ''),
                    make_option(c('--rat'), help = 'Name of the ratio matrix output file', default = ''),
                    make_option(c('--cov'), help = 'Name of the coverage matrix output file', default = ''),
                    make_option(c('--av'), help = 'Name of the vcf file for annotation in annovar', default = ''),
                    make_option(c('--percSamples'), default = 0.5, help = '% of samples editing site is required to validate across'),
                    make_option(c('--minER'), default = 0.1, help = 'minimum editing ratio'))

option.parser <- OptionParser(option_list = option_list)
opt <- parse_args(option.parser)

inDir <- opt$inDir
coverage_df_out <- opt$cov
ratio_df_out <- opt$rat
annovar_out <- opt$av
perc_samples <- opt$percSamples
min_edrate <- opt$minER

#aggregating across samples
#inDir <- "/sc/arion/projects/ad-omics/flora/cohorts/Navarro_stim/editing-pipeline/editing-pipeline/jacusa-pipeline/snakemake-workspace/"
#path = inDir
#print(path)

#message("inDirectory") 
#print(inDir)


files <- list.files(path = inDir, pattern = "*.filt$",full.names = TRUE, recursive = TRUE)

#files <- head(files,10)

message(" * found ", length(files) , " jacusa files")
#print(files)
#stop()

message(" * reading files")
sample_ids <- gsub(".filt", "", basename(files) )

data <- files %>% map(read_tsv, col_types = "cnnnnn")
names(data) <- sample_ids

message(" * merging files!")

#aggDF <- reduce(data, full_join, by = "ESid")

# get list of  all sites but remove singletons!
all_sites <- map(data, "ESid") %>% unlist() 
all_sites <- all_sites[duplicated(all_sites)]

message(" * ", length(all_sites), " unique sites found")

# coverage matrix - total site coverage in each sample 
joint_sites <- map2(data, sample_ids, ~{
  print(.y)
  d <- column_to_rownames(.x, "ESid")
  return(d[all_sites,])
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

print(head(ratio_df))

message(" * ", nrow(coverage_df), " sites found in total")
#save.image("debug.RData")
#stop()
#cohort level-filtering: editing sites must validate across 50% of samples and have editing efficiency of at least 10% with default settings
sample_n <- ceiling(length(files) * perc_samples) #where the --percSamples flag comes in

coverage_df_filt <- coverage_df[which(rowSums(!is.na(coverage_df)) >= sample_n),]

message( " * ", nrow(coverage_df_filt), " sites are found in at least ", perc_samples * 100, "% of samples (", sample_n, ")" )

ratio_df_filt <- ratio_df[which(rowMeans(ratio_df, na.rm = TRUE) >= min_edrate),] #where the --minER flag comes in

message( " * ", nrow(ratio_df_filt), " sites have at least ", min_edrate * 100, "% mean editing rate" )

clean_sites <- intersect(row.names(coverage_df_filt), row.names(ratio_df_filt) )

clean_sites <- clean_sites[order(clean_sites) ]

message( " * keeping ", length(clean_sites), " editing sites total" )

coverage_df_final <- coverage_df[clean_sites,] %>% rownames_to_column("ESid")

ratio_df_final <- ratio_df[clean_sites,] %>% rownames_to_column("ESid")

#coverage_df_final <- arrange(coverage_df_final, ESid)
#ratio_df_final <- arrange(ratio_df_final, ESid)

write_tsv(coverage_df_final, file = coverage_df_out)
write_tsv(ratio_df_final, file = ratio_df_out)

#formatting for annovar
annovar <- tibble(ESid = ratio_df_final$ESid) %>%
    tidyr::separate(col = ESid, into = c("Chr", "Stop", "Ref", "Alt"), sep = ":") %>%
    mutate(Start = Stop) %>%
    select(Chr, Start, Stop, Ref, Alt)

#annovar <- data.frame(ratio_df_final$ESid, do.call(rbind, strsplit(ratio_dffinal$ESid, split = ":", fixed = TRUE)))
#annovar <- annovar[,c(2,3,3,4,5)]
#colnames(annovar) <- c("Chr", "Start", "Stop", "Ref", "Alt")

message(" * writing out in VCF format for Annovar")
print(head(annovar) )

write_tsv(annovar, file = annovar_out, col_names = FALSE)

#, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)



