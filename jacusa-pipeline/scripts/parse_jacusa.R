#!/usr/bin/env Rscript

library(optparse)
library(splitstackshape)
library(dplyr)
library(tidyverse)

option_list <- list(make_option(c('--input'), help = '', default = ''),
                    make_option(c('--output'), help = '', default = ''),
                    make_option(c('--minRate'), default = 0.01, help = 'minimum editing rate'),
                    make_option(c('--altDepth'), default = 3, help = 'required read coverage for alternative allele'),
                    make_option(c('--siteDepth'), default = 10, help = 'required read coverage of the locus'))

option.parser <- OptionParser(option_list = option_list)
opt <- parse_args(option.parser)

input <- opt$input
output <- opt$output
altDepth <- opt$altDepth
siteDepth <- opt$siteDepth
minRate <- opt$minRate
sampleFile <- basename(input)
sampleID <- gsub(".out", "", sampleFile)

message("Loading Jacusa2 output from ",input)
if (!file.exists(input)) stop("File ",input," does not exist")

# JH refactoring into tidyverse
df <- read_tsv(input, skip = 1)

# split base vector into 4
message(" * splitting counts by base")

df <- tidyr::separate(df, col = bases11, into = c("A","C","G","T"), sep = ",", convert = TRUE)

df <- select(df, chr =  `#contig`, pos = end, score, ref, "A", "C", "G", "T") 

# sum coverage across bases to get total
df$total_cov <- rowSums(df[,5:8])

df <- df %>% 
    mutate(ref_cov = case_when(
        ref == "A" ~ .[[5]], 
        ref == "C" ~ .[[6]], 
        ref == "G" ~ .[[7]], 
        ref == "T" ~ .[[8]]
        )) %>%
    mutate(chrpos = paste0(chr, ":", pos) ) %>%
    select(chrpos, score, ref, "A", "C", "G", "T", ref_cov, total_cov)

message(" * ", nrow(df), " sites loaded")
# pivot
df_long <- df %>% pivot_longer(!c(chrpos, ref, score, ref_cov, total_cov), names_to = "alt", values_to = "alt_cov") 
# remove rows where ref & alt are the same, remove alt with 0 coverage
df_long <- filter(df_long, ref != alt & alt_cov > 0)
# sort and remove multiallelic sites
df_long2 <- df_long[order(df_long$chrpos, -abs(df_long$"alt_cov") ), ]
multiallelic <- df_long$chrpos[duplicated(df_long$chrpos)]
df_long3 <- df_long2[!df_long2$chrpos %in% multiallelic,]

#require editing site to be covered by >=10 (default) reads and >=3 (default) reads covering the edited allele
df_filt <- 
    df_long3 %>%
    filter(total_cov >= siteDepth & alt_cov >= altDepth) %>%
    mutate( 
        ESid = paste0(chrpos, ":", ref, ":", alt),
        edit_rate = alt_cov / total_cov 
        ) %>% 
    filter(edit_rate >= minRate) %>%
    separate(col = chrpos, into = c("chr", "pos"), sep = ":") %>%
    select(ESid, chr, pos, score, total_cov, ref_cov, alt_cov, edit_rate)

message(" * ", nrow(df_filt), " sites pass thresholds")

write_tsv(df_filt, file = output)
