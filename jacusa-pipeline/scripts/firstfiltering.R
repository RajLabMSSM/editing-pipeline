#!/usr/bin/env Rscript

library(optparse)
library(splitstackshape)
library(dplyr)
library(tidyverse)

option_list <- list(make_option(c('--input'), help = '', default = ''),
                    make_option(c('--output'), help = '', default = ''),
                    make_option(c('--altDepth'), default = 2, help = 'required read coverage for alternative allele'),
                    make_option(c('--siteDepth'), default = 10, help = 'required read coverage of the locus'))

option.parser <- OptionParser(option_list = option_list)
opt <- parse_args(option.parser)

input <- opt$input
output <- opt$output
altDepth <- opt$altDepth
siteDepth <- opt$siteDepth
sampleFile <- basename(input)
sampleID <- gsub(".out", "", sampleFile)

cat("Loading Jacusa2 output from",input,"\n")
if (!file.exists(input)) stop("File ",input," does not exist")
out <- read.delim(input, header = T, sep = "\t")

OUT <- cSplit(out, "bases11",",", direction = "wide") #split base vector into four columns
OUT <- OUT[,-c(2,4,6,7,8)] #"end" pos is true location so drop "start", "name" was a place holder used by Jacusa2, strand&info&filter are not used downstream
colnames(OUT) <- c("chr", "pos", "score", "ref", "A", "C", "G", "T") #assigning new column names to make obvious the order of the bases vector
OUT$totcov <- rowSums(OUT[,c(5:8)]) #creating total coverage column for easy filtering downstream
DF <- OUT %>% mutate(refcov = case_when(ref == "T" ~ .[[8]], ref == "A" ~ .[[5]], ref == "G" ~ .[[7]], ref == "C" ~ .[[6]]), chrpos = paste0(chr, ":", pos)) #creating chr:pos tag & ref cov columns
DF <- DF[,-c(1:2)] #dropping chr & pos separated
DF <- DF[,c(9,2,1,7,8,3,4,5,6)] #tidying column order
DFnew <- DF %>% pivot_longer(!c(chrpos, ref, score, refcov, totcov), names_to = "alt", values_to = "cov") #converting from wide to long format
DFnew <- DFnew[which(DFnew$ref != DFnew$alt & DFnew$totcov >= siteDepth & DFnew$cov > altDepth, ),] #remove rows where ref & alt are the same, remove alt with 0 coverage, require editing site to be covered by >=10 (default) reads and >=2 (defaul) reads covering the edited allele
DFnewer <- DFnew[!duplicated(DFnew$chrpos),] #drop multiallelic sites
df <- DFnewer %>% mutate(ESid = paste0(chrpos, ":", ref, ":", alt)) #creating ESid for sample-level aggregation
colnames(df) <- c("chrpos", "ref", "score", "totcov", "refcov", "alt", sampleID, "ESid")
colnames(df) <- gsub(pattern = ".out", replacement = "", colnames(df))
temp <- df[,c(8,2,5,6,7,4,3)]

write.table(temp, file = output, quote = F, sep = "\t", row.names = F, col.names = T)
