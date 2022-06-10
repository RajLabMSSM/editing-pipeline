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
covMatOut <- opt$cov
ratMatOut <- opt$rat
avOut <- opt$av
X <- opt$percSamples
Y <- opt$minER

#aggregating across samples
files <- dir(path = inDir, pattern = "*.filt")
data <- files %>% map(read_tsv)
aggDF <- reduce(data, full_join, by = "ESid")
covMat <- map(data, ~{
  d <- select(.x, ESid, totcov)
  colnames(d)[2] <- colnames(.x)[5]
  return(d)
}) %>% reduce(full_join, by = "ESid")

ratioMat <- map(data, ~{
  .x$edrat <- .x[,5]/.x[,6]
  d <- select(.x, ESid, edrat)
  colnames(d)[2] <- colnames(.x)[5]
  return(d)
}) %>% reduce(full_join, by = "ESid")

#cohort level-filtering: editing sites must validate across 50% of samples and have editing efficiency of at least 10% with default settings
N <- ceiling(length(files) * X) #where the --percSamples flag comes in
covMat <- covMat[which(rowSums(!is.na(covMat[,-1])) >= N),]
ratioMat <- ratioMat[which(rowMeans(ratioMat[,-1]) >= Y),] #where the --minER flag comes in

covMatfinal <- subset(covMat, (covMat$ESid %in% intersect(covMat$ESid, ratioMat$ESid)))
ratioMatfinal <- subset(ratioMat, (ratioMat$ESid %in% intersect(covMat$ESid, ratioMat$ESid)))

write.table(ratioMatfinal, file = ratMatOut, quote = F, sep = "\t", row.names = F, col.names = T)
write.table(covMatfinal, file = covMatOut, quote = F, sep = "\t", row.names = F, col.names = T)

#formatting for annovar
annovar <- data.frame(ratioMatfinal$ESid, do.call(rbind, strsplit(ratioMatfinal$ESid, split = ":", fixed = TRUE)))
annovar <- annovar[,c(2,3,3,4,5)]
colnames(annovar) <- c("Chr", "Start", "Stop", "Ref", "Alt")

write.table(annovar, file = avOut, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)



