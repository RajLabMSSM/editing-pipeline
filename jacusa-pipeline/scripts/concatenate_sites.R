library(optparse)
library(dplyr)
library(tidyverse)

option_list <- list(make_option(c('--deNovoIN'), help = 'Formatted and filtered de novo called editing sites file', default = ''),
                    make_option(c('--knownIN'), help = 'Formatted and filtered known editing sites call file', default = ''),
                    make_option(c('--output'), help = '', default = ''))

option.parser <- OptionParser(option_list = option_list)
opt <- parse_args(option.parser)

deNovo <- opt$deNovoIN
knownES <- opt$knownIN
output <- opt$output

message("Loading Jacusa2 de novo output from ",deNovo)
if (!file.exists(deNovo)) stop("File ",deNovo," does not exist")
message("Loading Jacusa2 known sites output from ",knownES)
if (!file.exists(knownES)) stop("File ",knownES," does not exist")

df_denovo <- read_tsv(deNovo)
df_known <- read_tsv(knownES)
df_all <- dplyr::union(df_denovo, df_known)

write_tsv(df_all, file = output)
