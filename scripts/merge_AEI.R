# merge all Alu Editing Index outputs together
# Jack Humphrey 2021

library(tidyverse)
library(optparse)

option_list <- list(
    make_option(c('--inFolder', '-i' ), help='The full path to the folder that contains the subfolders', default = "")
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)

inFolder <- opt$inFolder

all_files <- list.files(inFolder, pattern = "EditingIndex.csv", recursive = TRUE, full.names = TRUE)

message(paste0(" * found ", length(all_files), " AEI files" ))

all_res <- purrr::map_df(all_files, read_csv)

outfile <- paste0(inFolder, "/all_aei_results.tsv")

write_tsv(all_res, path = outfile)
