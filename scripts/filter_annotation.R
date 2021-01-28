library(optparse)
library(tidyverse)
library(rtracklayer)


option_list <- list(
    make_option(c('--input'), help='', default = "example"),
        make_option(c('--output'), help='', default = "example")
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)


input <- opt$input
output <- opt$output



gtf_file <- "/sc/hydra/projects/ad-omics/data/references/hg38_reference/GENCODE/gencode.v32.annotation.gtf.gz"

gtf <- import.gff(gtf_file)

gene_table <- data.frame(gene_id = gtf$gene_id, strand = strand(gtf), stringsAsFactors = FALSE )

gene_table$gene_id <- str_split_fixed( gene_table$gene_id, "\\.", 2)[,1]

annotations <- read_tsv(input)

anno_df <- left_join(annotations, gene_table, by = c("ANN[*].GENEID" = "gene_id") )

anno_df <- 
    mutate(anno_df, 
    verdict = case_when( 
        REF == "G" & ALT == "A" & strand == "-" ~ "PASS", 
        REF == "C" & ALT == "T" & strand == "+" ~ "PASS",
        TRUE ~ "FAIL")
    )

anno_df <- filter(anno_df, verdict == "PASS")

write_tsv(anno_df, path = output)
