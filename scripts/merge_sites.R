library(tidyverse)

library(optparse)


option_list <- list(
    make_option(c('--dataCode'), help='', default = "example"),
        make_option(c('--missingness'), help='', default = 0.8)
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)


dataCode <- opt$dataCode
missingness <- opt$missingness

outRData <- paste0(dataCode, "_merged_sites.RData")
outVCF <- paste0(dataCode, "_merged_sites.vcf")

# -o specifies output file
# positional arguments for each input sites.bed file

sites <- list.files(pattern = "sites.snp_filtered.bed", recursive = TRUE, full.names = TRUE)

all_sites <- purrr::map(sites, ~{
                readr::read_tsv(.x, col_names = c("chr","start","end", "info", "score", "strand"), col_types = "cnncnc") %>%
                dplyr::distinct() %>%
                as.data.frame()  
                })

sample_names <- gsub(".sites.snp_filtered.bed", "", basename(sites))

info_df <- 
    purrr::map2( .x = all_sites, .y = sample_names, 
            ~{
            df = mutate(.x, index = paste0(chr,":", start,"-",end,":",strand) ) %>% 
            dplyr::select(index, info)
        names(df)[2] <- .y
        df 
    } ) %>% 
purrr::reduce( dplyr::left_join, by = "index" )

info_df <- column_to_rownames(info_df, var = "index")

# create matrix of editing rates
split_info <- function(info, i){
    unlist(stringr::str_split_fixed(info, "\\|", 3)[,i])
}

# select sites covered in X% of samples
#missingness <- 0.75

clean_sites <- rowSums( !is.na(info_df) ) >= floor( missingness * ncol(info_df) )

clean_df <- info_df[ clean_sites,]

# get out editing ratios
editing_df <- apply(clean_df, MARGIN = c(1,2), FUN = function(x) {as.numeric(split_info(x, i = 3))})

# calculate mean and max editing rates
mean_editing <- rowMeans(editing_df, na.rm=TRUE)
max_editing <- apply( editing_df, MAR = 1, FUN = function(x) max(x, na.rm=TRUE) )

# hypothesis - C>T editing in homeostatic immune cells is going to be low ~ 5-10% at best
#editing_df[ max_editing > 0.01 &  mean_editing < 0.25,]

# output a VCF for annotation
# VCF columns: #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  MSMD0170-01

# ATTENTION - THESE VCFs are off by 1 - fix this
createVCF <- function(df){
    index <- row.names(df)
    mean_editing <- rowMeans(df, na.rm=TRUE)
    max_editing <- apply( df, MAR = 1, FUN = function(x) max(x, na.rm=TRUE) )
    index_df <- as.data.frame(str_split_fixed(index, ":|-", 4))
    names(index_df) <- c("#CHROM", "START", "POS", "STRAND")
    index_df$ID <- index
    index_df$REF <- ifelse(index_df$STRAND == "+", "C", "G")
    index_df$ALT <- ifelse(index_df$STRAND == "+", "T", "A")
    # quality can be number of samples with non-missing entries
    index_df$QUAL <- rowSums( !is.na(df) )
    index_df$FILTER <- mean_editing
    index_df$INFO <- max_editing
    index_df$FORMAT <- "."

    index_df <- select(index_df, `#CHROM`, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT)
    #return(index_df)
    vcf <- cbind(index_df, df)
    return(vcf)
}

save(info_df, editing_df, file = outRData)

vcf <- createVCF(editing_df)

write_tsv(vcf, path = outVCF) 
