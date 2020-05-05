library(tidyverse)

library(optparse)

# -o specifies output file
# positional arguments for each input sites.bed file

sites <- list.files(pattern = "sites.snp_filtered.bed", recursive = TRUE, full.names = TRUE)

all_sites <- purrr::map(sites, ~{
                readr::read_tsv(.x, col_names = c("chr","start","end", "info", "score", "strand")) %>%
                dplyr::distinct() %>%
                as.data.frame()  
                })

sample_names <- gsub(".sites.snp_filtered.bed", "", basename(sites))

info_df <- 
    purrr::map2( .x = all_sites, .y = sample_names, 
            ~{
            df = mutate(.x, index = paste0(chr,":", start,"-",end,":",strand, sep = ":") ) %>% 
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
x <- 0.75

clean_sites <- rowSums( !is.na(info_df) ) >= floor( x * ncol(info_df) )

clean_df <- info_df[ clean_sites,]

# get out editing ratios
editing_df <- apply(clean_df, MARGIN = c(1,2), FUN = function(x) {as.numeric(split_info(x, i = 3))})

# calculate mean and max editing rates
mean_editing <- rowMeans(editing_df, na.rm=TRUE)
max_editing <- apply( editing_df, MAR = 1, FUN = function(x) max(x, na.rm=TRUE) )

# hypothesis - C>T editing in homeostatic immune cells is going to be low ~ 5-10% at best
editing_df[ max_editing > 0.01 &  mean_editing < 0.25,]

# output a VCF for annotation
# VCF columns: #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  MSMD0170-01

# ATTENTION - THESE VCFs are off by 1 - fix this
createVCF <- function(df){
    index <- row.names(df)
    index_df <- as.data.frame(str_split_fixed(index, ":|-", 4))
    names(index_df) <- c("#CHROM", "POS", "END", "STRAND")
    index_df$ID <- index
    index_df$REF <- ifelse(index_df$STRAND == "+", "C", "G")
    index_df$ALT <- ifelse(index_df$STRAND == "+", "T", "A")
    # quality can be number of samples with non-missing entries
    index_df$QUAL <- rowSums( !is.na(df) )
    index_df$FILTER <- "."
    index_df$INFO <- index
    index_df$FORMAT <- "."

    index_df <- select(index_df, `#CHROM`, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT)
    #return(index_df)
    vcf <- cbind(index_df, df)
    return(vcf)
}
