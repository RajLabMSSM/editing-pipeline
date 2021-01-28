library(tidyverse)

library(optparse)

# -o specifies output file
# positional arguments for each input sites.bed file

=======

option_list <- list(
    make_option(c('--prefix'), help = 'prefix to add to search path, if required', default = ''),
    make_option(c('--dataCode'), help='', default = "example"),
        make_option(c('--missingness'), help='', default = 0.8)
)

option.parser <- OptionParser(option_list=option_list)
opt <- parse_args(option.parser)


dataCode <- opt$dataCode
missingness <- opt$missingness
prefix <- opt$prefix
outRData <- paste0(dataCode, "_merged_sites.RData")
outVCF <- paste0(dataCode, "_merged_sites.vcf")

# -o specifies output file
# positional arguments for each input sites.bed file

# if prefix variable == "all" then use all files for merging
if( prefix == "all"){ prefix <- "" }

search_pattern <- paste0(prefix,".*sites.snp_filtered.bed")

message( " * Using search pattern:", search_pattern)

sites <- list.files(pattern = search_pattern, recursive = TRUE, full.names = TRUE)

message( " * ", length(sites), " files detected")
stopifnot(length(sites) > 0)
all_sites <- purrr::map(sites, ~{
                message(" * reading in ",.x)
                readr::read_tsv(.x, col_names = c("chr","start","end", "info", "score", "strand"), col_types = "cnncnc") %>%
                dplyr::distinct() %>%
                as.data.frame()  
                })

sample_names <- gsub(".sites.snp_filtered.bed", "", basename(sites))

info_df <- 
    purrr::map2( .x = all_sites, .y = sample_names, 
            ~{
            message(" * processing ", .y) 
            df <- .x %>%
            dplyr::mutate(index = paste0(chr,":", start,"-",end,":",strand) ) %>% 
            dplyr::mutate(info = paste0(info, "|", score) ) %>%
            dplyr::select(index, info)
        names(df)[2] <- .y
        df 
    } )

rm(all_sites)

# with lots of samples, full joining each is impossible
# work in batches instead?

# get all indexes out and create unique list
# for each sample leftjoin with the index list
all_index <- unique(unlist(map(info_df, "index")))
index_df <- data.frame(index = all_index, stringsAsFactors=FALSE)

# join each sample to the master list
all_index <- map( info_df, ~{ 
    message(" * joining ", colnames(.x)[2] )
    left_join(index_df, .x, by = "index") %>% 
    column_to_rownames(var = "index") 
})

# bind together
info_df <- do.call( cbind, all_index)


# info_df is huge now - only retain sites found in at least 25% samples for saving
sites_5 <- rowSums( !is.na(info_df) ) >= floor(0.25 * ncol(info_df) )
info_df <- info_df[ sites_5,]

# select sites covered in X% of samples
#missingness <- 0.5

clean_sites <- rowSums( !is.na(info_df) ) >= floor( missingness * ncol(info_df) )

clean_df <- info_df[ clean_sites,]


# create matrix of editing rates
split_info <- function(info, i){
    unlist(stringr::str_split_fixed(info, "\\|", 4)[,i])
}


# get out editing ratios
editing_df <- apply(clean_df, MARGIN = c(1,2), FUN = function(x) {as.numeric(split_info(x, i = 3))})

# get confidence scores for each site 
confidence_df <-  apply(clean_df, MARGIN = c(1,2), FUN = function(x) {as.numeric(split_info(x, i = 4))})

# hypothesis - C>T editing in homeostatic immune cells is going to be low ~ 5-10% at best
#editing_df[ max_editing > 0.01 &  mean_editing < 0.25,]

# output a VCF for annotation
# VCF columns: #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  MSMD0170-01

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

save(info_df, editing_df, confidence_df, file = outRData)

vcf <- createVCF(editing_df)

# add in VCF header so IGV recognises it
vcf_header <- "##fileformat=VCFv4.3"
writeLines(vcf_header, con = outVCF)
write_tsv(vcf, path = outVCF, append = TRUE, col_names = TRUE) 
