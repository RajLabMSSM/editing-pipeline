library(splitstackshape)
library(dplyr)
library(tidyverse)
library(plyr)
library(optparse)

##establishing options for Rscript
option_list <- list(make_option(c('--inFolder', '-i'), 
                                help = 'The full path to the directory containing Jacusa2 output files',
                                default = ""))
option.parser <- OptionParser(option_list = option_list)
opt <- parse_args(option.parser)
inFolder <- opt$inFolder
all_files <- list.files(inFolder, pattern = ".out", recursive = T, full.names = T)
message(paste0(" * found ", length(all_files), " Jacusa2 output files" ))

fileNames <- list.files(inFolder, pattern = ".out", recursive = T, full.names = T)
fileNumbers <- seq(fileNames)
message(paste0(" * found ", length(fileNames), " Jacusa2 output files "))

##for loop of filtering samples


#restructuring Jacusa2 output and light sample-level filtering
setwd("~/Downloads/denovofiles/")
extension <- "out"
fileNames <- Sys.glob(paste("*.", extension, sep = ""))
fileNumbers <- seq(fileNames)

for (fileNumber in fileNumbers){
  out <- read.delim(fileNames[fileNumber], header = T, sep = "\t") #read in Jacusa2 de novo caller output
  OUT <- cSplit(out, "bases11",",", direction = "wide") #split base vector into four columns
  OUT <- OUT[,-c(2,4,6,7,8)] #drop "start" because we only need "end" position which is true location, name column was a Jacusa2 place holder, strand is missing, info&filter are unimportant downstream
  colnames(OUT) <- c("chr", "pos", "score", "ref", "A", "C", "G", "T") #assigning new column names to make obvious the order of the bases vector
  OUT$totcov <- rowSums(OUT[,c(5:8)]) #creating total coverage column for easy filtering downstream
  DF <- OUT %>% mutate(refcov = case_when(ref == "T" ~ .[[8]], ref == "A" ~ .[[5]], ref == "G" ~ .[[7]], ref == "C" ~ .[[6]]), chrpos = paste0(chr, ":", pos)) #creating chr:pos tag & ref cov columns
  DF <- DF[,-c(1:2)] #dropping chr & pos separated
  DF <- DF[,c(9,2,1,7,8,3,4,5,6)] #tidying column order
  DFnew <- DF %>% pivot_longer(!c(chrpos, ref, score, refcov, totcov), names_to = "alt", values_to = "cov") #converting from wide to long format
  covThreshold <- 10
  DFnew <- DFnew[which(DFnew$ref != DFnew$alt & DFnew$totcov >= covThreshold & DFnew$cov > 2, ),] #remove rows where ref & alt are the same, require >10 reads coverage of the site, remove alt with 0 coverage, require alt to have >2 reads
  DFnewer <- DFnew[!duplicated(DFnew$chrpos),] #drop multiallelic sites
  df <- DFnewer %>% mutate(ESid = paste0(chrpos, ":", ref, ":", alt)) #creating ESid for sample-level aggregation
  colnames(df) <- c("chrpos", "ref", "score", "totcov", "refcov", "alt", fileNames[fileNumber], "ESid")
  colnames(df) <- gsub(pattern = ".out", replacement = "", colnames(df))
  temp <- df[,c(8,2,5,6,7,4,3)]
  
  #outfile <- paste0(inFolder, paste0(fileNames[fileNumber], ".Jacusa2filt"))
  #write.tsv(temp, path = outfile, quote = F, row.names = F, col.names = T)
  write.table(temp, file = paste0(fileNames[fileNumber], ".Jacusa2filt"), quote = F, sep = "\t", row.names = F, col.names = T)
}

#aggregating across samples
library(optparse)
library(purrr)
library(dplyr)
library(tidyverse)

files <- dir(path = , pattern = "*.Jacusa2filt")
data <- files %>% map(read_tsv) #%>% join_all(,by = "ESid", type = "full")
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

#cohort level-filtering: editing sites must validate across 50% of samples and have editing efficiency of at least 10%
N <- ceiling(length(fileNumbers) * 0.5) #set threshold for inter-donor validation
covMat <- covMat[which(rowSums(!is.na(covMat[,-1])) >= N),]
ratioMat <- ratioMat[which(rowMeans(ratioMat[,-1]) >= 0.1),] #set editing ratio threshold

covMatfinal <- subset(covMat, (covMat$ESid %in% intersect(covMat$ESid, ratioMat$ESid)))
ratioMatfinal <- subset(ratioMat, (ratioMat$ESid %in% intersect(covMat$ESid, ratioMat$ESid)))

save(covMatfinal, file = "~/Downloads/denovofiles/covMat.Rda") 
save(ratioMatfinal, file = "~/Downloads/denovofiles/ratioMat.Rda")










