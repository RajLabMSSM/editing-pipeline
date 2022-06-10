#!/usr/bin/env Rscript

library(optparse)
library(dplyr)

option_list <- list(make_option(c('--inAnno'), help = 'Name of the annovar annotated file with no common SNPs', default = ''),
                    make_option(c('--inRat'), help = 'Name of the ratio matrix file before SNP filtering', default = ''),
                    make_option(c('--inCov'), help = 'Name of the coverage matrix file before SNP filtering', default = ''),
                    make_option(c('--outRat'), help = 'Name of the ratio matrix final filtered output file', default = ''),
                    make_option(c('--outCov'), help = 'Name of the coverage matrix final filtered output file', default = ''),
                    make_option(c('--outAnno'), help = 'Name of the output file with annotations for editing sites', default = ''))

option.parser <- OptionParser(option_list = option_list)
opt <- parse_args(option.parser)

noCommon <- opt$inAnno
covMatRaw <- opt$inCov
ratMatRaw <- opt$inRat
annotationsOut <- opt$outAnno
covMatOutFilt <- opt$outCov
ratMatOutFilt <- opt$outRat

annoOut <- read.delim(noCommon, header = T, sep = "\t")
anno <- annoOut[,c(1:7,26)]
colnames(anno) <- c("Chr", "Start", "Stop", "Ref", "Alt", "Region", "Gene", "Repeat")
anno <- anno %>% mutate(ESid = paste0(Chr, ":", Start, ":", Ref, ":", Alt))
ESannotated <- anno[,c(9,6,7,8)]

covMatfinal <- read.delim(covMatRaw,row.names = F, col.names = T, sep = "\t")
ratioMatfinal <- read.delim(ratMatRaw, row.names = F, col.names = T, sep = "\t")

coverageMatrixFinal <- subset(covMatfinal, (covMatfinal$ESid %in% ESannotated$ESid))
ratioMatrixFinal <- subset(ratioMatfinal, (ratioMatfinal$ESid %in% ESannotated$ESid))
ESannotatedFinal <- subset(ESannotated, (ESannotated$ESid %in% intersect(ratioMatrixFinal$ESid, coverageMatrixFinal$ESid)))

write.table(coverageMatrixFinal, file = covMatOutFilt, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(ratioMatrixFinal, file = ratMatOutFilt, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(ESannotatedFinal, file = annotationsOut, append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

