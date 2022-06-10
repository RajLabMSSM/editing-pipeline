library(dplyr)

annoOut <- read.delim("~/Downloads/AnswerALS.myanno.hg38_multianno.txt.noCommon.txt", sep = "\t", header = T)
anno <- annoOut[,c(1:7,26)]
colnames(anno) <- c("Chr", "Start", "Stop", "Ref", "Alt", "Region", "Gene", "Repeat")
anno <- anno %>% mutate(ESid = paste0(Chr, ":", Start, ":", Ref, ":", Alt))
ESannotated <- anno[,c(9,6,7,8)]

load("~/Downloads/denovofiles/ratioMat.Rda")
load("~/Downloads/denovofiles/covMat.Rda")
coverageMatrixFinal <- subset(covMatfinal, (covMatfinal$ESid %in% ESannotated$ESid))
ratioMatrixFinal <- subset(ratioMatfinal, (ratioMatfinal$ESid %in% ESannotated$ESid))
ESannotatedFinal <- subset(ESannotated, (ESannotated$ESid %in% intersect(ratioMatrixFinal$ESid, coverageMatrixFinal$ESid)))
save(coverageMatrixFinal, file = "FinalCoverageMatrix.Rda")
save(ratioMatrixFinal, file = "FinalRatioMatrix.Rda")
save(ESannotatedFinal, file = "FinalAnnotationsMatrix.Rda")