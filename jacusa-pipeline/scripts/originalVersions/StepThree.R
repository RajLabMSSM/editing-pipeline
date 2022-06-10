load("~/Downloads/denovofiles/ratioMat.Rda")
load("~/Downloads/denovofiles/covMat.Rda")

annovar <- data.frame(ratioMat$ESid, do.call(rbind, strsplit(ratioMat$ESid, split = ":", fixed = TRUE)))
annovar <- annovar[,c(2,3,3,4,5)]
colnames(annovar) <- c("Chr", "Start", "Stop", "Ref", "Alt")

write.table(annovar, file = "~/Downloads/AnswerALStest.avinput",
            append = FALSE, quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE)

