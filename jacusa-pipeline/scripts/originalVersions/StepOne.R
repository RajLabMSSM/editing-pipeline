setwd("/sc/arion/projects/ad-omics/winston/answerALS/")

file_list=read.table("testsamples.txt", stringsAsFactors=FALSE)
files=basename(file_list$V1)
full_run_input <- NULL
full_run_input <- paste('#!/bin/bash\n', sep="")
for(item in 1:length(files)){
        file_name0=(files[item - 1])
        file_name1=(files[item])
        sh_file <- paste(file_name1,"_jacusa2.sh", sep="")
        full_run_input <- paste(full_run_input,'bsub < ', sh_file, ';\n sleep 1;\n', sep="")
        sh_runfile <- "run_all_jacusa2.sh"
        fileConn <- file(sh_file)

writeLines(paste('
#!/bin/bash
#BSUB -W 10:00
#BSUB -n 5
#BSUB -q express
#BSUB -P acc_ad-omics
#BSUB -J ',file_name1,'
#BSUB -cwd /sc/arion/projects/ad-omics/winston/answerALS
#BSUB -o /sc/arion/projects/ad-omics/winston/answerALS/jacusa2logfiles/',file_name1,'.log.out
#BSUB -e /sc/arion/projects/ad-omics/winston/answerALS/jacusa2logfiles/',file_name1,'.log.err
#BSUB -u winston.cuddleston@icahn.mssm.edu
#BSUB -R rusage[mem=14000]
#BSUB -R span[hosts=1]
#BSUB -L /bin/bash

ml jacusa2
ml bedtools

java -jar $JACUSA2_JAR call-1 -r /sc/arion/projects/ad-omics/winston/answerALS/',file_name1,'.out /sc/arion/projects/bigbrain/data/Answer_ALS/bams_from_RAPiD/',file_name1,'.bam -p 10 -a D,M,Y,E:file=/sc/arion/projects/ad-omics/winston/jacusa_pipeline/hg38-blacklist.v2_sort.bed:type=BED -s -m 20 -R /sc/arion/projects/ad-omics/data/references/GRCh38_references/GRCh38.primary_assembly.genome.fa -P FR-SECONDSTRAND -F 1024

sed -i '1d' ',file_name1,'.out

', sep=""), fileConn)
close(fileConn)
}

fileConn <- file(sh_runfile)
writeLines(paste(full_run_input, sep=""), fileConn)
close(fileConn)


