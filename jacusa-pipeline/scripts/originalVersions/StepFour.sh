#BSUB -W 10:00
#BSUB -n 5
#BSUB -q express
#BSUB -P acc_ad-omics
#BSUB -J annoALS
#BSUB -cwd /sc/arion/projects/ad-omics/winston/jacusa_pipeline/annovar
#BSUB -o /sc/arion/projects/ad-omics/winston/jacusa_pipeline/annovar/annovar.log.out
#BSUB -e /sc/arion/projects/ad-omics/winston/jacusa_pipeline/annovar/annovar.log.err
#BSUB -u winston.cuddleston@icahn.mssm.edu
#BSUB -R rusage[mem=14000]
#BSUB -R span[hosts=1]
#BSUB -L /bin/bash

ml annovar
ml bcftools

IN="/sc/arion/projects/ad-omics/winston/jacusa_pipeline/annovar/AnswerALStest.avinput"
OUT="/sc/arion/projects/ad-omics/winston/jacusa_pipeline/annovar/AnswerALS.myanno"

table_annovar.pl $IN /sc/arion/projects/ad-omics/winston/jacusa_pipeline/humandb/ -buildver hg38 -out $OUT -remove -protocol refGene,dbsnp153CommonSNV,gnomad30_genome,phastConsElements30way,rmsk,rediportal_012920 -operation g,f,f,r,r,f --argument ,,,'--colsWanted 5','--colsWanted 10&11&12', -nastring "." --otherinfo --thread 10 --maxGeneThread 10

awk 'BEGIN{OFS=FS="\t"}{if ( ($11=="." || $11=="dbSNP153CommonSNV") && ( $12<0.05 || $12=="AF")) print $0}' AnswerALS.myanno.hg38_multianno.txt > AnswerALS.myanno.hg38_multianno.txt.noCommon.txt
