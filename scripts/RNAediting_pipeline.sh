#!/bin/bash
#BSUB -J "MULTI"
#BSUB -P acc_buxbaj01a
#BSUB -q premium
#BSUB -o ./step1_%J
#BSUB -u enrico.mossotto@mssm.edu
#BSUB -L /bin/bash
#BSUB -W 10:00
#BSUB -n 40
#BSUB -R rusage[mem=4000]
#BSUB -R span[hosts=1]


. variables
INFILE=$BAM
OUTFILE=${sampleID}.rediout.gz
REF=/sc/arion/projects/PBG/REFERENCES/GRCh38/FASTA/GRCh38.chrom.fa
REFINDEX=/sc/arion/projects/PBG/REFERENCES/GRCh38/FASTA/GRCh38.chrom.fa.fai
PROC=40
HOMPOL=/sc/hydra/projects/buxbaj01a/enrico/homopolymeric_sites_hg38.txt
#minimum number of supporting reads
CHANGE_DP=3
#minimum depth at locus
LOCUS_DP=5

#### Sort and index BAM file

module load samtools

samtools sort -m 3G -o ${sampleID}_sorted.bam -@ $PROC $INFILE
samtools index ${sampleID}_sorted.bam -@ $PROC

module purge

module load reditools
### create temp directory
mkdir temp_cov
mkdir temp_results

######## REDITOOLS ################
## generate coverage file for multi-threading (~1h)

/hpc/packages/minerva-centos7/reditools/2.0/reditools2.0/extract_coverage_dynamic.sh ${sampleID}_sorted.bam temp_cov/ $REFINDEX

# Run parallel reditools
mpirun -np 40 /hpc/packages/minerva-centos7/reditools/2.0/reditools2.0/src/cineca/parallel_reditools.py -f ${sampleID}_sorted.bam -r $REF -S -s 2 -ss 5 -mrl 50 -q 10 -bq 20 -C -T 2 -m $HOMPOL -os 5 -t temp_results/ -Z $REFINDEX -G temp_cov/${sampleID}*.cov -D temp_cov/ -o ${sampleID}_out

# Merge outputs in a single file
module load htslib/1.7

/hpc/packages/minerva-centos7/reditools/2.0/reditools2.0/merge.sh temp_results $OUTFILE 40

bgzip -d  ${sampleID}.rediout.gz

#remove temp folders and files
rm -fr temp_cov
rm -fr temp_results
rm ${sampleID}_sorted.bam ${sampleID}_sorted.bam.bai

############ FILTERING and ANNOTATION#########
module unload reditools
module load annovar


#set a minimum depth for the locus of 5
awk 'BEGIN{OFS=FS="\t"}{if($5>5) print $0}' ${sampleID}.rediout  >${sampleID}_final_sorted_cov5.redi

#split multi editing events in multiple lines
awk 'BEGIN{OFS=FS="\t"} {if(length($8)>2) {split($8,a," "); for (i in a){$8=a[i];print $0}} else print $0}' ${sampleID}_final_sorted_cov5.redi > ${sampleID}_final_sorted_cov5_linesplit.redi

rm ${sampleID}_final_sorted_cov5.redi

#Format for annovar and reverse ref/alt for reverse transcripts
awk 'BEGIN{OFS=FS="\t"} {if(NR==1) true; else if (NR>1){ split($8,a,""); if ($4==0) {if ($3=="A") {$3="T"} else if ($3=="C") {$3="G"} else if ($3=="G") {$3="C"} else if ($3=="T") {$3="A"}; if (a[2]=="A") {a[2]="T"} else if (a[2]=="C") {a[2]="G"} else if (a[2]=="G") {a[2]="C"} else if (a[2]=="T") {a[2]="A"}}; print $1,$2,$2,$3,a[2],$4,$5,$6,$7,$8}}' ${sampleID}_final_sorted_cov5_linesplit.redi > ${sampleID}_annovar.input


rm ${sampleID}_final_sorted_cov5_linesplit.redi 

### Filter low quality variants
module load python

python /sc/arion/projects/buxbaj01a/enrico/PIPELINE/filter_DP.py ${sampleID}_annovar.input $CHANGE_DP $LOCUS_DP > ${sampleID}_annovar_filtered.input

### Filter out known genomic position (if vcf available)
#format vcf into annovar db
zcat $VCF >tmp1
awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$2,$4,$5}' tmp1 >genomic_hits
rm tmp1
#annotate variants by genomics hits
/sc/arion/projects/buxbaj01a/resources/annovar/annotate_variation.pl --filter --dbtype generic --genericdbfile genomic_hits -build hg38 ${sampleID}_annovar_filtered.input .

# mantain only variants not matching genomic data
mv ${sampleID}_annovar_filtered.input.hg38_generic_filtered ${sampleID}_annovar_filtered_nogen.input 

#clean up
rm genomic_hits ${sampleID}_annovar_filtered.input ${sampleID}_annovar_filtered.input.invalid_input ${sampleID}_annovar_filtered.input.hg38_generic_dropped ${sampleID}_annovar_filtered.input.log


###Run annovar

/sc/arion/projects/buxbaj01a/resources/annovar/table_annovar.pl ${sampleID}_annovar_filtered_nogen.input /sc/arion/projects/buxbaj01a/resources/annovar/humandb/ -buildver hg38 -out ${sampleID} -remove -protocol refGene,dbsnp153CommonSNV,gnomad30_genome,phastConsElements30way,rmsk,rediportal_012920 -operation g,f,f,r,r,f --argument ,,,'--colsWanted 5','--colsWanted 10&11&12', -nastring "." --otherinfo --thread $PROC --maxgenethread $PROC

gzip -f ${sampleID}.hg38_multianno.txt

