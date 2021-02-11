import os
import pandas as pd

reditools = "/hpc/packages/minerva-centos7/reditools/2.0/reditools2.0/"
KNOWN_SNPS_BED = "/sc/arion/projects/breen_lab/AEI/SNPs/HomoSapiens/ucscHg38CommonGenomicSNPs150.bed.gz"
INFOLDER =  "/sc/arion/projects/als-omics/microglia_stimulated/editing/reditools/test/"
OUTFOLDER = "/sc/arion/projects/als-omics/microglia_stimulated/editing/reditools/"
SAMPLES = ["test"]

# needs to ideally be in sample metadata when you have mixed stranding between batches
STRAND = 0

CHANGE_DP = 3        #minimum number of supporting reads
LOCUS_DP=5 #minimum depth at a locus
# homopolymeric sites in the hg38 genome build are recorded by reditools when you run the -c option. This output is now saved for use.
rule all:
    input:
        expand(OUTFOLDER + "{SAMPLE}/{SAMPLE}.reditools.filtered.gz", SAMPLE = SAMPLES)

rule reditools:
    input: INFOLDER + "{SAMPLE}.bam"
    output: OUTFOLDER + "{SAMPLE}/{SAMPLE}.rediout.gz"
    params:
        REF = "/sc/arion/projects/ad-omics/data/references//hg38_reference/hg38.fa",
        REFINDEX = "/sc/arion/projects/ad-omics/data/references//hg38_reference/hg38.standard_chr.fa.fai",
        #REFINDEX = "/sc/arion/projects/breen_lab/AEI/Genomes/HomoSapiens/ucscHg38Genome.fa.fai",
        threads = 4,
        HOMPOL = "/sc/arion/projects/ad-omics/data/references/editing/homopolymeric_sites_hg38.txt",
        
    shell:
        "module load reditools;"
        "mkdir -p temp_cov;"
        "mkdir -p temp_results;"
        "cd {OUTFOLDER}/{wildcards.SAMPLE}/;"

        ######## REDITOOLS ################
        ## generate coverage file for multi-threading (~1h)

        "{reditools}/extract_coverage_dynamic.sh {input} temp_cov/ {params.REFINDEX};"

        # Run parallel reditools
        "mpirun -np {params.threads} {reditools}/src/cineca/parallel_reditools.py -f {input} "
        " -r {params.REF} -S -s {STRAND} -ss 5 -mrl 50 -q 10 -bq 20 -C -T 2 -c -m {params.HOMPOL} "
        " -os 5 -t temp_results/  -Z {params.REFINDEX} -G temp_cov/{wildcards.SAMPLE}*.cov "
        " -D temp_cov/ -o {wildcards.SAMPLE}_out; "

        # Merge outputs in a single file
        "module load htslib/1.7; "

        "{reditools}/merge.sh temp_results {output} {params.threads}; "

        #"bgzip -d  {wildcards.SAMPLE}.rediout.gz;"

        #remove temp folders and files
        #"rm -fr temp_cov;"
        #"rm -fr temp_results"


############ FILTERING and ANNOTATION#########
rule filtering:
    input:
        OUTFOLDER + "{SAMPLE}/{SAMPLE}.rediout.gz" 
    output:
        OUTFOLDER + "{SAMPLE}/{SAMPLE}.reditools.filtered.gz"
    params:
        prefix = OUTFOLDER + "{SAMPLE}/{SAMPLE}",
        filter_awk = "scripts/filter_reditools.awk",
        filter_python = "scripts/filter_DP.py"
        dbsnp = "/sc/arion/projects/ad-omics/data/references//hg38_reference/dbSNP/hg38_common_snps_dbSNP_153.bed.gz"
    shell:
        "module load annovar;"
        "sh {params.filter_awk} {input} {params.prefix}_annovar.input; "
        ### Filter low quality variants
        "module load python;"
        "python {params.filter_python} {params.prefix}_annovar.input {CHANGE_DP} {LOCUS_DP} > {params.prefix}_annovar_filtered.input;"
        "ml bedtools;"
        "bedtools intersect -v -a {params.prefix}_annovar_filtered.input} -b {params.dbsnp} > {output}"

        ### Filter out known genomic position (if vcf available)
        #format vcf into annovar db
        #zcat $VCF >tmp1
        #awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$2,$4,$5}' tmp1 >genomic_hits
        #rm tmp1
        #annotate variants by genomics hits
        # filter out common SNPs?
        #"annotate_variation.pl --filter --dbtype bed --bedfile {KNOWN_SNPS_BED} -build hg38 {params.prefix}_annovar_filtered.input . ;"

        # mantain only variants not matching genomic data
        #"mv {params.prefix}_annovar_filtered.input.hg38_generic_filtered {output}" #{params.prefix}_annovar_filtered_nogen.input;" 

        #clean up
        "rm genomic_hits {params.prefix}_annovar_filtered.input {params.prefix}_annovar_filtered.input.invalid_input {params.prefix}_annovar_filtered.input.hg38_generic_dropped {params.prefix}_annovar_filtered.input.log"


        ###Run annovar
rule annovar:
        "table_annovar.pl {params.prefix}_annovar_filtered_nogen.input /sc/arion/projects/buxbaj01a/resources/annovar/humandb/ -buildver hg38 -out {params.prefix} "
        " -remove -protocol refGene,dbsnp153CommonSNV,gnomad30_genome,phastConsElements30way,rmsk,rediportal_012920 -operation g,f,f,r,r,f "\
        " --argument ,,,'--colsWanted 5','--colsWanted 10&11&12', -nastring \".\" --otherinfo --thread {params.threads} --maxgenethread {params.threads};"
        # gzip
        "gzip -f {params.prefix}.hg38_multianno.txt"
