import os
import pandas as pd

# reference data
reditools = "/hpc/packages/minerva-centos7/reditools/2.0/reditools2.0/"
KNOWN_SNPS_BED = "/sc/arion/projects/breen_lab/AEI/SNPs/HomoSapiens/ucscHg38CommonGenomicSNPs150.bed.gz"

# make a test dataset with config.yaml and samples.tsv for only 1 sample
# extract chr21 from a bam file using
# samtools view -bh {sample} chr21 > test.bam
# samtools index test.bam



# hard-coded input and output - needs to be put into separate config.yaml
INFOLDER =  "/sc/arion/projects/als-omics/microglia_stimulated/editing/reditools/test/"
OUTFOLDER = "/sc/arion/projects/als-omics/microglia_stimulated/editing/reditools/"
OUTFOLDER = "reditools/"
# needs separate samples.tsv
SAMPLES = ["test"]

# cohort or sample specific variables
# need to be in separate config or samples.TSV
# needs to ideally be in sample metadata when you have mixed stranding between batches

# have a strand column in samples.tsv so then the pipeline can adjust the strand parameter automatically
STRAND = 0
CHANGE_DP = 3 #minimum number of supporting reads
LOCUS_DP=5 #minimum depth at a locus - in paper they use 10. not sure what to do

# homopolymeric sites in the hg38 genome build are recorded by reditools when you run the -c option. This output is now saved for use.
rule all:
    input:
        expand(OUTFOLDER + "{SAMPLE}/{SAMPLE}.reditools.filtered.gz", SAMPLE = SAMPLES)

rule reditools:
    input: INFOLDER + "{SAMPLE}.bam"
    output: OUTFOLDER + "{SAMPLE}/reditools/{SAMPLE}/{SAMPLE}.rediout.gz"
    params:
        REF = "/sc/arion/projects/ad-omics/data/references//hg38_reference/hg38.fa",
        REFINDEX = "/sc/arion/projects/ad-omics/data/references/hg38_reference/hg38.standard_chr.fa.fai",
        #REFINDEX = "/sc/arion/projects/breen_lab/AEI/Genomes/HomoSapiens/ucscHg38Genome.fa.fai",
        threads = 4,
        HOMPOL = "/sc/arion/projects/ad-omics/data/references/editing/homopolymeric_sites_hg38.txt"
        
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
# rewrite AWK filters in python
# add to filter_DP to create just one big python script for filtering
# bedtools - understand what this is doing - set intersections
rule filtering:
    input:
        OUTFOLDER + "{SAMPLE}/reditools/{SAMPLE}/{SAMPLE}.rediout.gz"
        #OUTFOLDER + "{SAMPLE}/{SAMPLE}.rediout.gz" 
    output:
        OUTFOLDER + "{SAMPLE}/{SAMPLE}.reditools.filtered.gz"
    params:
        prefix = OUTFOLDER + "{SAMPLE}/{SAMPLE}",
        filter_awk = "scripts/filter_reditools.awk",
        filter_python = "scripts/filter_DP.py",
        dbsnp = "/sc/arion/projects/ad-omics/data/references//hg38_reference/dbSNP/hg38_common_snps_dbSNP_153.bed.gz"
    shell:
        "module load annovar;"
        "sh {params.filter_awk} {input} {params.prefix}_annovar.input; "
        ### Filter low quality variants
        "module load python;"
        "python {params.filter_python} {params.prefix}_annovar.input {CHANGE_DP} {LOCUS_DP} > {params.prefix}_annovar_filtered.input;"
        "ml bedtools;"
        "bedtools intersect -v -a {params.prefix}_annovar_filtered.input} -b {params.dbsnp} > {output}"
        #clean up
        # make sure to not delete the important parts - just the trash

        #"rm genomic_hits {params.prefix}_annovar_filtered.input "
        #" {params.prefix}_annovar_filtered.input.invalid_input {params.prefix}_annovar_filtered.input.hg38_generic_dropped "
        #"{params.prefix}_annovar_filtered.input.log"

###Run annovar
# need to add inputs and outputs , and change rule all to accept the output of this step
#rule annovar:
##        "table_annovar.pl {params.prefix}_annovar_filtered_nogen.input /sc/arion/projects/buxbaj01a/resources/annovar/humandb/ -buildver hg38 -out {params.prefix} "
##        " -remove -protocol refGene,dbsnp153CommonSNV,gnomad30_genome,phastConsElements30way,rmsk,rediportal_012920 -operation g,f,f,r,r,f "\
##        " --argument ,,,'--colsWanted 5','--colsWanted 10&11&12', -nastring \".\" --otherinfo --thread {params.threads} --maxgenethread {params.threads};"
##        # gzip
##        "gzip -f {params.prefix}.hg38_multianno.txt"a


## new rule - filter on annovar

# read Breen/Cuddleston plus other paper's methods to see which filters to apply

# another new rule - bring together all samples in cohorta

# penultimate rule - perform filtering across samples - look at Cuddleston paper or others for examples

# final rule - deal with missing data. go back and re-call those specific sites using samtools mpileup on the original BAM files - look at Gabay 2022 paper for inspiration
