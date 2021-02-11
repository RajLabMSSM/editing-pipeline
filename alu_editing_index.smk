# Part of RNA Editing Pipeline
# Jack Humphrey, Raj Lab
# 2021
refDir: "/sc/arion/projects/breen_lab/AEI/"


import pandas as pd
import os
# inFolder must be full path
#import yaml
##import json

#inFolder = config["inFolder"]
#genes_bed = config["genes_bed"]

#dataCode = config['dataCode']
#groups = config['groups']

# hardcoded for now
#refDir = '/sc/hydra/projects/PBG/REFERENCES/GRCh38/FASTA/'
#dbSNPDir = '/sc/hydra/projects/ad-omics/data/references/hg38_reference/dbSNP/'

# create nested dictionary for storing reference files
#genomePath = "/sc/hydra/projects/ad-omics/data/references/hg38_reference/hg38.fa"
#knownSNPPath = "/sc/hydra/projects/ad-omics/data/references/hg38_reference/dbSNP/dbSNP_GCF_000001405.38.bed"

#config = {'input_bam': {'class': 'File', 'path': 'NA'}, 'reference': {'class': 'File', 'path': genomePath}, 'known_snp': {'class': 'File', 'path': knownSNPPath} }

outDir = config["outDir"]

metadata =  pd.read_csv(config["metadata"], sep = "\t")

SAMPLES = metadata["sample"]
BAMS = metadata["bam_path"] + metadata["sample"] + ".bam"
metadata_dict = metadata.set_index('sample').T.to_dict()

paired = config["paired"]
stranded = config["stranded"]

rule all:
    input:
        expand( outDir + "{SAMPLE}/EditingIndex.csv", SAMPLE = SAMPLES)
        #expand(dataCode + "_{group}_merged_sites.annotated.filtered.txt", group = groups),
        #expand( "{sample}/{sample}.config.json", sample = samples),
        #expand( "{sample}/{sample}.sites.snp_filtered.bed", sample = samples)


rule AEI:
    input:
        #metadata_dict["{SAMPLE}"]["bam_path"] + "/{SAMPLE}.bam"
    output:
        outDir + "{SAMPLE}/EditingIndex.csv"
    params:
        #bam_dir = metadata_dict["{SAMPLE}"]["bam_path"],
        bam_suffix = ".bam",
        out_dir = outDir + "/{SAMPLE}",
        genes_expression = refDir + "/GenesExpression/HomoSapiens/ucscHg38GTExGeneExpression.bed.gz",
        refseq = refDir + "/RefSeqAnnotations/HomoSapiens/ucscHg38RefSeqCurated.bed.gz",
        snps = refDir + "/SNPs/HomoSapiens/ucscHg38CommonGenomicSNPs150.bed.gz",
        gf = refDir + "/Genomes/HomoSapiens/ucscHg38Genome.fa",
        rb = refDir + "/Regions/HomoSapiens/ucscHg38Alu.bed.gz"
    run:
        bam_dir = metadata_dict[wildcards.SAMPLE]["bam_path"]
        
        if stranded == True:
            strand_cmd = " --stranded"
        else:
            strand_cmd = ""  

        if paired == True: 
            paired_cmd = " --paired"
        else:
            paired_cmd = ""
      
        shell( "mkdir -p {params.out_dir}; \
                cd {params.out_dir}; \
                ml rnaeditingindexer; \
                ml samtools/1.9; \
                rm -rf {params.out_dir}/flags; \
                RNAEditingIndex \
                -d {bam_dir} \
                -f {params.bam_suffix} \
                -o {params.out_dir} \
                --log_path {params.out_dir} \
                --genes_expression {params.genes_expression} \
                --refseq {params.refseq} \
                --snps {params.snps} \
                -gf {params.gf} \
                -rb {params.rb} \
                --genome UserProvided \
                {paired_cmd} \
                {strand_cmd} "
        )

# this may have to change depending on sample naming
wildcard_constraints:
    SAMPLE = '[A-Za-z0-9_-]+'

