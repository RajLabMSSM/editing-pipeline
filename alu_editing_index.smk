# Part of RNA Editing Pipeline
# Jack Humphrey, Raj Lab
# 2021



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
refDir = config["refDir"]


metadata =  pd.read_csv(config["metadata"], sep = "\t")

SAMPLES = metadata["sample"]

metadata_dict = metadata.set_index('sample').T.to_dict()

print(metadata_dict)

rule all:
    input:
        expand( outDir + "{SAMPLE}/EditingIndex.csv", SAMPLE = SAMPLES)
        #expand(dataCode + "_{group}_merged_sites.annotated.filtered.txt", group = groups),
        #expand( "{sample}/{sample}.config.json", sample = samples),
        #expand( "{sample}/{sample}.sites.snp_filtered.bed", sample = samples)


rule AEI:
    input:
        metadata_dict[wildcards.SAMPLE]["bam_path"] + "/{SAMPLE}.bam"
    output:
        outDir + "{SAMPLE}/EditingIndex.csv"
    params:
        bam_dir = metadata_dict["{SAMPLE}"]["bam_path"],
        bam_suffix = ".bam",
        out_dir = outDir + "/{SAMPLE}",
        genes_expression = refDir + "/GenesExpression/HomoSapiens/ucscHg38GTExGeneExpression.bed.gz",
        refseq = refDir + "/RefSeqAnnotations/HomoSapiens/ucscHg38RefSeqCurated.bed.gz",
        snps = refDir + "/SNPs/HomoSapiens/ucscHg38CommonGenomicSNPs150.bed.gz",
        gf = refDir + "/Genomes/HomoSapiens/ucscHg38Genome.fa",
        rb = refDir + "/Regions/HomoSapiens/ucscHg38Alu.bed.gz"
    shell:
        "mkdir -p {params.out_dir}; "
        "cd {params.out_dir}; "
        "ml rnaeditingindexer; "
        "ml samtools/1.9; "
        "rm -r {params.out_dir}/flags; "
        "RNAEditingIndex "
        "-d $PWD/{params.bam_dir} "
        " -f {params.bam_suffix} "
        " -o $PWD/{params.out_dir} "
        " --log_path $PWD/{outDir} "
        "--genes_expression {params.genes_expression} "
        "--refseq {params.refseq} "
        "--snps {params.snps} "
        " -gf {params.gf} "
        " -rb {params.rb} "
        "--genome UserProvided " 
        "--paired_end " 
        "--stranded "
        "--verbose "

# this may have to change depending on sample naming
wildcard_constraints:
    SAMPLE = '[A-Za-z0-9_-]+'

