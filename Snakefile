
sailor = "/sc/hydra/projects/ad-omics/data/software/sailor/sailor-1.2.0.img"

import pandas as pd
import os
# inFolder must be full path
import yaml
import json

metadata =  "test_samples.tsv"

samples = pd.read_csv(metadata, sep = "\t")['sample']

inFolder = "/hpc/users/humphj04/pipelines/editing-pipeline/test/"

# create nested dictionary for storing reference files
#genomePath = "/sc/hydra/projects/ad-omics/data/references/hg38_reference/hg38.fa"
#knownSNPPath = "/sc/hydra/projects/ad-omics/data/references/hg38_reference/dbSNP/dbSNP_GCF_000001405.38.bed"

#config = {'input_bam': {'class': 'File', 'path': 'NA'}, 'reference': {'class': 'File', 'path': genomePath}, 'known_snp': {'class': 'File', 'path': knownSNPPath} }


rule all:
    input:
        expand( "{sample}/{sample}.config.yaml", sample = samples),
        expand( "{sample}/{sample}.fwd.sorted.rmdup.readfiltered.formatted.varfiltered.snpfiltered.ranked.conf", sample = samples)

rule writeConfig:
    input:
        #genome = genomePath,
        #snps = knownSNPPath,
        bam = inFolder + "{sample}.bam",
        bai = inFolder + "{sample}.bam.bai"
    output:
        #bam = "{sample}/input/{sample}.bam",
        #bai = "{sample}/input/{sample}.bam.bai",
        #genome = refFolder + "hg38.fa",
        #snps = refFolder + "dbSNP.bed",
        config = "{sample}.config.yaml"
    run:
        #os.symlink(input.bam, output.bam)
        #os.symlink(input.bai, output.bai)
        #os.symlink(input.genome, output.genome)
        #os.symlink(input.snps, output.snps)
        # create data structure and save as YAML
        config['input_bam']['path'] = input.bam
        #config['reference']['path'] = output.genome
        #config['known_snp']['path'] = output.snps
        stream = open(output.config, 'w')
        json.dump(config, stream)

wildcard_constraints:
    sample = '[A-Za-z0-9_]+'

rule SAILOR:
    input:
        config = "{sample}.config.yaml",
        genome = config['reference']['path'],
        snps = config['known_snp']['path'],
        bam = inFolder + "{sample}.bam",
        bai = inFolder + "{sample}.bam.bai"
    output:
        "{sample}/{sample}.config.yaml",
        "{sample}/{sample}.fwd.sorted.rmdup.readfiltered.formatted.varfiltered.snpfiltered.ranked.conf"
    params: 
        config = "{sample}.config.yaml"
    shell:
        "ml singularity/3.5.2;"
        #"cd {wildcards.sample};"
        "singularity run {sailor} {input.config}; "
        "mkdir {wildcards.sample};"
        "mv {wildcards.sample}.* {wildcards.sample}/"
        #"--alpha 0 --beta 0 --ct --dp DP4 "
        #"--edge_mutation 5 --edit_fraction 0.01 " 
        #" --fwd_is_reverse --input_bam {input.bam} " 
        #" --junction_overhang 10 --keep_all_edited --known_snp {input.snps} "
        #" --min_variant_coverage 5 "
        #" --non_ag 1 --reference {input.genome} --rev_is_reverse"
        #" [--reverse_stranded_library] [--single_end] --skip_duplicate_removal "
