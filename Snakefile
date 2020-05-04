# RNA Editing Pipeline
# Jack Humphrey, Raj Lab
# 2020

# wraps around the SAILOR pipeline, maintained by Brian Yee (https://github.com/YeoLab/sailor)

sailor = "/sc/hydra/projects/ad-omics/data/software/sailor/sailor-1.2.0.img"

import pandas as pd
import os
# inFolder must be full path
import yaml
import json

metadata =  pd.read_csv(config["samples"], sep = "\t")

samples = metadata['sample']

metadata_dict = metadata.set_index('sample').T.to_dict()

inFolder = config["inFolder"]

# hardcoded for now
refDir = '/sc/hydra/projects/PBG/REFERENCES/GRCh38/FASTA/'
dbSNPDir = '/sc/hydra/projects/ad-omics/data/references/hg38_reference/dbSNP/'

# create nested dictionary for storing reference files
#genomePath = "/sc/hydra/projects/ad-omics/data/references/hg38_reference/hg38.fa"
#knownSNPPath = "/sc/hydra/projects/ad-omics/data/references/hg38_reference/dbSNP/dbSNP_GCF_000001405.38.bed"

#config = {'input_bam': {'class': 'File', 'path': 'NA'}, 'reference': {'class': 'File', 'path': genomePath}, 'known_snp': {'class': 'File', 'path': knownSNPPath} }


rule all:
    input:
        expand( "{sample}/{sample}.config.json", sample = samples),
        expand( "{sample}/{sample}.fwd.sorted.rmdup.readfiltered.formatted.varfiltered.snpfiltered.ranked.conf", sample = samples)

rule writeConfig:
    output:
        config = "{sample}.config.json"
    run:
        # subset metadata - get the corresponding BAM file
        bam = metadata_dict[wildcards.sample]['bam']
        # create data structure and save as JSON
        config['input_bam']['path'] = bam
        stream = open(output.config, 'w')
        json.dump(config, stream)

# this may have to change depending on sample naming
wildcard_constraints:
    sample = '[A-Za-z0-9_-]+'

rule SAILOR:
    input:
        config = "{sample}.config.json",
    output:
        "{sample}/{sample}.config.json",
        "{sample}/{sample}.fwd.sorted.rmdup.readfiltered.formatted.varfiltered.snpfiltered.ranked.conf"
    run:
        bam = metadata_dict[wildcards.sample]['bam']
        bamDir = os.path.dirname(bam)
        shell("ml singularity/3.5.2; singularity run --bind {bamDir} --bind {refDir} --bind {dbSNPDir} {sailor} {input.config}")
        shell("mkdir -p {wildcards.sample}")
        shell("mv {wildcards.sample}.* {wildcards.sample}/")
