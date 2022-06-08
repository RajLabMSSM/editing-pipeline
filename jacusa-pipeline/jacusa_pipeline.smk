#RNA Editing Pipline
#Winston Cuddleston, Raj Lab
#2022

import pandas as pd
import os

#hardcoded for now
refDir = '/sc/arion/projects/ad-omics/data/references/GRCh38_references/'
editingRefDir = '/sc/arion/projects/ad-omics/data/references/editing/'

rule all:

rule Jacusa:
	input:
		bam_dir = metadata_dict[wildcards.SAMPLE]["bam_path"]
		full_bam = bam_dir + "/" + wildcards.SAMPLE + ".bam"
		lib = metadata_dict[wildcards.SAMPLE]["library"]
		blacklist = editingRefDir + "hg38-blacklist.v2_sort.bed"
		fastaref = refDir + "GRCh38.primary_assembly.genome.fa"
		
	output:
		"/sc/arion/projects/ad-omics/winston/editing_with_snakemake/{sample}.out"
	run:
		shell("ml jacusa2; \
			   java -jar $JACUSA2_JAR call-1 -r {output} {input.full_bam} \
			   -p 10 -a D,M,Y,E:file={input.blacklist}:type=BED -s -m 20 \
			   -R {input.fastaref} -P {input.lib} -F 1024"
			   )

rule firstFiltering:
	input:
		"/sc/arion/projects/ad-omics/winston/editing_with_snakemake/{sample}.out"
	output:
		"/sc/arion/projects/ad-omics/winston/editing_with_snakemake/{sample}.filtered.out"
	run:
		shell("sed -i '1d' {sample}.out")
		shell("ml R/3.6.0;")
		shell("Rscript ")
