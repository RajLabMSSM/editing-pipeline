#RNA Editing Pipline
#Winston Cuddleston, Raj Lab
#2022

import pandas as pd
import os

#references
refDir = config["refDir"]
editingRefDir = config["editingRefDir"]
humandbDir = config["humandbDir"]

#data
projectDir = config["projectDir"]
metadata =  pd.read_csv(config["sample"], sep = "\t")
samples = metadata['sample']
metadata_dict = metadata.set_index('sample').T.to_dict()

rule all:
	input:
		projectDir + "ESannotations.txt",
		projectDir + "filtCoverageMatrix.txt",
		projectDir + "filtRatioMatrix.txt"

rule jacusa:
	input:
		bam_dir = metadata_dict[wildcards.SAMPLE]["bam_path"],
		full_bam = bam_dir + "/" + {wildcards.sample} + ".bam",
		lib = metadata_dict[wildcards.SAMPLE]["library"],
		blacklist = editingRefDir + "hg38-blacklist.v2_sort.bed",
		fastaref = refDir + "GRCh38.primary_assembly.genome.fa"
		
	output:
		outFile = projectDir + {wildcards.sample} + ".out"
	shell:
		"ml jacusa2;"
		"java -jar $JACUSA2_JAR call-1"
		"-r {output.outFile} {input.full_bam}"
		"-p 10 -a D,M,Y,E:file={input.blacklist}:type=BED"
		"-s -m 20 -R {input.fastaref} -P {input.lib} -F 1024"

rule firstFiltering:
	input:
		inFile = projectDir + {wildcards.sample} + ".out"
	output:
		outFile = projectDir + {wildcards.sample} + ".filt"
	params:
		script = "scripts/firstfiltering.R"
		
	shell:
		"sed -i '1d' {input.inFile};"
		"ml R/4.0.2;"
		"Rscript {params.script}"
		"--input {input.inFile}"
		"--output {output.outFile}"
		
rule mergeSecondFiltering:
	input:
		expand(projectDir + "{sample}.filt", sample = samples)
	output:
		covMat = projectDir + "coverageMatrix.txt",
		ratioMat = projectDir + "ratioMatrix.txt",
		av = projectDir + "avinput.txt"
	params:
		script = "scripts/mergeSecondFiltering.R",
		inDir = projectDir,
		perc = 0.5,
		er = 0.1
	shell:
		"ml R/4.0.2;"
		"Rscript {params.script}"
		"--inDir {params.inDir}"
		"--rat {output.ratioMat}"
		"--cov {output.covMat}"
		"--av {output.av}"
		"--percSamples {params.perc}"
		"--minER {params.er}"

rule annovar:
	input:
		av = projectDir + "avinput.txt"
	output:
		anno = projectDir +  "myanno.hg38_multianno.txt"
	params:
		outfile = projectDir + "myanno",
		humandb = humandbDir
	shell:
		"ml annovar;"
		"table_annovar.pl {input.avIN} {params.humandb} -buildver hg38 \
		-out {params.outfile} -remove -protocol refGene,dbsnp153CommonSNV,\
		gnomad30_genome,phastConsElements30way,rmsk,rediportal_012920\
		-operation g,f,f,r,r,f --argument ,,,'--colsWanted 5','--colsWanted 10&11&12',\
		-nastring "." --otherinfo --thread 10 --maxGeneThread 10"
		
rule dropSNPs:
	input:
		anno = projectDir + "myanno.hg38_multianno.txt"
	output:
		filtanno = projectDir + "myanno.hg38_multianno.txt.noCommon.txt"
	shell:	
		"awk 'BEGIN{OFS=FS="\t"}{if ( ($11=="." || $11=="dbSNP153CommonSNV") \
		&& ( $12<0.05 || $12=="AF")) print $0}' {input.anno} > {output.filtanno}"

rule annovarFiltering:
	input:
		filtanno = projectDir + "myanno.hg38_multianno.txt.noCommon.txt",
		covMat = projectDir + "coverageMatrix.txt",
		ratioMat = projectDir + "ratioMatrix.txt"
	output:
		ESanno = projectDir + "ESannotations.txt",
		filtCovMat = projectDir + "filtCoverageMatrix.txt",
		filtRatioMat = projectDir + "filtRatioMatrix.txt"
	params:
		script = "scripts/annovarfiltering.R"
	shell:
		"ml R/4.0.2;"
		"Rscript {params.script}"
		"--inAnno {input.filtanno}"
		"--inRat {input.ratioMat}"
		"--inCov {input.covMat}"
		"--outAnno {output.ESanno}"
		"--outRat {output.filtRatMat}"
		"--outCov {output.filtCovMat}"
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		