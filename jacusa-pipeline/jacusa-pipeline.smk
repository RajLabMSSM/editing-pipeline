#RNA Editing Pipline
#Winston Cuddleston, Flora Seo, Jack Humphrey, Raj Lab
#2022

R_VERSION = "R/4.0.3"
jacusa_threads = 10
annovar_threads = 12
import pandas as pd
import os

#references
refDir = config["refDir"]
editingRefDir = config["editingRefDir"]
humandbDir = config["humandbDir"]

#data
projectDir = config["projectDir"]
metadata =  pd.read_csv(config["metadata"], sep = "\t")
samples = metadata['sample']
metadata_dict = metadata.set_index('sample').T.to_dict()

rule all:
	input:
		projectDir + "ESannotations.txt",
		projectDir + "filtCoverageMatrix.txt",
		projectDir + "filtRatioMatrix.txt"

rule jacusa:
	input:
		blacklist = editingRefDir + "hg38-blacklist.v2_sort.bed",
		fastaref = refDir + "GRCh38.primary_assembly.genome.fa"
	output:
		outFile = projectDir + "{sample}/{sample}.out"
	run:
		bam_dir = metadata_dict[wildcards.sample]["bam_path"]
		full_bam = bam_dir + "/" + wildcards.sample + ".bam"
		lib = metadata_dict[wildcards.sample]["library"]
		shell("ml jacusa2; \
		       java -jar $JACUSA2_JAR call-1 \
		       -r {output.outFile} {full_bam} \
		       -p {jacusa_threads} -a D,M,Y,E:file={input.blacklist}:type=BED \
		       -s -m 20 -R {input.fastaref} -P {lib} -F 1024"
		       )

# format jacusa outputs for R
# split multi-allelic sites
# light filtering - total coverage at least 10 reads
# at least 2 edited reads
rule firstFiltering:
	input:
		inFile = projectDir + "{sample}/{sample}.out"
	output:
		outFile = projectDir + "{sample}/{sample}.filt"
	params:
		script = "scripts/firstfiltering.R"
		
	shell:
		# sed -i (work inplace) # find all line that contains 'id' and delete that line. 
		# without this line, running inputfile through Rscript works
		"sed -i '1d' {input.inFile};"
		"ml {R_VERSION};"
		"Rscript {params.script}"
		#adding blank space " --input.."
		" --input {input.inFile}"
		" --output {output.outFile};"
		"rm '*.out.filtered'"	

# merge files together
rule mergeSecondFiltering:
	input:
		expand(projectDir + "{sample}/{sample}.filt", sample = samples)
	output:
		covMat = projectDir + "coverageMatrix.txt",
		ratioMat = projectDir + "ratioMatrix.txt",
		av = projectDir + "avinput.txt"
	params:
		# changing mergeSecondFiltering.R to
		# secondfiltering.R
		script = "scripts/secondfiltering.R",
		# trying-> giving a full path
		inDir = projectDir,
		perc = 0.5,
		er = 0.1
	shell:
		"ml {R_VERSION};"
		"Rscript {params.script}"
		# adding blank space " --inDir.."
		" --inDir {params.inDir}"
		" --rat {output.ratioMat}"
		" --cov {output.covMat}"
		" --av {output.av}"
		" --percSamples {params.perc}"
		" --minER {params.er}"

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
		"table_annovar.pl {input.av} {params.humandb} -buildver hg38"
		" -out {params.outfile} -remove -protocol refGene,dbsnp153CommonSNV,"
		"gnomad30_genome,phastConsElements30way,rmsk,rediportal_012920"
		" -operation g,f,f,r,r,f"
		#" --argument ,,, \'--colsWanted 5\',\'--colsWanted 10&11&12\',"
		# the thread was changed from 10 to 4 
		" -nastring \'.\' --otherinfo --thread {annovar_threads} --maxGeneThread {annovar_threads}"


		# ml annovar
		# table_annovar.pl /sc/arion/projects/ad-omics/flora/cohorts/Navarro_stim/editing-pipeline/editing-pipeline/jacusa-pipeline/snakemake-workspace/avinput.txt /sc/arion/projects/ad-omics/data/references/editing/humandb/ -buildver hg38 -out /sc/arion/projects/ad-omics/flora/cohorts/Navarro_stim/editing-pipeline/editing-pipeline/jacusa-pipeline/snakemake-workspace/myanno -remove -protocol refGene,dbsnp153CommonSNV,gnomad30_genome,phastConsElements30way,rmsk,rediportal_102920 -operation g,f,f,r,r,f

		# The following ( last  argument) causes error. 
		# unrecognized argument
		# --argument ,,, '--colsWanted 5','--colsWanted 10&11&12',

rule annovarFiltering:
	input:
		filtanno = projectDir + "myanno.hg38_multianno.txt",
		covMat = projectDir + "coverageMatrix.txt",
		ratioMat = projectDir + "ratioMatrix.txt"
	output:
		ESanno = projectDir + "ESannotations.txt",
		filtCovMat = projectDir + "filtCoverageMatrix.txt",
		filtRatioMat = projectDir + "filtRatioMatrix.txt"
	params:
		script = "scripts/annovarfiltering.R"
	shell:
		"ml {R_VERSION};"
		"Rscript {params.script}"
		" --inAnno {input.filtanno}"
		" --inRat {input.ratioMat}"
		" --inCov {input.covMat}"
		" --outAnno {output.ESanno}"
		" --outRat {output.filtRatioMat}"
		" --outCov {output.filtCovMat}"
