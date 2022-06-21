#RNA Editing Pipline
#Winston Cuddleston, Flora Seo, Jack Humphrey, Raj Lab
#2022

R_VERSION = "R/4.0.3"
jacusa_threads = 40
annovar_threads = 8
import pandas as pd
import os

#references
refDir = config["refDir"]
editingRefDir = config["editingRefDir"]
humandbDir = config["humandbDir"]
gencode = "/sc/arion/projects/ad-omics/data/references//hg38_reference/GENCODE/gencode.v38.primary_assembly/gencode.v38.primary_assembly.gene_strand.tsv.gz" # just gene_id, gene_name, and strand
#data
projectDir = config["projectDir"]
metadata =  pd.read_csv(config["metadata"], sep = "\t")
samples = metadata['sample']
metadata_dict = metadata.set_index('sample').T.to_dict()

rule all:
    input:
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
        shell("rm {projectDir}/{wildcards.sample}/{wildcards.sample}.out.filtered")

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
        #"sed -i '1d' {input.inFile};"
        "ml {R_VERSION};"
        "Rscript {params.script}"
        #adding blank space " --input.."
        " --input {input.inFile}"
        " --output {output.outFile}"

# merge files together
rule mergeSecondFiltering:
    input:
        expand(projectDir + "{sample}/{sample}.filt", sample = samples)
    output:
        covMat = projectDir + "coverageMatrix.txt",
        ratioMat = projectDir + "ratioMatrix.txt",
        av = projectDir + "avinput.txt"
    params:
        script = "scripts/secondfiltering.R",
        inDir = projectDir,
        perc = 0.5,
        er = 0.1
    shell:
        "ml {R_VERSION};"
        "Rscript {params.script}"
        " --inDir {params.inDir}"
        " --rat {output.ratioMat}"
        " --cov {output.covMat}"
        " --av {output.av}"
        " --percSamples {params.perc}"
        " --minER {params.er}"

# annotate variants
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
        #"-dbtype ensGene " # use Ensembl gene annotation
        " -out {params.outfile} -remove -protocol refGene,dbsnp153CommonSNV,"
        "gnomad30_genome,phastConsElements30way,rmsk,rediportal_012920"
        " -operation g,f,f,r,r,f "
        " -arg '--neargene 10000',,,,, " # flank nearby genes by 10kb
        #" --argument ,,, \'--colsWanted 5\',\'--colsWanted 10&11&12\',"
        # the thread was changed from 10 to 4 
        " -nastring \'.\' --otherinfo --thread {annovar_threads} --maxGeneThread {annovar_threads}"


        # ml annovar
        # table_annovar.pl /sc/arion/projects/ad-omics/flora/cohorts/Navarro_stim/editing-pipeline/editing-pipeline/jacusa-pipeline/snakemake-workspace/avinput.txt /sc/arion/projects/ad-omics/data/references/editing/humandb/ -buildver hg38 -out /sc/arion/projects/ad-omics/flora/cohorts/Navarro_stim/editing-pipeline/editing-pipeline/jacusa-pipeline/snakemake-workspace/myanno -remove -protocol refGene,dbsnp153CommonSNV,gnomad30_genome,phastConsElements30way,rmsk,rediportal_102920 -operation g,f,f,r,r,f

        # The following ( last  argument) causes error. 
        # unrecognized argument
        # --argument ,,, '--colsWanted 5','--colsWanted 10&11&12',

# more filtering
rule annovarFiltering:
    input:
        filtanno = projectDir + "myanno.hg38_multianno.txt",
        covMat = projectDir + "coverageMatrix.txt",
        ratioMat = projectDir + "ratioMatrix.txt"
    output:
        ESanno = projectDir + "editing_annotation.txt",
        filtCovMat = projectDir + "filtCoverageMatrix.txt",
        filtRatioMat = projectDir + "filtRatioMatrix.txt",
        outBed = projectDir + "editing_sites.bed"
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
        " --outBed {output.outBed}"
        " --gencode {gencode}"
