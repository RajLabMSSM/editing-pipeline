#RNA Editing Pipline
#Winston Cuddleston, Flora Seo, Jack Humphrey, Raj Lab
#2022

R_VERSION="R/4.0.3"
jacusa_threads=40
jacusa_multi_threads=10
annovar_threads=8
import pandas as pd
import os

# references
refDir=config["refDir"]
editingRefDir=config["editingRefDir"]
humandbDir=config["humandbDir"]
gencodeDir=config["gencodeDir"]
metaDF=config["metadata"]


# data
projectDir=config["projectDir"]
metadata= pd.read_csv(config["metadata"], sep="\t")
samples=metadata['sample']
metadata_dict=metadata.set_index('sample').T.to_dict()

# variables - put these in config, experiment with relaxing
min_coverage=config["min_coverage"]
min_edit_rate=config["min_edit_rate"]
filter_missing_samples=config["sample_missingness"]
filter_missing_sites=config["site_missingness"]

chromosomes=["chr" + str(i) for i in range(1,23) ] + ["chrX", "chrY"]

rule all:
    input:
        projectDir + "all_sites_pileup_coverage.tsv.gz",
        projectDir + "all_sites_pileup_editing.tsv.gz",
        projectDir + "all_sites_pileup_annotation.tsv.gz"

# de novo calling of RNA editing genome-wide
# remove blacklist regions
# at some point offer ability to filter on VCF        
rule jacusa:
    input:
        blacklist = editingRefDir + "hg38-blacklist.v2_sort.bed",
        fastaref = refDir + "GRCh38.primary_assembly.genome.fa"
    output:
        outFile = projectDir + "{sample}/{sample}.out"
    run:
        bam_file = os.path.join(metadata_dict[wildcards.sample]["bam_path"], wildcards.sample + ".bam")
        lib = metadata_dict[wildcards.sample]["library"]
        shell("ml jacusa2; \
               java -jar $JACUSA2_JAR call-1 \
               -r {output.outFile} {bam_file} \
               -p {jacusa_threads} -a D,M,Y,E:file={input.blacklist}:type=BED \
               -s -m 20 -R {input.fastaref} -P {lib} -F 1024"
               )
        shell("rm {projectDir}/{wildcards.sample}/{wildcards.sample}.out.filtered")

# supervised calling of RNA editing sites from REDIportal
rule known_sites:
    input:
        bed=editingRefDir + "REDIsites.bed"
    output:
        outFile=projectDir + "{sample}/{sample}.knownES.out"
    run:
        bam_file=os.path.join(metadata_dict[wildcards.sample]["bam_path"], wildcards.sample + ".bam")
        lib=metadata_dict[wildcards.sample]["library"]
        shell("ml jacusa2; \
        java -jar $JACUSA2_JAR call-1 \
        -r {output.outFile} -b {input.bed} {bam_file} \
        -p {jacusa_multi_threads} \
        -a D,Y -P {lib} -F 1024 \
        -A -s -m 20")
        shell("rm {projectDir}/{wildcards.sample}/{wildcards.sample}.knownES.out.filtered")

# format jacusa outputs for R
# first known sites
# split multi-allelic sites
# light filtering - less stringent than de novo - total coverage at least 5 reads
# at least 3 edited reads
# minimum editing rate > 0.01
rule parse_known:
    input:
        inFile=projectDir + "{sample}/{sample}.knownES.out"
    output:
        outFile=projectDir + "{sample}/{sample}.knownES.filt"
    params:
        script="scripts/parse_knownSites.R"
    shell:
        "ml {R_VERSION};"
        "Rscript {params.script} "
        "--input {input.inFile} "
        "--output {output.outFile}"

# then de novo sites
# split multi-allelic sites
# light filtering - more stringent than known - total coverage at least 10 reads
# at least 3 edited reads
# minimum editing rate > 0.01
rule parse_jacusa:
    input:
        inFile=projectDir + "{sample}/{sample}.out"
    output:
        outFile=projectDir + "{sample}/{sample}.filt"
    params:
        script="scripts/parse_jacusa.R"
    shell:
        "ml {R_VERSION};"
        "Rscript {params.script} "
        "--input {input.inFile} "
        "--output {output.outFile}"

# concatenate de novo sites and known sites
rule cat_all_sites:
    input:
        deNovo=projectDir + "{sample}/{sample}.filt",
        known=projectDir + "{sample}/{sample}.knownES.filt"
    output:
        outFile=projectDir + "{sample}/{sample}.filtAll"
    params:
        script="scripts/concatenate_sites.R"
    shell:
        "ml {R_VERSION};"
        "Rscript {params.script} "
        "--deNovoIN {input.deNovo} "
        "--knownIN {input.known} "
        "--output {output.outFile}"

# merge files together
rule merge_jacusa:
    input:
        expand(projectDir + "{sample}/{sample}.filtAll", sample=samples)
    output:
        projectDir + "merge/{chrom}_coverage.tsv.gz",
        projectDir + "merge/{chrom}_ratio.tsv.gz"
    params:
        script="scripts/merge_jacusa.R" 
    shell:
        "ml {R_VERSION};"
        "Rscript {params.script} "
        "--inDir {projectDir} "
        "--chr {wildcards.chrom}"

# apply cohort level filters
rule filter_cohort:
    input:
        projectDir + "merge/{chrom}_coverage.tsv.gz",
        projectDir + "merge/{chrom}_ratio.tsv.gz"
    output:
        covMat=projectDir + "filter/{chrom}_coverage.tsv.gz",
        ratMat=projectDir + "filter/{chrom}_ratio.tsv.gz"
    params:
        script="scripts/filter_cohort.R",
        N=min_coverage,
        er=min_edit_rate
    shell:
        "ml {R_VERSION};"
        "Rscript {params.script} "
        "--inDir {projectDir} "
        "--chr {wildcards.chrom} "
        "--minSamples {params.N} "
        "--minER {params.er} "
        "--group {metaDF} "

# concatenate filtered chunks together
# write out VCF format for ANNOVAR
rule concatenate_cohort:
    input:
        expand(projectDir + "filter/{chrom}_coverage.tsv.gz", chrom=chromosomes),
        expand(projectDir + "filter/{chrom}_ratio.tsv.gz", chrom=chromosomes)
    output:
        cov=projectDir + "all_sites_coverage.tsv.gz",
        rat=projectDir + "all_sites_ratio.tsv.gz",
        av=projectDir + "avinput.txt" 
    params:
        script="scripts/concatenate_cohort.R"
    shell:
        "ml {R_VERSION};"
        "Rscript {params.script} "
        "--inDir {projectDir}"
        
# annotate variants
# look into using GENCODE instead of refGene at some point
rule annovar:
    input:
        av=projectDir + "avinput.txt"
    output:
        anno=projectDir +  "myanno.hg38_multianno.txt"
    params:
        outfile=projectDir + "myanno",
        humandb=humandbDir
    shell:
        "ml annovar;"
        "table_annovar.pl {input.av} {params.humandb} -buildver hg38"
        " -out {params.outfile} -remove -protocol refGene,dbsnp153CommonSNV,"
        "gnomad30_genome,phastConsElements30way,rmsk,rediportal_012920"
        " -operation g,f,f,r,r,f "
        " -arg '--neargene 10000',,,,, " # flank nearby genes by 10kb
        " -nastring \'.\' --otherinfo --thread {annovar_threads} --maxGeneThread {annovar_threads}" # the thread was changed from 10 to 4 

# more filtering
# don't need to output matrices at this point, only annotation and BED files needed
rule filter_annovar:
    input:
        filtanno=projectDir + "myanno.hg38_multianno.txt",
        cov=projectDir + "all_sites_coverage.tsv.gz",
        rat=projectDir + "all_sites_ratio.tsv.gz",
        gencode=gencodeDir + "gencode.v38.primary_assembly.gene_strand.tsv.gz"
    output:
        ESanno=projectDir + "editing_annotation.txt",
        outBed=projectDir + "editing_sites.bed"
    params:
        script="scripts/filter_annovar.R"
    shell:
        "ml {R_VERSION};"
        "Rscript {params.script} "
        "--inAnno {input.filtanno} "
        "--inRat {input.rat} "
        "--inCov {input.cov} "
        "--outAnno {output.ESanno} "
        "--outBed {output.outBed} "
        "--gencode {input.gencode}"

# get counts of each filtered site in each sample
# run jacusa with set of sites to fill in missing data
rule jacusa_pileup:
    input:
        bed=projectDir + "editing_sites.bed"
    output:
        outFile=projectDir + "{sample}/{sample}.pileup.txt"
    run:
        bam_file = os.path.join(metadata_dict[wildcards.sample]["bam_path"], wildcards.sample + ".bam")
        lib = metadata_dict[wildcards.sample]["library"]
        shell("ml jacusa2; java -jar $JACUSA2_JAR call-1 \
            -r {output.outFile} -b {input.bed} {bam_file} \
            -p {jacusa_multi_threads} \
            -a D,Y -P {lib} -F 1024 \
            -A -s -m 20")
        shell("rm {projectDir}/{wildcards.sample}/{wildcards.sample}.pileup.txt.filtered")    

# format final outputs for R
# split multi-allelic sites
# create DTU file here?
rule parse_pileup:
    input:
        pu=projectDir + "{sample}/{sample}.pileup.txt"
    output:
        parsed=projectDir + "{sample}/{sample}.pileup.parsed.txt"
    params:
        script="scripts/parse_pileup.R"
    shell:
        "ml {R_VERSION};"
        "Rscript {params.script} "
        "--input {input.pu} "
        "--output {output.parsed}"

# merge pileups
# final annotation - flip allele orientation based on gene
# any final filtering? sample missingness?
rule merge_pileup:
    input:
        expand(projectDir + "{sample}/{sample}.pileup.parsed.txt", sample=samples),
        anno=projectDir + "editing_annotation.txt"
    output:
        projectDir + "all_sites_pileup_coverage.tsv.gz",
        projectDir + "all_sites_pileup_editing.tsv.gz",
        projectDir + "all_sites_pileup_annotation.tsv.gz",
        projectDir + "all_sites_pileup_dtu.tsv.gz"
    params:
        script="scripts/merge_pileup.R" 
    shell:
        "ml {R_VERSION};"
        "Rscript {params.script} "
        "--inDir {projectDir} "
        "--anno {input.anno} "
        "--missing_sites {filter_missing_sites} "
        "--missing_samples {filter_missing_samples} "
