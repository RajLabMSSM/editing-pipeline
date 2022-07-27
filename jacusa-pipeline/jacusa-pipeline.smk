#RNA Editing Pipline
#Winston Cuddleston, Flora Seo, Jack Humphrey, Raj Lab
#2022

R_VERSION = "R/4.0.3"
jacusa_threads = 40
jacusa_multi_threads = 10
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

# variables - put these in config, experiment with relaxing
min_coverage = config["min_coverage"]
min_edit_rate = config["min_edit_rate"]

chromosomes = ["chr" + str(i) for i in range(1,23) ] + ["chrX", "chrY"]

rule all:
    input:
        projectDir + "all_sites_pileup_coverage.tsv.gz",
        projectDir + "all_sites_pileup_editing.tsv.gz",
        projectDir + "all_sites_pileup_annotation.tsv.gz"
        #projectDir + "filtCoverageMatrix.txt",
        #projectDir + "filtRatioMatrix.txt" #,
        #expand(projectDir + "{sample}/jacusa_multi_pileup.txt", sample = samples)

# call RNA editing genome-wide
# remove blacklist regions
# at some point offer ability to filter on VCF
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
# minimum editing rate > 0.01
rule parse_jacusa:
    input:
        inFile = projectDir + "{sample}/{sample}.out"
    output:
        outFile = projectDir + "{sample}/{sample}.filt"
    params:
        script = "scripts/parse_jacusa.R"
        
    shell:
        "ml {R_VERSION};"
        "Rscript {params.script}"
        " --input {input.inFile}"
        " --output {output.outFile}"

# merge files together 
rule merge_jacusa:
    input:
        expand(projectDir + "{sample}/{sample}.filt", sample = samples)
    output:
        projectDir + "merge/{chrom}_coverage.tsv.gz",
        projectDir + "merge/{chrom}_ratio.tsv.gz"
    params:
        script = "scripts/merge_jacusa.R" 
    shell:
        "ml {R_VERSION};"
        "Rscript {params.script} --inDir {projectDir} --chr {wildcards.chrom} "


# apply cohort level filters
rule filter_cohort:
    input:
        projectDir + "merge/{chrom}_coverage.tsv.gz",
        projectDir + "merge/{chrom}_ratio.tsv.gz"
    output:
        covMat = projectDir + "filter/{chrom}_coverage.tsv.gz",
        ratMat = projectDir + "filter/{chrom}_ratio.tsv.gz"
    params:
        script = "scripts/filter_cohort.R",
        perc = min_coverage,
        er = min_edit_rate
    shell:
        "ml {R_VERSION};"
        "Rscript {params.script}"
        " --inDir {projectDir}"
        " --chr {wildcards.chrom} "
        " --percSamples {params.perc}"
        " --minER {params.er}"

# concatenate filtered chunks together
# write out VCF format for ANNOVAR
rule concatenate_cohort:
    input:
        expand(projectDir + "filter/{chrom}_coverage.tsv.gz", chrom = chromosomes),
        expand(projectDir + "filter/{chrom}_ratio.tsv.gz", chrom = chromosomes)
    output:
        cov = projectDir + "all_sites_coverage.tsv.gz",
        rat = projectDir + "all_sites_ratio.tsv.gz",
        av = projectDir + "avinput.txt" 
    params:
        script = "scripts/concatenate_cohort.R"
    shell:
        "ml {R_VERSION};"
        "Rscript {params.script} "
        " --inDir {projectDir} "
        
# annotate variants
# look into using GENCODE instead of refGene at some point
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

# more filtering
# don't need to output matrices at this point, only annotation and BED files needed
rule filter_annovar:
    input:
        filtanno = projectDir + "myanno.hg38_multianno.txt",
        cov = projectDir + "all_sites_coverage.tsv.gz",
        rat = projectDir + "all_sites_ratio.tsv.gz"
    output:
        ESanno = projectDir + "editing_annotation.txt",
        filtCovMat = projectDir + "filtCoverageMatrix.txt",
        filtRatioMat = projectDir + "filtRatioMatrix.txt",
        outBed = projectDir + "editing_sites.bed"
    params:
        script = "scripts/filter_annovar.R"
    shell:
        "ml {R_VERSION};"
        "Rscript {params.script}"
        " --inAnno {input.filtanno}"
        " --inRat {input.rat}"
        " --inCov {input.cov}"
        " --outAnno {output.ESanno}"
        " --outRat {output.filtRatioMat}"
        " --outCov {output.filtCovMat}"
        " --outBed {output.outBed}"
        " --gencode {gencode}"

# get counts of each filtered site in each sample

# run jacusa with set of sites to fill in missing data
rule jacusa_pileup:
    input:
        bed = projectDir + "editing_sites.bed"
    output:
        pu = projectDir + "{sample}/{sample}_pileup.txt",
        parsed = projectDir + "{sample}/{sample}_pileup_parsed.txt"
    params:
        script = "scripts/parse_pileup.R"
    run:
        bam_file = os.path.join(metadata_dict[wildcards.sample]["bam_path"], wildcards.sample + ".bam")
        
        shell("ml jacusa2; java -jar $JACUSA2_JAR call-1 \
            -r {output.pu} -b {input.bed} {bam_file} \
            -p {jacusa_multi_threads}\
            -a D -F 1024 \
            -A; \
            ml {R_VERSION};\
            Rscript {params.script} --input {output.pu} --output {output.parsed}")

# merge pileups
# final annotation - flip allele orientation based on gene
# any final filtering? sample missingness?
rule merge_pileup:
    input:
        expand(projectDir + "{sample}/{sample}_pileup_parsed.txt", sample = samples),
        anno = projectDir + "editing_annotation.txt"
    output:
        projectDir + "all_sites_pileup_coverage.tsv.gz",
        projectDir + "all_sites_pileup_editing.tsv.gz",
        projectDir + "all_sites_pileup_annotation.tsv.gz"
    params:
        script = "scripts/merge_pileup.R" 
    shell:
        "ml {R_VERSION};"
        "Rscript {params.script} --inDir {projectDir} --anno {input.anno} "
