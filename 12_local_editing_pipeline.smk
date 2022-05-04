# This is a local pipeline written for debug purposes
# This is a copy of a local pipeline to debug the rule 1 by seperating rule 1 to rule1_1 and rule1_2

# Error report
# 12_local_editing_pipeline.smk
# reditools_1 : works fine (input SAMPLE.bam)
# reditools_2 : doesn't work in pipeline (with snakemake -s.. command)
#	      : however, the code works fine outside of the pipe (input SAMPLE.bam, output SAMPLE.rediout.gz)
#	      : however, when feeding the output gz file from the command outside of the terminal into the same pipeline to proceed on rule 3 and 4, the snakemake recognizes the file to be incompelte. won'r run
#	      : therefore, the pipe with four rules were divided into two, 12_ containing rule reditools_1 and reditools_2, and 34_ containing rule filtering and annovar.

# 34_local_editing_pipeline.smk 
# filtering and annovar   : in 34_ pipeline. Feeding the SAMPLE.rediout.gz created outside the pipe(reditools_2) to rule filtering and annovar in 34_ pipeline works fine (input SAMPLE.rediout.gz, output SAMPLE..annovar.csv)


import os
import pandas as pd

reditools = config["reditools"]
INFOLDER = config["INFOLDER"]
OUTFOLDER = config["OUTFOLDER"]

threads = config["threads"]
STRAND = config["STRAND"]
CHANGE_DP = config["CHANGE_DP"]
LOCUS_DP = config["LOCUS_DP"]

#for processing multiple samples
##metadata = pd.read_csv(config["metadata"], sep="\t")
##SAMPLES =  metadata["SAMPLE"]
##BAMS = metadata["bam_path"] + metadata["SAMPLE"] + ".bam"
##metadata_dict = metadata.set_index('SAMPLE').T.to_dict()

# fix
SAMPLES = ["new_test"]

rule all:
	input:
		expand(OUTFOLDER + "{SAMPLE}.rediout.gz", SAMPLE = SAMPLES)
		#expand(OUTFOLDER + "{SAMPLE}.annovar.hg38_multianno.csv", SAMPLE = SAMPLES)

rule reditools_1:
	input: INFOLDER + "{SAMPLE}.bam"
	output: OUTFOLDER + "{SAMPLE}.rediout_cov"
	# ouptut of reditools_1 is not used as input of reditools_2, it is just assigned for the form of the snakemake pipe. reditools_1 creats multiple files under new_temp_cov, and those are used in reditools_2, but it is assigned in the code line, not as input ( input is still the bam file)
	params:
		REFINDEX = "/sc/arion/projects/ad-omics/data/references/hg38_reference/hg38.standard_chr.fa.fai"
	shell:
		"module load reditools;"
		"{reditools}/extract_coverage_dynamic.sh {input} `pwd`/new_temp_cov/ {params.REFINDEX} > {output};"
		#"mkdir -p new_temp_cov;"
		#"mkdir -p new_temp_results;"

		# with new test bam
		# /hpc/packages/minerva-centos7/reditools/2.0/reditools2.0/extract_coverage_dynamic.sh test/new_test.bam new_temp_cov/ /sc/arion/projects/ad-omics/data/references/hg38_reference/hg38.standard_chr.fa.fai
		
		# /hpc/packages/minerva-centos7/reditools/2.0/reditools2.0/extract_coverage_dynamic.sh `pwd`/test/test.bam `pwd`/temp_cov/ /sc/arion/projects/ad-omics/data/references/hg38_reference/hg38.standard_chr.fa.fai
	                                 
# rule reditools_2 stops to continue 

rule reditools_2:
	input: INFOLDER + "{SAMPLE}.bam"
	output: OUTFOLDER + "{SAMPLE}.rediout.gz"
	params:
		REFINDEX = "/sc/arion/projects/ad-omics/data/references/hg38_reference/hg38.standard_chr.fa.fai",
                REF = "/sc/arion/projects/ad-omics/data/references/hg38_reference/hg38.fa"
	shell:
		"module load reditools;"
                "module load htslib/1.7;" # not needed 
                "module load samtools;" # not needed
                "mpirun -np {threads} {reditools}/src/cineca/parallel_reditools.py -f {input} -r {params.REF} -S -s {STRAND} -ss {LOCUS_DP} -mrl 50 -q 10 -bq 20 -C -T 2 -c -m /sc/arion/projects/ad-omics/data/references/editing/homopolymeric_sites_hg38.txt -os {LOCUS_DP} -t `pwd`/new_temp_results/ -Z {params.REFINDEX} -G `pwd`/new_temp_cov/{wildcards.SAMPLE}.cov -D `pwd`/new_temp_cov/ -o {wildcards.SAMPLE}_out;"
		"{reditools}/merge.sh `pwd`/new_temp_results/ {output} {threads};"

                # this param is not used bc it has error calling it
                # HOMPOL = "/sc/arion/projects/ad-omics/data/references/editing/homopolymeric_sites_hg38.txt"
		
		# with new test bam (8 nodes, chr12, smaller)
		# mpirun -np 8 /hpc/packages/minerva-centos7/reditools/2.0/reditools2.0/src/cineca/parallel_reditools.py -f `pwd`/test/new_test.bam -r /sc/arion/projects/ad-omics/data/references/hg38_reference/hg38.fa -S -s 0 -ss 5 mrl 50 -q 10 -bq 20 -C -T 2 -c -m /sc/arion/projects/ad-omics/data/references/editing/homopolymeric_sites_hg38.txt -os 5 -t `pwd`/new_temp_results/ -Z /sc/arion/projects/ad-omics/data/references/hg38_reference/hg38.standard_chr.fa.fai -G `pwd`/new_temp_cov/new_test.cov -D `pwd`/new_temp_cov/ -o new_test_out 
		
		# with old test bam (8 nodes, chr21, bigger)
		# mpirun -np 8 /hpc/packages/minerva-centos7/reditools/2.0/reditools2.0/src/cineca/parallel_reditools.py -f test/test.bam -r /sc/arion/projects/ad-omics/data/references/hg38_reference/hg38.fa -S -s 0 -ss 5 mrl 50 -q 10 -bq 20 -C -T 2 -c -m /sc/arion/projects/ad-omics/data/references/editing/homopolymeric_sites_hg38.txt -os 5 -t temp_results/ -Z /sc/arion/projects/ad-omics/data/references/hg38_reference/hg38.standard_chr.fa.fai -G temp_cov/test.cov -D temp_cov/ -o test_out 

		 # translation
		 # mpirun -np 8 /hpc/packages/minerva-centos7/reditools/2.0/reditools2.0/src/cineca/parallel_reditools.py -f(FILE) test/test.bam -r(REFERENCE) /sc/arion/projects/ad-omics/data/references/hg38_reference/hg38.fa -S(STRICT_ACTIVATE) -s(STRAND) 0 -ss(SPLICING_SPAN) 5 -mrl(MIN_READ_LENGTH) 50 -q(MIN_READ_QUALITY) 10 -bq(MIN_BASE_QUALITY) 20 -C(STRAND_CORRECTION) -T(STRAND_CONFIDENCE) 2 -c(CREATE_OMOPOLYMERIC_FILE) -m(OMOPOLYMERIC_FILE) /sc/arion/projects/ad-omics/data/references/editing/homopolymeric_sites_hg38.txt -os(OMOPOLYMERIC_SPAN) 5 -t(TEMP_DIR) `pwd`/temp_results/ -Z(CREATE_CHR_SIZE_FILE) /sc/arion/projects/ad-omics/data/references/hg38_reference/hg38.standard_chr.fa.fai -G(COVERAGE_FILE) `pwd`/temp_cov/test(*.cov).cov -D(COVERAGE_DIR) temp_cov/ -o test_out 

		#"module load htslib/1.7;"
		#"module load samtools;"
		#"{reditools}/merge.sh `pwd`/new_temp_results/ {output} {threads};"
		#"bgzip -d {output};"
		#"rm -fr temp_cov;"
		#"rm -fr temp_results;"

#rule filtering:
#	input: OUTFOLDER + "{SAMPLE}.rediout.gz"
#	output: OUTFOLDER + "{SAMPLE}.reditools.filtered.gz"
#	params:
#		perfix = OUTFOLDER + "{SAMPLE}",
#		filter_awk = "/sc/arion/projects/ad-omics/flora/cohorts/Efthymiou_stim/debug/scripts/filter_reditools.awk",
#		filter_python = "/sc/arion/projects/ad-omics/flora/cohorts/Efthymiou_stim/debug/scripts/filter_DP.py",
#		dbsnp = "/sc/arion/projects/ad-omics/data/references/hg38_reference/dbSNP/hg38_common_snps_dbSNP_153.bed.gz"
#	shell: 
#		"module load annovar;"
#		"sh {params.filter_awk} {input} {params.perfix}_annovar.input;"
#	
#	# "sh scripts/filter_reditools.awk test/test.reditout.gz test/output/test_annovar.input" 
#		"module load python;"
#		"python {params.filter_python} {params.perfix}_annovar.input {CHANGE_DP} {LOCUS_DP} > {params.perfix}_annovar_filtered.input;"
#		# "python scripts/filter_DP.py test/output/test_annovar.input 3 5 > test/output/test_annovar_filtered.input"
#
#		"ml bedtools;"
#		"bedtools intersect -v -a {params.perfix}_annovar_filtered.input -b {params.dbsnp} > {params.perfix}.reditools.filtered.gz;"
#		# "bedtools intersect -v -a test/output/test_annovar_filtered.input -b /sc/arion/projects/ad-omics/data/references/hg38_reference/dbSNP/hg38_common_snps_dbSNP_153.bed.gz > test/output/test.reditools.filtered.gz"

#rule annovar: 
#	input: OUTFOLDER + "{SAMPLE}.reditools.filtered.gz"
#	output: OUTFOLDER + "{SAMPLE}.annovar.hg38_multianno.csv"
#	params:
#		perfix = OUTFOLDER + "{SAMPLE}",
#		humandb = "/sc/arion/projects/buxbaj01a/resources/annovar/humandb",
#		threads = 4
#	shell:
#		"module load annovar;"
#		#"table_annovar.pl {input} {params.humandb} -buildver hg38 -remove -protocal refGene,dbsnp153CommonSNV,gnomad30_genome,rmsk,rediportal_012920 -operation g,f,f,r,f -nastring . --otherinfo --thread {threads} --maxgenethread {threads} -csvout -outfile {params.perfix}.annovar;"
#		"table_annovar.pl {input} {params.humandb}/ --buildver hg38 --out --remove --protocol refGene,dbsnp153CommonSNV,gnomad30_genome,rmsk,rediportal_012920 --operation g,f,f,r,f --nastring . --otherinfo --thread {params.threads} --maxgenethread {params.threads} --csvout --outfile {params.perfix}.annovar;"
#		# "table_annovar.pl test/output/test.reditools.filtered.gz /sc/arion/projects/buxbaj01a/resources/annovar/humandb -buildver hg38 -remove -protocal refGene,dbsnp153CommonSNV,gnomad30_genome,rmsk,rediportal_012920 -operation g,f,f,r,f -nastring . --otherinfo --thread 4 --maxgenethread 4 -csvout -outfile test/output/test.annovar"
