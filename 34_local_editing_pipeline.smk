# This is a local pipeline written for debug purposes
# This is a pipeline that contains rule filtering and annovar, (3rd and 4th rules)
# This was made to be used with the complete intput file 'SAMPLE.rediout.ga' made from rule 1 and rule 2 in pipie 12_local..smk. 
# adding change to see if git add works
# Error report
# # 12_local_editing_pipeline.smk
# # reditools_1 : works fine (input SAMPLE.bam)
# # reditools_2 : doesn't work in pipeline (with snakemake -s.. command)
# #             : however, the code works fine outside of the pipe (input SAMPLE.bam, output SAMPLE.rediout.gz)
# #             : however, when feeding the output gz file from the command outside of the terminal into the same pipeline to proceed on rule 3 and 4, the snakemake recognizes the file to be incompelte. won'r run
# #             : therefore, the pipe with four rules were divided into two, 12_ containing rule reditools_1 and reditools_2, and 34_ containing rule filtering and annovar.
#
# # 34_local_editing_pipeline.smk 
# # filtering and annovar   : in 34_ pipeline. Feeding the SAMPLE.rediout.gz created outside the pipe(reditools_2) to rule filtering and annovar in 34_ pipeline works fine (input SAMPLE.rediout.gz, output SAMPLE.annovar.hg38_multianno.csv)

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
		#expand(OUTFOLDER + "{SAMPLE}.rediout.gz", SAMPLE = SAMPLES)
		expand(OUTFOLDER + "{SAMPLE}.annovar.hg38_multianno.csv", SAMPLE = SAMPLES)

rule filtering:
	input: OUTFOLDER + "{SAMPLE}.rediout.gz"
	output: OUTFOLDER + "{SAMPLE}.reditools.filtered.gz"
	params:
		perfix = OUTFOLDER + "{SAMPLE}",
		filter_awk = "/sc/arion/projects/ad-omics/flora/cohorts/Efthymiou_stim/debug/scripts/filter_reditools.awk",
		filter_python = "/sc/arion/projects/ad-omics/flora/cohorts/Efthymiou_stim/debug/scripts/filter_DP.py",
		dbsnp = "/sc/arion/projects/ad-omics/data/references/hg38_reference/dbSNP/hg38_common_snps_dbSNP_153.bed.gz"
	shell: 
		"module load annovar;"
		"sh {params.filter_awk} {input} {params.perfix}_annovar.input;"
	
	# "sh scripts/filter_reditools.awk test/test.reditout.gz test/output/test_annovar.input" 
		"module load python;"
		"python {params.filter_python} {params.perfix}_annovar.input {CHANGE_DP} {LOCUS_DP} > {params.perfix}_annovar_filtered.input;"
		# "python scripts/filter_DP.py test/output/test_annovar.input 3 5 > test/output/test_annovar_filtered.input"

		"ml bedtools;"
		"bedtools intersect -v -a {params.perfix}_annovar_filtered.input -b {params.dbsnp} > {params.perfix}.reditools.filtered.gz;"
		# "bedtools intersect -v -a test/output/test_annovar_filtered.input -b /sc/arion/projects/ad-omics/data/references/hg38_reference/dbSNP/hg38_common_snps_dbSNP_153.bed.gz > test/output/test.reditools.filtered.gz"

rule annovar: 
	input: OUTFOLDER + "{SAMPLE}.reditools.filtered.gz"
	output: OUTFOLDER + "{SAMPLE}.annovar.hg38_multianno.csv"
	params:
		perfix = OUTFOLDER + "{SAMPLE}",
		humandb = "/sc/arion/projects/buxbaj01a/resources/annovar/humandb",
		threads = 4
	shell:
		"module load annovar;"
		#"table_annovar.pl {input} {params.humandb} -buildver hg38 -remove -protocal refGene,dbsnp153CommonSNV,gnomad30_genome,rmsk,rediportal_012920 -operation g,f,f,r,f -nastring . --otherinfo --thread {threads} --maxgenethread {threads} -csvout -outfile {params.perfix}.annovar;"
		"table_annovar.pl {input} {params.humandb}/ --buildver hg38 --out --remove --protocol refGene,dbsnp153CommonSNV,gnomad30_genome,rmsk,rediportal_012920 --operation g,f,f,r,f --nastring . --otherinfo --thread {params.threads} --maxgenethread {params.threads} --csvout --outfile {params.perfix}.annovar;"
		# "table_annovar.pl test/output/test.reditools.filtered.gz /sc/arion/projects/buxbaj01a/resources/annovar/humandb -buildver hg38 -remove -protocal refGene,dbsnp153CommonSNV,gnomad30_genome,rmsk,rediportal_012920 -operation g,f,f,r,f -nastring . --otherinfo --thread 4 --maxgenethread 4 -csvout -outfile test/output/test.annovar"
