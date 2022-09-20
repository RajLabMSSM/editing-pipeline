#!/bin/sh
# properties = {"type": "single", "rule": "jacusa", "local": false, "input": ["/sc/arion/projects/ad-omics/data/references/editing/hg38-blacklist.v2_sort.bed", "/sc/arion/projects/ad-omics/data/references/GRCh38_references/GRCh38.primary_assembly.genome.fa"], "output": ["/sc/arion/projects/ad-omics/winston/bigbrain/editing/AnswerALS/Run6/CASE-NEUPG593UTV-7624-T_P68/CASE-NEUPG593UTV-7624-T_P68.out"], "wildcards": {"sample": "CASE-NEUPG593UTV-7624-T_P68"}, "params": {}, "log": [], "threads": 1, "resources": {"tmpdir": "/tmp"}, "jobid": 711, "cluster": {"queue": "express", "cores": 4, "mem": 3750, "time": "720", "name": "$(basename $(pwd)):jacusa:sample=CASE-NEUPG593UTV-7624-T_P68", "output": "logs/jacusa:sample=CASE-NEUPG593UTV-7624-T_P68.stdout", "error": "logs/jacusa:sample=CASE-NEUPG593UTV-7624-T_P68.stderr", "himem": ""}}
 cd /sc/arion/projects/ad-omics/winston/editing-pipeline/jacusa-pipeline && \
/sc/arion/projects/als-omics/conda/envs/snakemake/bin/python3.6 \
-m snakemake /sc/arion/projects/ad-omics/winston/bigbrain/editing/AnswerALS/Run6/CASE-NEUPG593UTV-7624-T_P68/CASE-NEUPG593UTV-7624-T_P68.out --snakefile /sc/arion/projects/ad-omics/winston/editing-pipeline/jacusa-pipeline/jacusa-pipeline-winston.smk \
--force --cores all --keep-target-files --keep-remote --max-inventory-time 0 \
--wait-for-files /sc/arion/projects/ad-omics/winston/editing-pipeline/jacusa-pipeline/.snakemake/tmp.flaj49p6 /sc/arion/projects/ad-omics/data/references/editing/hg38-blacklist.v2_sort.bed /sc/arion/projects/ad-omics/data/references/GRCh38_references/GRCh38.primary_assembly.genome.fa --latency-wait 30 \
 --attempt 1 --force-use-threads --scheduler ilp \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
 --configfiles /sc/arion/projects/ad-omics/winston/editing-pipeline/jacusa-pipeline/config.yaml --config mode=Null  --allowed-rules jacusa --nocolor --notemp --no-hooks --nolock --scheduler-solver-path /sc/arion/projects/als-omics/conda/envs/snakemake/bin \
--mode 2  --use-conda  --conda-base-path /hpc/packages/minerva-centos7/anaconda3/2020.11  --default-resources "tmpdir=system_tmpdir" 
