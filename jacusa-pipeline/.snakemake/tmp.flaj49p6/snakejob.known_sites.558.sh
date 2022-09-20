#!/bin/sh
# properties = {"type": "single", "rule": "known_sites", "local": false, "input": ["/sc/arion/projects/ad-omics/data/references/editing/REDIsites.bed"], "output": ["/sc/arion/projects/ad-omics/winston/bigbrain/editing/AnswerALS/Run6/CASE-NEUKP160HX8-6874-T_P45/CASE-NEUKP160HX8-6874-T_P45.knownES.out"], "wildcards": {"sample": "CASE-NEUKP160HX8-6874-T_P45"}, "params": {}, "log": [], "threads": 1, "resources": {"tmpdir": "/tmp"}, "jobid": 558, "cluster": {"queue": "express", "cores": 8, "mem": 10000, "time": "720", "name": "$(basename $(pwd)):known_sites:sample=CASE-NEUKP160HX8-6874-T_P45", "output": "logs/known_sites:sample=CASE-NEUKP160HX8-6874-T_P45.stdout", "error": "logs/known_sites:sample=CASE-NEUKP160HX8-6874-T_P45.stderr", "himem": ""}}
 cd /sc/arion/projects/ad-omics/winston/editing-pipeline/jacusa-pipeline && \
/sc/arion/projects/als-omics/conda/envs/snakemake/bin/python3.6 \
-m snakemake /sc/arion/projects/ad-omics/winston/bigbrain/editing/AnswerALS/Run6/CASE-NEUKP160HX8-6874-T_P45/CASE-NEUKP160HX8-6874-T_P45.knownES.out --snakefile /sc/arion/projects/ad-omics/winston/editing-pipeline/jacusa-pipeline/jacusa-pipeline-winston.smk \
--force --cores all --keep-target-files --keep-remote --max-inventory-time 0 \
--wait-for-files /sc/arion/projects/ad-omics/winston/editing-pipeline/jacusa-pipeline/.snakemake/tmp.flaj49p6 /sc/arion/projects/ad-omics/data/references/editing/REDIsites.bed --latency-wait 30 \
 --attempt 1 --force-use-threads --scheduler ilp \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
 --configfiles /sc/arion/projects/ad-omics/winston/editing-pipeline/jacusa-pipeline/config.yaml --config mode=Null  --allowed-rules known_sites --nocolor --notemp --no-hooks --nolock --scheduler-solver-path /sc/arion/projects/als-omics/conda/envs/snakemake/bin \
--mode 2  --use-conda  --conda-base-path /hpc/packages/minerva-centos7/anaconda3/2020.11  --default-resources "tmpdir=system_tmpdir" 

