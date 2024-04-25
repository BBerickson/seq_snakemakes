#!/bin/sh
# properties = {"type": "single", "rule": "fastqc", "local": false, "input": ["Temp4/fastqs/PP1H66KpDoxS2Ser5_CKDL230028616-1A_HCLMNDSX7_L4_1.fq.gz"], "output": ["Temp4/fastqc/PP1H66KpDoxS2Ser5_CKDL230028616-1A_HCLMNDSX7_L4_1_fastqc.zip", "Temp4/fastqc/PP1H66KpDoxS2Ser5_CKDL230028616-1A_HCLMNDSX7_L4_1_fastqc.html"], "wildcards": {"fastq": "PP1H66KpDoxS2Ser5_CKDL230028616-1A_HCLMNDSX7_L4_1"}, "params": {"job_name": "PP1H66KpDoxS2Ser5_CKDL230028616-1A_HCLMNDSX7_L4_1_fastqc", "memory": 20, "out": "Temp4/fastqc"}, "log": ["Temp4/logs/PP1H66KpDoxS2Ser5_CKDL230028616-1A_HCLMNDSX7_L4_1_fastqc.out", "Temp4/logs/PP1H66KpDoxS2Ser5_CKDL230028616-1A_HCLMNDSX7_L4_1_fastqc.err"], "threads": 6, "resources": {"tmpdir": "/tmp"}, "jobid": 8, "cluster": {}}
 cd /beevol/home/erickson/Ben_pipelines/snakemake_pipelines && \
/cluster/software/modules-python/python/3.8.5/bin/python3 \
-m snakemake Temp4/fastqc/PP1H66KpDoxS2Ser5_CKDL230028616-1A_HCLMNDSX7_L4_1_fastqc.zip --snakefile /beevol/home/erickson/Ben_pipelines/snakemake_pipelines/pipelines/ChIPseq.snake \
--force --cores all --keep-target-files --keep-remote --max-inventory-time 0 \
--wait-for-files '/beevol/home/erickson/Ben_pipelines/snakemake_pipelines/.snakemake/tmp.0i3ompat' 'Temp4/fastqs/PP1H66KpDoxS2Ser5_CKDL230028616-1A_HCLMNDSX7_L4_1.fq.gz' --latency-wait 5 \
 --attempt 1 --force-use-threads --scheduler greedy \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
 --configfiles /beevol/home/erickson/Ben_pipelines/snakemake_pipelines/samples_ChIPseq3.yaml /beevol/home/erickson/Ben_pipelines/snakemake_pipelines/pipelines/ref/hg38.yaml  --allowed-rules fastqc --nocolor --notemp --no-hooks --nolock --scheduler-solver-path /cluster/software/modules-python/python/3.8.5/bin \
--mode 2  --local-groupid 8  --default-resources "tmpdir=system_tmpdir" 

