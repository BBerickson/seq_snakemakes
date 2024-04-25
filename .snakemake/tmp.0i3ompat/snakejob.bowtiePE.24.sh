#!/bin/sh
# properties = {"type": "single", "rule": "bowtiePE", "local": false, "input": ["Temp4/PP1H66KpDoxS2Ser5_R1.fastq.gz", "Temp4/PP1H66KpDoxS2Ser5_R2.fastq.gz"], "output": ["Temp4/PP1H66KpDoxS2Ser5_hg38.bam", "Temp4/PP1H66KpDoxS2Ser5_hg38.bam.bai", "Temp4/bams/PP1H66KpDoxS2Ser5_hg38_bowtie_stats.txt"], "wildcards": {"sample": "PP1H66KpDoxS2Ser5"}, "params": {"job_name": "PP1H66KpDoxS2Ser5_bowtiePE", "memory": 64, "idx": "/beevol/home/erickson/ref/hg38/hg38", "args": " --local --no-mixed --no-discordant ", "args2": " -bF4q30 ", "sortname": "Temp4/PP1H66KpDoxS2Ser5.temp"}, "log": ["Temp4/logs/PP1H66KpDoxS2Ser5_bowtie.out", "Temp4/logs/PP1H66KpDoxS2Ser5_bowtie.err"], "threads": 12, "resources": {"tmpdir": "/tmp"}, "jobid": 24, "cluster": {}}
 cd /beevol/home/erickson/Ben_pipelines/snakemake_pipelines && \
/cluster/software/modules-python/python/3.8.5/bin/python3 \
-m snakemake Temp4/bams/PP1H66KpDoxS2Ser5_hg38_bowtie_stats.txt --snakefile /beevol/home/erickson/Ben_pipelines/snakemake_pipelines/pipelines/ChIPseq.snake \
--force --cores all --keep-target-files --keep-remote --max-inventory-time 0 \
--wait-for-files '/beevol/home/erickson/Ben_pipelines/snakemake_pipelines/.snakemake/tmp.0i3ompat' 'Temp4/PP1H66KpDoxS2Ser5_R1.fastq.gz' 'Temp4/PP1H66KpDoxS2Ser5_R2.fastq.gz' --latency-wait 5 \
 --attempt 1 --force-use-threads --scheduler greedy \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
 --configfiles /beevol/home/erickson/Ben_pipelines/snakemake_pipelines/samples_ChIPseq3.yaml /beevol/home/erickson/Ben_pipelines/snakemake_pipelines/pipelines/ref/hg38.yaml  --allowed-rules bowtiePE --nocolor --notemp --no-hooks --nolock --scheduler-solver-path /cluster/software/modules-python/python/3.8.5/bin \
--mode 2  --local-groupid 24  --default-resources "tmpdir=system_tmpdir" 

