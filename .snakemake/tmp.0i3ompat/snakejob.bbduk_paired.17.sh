#!/bin/sh
# properties = {"type": "single", "rule": "bbduk_paired", "local": false, "input": ["Temp4/PP1H66KpDS2Ser53_R1_clumpify.fastq.gz", "Temp4/PP1H66KpDS2Ser53_R2_clumpify.fastq.gz"], "output": ["Temp4/PP1H66KpDS2Ser53_R1.fastq.gz", "Temp4/PP1H66KpDS2Ser53_R2.fastq.gz", "Temp4/bbduk/PP1H66KpDS2Ser53_bbduk.log", "Temp4/bbduk/PP1H66KpDS2Ser53_bbduk_KTrimmed.log"], "wildcards": {"sample": "PP1H66KpDS2Ser53"}, "params": {"job_name": "PP1H66KpDS2Ser53_bbduk", "memory": 120, "args": " overwrite=t ftm=5 minlen=50 ktrim=r k=23 hdist=1 mink=11 tpe tbo ecco=t ", "bar": ""}, "log": ["Temp4/logs/PP1H66KpDS2Ser53_bbduk.out", "Temp4/logs/PP1H66KpDS2Ser53_bbduk.err"], "threads": 12, "resources": {"tmpdir": "/tmp"}, "jobid": 17, "cluster": {}}
 cd /beevol/home/erickson/Ben_pipelines/snakemake_pipelines && \
/cluster/software/modules-python/python/3.8.5/bin/python3 \
-m snakemake Temp4/bbduk/PP1H66KpDS2Ser53_bbduk.log --snakefile /beevol/home/erickson/Ben_pipelines/snakemake_pipelines/pipelines/ChIPseq.snake \
--force --cores all --keep-target-files --keep-remote --max-inventory-time 0 \
--wait-for-files '/beevol/home/erickson/Ben_pipelines/snakemake_pipelines/.snakemake/tmp.0i3ompat' 'Temp4/PP1H66KpDS2Ser53_R1_clumpify.fastq.gz' 'Temp4/PP1H66KpDS2Ser53_R2_clumpify.fastq.gz' --latency-wait 5 \
 --attempt 1 --force-use-threads --scheduler greedy \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
 --configfiles /beevol/home/erickson/Ben_pipelines/snakemake_pipelines/samples_ChIPseq3.yaml /beevol/home/erickson/Ben_pipelines/snakemake_pipelines/pipelines/ref/hg38.yaml  --allowed-rules bbduk_paired --nocolor --notemp --no-hooks --nolock --scheduler-solver-path /cluster/software/modules-python/python/3.8.5/bin \
--mode 2  --local-groupid 17  --default-resources "tmpdir=system_tmpdir" 

