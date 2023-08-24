#!/bin/sh
# properties = {"type": "single", "rule": "create_beds", "local": false, "input": ["Xrn2/bams/293flpin_Xrn2_MT_R2_hg38_220215_220218_dedup_filtered_gene.bam"], "output": ["Xrn2/beds/293flpin_Xrn2_MT_R2_hg38_220215_220218.3bed.gz", "Xrn2/beds/293flpin_Xrn2_MT_R2_hg38_220215_220218_3bed.txt"], "wildcards": {"sample": "293flpin_Xrn2_MT_R2"}, "params": {"job_name": "293flpin_Xrn2_MT_R2_create_beds", "memory": 96, "genes": "/beevol/home/erickson/ref/hg38/ref/hg38_4_-5kb_+5kb.bed", "mask": "/beevol/home/erickson/ref/hg38/ref/hg38_snoRNAs.bed"}, "log": ["Xrn2/logs/293flpin_Xrn2_MT_R2_beds.out", "Xrn2/logs/293flpin_Xrn2_MT_R2_beds.err"], "threads": 12, "resources": {"tmpdir": "/tmp"}, "jobid": 9, "cluster": {}}
 cd /beevol/home/erickson/Ben_pipelines/snakemake_pipelines/NETseq_UMI_findPauses && \
/cluster/software/modules-python/python/3.8.5/bin/python3 \
-m snakemake Xrn2/beds/293flpin_Xrn2_MT_R2_hg38_220215_220218.3bed.gz --snakefile /beevol/home/erickson/Ben_pipelines/snakemake_pipelines/NETseq_UMI_findPauses/pipelines/mNETseq_pauses.snake \
--force --cores all --keep-target-files --keep-remote --max-inventory-time 0 \
--wait-for-files '/beevol/home/erickson/Ben_pipelines/snakemake_pipelines/NETseq_UMI_findPauses/.snakemake/tmp.wl2jkw3l' 'Xrn2/bams/293flpin_Xrn2_MT_R2_hg38_220215_220218_dedup_filtered_gene.bam' --latency-wait 5 \
 --attempt 1 --force-use-threads --scheduler greedy \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
 --configfiles /beevol/home/erickson/Ben_pipelines/snakemake_pipelines/NETseq_UMI_findPauses/samples.yaml /beevol/home/erickson/Ben_pipelines/snakemake_pipelines/NETseq_UMI_findPauses/pipelines/ref/hg38.yaml  --allowed-rules create_beds --nocolor --notemp --no-hooks --nolock --scheduler-solver-path /cluster/software/modules-python/python/3.8.5/bin \
--mode 2  --local-groupid 9  --default-resources "tmpdir=system_tmpdir" 

