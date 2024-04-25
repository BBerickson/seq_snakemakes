#!/bin/sh
# properties = {"type": "single", "rule": "bamCoverage", "local": false, "input": ["Temp_NET/bams_sub/PP1H66KmDS2Ser53_hg38_230426.bam", "Temp_NET/bams_sub/PP1H66KpDS2Ser53_hg38_230426.bam", "Temp_NET/stats/Temp_NET_results.tsv"], "output": ["Temp_NET/bw/PP1H66KSer5-r3_hg38_230426_norm_3end_subsample_fw.bw", "Temp_NET/bw/PP1H66KSer5-r3_hg38_230426_norm_3end_subsample_rev.bw", "Temp_NET/bw/PP1H66KSer5-r3_bamCoverage_stats.txt"], "wildcards": {"group": "PP1H66KSer5-r3"}, "params": {"job_name": "PP1H66KSer5-r3_bamCoverage", "memory": 20, "args": " --binSize 1 --normalizeUsing None --skipNonCoveredRegions --Offset -1 ", "argsf": " reverse ", "argsr": " forward ", "scale": 1, "color": "0,0,0", "url_fw": "http://amc-sandbox.ucdenver.edu/User7/230426_Temp_NET/PP1H66KSer5-r3_hg38_230426_norm_3end_subsample_fw.bw", "url_rev": "http://amc-sandbox.ucdenver.edu/User7/230426_Temp_NET/PP1H66KSer5-r3_hg38_230426_norm_3end_subsample_rev.bw"}, "log": ["Temp_NET/logs/PP1H66KSer5-r3_bamCoverage.out", "Temp_NET/logs/PP1H66KSer5-r3_bamCoverage.err"], "threads": 12, "resources": {"tmpdir": "/tmp"}, "jobid": 50, "cluster": {}}
 cd /beevol/home/erickson/Ben_pipelines/snakemake_pipelines && \
/cluster/software/modules-python/python/3.8.5/bin/python3 \
-m snakemake Temp_NET/bw/PP1H66KSer5-r3_bamCoverage_stats.txt --snakefile /beevol/home/erickson/Ben_pipelines/snakemake_pipelines/pipelines/NETseq_spikein.snake \
--force --cores all --keep-target-files --keep-remote --max-inventory-time 0 \
--wait-for-files '/beevol/home/erickson/Ben_pipelines/snakemake_pipelines/.snakemake/tmp.19qcpfb0' 'Temp_NET/bams_sub/PP1H66KmDS2Ser53_hg38_230426.bam' 'Temp_NET/bams_sub/PP1H66KpDS2Ser53_hg38_230426.bam' 'Temp_NET/stats/Temp_NET_results.tsv' --latency-wait 5 \
 --attempt 1 --force-use-threads --scheduler greedy \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
 --configfiles /beevol/home/erickson/Ben_pipelines/snakemake_pipelines/samples_NETseq.yaml /beevol/home/erickson/Ben_pipelines/snakemake_pipelines/pipelines/ref/hg38_mm10.yaml  --allowed-rules bamCoverage --nocolor --notemp --no-hooks --nolock --scheduler-solver-path /cluster/software/modules-python/python/3.8.5/bin \
--mode 2  --local-groupid 50  --default-resources "tmpdir=system_tmpdir" 

