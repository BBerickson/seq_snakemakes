# in project folder sync pipeline


rsync -artuv ~/Ben_pipelines/snakemake_pipelines/pipelines .
cp -r ~/Ben_pipelines/snakemake_pipelines/misc/bw_to_matrix/* .