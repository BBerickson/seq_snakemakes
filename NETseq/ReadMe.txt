# in project folder sync NETseq pipeline

rsync -artuv ~/Ben_pipelines/snakemake_pipelines/NETseq/* .
rsync -artuv ~/Ben_pipelines/snakemake_pipelines/pipelines .
rsync -artuv ~/Ben_pipelines/snakemake_pipelines/src .
