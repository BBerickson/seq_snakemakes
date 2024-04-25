# in project folder sync NETseq pipeline with find pauses addon

rsync -artuv ~/Ben_pipelines/snakemake_pipelines/NETseq_UMI/* .
rsync -artuv ~/Ben_pipelines/snakemake_pipelines/NETseq_UMI_findPauses/* .
rsync -artuv ~/Ben_pipelines/snakemake_pipelines/pipelines .
rsync -artuv ~/Ben_pipelines/snakemake_pipelines/src .
