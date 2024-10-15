# in project folder sync NETseq pipeline

rsync -artuv ~/Ben_pipelines/snakemake_pipelines/NETseq_UMI/* .
rsync -artuv ~/Ben_pipelines/snakemake_pipelines/pipelines .
rsync -artuv ~/Ben_pipelines/snakemake_pipelines/src .
cp ~/Ben_pipelines/snakemake_pipelines/NETseq_UMI/src/Rmds/qc-template.Rmd src/Rmds/