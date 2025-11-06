# in project folder sync ChIPseq pipeline

rsync -artuv /beevol/home/erickson/Ben_pipelines/snakemake_pipelines/ChIPseq/* .
rsync -artuv /beevol/home/erickson/Ben_pipelines/snakemake_pipelines/pipelines .

######
wget https://github.com/BBerickson/seq_snakemakes/blob/main/ChIPseq/setup_pipeline.sh
bash setup_pipeline.sh
