# in project folder sync NETseq pipeline

rsync -artuv /beevol/home/erickson/Ben_pipelines/snakemake_pipelines/NETseq/* .
wget https://raw.githubusercontent.com/BBerickson/seq_snakemakes/main/NETseq/setup_pipeline.sh
# run script with one of the profiles: Bodhi | Alpine
bash setup_pipeline.sh Bodhi

# edit NETseq_samples.yaml with your sample information
# run script by submitting to Bodhi (lsf) or Alpine (slerm) 
