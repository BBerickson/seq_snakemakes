# in project folder get ChIPseq setup script

wget https://raw.githubusercontent.com/BBerickson/seq_snakemakes/main/ChIPseq/setup_pipeline.sh
# run script with one of the profiles: Bodhi | Alpine
bash setup_pipeline.sh Bodhi

# edit ChIPseq_samples.yaml with your sample information
# run script by submitting to Bodhi (lsf) or Alpine (slerm) 

