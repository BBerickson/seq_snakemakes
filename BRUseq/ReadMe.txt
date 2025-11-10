# in project folder sync BRUseq pipeline

wget https://raw.githubusercontent.com/BBerickson/seq_snakemakes/main/BRUseq/setup_pipeline.sh
# run script with one of the profiles: Bodhi | Alpine
bash setup_pipeline.sh Bodhi

# edit BRUseq_samples.yaml with your sample information
# run script by submitting to Bodhi (lsf) or Alpine (slerm) 

