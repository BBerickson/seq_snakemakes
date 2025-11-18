# in project folder sync ChIPseq_duel_IP pipeline

wget https://raw.githubusercontent.com/BBerickson/seq_snakemakes/main/ChIPseq_duel_IP/setup_pipeline.sh
# run script with one of the profiles: Bodhi | Alpine
bash setup_pipeline.sh Bodhi

# edit ChIPseq_duel_IP_samples.yaml with your sample information
# run script by submitting to Bodhi (lsf) or Alpine (slerm) 

!!! Genome must point to ref yaml file with 2 genomes !!!
!!! see workflow/ref/hg38_polIII_SEQL_hg38.yaml for example !!!