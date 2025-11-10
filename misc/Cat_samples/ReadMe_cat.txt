# in project folder sync pipeline

wget https://raw.githubusercontent.com/BBerickson/seq_snakemakes/main/misc/Cat_samples/setup_pipeline.sh
# run script with one of the profiles: Bodhi | Alpine
bash setup_pipeline.sh Bodhi

# edit samples_cat.yaml with your sample information
# run script by submitting to Bodhi (lsf) or Alpine (slerm) 

