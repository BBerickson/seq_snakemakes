# in project folder sync rMATS add on pipeline

wget https://raw.githubusercontent.com/BBerickson/seq_snakemakes/main/misc/RNAseq_rMATS/setup_rMATS_pipeline.sh
# run script with one of the profiles: Bodhi | Alpine
bash setup_rMATS_pipeline.sh Bodhi

# edit samples_rMATS.yaml and run_rMATS.sh with your sample information
# run script by submitting to Bodhi (lsf) or Alpine (slerm) after running main alignment script

# the additinal scripts are not set up for Alpine yet

