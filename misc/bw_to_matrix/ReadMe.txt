# in project folder sync pipeline

wget https://raw.githubusercontent.com/BBerickson/seq_snakemakes/main/misc/bw_to_matrix/setup_pipeline.sh
# run script with one of the profiles: Bodhi | Alpine
bash setup_pipeline.sh Bodhi

# edit BW_samples.yaml with your sample information
# run script by submitting to Bodhi (lsf) or Alpine (slerm) 

### Pipeline Overview ###
1. Coverage Track Generation

2a bigWig URLs (Bodhi option) — Copy bigWig to sandbox and generate UCSC URLs

3. Matrix Computation

3a deepTools computeMatrix — Create sense/antisense stranded matrix files
3b BW URLs (Bodhi option) — Copy matrix to sandbox with URLs for bentools
