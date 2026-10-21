# in project folder sync the RNA-seq + translational-efficiency (TE) pipeline
# This pulls BOTH the RNA-seq arm (RNAseq.smk/star.smk) and the TE join
# (RiboSeq_TE.smk) into one project folder -- Phase 2 always needs a
# completed, paired RNA-seq run, so this one script sets both up together.
# Run RiboSeq.smk or RiboSeq_UMI.smk SEPARATELY (own setup_pipeline.sh, own
# project folder) for the Ribo-seq arm -- see RiboSeq/ReadMe.txt or
# RiboSeq_UMI/ReadMe.txt. This folder does not set up Ribo-seq itself.

wget https://raw.githubusercontent.com/BBerickson/seq_snakemakes/main/RiboSeq_TE/setup_pipeline.sh
# run script with one of the profiles: Bodhi
bash setup_pipeline.sh Bodhi

# edit RNAseq_samples.yaml with your paired RNA-seq sample information and
# run one of run_RNAseq_Hisat_bodhi.sh / run_RNAseq_STAR_bodhi.sh


### Pipeline Overview -- Phase 2: Translational Efficiency ###

Joins a COMPLETED Ribo-seq run (RiboSeq.smk/RiboSeq_UMI.smk, own project
folder -- see RiboSeq/ReadMe.txt or RiboSeq_UMI/ReadMe.txt) with the
COMPLETED, PAIRED RNA-seq run set up by this same folder, into anota2seq
(primary, needs >=3 replicates/group) + xtail (secondary, no replicate
minimum) translational-efficiency analysis. Does no alignment/quantification
of its own -- only reads both runs' completed featureCounts outputs. Works
identically whether the Ribo-seq side used UMIs or not -- RiboSeq_TE.smk
doesn't care.

Once BOTH this RNA-seq run and the separately-run Ribo-seq run are done,
edit RiboSeq_TE_samples.yaml (RIBO_PROJ/RNA_PROJ paths, CONTRASTS by sample
nickname, CMD_PARAMS.use_psite_counts if the Ribo-seq run used
enable_psite_cds_counts: true), then either:
  - run run_RiboSeq_TE_bodhi.sh directly, once both arms are already done, or
  - run run_RiboSeq_TE_pipeline_bodhi.sh --ribo <ribo_project_dir> --rna . --te .
    to launch the Ribo-seq run + this RNA-seq run + the TE join in one
    command -- the two independent arms run in parallel, TE is deferred via
    a SLURM job dependency until both finish successfully.

GO enrichment on TE-changed genes (CMD_PARAMS.run_go_enrichment, off by
default) needs a species-specific org.*.eg.db annotation package installed
OUTSIDE the container first -- see workflow/scripts/install_go_annotation_db.sh
and ref/<your_genome>/go_annotation/README.txt (hg38/human already set up
in this repo's own test project; same script installs any other organism,
e.g. mouse, into the same shared location, no container rebuild needed).
