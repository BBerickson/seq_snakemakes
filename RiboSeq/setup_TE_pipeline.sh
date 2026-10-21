#!/bin/bash
# setup_TE_pipeline.sh - Download specific files from GitHub
# Phase 2 (optional): creates sibling rnaseq/ and te/ subfolders next to this
# already-set-up Ribo-seq project folder, pulling the paired RNA-seq arm
# (RNAseq.smk/star.smk) into rnaseq/ and the TE join (RiboSeq_TE.smk) into
# te/. NOT the same folder as this Ribo-seq run -- Snakemake locks its whole
# working directory for a run, so two snakemake invocations sharing one CWD
# (Ribo-seq + RNA-seq) collide the moment both run at once, which is
# exactly what run_RiboSeq_TE_pipeline_bodhi.sh does (submits both in
# parallel). Run this after setup_pipeline.sh (Phase 1) has already been
# run here. See ReadMe.txt's Phase 2 section for the full walkthrough.

GITHUB_USER="BBerickson"
REPO="seq_snakemakes"
BRANCH="main"

show_usage() {
    echo "Usage: $0 <profile_type>"
    echo ""
    echo "Available profiles:"
    echo "  Bodhi      - SLURM"
    echo ""
    echo "Example: $0 Bodhi"
    exit 1
}

if [ $# -eq 0 ]; then
    show_usage
fi

# This script only adds Phase 2's own files (RNA-seq arm + TE join) --
# RiboSeq_samples.yaml/ReadMe.txt/run_RiboSeq_bodhi.sh are Phase 1's
# responsibility (setup_pipeline.sh), so check they're already here rather
# than silently leaving an incomplete folder if someone runs this first.
RIBO_CONFIG_FILE="RiboSeq_samples.yaml"
if [ ! -f "$RIBO_CONFIG_FILE" ]; then
    echo "Error: $RIBO_CONFIG_FILE not found in this folder." >&2
    echo "Run setup_pipeline.sh (Phase 1) here first: bash setup_pipeline.sh Bodhi" >&2
    exit 1
fi

# RNA-seq arm (rnaseq/) -- same files RNAseq/setup_pipeline.sh itself pulls
RNA_SMK_FILE="RNAseq.smk"
RNA_SMK_FILE2="star.smk"
RNA_CONFIG_FILE="RNAseq_samples.yaml"
RNA_MTX_CONFIG_FILE="Stranded_matrix.yaml"

# TE join (te/)
TE_SMK_FILE="RiboSeq_TE.smk"
TE_CONFIG_FILE="RiboSeq_TE_samples.yaml"

PIPELINE_TYPE="$1"

case "$PIPELINE_TYPE" in
    Bodhi)
        RNA_SUBMIT_GLOB1="run_RNAseq_Hisat_bodhi.sh"
        RNA_SUBMIT_GLOB2="run_RNAseq_STAR_bodhi.sh"
        TE_SUBMIT_GLOB="run_RiboSeq_TE_bodhi.sh"
        PIPELINE_SUBMIT_GLOB="run_RiboSeq_TE_pipeline_bodhi.sh"
        PROFILES="Bodhi"
        ;;
    *)
        echo "Error: Unknown pipeline type '$PIPELINE_TYPE'"
        show_usage
        ;;
esac


# Download entire repo as tarball, extract only what we need
TEMP_DIR=$(mktemp -d)
echo "Downloading pipeline repository..."

wget -q "https://github.com/$GITHUB_USER/$REPO/archive/refs/heads/$BRANCH.tar.gz" -O "$TEMP_DIR/repo.tar.gz"

if [ $? -ne 0 ]; then
    echo "Error: Failed to download repository"
    rm -rf "$TEMP_DIR"
    exit 1
fi

echo "Extracting files..."
tar -xzf "$TEMP_DIR/repo.tar.gz" -C "$TEMP_DIR"
if [ $? -ne 0 ]; then
    echo "Error: Failed to extract repository"
    rm -rf "$TEMP_DIR"
    exit 1
fi

EXTRACT_DIR="$TEMP_DIR/$REPO-$BRANCH"

# Each arm gets its own complete, self-contained workflow/ tree (same
# convention every other pipeline folder in this repo already uses) -- not
# shared with the Ribo-seq run already set up flat in this folder, and not
# shared between rnaseq/ and te/ either, so each one's Snakemake invocation
# locks only its own subfolder.
mkdir -p rnaseq/workflow/{profiles,ref,rules,scripts}
mkdir -p te/workflow/{profiles,ref,rules,scripts}

# RNA-seq arm
rsync -a "$EXTRACT_DIR/workflow/$RNA_SMK_FILE" rnaseq/workflow/
rsync -a "$EXTRACT_DIR/workflow/$RNA_SMK_FILE2" rnaseq/workflow/
if [ ! -f "rnaseq/$RNA_CONFIG_FILE" ]; then
    rsync -a "$EXTRACT_DIR/workflow/configs/$RNA_CONFIG_FILE" rnaseq/
fi
rsync -a "$EXTRACT_DIR/workflow/"*_matrix.smk rnaseq/workflow/
rsync -a "$EXTRACT_DIR/workflow/configs/$RNA_MTX_CONFIG_FILE" rnaseq/
rsync -a "$EXTRACT_DIR/workflow/configs/"*_matrix.yaml rnaseq/workflow/
rsync -a "$EXTRACT_DIR/workflow/submit_scripts/$RNA_SUBMIT_GLOB1" rnaseq/ 2>/dev/null
rsync -a "$EXTRACT_DIR/workflow/submit_scripts/$RNA_SUBMIT_GLOB2" rnaseq/ 2>/dev/null
rsync -a "$EXTRACT_DIR/workflow/rules" rnaseq/workflow/
rsync -a "$EXTRACT_DIR/workflow/configs/ref" rnaseq/workflow/
rsync -a "$EXTRACT_DIR/workflow/profiles/$PROFILES" rnaseq/workflow/profiles/
rsync -a "$EXTRACT_DIR/workflow/scripts" rnaseq/workflow/
rsync -a "$EXTRACT_DIR/workflow/profiles/$PROFILES/ref/"* rnaseq/workflow/ref/

# TE join -- reads RIBO_PROJ (this Ribo-seq run, one level up) and RNA_PROJ
# (../rnaseq/<proj>, once run). NOT copied conditionally like
# $RNA_CONFIG_FILE above -- RIBO_PROJ/RNA_PROJ/CONTRASTS always need
# hand-editing per project anyway.
rsync -a "$EXTRACT_DIR/workflow/$TE_SMK_FILE" te/workflow/
if [ ! -f "te/$TE_CONFIG_FILE" ]; then
    rsync -a "$EXTRACT_DIR/workflow/configs/$TE_CONFIG_FILE" te/
fi
rsync -a "$EXTRACT_DIR/workflow/submit_scripts/$TE_SUBMIT_GLOB" te/ 2>/dev/null
rsync -a "$EXTRACT_DIR/workflow/rules" te/workflow/
rsync -a "$EXTRACT_DIR/workflow/configs/ref" te/workflow/
rsync -a "$EXTRACT_DIR/workflow/profiles/$PROFILES" te/workflow/profiles/
rsync -a "$EXTRACT_DIR/workflow/scripts" te/workflow/
rsync -a "$EXTRACT_DIR/workflow/profiles/$PROFILES/ref/"* te/workflow/ref/

# Orchestrator script (optional) -- lives at the TOP level, alongside this
# already-set-up Ribo-seq run, since it launches all three arms by
# directory path (--ribo . --rna ./rnaseq --te ./te)
rsync -a "$EXTRACT_DIR/workflow/submit_scripts/$PIPELINE_SUBMIT_GLOB" . 2>/dev/null
# Cleanup
rm -rf "$TEMP_DIR"

echo ""
echo "✓ $PIPELINE_TYPE TE pipeline setup complete!"
echo ""
echo "!!!STAR is not set up for mixed genome alignment!!!"
echo "Files created:"
echo "  - rnaseq/workflow/$RNA_SMK_FILE"
echo "  - rnaseq/workflow/$RNA_SMK_FILE2"
echo "  - rnaseq/$RNA_CONFIG_FILE"
echo "  - te/workflow/$TE_SMK_FILE"
echo "  - te/$TE_CONFIG_FILE"
echo "  - $PIPELINE_SUBMIT_GLOB"
echo ""
echo "Next steps:"
echo "  1. cd rnaseq, edit $RNA_CONFIG_FILE with your paired RNA-seq sample"
echo "     information, then run one of its submit scripts ($RNA_SUBMIT_GLOB1 /"
echo "     $RNA_SUBMIT_GLOB2)"
echo "  2. Once both this Ribo-seq run and that RNA-seq run are done, cd te and"
echo "     edit $TE_CONFIG_FILE (RIBO_PROJ -- one level up at ../<this Ribo-seq"
echo "     run's PROJ> -- and RNA_PROJ -- ../rnaseq/<its PROJ> --, plus CONTRASTS"
echo "     by sample nickname), then run $TE_SUBMIT_GLOB from inside te/ -- or"
echo "     from this top folder use"
echo "     ./$PIPELINE_SUBMIT_GLOB --ribo . --rna ./rnaseq --te ./te to launch the"
echo "     Ribo-seq run + the RNA-seq run + the TE join in one command (the two"
echo "     independent arms run in parallel in their own subfolders, TE deferred"
echo "     until both finish)"
echo "  3. For GO enrichment (off by default), see"
echo "     te/workflow/scripts/install_go_annotation_db.sh and"
echo "     ref/<your_genome>/go_annotation/README.txt"
