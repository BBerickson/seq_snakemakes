#!/bin/bash
# setup_pipeline.sh - Download specific files from GitHub

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

SMK_FILE="RiboSeq_UMI.smk"
CONFIG_FILE="RiboSeq_UMI_samples.yaml"
TE_SMK_FILE="RiboSeq_TE.smk"
TE_CONFIG_FILE="RiboSeq_TE_samples.yaml"

PIPELINE_TYPE="$1"

case "$PIPELINE_TYPE" in
    Bodhi)
        SUBMIT_GLOB="run_RiboSeq_UMI_bodhi.sh"
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

# Create directory structure
mkdir -p workflow/{profiles,ref,rules,scripts}

# Copy pipeline-specific files -- Phase 1 (alignment + QC)
rsync -a "$EXTRACT_DIR/RiboSeq_UMI/ReadMe.txt" .
rsync -a "$EXTRACT_DIR/workflow/$SMK_FILE" workflow/
if [ ! -f "$CONFIG_FILE" ]; then
    rsync -a "$EXTRACT_DIR/workflow/configs/$CONFIG_FILE" .
fi
rsync -a "$EXTRACT_DIR/workflow/submit_scripts/$SUBMIT_GLOB" . 2>/dev/null

# Phase 2 (translational-efficiency join against a completed RNAseq.smk run)
# -- run as a second, separate snakemake workflow against the same project
# folder, once both this Ribo-seq run and the paired RNA-seq run are done
# (see ReadMe.txt). NOT copied conditionally like $CONFIG_FILE above --
# RIBO_PROJ/RNA_PROJ/CONTRASTS always need hand-editing per project anyway.
rsync -a "$EXTRACT_DIR/workflow/$TE_SMK_FILE" workflow/
if [ ! -f "$TE_CONFIG_FILE" ]; then
    rsync -a "$EXTRACT_DIR/workflow/configs/$TE_CONFIG_FILE" .
fi
rsync -a "$EXTRACT_DIR/workflow/submit_scripts/$TE_SUBMIT_GLOB" . 2>/dev/null
rsync -a "$EXTRACT_DIR/workflow/submit_scripts/$PIPELINE_SUBMIT_GLOB" . 2>/dev/null

# Copy shared directories
rsync -a "$EXTRACT_DIR/workflow/rules" workflow/
rsync -a "$EXTRACT_DIR/workflow/configs/ref" workflow/
rsync -a "$EXTRACT_DIR/workflow/profiles/$PROFILES" workflow/profiles/
rsync -a "$EXTRACT_DIR/workflow/scripts" workflow/

# Copy profile files
rsync -a "$EXTRACT_DIR/workflow/profiles/$PROFILES/ref/"* workflow/ref/
# Cleanup
rm -rf "$TEMP_DIR"

echo ""
echo "✓ $PIPELINE_TYPE pipeline setup complete!"
echo ""
echo "Files created:"
echo "  - workflow/$SMK_FILE"
echo "  - $CONFIG_FILE"
echo "  - workflow/$TE_SMK_FILE"
echo "  - $TE_CONFIG_FILE"
echo "  - workflow/rules/"
echo "  - workflow/ref/"
echo "  - workflow/profiles/"
echo "  - workflow/scripts/ (includes install_go_annotation_db.sh for GO enrichment)"
echo ""
echo "Next steps:"
echo "  1. Edit $CONFIG_FILE with your sample information and run its submit script"
echo "  2. Separately, run RNAseq.smk (its own setup_pipeline.sh, in its own project"
echo "     folder) against the matched RNA-seq fastqs, if you want translational"
echo "     efficiency (Phase 2) -- not needed for Ribo-seq alignment/QC alone"
echo "  3. Edit $TE_CONFIG_FILE to point at both completed runs, then run its submit"
echo "     script -- or use $PIPELINE_SUBMIT_GLOB to launch all three stages in one"
echo "     command (Ribo-seq + RNA-seq in parallel, TE deferred until both finish)"
echo "  4. For GO enrichment (off by default), see"
echo "     workflow/scripts/install_go_annotation_db.sh and"
echo "     ref/<your_genome>/go_annotation/README.txt"
