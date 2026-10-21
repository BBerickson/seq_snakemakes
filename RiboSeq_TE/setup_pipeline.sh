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

# RNA-seq arm -- Phase 2 always needs a completed, paired RNA-seq run, so
# this script pulls the same files RNAseq/setup_pipeline.sh itself pulls,
# bundled here rather than requiring a second, separate setup step.
RNA_SMK_FILE="RNAseq.smk"
RNA_SMK_FILE2="star.smk"
RNA_CONFIG_FILE="RNAseq_samples.yaml"
RNA_MTX_CONFIG_FILE="Stranded_matrix.yaml"

# TE join itself
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

# Create directory structure
mkdir -p workflow/{profiles,ref,rules,scripts}

# RNA-seq arm -- same files RNAseq/setup_pipeline.sh itself pulls
rsync -a "$EXTRACT_DIR/RiboSeq_TE/ReadMe.txt" .
rsync -a "$EXTRACT_DIR/workflow/$RNA_SMK_FILE" workflow/
rsync -a "$EXTRACT_DIR/workflow/$RNA_SMK_FILE2" workflow/
if [ ! -f "$RNA_CONFIG_FILE" ]; then
    rsync -a "$EXTRACT_DIR/workflow/configs/$RNA_CONFIG_FILE" .
fi
rsync -a "$EXTRACT_DIR/workflow/"*_matrix.smk workflow/
rsync -a "$EXTRACT_DIR/workflow/configs/$RNA_MTX_CONFIG_FILE" .
rsync -a "$EXTRACT_DIR/workflow/configs/"*_matrix.yaml workflow/
rsync -a "$EXTRACT_DIR/workflow/submit_scripts/$RNA_SUBMIT_GLOB1" . 2>/dev/null
rsync -a "$EXTRACT_DIR/workflow/submit_scripts/$RNA_SUBMIT_GLOB2" . 2>/dev/null

# TE join -- reads this completed RNA-seq run plus a separately-run,
# completed RiboSeq.smk/RiboSeq_UMI.smk run (its own setup_pipeline.sh, own
# project folder). NOT copied conditionally like $RNA_CONFIG_FILE above --
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
echo "!!!STAR is not set up for mixed genome alignment!!!"
echo "Files created:"
echo "  - workflow/$RNA_SMK_FILE"
echo "  - workflow/$RNA_SMK_FILE2"
echo "  - $RNA_CONFIG_FILE"
echo "  - workflow/$TE_SMK_FILE"
echo "  - $TE_CONFIG_FILE"
echo "  - workflow/rules/"
echo "  - workflow/ref/"
echo "  - workflow/profiles/"
echo "  - workflow/scripts/ (includes install_go_annotation_db.sh for GO enrichment)"
echo ""
echo "Next steps:"
echo "  1. Edit $RNA_CONFIG_FILE with your paired RNA-seq sample information, then"
echo "     run one of its submit scripts ($RNA_SUBMIT_GLOB1 / $RNA_SUBMIT_GLOB2)"
echo "  2. Separately, run RiboSeq.smk or RiboSeq_UMI.smk to completion (its own"
echo "     setup_pipeline.sh, own project folder) against the matched Ribo-seq fastqs"
echo "  3. Once both arms are done, edit $TE_CONFIG_FILE (RIBO_PROJ/RNA_PROJ paths,"
echo "     CONTRASTS by sample nickname), then run $TE_SUBMIT_GLOB -- or use"
echo "     $PIPELINE_SUBMIT_GLOB --ribo <ribo_project_dir> --rna . --te . to launch"
echo "     the Ribo-seq run + this RNA-seq run + the TE join in one command (the"
echo "     two independent arms run in parallel, TE deferred until both finish)"
echo "  4. For GO enrichment (off by default), see"
echo "     workflow/scripts/install_go_annotation_db.sh and"
echo "     ref/<your_genome>/go_annotation/README.txt"
