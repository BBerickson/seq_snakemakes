#!/bin/bash
# setup_TE_pipeline.sh - Download specific files from GitHub
# Alternative to setup_pipeline.sh for when you already know you want
# translational efficiency (Phase 2): sets up all three arms from scratch,
# in one call, in a fresh/empty project folder -- riboseq/ (Phase 1) and
# rnaseq/ each get their own complete, self-contained workflow/ tree; the
# TE join (RiboSeq_TE.smk) lives at this TOP level. Separate subfolders for
# riboseq/rnaseq, not this same folder -- Snakemake locks its whole working
# directory for a run, so two snakemake invocations sharing one CWD would
# collide the moment both run at once, which is exactly what
# run_RiboSeq_TE_pipeline_bodhi.sh does (submits both in parallel). See
# ReadMe.txt's Phase 2 section for the full walkthrough.
#
# If you already have a completed Ribo-seq run elsewhere (set up via this
# folder's own setup_pipeline.sh) and just want to add TE later, don't run
# this script -- just point RiboSeq_TE_samples.yaml's RIBO_PROJ directly at
# that existing run instead (any path works, nothing here requires it).

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

# Ribo-seq arm (riboseq/) -- same files this folder's own setup_pipeline.sh
# pulls, just targeting a subfolder instead of landing flat here.
RIBO_SMK_FILE="RiboSeq.smk"
RIBO_CONFIG_FILE="RiboSeq_samples.yaml"

# RNA-seq arm (rnaseq/) -- same files RNAseq/setup_pipeline.sh itself pulls
RNA_SMK_FILE="RNAseq.smk"
RNA_SMK_FILE2="star.smk"
RNA_CONFIG_FILE="RNAseq_samples.yaml"
RNA_MTX_CONFIG_FILE="Stranded_matrix.yaml"

# TE join -- this top level
TE_SMK_FILE="RiboSeq_TE.smk"
TE_CONFIG_FILE="RiboSeq_TE_samples.yaml"

PIPELINE_TYPE="$1"

case "$PIPELINE_TYPE" in
    Bodhi)
        RIBO_SUBMIT_GLOB="run_RiboSeq_bodhi.sh"
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
# convention every other pipeline folder in this repo already uses) --
# riboseq/ and rnaseq/ lock only their own subfolder when running, and
# this top level (TE) is never run concurrently with either (it's deferred
# until both finish), so it doesn't need its own subfolder.
mkdir -p riboseq/workflow/{profiles,ref,rules,scripts}
mkdir -p rnaseq/workflow/{profiles,ref,rules,scripts}
mkdir -p workflow/{profiles,ref,rules,scripts}

# Ribo-seq arm
rsync -a "$EXTRACT_DIR/RiboSeq/ReadMe.txt" riboseq/
rsync -a "$EXTRACT_DIR/workflow/$RIBO_SMK_FILE" riboseq/workflow/
if [ ! -f "riboseq/$RIBO_CONFIG_FILE" ]; then
    rsync -a "$EXTRACT_DIR/workflow/configs/$RIBO_CONFIG_FILE" riboseq/
fi
rsync -a "$EXTRACT_DIR/workflow/submit_scripts/$RIBO_SUBMIT_GLOB" riboseq/ 2>/dev/null
rsync -a "$EXTRACT_DIR/workflow/rules" riboseq/workflow/
rsync -a "$EXTRACT_DIR/workflow/configs/ref" riboseq/workflow/
rsync -a "$EXTRACT_DIR/workflow/profiles/$PROFILES" riboseq/workflow/profiles/
rsync -a "$EXTRACT_DIR/workflow/scripts" riboseq/workflow/
rsync -a "$EXTRACT_DIR/workflow/profiles/$PROFILES/ref/"* riboseq/workflow/ref/

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

# TE join -- this top level. RIBO_PROJ/RNA_PROJ/CONTRASTS always need
# hand-editing per project anyway, so TE_CONFIG_FILE isn't copied
# conditionally like the other two configs above.
rsync -a "$EXTRACT_DIR/workflow/$TE_SMK_FILE" workflow/
if [ ! -f "$TE_CONFIG_FILE" ]; then
    rsync -a "$EXTRACT_DIR/workflow/configs/$TE_CONFIG_FILE" .
fi
rsync -a "$EXTRACT_DIR/workflow/submit_scripts/$TE_SUBMIT_GLOB" . 2>/dev/null
rsync -a "$EXTRACT_DIR/workflow/submit_scripts/$PIPELINE_SUBMIT_GLOB" . 2>/dev/null
rsync -a "$EXTRACT_DIR/workflow/rules" workflow/
rsync -a "$EXTRACT_DIR/workflow/configs/ref" workflow/
rsync -a "$EXTRACT_DIR/workflow/profiles/$PROFILES" workflow/profiles/
rsync -a "$EXTRACT_DIR/workflow/scripts" workflow/
rsync -a "$EXTRACT_DIR/workflow/profiles/$PROFILES/ref/"* workflow/ref/
# Cleanup
rm -rf "$TEMP_DIR"

echo ""
echo "✓ $PIPELINE_TYPE TE pipeline setup complete!"
echo ""
echo "!!!STAR is not set up for mixed genome alignment!!!"
echo "Files created:"
echo "  - riboseq/workflow/$RIBO_SMK_FILE"
echo "  - riboseq/$RIBO_CONFIG_FILE"
echo "  - rnaseq/workflow/$RNA_SMK_FILE"
echo "  - rnaseq/workflow/$RNA_SMK_FILE2"
echo "  - rnaseq/$RNA_CONFIG_FILE"
echo "  - workflow/$TE_SMK_FILE"
echo "  - $TE_CONFIG_FILE"
echo "  - $PIPELINE_SUBMIT_GLOB"
echo ""
echo "Next steps:"
echo "  1. cd riboseq, edit $RIBO_CONFIG_FILE with your sample information, and"
echo "     run $RIBO_SUBMIT_GLOB"
echo "  2. cd rnaseq, edit $RNA_CONFIG_FILE with your paired RNA-seq sample"
echo "     information, then run one of its submit scripts ($RNA_SUBMIT_GLOB1 /"
echo "     $RNA_SUBMIT_GLOB2)"
echo "  3. Once both arms are done, edit $TE_CONFIG_FILE here (RIBO_PROJ --"
echo "     riboseq/<its PROJ> -- and RNA_PROJ -- rnaseq/<its PROJ> --, plus"
echo "     CONTRASTS by sample nickname), then run $TE_SUBMIT_GLOB -- or use"
echo "     ./$PIPELINE_SUBMIT_GLOB --ribo ./riboseq --rna ./rnaseq --te . to"
echo "     launch all three stages in one command (the two independent arms run"
echo "     in parallel in their own subfolders, TE deferred until both finish)"
echo "  4. For GO enrichment (off by default), see"
echo "     workflow/scripts/install_go_annotation_db.sh and"
echo "     ref/<your_genome>/go_annotation/README.txt"
