#!/bin/bash
# setup_pipeline.sh - Download specific files from GitHub

GITHUB_USER="BBerickson"
REPO="seq_snakemakes"
BRANCH="main"

show_usage() {
    echo "Usage: $0 <profile_type>"
    echo ""
    echo "Available profiles:"
    echo "  Bodhi      - LSF"
    echo "  Alpine     - SLURM"
    echo ""
    echo "Example: $0 Bodhi"
    exit 1
}

if [ $# -eq 0 ]; then
    show_usage
fi

CONFIG_FILE="BW_samples.yaml"

PIPELINE_TYPE="$1"

case "$PIPELINE_TYPE" in
    Bodhi)
        SUBMIT_GLOB="run_BW_bodhi.sh"
        PROFILES="Bodhi"
        ;;
    Alpine)
        SUBMIT_GLOB="run_BW_alpine.sh"
        PROFILES="Alpine"
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
mkdir -p workflow/{profiles,ref,rules,Rmds,scripts}

# Copy pipeline-specific files
rsync -a "$EXTRACT_DIR/misc/bw_to_matrix/ReadMe.txt" .
if [ ! -f "$CONFIG_FILE" ]; then
    rsync -a "$EXTRACT_DIR/workflow/configs/$CONFIG_FILE" .
fi
rsync -a "$EXTRACT_DIR/workflow/"*_matrix.smk workflow/
rsync -a "$EXTRACT_DIR/workflow/configs/"*_matrix.yaml .

rsync -a "$EXTRACT_DIR/workflow/submit_scripts/$SUBMIT_GLOB" . 2>/dev/null

# Copy shared directories
rsync -a "$EXTRACT_DIR/workflow/rules" workflow/
rsync -a "$EXTRACT_DIR/workflow/configs/ref" workflow/
rsync -a "$EXTRACT_DIR/workflow/profiles/$PROFILES" workflow/profiles/
rsync -a "$EXTRACT_DIR/workflow/Rmds" workflow/
rsync -a "$EXTRACT_DIR/workflow/scripts" workflow/

# Copy profile files
rsync -a "$EXTRACT_DIR/workflow/profiles/$PROFILES/ref/"* workflow/ref/
# Cleanup
rm -rf "$TEMP_DIR"

echo ""
echo "✓ $PIPELINE_TYPE pipeline setup complete!"
echo ""
echo "Files created:"
echo "  - $CONFIG_FILE"
echo "  - workflow/rules/"
echo "  - workflow/ref/"
echo "  - workflow/profiles/"
echo "  - workflow/Rmd/"
echo "  - workflow/scripts/"
echo ""
echo "Next steps:"
echo "  1. Edit $CONFIG_FILE with your sample information"
echo "  2. Edit and run the submit script"


