#!/usr/bin/env bash
# Installs a Bioconductor org.*.eg.db annotation package into a shared R
# library directory, OUTSIDE the Ribo-seq Singularity container -- see
# ref/hg38_p14/go_annotation/README.txt for why (container stays
# genome-agnostic + stable; adding a new organism never needs a rebuild).
#
# This is a one-time SETUP step and needs internet access from wherever you
# run it (confirmed both the Bodhi login node and its compute nodes could
# reach bioconductor.org as of 2026-10-16) -- the actual RiboSeq_TE.smk
# pipeline never needs network access at runtime, it only reads whatever
# this script already installed.
#
# Usage:
#   ./install_go_annotation_db.sh <org_db_package> <lib_dir> [container]
#
# Example (human, matching this repo's existing hg38_p14 setup):
#   ./install_go_annotation_db.sh org.Hs.eg.db /beevol/home/erickson/ref/hg38_p14/go_annotation
#
# Example (mouse, same shared directory -- org.*.eg.db is per-organism, not
# per-genome-build, so one directory can hold every species you need):
#   ./install_go_annotation_db.sh org.Mm.eg.db /beevol/home/erickson/ref/hg38_p14/go_annotation
#
# Other common package names: org.Rn.eg.db (rat), org.Dm.eg.db (fly),
# org.Dr.eg.db (zebrafish), org.Sc.sgd.db (yeast). Full list:
# https://bioconductor.org/packages/release/BiocViews.html#___OrgDb

set -o nounset -o pipefail -o errexit

CONTAINER="${3:-docker://bberickson/riboseq:1.0.0}"

if [[ $# -lt 2 ]]; then
    echo "Usage: $0 <org_db_package> <lib_dir> [container]" >&2
    echo "Example: $0 org.Hs.eg.db /beevol/home/erickson/ref/hg38_p14/go_annotation" >&2
    exit 1
fi

ORGDB_PACKAGE="$1"
LIB_DIR="$2"

mkdir -p "$LIB_DIR"

echo "Installing $ORGDB_PACKAGE into $LIB_DIR (via $CONTAINER)..."
singularity exec --bind /beevol/home "$CONTAINER" Rscript -e "
  BiocManager::install('${ORGDB_PACKAGE}', lib = '${LIB_DIR}', update = FALSE, ask = FALSE, dependencies = TRUE)
"

echo "Verifying it loads standalone from $LIB_DIR..."
singularity exec --bind /beevol/home "$CONTAINER" Rscript -e "
  .libPaths(c('${LIB_DIR}', .libPaths()))
  library(${ORGDB_PACKAGE})
  cat('OK --', '${ORGDB_PACKAGE}', 'version', as.character(packageVersion('${ORGDB_PACKAGE}')), 'loaded from', find.package('${ORGDB_PACKAGE}'), '\n')
"

echo ""
echo "Done. To use this in a RiboSeq_TE run, set in its RiboSeq_TE_samples.yaml:"
echo "  GO_ORGANISM_DB: \"${ORGDB_PACKAGE}\""
echo "  GO_ANNOTATION_LIB_DIR: \"${LIB_DIR}\""
