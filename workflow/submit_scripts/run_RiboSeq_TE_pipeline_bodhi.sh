#!/usr/bin/env bash
# Lightweight orchestrator for the full Ribo-seq -> RNA-seq -> TE pipeline,
# chaining the three ALREADY-WORKING, independently-maintained driver
# scripts (run_RiboSeq[_UMI]_bodhi.sh / run_RNAseq_*_bodhi.sh /
# run_RiboSeq_TE_bodhi.sh) via SLURM job dependencies, rather than merging
# them into one Snakemake DAG. A true merge isn't worth it here: every
# pipeline in this repo (BRUseq, RNAseq, RiboSeq, ...) is a self-contained
# .smk file with its own module-level config variables (PROJ, SAMPLES,
# CMD_PARAMS, ...) -- RiboSeq.smk and RNAseq.smk both define the same
# names, so combining them into one DAG would need Snakemake's module/
# namespacing system just for this one pipeline, making it the one outlier
# in the whole repo. This script only touches orchestration, not rule logic.
#
# The Ribo-seq and RNA-seq arms are fully independent (different samples,
# different outputs, no shared files) and are submitted IN PARALLEL. The TE
# arm is submitted with --dependency=afterok on both, so it waits until
# both finish successfully before SLURM starts it -- if either arm fails,
# the TE job is never started (SLURM cancels it, per afterok semantics).
# TE is skipped entirely if no RNA-seq arm is given (Ribo-seq-only use case
# needs nothing beyond the existing run_RiboSeq[_UMI]_bodhi.sh script).
#
# Usage:
#   ./run_RiboSeq_TE_pipeline_bodhi.sh --ribo <ribo_project_dir> \
#       [--rna <rna_project_dir> --te <te_project_dir>] \
#       [--ribo-script run_RiboSeq_UMI_bodhi.sh] \
#       [--rna-script run_RNAseq_Hisat_bodhi.sh]
#
# Each <..._dir> is a project directory already set up exactly as if you
# were about to run its own driver script by hand (samples.yaml, workflow/,
# logs/ present) -- this script does not create or validate that setup, it
# only submits the sbatch jobs and wires the dependency between them.

set -o nounset -o pipefail -o errexit

usage() {
    echo "Usage: $0 --ribo <ribo_project_dir> [--rna <rna_project_dir> --te <te_project_dir>] [--ribo-script <name>] [--rna-script <name>]" >&2
    exit 1
}

RIBO_DIR=""
RNA_DIR=""
TE_DIR=""
RIBO_SCRIPT="run_RiboSeq_bodhi.sh"
RNA_SCRIPT="run_RNAseq_STAR_bodhi.sh"
TE_SCRIPT="run_RiboSeq_TE_bodhi.sh"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --ribo)        RIBO_DIR="$2"; shift 2 ;;
        --rna)         RNA_DIR="$2"; shift 2 ;;
        --te)          TE_DIR="$2"; shift 2 ;;
        --ribo-script) RIBO_SCRIPT="$2"; shift 2 ;;
        --rna-script)  RNA_SCRIPT="$2"; shift 2 ;;
        --te-script)   TE_SCRIPT="$2"; shift 2 ;;
        *) usage ;;
    esac
done

[[ -z "$RIBO_DIR" ]] && usage
if [[ -n "$TE_DIR" && -z "$RNA_DIR" ]]; then
    echo "Error: --te requires --rna (TE analysis needs a completed RNA-seq arm)." >&2
    exit 1
fi

echo "Submitting Ribo-seq arm ($RIBO_DIR/$RIBO_SCRIPT)..."
RIBO_JOBID=$(cd "$RIBO_DIR" && sbatch --parsable "$RIBO_SCRIPT")
echo "  Ribo-seq driver job: $RIBO_JOBID"

if [[ -n "$RNA_DIR" ]]; then
    echo "Submitting RNA-seq arm ($RNA_DIR/$RNA_SCRIPT), in parallel with Ribo-seq..."
    RNA_JOBID=$(cd "$RNA_DIR" && sbatch --parsable "$RNA_SCRIPT")
    echo "  RNA-seq driver job: $RNA_JOBID"

    if [[ -n "$TE_DIR" ]]; then
        echo "Submitting TE arm ($TE_DIR/$TE_SCRIPT), deferred until both arms complete (--dependency=afterok:$RIBO_JOBID:$RNA_JOBID)..."
        TE_JOBID=$(cd "$TE_DIR" && sbatch --parsable --dependency="afterok:${RIBO_JOBID}:${RNA_JOBID}" "$TE_SCRIPT")
        echo "  TE driver job: $TE_JOBID (waits on $RIBO_JOBID and $RNA_JOBID)"
    fi
fi
