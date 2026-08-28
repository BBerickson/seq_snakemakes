#!/usr/bin/env bash

#SBATCH --job-name=ChIPseq
#SBATCH --output=logs/ChIPseq_%j.out
#SBATCH --error=logs/ChIPseq_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition=normal 
#SBATCH --qos=normal

set -o nounset -o pipefail -o errexit -x

mkdir -p logs

# Configuration
PROFILE="workflow/profiles/Bodhi"
LSF_CONFIG="workflow/profiles/Bodhi/Bodhi_config.yaml"
SSH_KEY_DIR="${HOME}/.ssh"
SAMPLES_FILE="ChIPseq_dual_IP_samples.yaml"
MATRIX_FILE="UnStranded_matrix.yaml"
SINGULARITY_PREFIX="/beevol/home/${USER}/.singularity_cache"

snakemake \
    --profile ${PROFILE} \
    --snakefile workflow/ChIPseq_dual_IP.smk \
    --configfile ${SAMPLES_FILE} ${LSF_CONFIG} \
    --singularity-prefix "${SINGULARITY_PREFIX}" \
    --config SSH_KEY_DIR="${SSH_KEY_DIR}"

snakemake \
    --profile ${PROFILE} \
    --snakefile workflow/UnStranded_matrix.smk \
    --configfile ${SAMPLES_FILE} ${LSF_CONFIG} ${MATRIX_FILE} \
    --singularity-prefix "${SINGULARITY_PREFIX}" \
    --config SSH_KEY_DIR="${SSH_KEY_DIR}" \
    --rerun-triggers mtime input code

SAMPLES_FILE2="ChIPseq_dual_IP_tRNA_samples.yaml"
MATRIX_FILE="workflow/small_matrix.yaml"

# SAMPLES_FILE2 swaps primary/spike roles to build matrix files "of the spike
# alignment" using its own scalefactor (see workflow/ref/readme.txt). 
# make a dummy results.tsv to let run trigger
SWAPPED_PRIMARY=$(python3 -c "
import yaml
cfg = yaml.safe_load(open('${SAMPLES_FILE2}'))
print(next(e['name'] for e in cfg['INDEXES'] if e.get('role') == 'primary'))
")
PROJ_NAME=$(python3 -c "
import yaml
cfg = {}
for f in ['${SAMPLES_FILE}', '${SAMPLES_FILE2}']:
    cfg.update(yaml.safe_load(open(f)))
print(cfg['PROJ'])
")
MARKER="${PROJ_NAME}/report/${PROJ_NAME}_${SWAPPED_PRIMARY}_results.tsv"
if [ ! -f "${MARKER}" ]; then
    mkdir -p "$(dirname "${MARKER}")"
    touch "${MARKER}"
fi

snakemake \
    --profile ${PROFILE} \
    --snakefile workflow/UnStranded_matrix.smk \
    --configfile ${SAMPLES_FILE} ${SAMPLES_FILE2} ${LSF_CONFIG} ${MATRIX_FILE} \
    --singularity-prefix "${SINGULARITY_PREFIX}" \
    --config SSH_KEY_DIR="${SSH_KEY_DIR}" \
    --rerun-triggers mtime input code
    