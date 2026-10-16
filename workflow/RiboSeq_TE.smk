# ===== Snake file for Ribo-seq / RNA-seq translational-efficiency (TE) join ===
# Phase 2 of the Ribo-seq pipeline work. Unlike every other top-level smk file
# in this repo, this one does no alignment/quantification of its own -- it
# only reads two ALREADY-COMPLETED projects' featureCounts outputs
# (RiboSeq.smk/RiboSeq_UMI.smk's and RNAseq.smk's), so there's no
# load_index_config/process_samples/FASTQS boilerplate here.

shell.executable("/bin/bash")
shell.prefix("[ -f ~/.bash_profile ] && source ~/.bash_profile; set -o nounset -o pipefail -o errexit -x; ")

import os

# Reused for memory_estimator only -- no index/sample-processing helpers
# apply here.
include: workflow.source_path("scripts/funs.py")

SINGULARITY = config.get("CONTAINER", "")
USE_SINGULARITY = bool(SINGULARITY and SINGULARITY.strip())
if USE_SINGULARITY:
    singularity:
        config["CONTAINER"]

PROJ        = config.get("PROJ")
RIBO_PROJ   = config.get("RIBO_PROJ")
RIBO_INDEX  = config.get("RIBO_INDEX")
RNA_PROJ    = config.get("RNA_PROJ")
RNA_INDEX   = config.get("RNA_INDEX")
CONTRASTS   = config.get("CONTRASTS")
CMD_PARAMS  = config.get("CMD_PARAMS", {})
USER        = config.get("USER")
RUN_GO      = config.get("run_go_enrichment", False)

os.makedirs(PROJ + "/report", exist_ok=True)

CONTRAST_NAMES = [c["name"] for c in CONTRASTS]
CONTRAST_BY_NAME = {c["name"]: c for c in CONTRASTS}

# CONTRASTS' group1/group2 are nicknames shared between the two completed
# runs -- but each run's OWN fastq-prefix naming is independent (e.g.
# "ribo_untr_a" vs "mrna_untr_a" for the same nickname), since
# featureCounts files are named after the fastq prefix, not the nickname.
# Translate through each project's own samples.yaml before constructing
# featureCounts paths.
RIBO_NICKNAME_MAP = load_nickname_map(config.get("RIBO_SAMPLES_YAML"))
RNA_NICKNAME_MAP  = load_nickname_map(config.get("RNA_SAMPLES_YAML"))

wildcard_constraints:
    contrast = r"[a-zA-Z0-9_\-]+"

rule all:
    input:
        expand(PROJ + "/report/{contrast}_TE_results.tsv", contrast=CONTRAST_NAMES),
        expand(PROJ + "/report/{contrast}_MA_plot.png", contrast=CONTRAST_NAMES),
        expand(PROJ + "/report/{contrast}_volcano_plot.png", contrast=CONTRAST_NAMES),
        expand(PROJ + "/report/{contrast}_4quadrant_plot.png", contrast=CONTRAST_NAMES),
        [] if not RUN_GO else expand(PROJ + "/report/{contrast}_GO_enrichment.tsv", contrast=CONTRAST_NAMES)


rule te_analysis:
    input:
        # Optionally uses rules/03e_riboseq_psite_cds_counts.snake's
        # P-site-restricted, CDS-only counts instead of plain featureCounts
        # (same 9-column schema, drop-in) -- off by default since it depends
        # on RIBO_PROJ having been run with enable_psite_cds_counts: true.
        ribo_fc = lambda wildcards: expand(
            RIBO_PROJ + "/counts/{sample}_aligned_{index}_" +
                ("psite_cds_featureCounts.tsv" if CMD_PARAMS.get("use_psite_counts") else "featureCounts.tsv"),
            sample = [RIBO_NICKNAME_MAP[n] for n in
                      CONTRAST_BY_NAME[wildcards.contrast]["group1"] + CONTRAST_BY_NAME[wildcards.contrast]["group2"]],
            index = RIBO_INDEX
        ),
        rna_fc = lambda wildcards: expand(
            RNA_PROJ + "/counts/{sample}_aligned_{index}_featureCounts.tsv",
            sample = [RNA_NICKNAME_MAP[n] for n in
                      CONTRAST_BY_NAME[wildcards.contrast]["group1"] + CONTRAST_BY_NAME[wildcards.contrast]["group2"]],
            index = RNA_INDEX
        )
    output:
        results = PROJ + "/report/{contrast}_TE_results.tsv",
        ma       = PROJ + "/report/{contrast}_MA_plot.png",
        volcano  = PROJ + "/report/{contrast}_volcano_plot.png",
        quad     = PROJ + "/report/{contrast}_4quadrant_plot.png"
    params:
        job_name = "{contrast}_te_analysis",
        script   = "workflow/scripts/te_analysis.R",
        group1   = lambda wildcards: " ".join(CONTRAST_BY_NAME[wildcards.contrast]["group1"]),
        group2   = lambda wildcards: " ".join(CONTRAST_BY_NAME[wildcards.contrast]["group2"]),
        filter_biotype = CMD_PARAMS.get("filter_biotype", ""),
        outdir   = PROJ + "/report"
    resources:
        # min bumped 2->16: OOM-killed at ~1.3GB against the 2GB floor --
        # xtail's internal permutation-based significance testing scales
        # with gene/permutation count, not featureCounts input file size
        # (a few MB here), so the file-size estimate was never going to be
        # representative. Same pattern as sortmerna/riboWaltz elsewhere in
        # this pipeline.
        memory = lambda wildcards, input: memory_estimator(input, 1, 16, 32)
    threads:
        4
    shell:
        """
        mkdir -p {params.outdir}
        Rscript {params.script} \
            "{input.ribo_fc}" \
            "{input.rna_fc}" \
            "{params.group1}" \
            "{params.group2}" \
            "{wildcards.contrast}" \
            "{params.filter_biotype}" \
            "{threads}" \
            "{params.outdir}"
        """


rule te_go_enrichment:
    input:
        PROJ + "/report/{contrast}_TE_results.tsv"
    output:
        PROJ + "/report/{contrast}_GO_enrichment.tsv"
    params:
        job_name = "{contrast}_go_enrichment",
        script   = "workflow/scripts/te_go_enrichment.R",
        orgdb    = config.get("GO_ORGANISM_DB", "org.Hs.eg.db")
    resources:
        memory = lambda wildcards, input: memory_estimator(input, 1, 2, 8)
    threads:
        1
    shell:
        """
        Rscript {params.script} "{input}" "{params.orgdb}" "{output}"
        """
