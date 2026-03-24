# How to Make INDEX Ref Files for Snakemake Pipeline

## Standard Alignment with Ref Files
1. See hg38.yaml for an example.
    a. INDEX_MAP must match the prefix of the index files
    b. INDEXES must match the yaml file name (e.g., "hg38", "hg38.yaml")

## Standard Alignment with Subset Ref Files
1. See hg38_tRNA.yaml for an example.
    a. GENELIST value will be appended to output file names (e.g., "tRNA")

## Standard Alignment (no ref files)
1. See ecoli.yaml for an example.
    a. Ref bed files can be left empty

## Mixed, Sequential, and Dual Alignment with Ref Files
1. See hg38_mm39.yaml for an example.
    a. INDEX_MAP must match the prefix of the indexes
    b. INDEXES must match the names of other existing ref yaml files