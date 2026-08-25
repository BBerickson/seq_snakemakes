# How to Make INDEX Ref Files for Snakemake Pipeline

Each ref yaml here describes ONE genome index's biological metadata (GTF,
MASK, FA_SAMPLE, MY_REF, FW_REF/REV_REF, PI_REF, etc.). A project's
samples.yaml selects which ref yaml(s) to use and how they relate to each
other via ALIGNMENT_MODE and an INDEXES list of {name, role} pairs -- it no
longer points at a single combo "GENOME" yaml.

    ALIGNMENT_MODE: single        # single | split | sequential
    INDEXES:
      - name: hg38                # must match a ref yaml file name (hg38.yaml)
        role: primary              # optional, defaults to "primary" for 1 entry

- ALIGNMENT_MODE: single   -- one real genome, one index, one alignment pass.
- ALIGNMENT_MODE: split    -- one pre-built combined index (e.g. hg38+mm39),
  aligned once, then split post-hoc into per-genome bams by contig name.
  Requires exactly 2 INDEXES entries (role: primary + role: spike) plus the
  combined index's own INDEX_PATH/INDEX_MAP set directly in samples.yaml
  (next to INDEXES), since that combined index is specific to the project,
  not reusable per-genome metadata:

    ALIGNMENT_MODE: split
    INDEXES:
      - name: hg38
        role: primary
      - name: mm39
        role: spike
    INDEX_PATH: "/path/to/combined/index/"
    INDEX_MAP:  "hg38_mm39"

- ALIGNMENT_MODE: sequential -- two real, separate indexes, aligned in two
  passes (e.g. ChIPseq_dual_IP.smk). Also requires exactly 2 INDEXES entries
  (role: primary + role: spike); the spike-role index is aligned FIRST, the
  primary-role index SECOND and is the reported/results index. Each index's
  own INDEX_PATH/INDEX_MAP come from its own ref yaml (no extra top-level
  keys needed, unlike split mode).

## Standard alignment
See hg38.yaml for an example.
    a. INDEX_MAP must match the prefix of the index files
    b. Give the ref yaml the same name you'll reference from samples.yaml's
       INDEXES list (e.g. hg38.yaml -> name: hg38)

## Standard alignment with subset ref/region files (gene-list variants)
See hg38_tRNA.yaml vs. hg38.yaml for an example: same INDEX_PATH/INDEX_MAP
(same real alignment index) but different MASK/MY_REF/FW_REF/REV_REF and a
GENELIST tag, which gets appended to output file names (e.g. "tRNA"). Select
a variant by using its yaml's name (e.g. name: hg38_tRNA) instead of the
base genome's name for whichever INDEXES entry it applies to -- this works
identically for a single-mode primary index or a sequential-mode primary/
spike index (e.g. hg38_polIII's dual-IP with a tRNA-only primary variant).

To reuse ONE existing alignment to build matrix files for TWO (or more) gene
lists (e.g. whole-genome hg38 and tRNA-only hg38_tRNA) in a SINGLE run, add
a GENELISTS list at the top level of the matrix pipeline's samples.yaml,
alongside INDEXES -- Pattern A (alignment) files ignore this key entirely,
so it has no effect on ChIPseq.smk/RNAseq.smk/etc., only on
Stranded_matrix.smk/UnStranded_matrix.smk/bidirectional_matrix.smk:

    ALIGNMENT_MODE: single
    INDEXES:
      - name: hg38
    GENELISTS:
      - hg38
      - hg38_tRNA

Every GENELISTS entry must resolve to the SAME real alignment index as the
primary INDEXES entry (checked at startup), and every entry's GENELIST tag
must be unique (also checked at startup) since that tag is what keeps their
matrix/heatmap output filenames from colliding (e.g. "..._5_1.0k_2.0k_25bin_
norm_CPM_sense_matrix.gz" vs. "..._5_1.0k_2.0k_25bin_tRNA_norm_CPM_sense_
matrix.gz"). Bigwigs are NOT duplicated per gene list -- they're genome-wide
and shared across every GENELISTS entry, since MASK (used at the bamCoverage/
bigwig step) always comes from the primary INDEXES entry's own yaml, never
from a GENELISTS variant; only the per-region MY_REF/FW_REF/REV_REF/PI_REF
region files used at the computeMatrix step vary per gene list.

If you'd rather build gene lists one at a time instead (e.g. from different
projects, or just to keep runs smaller), omit GENELISTS and swap which yaml
the primary INDEXES entry names instead (name: hg38 vs. name: hg38_tRNA),
running the matrix pipeline once per gene list -- both approaches produce
the same output files, just via one run or several.

## Standard alignment (no ref files)
See ecoli.yaml for an example.
    a. Ref bed files can be left empty

## Scaling-factor region set (SCALING_REF)
CMD_PARAMS.scalingFactors in samples.yaml controls multiBamSummary's mode +
flags only (e.g. "BED-file --centerReads" or "bins --binSize 100
--centerReads") -- it never holds a path. In "BED-file" mode, the region
file comes from a ref yaml's SCALING_REF (see hg38_polIII.yaml), looked up
by workflow/rules/03b_results.snake per index, but preferring the spike
index's SCALING_REF over the index actually being summarized -- so a
dual-IP/split setup's primary index reuses the SAME spike-defined region
set as the spike itself, giving both alignments one shared scaling factor.
With no spike role (single mode), each index just uses its own SCALING_REF.
Leave scalingFactors empty (the default in every ref yaml here except
hg38_polIII.yaml) to skip this step entirely.

## Naming
A ref yaml's filename is what samples.yaml's INDEXES[].name looks it up by
(hg38_tRNA.yaml -> name: hg38_tRNA) -- but the ref yaml's own INDEXES: field
is what actually gets used as the index label in bam/bw/matrix filenames
and everywhere else downstream. For a plain genome yaml these are the same
value (hg38.yaml's INDEXES: "hg38"). For a gene-list variant they are NOT
the same on purpose: hg38_tRNA.yaml's own INDEXES: field is "hg38", not
"hg38_tRNA" -- that's what lets name: hg38_tRNA resolve to the SAME
alignment's bam/bw files as name: hg38, instead of looking for
nonexistent "...aligned_hg38_tRNA_..." files. When adding a new plain
genome yaml, keep its internal INDEXES: value equal to its filename; when
adding a gene-list variant, set its internal INDEXES: value to match the
base genome yaml it shares an alignment with.
