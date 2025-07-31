#! /usr/bin/env bash

#BSUB -J ggsashimi[1-2]
#BSUB -e logs/ggsashimi.%J.%I.err
#BSUB -o logs/ggsashimi.%J.%I.out
#BSUB -R "rusage[mem=4] span[hosts=1]"

mkdir -p logs

UPDWOWN=(more_inclusion
more_skipping)

updown=${UPDWOWN[$(($LSB_JOBINDEX - 1))]}
mkdir -p $updown

gtf=/beevol/home/erickson/ref/hg38_p14/base_ref/hg38_chrfix_GRCh38.113.gtf

# need to run makeSashimi() from rMATS_plots.R first
tc="treatment:control" # same as GROUPS_COMP: in samples_rAMTS.yaml
PROJ="180629_SY351" # same as PROJ: in samples.yaml

while IFS=$'\t' read -r col1 col2; do
    python /beevol/home/erickson/ggsashimi/ggsashimi_BE.py \
        -b $PROJ/rmats/$tc\_bams.tsv  \
        -c $col1 \
        -g $gtf \
        -P $PROJ/rmats/$tc\_palette.txt \
        -o $updown/$col2 \
        -M 10 -C 3 -O 3 \
        --fix-y-scale \
        --shrink --alpha 0.25 \
        --base-size=20 --ann-height=4 \
        --height=3 --width=18 \
        --strand MATE1_SENSE \
        --out-strand plus \
        --out-format png
done < $PROJ/rmats/$tc\_$updown\_pos.txt

while IFS=$'\t' read -r col1 col2; do
    python /beevol/home/erickson/ggsashimi/ggsashimi_BE.py \
        -b $PROJ/rmats/$tc\_bams.tsv  \
        -c $col1 \
        -g $gtf \
        -P $PROJ/rmats/$tc\_palette.txt \
        -o $updown/$col2 \
        -M 10 -C 3 -O 3 \
        --fix-y-scale \
        --shrink --alpha 0.25 \
        --base-size=20 --ann-height=4 \
        --height=3 --width=18 \
        --strand MATE1_SENSE \
        --out-strand minus \
        --out-format png
done < $PROJ/rmats/$tc\_$updown\_neg.txt


Rscript /beevol/home/erickson/Ben_pipelines/R_scripts/pngToPDF.R $updown $tc\_$updown.pdf
rm -r $updown

#### tab in bams.tsv can be a bit sensitive ###

# usage: ggsashimi.py [-h] -b BAM -c COORDINATES [-o OUT_PREFIX] [-S OUT_STRAND]
#                     [-M MIN_COVERAGE] [-j JUNCTIONS_BED] [-g GTF] [-s STRAND]
#                     [--shrink] [-O OVERLAY] [-A AGGR] [-C COLOR_FACTOR]
#                     [--alpha ALPHA] [-P PALETTE] [-L LABELS] [--fix-y-scale]
#                     [--height HEIGHT] [--ann-height ANN_HEIGHT]
#                     [--width WIDTH] [--base-size BASE_SIZE] [-F OUT_FORMAT]
#                     [-R OUT_RESOLUTION] [--debug-info] [--version]
# 
# Create sashimi plot for a given genomic region
# 
# optional arguments:
#   -h, --help            show this help message and exit
#   -b BAM, --bam BAM     Individual bam file or file with a list of bam files.
#                         In the case of a list of files the format is tsv:
#                         1col: id for bam file, 2col: path of bam file, 3+col:
#                         additional columns
#   -c COORDINATES, --coordinates COORDINATES
#                         Genomic region. Format: chr:start-end. Remember that
#                         bam coordinates are 0-based
#   -o OUT_PREFIX, --out-prefix OUT_PREFIX
#                         Prefix for plot file name [default=sashimi]
#   -S OUT_STRAND, --out-strand OUT_STRAND
#                         Only for --strand other than 'NONE'. Choose which
#                         signal strand to plot: <both> <plus> <minus>
#                         [default=both]
#   -M MIN_COVERAGE, --min-coverage MIN_COVERAGE
#                         Minimum number of reads supporting a junction to be
#                         drawn [default=1]
#   -j JUNCTIONS_BED, --junctions-bed JUNCTIONS_BED
#                         Junction BED file name [default=no junction file]
#   -g GTF, --gtf GTF     Gtf file with annotation (only exons is enough)
#   -s STRAND, --strand STRAND
#                         Strand specificity: <NONE> <SENSE> <ANTISENSE>
#                         <MATE1_SENSE> <MATE2_SENSE> [default=NONE]
#   --shrink              Shrink the junctions by a factor for nicer display
#                         [default=False]
#   -O OVERLAY, --overlay OVERLAY
#                         Index of column with overlay levels (1-based)
#   -A AGGR, --aggr AGGR  Aggregate function for overlay: <mean> <median>
#                         <mean_j> <median_j>. Use mean_j | median_j to keep
#                         density overlay but aggregate junction counts
#                         [default=no aggregation]
#   -C COLOR_FACTOR, --color-factor COLOR_FACTOR
#                         Index of column with color levels (1-based)
#   --alpha ALPHA         Transparency level for density histogram [default=0.5]
#   -P PALETTE, --palette PALETTE
#                         Color palette file. tsv file with >=1 columns, where
#                         the color is the first column. Both R color names and
#                         hexadecimal values are valid
#   -L LABELS, --labels LABELS
#                         Index of column with labels (1-based) [default=1]
#   --fix-y-scale         Fix y-scale across individual signal plots
#                         [default=False]
#   --height HEIGHT       Height of the individual signal plot in inches
#                         [default=2]
#   --ann-height ANN_HEIGHT
#                         Height of annotation plot in inches [default=1.5]
#   --width WIDTH         Width of the plot in inches [default=10]
#   --base-size BASE_SIZE
#                         Base font size of the plot in pch [default=14]
#   -F OUT_FORMAT, --out-format OUT_FORMAT
#                         Output file format: <pdf> <svg> <png> <jpeg> <tiff>
#                         [default=pdf]
#   -R OUT_RESOLUTION, --out-resolution OUT_RESOLUTION
#                         Output file resolution in PPI (pixels per inch).
#                         Applies only to raster output formats [default=300]
#   --debug-info          Show several system information useful for debugging
#                         purposes [default=None]
#   --version             show program's version number and exit
