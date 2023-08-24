#! /usr/bin/env bash

FUNS_DIR=$(dirname $(realpath ${BASH_SOURCE[0]}))


# Count reads and normalize
intersect_reads() {
    local reads=$1
    local region=$2
    local threads=$3
    local out=$4
    local out_N=$5
    local flags=$6

    local tmp_dir=$(dirname $out)
    local tmp_1=$(mktemp -p $tmp_dir tmp_bed.XXXXX)
    local tmp_2=$(mktemp -p $tmp_dir tmp_bed.XXXXX)

    zcat $reads \
        > $tmp_1

    cat $tmp_1 \
        | bedtools intersect -sorted -a $region -b - $flags \
        | sort -S1G --parallel=$threads -k1,1 -k2,2n \
        > $tmp_2

    $FUNS_DIR/norm_bed $tmp_1 $tmp_2 7 -len \
        | sort -S1G --parallel=$threads -k1,1 -k2,2n \
        | pigz -p $threads \
        > $out_N

    cat $tmp_2 \
        | pigz -p $threads \
        > $out

    rm $tmp_1 $tmp_2
}
