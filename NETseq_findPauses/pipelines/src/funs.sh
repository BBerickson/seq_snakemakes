#! /usr/bin/env bash

FUNS_DIR=$(dirname "$(realpath "${BASH_SOURCE[0]}")")


# Count reads and normalize
intersect_reads() {
    local reads="$1"
    local region="$2"
    local threads="$3"
    local out="$4"
    local out_N="$5"
    local flags=$6

    local tmp_dir=$(dirname "$out")
    local tmp_1=$(mktemp -p $tmp_dir tmp_bed.XXXXX)
    local tmp_2=$(mktemp -p $tmp_dir tmp_bed.XXXXX)

    zcat "$reads" \
        > $tmp_1

    cat $tmp_1 \
        | bedtools intersect -sorted -a "$region" -b - $flags \
        | sort -S1G --parallel="$threads" -k1,1 -k2,2n \
        > $tmp_2

    "$FUNS_DIR/norm_bed" $tmp_1 $tmp_2 7 -len \
        | sort -S1G --parallel="$threads" -k1,1 -k2,2n \
        | pigz -p "$threads" \
        > "$out_N"

    cat $tmp_2 \
        | pigz -p "$threads" \
        > "$out"

    rm $tmp_1 $tmp_2
}


# Create bigwigs
create_bigwig() {
    local reads="$1"
    local strand="$2"
    local mtplyr="$3"
    local bg="$4"
    local bg_N="$5"
    local bw="$6"
    local chroms="$7"
    local threads="$8"

    tmp_1=$(mktemp tmp.XXXXX)
    tmp_2=$(mktemp tmp.XXXXX)
    tmp_3=$(mktemp tmp.XXXXX)

    zcat -f "$reads" \
        > $tmp_1

    cat $tmp_1 \
        | sort -S1G --parallel="$threads" -k1,1 -k2,2n \
        | bedtools genomecov -bg $strand -i - -g "$chroms" \
        > $tmp_2

    $FUNS_DIR/norm_bed $tmp_1 $tmp_2 4 \
        | sort -S1G --parallel="$threads" -k1,1 -k2,2n \
        | awk -v mtplyr="$mtplyr" -v OFS="\t" '{$4 = $4 * mtplyr; print}' \
        > $tmp_3
    
    bedGraphToBigWig \
        $tmp_3 \
        "$chroms" \
        "$bw"

    cat $tmp_2 \
        | pigz -p "$threads" \
        > "$bg"

    cat $tmp_3 \
        | pigz -p "$threads" \
        > "$bg_N"

    rm $tmp_1 $tmp_2 $tmp_3
}


# Write URLs
create_urls() {
    local paths=$1
    local ssh="$2"
    local ssh_key_dir="$3"
    local url="$4"
    local url_file="$5"
    local grp_re="${6:-}"

    local header="URL\tSAMPLE\tSTRAND"
    local cmd='echo -e $bw_url\\t$sam\\t$strand'

    # When running in docker container $HOME becomes workflow working directory
    # need to link .ssh folder to ensure ssh credentials are available within
    # container
    if [ ! -d "$HOME/.ssh" ] && [ ! -L "$HOME/.ssh" ]
    then
        ln -s "$ssh_key_dir" "$HOME"
    fi

    if [ ! -z "$grp_re" ]
    then
        local header="${header/STRAND/GROUP\\\tSTRAND}"
        local cmd="${cmd/strand/grp\\\\t\$strand}"
    fi

    if [ ! -s "$url_file" ]
    then
        echo -e "$header" \
            > "$url_file"
    fi

    # Create remote directory if it does not exist
    # add trailing slash to ssh path to ensure symlink is not overwritten
    # by rsync
    local ssh="$ssh/"
    local hostname=$(echo "$ssh" | cut -d ':' -f 1)
    local hostdir=$(echo "$ssh" | cut -d ':' -f 2-)

    ssh -o StrictHostKeyChecking=no "$hostname" "
    if [ ! -d \"\$(readlink -f $hostdir)\" ];
    then
        mkdir -p -m 755 $hostdir;
    fi
    "

    # Transfer bigwigs to amc-sandbox
    for path in $paths
    do
        local bw=$(basename "$path")
        local bw_url="$url/$bw"
        local strand=$(echo "$bw" | grep -oP "(?<=_)[a-zA-Z]+(?=.[a-z]+$)")
        local sffx="$strand"
        
        if [ ! -z "$grp_re" ]
        then
            local grp=$(echo "$bw" | grep -oP "(?<=-)""$grp_re""(?=_"$strand".[a-z]+$)")
            local sffx="$grp"
        fi
        
        local sam=$(echo "$bw" | grep -oP ".+(?<=[_-])""$sffx")

        rsync -e 'ssh -o StrictHostKeyChecking=no' --perms --chmod=ugo+r "$path" "$ssh"

        eval "$cmd" \
            >> "$url_file"
    done
}


