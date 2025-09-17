#!/usr/bin/env python3
# Modified from Anshul Kundaje's lab
# piped script to take multimappers and randomly select up to -k alignments
# requires a qname sorted file!!
import sys
import random
import argparse

def parse_args():
    '''
    Gives options
    '''
    parser = argparse.ArgumentParser(
        description='Randomly selects up to k alignments per read and discards the rest')
    parser.add_argument('-k', help='Maximum number of alignments to keep per read', required=True)
    parser.add_argument('--paired-end', dest='paired_ended',
                        action='store_true', help='Data is paired-end')
    parser.add_argument('--seed', type=int, help='Random seed for reproducibility')
    args = parser.parse_args()
    alignment_cutoff = int(args.k)
    paired_ended = args.paired_ended
    seed = args.seed
    return alignment_cutoff, paired_ended, seed

def process_read_group(reads, alignment_cutoff):
    '''
    Randomly select up to alignment_cutoff reads from the group
    '''
    if len(reads) <= alignment_cutoff:
        # Keep all reads if we have fewer than or equal to the cutoff
        return reads
    else:
        # Randomly sample alignment_cutoff reads
        return random.sample(reads, alignment_cutoff)

if __name__ == "__main__":
    '''
    Runs the random selection of multimapped reads
    '''
    alignment_cutoff, paired_ended, seed = parse_args()
    
    # Set random seed for reproducibility if provided
    if seed is not None:
        random.seed(seed)
    
    if paired_ended:
        alignment_cutoff = int(alignment_cutoff) * 2
    
    # Store each line in sam file as a list of reads
    current_reads = []
    current_qname = ''
    
    for line in sys.stdin:
        read_elems = line.strip().split('\t')
        
        # Pass through header lines
        if read_elems[0].startswith('@'):
            sys.stdout.write(line)
            continue
        
        # Keep taking lines that have the same qname
        if read_elems[0] == current_qname:
            # Add line to current reads
            current_reads.append(line)
        else:
            # Process the previous read group if it exists
            if len(current_reads) > 0:
                # Randomly select up to alignment_cutoff reads
                selected_reads = process_read_group(current_reads, alignment_cutoff)
                
                # Output the selected reads
                for read in selected_reads:
                    sys.stdout.write(read)
            
            # Start new read group
            current_reads = [line]
            current_qname = read_elems[0]
    
    # Process the final read group
    if len(current_reads) > 0:
        selected_reads = process_read_group(current_reads, alignment_cutoff)
        for read in selected_reads:
            sys.stdout.write(read)
