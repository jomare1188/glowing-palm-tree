#!/usr/bin/env python

import os
import sys
from Bio import SeqIO

def split_fasta(input_file, num_parts, output_dir):
    # Open input file and parse sequences
    sequences = []
    with open(input_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            sequences.append(record)
    
    # Determine the number of sequences per file
    num_seqs_per_file = len(sequences) // num_parts
    remainder = len(sequences) % num_parts
    
    # Create output files
    for i in range(1, num_parts + 1):
        output_file = os.path.join(output_dir, f"part_{i}_{os.path.splitext(os.path.basename(input_file))[0]}.fasta")
        out_handle = open(output_file, "w")
        
        # Write sequences to output file
        seq_count = 0
        for record in sequences[(i-1)*num_seqs_per_file:i*num_seqs_per_file]:
            SeqIO.write(record, out_handle, "fasta")
            seq_count += 1
        if remainder > 0 and i == num_parts:
            for record in sequences[-remainder:]:
                SeqIO.write(record, out_handle, "fasta")
                seq_count += 1
            remainder = 0
            
        out_handle.close()

if __name__ == '__main__':
    input_file = sys.argv[1]
    num_parts = int(sys.argv[2])
    output_dir = sys.argv[3]
    split_fasta(input_file, num_parts, output_dir)

