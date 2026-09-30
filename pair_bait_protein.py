#!/usr/bin/env python3
"""
make_pair_fastas.py
-------------------
Parses a multi-FASTA file, identifies a designated bait sequence (either via gene name keyword 
or a separate bait FASTA file), and generates paired/multimer FASTA files with colon separation (SEQ_A:SEQ_B) 
for AlphaFold-Multimer / ColabFold.

Usage:
    python3 make_pair_fastas.py
    python3 make_pair_fastas.py --input ./proteins/multi_fasta.fasta --bait SIX1 --output ./pair_fastas
    python3 make_pair_fastas.py --input ./proteins/multi_fasta.fasta --bait-fasta ./bait.fasta --output ./pair_fastas
    python3 make_pair_fastas.py --input ./proteins/multi_fasta.fasta --bait SIX1 --bait-copy 2 --target-copy 1
"""

import os
import re
import argparse

def parse_fasta(fasta_path):
    """Reads a FASTA file into a dict of {header: sequence}."""
    sequences = {}
    current_header = None
    current_seq = []

    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_header:
                    sequences[current_header] = ''.join(current_seq)
                current_header = line[1:]
                current_seq = []
            else:
                current_seq.append(line)
        if current_header:
            sequences[current_header] = ''.join(current_seq)

    return sequences

def extract_gene_name(header):
    """
    Extracts a clean, short gene/protein identifier from headers like:
    sp|Q9Y5J3|HEY1_HUMAN Hairy enhancer... -> HEY1_HUMAN
    """
    # Look for GN= Gene Name pattern first
    gn_match = re.search(r'GN=([A-Za-z0-9_-]+)', header)
    if gn_match:
        return gn_match.group(1)
    
    # Try parsing UniProt standard ID format (sp|ACC|NAME_SPECIES)
    parts = header.split('|')
    if len(parts) >= 3:
        return parts[2].split()[0]
    
    # Fallback to first contiguous string
    return header.split()[0].replace('/', '_').replace('\\', '_')

def main():
    parser = argparse.ArgumentParser(description="Generate paired FASTA files for ColabFold Multimer predictions.")
    parser.add_argument("-i", "--input", default="./proteins/multi_fasta.fasta", help="Path to input multi-FASTA file.")
    parser.add_argument("-o", "--output", default="./pair_fastas", help="Directory where paired FASTAs will be saved.")
    parser.add_argument("-b", "--bait", default="SIX1", help="Gene or header keyword for the bait protein (default: SIX1).")
    parser.add_argument("-bf", "--bait-fasta", default=None, help="Path to a separate single- or multi-FASTA file to use as bait.")
    parser.add_argument("--bait-copy", type=int, default=1, help="Number of bait protein copies in the multimer string (default: 1).")
    parser.add_argument("--target-copy", type=int, default=1, help="Number of target protein copies in the multimer string (default: 1).")

    args = parser.parse_args()

    if not os.path.exists(args.input):
        print(f"[ERROR] Input file not found: {args.input}")
        return

    # Automatically ensure the specified output directory exists
    os.makedirs(args.output, exist_ok=True)
    sequences = parse_fasta(args.input)

    # Locate bait protein
    bait_header = None
    bait_seq = None

    if args.bait_fasta:
        if not os.path.exists(args.bait_fasta):
            print(f"[ERROR] Bait FASTA file not found: {args.bait_fasta}")
            return
        bait_sequences = parse_fasta(args.bait_fasta)
        if not bait_sequences:
            print(f"[ERROR] No sequences found in bait FASTA file: {args.bait_fasta}")
            return
        # Use the first sequence in the separate bait FASTA file
        bait_header, bait_seq = next(iter(bait_sequences.items()))
    else:
        for header, seq in sequences.items():
            if args.bait.upper() in header.upper():
                bait_header = header
                bait_seq = seq
                break

    if not bait_header:
        print(f"[ERROR] Could not find any header matching bait keyword '{args.bait}' in {args.input}")
        print("Available headers in file:")
        for h in sequences.keys():
            print(f"  - {h}")
        return

    bait_gene = extract_gene_name(bait_header)
    print(f"==================================================")
    print(f"Bait Protein Identified : {bait_gene}")
    print(f"Bait Sequence Length     : {len(bait_seq)} aa")
    print(f"Bait Copy Count          : {args.bait_copy}")
    print(f"Target Copy Count        : {args.target_copy}")
    print(f"Output Directory        : {args.output}")
    print(f"==================================================\n")

    generated_count = 0

    for header, seq in sequences.items():
        if header == bait_header:
            continue  # Skip pairing bait with itself

        partner_gene = extract_gene_name(header)
        pair_id = f"{bait_gene}_X_{partner_gene}"
        output_filename = os.path.join(args.output, f"{pair_id}.fasta")

        # Construct colon-separated sequence based on requested copy counts
        chain_seqs = ([bait_seq] * args.bait_copy) + ([seq] * args.target_copy)
        multimer_seq_str = ":".join(chain_seqs)

        # Write paired format required by ColabFold
        with open(output_filename, 'w') as out_f:
            out_f.write(f">{pair_id}\n")
            out_f.write(f"{multimer_seq_str}\n")

        print(f"[CREATED] {output_filename}")
        print(f"          {bait_gene} ({len(bait_seq)} aa) x{args.bait_copy} : {partner_gene} ({len(seq)} aa) x{args.target_copy}")
        generated_count += 1

    print(f"\n[DONE] Successfully created {generated_count} paired FASTA files in '{args.output}'.")

if __name__ == "__main__":
    main()
