#!/usr/bin/env python3

import re
import sys
import argparse
from pathlib import Path

def reverse_complement(sequence):
    """Generate reverse complement of DNA sequence"""
    complement = str.maketrans('ATCG', 'TAGC')
    return sequence.translate(complement)[::-1]

def build_satellite_regex(monomer, min_copies):
    """Build regex patterns for all rotations of satellite monomer"""
    patterns = []
    # Generate all rotations of the monomer
    for i in range(len(monomer)):
        rotated = monomer[i:] + monomer[:i]
        patterns.append(f"({rotated}){{{min_copies},}}")
    return "|".join(patterns)

def analyze_sequence(seq_id, sequence, monomer, min_copies=1):
    """Analyze a single sequence for satellite content"""
    seq_len = len(sequence)
    
    # Build regex patterns for forward and reverse complement
    forward_pattern = re.compile(build_satellite_regex(monomer, min_copies))
    rc_monomer = reverse_complement(monomer)
    reverse_pattern = re.compile(build_satellite_regex(rc_monomer, min_copies))
    
    # Find matches
    forward_bp = 0
    reverse_bp = 0
    max_copies_forward = 0
    max_copies_reverse = 0
    
    # Analyze forward strand
    for match in forward_pattern.finditer(sequence):
        match_len = match.end() - match.start()
        forward_bp += match_len
        copies = match_len // len(monomer)
        max_copies_forward = max(max_copies_forward, copies)
    
    # Analyze reverse strand
    for match in reverse_pattern.finditer(sequence):
        match_len = match.end() - match.start()
        reverse_bp += match_len
        copies = match_len // len(monomer)
        max_copies_reverse = max(max_copies_reverse, copies)
    
    # Take the best result (forward or reverse)
    sat_block_bp = max(forward_bp, reverse_bp)
    max_copies = max(max_copies_forward, max_copies_reverse)
    p_sat_block = (sat_block_bp / seq_len) * 100 if seq_len > 0 else 0
    
    return {
        'name': seq_id,
        'seq_len': seq_len,
        'sat_block_bp': sat_block_bp,
        'p_sat_block': p_sat_block,
        'max_copies': max_copies
    }

def extract_read_name(header):
    """Extract just the read name from complex FASTA headers"""
    # Remove the '>' character first
    header = header[1:] if header.startswith('>') else header
    
    # Split by space and take the first part (the actual read ID)
    read_name = header.split()[0]
    
    # For Nanopore reads, sometimes there are additional suffixes like /1
    # Remove common suffixes
    if read_name.endswith('/1') or read_name.endswith('/2'):
        read_name = read_name[:-2]
    
    return read_name

def read_fasta(filename):
    """Simple FASTA reader"""
    sequences = {}
    current_id = None
    current_seq = []
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    sequences[current_id] = ''.join(current_seq)
                # Extract clean read name from header
                current_id = extract_read_name(line)
                current_seq = []
            else:
                current_seq.append(line)
        
        # Don't forget the last sequence
        if current_id:
            sequences[current_id] = ''.join(current_seq)
    
    return sequences

def main():
    parser = argparse.ArgumentParser(description='Simple satellite DNA detection')
    parser.add_argument('--input_file', required=True, help='Input FASTA file')
    parser.add_argument('--satellite', required=True, help='Satellite monomer sequence')
    parser.add_argument('--min_copies', type=int, default=1, help='Minimum number of tandem copies')
    parser.add_argument('--output', help='Output file (default: stdout)')
    
    args = parser.parse_args()
    
    # Read sequences
    sequences = read_fasta(args.input_file)
    
    # Prepare output
    output_file = open(args.output, 'w') if args.output else sys.stdout
    
    # Write header
    header = "name$\tseq_len\tsat_block_bp\tP_sat_block\tmax_copy"
    print(header, file=output_file)
    
    # Analyze each sequence
    for seq_id, sequence in sequences.items():
        result = analyze_sequence(seq_id, sequence, args.satellite, args.min_copies)
        
        # Format output to match original
        output_line = f"{result['name']}\t{result['seq_len']}\t{result['sat_block_bp']}\t{result['p_sat_block']:.1f}\t{result['max_copies']}"
        print(output_line, file=output_file)
    
    if args.output:
        output_file.close()

if __name__ == '__main__':
    main()
