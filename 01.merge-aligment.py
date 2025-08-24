#!/usr/bin/env python3
"""
Script: 01.merge_alignment.py
Author: Liqiang Hou
Date: 2025.08

Description:
    Concatenate multiple sequence alignments (MSA) from different genes into a single supermatrix.
    This script reads multiple FASTA alignment files and combines sequences with the same ID,
    filling gaps with dashes (-) for missing sequences.

Usage:
    Method 1 - Command line arguments (recommended):
        python 01.merge_alignment.py -i <input_dir> -o <output.fasta>
    
    Method 2 - Edit paths in script:
        Edit input_directory and output_file_path in __main__ section, then run:
        python 01.merge_alignment.py

Arguments:
    -i, --input   : Directory containing FASTA alignment files to concatenate
    -o, --output  : Output file path for concatenated alignment (default: concatenated_sequences.fasta)

Input:
    - Directory containing multiple FASTA format alignment files
    - Each file should contain aligned sequences (same length within file)
    - Sequence IDs should be consistent across files

Output:
    - Single FASTA file with concatenated sequences
    - Missing sequences are filled with gaps (-)
    - Sequences are wrapped at 60 characters per line

Example:
    Input files:
        gene1.fasta:
            >sp1
            ATCG
            >sp2
            ACCG
        
        gene2.fasta:
            >sp1
            GGAA
            >sp3
            GTAT
    
    Output (concatenated_sequences.fasta):
        >sp1
        ATCGGGAA
        >sp2
        ACCG----
        >sp3
        ----GTAT

Requirements:
    - Python 3.6+
    - No external dependencies required

Notes:
    - Assumes all sequences within a single file have the same length (properly aligned)
    - IDs are sorted numerically if they are numbers, otherwise alphabetically
    - Memory efficient: processes files sequentially
"""

import os
import argparse
import sys

def read_fasta(file_path):
    """
    Read a FASTA file and return a dictionary {ID: sequence}
    
    Args:
        file_path (str): Path to FASTA file
    
    Returns:
        dict: Dictionary mapping sequence IDs to sequences
    """
    sequences = {}
    with open(file_path, 'r') as file:
        current_id = None
        for line in file:
            line = line.strip()
            if line.startswith(">"):  # Sequence ID line
                current_id = line[1:]  # Remove ">" to extract ID
                sequences[current_id] = ""  # Initialize sequence for this ID as empty
            else:
                if current_id:
                    sequences[current_id] += line  # Append sequence data
    return sequences

def get_all_ids_and_lengths(input_dir):
    """
    Get the set of all IDs across files and the sequence length in each file
    
    Args:
        input_dir (str): Directory containing FASTA files
    
    Returns:
        tuple: (sorted list of all IDs, dictionary of file lengths)
    """
    all_ids = set()
    file_lengths = {}

    for file_name in os.listdir(input_dir):
        if not file_name.endswith(('.fasta', '.fa', '.fas', '.fna', '.faa')):
            continue  # Skip non-FASTA files
            
        file_path = os.path.join(input_dir, file_name)
        if os.path.isfile(file_path):
            sequences = read_fasta(file_path)
            all_ids.update(sequences.keys())  # Collect all IDs
            if sequences:
                # Assume equal length sequences in a file (proper alignment)
                file_lengths[file_name] = len(next(iter(sequences.values())))
            else:
                file_lengths[file_name] = 0  # Length is 0 for empty files

    # Sort IDs: numerically if possible, otherwise alphabetically
    def sort_key(x):
        try:
            return (0, int(x))  # Numeric IDs come first
        except ValueError:
            return (1, x)  # Then alphabetic IDs
    
    return sorted(all_ids, key=sort_key), file_lengths

def concatenate_sequences_stream(input_dir, output_file):
    """
    Process sequences file-by-file, concatenate sequences with the same ID,
    and write to file to reduce memory usage
    
    Args:
        input_dir (str): Directory containing input FASTA files
        output_file (str): Path to output concatenated FASTA file
    
    Returns:
        tuple: (number of sequences, total alignment length)
    """
    # Get all IDs and per-file sequence lengths
    all_ids, file_lengths = get_all_ids_and_lengths(input_dir)
    
    if not all_ids:
        print("Warning: No sequences found in input directory!")
        return 0, 0
    
    print(f"Found {len(all_ids)} unique sequence IDs across {len(file_lengths)} files")

    # Open output file
    with open(output_file, 'w') as out_file:
        # Initialize concatenation result for each ID as an empty string
        concatenated_sequences = {seq_id: "" for seq_id in all_ids}

        # Process files in sorted order for reproducibility
        sorted_files = sorted([f for f in os.listdir(input_dir) 
                              if f.endswith(('.fasta', '.fa', '.fas', '.fna', '.faa'))])
        
        for file_name in sorted_files:
            file_path = os.path.join(input_dir, file_name)
            if os.path.isfile(file_path):
                # Read sequences from the current file
                sequences = read_fasta(file_path)
                seq_length = file_lengths.get(file_name, 0)  # Sequence length in the current file
                
                print(f"Processing {file_name}: {len(sequences)} sequences, length {seq_length}")

                # Update concatenated result for each ID
                for seq_id in all_ids:
                    if seq_id in sequences:
                        concatenated_sequences[seq_id] += sequences[seq_id]
                    else:
                        # Fill with '-' when sequence is missing in this file
                        concatenated_sequences[seq_id] += "-" * seq_length

        # Write final results to file
        total_length = 0
        for seq_id, concatenated_seq in concatenated_sequences.items():
            out_file.write(f">{seq_id}\n")
            # Write 60 characters per line for readability
            for i in range(0, len(concatenated_seq), 60):
                out_file.write(concatenated_seq[i:i+60] + '\n')
            if total_length == 0:
                total_length = len(concatenated_seq)
        
        return len(all_ids), total_length

def parse_arguments():
    """
    Parse command line arguments
    
    Returns:
        argparse.Namespace: Parsed arguments
    """
    parser = argparse.ArgumentParser(
        description='Concatenate multiple sequence alignments into a supermatrix',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example:
    python merge_alignment.py -i gene_alignments/ -o supermatrix.fasta
    
This will concatenate all FASTA files in gene_alignments/ directory
and save the result to supermatrix.fasta
        """
    )
    
    parser.add_argument('-i', '--input', 
                       help='Input directory containing FASTA alignment files',
                       default=None)
    parser.add_argument('-o', '--output', 
                       help='Output concatenated FASTA file',
                       default='concatenated_sequences.fasta')
    
    return parser.parse_args()

if __name__ == "__main__":
    # Parse command line arguments
    args = parse_arguments()
    
    # Use command line arguments if provided, otherwise use default paths
    if args.input:
        input_directory = args.input
        output_file_path = args.output
    else:
        # Default paths - edit these if not using command line arguments
        input_directory = "path_to_your_sequence_files"  # Change this to your input directory
        output_file_path = "concatenated_sequences.fasta"  # Change this to your desired output
        
        print("No command line arguments provided. Using paths defined in script.")
        print("Edit the paths in __main__ section or use command line arguments:")
        print("  python merge_alignment.py -i <input_dir> -o <output.fasta>")
    
    # Check if input directory exists
    if not os.path.exists(input_directory):
        print(f"Error: Input directory '{input_directory}' does not exist!")
        sys.exit(1)
    
    if not os.path.isdir(input_directory):
        print(f"Error: '{input_directory}' is not a directory!")
        sys.exit(1)
    
    # Execute sequence concatenation
    print(f"Input directory: {input_directory}")
    print(f"Output file: {output_file_path}")
    print("Starting concatenation...")
    
    num_seqs, total_length = concatenate_sequences_stream(input_directory, output_file_path)
    
    if num_seqs > 0:
        print(f"\n✓ Concatenation completed successfully!")
        print(f"  - Number of sequences: {num_seqs}")
        print(f"  - Total alignment length: {total_length} bp")
        print(f"  - Results saved to: {output_file_path}")
    else:
        print("\n✗ No sequences were concatenated. Please check your input files.")