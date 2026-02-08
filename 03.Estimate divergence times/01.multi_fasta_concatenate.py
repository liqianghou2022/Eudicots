#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os

def read_fasta(file_path):
    """
    Read a FASTA file and return a dictionary {ID: sequence}
    """
    sequences = {}
    with open(file_path, 'r') as file:
        current_id = None
        for line in file:
            line = line.strip()
            if line.startswith(">"):  # Sequence ID line
                current_id = line[1:]  # Remove ">" to extract ID
                sequences[current_id] = ""  # Initialize sequence for this ID
            else:
                if current_id:
                    sequences[current_id] += line  # Append sequence data
    return sequences

def get_all_ids_and_lengths(input_dir):
    """
    Collect all sequence IDs across files and record the sequence length of each file
    """
    all_ids = set()
    file_lengths = {}

    for file_name in os.listdir(input_dir):
        file_path = os.path.join(input_dir, file_name)
        if os.path.isfile(file_path):
            sequences = read_fasta(file_path)
            all_ids.update(sequences.keys())  # Collect all IDs
            if sequences:
                file_lengths[file_name] = len(next(iter(sequences.values())))  # Assume equal length per file
            else:
                file_lengths[file_name] = 0  # Empty file has length 0

    return sorted(all_ids, key=lambda x: int(x)), file_lengths

def concatenate_sequences_stream(input_dir, output_file):
    """
    Process files sequentially and concatenate sequences with the same ID,
    writing results to file to reduce memory usage
    """
    # Get all IDs and sequence lengths for each file
    all_ids, file_lengths = get_all_ids_and_lengths(input_dir)

    # Open output file
    with open(output_file, 'w') as out_file:
        # Initialize concatenated sequence for each ID
        concatenated_sequences = {seq_id: "" for seq_id in all_ids}

        # Process files in order
        for file_name in os.listdir(input_dir):
            file_path = os.path.join(input_dir, file_name)
            if os.path.isfile(file_path):
                # Read sequences from current file
                sequences = read_fasta(file_path)
                seq_length = file_lengths[file_name]  # Sequence length in this file

                # Update concatenated sequences for each ID
                for seq_id in all_ids:
                    if seq_id in sequences:
                        concatenated_sequences[seq_id] += sequences[seq_id]
                    else:
                        concatenated_sequences[seq_id] += "-" * seq_length  # Fill missing with '-'

        # Write final concatenated sequences to file
        for seq_id, concatenated_seq in concatenated_sequences.items():
            out_file.write(f">{seq_id}\n")
            # Write 60 characters per line
            for i in range(0, len(concatenated_seq), 60):
                out_file.write(concatenated_seq[i:i+60] + '\n')

if __name__ == "__main__":
    # Input directory containing all sequence files
    input_directory = "./00.trimal"
    # Output result file
    output_file_path = "concatenated_sequences.fasta"

    # Execute sequence concatenation
    concatenate_sequences_stream(input_directory, output_file_path)
    print(f"All sequences have been successfully concatenated. Output saved to {output_file_path}")
