#!/usr/bin/env python3
"""
Script: 03.extract_sequences_from_csv.py
Author: Liqiang Hou
Date: 2025.08

Description:
    Extract sequences from a FASTA file based on gene IDs in a CSV file.
    Each CSV row becomes a separate FASTA file with sequences numbered by column position.

Usage:
    python 03.extract_sequences_from_csv.py
    
    Edit pep_file and csv_file variables below before running.

Input:
    - all.pep: Master FASTA file with all sequences
    - CSV file: Each row contains gene IDs (comma-separated)

Output:
    - Creates directories 1/, 2/, 3/... for each row
    - Each contains a FASTA file with sequences renamed to column numbers

Example:
    CSV row: gene1,gene2,gene3
    Output: 1/1.fa containing sequences >1, >2, >3

Requirements:
    - Biopython (pip install biopython)
"""

import csv
import os
from Bio import SeqIO

# === EDIT THESE PATHS ===
pep_file = "all.pep"        # Your FASTA file
csv_file = "your_file.csv"  # Your CSV file with gene IDs
# ========================

# Read all sequences into memory
sequences = {record.id: record for record in SeqIO.parse(pep_file, "fasta")}
print(f"Loaded {len(sequences)} sequences")

# Process CSV file row by row
with open(csv_file, newline='', encoding='utf-8') as csvfile:
    reader = csv.reader(csvfile)
    
    for row_index, row in enumerate(reader, start=1):
        output_seqs = []
        
        # Extract sequences for this row
        for col_index, gene_id in enumerate(row, start=1):
            gene_id = gene_id.strip()
            if gene_id and gene_id in sequences:
                seq_record = sequences[gene_id].copy()
                seq_record.id = str(col_index)  # Rename to column number
                seq_record.description = ""
                output_seqs.append(seq_record)
        
        # Write output if sequences found
        if output_seqs:
            dir_name = f"{row_index}"
            os.makedirs(dir_name, exist_ok=True)
            
            output_file = os.path.join(dir_name, f"{row_index}.fa")
            SeqIO.write(output_seqs, output_file, "fasta")
            print(f"Row {row_index}: {len(output_seqs)} sequences -> {output_file}")

print("Extraction completed!")