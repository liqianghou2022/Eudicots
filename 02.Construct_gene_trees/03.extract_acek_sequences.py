"""
Script: Extract ACEK-corresponding sequences from species-specific sequences

Description
-----------
Based on the filtered syntenic gene matrix ("Species-one-gene.csv") and the 
sequence files of all species (either CDS or protein sequences), this script 
retrieves sequences corresponding to ACEK genes from all species and saves them 
as individual FASTA files.

For each ACEK gene:
- A FASTA file is created named after the first gene in the corresponding row.
- Each sequence within the file corresponds to a species, with the sequence ID 
  replaced by the species/column identifier from the matrix.

Input
-----
1. Filtered syntenic gene matrix
   - File: Species-one-gene.csv
   - No header
   - Rows: ACEK genes
   - Columns: one gene per species (selected representative after filtering)

2. Species sequences
   - File: all.clean.cds (FASTA format, can also be protein sequences)
   - Sequence IDs correspond to gene names used in the CSV

Output
------
- One FASTA file per ACEK gene
- File name: ACEK gene ID + ".fa"
- Each FASTA file contains sequences of the corresponding genes across all species
- Sequence IDs within the FASTA are replaced by the column name (species identifier)

Note
----
If a gene in the CSV is not found in the sequence file, a warning is printed.
"""
import os
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# --- 1. Read CSV file (no headers) ---
df = pd.read_csv("Species-one-gene.csv", header=None, dtype=str)
# Assign column names as 1, 2, ..., n
df.columns = [str(i+1) for i in range(df.shape[1])]

# --- 2. Load all sequences into a dictionary ---
fasta_file = "all.clean.cds"
seq_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

# --- 3. Set up output directory ---
out_dir = "output_fasta-fix"
os.makedirs(out_dir, exist_ok=True)

# --- 4. Iterate over rows and extract sequences ---
for i, row in df.iterrows():
    genes = row.dropna()  # Non-empty genes
    if genes.empty:
        continue

    # Output FASTA file = first ACEK gene ID + ".fa"
    out_fa = os.path.join(out_dir, f"{genes.iloc[0]}.fa")

    records = []
    for col_name, g in genes.items():
        if g in seq_dict:
            # Create SeqRecord with ID replaced by column name (species ID)
            new_record = SeqRecord(seq_dict[g].seq, id=col_name, description="")
            records.append(new_record)
        else:
            print(f"Warning: gene {g} not found in fasta")

    if records:
        SeqIO.write(records, out_fa, "fasta")

print(f"Extraction completed! All files are saved in: {out_dir}")
