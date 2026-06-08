#!/usr/bin/env python3
from Bio import SeqIO
import glob

# All FASTA files, sorted by filename to ensure consistent order
fa_files = sorted(glob.glob("*.fas"))

# Collect all species names (assuming each FASTA has at least one sequence)
species_set = set()
for f in fa_files:
    for rec in SeqIO.parse(f, "fasta"):
        sp = rec.id.split("_")[0]  # Use the first field as the species name
        species_set.add(sp)

species_list = sorted(list(species_set))  # Alphabetical order (or custom order if needed)

# Initialize the super-sequence dictionary for each species
superseq = {sp: "" for sp in species_list}

# Iterate over each FASTA file and concatenate sequences by species
for f in fa_files:
    # Store sequences for each species in the current file
    seq_dict = {}
    for rec in SeqIO.parse(f, "fasta"):
        sp = rec.id.split("_")[0]
        seq_dict[sp] = str(rec.seq).upper()  # Convert to uppercase

    # Sequence length for the current file
    seqlen = max(len(s) for s in seq_dict.values()) if seq_dict else 0

    # Concatenate to the super-sequence, fill missing species with '-'
    for sp in species_list:
        seq = seq_dict.get(sp, "-" * seqlen)
        superseq[sp] += seq

# Output the final supermatrix FASTA
with open("supermatrix.fa", "w") as out:
    for sp in species_list:
        out.write(f">{sp}\n{superseq[sp]}\n")

