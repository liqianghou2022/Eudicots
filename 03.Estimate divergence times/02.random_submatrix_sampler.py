import random
import os
from Bio import SeqIO

# Input file and parameters
input_fasta = "concatenated_sequences.fasta"
output_dir = "./random-100-2"  # Output directory
subseq_length = 8000000
num_files = 100

# Create output directory if it does not exist
os.makedirs(output_dir, exist_ok=True)

# Read sequences and check that all have the same length
records = list(SeqIO.parse(input_fasta, "fasta"))
seq_lengths = {len(str(record.seq)) for record in records}
if len(seq_lengths) != 1:
    raise ValueError("All sequences must have the same length. Please check the input alignment!")

full_length = seq_lengths.pop()
if full_length < subseq_length:
    raise ValueError("The target subsequence length is greater than the full sequence length!")

for idx in range(1, num_files + 1):
    # Randomly select an interval
    start = random.randint(0, full_length - subseq_length)
    end = start + subseq_length
    output_fasta = os.path.join(output_dir, f"random_{idx:03d}.fasta")
    print(f"[{idx:03d}] Extracted interval: {start+1}-{end}, output file: {output_fasta}")
    with open(output_fasta, 'w') as out_f:
        for record in records:
            subseq = str(record.seq)[start:end]
            out_f.write(f">{record.id}\n")
            out_f.write(subseq + '\n')  # Write sequence in a single line

print(f"All done! {num_files} FASTA files have been written to {output_dir}/")
