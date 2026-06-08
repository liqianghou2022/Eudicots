#!/bin/bash

# Set the directory containing gene tree files
input_dir="./cds-fa"

# Set the output file name
output_file="combine.tsv"

# Empty the output file if it already exists, or create a new one
> "$output_file"

# Iterate over all .treefile files
for tree_file in "$input_dir"/*.treefile; do
    # Get the file name without the directory path
    filename=$(basename "$tree_file")

    # Optionally remove the extension to use the file name as the gene ID
    id="${filename%.treefile}"

    # Read the tree content; assumes each file contains one Newick tree on a single line
    tree=$(cat "$tree_file")

    # Write the gene ID and tree to the combined output file
    echo -e "${id}\t${tree}" >> "$output_file"
done