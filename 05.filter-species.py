#!/usr/bin/env python3
"""
Script: 05.filter_trees_by_species.py
Author: Liqiang Hou
Date: 2025.08

Description:
    Filter phylogenetic trees based on the minimum number of leaf nodes (species/taxa).
    This script reads Newick format trees and keeps only those with at least a 
    specified number of terminal nodes.

Usage:
    python 05.filter_trees_by_species.py -i input.nwk -o output.nwk -t 10

Arguments:
    -i, --input     : Input file containing Newick trees (one tree per line)
    -o, --output    : Output file for filtered trees
    -t, --threshold : Minimum number of leaf nodes (species) required

Input:
    - Newick format tree file
    - Each line should contain one complete tree
    - Trees can have branch lengths and support values

Output:
    - Filtered Newick file containing only trees with >= threshold species
    - Original tree format and branch lengths are preserved

Example:
    # Keep only trees with at least 10 species
    python 05.filter_trees_by_species.py -i all_trees.nwk -o filtered_trees.nwk -t 10
    
    # Keep only trees with at least 50 species
    python 05.filter_trees_by_species.py -i gene_trees.nwk -o large_trees.nwk -t 50

Tree Examples:
    Input trees:
        ((A,B),(C,D));              # 4 species - kept if threshold <= 4
        ((A,B),C);                  # 3 species - kept if threshold <= 3
        (A,B);                      # 2 species - kept if threshold <= 2

Filtering Logic:
    - Counts the number of leaf nodes (terminal taxa) in each tree
    - Trees with number_of_leaves >= threshold are kept
    - Trees with fewer leaves are discarded
    - Empty lines are automatically skipped

Requirements:
    - Python 3.6+
    - ete3 (pip install ete3)

Notes:
    - This filter is useful for removing gene trees with too few taxa
    - Small trees may not be informative for phylogenetic analyses
    - Trees that cannot be parsed are skipped with a warning message
    - The script preserves the original Newick format including branch lengths
"""

from ete3 import Tree
import argparse

def main():
    parser = argparse.ArgumentParser(description="Filter Newick trees by minimum number of leaf species.")
    parser.add_argument("-i", "--input", required=True, help="Input file with Newick trees (one per line)")
    parser.add_argument("-o", "--output", required=True, help="Output file for filtered trees")
    parser.add_argument("-t", "--threshold", type=int, required=True, help="Minimum leaf (species) count threshold")
    args = parser.parse_args()

    input_file = args.input
    output_filtered_file = args.output
    minimum_species_threshold = args.threshold

    # Read all gene trees
    with open(input_file, "r") as infile:
        trees = infile.readlines()

    # Initialize filtered results
    filtered_trees = []

    # Iterate over each gene tree
    for newick in trees:
        newick = newick.strip()
        if not newick:  # Skip empty lines
            continue

        try:
            # Parse tree
            tree = Tree(newick, format=1)  # Read Newick-format tree
            # Get number of leaves (species)
            num_species = len(tree.get_leaf_names())

            # Keep trees meeting the threshold
            if num_species >= minimum_species_threshold:
                filtered_trees.append(newick + "\n")

        except Exception as e:
            print(f"Failed to parse tree: {newick[:30]}... Error: {e}")

    # Write filtered results to file
    with open(output_filtered_file, "w") as outfile:
        outfile.writelines(filtered_trees)

    print(f"Filtering completed. Kept {len(filtered_trees)} trees. Results written to {output_filtered_file}")

if __name__ == "__main__":
    main()