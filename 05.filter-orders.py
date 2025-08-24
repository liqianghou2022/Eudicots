#!/usr/bin/env python3
"""
Script: 05.filter-orders
Author: Liqiang Hou
Date: 2025.08

Description:
    Filter phylogenetic trees based on taxonomic group coverage. This script reads
    a species-to-group mapping from a CSV file and filters trees to keep only those
    that contain species from at least a specified number of different groups.
    Branch lengths are preserved in the output.

Usage:
    python 05.filter-orders -i input.nwk -c mapping.csv -o output.nwk [-g 30]

Arguments:
    -i, --input  : Input Newick file containing one or more phylogenetic trees
    -c, --csv    : CSV mapping file with species-to-group assignments (no header)
    -o, --output : Output Newick file for filtered trees
    -g, --groups : Minimum number of groups required (default: 30)

Input Files:
    1. Newick file: Standard Newick format with one or more trees
    2. CSV mapping file: Two columns without header
       Column 1: Species/sequence ID
       Column 2: Group/category name
       
       Example:
       species1,GroupA
       species2,GroupA
       species3,GroupB
       species4,GroupC

Output:
    - Filtered Newick file containing only trees that meet the group coverage criterion
    - Original branch lengths and tree topology are preserved
    - Console output showing number of trees retained

Example:
    # Keep trees with at least 30 different groups represented
    python filter_trees_by_groups.py -i trees.nwk -c species_groups.csv -o filtered.nwk -g 30
    
    # Keep trees with at least 10 different groups
    python filter_trees_by_groups.py -i trees.nwk -c species_groups.csv -o filtered.nwk -g 10

Filtering Logic:
    - A tree passes if it contains species from >= specified number of groups
    - A group is considered "present" if at least one species from that group is in the tree
    - Trees with fewer groups than the threshold are discarded

Requirements:
    - Python 3.6+
    - Biopython (pip install biopython)
    - pandas (pip install pandas)

Notes:
    - Species IDs in the tree must match exactly with IDs in the CSV file
    - Group names are case-sensitive
    - Empty trees or trees with no matching species are automatically filtered out
    - This is useful for ensuring phylogenetic diversity in tree selection
"""

from Bio import Phylo
import pandas as pd
import argparse

def load_groups_from_csv(csv_file):
    """
    Load species-to-group mapping from a CSV file.
    :param csv_file: Path to CSV file (two columns: id, group)
    :return: Dict mapping group -> list of species IDs (as strings)
    """
    df = pd.read_csv(csv_file, header=None, names=['id', 'group'])
    groups = df.groupby('group')['id'].apply(lambda x: list(map(str, x))).to_dict()
    return groups

def check_tree_for_groups(tree, groups, required_groups_count=10):
    """
    Check whether the tree contains at least a specified number of groups.
    :param tree: Phylogenetic tree
    :param groups: Dict of group -> list of species IDs
    :param required_groups_count: Minimum number of groups that must be present
    :return: True if condition is met, else False
    """
    terminal_names = {term.name for term in tree.get_terminals()}
    groups_covered = 0

    for group_nodes in groups.values():
        if any(node in terminal_names for node in group_nodes):
            groups_covered += 1
        if groups_covered >= required_groups_count:
            return True
    return False

def filter_trees_by_groups(file_path, groups, output_file, required_groups_count=10):
    """
    Filter trees by group coverage, preserving original branch lengths.
    :param file_path: Input Newick file path
    :param groups: Dict of group -> list of species IDs
    :param output_file: Output Newick file path
    :param required_groups_count: Minimum number of groups required
    :return: Number of trees that meet the condition
    """
    with open(file_path, 'r') as f:
        trees = list(Phylo.parse(f, 'newick'))  # Read all trees

    # Keep trees that satisfy the group coverage criterion (branch lengths preserved)
    valid_trees = [tree for tree in trees if check_tree_for_groups(tree, groups, required_groups_count)]

    # Write filtered trees to output, preserving branch lengths
    with open(output_file, 'w') as out_f:
        Phylo.write(valid_trees, out_f, 'newick')

    return len(valid_trees)

# Main
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter Newick trees by the number of groups present.")
    parser.add_argument("-i", "--input", required=True, help="Input Newick file")
    parser.add_argument("-c", "--csv", required=True, help="CSV mapping file: two columns [id, group], no header")
    parser.add_argument("-o", "--output", required=True, help="Output Newick file")
    parser.add_argument("-g", "--groups", type=int, default=30, help="Minimum number of groups required (default: 30)")
    args = parser.parse_args()

    csv_file = args.csv
    newick_file = args.input
    output_file = args.output
    required_groups = args.groups

    # Load group mapping
    groups = load_groups_from_csv(csv_file)

    # Filter trees by group coverage (branch lengths preserved)
    count = filter_trees_by_groups(newick_file, groups, output_file, required_groups)

    # Report
    print(f"Number of trees after filtering: {count}")
    print(f"Filtered trees saved to '{output_file}'.")