#!/usr/bin/env python3
"""
Script: 04.filter_trees_by_support.py
Author: Liqiang Hou
Date: 2025.08

Description:
    Filter phylogenetic trees in Newick format based on bootstrap support values 
    and branch length thresholds. Trees with any support value or branch length 
    below the specified thresholds will be filtered out.

Usage:
    python 04.filter_trees_by_support.py -i input.nwk -o output.nwk [-t 0.7] [-b 0.01]

Arguments:
    -i, --infile   : Input Newick file containing one or more trees
    -o, --outfile  : Output file for filtered trees
    -t, --threshold_support : Minimum support value threshold (default: 0.7)
    -b, --threshold_branch  : Minimum branch length threshold (default: 0.01)

Input:
    - Newick format tree file (.nwk)
    - Each line should contain one tree
    - Support values should appear after closing parentheses and before colons
    - Branch lengths should appear after colons

Output:
    - Filtered Newick file containing only trees that pass both thresholds
    - Summary statistics printed to console

Example:
    # Filter trees keeping only those with all support values >= 0.8 and branch lengths >= 0.02
    python 04.filter_trees_by_support.py -i trees.nwk -o filtered_trees.nwk -t 0.8 -b 0.02
    
    # Use default thresholds (support >= 0.7, branch length >= 0.01)
    python 04.filter_trees_by_support.py -i trees.nwk -o filtered_trees.nwk

Tree Format Example:
    ((A:0.1,B:0.2)0.95:0.05,(C:0.15,D:0.18)0.88:0.03)1.0:0.0;
    
    Where:
    - 0.95, 0.88, 1.0 are support values
    - 0.05, 0.03, 0.0 are internal branch lengths

Filtering Logic:
    - A tree passes if ALL support values >= threshold_support
    - AND ALL internal branch lengths >= threshold_branch
    - Trees missing support values or branch lengths are discarded

Requirements:
    - Python 3.6+
    - No external dependencies

Notes:
    - Support values can be bootstrap values (0-100) or posterior probabilities (0-1)
    - The script checks internal node support and branch lengths only
    - Terminal branch lengths are not considered in filtering
"""

import argparse
import re

def extract_supports(newick):
    """Extract support values (numbers after internal node parentheses, right before the colon)"""
    return [float(s) for s in re.findall(r'\)\s*([\d\.eE+-]+)\s*:', newick)]

def extract_branch_lengths(newick):
    """Extract all internal node branch lengths (numbers after the colon following the parentheses)"""
    # This regex matches the branch_length in the pattern `)support:branch_length`
    return [float(s) for s in re.findall(r'\)\s*[\d\.eE+-]+\s*:\s*([\d\.eE+-]+)', newick)]

def main():
    parser = argparse.ArgumentParser(description="Filter trees by support and branch length.")
    parser.add_argument("-i", "--infile", required=True, help="Input nwk file")
    parser.add_argument("-o", "--outfile", required=True, help="Output filtered nwk file")
    parser.add_argument("-t", "--threshold_support", type=float, default=0.7, help="Support threshold, default 0.7")
    parser.add_argument("-b", "--threshold_branch", type=float, default=0.01, help="Branch length threshold, default 0.01")
    args = parser.parse_args()

    n_total = 0
    n_passed = 0

    with open(args.infile) as fin, open(args.outfile, "w") as fout:
        for line in fin:
            line = line.strip()
            if not line:
                continue
            n_total += 1

            supports = extract_supports(line)
            branches = extract_branch_lengths(line)

            # Filter out if any condition is not met
            if not supports or not branches:
                continue  # No support or branch length; discard

            if min(supports) >= args.threshold_support and min(branches) >= args.threshold_branch:
                fout.write(line + "\n")
                n_passed += 1

    print(f"Total trees:    {n_total}")
    print(f"Trees passed:   {n_passed}")
    print(f"Support cutoff: {args.threshold_support}")
    print(f"Branch length cutoff: {args.threshold_branch}")

if __name__ == "__main__":
    main()