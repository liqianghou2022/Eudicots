#!/usr/bin/env python3
"""
Script: 04.wgt_tree_stat.py
Author: Liqiang Hou
Date: 2024.12

Description:
    Analyze gene trees to determine whether two species share a Whole Genome Triplication (WGT) event
    or underwent independent WGT events. The script examines the monophyly of gene copies from each
    species to infer the evolutionary scenario.

Usage:
    python 04.wgt_tree_stat.py
    
    Run this script in a directory containing .nwk files (Newick format gene trees).
    Edit the A and B sets below to match your species' gene copy names.

Input:
    - Multiple .nwk files containing gene trees in Newick format
    - Each file can contain one or more trees separated by semicolons

Output:
    - WGT_support_summary.txt: Summary table with support statistics for each file

Interpretation:
    - Shared WGT: Both species' copies are intermingled (neither forms monophyletic group)
                  Expected pattern: ((A1,B1),(A2,B2),(A3,B3))
    - Independent WGT: At least one species' copies form monophyletic group
                       Expected pattern: ((A1,A2,A3),(B1,B2,B3))
    - Higher Shared_Ratio (>0.5) suggests shared WGT event
    - Higher NonShared_Ratio (>0.5) suggests independent WGT events

Requirements:
    - Python 3.6+
    - ete3 library (pip install ete3)
    
Example gene naming:
    A = {"1", "2", "3"}  # Arabidopsis thaliana copies
    B = {"4", "5", "6"}  # Oryza sativa copies
"""

from ete3 import Tree
import glob

# === Parameter Settings ===
# Define gene copy sets for each species (3 copies per species from WGT)
A = {"1", "2", "3"}  # Replace with species A gene copy names (e.g., {"AtA1", "AtA2", "AtA3"})
B = {"4", "5", "6"}  # Replace with species B gene copy names (e.g., {"OsB1", "OsB2", "OsB3"})

# Output file name
output_file = "WGT_support_summary.txt"

# Process all .nwk files and write results
with open(output_file, "w") as out:
    # Write header
    out.write("File\tTotal_Trees\tNonShared_Count\tNonShared_Ratio\tShared_Count\tShared_Ratio\n")

    # Process each .nwk file in current directory
    for nwk_file in glob.glob("*.nwk"):
        total = 0
        support_shared = 0      # Supports shared WGT (both A and B are non-monophyletic)
        support_non_shared = 0  # Does not support shared WGT (A or B is monophyletic)

        # Read and parse trees from file
        with open(nwk_file) as f:
            nwk = f.read().strip()
            # Ensure proper formatting
            if not nwk.endswith(";"):
                nwk += ";"
            # Split multiple trees if present
            trees = [x.strip() + ";" for x in nwk.split(";") if x.strip()]

        # Analyze each tree
        for t_str in trees:
            try:
                # Parse tree with format 1 (flexible with support values)
                t = Tree(t_str, format=1)
            except:
                # Skip malformed trees
                continue

            total += 1
            
            # Get leaf names and find species copies present in this tree
            leaves = set(l.name for l in t.iter_leaves())
            a_here = list(A & leaves)  # Species A copies in this tree
            b_here = list(B & leaves)  # Species B copies in this tree

            # Check monophyly for each species
            # If only one copy present, consider as monophyletic (True)
            monoA = t.check_monophyly(values=a_here, target_attr="name")[0] if len(a_here) > 1 else True
            monoB = t.check_monophyly(values=b_here, target_attr="name")[0] if len(b_here) > 1 else True

            # Classify tree topology
            if monoA or monoB:
                # At least one species is monophyletic -> supports independent WGT
                support_non_shared += 1
            else:
                # Both species are non-monophyletic -> supports shared WGT
                support_shared += 1

        # Calculate ratios
        if total > 0:
            ratio_non_shared = support_non_shared / total
            ratio_shared = support_shared / total
        else:
            ratio_non_shared = ratio_shared = 0

        # Write results for this file
        out.write(f"{nwk_file}\t{total}\t{support_non_shared}\t{ratio_non_shared:.4f}\t{support_shared}\t{ratio_shared:.4f}\n")

# Print completion message
print(f"Analysis completed. Results saved to {output_file}")
print(f"Interpretation: Higher Shared_Ratio (>0.5) suggests shared WGT event")
print(f"                Higher NonShared_Ratio (>0.5) suggests independent WGT events")