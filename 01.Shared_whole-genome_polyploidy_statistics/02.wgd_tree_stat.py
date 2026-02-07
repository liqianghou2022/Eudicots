#!/usr/bin/env python3
"""
Script: 04.wgd_tree_stat.py
Author: Liqiang Hou
Date: 2024.12

Description:
    Analyze gene trees to determine whether two species share a Whole Genome Duplication (WGD) event
    or underwent independent WGD events. The script examines the monophyly of gene copies from each
    species to infer the evolutionary scenario.

Usage:
    python 04.wgd_tree_stat.py
    
    Run this script in a directory containing .nwk files (Newick format gene trees).
    Edit the A and B sets below to match your species' gene copy names.

Input:
    - Multiple .nwk files containing gene trees in Newick format
    - Each file can contain one or more trees separated by semicolons
    - Each species should have 2 copies per gene (from WGD)

Output:
    - WGD_support_summary.txt: Summary table with support statistics for each file

Interpretation:
    - Independent WGD: Both species' copies form monophyletic groups
                       Expected pattern: ((A1,A2),(B1,B2))
                       Indicates WGD occurred after speciation
    
    - Shared WGD: Neither species' copies form monophyletic groups
                  Expected pattern: ((A1,B1),(A2,B2))
                  Indicates WGD occurred before speciation
    
    - Uncertain: Only one species forms monophyletic group
                 May indicate gene loss or complex evolutionary history
    
    Decision criteria:
    - Shared_Ratio > 0.5: Evidence supports shared WGD event
    - Ind_Ratio > 0.5: Evidence supports independent WGD events
    - High Uncertain ratio: Results inconclusive, needs further investigation

Requirements:
    - Python 3.6+
    - ete3 library (pip install ete3)

Example gene naming:
    A = {"1", "2"}  # Arabidopsis thaliana copies
    B = {"3", "4"}  # Oryza sativa copies
"""

from ete3 import Tree
import glob

# === Parameter Settings ===
# Define gene copy sets for each species (2 copies per species from WGD)
A = {"1", "2"}  # Replace with species A gene copy names (e.g., {"speciesA_copy1", "speciesA_copy2"})
B = {"3", "4"}  # Replace with species B gene copy names (e.g., {"speciesB_copy1", "speciesB_copy2"})

# Output file name
output_file = "WGD_support_summary.txt"

# Process all .nwk files and write results
with open(output_file, "w") as out:
    # Write header
    out.write("File\tTotal\tIndependent\tShared\tUncertain\tInd_Ratio\tShared_Ratio\n")
    
    # Process each .nwk file in current directory
    for nwk_file in glob.glob("*.nwk"):
        total = independent = shared = uncertain = 0
        
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
                
                # Get leaf names and find species copies present in this tree
                leaves = set(l.name for l in t.iter_leaves())
                a_here = list(A & leaves)  # Species A copies in this tree
                b_here = list(B & leaves)  # Species B copies in this tree
                
                # Skip incomplete trees (need both copies for meaningful analysis)
                if len(a_here) < 2 or len(b_here) < 2:
                    continue
                
                total += 1
                
                # Check monophyly for each species
                monoA = t.check_monophyly(values=a_here, target_attr="name")[0]
                monoB = t.check_monophyly(values=b_here, target_attr="name")[0]
                
                # Classify tree topology based on monophyly patterns
                if monoA and monoB:
                    # Both species form monophyletic groups -> independent WGD
                    # Pattern: ((A1,A2),(B1,B2))
                    independent += 1
                elif not monoA and not monoB:
                    # Neither species forms monophyletic group -> shared WGD
                    # Pattern: ((A1,B1),(A2,B2))
                    shared += 1
                else:
                    # Mixed pattern -> uncertain (could be gene loss or other factors)
                    uncertain += 1
            except:
                # Skip malformed trees
                continue
        
        # Calculate ratios for this file
        if total > 0:
            ind_ratio = independent / total
            shared_ratio = shared / total
            out.write(f"{nwk_file}\t{total}\t{independent}\t{shared}\t{uncertain}\t{ind_ratio:.4f}\t{shared_ratio:.4f}\n")

# Print completion message and interpretation guide
print(f"Analysis completed. Results saved to {output_file}")
print("\nInterpretation guide:")
print("  - Independent: Both species form monophyletic groups (WGD after speciation)")
print("  - Shared: Neither species forms monophyletic groups (WGD before speciation)")
print("  - Uncertain: Mixed pattern (may indicate gene loss)")
print("\nDecision criteria:")
print("  - Shared_Ratio > 0.5: Evidence supports shared WGD")
print("  - Ind_Ratio > 0.5: Evidence supports independent WGD")