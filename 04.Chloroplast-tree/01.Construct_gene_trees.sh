Script: Gene Tree Construction for Syntenic Genes

Workflow Description
--------------------
All syntenic genes were used to construct individual gene trees through a three-step pipeline:

1. Sequence Alignment
   - Align sequences for each gene using MAFFT with automatic strategy selection.
   - Parallel execution accelerates alignment for all FASTA files.

2. Alignment Trimming
   - Remove low-information or gap-rich regions using trimAl (automated mode).
   - Retains only well-aligned, informative sites for phylogenetic inference.

3. Gene Tree Reconstruction
   - Infer maximum-likelihood gene trees using IQ-TREE.
   - Automatic model selection (-m MFP) and 1000 ultrafast bootstrap replicates (-bb 1000).
   - Parallel execution for efficiency.

# === Step 1: Multiple sequence alignment with MAFFT ===
for f in *.fa; do
    # Align each FASTA file quietly using MAFFT --auto
    echo "mafft --auto --quiet $f > $f.fas"
done > 01.mafft-00.sh

# Run all alignment jobs in parallel (40 cores)
parallel -j 40 < 01.mafft-00.sh


# === Step 2: Trim low-quality regions with trimAl ===
for f in *.fas; do
    # Automatically remove gap-rich or low-information regions
    echo "trimal -in $f -out $f.best.fas -automated1"
done > 02.trimal-00.sh

# Run all trimming jobs in parallel (40 cores)
parallel -j 40 < 02.trimal-00.sh


# Obtain the concatenated sequence
python 01.chuanlian.py


#  Build the tree using iqtree
/home/houlq2025/software/iqtree-3.0.1-Linux-intel/bin/iqtree3 -s *.fas -pre supermatrix.triml -bb 1000 -m MFP -quiet -redo -T AUTO
