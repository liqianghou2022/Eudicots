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


# === Step 3: Construct gene trees with IQ-TREE ===
for f in *.best.fas; do
    # Build maximum-likelihood tree with automatic model selection
    # 1000 ultrafast bootstrap replicates, quiet mode, redo if needed, auto threads
    echo "iqtree3 -s $f -pre ${f}.tre -bb 1000 -m MFP -quiet -redo -T AUTO"
done > 03.iqtree-00.sh

# Run all tree reconstruction jobs in parallel (40 cores)
parallel -j 40 < 03.iqtree-00.sh


# === Step 4: Get all gene trees
cat *.treefile > 144-sp.nwk


# === Remove leaves (tips) exhibiting bootstrap support values less than 50 in the gene trees.
nw_ed 144-sp.nwk 'i & b<=50' o > 144-sp-BS50.nwk

# === Potential long-branch attraction artifacts due to sampling biases were removed with TreeShrink
run_treeshrink.py -t 144-sp-BS50.nwk -o 144-sp-BS50-ts.nwk -m per-gene -c

