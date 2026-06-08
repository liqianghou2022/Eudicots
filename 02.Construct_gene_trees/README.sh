Script: Gene Tree Construction for Syntenic Genes

## Workflow Description

All syntenic genes were used to construct individual gene trees and to infer collinearity-based species trees through a twelve-step pipeline. First, syntenic gene sets corresponding to ACEK genes were identified from cross-species collinearity analyses. Initial gene trees were then constructed to identify abnormally long branches and potential paralogous or erroneous sequences. After long-branch pruning and taxon-coverage filtering, the retained species/gene IDs were extracted for each gene tree. These IDs were subsequently used to retrieve the corresponding CDS and protein sequences from the original filtered gene sets. The cleaned gene sets were then re-aligned, re-trimmed and re-estimated to generate the final gene trees used for downstream species-tree inference.

1. Syntenic Gene Set Construction

   * Syntenic gene sets corresponding to ACEK genes were obtained for each species using the `-a` parameter of WGDI.
   * The resulting WGDI output CSV files were merged into a single syntenic gene matrix using `01.merge_wgdi_csv.py`.
   * The merged syntenic gene matrix was simplified using 02.filter_syntenic_matrix_by_species.py, with a species-to-column mapping file, mapping.csv.
   * This filtering step reduced redundancy caused by recent whole-genome duplication events, following the strategy described in the Materials and Methods and Supplementary Materials.
   * Based on the simplified syntenic matrix, CDS and protein sequence sets for each syntenic gene were extracted using `03.extract_acek_sequences.py`.

2. Protein Sequence Alignment

   * Protein sequences for each syntenic gene were aligned using MAFFT with automatic strategy selection.
   * All alignment jobs were executed in parallel to accelerate the processing of multiple FASTA files.

3. Codon Alignment Generation

   * The pre-aligned protein sequence alignments were used to guide the generation of corresponding CDS codon alignments.
   * This step was performed using `pxaa2cdn` from the phyx package.
   * CDS sequences were aligned according to the protein alignment framework, thereby preserving codon structure.

4. Alignment Trimming

   * Codon alignments were trimmed using `pxclsq` from the phyx package.
   * Low-information or gap-rich regions were removed with a trimming threshold of 0.1.
   * The resulting cleaned alignments retained well-aligned and phylogenetically informative sites for downstream phylogenetic inference.

5. Initial Gene Tree Reconstruction

   * Initial maximum-likelihood gene trees were inferred using IQ-TREE.
   * For each trimmed codon alignment, the best-fit substitution model was selected automatically using ModelFinder (`-m MFP`).
   * Branch support was assessed with 1000 ultrafast bootstrap replicates (`-bb 1000`).
   * Tree reconstruction jobs were executed in parallel for computational efficiency.

6. Gene Tree Collection

   * All resulting gene tree files were collected and merged into a two-column table containing gene IDs and corresponding Newick gene trees.
   * This step was performed by running `06.combined.sh`.

7. Long-Branch Pruning and Removal of Potential Paralogous Sequences

   * The initial gene trees were screened for abnormally long branches using a branch-length threshold of 0.8.
   * This step was performed using `06.prune_long_branch_trees`.
   * Trees without long branches were written to `01.correct_trees.tsv`.
   * Trees containing long branches were written to `02.wrong_trees.tsv`.
   * For trees with long branches, potentially problematic tips or clades were pruned, and the resulting trees were written to `03.pruned_trees.tsv`.
   * A detailed pruning report was generated as `04.pruned_trees.report.tsv`.
   * Trees without long branches and pruned trees were combined into `combine.new.tsv`.
   * This step was used to identify potential paralogous, erroneous or highly divergent sequences.
   * The pruned trees were used only to determine which sequences should be retained or removed, and were not treated as the final gene trees for species-tree inference.

8. Filtering by Order Coverage

   * The long-branch-filtered gene tree set, `combine.new.tsv`, was further filtered according to the number of represented taxonomic orders.
   * Only gene trees containing species from at least 22 orders were retained.
   * This filtering was performed using `08.filter-orders` with the species-to-order mapping file `mapping-5.csv`.
   * The filtered output was saved as `combine-22order.tsv`.

9. Filtering by Species Coverage

   * The order-filtered gene trees were further filtered according to species coverage.
   * Only gene trees containing at least 73 species were retained.
   * This filtering was performed using `09.filter-species`.
   * The filtered output was saved as `combine-22order-73.tsv`.

10. Extraction of Cleaned Gene Sets and Reconstruction of Final Gene Trees

* The retained leaves in each filtered gene tree were extracted using `10.extract_tree_leaves.py`.
* This step generated an `id` directory containing one ID file for each syntenic gene.
* Each ID file recorded the retained species/tip names and their corresponding gene IDs in the gene set.
* These ID files were used to extract the corresponding CDS sequences from `filtered_cds` and protein sequences from `filtered_pep` using `seqkit`.
* The extracted sequences represented the cleaned gene sets after removing long-branch-associated potential paralogs and applying taxon-coverage filters.
* Protein sequence alignment, codon alignment generation, alignment trimming, maximum-likelihood tree reconstruction and gene tree collection were then repeated following Steps 2–6.
* This procedure generated the final rebuilt gene tree set, `combine-rebuilt.tsv`.
* Therefore, the final gene trees used for downstream species-tree inference were re-estimated from cleaned alignments, rather than obtained by simply pruning tips from the initial gene trees.

11. Collinearity-Based Species-Tree Inference

* The rebuilt gene trees were extracted from the two-column table.
* The Newick trees were saved as `combine-rebuilt.nwk` and used as input for ASTRAL 5.6.9.
* A collinearity-based species tree was inferred from the cleaned syntenic gene trees.
* Branch annotations were generated using the `-t 2` option, which reports detailed quartet-based support information.

12. Constrained Species-Tree Inference

* A final merged phylogenetic tree was inferred using ASTRAL 5.6.9 under a constraint-tree framework.
* The constraint tree, `species_constraint.tre`, was constructed primarily based on shared chromosomal fusion events and shared whole-genome polyploidization events.
* The cleaned collinearity-based gene trees were used as supplementary phylogenetic evidence under the constraint imposed by the karyotype-evolution framework.
* The `-j` parameter was used to enforce the constraint tree during species-tree inference.
* The `-t 2` option was used to output detailed quartet-based branch annotations.
Script: Gene Tree Construction for Syntenic Genes


# === Step 1: Construction of syntenic gene sets ===
# First, syntenic gene sets corresponding to ACEK genes were obtained for each species
# using the "-a" parameter of WGDI:
# https://github.com/SunPengChuan/wgdi

# Merge all WGDI output CSV files into a single syntenic gene matrix，命名为merge.all.csv
python 01.merge_wgdi_csv.py

# To simplify the overall merged CSV file above, a corresponding species mapping.csv file for each column needs to be provided
# The filtering strategy is described in the Materials and Methods and Supplementary Materials.
python 02.filter_syntenic_matrix_by_species.py

# Extract syntenic CDS and protein sequence sets according to the simplified syntenic matrix
python 03.extract_acek_sequences.py

# === Step 2: Multiple sequence alignment with MAFFT ===
for f in pep-fa/*.fa; do
    # Align each FASTA file quietly using MAFFT --auto
    echo "mafft --auto --quiet $f > $f.fas"
done > 02.mafft-00.sh

# Run all alignment jobs in parallel (40 cores)
parallel -j 40 < 02.mafft-00.sh


# === Step 3: Using the pre-aligned protein sequence alignment to guide the generation of a corresponding CDS codon alignment ===
for f in pep-fa/*.fas                      
do                                                                                                                                                   
    base=$(basename $f .fas)
    echo "~/phyx/src/pxaa2cdn -a pep-fa/${base}.fas -n cds-fa/${base} -s -o cds-fa/${base}.aln"
done > 03.run_pxaa2cdn.sh

# Run all guiding jobs in parallel (40 cores)
parallel -j 40 < 03.run_pxaa2cdn.sh


# === Step 4: Trimming was performed using pxclsq.
for f in cds-fa/*.aln; do    
echo "~/phyx/src/pxclsq -s $f -o $f.cln -p 0.1"                                                                               
done >> 04.pxclsq.sh

# Run all pxclsq trimming jobs in parallel (40 cores)
parallel -j 40 < 04.pxclsq.sh

# === Step 5: Construct gene trees with IQ-TREE ===
for f in cds-fa/*.cln; do
    # Build maximum-likelihood tree with automatic model selection
    # 1000 ultrafast bootstrap replicates, quiet mode, redo if needed, auto threads
    echo "iqtree3 -s $f -pre ${f}.tre -bb 1000 -m MFP -quiet -redo -T AUTO"
done > 05.iqtree.sh

# Run all tree reconstruction jobs in parallel (40 cores)
parallel -j 10 < 05.iqtree.sh


# === Step 6: Get all gene trees
# Run combined.sh to merge all the gene trees
sh 06.combined.sh

# === Step 7: prune long branch nodes in every gene trees 
python 07.prune_long_branch_trees --trees combine.tsv --threshold 0.8 --correct 01.correct_trees.tsv --wrong 02.wrong_trees.tsv --pruned 03.pruned_trees.tsv --report 04.pruned_trees.report.tsv
# Obtain the new gene tree after removing all paralogous genes
cat 01.correct_trees.tsv 03.pruned_trees.tsv > combine.new.tsv

# === Step 8: Filter gene trees by order coverage ===
python 08.filter-orders -m mapping-5.csv -i combine.new.tsv -o combine-22order.tsv -n 22 -d tab

# === Step 9: Filter gene trees by species coverage ===
python 09.filter-species -i combine-22order.tsv -o combine-22order-73.tsv -n 73

# === Step 10: Extract cleaned gene sets after long-branch pruning and taxon-coverage filtering ===
# Extract retained tip names from each filtered gene tree.
# The output id directory contains one ID file for each gene tree.
# Each ID file links retained species/tip names to the corresponding gene IDs in the gene set.
python 10.extract_tree_leaves.py combine-22order-73.tsv

# Generate seqkit commands to extract the cleaned CDS sequence sets
for f in id/*
do
    base=$(basename "$f")
    echo "seqkit grep -f $f filtered_cds/${base}.fa > cds-fa/$base.fa"
done > seqkit_seq-cds.sh

# Extract the cleaned CDS sequence sets
parallel -j 10 < seqkit_seq-cds.sh

# Generate seqkit commands to extract the cleaned protein sequence sets
for f in id/*
do
    base=$(basename "$f")
    echo "seqkit grep -f $f filtered_pep/${base}.fa > pep-fa/$base.fa"
done > seqkit_seq-pep.sh

# Extract the cleaned protein sequence sets
parallel -j 10 < seqkit_seq-pep.sh

# Reconstruct final gene trees using the cleaned gene sets.
# Steps 1, 2, 3, 4 and 5 are repeated to generate the final rebuilt gene tree set:
# combine-rebuilt.tsv

# === Step 11: Infer the collinearity-based species tree using ASTRAL 5.6.9 ===
cut -f 2 combine-rebuilt.tsv > combine-rebuilt.nwk
java -jar ~/Constrained-search-master/astral.5.6.9.jar -i combine-rebuilt.nwk -o combine-rebuilt.nwk.tre -t 2

# === Step 12: Infer the final merged phylogenetic tree using ASTRAL 5.6.9 with a constraint tree ===
# The constraint tree was constructed primarily from shared chromosomal fusion events and shared whole-genome polyploidization events, while collinearity-based gene trees were used as supplementary evidence.
java -jar ~/Constrained-search-master/astral.5.6.9.jar -j species_constraint.tre -i combine-rebuilt.nwk -o merged-phylogenetic.tre -t 2