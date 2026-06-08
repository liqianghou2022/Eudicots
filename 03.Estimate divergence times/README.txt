## 01.multi_fasta_concatenate.py: Concatenation of Cleaned CDS Alignments

The script `01.multi_fasta_concatenate.py` was used to concatenate gene-wise CDS alignments from the `cds-fa-rebuilt-aln-cln` directory. The resulting concatenated alignment matrix was used for downstream branch-length estimation.

Usage:

```bash
python 01.multi_fasta_concatenate.py
```

## 02.random_submatrix_sampler.py: Random Subsampling of the Concatenated Alignment Matrix

The script `02.random_submatrix_sampler.py` was used to randomly sample submatrices from the concatenated FASTA alignment. A total of 1,000 random subsampled matrices were generated and subsequently used to estimate branch lengths on the fixed topology of `merged-phylogenetic.tre`.

Usage:

```bash
python 02.random_submatrix_sampler.py
```

## Branch-Length Estimation for 1,000 Randomly Sampled Matrices

IQ-TREE was used to estimate branch lengths for each randomly sampled FASTA matrix while fixing the topology to `merged-phylogenetic.tre`. The GTR+G substitution model was applied to DNA alignments.

```bash
# Generate IQ-TREE commands for all sampled FASTA files to estimate branch lengths
for f in random-1000/*.fasta; do
    echo "/groups/lzu_public/home/houlq2025/software/iqtree-3.0.1-Linux-intel/bin/iqtree3 -s $f -st DNA -m GTR+G -te merged-phylogenetic.tre -T AUTO"
done > run.iqtree.sh

sh run.iqtree.sh
```

## Divergence-Time Estimation Using TreePL

Divergence times were estimated using `03.TreePL_time.sh`. Cross-validation was first performed to identify the optimal smoothing parameter.

The cross-validation settings were:

```text
cv = true
cvstart = 0.0001
cvstop = 1000
cvmultstep = 10
randomcv = true
cviter = 5
```

After the optimal smoothing value was determined based on the cross-validation results, the final divergence-time estimation was performed with cross-validation disabled:

```text
smooth = best_smoothing_value
cv = false
```