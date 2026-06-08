"""
Script: Species-level filtering of syntenic gene matrix for recent WGD lineages

Description
-----------
Based on the previously constructed syntenic gene matrix, we further selected
syntenic genes from species that experienced recent whole-genome duplication (WGD)
events. This filtering step reduces matrix redundancy and substantially decreases
computational demands in downstream analyses.

For species with multiple retained syntenic gene copies resulting from recent WGDs,
one representative gene is selected per ACEK gene. When multiple valid candidates
are available, one gene is randomly chosen. Species without recent WGDs typically
retain a single syntenic copy and are kept unchanged.

Detailed biological rationale and validation of this filtering strategy are provided
in Supplementary Figure 16 and are not repeated here.

Input
-----
1. Syntenic gene matrix (CSV format)
   - File: merge.all.csv
   - No header
   - Rows correspond to ACEK genes
   - Columns correspond to syntenic genes from multiple species
   - All values are read as strings

2. Column-to-species mapping file (CSV format)
   - File: mapping.csv
   - Two columns (no header):
       * Column 1: column index in the syntenic matrix
       * Column 2: species or species-group identifier
   - Columns sharing the same identifier belong to the same species

Output
------
- Filtered syntenic gene matrix (CSV format)
  - File: Species-one-gene.csv
  - One column per species
  - Each column contains a single representative syntenic gene per ACEK gene
  - Missing values indicate gene loss or absence of detectable synteny
"""
import pandas as pd
import numpy as np

# Read the syntenic gene matrix (all values as strings)
data = pd.read_csv('merge.all.csv', header=None, dtype=str)

# Assign sequential column indices starting from 1
data.columns = range(1, data.shape[1] + 1)

# Read the column-to-species mapping file
mapping = pd.read_csv('mapping.csv', header=None, names=['col', 'type'])

# Store one representative column per species
cols_dict = {}

# Process columns by species (or species group)
for t, group in mapping.groupby('type'):
    cols = group['col'].tolist()
    sub = data[cols]

    # Single-copy species: keep directly
    if len(cols) == 1:
        chosen = sub.iloc[:, 0]
    else:
        # Multi-copy species (recent WGD): select one representative per row
        def choose_value(row):
            values = [v for v in row if pd.notna(v) and v != '']
            if len(values) == 0:
                return np.nan
            elif len(values) == 1:
                return values[0]
            else:
                return np.random.choice(values)

        chosen = sub.apply(choose_value, axis=1)

    cols_dict[t] = chosen

# Merge all species-level columns
result = pd.concat(cols_dict, axis=1)

# Write the filtered syntenic gene matrix
result.to_csv('Species-one-gene.csv', index=False)
