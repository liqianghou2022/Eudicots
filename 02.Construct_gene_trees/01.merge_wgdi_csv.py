"""
Script: merge WGDI syntenic gene CSV files by ACEK gene IDs

Description
-----------
After obtaining syntenic gene sets between the ancestral core eudicot karyotype (ACEK)
and each extant species using WGDI (https://github.com/SunPengChuan/wgdi), one CSV file
is generated per species. Therefore, the total number of CSV files corresponds to the
number of species analyzed.

For each CSV file:
- The first column contains ACEK gene IDs.
- The remaining columns contain syntenic genes from the corresponding species.
- Each column represents one syntenic gene copy identified by collinearity analysis.
- If multiple columns are present, this indicates that the species experienced
  additional whole-genome duplication (WGD) events after the γ (gamma) hexaploidization.

This script merges all species-specific CSV files into a single matrix using the
ACEK gene IDs (first column) as the common key. The merged matrix integrates syntenic
gene information across all species, while missing values indicate gene loss or the
absence of detectable collinearity relationships.

Input
-----
- Multiple CSV files (*.csv)
- No header row
- First column: ACEK gene IDs
- Remaining columns: syntenic genes from one species

Output
------
- One merged CSV file
- Rows correspond to ACEK gene IDs
- Columns correspond to syntenic gene copies across all species
- Missing values indicate lineage-specific gene loss or absence of synteny
"""

import pandas as pd
import glob

def merge_csv_by_first_column(csv_files, output_file):
    # Initialize an empty DataFrame
    merged_df = None

    for file in csv_files:
        # Read the current CSV file
        df = pd.read_csv(file, header=None, index_col=0)

        # If merged_df is empty, assign directly; otherwise, merge by the first column (index)
        if merged_df is None:
            merged_df = df
        else:
            # Merge using the first column (ACEK gene IDs) as the index
            merged_df = pd.merge(
                merged_df,
                df,
                left_index=True,
                right_index=True,
                how='outer'
            )

    # Reset index so that the ACEK gene IDs become the first column again
    merged_df.reset_index(drop=False, inplace=True)

    # Write the merged DataFrame to a new CSV file
    merged_df.to_csv(output_file, index=False, header=False)
    print(f"Merging completed. Output file written to: {output_file}")

# Path to the directory containing species-specific CSV files
csv_files = glob.glob(
    './*.csv'
)

# Path and filename of the merged output CSV file
output_file = (
    './merge.all.csv'
)

merge_csv_by_first_column(csv_files, output_file)
