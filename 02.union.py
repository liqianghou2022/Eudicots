#!/usr/bin/env python3
# Author: Liqiang Hou 2025.08
# usage: python 02.union.py -i "/path/to/csv_dir/*.csv" -o /path/to/output.csv  (or edit the variables below)

import pandas as pd
import glob

def merge_csv_by_first_column(csv_files, output_file):
    # Initialize an empty DataFrame
    merged_df = None

    for file in csv_files:
        # Read the current CSV file
        df = pd.read_csv(file, header=None, index_col=0)
        # Check if merged_df is empty; if so, assign directly; otherwise, merge
        if merged_df is None:
            merged_df = df
        else:
            # Merge by the first column (used as index here)
            merged_df = pd.merge(merged_df, df, left_index=True, right_index=True, how='outer')

    # Reset the index because the index (first column) is part of the data
    merged_df.reset_index(drop=False, inplace=True)

    # Write the merged DataFrame to a new CSV file
    merged_df.to_csv(output_file, index=False, header=False)
    print(f"Merging completed. Output written to: {output_file}")

# Change to the folder path pattern containing your CSV files (glob pattern supported)
csv_files = glob.glob('/path/to/csv_dir/*.csv')
# Change to your desired output file path and name
output_file = '/path/to/output/merged.csv'

merge_csv_by_first_column(csv_files, output_file)