#!/usr/bin/env python3
"""
Script: 02.merge_csv_union.py
Author: Liqiang Hou
Date: 2025.08

Description:
    Merge multiple CSV files using the first column as key (outer join).
    Missing values are left empty.

Usage:
    python 02.merge_csv_union.py -i "/path/to/csv_dir/*.csv" -o /path/to/output.csv
    
    Or edit the csv_files_pattern and output_file variables below and run:
    python 02.merge_csv_union.py

Input:
    - Multiple CSV files without headers
    - First column is used as the merge key

Output:
    - Single merged CSV file
    - All rows from all files are kept (union)

Example:
    file1.csv: gene1,10,20
               gene2,30,40
    file2.csv: gene1,15,25
               gene3,35,45
    Result:    gene1,10,20,15,25
               gene2,30,40,,
               gene3,,,35,45

Requirements:
    - pandas (pip install pandas)
"""

import pandas as pd
import glob
import sys

# Parse command line arguments if provided
if len(sys.argv) > 1 and '-i' in sys.argv:
    idx = sys.argv.index('-i')
    csv_files_pattern = sys.argv[idx + 1]
else:
    # Edit this path if not using command line
    csv_files_pattern = '/path/to/csv_dir/*.csv'

if len(sys.argv) > 1 and '-o' in sys.argv:
    idx = sys.argv.index('-o')
    output_file = sys.argv[idx + 1]
else:
    # Edit this path if not using command line
    output_file = '/path/to/output/merged.csv'

def merge_csv_by_first_column(csv_files, output_file):
    """Merge CSV files using first column as key"""
    if not csv_files:
        print("Error: No CSV files found!")
        return
    
    print(f"Merging {len(csv_files)} files...")
    
    # Initialize empty DataFrame
    merged_df = None
    
    for file in csv_files:
        # Read CSV using first column as index
        df = pd.read_csv(file, header=None, index_col=0)
        
        if merged_df is None:
            merged_df = df
        else:
            # Merge with outer join (keeps all rows)
            merged_df = pd.merge(merged_df, df, left_index=True, right_index=True, how='outer')
    
    # Reset index to make first column a regular column
    merged_df.reset_index(drop=False, inplace=True)
    
    # Save to CSV without headers
    merged_df.to_csv(output_file, index=False, header=False)
    print(f"Done! Output saved to: {output_file}")

# Get list of CSV files
csv_files = glob.glob(csv_files_pattern)

# Run merge
merge_csv_by_first_column(csv_files, output_file)