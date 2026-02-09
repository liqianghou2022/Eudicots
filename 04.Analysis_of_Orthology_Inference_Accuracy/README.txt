01.run_species-genes.py	This script searches BLAST result files for gene IDs linked to each standard gene and outputs per-gene summary CSV files.
02.count_num		It summarizes multiple CSV files by counting how many species have a “yes” and reports duplicates.
03.transposed.py		It transposes a CSV file (rows ↔ columns) and saves the result as a new CSV.
04.rm_unnamed.sh		It deletes Unnamed: columns from the CSV, skipping the header, and overwrites the file.
05.1_0.py			This script converts all non-missing values in the CSV to 1 and missing values to 0, then saves it as a new file.
06.rm_dup.sh		Merges rows by the first column taking the max per numeric field, then lists columns that sum to zero.
07.rm.py			Maps filtered IDs to their statistics from the original file and outputs a two-row CSV, leaving missing values blank.