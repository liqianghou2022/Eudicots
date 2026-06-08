01.filter_alignment_rate.sh	Filters all .txt files, keeping lines where the third column ≥ 90, and saves them as _filtered.txt.
02.group_genes.py	                This script extracts gene IDs from a CSV file, groups them by prefix, allows manual editing of group mappings, and outputs the final grouped gene list.
03.add_underscore	                Replaces spaces with underscores in us-1kp-sp.csv in place.
04.run_species_genes	Processes gene collinearity and BLAST results to calculate group-level correctness and generate summary and detailed output files.