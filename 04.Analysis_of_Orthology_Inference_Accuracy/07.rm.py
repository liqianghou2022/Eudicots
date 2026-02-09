import csv, sys

infile  = "2.csv"   # Replace with your input file
outfile = "acek_mapped.csv"

with open(infile, "r", encoding="utf-8") as f:
    lines = [line.rstrip("\n") for line in f if line.strip() != ""]

if len(lines) < 3:
    sys.exit("Input file must contain at least three lines: line 1 = IDs, line 2 = statistics, line 3 = filtered IDs")

# Split fields and trim whitespace
def split_clean(s):
    # Use csv module for robustness (also handles quoted fields)
    return [cell.strip() for cell in next(csv.reader([s]))]

ids_all   = split_clean(lines[0])
stats_all = split_clean(lines[1])
ids_keep  = split_clean(lines[2])

if len(ids_all) != len(stats_all):
    sys.exit(f"Number of IDs in line 1 ({len(ids_all)}) does not match number of statistics in line 2 ({len(stats_all)})")

# Build ID -> statistic mapping (keep the first occurrence if duplicated)
id2stat = {}
for i, gid in enumerate(ids_all):
    if gid not in id2stat:
        id2stat[gid] = stats_all[i] if i < len(stats_all) else ""

# Map IDs in line 3 to statistics; fill missing values with empty or 0
mapped_stats = []
missing = []
for gid in ids_keep:
    if gid in id2stat:
        mapped_stats.append(id2stat[gid])
    else:
        mapped_stats.append("")   # Change to "0" if preferred
        missing.append(gid)

# Write two-line CSV: first row = filtered IDs; second row = corresponding statistics
with open(outfile, "w", newline="", encoding="utf-8") as out:
    w = csv.writer(out)
    w.writerow(ids_keep)
    w.writerow(mapped_stats)

print(f"Wrote: {outfile}")
if missing:
    print("Warning: the following IDs were not found in line 1; values were left empty:")
    for m in missing:
        print("  -", m)