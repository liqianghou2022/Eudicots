#!/usr/bin/env python3
import csv
import subprocess
import shlex
from pathlib import Path
import argparse
import sys

# Configuration (can be overridden by command-line arguments)
MAP_TXT = "222-ACEK.blast.txt"            # Text file: each line like std_gene<TAB>gXXXX (optional validation)
IDS_CSV = "37-ACEK-202.id.csv"            # First column: standard gene; remaining columns: corresponding gene IDs
SPECIES2GENES_CSV = "species-gene.csv"    # Each row: species_code, gene1, gene2, ...
BLAST_DIR = "00.all.blast"
BLAST_GLOB = "*_blast_filtered.txt"       # Glob pattern for BLAST input files
OUTPUT_DIR = "02.result"
STD_GENE_PREFIX = "acek_cja_vvi"          # Only process standard genes with this prefix (empty = process all)

def parse_args():
    ap = argparse.ArgumentParser(description="Search BLAST results and write per-standard-gene CSV files.")
    ap.add_argument("--blast-dir", default=BLAST_DIR, help="Directory containing BLAST files")
    ap.add_argument("--glob", default=BLAST_GLOB, help="Glob pattern of BLAST files in --blast-dir")
    ap.add_argument("--output-dir", default=OUTPUT_DIR, help="Directory for output CSV files")
    ap.add_argument("--map-txt", default=MAP_TXT, help="TXT file for optional standard-gene validation")
    ap.add_argument("--ids-csv", default=IDS_CSV)
    ap.add_argument("--species2genes-csv", default=SPECIES2GENES_CSV)
    ap.add_argument("--std-prefix", default=STD_GENE_PREFIX,
                    help="Only process standard genes starting with this prefix (empty for all)")
    return ap.parse_args()

def grep_lines_fixed(pattern, file):
    cmd = f"grep -F {shlex.quote(pattern)} {shlex.quote(str(file))}"
    try:
        out = subprocess.check_output(cmd, shell=True, stderr=subprocess.DEVNULL)
        return out.decode("utf-8", errors="ignore").splitlines()
    except subprocess.CalledProcessError:
        return []

# 0) Load species → genes mapping and build gene → species map (first occurrence wins)
def load_gene2species(csv_path):
    gene2species = {}
    with open(csv_path, newline="", encoding="utf-8") as fin:
        rdr = csv.reader(fin)
        for row in rdr:
            if not row:
                continue
            sp = row[0].strip()
            for gid in row[1:]:
                gid = gid.strip()
                if gid and gid not in gene2species:
                    gene2species[gid] = sp
    return gene2species

# 1) Load standard genes and their associated gene IDs from CSV
def load_std_gene_rows(ids_csv, prefix=None):
    result = {}  # std_gene -> [gene_ids...]
    with open(ids_csv, newline="", encoding="utf-8") as fin:
        rdr = csv.reader(fin)
        for row in rdr:
            if not row:
                continue
            std_gene = row[0].strip()
            if not std_gene:
                continue
            if prefix and not std_gene.startswith(prefix):
                continue
            gene_ids = [x.strip() for x in row[1:] if x and x.strip()]
            result[std_gene] = gene_ids
    return result

# 2) Optional: check whether standard genes are present in map.txt (no interruption if missing)
def read_std_in_map_txt(map_txt):
    present = set()
    try:
        with open(map_txt, encoding="utf-8") as f:
            for line in f:
                line = line.rstrip("\n")
                if not line or line.startswith("#"):
                    continue
                # Allow a single column; if TAB-separated, use the first field
                parts = line.split("\t")
                if parts:
                    present.add(parts[0])
    except FileNotFoundError:
        pass
    return present

def main():
    args = parse_args()

    blast_dir = Path(args.blast_dir)
    if not blast_dir.is_dir():
        raise SystemExit(f"ERROR: BLAST_DIR not found: {blast_dir}")

    # List all BLAST files to be searched (using glob pattern)
    blast_files = sorted(p for p in blast_dir.glob(args.glob) if p.is_file())
    if not blast_files:
        raise SystemExit(f"ERROR: no {args.glob} found in {blast_dir}")

    # Create output directory if needed
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    gene2species = load_gene2species(args.species2genes_csv)
    std_rows = load_std_gene_rows(args.ids_csv, prefix=(args.std_prefix or None))
    std_in_map = read_std_in_map_txt(args.map_txt)

    print(f"Total standard genes to process: {len(std_rows)}")
    print(f"BLAST files to search: {len(blast_files)}")
    print(f"Output directory: {out_dir.resolve()}")

    for std_gene, row_gene_ids in std_rows.items():
        out_file = out_dir / f"{std_gene}.final.csv"
        with open(out_file, "w", newline="", encoding="utf-8") as fout:
            w = csv.writer(fout)
            w.writerow(["gene_id", "species_code", "found_in_blast", "total_hits", "files"])

            for gid in row_gene_ids:
                sp = gene2species.get(gid, "")
                total_hits = 0
                hit_files = []
                if sp:
                    for bf in blast_files:
                        lines = grep_lines_fixed(sp, bf)
                        if not lines:
                            continue
                        hits = [ln for ln in lines if gid in ln]
                        if hits:
                            total_hits += len(hits)
                            hit_files.append(bf.name)
                w.writerow([
                    gid,
                    sp,
                    "yes" if total_hits > 0 else "no",
                    total_hits,
                    ";".join(hit_files)
                ])

        # Optional warning if standard gene not found in map.txt
        if std_in_map and std_gene not in std_in_map:
            print(f"Warning: {std_gene} not present in map.txt (validation skipped).")

        print(f"Wrote {out_file}")

    print("All done.")

if __name__ == "__main__":
    main()
