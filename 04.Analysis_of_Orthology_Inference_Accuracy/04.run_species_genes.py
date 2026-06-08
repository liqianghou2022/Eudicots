#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import csv
import os
from collections import defaultdict

# ================== Path configuration ==================
MATRIX_FILE = "ACEK-222.id.csv"
GROUP_FILE = "00.grouped_genes_final.txt"   # Provides group names and gene mappings
SPECIES_MAP_FILE = "us-1kp-sp.csv"          # Three columns: full species name, species code, group name
GENE_TO_G_FILE = "222-ACEK.blast.txt"
BLAST_DIR = "/data/01/user257/project/2025.houliqiang/03.ACEK-gene.blast-new/01.202-gene/00.all.blast"
BLAST_SUFFIX = "_blast_filtered.txt"
SUMMARY_OUT = "summary_correct_rate.csv"
DETAIL_OUT = "detailed_all_species.txt"
# ========================================================

def load_gene_to_group(group_file):
    """
    Load mapping from gene ID -> group name from the grouping file.
    Also return the set of all group names (used as valid prefixes).
    """
    gene_to_group = {}
    all_groups = set()
    current_group = None
    with open(group_file, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("组:"):
                parts = line.split()
                if len(parts) >= 2:
                    current_group = parts[1]
                    all_groups.add(current_group)
            elif current_group is not None and not line.startswith("="):
                gene_id = line.strip()
                if gene_id:
                    gene_to_group[gene_id] = current_group
    return gene_to_group, all_groups

def load_species_mapping(map_file):
    """
    Return:
    - name_to_code: {full species name: species code}
    - code_to_expected_group: {species code: expected group name} (consistent with grouping file)
    """
    name_to_code = {}
    code_to_expected_group = {}
    with open(map_file, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line or ',' not in line:
                continue
            parts = line.split(',')
            if len(parts) >= 3:
                name = parts[0].strip()
                code = parts[1].strip()
                group = parts[2].strip()
                name_to_code[name] = code
                code_to_expected_group[code] = group
    return name_to_code, code_to_expected_group

def load_gene_to_g(map_file):
    gene_to_g = {}
    with open(map_file, 'r', encoding='utf-8') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                gene_to_g[parts[0]] = parts[1]
    return gene_to_g

def parse_matrix_with_groups(matrix_file, gene_to_group, name_to_code):
    """
    Parse the matrix:
    - Get species codes from column headers.
    - Map each gene ID to its group using the grouping file.
    Return: {query_gene: {group_name: (species_code, {set of gene IDs})}}
    """
    with open(matrix_file, 'r', encoding='utf-8') as f:
        lines = [line.rstrip('\n') for line in f if line.strip()]
    if not lines:
        return {}

    header = lines[0].split(',')
    col_code = []
    for idx, sp_name in enumerate(header):
        sp_name = sp_name.strip()
        if idx == 0:
            continue
        code = name_to_code.get(sp_name, sp_name)
        col_code.append(code)

    result = {}
    for line in lines[1:]:
        fields = line.split(',')
        if not fields:
            continue
        our_gene = fields[0].strip()
        if not our_gene:
            continue
        gene_data = {}
        for i, cell in enumerate(fields[1:]):
            if i >= len(col_code):
                break
            cell = cell.strip()
            if not cell:
                continue
            code = col_code[i]
            group = gene_to_group.get(cell, '__unmapped__')
            if group not in gene_data:
                gene_data[group] = (code, set())
            gene_data[group][1].add(cell)
        if gene_data:
            result[our_gene] = gene_data
    return result

def parse_blast_file_for_query(blast_path, code_to_expected_group, all_groups):
    """
    Parse BLAST file.
    Determine the group of each gene by its prefix (using all_groups as prefix set).
    Only keep entries whose group matches the expected group.

    Return: {group_name: {set of gene IDs}}
    """
    blast_data = defaultdict(set)
    if not os.path.exists(blast_path):
        return blast_data
    with open(blast_path, 'r', encoding='utf-8') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 2:
                continue
            species_code = parts[0]
            gene_id = parts[1]
            if species_code not in code_to_expected_group:
                continue
            expected_group = code_to_expected_group[species_code]

            # Use longest prefix match to determine which group the gene belongs to
            matched_group = None
            max_len = 0
            for g in all_groups:
                if gene_id.startswith(g) and len(g) > max_len:
                    matched_group = g
                    max_len = len(g)

            if matched_group == expected_group:
                blast_data[expected_group].add(gene_id)
            # If matched group is not the expected one, ignore it
    return dict(blast_data)

def main():
    print("Loading gene grouping file...")
    gene_to_group, all_groups = load_gene_to_group(GROUP_FILE)
    print(f"Loaded {len(gene_to_group)} genes and {len(all_groups)} groups.")

    print("Loading species mapping file...")
    name_to_code, code_to_expected_group = load_species_mapping(SPECIES_MAP_FILE)

    print("Loading gene correspondence file...")
    gene_to_g = load_gene_to_g(GENE_TO_G_FILE)

    print("Parsing collinearity matrix...")
    matrix_dict = parse_matrix_with_groups(MATRIX_FILE, gene_to_group, name_to_code)

    summary_rows = []
    detail_lines = []

    for our_gene, matrix_data in matrix_dict.items():
        if our_gene not in gene_to_g:
            continue
        g_number = gene_to_g[our_gene]
        blast_file = os.path.join(BLAST_DIR, f"{g_number}{BLAST_SUFFIX}")
        blast_data = parse_blast_file_for_query(blast_file, code_to_expected_group, all_groups)

        yes_count = 0
        no1_count = 0
        no2_count = 0
        species_results = []

        for group, (code, mat_genes) in matrix_data.items():
            if group == '__unmapped__':
                continue
            if group in blast_data:
                blast_genes = blast_data[group]
                if mat_genes & blast_genes:
                    status = "yes"
                    yes_count += 1
                else:
                    status = "no1"
                    no1_count += 1
                gene_items = []
                for gene in sorted(mat_genes):
                    tag = "yes" if gene in blast_genes else "no1"
                    gene_items.append(f"{gene}[{tag}]")
                species_results.append(f"{group}=>{'|'.join(gene_items)}")
            else:
                status = "no2"
                no2_count += 1
                gene_items = [f"{gene}[no2]" for gene in sorted(mat_genes)]
                species_results.append(f"{group}=>{'|'.join(gene_items)}")

        total_valid = yes_count + no1_count
        rate = (yes_count / total_valid * 100) if total_valid > 0 else 0.0

        summary_rows.append({
            "gene_us": our_gene,
            "gene_other": g_number,
            "total_groups": len(matrix_data),
            "valid_groups": total_valid,
            "yes": yes_count,
            "no1": no1_count,
            "no2": no2_count,
            "correct_rate": f"{rate:.2f}%"
        })

        detail_line = f"{our_gene}," + ",".join(species_results)
        detail_lines.append(detail_line)

    with open(SUMMARY_OUT, 'w', encoding='utf-8', newline='') as f:
        fieldnames = ["gene_us", "gene_other", "total_groups", "valid_groups",
                      "yes", "no1", "no2", "correct_rate"]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(summary_rows)
    print(f"Summary saved to {SUMMARY_OUT}")

    with open(DETAIL_OUT, 'w', encoding='utf-8') as f:
        f.write("\n".join(detail_lines))
    print(f"Detailed results saved to {DETAIL_OUT}")

if __name__ == "__main__":
    main()