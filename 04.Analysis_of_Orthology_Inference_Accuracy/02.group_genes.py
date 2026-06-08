import csv
import re
from collections import defaultdict

# ========== Configuration ==========
CSV_FILE = "ACEK-222.id.csv"          # Input file
PREFIX_LIST_FILE = "00.prefix_list.csv"  # Step 1 output (for manual editing)
MAPPING_FILE = "00.group_mapping.csv"    # Mapping file after editing
OUTPUT_FILE = "00.grouped_genes_final.txt"  # Final grouped result
# ================================

PATTERN = re.compile(r'^(.*?)([gG]\d+)$')

def extract_gene_ids(file_path):
    """Read CSV file, skip header, and collect all non-empty gene IDs"""
    genes = []
    with open(file_path, 'r', encoding='utf-8') as f:
        reader = csv.reader(f)
        first_row = True
        for row in reader:
            if first_row:
                first_row = False
                # If the first cell does not end with g/G + digits, treat as header and skip
                if row and not re.search(r'[gG]\d+$', row[0].strip()):
                    continue
            for cell in row:
                cell = cell.strip()
                if cell:
                    genes.append(cell)
    return genes

def generate_prefix_list(genes):
    """Extract all unique prefixes and provide suggested group names"""
    prefixes = defaultdict(set)
    unmatched = []
    for gid in genes:
        m = PATTERN.match(gid)
        if m:
            prefix = m.group(1)
            prefixes[prefix].add(gid)
        else:
            unmatched.append(gid)
    return prefixes, unmatched

def suggest_group_name(prefix):
    """Automatically suggest a group name by removing trailing numbers and optional lowercase letters"""
    return re.sub(r'[\d_]+[a-z]?$', '', prefix)

def write_prefix_list(prefixes, unmatched, output_file):
    """Generate a CSV file with three columns: prefix, suggested group name, and example gene"""
    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['prefix', 'suggested_group', 'example_gene'])
        for prefix in sorted(prefixes.keys()):
            suggestion = suggest_group_name(prefix)
            example = list(prefixes[prefix])[0]
            writer.writerow([prefix, suggestion, example])
        for gid in unmatched:
            writer.writerow(['', '', gid])

def read_mapping(mapping_file):
    """Read the user-edited mapping file and return a dictionary {prefix: final_group_name}"""
    mapping = {}
    with open(mapping_file, 'r', encoding='utf-8') as f:
        reader = csv.reader(f)
        next(reader)  # Skip header
        for row in reader:
            if len(row) >= 2 and row[0].strip() and row[1].strip():
                prefix = row[0].strip()
                group = row[1].strip()
                mapping[prefix] = group
    return mapping

def group_genes_with_mapping(genes, mapping):
    """Group gene IDs based on the mapping table"""
    groups = defaultdict(list)
    unmatched = []
    for gid in genes:
        m = PATTERN.match(gid)
        if m:
            prefix = m.group(1)
            group_name = mapping.get(prefix, prefix)  # Use prefix itself if not mapped
            groups[group_name].append(gid)
        else:
            unmatched.append(gid)
    return groups, unmatched

def print_groups(groups, unmatched, output_file=None):
    """Print grouped results and optionally save to file"""
    lines = ["===== Final Grouping Result =====\n"]
    for group in sorted(groups.keys()):
        lines.append(f"Group: {group}  (Total {len(groups[group])} genes)")
        for g in sorted(groups[group]):
            lines.append(f"  {g}")
        lines.append("")
    if unmatched:
        lines.append("===== Unmatched Gene IDs =====\n")
        for g in sorted(unmatched):
            lines.append(f"  {g}")
    output_text = "\n".join(lines)
    print(output_text)
    if output_file:
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(output_text)
        print(f"\nResult saved to: {output_file}")

# ========== Main Workflow ==========
if __name__ == "__main__":
    print("Reading CSV file...")
    genes = extract_gene_ids(CSV_FILE)
    print(f"Total {len(genes)} gene IDs loaded\n")

    # Step 1: Generate prefix list
    print("Step 1: Generating prefix list for manual editing...")
    prefixes, unmatched = generate_prefix_list(genes)
    write_prefix_list(prefixes, unmatched, PREFIX_LIST_FILE)
    print(f"File '{PREFIX_LIST_FILE}' has been generated. Please open and edit the 'suggested_group' column.")
    print("Assign the same group name to prefixes you want to merge (e.g., set both FvirAHapA and FvirBHapB to Fvir).")
    print(f"After editing, save it as '{MAPPING_FILE}' in the same directory.")
    print("\nPress Enter to continue after finishing...")
    input()

    # Step 2: Read mapping and group genes
    try:
        mapping = read_mapping(MAPPING_FILE)
    except FileNotFoundError:
        print(f"Error: mapping file '{MAPPING_FILE}' not found. Please make sure it is created as instructed.")
        exit(1)

    print("Grouping genes based on mapping file...")
    groups, unmatched = group_genes_with_mapping(genes, mapping)
    print_groups(groups, unmatched, OUTPUT_FILE)