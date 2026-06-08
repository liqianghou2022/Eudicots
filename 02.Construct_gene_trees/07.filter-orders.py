#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Filter a two-column gene tree file based on taxonomic group coverage.

Input file format:
    Column 1: gene_id
    Column 2: newick_tree

Mapping file format:
    No header
    Column 1: species ID
    Column 2: group name

Filtering logic:
    A tree is retained if it contains species from at least the specified number
    of groups. A group is considered "present" if at least one species ID from
    that group occurs among the leaf nodes of the tree.

Output file format:
    Column 1: gene_id
    Column 2: newick_tree
"""

import argparse
import csv
import sys
from io import StringIO

import pandas as pd
from Bio import Phylo


def load_groups_from_csv(csv_file):
    """
    Load species-to-group information from a CSV file without a header.

    Column 1: species ID
    Column 2: group name

    Returns:
        {
            "group1": ["2", "3", "4"],
            "group2": ["10", "11"]
        }
    """
    df = pd.read_csv(csv_file, header=None, names=["id", "group"], dtype=str)

    # Strip leading and trailing whitespace to avoid matching failures
    df["id"] = df["id"].astype(str).str.strip()
    df["group"] = df["group"].astype(str).str.strip()

    # Remove rows with missing or empty values
    df = df.dropna(subset=["id", "group"])
    df = df[(df["id"] != "") & (df["group"] != "")]

    groups = df.groupby("group")["id"].apply(list).to_dict()
    return groups


def parse_tree_from_string(tree_str):
    """
    Parse a Newick tree from a string.
    """
    handle = StringIO(tree_str.strip())
    tree = Phylo.read(handle, "newick")
    return tree


def check_tree_for_groups(tree, groups, required_groups_count=10):
    """
    Check whether a tree contains at least the specified number of groups.

    The logic is consistent with the original script:
    a group is counted as covered if any species from that group occurs in the tree.
    """
    terminal_names = {str(term.name).strip() for term in tree.get_terminals() if term.name is not None}
    groups_covered = 0

    for group_nodes in groups.values():
        if any(str(node).strip() in terminal_names for node in group_nodes):
            groups_covered += 1
        if groups_covered >= required_groups_count:
            return True

    return False


def detect_delimiter(file_path, user_delimiter=None):
    """
    Detect the input delimiter automatically or use the delimiter specified by the user.

    Supported delimiters:
        tab / comma / semicolon / pipe / space
    """
    if user_delimiter:
        mapping = {
            "tab": "\t",
            "comma": ",",
            "semicolon": ";",
            "pipe": "|",
            "space": " ",
        }
        if user_delimiter not in mapping:
            raise ValueError(f"Unsupported delimiter type: {user_delimiter}")
        return mapping[user_delimiter]

    # Automatically detect the delimiter
    with open(file_path, "r", encoding="utf-8") as f:
        sample = f.read(4096)

    try:
        dialect = csv.Sniffer().sniff(sample, delimiters="\t,;| ")
        return dialect.delimiter
    except csv.Error:
        # Use tab as the default delimiter
        return "\t"


def filter_trees_by_groups(
    input_file,
    groups,
    output_file,
    required_groups_count=10,
    delimiter=None,
    skip_bad_trees=False
):
    """
    Filter a two-column gene tree file based on group coverage.

    Parameters:
        input_file: two-column input file; column 1 is gene_id and column 2 is tree
        groups: dictionary of group -> [species IDs]
        output_file: output file path
        required_groups_count: minimum number of groups required for a tree to be retained
        delimiter: input and output delimiter
        skip_bad_trees: if True, malformed trees are skipped; otherwise, the script exits with an error

    Returns:
        total_count: total number of input trees
        valid_count: number of retained trees
        bad_count: number of trees that failed to parse or were skipped
    """
    total_count = 0
    valid_count = 0
    bad_count = 0

    with open(input_file, "r", encoding="utf-8") as fin, \
         open(output_file, "w", encoding="utf-8", newline="") as fout:

        reader = csv.reader(fin, delimiter=delimiter)
        writer = csv.writer(fout, delimiter=delimiter)

        for line_num, row in enumerate(reader, start=1):
            if not row:
                continue

            total_count += 1

            if len(row) < 2:
                msg = f"[Line {line_num}] Fewer than two columns; unable to read gene_id and tree. Content: {row}"
                if skip_bad_trees:
                    print("Warning:", msg, file=sys.stderr)
                    bad_count += 1
                    continue
                raise ValueError(msg)

            gene_id = row[0].strip()
            tree_str = row[1].strip()

            if gene_id == "" or tree_str == "":
                msg = f"[Line {line_num}] gene_id or tree is empty."
                if skip_bad_trees:
                    print("Warning:", msg, file=sys.stderr)
                    bad_count += 1
                    continue
                raise ValueError(msg)

            try:
                tree = parse_tree_from_string(tree_str)
            except Exception as e:
                msg = f"[Line {line_num}] Failed to parse tree for gene_id={gene_id}. Error: {e}"
                if skip_bad_trees:
                    print("Warning:", msg, file=sys.stderr)
                    bad_count += 1
                    continue
                raise ValueError(msg)

            if check_tree_for_groups(tree, groups, required_groups_count):
                # Write back the original tree string to preserve the input format as much as possible
                writer.writerow([gene_id, tree_str])
                valid_count += 1

    return total_count, valid_count, bad_count


def main():
    parser = argparse.ArgumentParser(
        description="Filter a two-column gene tree file based on group coverage "
                    "(column 1: gene_id; column 2: newick_tree)."
    )

    parser.add_argument(
        "-m", "--mapping",
        required=True,
        help="Path to the mapping file; no header, column 1 is species ID and column 2 is group."
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Path to the input gene tree file; two columns: gene_id and tree."
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Path to the output file."
    )
    parser.add_argument(
        "-n", "--required-groups",
        type=int,
        required=True,
        help="Minimum number of groups required for a tree to be retained."
    )
    parser.add_argument(
        "-d", "--delimiter",
        choices=["tab", "comma", "semicolon", "pipe", "space"],
        default=None,
        help="Delimiter of the input/output file. If not specified, it will be detected automatically."
    )
    parser.add_argument(
        "--skip-bad-trees",
        action="store_true",
        help="Skip malformed or unparseable trees and continue processing."
    )

    args = parser.parse_args()

    groups = load_groups_from_csv(args.mapping)

    if len(groups) == 0:
        raise ValueError("No valid groups were read from the mapping file.")

    delimiter = detect_delimiter(args.input, args.delimiter)

    total_count, valid_count, bad_count = filter_trees_by_groups(
        input_file=args.input,
        groups=groups,
        output_file=args.output,
        required_groups_count=args.required_groups,
        delimiter=delimiter,
        skip_bad_trees=args.skip_bad_trees
    )

    print(f"Number of groups in the mapping file: {len(groups)}")
    print(f"Total input trees: {total_count}")
    print(f"Number of retained trees: {valid_count}")
    print(f"Number of failed/skipped trees: {bad_count}")
    print(f"Output file: {args.output}")


if __name__ == "__main__":
    main()