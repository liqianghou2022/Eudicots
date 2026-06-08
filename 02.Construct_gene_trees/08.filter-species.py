#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from ete3 import Tree
import argparse
import csv
import sys


def get_delimiter(delim_name):
    mapping = {
        "tab": "\t",
        "comma": ",",
        "space": " ",
        "semicolon": ";",
        "pipe": "|"
    }
    return mapping[delim_name]


def count_filtered_leaves(tree, excluded_nodes):
    """
    Count the number of leaf nodes included in the filtering criterion.

    Leaf names listed in excluded_nodes are ignored and not counted.
    """
    leaf_names = tree.get_leaf_names()
    kept = [x for x in leaf_names if x not in excluded_nodes]
    return len(kept)


def main():
    parser = argparse.ArgumentParser(
        description="Filter a two-column gene tree file by the number of leaf nodes "
                    "(column 1: gene_id; column 2: gene_tree)."
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Path to the input file; two columns: gene_id and tree."
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Path to the output file containing gene trees that pass the filter."
    )
    parser.add_argument(
        "-n", "--minimum-species-threshold",
        type=int,
        required=True,
        help="Minimum number of retained leaf nodes required after excluding specified nodes."
    )
    parser.add_argument(
        "-d", "--delimiter",
        choices=["tab", "comma", "space", "semicolon", "pipe"],
        default="tab",
        help="Delimiter used in the input/output files. Default: tab."
    )
    parser.add_argument(
        "--exclude",
        default="2,3,4",
        help="Comma-separated leaf names to exclude from the filtering count. Default: 2,3,4."
    )
    parser.add_argument(
        "--skip-bad-trees",
        action="store_true",
        help="Skip trees that fail to parse instead of terminating with an error."
    )

    args = parser.parse_args()

    delimiter = get_delimiter(args.delimiter)
    excluded_nodes = {x.strip() for x in args.exclude.split(",") if x.strip()}

    total_trees = 0
    kept_trees = 0
    bad_trees = 0

    with open(args.input, "r", encoding="utf-8") as infile,          open(args.output, "w", encoding="utf-8", newline="") as outfile:

        reader = csv.reader(infile, delimiter=delimiter)
        writer = csv.writer(outfile, delimiter=delimiter)

        for line_num, row in enumerate(reader, start=1):
            if not row:
                continue

            total_trees += 1

            if len(row) < 2:
                msg = f"[Line {line_num}] Fewer than two columns; unable to read gene_id and tree: {row}"
                if args.skip_bad_trees:
                    print("Warning:", msg, file=sys.stderr)
                    bad_trees += 1
                    continue
                else:
                    raise ValueError(msg)

            gene_id = row[0].strip()
            tree_str = row[1].strip()

            if gene_id == "" or tree_str == "":
                msg = f"[Line {line_num}] gene_id or tree is empty."
                if args.skip_bad_trees:
                    print("Warning:", msg, file=sys.stderr)
                    bad_trees += 1
                    continue
                else:
                    raise ValueError(msg)

            try:
                tree = Tree(tree_str, format=1)
                num_species = count_filtered_leaves(tree, excluded_nodes)

                if num_species >= args.minimum_species_threshold:
                    writer.writerow([gene_id, tree_str])
                    kept_trees += 1

            except Exception as e:
                msg = f"[Line {line_num}] Failed to parse tree for gene_id={gene_id}. Error: {e}"
                if args.skip_bad_trees:
                    print("Warning:", msg, file=sys.stderr)
                    bad_trees += 1
                    continue
                else:
                    raise ValueError(msg)

    print(f"Total trees: {total_trees}")
    print(f"Retained trees: {kept_trees}")
    print(f"Failed tree parses: {bad_trees}")
    print(f"Excluded nodes: {','.join(sorted(excluded_nodes))}")
    print(f"Filtering threshold: {args.minimum_species_threshold}")
    print(f"Output file: {args.output}")


if __name__ == "__main__":
    main()
