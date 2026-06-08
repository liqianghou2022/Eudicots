#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import sys
from ete3 import Tree


def safe_filename(name):
    """Convert a gene name into a safe filename."""
    name = name.strip()
    name = re.sub(r'[\\/:"*?<>|]+', "_", name)
    return name


def main():
    if len(sys.argv) != 2:
        print(f"Usage: python {sys.argv[0]} input.tsv")
        sys.exit(1)

    input_file = sys.argv[1]
    outdir = "id"
    os.makedirs(outdir, exist_ok=True)

    total = 0
    success = 0
    failed = 0

    with open(input_file, "r", encoding="utf-8") as f:
        for line_num, line in enumerate(f, start=1):
            line = line.rstrip("\n")
            if not line.strip():
                continue

            parts = line.split("\t")
            if len(parts) < 2:
                print(f"[Warning] Line {line_num}: less than 2 columns, skipped", file=sys.stderr)
                failed += 1
                continue

            gene_id = parts[0].strip()
            tree_str = parts[1].strip()

            if not gene_id or not tree_str:
                print(f"[Warning] Line {line_num}: empty gene_id or tree, skipped", file=sys.stderr)
                failed += 1
                continue

            try:
                # format=1 is more suitable for trees containing internal node support values and branch lengths.
                try:
                    t = Tree(tree_str, format=1)
                except Exception:
                    t = Tree(tree_str, format=0)

                leaf_names = [leaf.name.strip() for leaf in t.iter_leaves()]

                outname = safe_filename(gene_id)
                outfile = os.path.join(outdir, outname)

                with open(outfile, "w", encoding="utf-8") as out:
                    for leaf in leaf_names:
                        out.write(leaf + "\n")

                success += 1

            except Exception as e:
                print(f"[Warning] Line {line_num}, gene {gene_id}: parse failed: {e}", file=sys.stderr)
                failed += 1

            total += 1

    print("Done.")
    print(f"Total processed lines: {total}")
    print(f"Successful: {success}")
    print(f"Failed: {failed}")
    print(f"Output directory: {outdir}")


if __name__ == "__main__":
    main()