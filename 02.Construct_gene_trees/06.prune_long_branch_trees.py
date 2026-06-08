#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Long-branch detection and pruning for gene trees.

Description
-----------
This script processes a collection of Newick gene trees and detects branches
whose branch lengths exceed a user-defined threshold. For decision-making, each
tree is internally rooted using predefined outgroups. By default, the outgroups
are tip names "2", "3", and "4".

The script classifies trees into two categories:

1. Correct trees:
   Trees without long branches after internal outgroup-based rooting.

2. Wrong trees:
   Trees containing at least one long branch.

For wrong trees, the script iteratively determines which leaves should be
removed according to a majority-rule strategy:

- If a long branch subtends more than half of all leaves in the tree, the long
  branch is retained and the other leaves are removed.
- Otherwise, all descendant leaves of the long branch are removed.

Important
---------
The outgroup-based rooting is used only for internal decision-making.
Final output trees remain based on the original unrooted tree structure, with
only the selected leaves removed.

Input format
------------
The input file specified by --trees supports two formats:

1. Single-column format:
   Each line contains one Newick tree.

   Example:
   ((1:0.1,2:0.2),(3:0.3,4:0.4));

2. Two-column tab-delimited format:
   The first column is a gene ID, and the second column is the Newick tree.

   Example:
   gene001    ((1:0.1,2:0.2),(3:0.3,4:0.4));

Usage examples
--------------
Basic usage:

    python prune_long_branch_trees.py \\
        --trees all_gene_trees.tsv \\
        --threshold 1.0

Specify output files:

    python prune_long_branch_trees.py \\
        --trees all_gene_trees.tsv \\
        --threshold 1.0 \\
        --correct 01.correct_trees.tsv \\
        --wrong 02.wrong_trees.tsv \\
        --pruned 03.pruned_trees.tsv \\
        --report 04.pruned_trees.report.tsv

Change the long-branch threshold:

    python prune_long_branch_trees.py \\
        --trees all_gene_trees.tsv \\
        --threshold 2.0

Requirements
------------
Python packages:
    ete3

Install ete3 if needed:

    conda install -c etetoolkit ete3

Outputs
-------
1. Correct tree file:
   Trees without long branches.

2. Wrong tree file:
   Original unrooted trees containing long branches.

3. Pruned tree file:
   Original unrooted trees after removing selected long-branch leaves.

4. Report file:
   Detailed pruning decisions for each tree.

5. Summary file:
   Overall counts of total trees, correct trees, wrong trees, and parse failures.
"""

import argparse
from itertools import combinations
from ete3 import Tree


# Default outgroups used only for internal rooting and decision-making.
OUTGROUPS = {"2", "3", "4"}


def parse_line_auto(line, line_no):
    """
    Automatically detect the input tree format.

    Supported formats:
    1. Single-column format:
       The entire line is treated as a Newick tree.

    2. Two-column tab-delimited format:
       gene_id <TAB> tree

    Parameters
    ----------
    line : str
        One input line.
    line_no : int
        Line number in the input file.

    Returns
    -------
    tuple
        (gene_id_or_none, tree_string)
    """
    raw = line.rstrip("\n")
    if not raw.strip():
        return None, None

    if "\t" in raw:
        parts = raw.split("\t")
        if len(parts) != 2:
            raise ValueError(f"Line {line_no} is not a valid two-column tab-delimited format.")
        gene_id = parts[0].strip()
        tree_str = parts[1].strip()
        if gene_id == "" or tree_str == "":
            raise ValueError(f"Line {line_no} has an empty gene_id or tree field.")
        return gene_id, tree_str

    return None, raw.strip()


def write_output_line(fh, gene_id, tree_str):
    """
    Write one output line.

    If the original input has a gene ID, output:
        gene_id <TAB> tree

    Otherwise, output:
        tree
    """
    if gene_id is None:
        fh.write(tree_str + "\n")
    else:
        fh.write(f"{gene_id}\t{tree_str}\n")


def get_present_outgroup_leaves(t, outgroups=OUTGROUPS):
    """
    Return outgroup leaves that are present in the current tree.
    """
    return [leaf for leaf in t.iter_leaves() if leaf.name in outgroups]


def get_subtree_leaf_count_for_nodes(t, leaves):
    """
    Return the number of leaves subtended by the MRCA of the given leaves.

    Parameters
    ----------
    t : ete3.Tree
        Input tree.
    leaves : list
        A list of ete3 leaf nodes.

    Returns
    -------
    tuple
        (number_of_leaves_under_MRCA, MRCA_node)
    """
    if len(leaves) == 1:
        return 1, leaves[0]
    mrca = t.get_common_ancestor(leaves)
    return len(mrca.get_leaves()), mrca


def choose_shortest_branch_leaf(leaves):
    """
    Choose the leaf with the shortest terminal branch length.
    """
    best = None
    best_dist = None

    for leaf in leaves:
        try:
            d = float(leaf.dist)
        except Exception:
            d = float("inf")

        if best is None or d < best_dist:
            best = leaf
            best_dist = d

    return best


def reroot_tree_majority_outgroups(t, outgroups=OUTGROUPS):
    """
    Root a tree using the available outgroups.

    Rooting rules:
    - 0 outgroups:
      Do not root the tree.

    - 1 outgroup:
      Root the tree using this single outgroup leaf.

    - 2 outgroups:
      If the two outgroups are monophyletic, root using their MRCA.
      Otherwise, root using the outgroup leaf with the shorter branch length.

    - 3 outgroups:
      If all three outgroups are monophyletic, root using their MRCA.
      If any monophyletic pair exists, root using the tightest monophyletic pair.
      If all three outgroups are separated, root using the outgroup leaf with the
      shortest branch length.
    """
    present = get_present_outgroup_leaves(t, outgroups)
    present_names = [x.name for x in present]

    if len(present) == 0:
        return t, [], "no_outgroup"

    if len(present) == 1:
        try:
            t.set_outgroup(present[0])
            return t, present_names, "rooted_1_outgroup"
        except Exception:
            return t, present_names, "root_failed_1"

    if len(present) == 2:
        a, b = present[0], present[1]
        try:
            count, mrca = get_subtree_leaf_count_for_nodes(t, [a, b])

            if count == 2:
                t.set_outgroup(mrca)
                return t, present_names, "rooted_2_outgroups_monophyletic"

            chosen = choose_shortest_branch_leaf([a, b])
            t.set_outgroup(chosen)
            return t, present_names, f"rooted_shorter_single_{chosen.name}"

        except Exception:
            return t, present_names, "root_failed_2"

    # Three outgroups are present.
    try:
        all_count, all_mrca = get_subtree_leaf_count_for_nodes(t, present)
        if all_count == 3:
            t.set_outgroup(all_mrca)
            return t, present_names, "rooted_3_outgroups_monophyletic"
    except Exception:
        pass

    monophyletic_pairs = []

    for a, b in combinations(present, 2):
        try:
            count, mrca = get_subtree_leaf_count_for_nodes(t, [a, b])
            if count == 2:
                monophyletic_pairs.append((a, b, mrca))
        except Exception:
            continue

    if len(monophyletic_pairs) > 0:
        best_pair = None
        best_score = None
        best_mrca = None

        for a, b, mrca in monophyletic_pairs:
            try:
                da = float(a.dist)
            except Exception:
                da = float("inf")

            try:
                db = float(b.dist)
            except Exception:
                db = float("inf")

            score = da + db

            if best_score is None or score < best_score:
                best_score = score
                best_pair = (a.name, b.name)
                best_mrca = mrca

        try:
            t.set_outgroup(best_mrca)
            return t, present_names, f"rooted_majority_pair_{best_pair[0]}_{best_pair[1]}"
        except Exception:
            return t, present_names, "root_failed_majority_pair"

    chosen = choose_shortest_branch_leaf(present)

    try:
        t.set_outgroup(chosen)
        return t, present_names, f"rooted_shortest_single_{chosen.name}"
    except Exception:
        return t, present_names, "root_failed_3_separate"


def cleanup_upward(node):
    """
    Clean up tree structure upward after deleting a leaf or node.

    This function removes empty nodes and collapses internal nodes with only
    one child. When a single-child internal node is collapsed, its branch length
    is added to the child branch length when possible.
    """
    while node is not None:
        if node.is_root():
            break

        children = node.get_children()

        if len(children) == 0:
            parent = node.up
            node.detach()
            node = parent

        elif len(children) == 1:
            child = children[0]
            parent = node.up

            try:
                child.dist += node.dist
            except Exception:
                pass

            child.detach()
            parent.add_child(child)
            node.detach()
            node = parent

        else:
            break


def cleanup_root(tree):
    """
    Clean up the root if it has only one child.

    This avoids producing trees with unnecessary single-child root structures
    after pruning.
    """
    while True:
        children = tree.get_children()
        if len(children) != 1:
            break

        child = children[0]

        tree.name = child.name

        try:
            tree.dist = child.dist
        except Exception:
            pass

        try:
            tree.support = child.support
        except Exception:
            pass

        grandchildren = child.get_children()
        child.detach()

        for gc in grandchildren:
            gc.detach()
            tree.add_child(gc)


def prune_leaves_and_cleanup(t, leaf_names_to_remove):
    """
    Remove selected leaves from a tree and clean up redundant internal nodes.
    """
    remove_set = set(leaf_names_to_remove)
    targets = [leaf for leaf in t.iter_leaves() if leaf.name in remove_set]

    for leaf in targets:
        parent = leaf.up
        leaf.detach()

        if parent is not None:
            cleanup_upward(parent)

    cleanup_root(t)

    return t


def collect_long_branch_pruning_targets_majority_rule(t, threshold):
    """
    Identify leaves to remove according to long-branch and majority-rule logic.

    For each branch with dist > threshold:

    - If the number of descendant leaves under this long branch is greater than
      half of all leaves in the tree, keep this long-branch subtree and remove
      all other leaves.

    - Otherwise, remove all descendant leaves under this long branch.

    Both internal branches and terminal branches are checked because node.dist
    is evaluated for every non-root node.
    """
    all_leaves = [leaf.name for leaf in t.iter_leaves()]
    all_leaf_set = set(all_leaves)
    total_leaves = len(all_leaves)

    remove_leaves = set()
    decision_records = []

    for node in t.traverse("preorder"):
        if node.is_root():
            continue

        try:
            dist = float(node.dist)
        except Exception:
            continue

        if dist <= threshold:
            continue

        desc_leaves = set(leaf.name for leaf in node.iter_leaves())
        desc_count = len(desc_leaves)

        if desc_count > total_leaves / 2:
            to_remove = all_leaf_set - desc_leaves
            strategy = "keep_long_branch_remove_others"
        else:
            to_remove = desc_leaves
            strategy = "remove_long_branch_descendants"

        remove_leaves |= to_remove

        decision_records.append({
            "node_name": node.name if node.name else "NA",
            "dist": dist,
            "desc_count": desc_count,
            "total_leaves": total_leaves,
            "strategy": strategy,
            "affected_leaf_count": len(to_remove),
            "affected_leaves": sorted(to_remove),
        })

    return remove_leaves, decision_records


def has_any_long_branch(t, threshold):
    """
    Check whether a tree contains any branch longer than the threshold.
    """
    for node in t.traverse("preorder"):
        if node.is_root():
            continue

        try:
            if float(node.dist) > threshold:
                return True
        except Exception:
            continue

    return False


def iterative_root_and_collect_removals(t_raw, threshold, max_rounds=100):
    """
    Iteratively root the working tree and collect leaves that should be removed.

    Important:
    - Rooting and long-branch decisions are made only on a working copy.
    - The original unrooted tree is not modified in this function.
    - The final removal set is later applied to a copy of the original unrooted tree.

    Returns
    -------
    tuple
        all_removed_leaves :
            A set of leaf names to remove.

        all_decisions :
            Detailed decision records for each iteration.

        round_count :
            Number of pruning iterations.

        final_has_long_branch :
            Whether the final working tree still contains long branches.

        root_history :
            Rooting status for each iteration.
    """
    t_work = t_raw.copy(method="deepcopy")

    all_removed = set()
    all_decisions = []
    root_history = []
    round_count = 0

    for r in range(1, max_rounds + 1):
        t_work, used_outgroups, root_status = reroot_tree_majority_outgroups(t_work, OUTGROUPS)

        root_history.append({
            "round": r,
            "used_outgroups": used_outgroups,
            "root_status": root_status
        })

        remove_leaves, decision_records = collect_long_branch_pruning_targets_majority_rule(
            t_work,
            threshold
        )

        if len(decision_records) == 0:
            final_has_long = has_any_long_branch(t_work, threshold)
            return all_removed, all_decisions, round_count, final_has_long, root_history

        round_count += 1

        all_decisions.append({
            "round": r,
            "decision_records": decision_records
        })

        new_to_remove = set(remove_leaves) - all_removed

        if len(new_to_remove) == 0:
            final_has_long = has_any_long_branch(t_work, threshold)
            return all_removed, all_decisions, round_count, final_has_long, root_history

        all_removed |= set(remove_leaves)
        t_work = prune_leaves_and_cleanup(t_work, remove_leaves)

    final_has_long = has_any_long_branch(t_work, threshold)

    return all_removed, all_decisions, round_count, final_has_long, root_history


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Detect and prune long branches in gene trees. "
            "Trees are internally rooted using outgroups 2/3/4 for decision-making, "
            "but final outputs remain based on the original unrooted trees."
        )
    )

    parser.add_argument(
        "--trees",
        required=True,
        help="Input tree file. Supports single-column Newick trees or two-column tab-delimited gene_id + tree format."
    )

    parser.add_argument(
        "--threshold",
        type=float,
        default=1.0,
        help="Branch-length threshold. Branches with dist > threshold are treated as long branches. Default: 1.0."
    )

    parser.add_argument(
        "--correct",
        default="01.correct_trees.tsv",
        help="Output file for trees without long branches. Default: 01.correct_trees.tsv."
    )

    parser.add_argument(
        "--wrong",
        default="02.wrong_trees.tsv",
        help="Output file for original unrooted trees containing long branches. Default: 02.wrong_trees.tsv."
    )

    parser.add_argument(
        "--pruned",
        default="03.pruned_trees.tsv",
        help="Output file for pruned trees based on the original unrooted trees. Default: 03.pruned_trees.tsv."
    )

    parser.add_argument(
        "--report",
        default="04.pruned_trees.report.tsv",
        help="Output detailed report file. Default: 04.pruned_trees.report.tsv."
    )

    args = parser.parse_args()

    total = 0
    parse_failed = 0
    correct_count = 0
    wrong_count = 0

    with open(args.trees, "r", encoding="utf-8") as fin, \
         open(args.correct, "w", encoding="utf-8") as fout_correct, \
         open(args.wrong, "w", encoding="utf-8") as fout_wrong, \
         open(args.pruned, "w", encoding="utf-8") as fout_pruned, \
         open(args.report, "w", encoding="utf-8") as freport:

        freport.write(
            "line_no\tgene_id\tstatus\tthreshold\t"
            "initial_long_branch_node_count\titerations\tfinal_has_long_branch\t"
            "rooting_summary\tremoved_leaf_count\tremoved_leaf_names\tdecision_summary\n"
        )

        for idx, line in enumerate(fin, start=1):
            if not line.strip():
                continue

            total += 1

            try:
                gene_id, tree_str = parse_line_auto(line, idx)
            except Exception as e:
                parse_failed += 1

                freport.write(
                    f"{idx}\tNA\tparse_failed\t{args.threshold}\t"
                    f"NA\tNA\tNA\tNA\tNA\tNA\t{str(e)}\n"
                )

                continue

            gene_id_out = gene_id if gene_id is not None else f"LINE_{idx}"

            try:
                t_raw = Tree(tree_str, format=1)
            except Exception as e:
                parse_failed += 1

                freport.write(
                    f"{idx}\t{gene_id_out}\tparse_failed\t{args.threshold}\t"
                    f"NA\tNA\tNA\tNA\tNA\tNA\tTree parsing failed: {str(e)}\n"
                )

                continue

            # Initial long-branch detection is performed on a temporary rooted tree.
            # However, all correct/wrong tree outputs keep the original unrooted form.
            t_init = t_raw.copy(method="deepcopy")
            t_init, init_outgroups, init_root_status = reroot_tree_majority_outgroups(t_init, OUTGROUPS)

            initial_remove_leaves, initial_decisions = collect_long_branch_pruning_targets_majority_rule(
                t_init,
                args.threshold
            )

            initial_long_branch_count = len(initial_decisions)

            rooting_summary = (
                f"init:{init_root_status}"
                f"[{','.join(init_outgroups) if init_outgroups else 'NA'}]"
            )

            # Correct tree: output the original unrooted tree.
            if initial_long_branch_count == 0:
                correct_count += 1

                correct_tree_str = t_raw.write(format=1).strip()
                write_output_line(fout_correct, gene_id, correct_tree_str)

                freport.write(
                    f"{idx}\t{gene_id_out}\tcorrect\t{args.threshold}\t"
                    f"0\t0\tFalse\t{rooting_summary}\t0\tNA\tNA\n"
                )

                continue

            # Wrong tree: output the original unrooted tree.
            wrong_count += 1

            wrong_tree_str = t_raw.write(format=1).strip()
            write_output_line(fout_wrong, gene_id, wrong_tree_str)

            # Internally root the tree, determine leaves to remove, and iterate.
            all_removed, all_decisions, round_count, final_has_long_branch, root_history = (
                iterative_root_and_collect_removals(t_raw, args.threshold)
            )

            # Apply the final removal set to the original unrooted tree for output.
            t_pruned_raw = t_raw.copy(method="deepcopy")

            if len(all_removed) > 0:
                t_pruned_raw = prune_leaves_and_cleanup(t_pruned_raw, all_removed)

            pruned_tree_str = t_pruned_raw.write(format=1).strip()
            write_output_line(fout_pruned, gene_id, pruned_tree_str)

            root_parts = []

            for rh in root_history:
                root_parts.append(
                    f"round{rh['round']}:{rh['root_status']}"
                    f"[{','.join(rh['used_outgroups']) if rh['used_outgroups'] else 'NA'}]"
                )

            decision_summary_parts = []

            for round_info in all_decisions:
                r = round_info["round"]

                for d in round_info["decision_records"]:
                    decision_summary_parts.append(
                        f"round{r}:{d['strategy']}|dist={d['dist']}|desc={d['desc_count']}/{d['total_leaves']}"
                    )

            freport.write(
                f"{idx}\t{gene_id_out}\twrong\t{args.threshold}\t"
                f"{initial_long_branch_count}\t{round_count}\t{final_has_long_branch}\t"
                f"{';'.join(root_parts) if root_parts else rooting_summary}\t"
                f"{len(all_removed)}\t"
                f"{','.join(sorted(all_removed)) if all_removed else 'NA'}\t"
                f"{';'.join(decision_summary_parts) if decision_summary_parts else 'NA'}\n"
            )

    summary_file = args.report + ".summary.txt"

    with open(summary_file, "w", encoding="utf-8") as fsum:
        fsum.write(f"Total trees\t{total}\n")
        fsum.write(f"Parse failed\t{parse_failed}\n")
        fsum.write(f"Correct\t{correct_count}\n")
        fsum.write(f"Wrong\t{wrong_count}\n")
        fsum.write(f"Threshold\t{args.threshold}\n")
        fsum.write(
            "Rule\tRoot by outgroups 2/3/4 only for internal decision-making; "
            "final outputs remain based on original unrooted trees. "
            "If two outgroups are not monophyletic, use the shorter-branch outgroup to root. "
            "If three outgroups are all separate with no monophyletic pair, use the shortest-branch outgroup to root. "
            "If three outgroups are not monophyletic but one pair is monophyletic, use that majority pair. "
            "If no outgroups exist, prune directly. "
            "Iterate pruning until no new long branches remain in the rooted working tree.\n"
        )

    print("Done.")
    print(f"Input:        {args.trees}")
    print(f"Threshold:    {args.threshold}")
    print(f"Correct:      {args.correct}")
    print(f"Wrong:        {args.wrong}")
    print(f"Pruned:       {args.pruned}")
    print(f"Report:       {args.report}")
    print(f"Summary:      {summary_file}")


if __name__ == "__main__":
    main()