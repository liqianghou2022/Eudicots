#!/usr/bin/env python3
"""
Script: 06.prune_tree_nodes.py
Author: Liqiang Hou
Date: 2025.08

Description:
    Remove specified nodes from phylogenetic trees while maintaining tree structure.
    When a node is removed, its siblings are reconnected to preserve tree topology.
    Branch lengths are adjusted accordingly.

Usage:
    python 06.prune_tree_nodes.py
    
    Before running, edit the following variables in the script:
    - input_file: Path to input Newick file
    - output_file: Path to output Newick file  
    - nodes_to_remove: Set of node names to remove

Input:
    - Newick format tree file (one or more trees, one per line)
    - Trees should have node names in format 1 (with names and branch lengths)

Output:
    - Pruned Newick file with specified nodes removed
    - Branch lengths are preserved and adjusted when nodes are deleted
    - Format: Newick with internal node names and branch lengths

Algorithm:
    1. For each node to be removed:
       - If node has one sibling: sibling is connected to grandparent
       - Branch length is adjusted (parent.dist + sibling.dist)
       - If no grandparent exists, sibling becomes new root
    2. Tree topology is maintained after node removal

Example:
    Original tree: ((A:0.1,B:0.2)node1:0.3,(C:0.4,D:0.5)node2:0.6)root;
    Remove node B: ((A:0.1)node1:0.3,(C:0.4,D:0.5)node2:0.6)root;
    
    The branch length from A to node1 remains 0.1
    Tree structure is preserved

Requirements:
    - Python 3.6+
    - ete3 (pip install ete3)

Notes:
    - Node names are case-sensitive
    - Nodes not found in the tree are silently skipped
    - When removing a node with one sibling, the sibling inherits the position
    - Branch lengths are summed when reconnecting nodes
    - Output uses format 5: includes both node names and branch lengths
    - Branch lengths use 10 decimal places to avoid scientific notation
"""

from ete3 import Tree

# 输入和输出文件
input_file = "all_nwk.tre"   # 你的 Newick 文件
output_file = "all_prune_nwk.tre"  # 处理后的 Newick 文件

# 需要删除的节点列表
nodes_to_remove = {"Node1","Node2","Node4".....}   #If you want to delete nodes:"Node1","Node2","Node4".....

# 读取文件，可能包含多个树
with open(input_file, "r") as f:
    tree_strings = f.read().strip().split("\n")  # 逐行读取每棵树

cleaned_trees = []

# 逐棵树进行处理
for tree_str in tree_strings:
    t = Tree(tree_str, format=1)  # 解析树

    # 删除指定的节点并重新连接树枝
    for node_name in nodes_to_remove:
        nodes = t.search_nodes(name=node_name)
        for node in nodes:
            parent = node.up  # 获取父节点
            if parent:
                siblings = node.get_sisters()  # 获取兄弟节点
                if len(siblings) == 1:
                    # 只有一个兄弟节点，则直接连到祖先节点，防止孤立
                    grandparent = parent.up
                    if grandparent:
                        siblings[0].detach()
                        grandparent.add_child(siblings[0], dist=parent.dist + siblings[0].dist)
                    else:
                        siblings[0].detach()
                        t = siblings[0]  # 重新设置树的根
                node.delete()

    # 以标准格式输出，并保留枝长（避免科学计数法）
    cleaned_trees.append(t.write(format=5, dist_formatter="%.10f"))

# 将处理后的所有树写入新文件
with open(output_file, "w") as outf:
    outf.write("\n".join(cleaned_trees) + "\n")

print("处理完成! 新文件保存为:", output_file)