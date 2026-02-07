01.merge_wgdi_csv.py – Merge syntenic gene sets from all species.

02.filter_syntenic_matrix_by_species.py – Simplify the original syntenic matrix by keeping only one syntenic gene per species.

03.extract_acek_sequences.py – Retrieve sequences of all species for each syntenic gene.

04.Construct_gene_trees.sh – Build gene trees from syntenic genes.

05.filter-orders.py – Filter gene tree sets based on the number of target taxa to generate a new set of trees.

05.filter-species.py (optional) – Filter gene tree sets based on the number of species to generate a new set of trees.

06.prune_tree_nodes.py (optional) – Remove problematic nodes from gene trees if necessary.