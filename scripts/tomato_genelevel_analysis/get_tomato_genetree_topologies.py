# -*- coding: utf-8 -*-
"""
Created on Tue Nov 17 11:44:49 2020

@author: Mark
"""

#This script was used to parse the coding sequence tree topology for each gene 
#in each trio.

from ete3 import Tree

treelist = []
treenames = []

with open(
"C:/Users/18126/OneDrive - Indiana University/Projects/intro_gene_expression/datasets/pease_tomato/james_trees_pruned.trees.txt") as genetrees:
    
    for tree in genetrees: 
        treelist.append(Tree(tree))
        
with open(
"C:/Users/18126/OneDrive - Indiana University/Projects/intro_gene_expression/datasets/pease_tomato/james_trees_pruned.names.txt") as names:
    
    for name in names: 
        treenames.append(name.strip())
        
intro_trio_topologies = []
cont_trio_topologies = []
intro_trio_treenames = []
cont_trio_treenames = []

for i in range(len(treelist)):
    
    
    intro_P1P2_dist = treelist[i].get_distance("LA3778", "LA1777", topology_only = True)
    intro_P1P3_dist = treelist[i].get_distance("LA3778", "LA1316", topology_only = True)
    intro_P2P3_dist = treelist[i].get_distance("LA1777", "LA1316", topology_only = True)
    cont_P1P2_dist = treelist[i].get_distance("LA1589", "LA3475", topology_only = True)
    cont_P1P3_dist = treelist[i].get_distance("LA1589", "LA0716", topology_only = True)
    cont_P2P3_dist = treelist[i].get_distance("LA3475", "LA0716", topology_only = True)
    
    if intro_P1P2_dist < intro_P1P3_dist and intro_P1P2_dist < intro_P2P3_dist:
        intro_trio_topologies.append("P1P2")
        intro_trio_treenames.append(treenames[i])
    elif intro_P1P3_dist < intro_P1P2_dist and intro_P1P3_dist < intro_P2P3_dist:
        intro_trio_topologies.append("P1P3")
        intro_trio_treenames.append(treenames[i])
    elif intro_P2P3_dist < intro_P1P2_dist and intro_P2P3_dist < intro_P1P3_dist:
        intro_trio_topologies.append("P2P3")
        intro_trio_treenames.append(treenames[i])
        
    if cont_P1P2_dist < cont_P1P3_dist and cont_P1P2_dist < cont_P2P3_dist:
        cont_trio_topologies.append("P1P2")
        cont_trio_treenames.append(treenames[i])
    elif cont_P1P3_dist < cont_P1P2_dist and cont_P1P3_dist < cont_P2P3_dist:
        cont_trio_topologies.append("P1P3")
        cont_trio_treenames.append(treenames[i])
    elif cont_P2P3_dist < cont_P1P2_dist and cont_P2P3_dist < cont_P1P3_dist:
        cont_trio_topologies.append("P2P3")
        cont_trio_treenames.append(treenames[i])
        
        
with open(
"C:/Users/18126/OneDrive - Indiana University/Projects/intro_gene_expression/datasets/pease_tomato/pease_intro_trio_tree_topologies.csv", "w") as intro_topologies:
    intro_topologies.write("gene_name, topology\n")
    for i in range(len(intro_trio_topologies)):
        line = ", ".join([intro_trio_treenames[i], intro_trio_topologies[i]])
        intro_topologies.write(line + "\n")
        
with open(
"C:/Users/18126/OneDrive - Indiana University/Projects/intro_gene_expression/datasets/pease_tomato/pease_cont_trio_tree_topologies.csv", "w") as cont_topologies:
    cont_topologies.write("gene_name, topology\n")
    for i in range(len(cont_trio_topologies)):
        line = ", ".join([cont_trio_treenames[i], cont_trio_topologies[i]])
        cont_topologies.write(line + "\n")

with open(
"C:/Users/18126/OneDrive - Indiana University/Projects/intro_gene_expression/datasets/pease_tomato/pease_cont_trio_tree_divergences.csv", "w") as cont_divergences:
    cont_divergences.write("gene_name, P1P2, P1P3, P2P3\n")
    for i in range(len(cont_trio_divergences)):
        line = ", ".join([cont_trio_treenames[i], str(cont_trio_divergences[i][0]),
                          str(cont_trio_divergences[i][1]), str(cont_trio_divergences[i][2])])
        print(line)
        cont_divergences.write(line + "\n")
