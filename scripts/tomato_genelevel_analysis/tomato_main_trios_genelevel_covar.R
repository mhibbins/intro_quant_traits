remove(list=ls())
library(ggplot2)
library(ape)
library(tidyverse)
library(gridExtra)
                            ### Load datasets ### 

#This script performs the gene-level analysis and generates the chi-squared
#tables reports in Figure 5 and Supplementary Tables 1 and 2.

#Gene expression data

ovule_exp_cpm <- read.csv(
  "C:/Users/18126/OneDrive - Indiana University/Projects/intro_gene_expression/datasets/ovule_exp_cpm.csv")

#Gene trees from Pease et al. 

tree_names <- read.table(
  "C:/Users/18126/OneDrive - Indiana University/Projects/intro_gene_expression/datasets/pease_tomato/james_trees_pruned.names.txt", quote="\"", comment.char="")

tree_names <- as.character(tree_names[,1])

intro_tree_topologies <- read.csv(
  "C:/Users/18126/OneDrive - Indiana University/Projects/intro_gene_expression/datasets/pease_tomato/pease_intro_trio_tree_topologies.csv")
cont_tree_topologies <- read.csv(
  "C:/Users/18126/OneDrive - Indiana University/Projects/intro_gene_expression/datasets/pease_tomato/pease_cont_trio_tree_topologies.csv")

                            ### Data tidying ### 

#Calculate means from individual replicates in expression data

LA3475 <- rowMeans(ovule_exp_cpm[,c("LA101A", "LA101B","LA101.C.A", "LA101.C.B")])
LA1316 <- rowMeans(ovule_exp_cpm[,c("LA1316.4A", "LA1316.4B", "LA1316.9")])
LA1589 <- rowMeans(ovule_exp_cpm[,c("LA1589.1", "LA1589.7")])
LA1777 <- ovule_exp_cpm[,"LA1777.21"]
LA2172 <- rowMeans(ovule_exp_cpm[,c("LA2172.19", "LA2172.9")])
LA3778 <- rowMeans(ovule_exp_cpm[,c("LA3778.226", "LA3778.8")])
LA4117 <- rowMeans(ovule_exp_cpm[,c("LA4117A.33", "LA4117A.38")])
LA0716 <- rowMeans(ovule_exp_cpm[,c("LA716A.1", "LA716A.2", "LA716.B", "LA716.C_S2")])

#Make tidy dataframe for gene expression data

tidy_ovule_exp_cpm <- as.data.frame(cbind("Gene" = ovule_exp_cpm$Gene, LA3475, LA1316, 
                                          LA1589, LA1777, LA2172, LA3778, LA4117, LA0716))

#Get overlapping set of gene tree topologies and expression orthologs

tidy_ovule_exp_cpm$Gene <- gsub('.{4}$', '', tidy_ovule_exp_cpm$Gene)
intro_tree_topologies$gene_name <- gsub('.{2}$', '', intro_tree_topologies$gene_name)
cont_tree_topologies$gene_name <- gsub('.{2}$', '', cont_tree_topologies$gene_name)

topologies_overlap <- intersect(cont_tree_topologies$gene_name, 
                                intro_tree_topologies$gene_name)
tidy_ovule_exp_cpm <- subset(tidy_ovule_exp_cpm,
                             Gene %in% topologies_overlap)
intro_tree_topologies <- subset(intro_tree_topologies, 
                                gene_name %in% tidy_ovule_exp_cpm$Gene)
cont_tree_topologies <- subset(cont_tree_topologies, 
                                gene_name %in% tidy_ovule_exp_cpm$Gene)

                ### Get expression tree for each ortholog ###

get_exp_tree <- function(gene, P1, P2, P3) {
  
  trio <- as.numeric(gene[c(P1, P2, P3)])
  dP1P2 <- abs(trio[1] - trio[2])
  dP1P3 <- abs(trio[1] - trio[3])
  dP2P3 <- abs(trio[2] - trio[3])
  
  if(dP1P2 < dP1P3 && dP1P2 < dP2P3) {return("P1P2")}
  else if(dP2P3 < dP1P2 && dP2P3 < dP1P3) {return("P2P3")}
  else if(dP1P3 < dP1P2 && dP1P3 < dP2P3) {return("P1P3")}
}

intro_trio_min <- apply(tidy_ovule_exp_cpm, 1, get_exp_tree,
                              "LA3778", "LA1777", "LA1316")
control_trio_min <- apply(tidy_ovule_exp_cpm, 1, get_exp_tree,
                                "LA1589", "LA3475", "LA0716")


                          ### Get results ### 

gene_name <- cont_tree_topologies$gene_name
cont_topology <- cont_tree_topologies$topology
intro_topology <- intro_tree_topologies$topology

results_df <- data.frame(cbind(gene_name, cont_topology, intro_topology,
                    control_trio_min, intro_trio_min))

cont_table <- table(results_df[ , c("cont_topology", "control_trio_min")])
intro_table <- table(results_df[ , c("intro_topology", "intro_trio_min")])

cont_result <- chisq.test(cont_table)
intro_result <- chisq.test(intro_table)

write.csv(results_df,
"C:/Users/18126/OneDrive - Indiana University/Projects/intro_gene_expression/results/gene_expression_trees.csv")





