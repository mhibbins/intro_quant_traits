# -*- coding: utf-8 -*-
"""
Created on Fri Nov 13 12:38:38 2020

@author: Mark
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

                      ### Parse the results file ### 

linelist = []

with open(
"C:/Users/18126/OneDrive - Indiana University/Projects/intro_gene_expression/Revision/Results/tomato_power_analysis_nointro_results.txt") as results_file:
    for line in results_file:
        linelist.append(line.strip())
        
        
pvals = []
signs = []

def isfloat(value):
    try:
        float(value)
        return True
    except ValueError:
        return False

for i in range(len(linelist)):
    
    split = linelist[i].split()
    if len(split) > 1:
    
        if split[0] == "[1]":
            pvals.append([])
            signs.append([])

        if isfloat(split[1]):
            pvals[-1].append(split[1])
            signs[-1].append(split[2])
            pvals[-1].append(split[3])
        else:
            signs[-1].append(split[1])
            pvals[-1].append(split[2])
            signs[-1].append(split[3])
        
              ### Add long format simulation parameters ### 
        
pvals = [item for sublist in pvals for item in sublist]
signs = [item for sublist in signs for item in sublist]


n_gene = [5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000, 5000,
           10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000, 10000,
           15000, 15000, 15000, 15000, 15000, 15000, 15000, 15000, 15000]
n_genes = [n_gene for i in range(100)]
n_genes = [item for sublist in n_genes for item in sublist]


intbranch = [1, 0.5, 0.1, 1, 0.5, 0.1, 1, 0.5, 0.1,
             1, 0.5, 0.1, 1, 0.5, 0.1, 1, 0.5, 0.1,
             1, 0.5, 0.1, 1, 0.5, 0.1, 1, 0.5, 0.1]
intbranches = [intbranch for i in range(100)]
intbranches = [item for sublist in intbranches for item in sublist]

                     ### Set up power matrices ### 

#rows are delta, columns are timing 

power_matrix = [[0,0,0],[0,0,0],[0,0,0]]


for i in range(len(pvals)):
    if float(pvals[i]) < 0.05 and signs[i] == "Yes": #sign is backwards
        if n_genes[i] == 5000:
            if intbranches[i] == 1:
                power_matrix[0][0] += 1
            elif intbranches[i] == 0.5:
                power_matrix[0][1] += 1
            elif intbranches[i] == 0.1:
                power_matrix[0][2] += 1
        elif n_genes[i] == 10000:
            if intbranches[i] == 1:
                power_matrix[1][0] += 1
            elif intbranches[i] == 0.5:
                power_matrix[1][1] += 1
            elif intbranches[i] == 0.1:
                power_matrix[1][2] += 1
        elif n_genes[i] == 15000:
            if intbranches[i] == 1:
                power_matrix[2][0] += 1
            elif intbranches[i] == 0.5:
                power_matrix[2][1] += 1
            elif intbranches[i] == 0.1:
                power_matrix[2][2] += 1
        
                    
for i in range(len(power_matrix)):
    for j in range(len(power_matrix)):
        power_matrix[i][j] = power_matrix[i][j]/100
 
    
                             ### Plot ### 
xticks = [5000, 10000, 15000]
yticks = [1, 0.5, 0.1]

heatmap = sns.heatmap(power_matrix, xticklabels=xticks,
                  yticklabels=yticks, annot=True, fmt="g", 
                  cmap = ["lightyellow"], linecolor = "k",
                  cbar = False)

plt.yticks(rotation=0)
plt.xlabel("Number of traits", fontsize=20)
plt.ylabel("$t_{2} - t_{1}$", fontsize=20)
plt.title("No introgression", fontsize=20)
#plt.gca().invert_xaxis()
heatmap.set_xticklabels(heatmap.get_xmajorticklabels(), fontsize = 15)
heatmap.set_yticklabels(heatmap.get_ymajorticklabels(), fontsize = 15)
plt.show()


        
