# -*- coding: utf-8 -*-
"""
Created on Fri Nov 13 12:38:38 2020

@author: Mark
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

#This script was used to generate the individual 3x3 plots in Figure 3 and Supplementary Figure 1. 
#Change the input file path and corresponding arguments to seaborn to generate the plot
for each set of conditions. 

                      ### Parse the results file ### 

linelist = []

with open(
"C:/Users/18126/OneDrive - Indiana University/Projects/intro_gene_expression/results/tomato_power_analysis_v2_10kgenes_results.txt") as results_file:
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


delta = [0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01,
         0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
         0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
deltas = [delta for i in range(100)]
deltas = [item for sublist in deltas for item in sublist]

timing = [0.5, 0.5, 0.5, 0.25, 0.25, 0.25, 0.1, 0.1, 0.1,
           0.5, 0.5, 0.5, 0.25, 0.25, 0.25, 0.1, 0.1, 0.1,
           0.5, 0.5, 0.5, 0.25, 0.25, 0.25, 0.1, 0.1, 0.1]
timings = [timing for i in range(100)]
timings = [item for sublist in timings for item in sublist]

intbranch = [1, 0.5, 0.1, 1, 0.5, 0.1, 1, 0.5, 0.1,
             1, 0.5, 0.1, 1, 0.5, 0.1, 1, 0.5, 0.1,
             1, 0.5, 0.1, 1, 0.5, 0.1, 1, 0.5, 0.1]
intbranches = [intbranch for i in range(100)]
intbranches = [item for sublist in intbranches for item in sublist]

                     ### Set up power matrices ### 

#rows are delta, columns are timing 

intbranch_1_power_matrix = [[0,0,0],[0,0,0],[0,0,0]]
intbranch_05_power_matrix = [[0,0,0],[0,0,0],[0,0,0]]
intbranch_01_power_matrix = [[0,0,0],[0,0,0],[0,0,0]]

for i in range(len(pvals)):
    if float(pvals[i]) < 0.05 and signs[i] == "Yes": #sign is backwards
        if intbranches[i] == 1:
            if deltas[i] == 0.01:
                if timings[i] == 0.5:
                    intbranch_1_power_matrix[0][0] += 1
                elif timings[i] == 0.25:
                    intbranch_1_power_matrix[0][1] += 1
                elif timings[i] == 0.1:
                    intbranch_1_power_matrix[0][2] += 1
            elif deltas[i] == 0.05:
                if timings[i] == 0.5:
                    intbranch_1_power_matrix[1][0] += 1
                elif timings[i] == 0.25:
                    intbranch_1_power_matrix[1][1] += 1
                elif timings[i] == 0.1:
                    intbranch_1_power_matrix[1][2] += 1
            elif deltas[i] == 0.1:
                if timings[i] == 0.5:
                    intbranch_1_power_matrix[2][0] += 1
                elif timings[i] == 0.25:
                    intbranch_1_power_matrix[2][1] += 1
                elif timings[i] == 0.1:
                    intbranch_1_power_matrix[2][2] += 1
        elif intbranches[i] == 0.5:
            if deltas[i] == 0.01:
                if timings[i] == 0.5:
                    intbranch_05_power_matrix[0][0] += 1
                elif timings[i] == 0.25:
                    intbranch_05_power_matrix[0][1] += 1
                elif timings[i] == 0.1:
                    intbranch_05_power_matrix[0][2] += 1
            elif deltas[i] == 0.05:
                if timings[i] == 0.5:
                    intbranch_05_power_matrix[1][0] += 1
                elif timings[i] == 0.25:
                    intbranch_05_power_matrix[1][1] += 1
                elif timings[i] == 0.1:
                    intbranch_05_power_matrix[1][2] += 1
            elif deltas[i] == 0.1:
                if timings[i] == 0.5:
                    intbranch_05_power_matrix[2][0] += 1
                elif timings[i] == 0.25:
                    intbranch_05_power_matrix[2][1] += 1
                elif timings[i] == 0.1:
                    intbranch_05_power_matrix[2][2] += 1
        elif intbranches[i] == 0.1:
            if deltas[i] == 0.01:
                if timings[i] == 0.5:
                    intbranch_01_power_matrix[0][0] += 1
                elif timings[i] == 0.25:
                    intbranch_01_power_matrix[0][1] += 1
                elif timings[i] == 0.1:
                    intbranch_01_power_matrix[0][2] += 1
            elif deltas[i] == 0.05:
                if timings[i] == 0.5:
                    intbranch_01_power_matrix[1][0] += 1
                elif timings[i] == 0.25:
                    intbranch_01_power_matrix[1][1] += 1
                elif timings[i] == 0.1:
                    intbranch_01_power_matrix[1][2] += 1
            elif deltas[i] == 0.1:
                if timings[i] == 0.5:
                    intbranch_01_power_matrix[2][0] += 1
                elif timings[i] == 0.25:
                    intbranch_01_power_matrix[2][1] += 1
                elif timings[i] == 0.1:
                    intbranch_01_power_matrix[2][2] += 1
                    
for i in range(len(intbranch_1_power_matrix)):
    for j in range(len(intbranch_1_power_matrix)):
        intbranch_1_power_matrix[i][j] = intbranch_1_power_matrix[i][j]/100
        intbranch_05_power_matrix[i][j] = intbranch_05_power_matrix[i][j]/100
        intbranch_01_power_matrix[i][j] = intbranch_01_power_matrix[i][j]/100

 
    
                             ### Plots ### 
xticks = [0.5, 0.25, 0.1]
yticks = [0.01, 0.05, 0.1]

heatmap = sns.heatmap(intbranch_1_power_matrix, xticklabels=xticks,
                  yticklabels=yticks, annot=True, fmt="g", 
                  cmap = "YlGnBu", cbar = False)
plt.yticks(rotation=0)
plt.xlabel("$t_{1} - t_{m}$", fontsize=20)
plt.ylabel("Rate of introgression", fontsize=20)
plt.title("$t_{2} - t_{1} = 1$", fontsize=20)
plt.gca().invert_xaxis()
heatmap.set_xticklabels(heatmap.get_xmajorticklabels(), fontsize = 15)
heatmap.set_yticklabels(heatmap.get_ymajorticklabels(), fontsize = 15)
plt.show()


        
