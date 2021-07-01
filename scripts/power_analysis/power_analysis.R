remove(list=ls())
library(tidyverse)
#library(MASS)
library(tmvtnorm)
library(matrixcalc)
library(corpcor)

# This script simulates a single replicate gene expression dataset 
# from the specified input parameters and tests for 
# significantly different expression covariance across 27 combinations of 
# delta, tm, and t1. Used to generate the results plotted in Figure 3
# and Supplementary Figure 1.

sigma2 = 1
n_genes = 15000

### Functions for expected covariances ### 

calc_cov_P1P2 <- function(t2, t1, tm, delta) {
  
  part1 = (1-delta)*(1 - exp(-(t2-t1)))*((exp(t2)*(t2-t1))/(exp(t2)-exp(t1)))
  part2 = (1-delta)*((1/3)*exp(-(t2-t1)))
  part3 = delta*((1/3)*exp(-(t2-tm)))
  P1P2_cov = part1 + part2 + part3
  
  return(P1P2_cov)
}

#Assuming P2/P3 introgression 

calc_cov_P2P3 <- function(t2, t1, tm, delta) { 
  
  part1 = delta*(1 - exp(-(t2-tm)))*((exp(t2)*(t2-tm))/(exp(t2)-exp(tm)))
  part2 = delta*((1/3)*exp(-(t2-tm)))
  part3 = (1-delta)*((1/3)*exp(-(t2-t1)))
  P2P3_cov = part1 + part2 + part3
  
  return(P2P3_cov)
  
}

calc_cov_P1P3 <- function(t2, t1, tm, delta) {
  
  part1 = (1-delta)*(1/3)*(exp(-(t2-t1)))
  part2 = delta*(1/3)*(exp(-(t2-tm)))
  P1P3_cov = part1 + part2
  
  return(P1P3_cov)
}


### Function to simulate expression data ### 

sim_expression_data <- function(t2, t1, delta, tm, sigma2,
                                n_genes) {
  
  #Calculate covariances 
  
  P1P2_covar <- calc_cov_P1P2(t2, t1, tm, delta)
  P2P3_covar <- calc_cov_P2P3(t2, t1, tm, delta)
  P1P3_covar <- calc_cov_P1P3(t2, t1, tm, delta)
  
  #Set up covariance matrix 
  
  cov_matrix <- c(t2, P1P2_covar, P1P3_covar, 
                  P1P2_covar, t2, P2P3_covar, 
                  P1P3_covar, P2P3_covar, t2)
  
  cov_matrix <- matrix((cov_matrix*sigma2), 3, 3)
  cov_matrix <- make.positive.definite(cov_matrix, tol=1e-3)
  
  #Simulate trait values 
  
  P1_traits <- vector()
  P2_traits <- vector()
  P3_traits <- vector() 
  
  for (i in 1:n_genes){
    
    trait_vector <- rtmvnorm(n = 1, c(0, 0, 0), cov_matrix)
    
    P1_traits <- append(P1_traits, trait_vector[1])
    P2_traits <- append(P2_traits, trait_vector[2])
    P3_traits <- append(P3_traits, trait_vector[3])
    
  }
  
  trait_matrix <- cbind("P1" = P1_traits, 
                        "P2" = P2_traits,
                        "P3" = P3_traits)
  
  return(trait_matrix)
  
}

### Functions to calculate asymmetry and pval from expression dataset ### 

calc_Q3 <- function(gene, P1, P2, P3) {
  
  trio <- as.numeric(gene[c(P1, P2, P3)])
  dP1P3 <- trio[1] - trio[3]
  dP2P3 <- trio[2] - trio[3]
  Q3 <- (abs(dP2P3) - abs(dP1P3))/(abs(dP2P3) + abs(dP1P3))
  
  return(Q3)
}

get_pval <- function(dataset) {
  
  Q3_vals <- apply(dataset, 1, calc_Q3, "P1", "P2", "P3")
  
  if (mean(Q3_vals) > 0) {
    expected_direction <- "No"
  } else {
    expected_direction <- "Yes"
  }
  
  
  pval <- t.test(x=Q3_vals, mu = 0)$p.value
  
  return(c(pval, expected_direction))
}

### Simulate datasets ### 

t2_vector <- rep(2, 27)
t1_vector <- rep(c(1, 1.5, 1.9), 9)
delta_vector <- rep(c(0.01, 0.05, 0.1), 
                    times = (c(9, 9, 9)))
tm_vector <- rep(c(0.5, 1, 1.4, 0.75, 1.25, 1.65, 0.9, 1.4, 1.8),
                 3)

results <- vector()

for (i in 1:length(t2_vector)){
  
  trait_matrix <- sim_expression_data(t2_vector[i], 
                                      t1_vector[i],
                                      delta_vector[i],
                                      tm_vector[i],
                                      sigma2, n_genes)
  
  result <- get_pval(trait_matrix)
  results <- append(results, result)
  
}

print(results)
 
