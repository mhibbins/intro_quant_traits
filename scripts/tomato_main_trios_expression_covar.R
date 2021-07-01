remove(list=ls())
library(tidyr)
library(dplyr)
library(ggplot2)
library(gridExtra)

#This script performs the analyses relating to Q3 in the tomato gene expression
#data, and plots the results. Generates the plots in Figure 4, panels B and C. 

                        ### Load expression dataset ### 

ovule_exp_cpm <- read.csv(
  "C:/Users/18126/OneDrive - Indiana University/Projects/intro_gene_expression/datasets/ovule_exp_cpm.csv")

                                ### Data tidying ### 

#Calculate means from individual replicates 

LA3475 <- rowMeans(ovule_exp_cpm[,c("LA101A", "LA101B","LA101.C.A", "LA101.C.B")])
LA1316 <- rowMeans(ovule_exp_cpm[,c("LA1316.4A", "LA1316.4B", "LA1316.9")])
LA1589 <- rowMeans(ovule_exp_cpm[,c("LA1589.1", "LA1589.7")])
LA1777 <- ovule_exp_cpm[,"LA1777.21"]
LA2172 <- rowMeans(ovule_exp_cpm[,c("LA2172.19", "LA2172.9")])
LA3778 <- rowMeans(ovule_exp_cpm[,c("LA3778.226", "LA3778.8")])
LA4117 <- rowMeans(ovule_exp_cpm[,c("LA4117A.33", "LA4117A.38")])
LA0716 <- rowMeans(ovule_exp_cpm[,c("LA716A.1", "LA716A.2", "LA716.B", "LA716.C_S2")])

#Make tidy dataframe

tidy_ovule_exp_cpm <- as.data.frame(cbind("Gene" = ovule_exp_cpm$Gene, LA3475, LA1316, 
                                          LA1589, LA1777, LA2172, LA3778, LA4117, LA0716))


           ### Calculate asymmetry stat across genes for both trios ### 

calc_CPM_asymmetry <- function(gene, P1, P2, P3) {
  
  trio <- as.numeric(gene[c(P1, P2, P3)])
  dP1P3 <- trio[1] - trio[3]
  dP2P3 <- trio[2] - trio[3]
  asym <- (abs(dP2P3) - abs(dP1P3))/(abs(dP2P3) + abs(dP1P3))
  
  return(asym)
}

intro_trio_asymmetry <- apply(tidy_ovule_exp_cpm, 1, calc_CPM_asymmetry,
                              "LA3778", "LA1777", "LA1316")
control_trio_asymmetry <- apply(tidy_ovule_exp_cpm, 1, calc_CPM_asymmetry,
                              "LA3475", "LA1589", "LA0716")

              ### Plot distributions of Q3 for both trios ### 

intro_trio_dist <- ggplot() + aes(intro_trio_asymmetry) + geom_density(size=1)
intro_trio_dist <- intro_trio_dist + theme_minimal(base_size = 20)
intro_trio_dist <- intro_trio_dist + geom_vline(xintercept = 0, 
                                                linetype="dashed",
                                                size=1)
intro_trio_dist <- intro_trio_dist + labs(x = "Q3")

cont_trio_dist <- ggplot() + aes(control_trio_asymmetry) + geom_density(size=1)
cont_trio_dist <- cont_trio_dist + theme_minimal(base_size = 20)
cont_trio_dist <- cont_trio_dist + geom_vline(xintercept = 0, 
                                                linetype="dashed",
                                                size=1)
cont_trio_dist <- cont_trio_dist + labs(x = "Q3")


              ### Bootstrap sign tests for significance ### 

intro_trio_successes = c(sum(intro_trio_asymmetry < 0), 
                         sum(intro_trio_asymmetry > 0)) #negatives, positives 
control_trio_successes = c(sum(control_trio_asymmetry < 0),
                           sum(control_trio_asymmetry > 0))

#Introgression trio 

intro_trio_sign_rank <- 0
boot_intro_diffs <- vector()

for (i in 1:10000) {
  
  intro_bootstrapped_signs <- sample(c(-1, 1), size = 14556, replace = TRUE)
  
  intro_bootstrap_successes <- c(sum(intro_bootstrapped_signs < 0), 
                           sum(intro_bootstrapped_signs > 0))
  
  diff <- intro_bootstrap_successes[1] - intro_bootstrap_successes[2]
  boot_intro_diffs <- append(boot_intro_diffs, diff)
  
  if ((intro_trio_successes[1]-intro_trio_successes[2]) < (
    (intro_bootstrap_successes[1]-intro_bootstrap_successes[2])
  )) {intro_trio_sign_rank <- intro_trio_sign_rank + 1}
  
}


#Control trio 

cont_trio_sign_rank <- 0
boot_cont_diffs <- vector()

for (i in 1:10000) {
  
  cont_bootstrapped_signs <- sample(c(-1, 1), size = 14556, replace = TRUE)
  
  cont_bootstrap_successes <- c(sum(cont_bootstrapped_signs < 0), 
                                 sum(cont_bootstrapped_signs > 0))
  
  diff <- cont_bootstrap_successes[1] - cont_bootstrap_successes[2]
  boot_cont_diffs <- append(boot_cont_diffs, diff)
  
  if ((control_trio_successes[2]-control_trio_successes[1]) < (
    (cont_bootstrap_successes[2]-cont_bootstrap_successes[1])
  )) {cont_trio_sign_rank <- cont_trio_sign_rank + 1}
  
}

intro_trio_sign_pval <- 1 - 2*abs(0.5-(intro_trio_sign_rank/10000))
cont_trio_sign_pval <- 1 - 2*abs(0.5-((cont_trio_sign_rank)/10000))

               ### Bootstrap mean tests for significance ### 

#Introgression trio 

intro_trio_mean_rank <- 0

for (i in 1:10000) {
  
  intro_bootstrapped_vals <- sample(intro_trio_asymmetry, size = 14556, replace = TRUE)
  
  if (mean(intro_bootstrapped_vals) >= 0)
   {intro_trio_mean_rank <- intro_trio_mean_rank + 1}
  
}


#Control trio 

cont_trio_mean_rank <- 0

for (i in 1:10000) {
  
  cont_bootstrapped_vals <- sample(control_trio_asymmetry, size = 14556, replace = TRUE)
  
  if (mean(cont_bootstrapped_vals) <= 0)
  {cont_trio_mean_rank <- cont_trio_mean_rank + 1}
  
  
}

intro_trio_mean_pval <- 1 - 2*abs(0.5-(intro_trio_mean_rank/10000))
cont_trio_mean_pval <- 1 - 2*abs(0.5-((cont_trio_mean_rank)/10000))

                       ### Plot mean asymmetries ### 

results_df <- as.data.frame(cbind("Introgression" = intro_trio_asymmetry, 
                                  "Control" = control_trio_asymmetry))

results_df_long <- gather(results_df, "Trio", "CPM_Asymmetry")
results_df_long$Trio <- factor(results_df_long$Trio, 
                               levels = c("Introgression", "Control"))

asymmetry_plot <- ggplot(results_df_long, aes(x = Trio,
                                              y = CPM_Asymmetry))
asymmetry_plot <- asymmetry_plot + stat_summary(fun.y = mean,
                                                geom = "point", 
                                                size = 4)
asymmetry_plot <- asymmetry_plot + stat_summary(fun.data = "mean_se", 
                                                geom = "errorbar",
                                                width=0.5,
                                                size = 1)
asymmetry_plot <- asymmetry_plot + geom_hline(aes(yintercept = 0),
                                              linetype = "dashed", 
                                              size = 1)
asymmetry_plot <- asymmetry_plot + theme_minimal(base_size = 30)
asymmetry_plot <- asymmetry_plot + labs(y = expression(paste(Q[3])))
asymmetry_plot <- asymmetry_plot + scale_x_discrete(
                                    labels = c("High", "Low"))
asymmetry_plot <- asymmetry_plot + coord_flip(ylim = c(-0.03, 0.03))

                        ### Plot sign counts ### 

obs_cont_diff <- control_trio_successes[1] - control_trio_successes[2]
obs_intro_diff <- intro_trio_successes[1] - intro_trio_successes[2]
boot_cont_diff <- mean(boot_cont_diffs)
boot_intro_diff <- mean(boot_intro_diffs)
boot_cont_se <- sd(boot_cont_diffs)#/sqrt(length(boot_cont_diffs))
boot_intro_se <- sd(boot_intro_diffs)#/sqrt(length(boot_intro_diffs))

lowtrio_counts_plot <- ggplot() + geom_density(aes(x=boot_cont_diffs),
                                               size = 2)
lowtrio_counts_plot <- lowtrio_counts_plot + scale_x_continuous(limits = c(-500, 500))
lowtrio_counts_plot <- lowtrio_counts_plot + theme_minimal(base_size = 20)
lowtrio_counts_plot <- lowtrio_counts_plot + theme(axis.title.x=element_blank(),
                                        plot.title = element_text(hjust = 0.5))
lowtrio_counts_plot <- lowtrio_counts_plot + labs(
  title = "Low trio")
lowtrio_counts_plot <- lowtrio_counts_plot + geom_vline(xintercept = obs_cont_diff,
                                                        linetype = "dashed", size = 1.5)

hightrio_counts_plot <- ggplot() + geom_density(aes(x=boot_intro_diffs),
                                               size = 2)
hightrio_counts_plot <- hightrio_counts_plot + theme_minimal(base_size = 20)
hightrio_counts_plot <- hightrio_counts_plot + labs(
  x = expression(paste("Negative ", Q[3], " values - positive ", 
                       Q[3], " values")),
  title = "High trio")
hightrio_counts_plot <- hightrio_counts_plot + theme(
  plot.title = element_text(hjust = 0.5))
hightrio_counts_plot <- hightrio_counts_plot + geom_vline(xintercept = obs_intro_diff,
                                                        linetype = "dashed", size = 1.5)

counts_plot <- grid.arrange(lowtrio_counts_plot, hightrio_counts_plot,
                            nrow = 2)

