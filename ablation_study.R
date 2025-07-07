
library(devtools)
#install_github("allenkei/CPDstergm", force = TRUE)
library(CPDstergm)



######################
# Simulation 1 (SBM) #
######################


network_stats=c("edges", "mutual")


set.seed(1)
SBM_list <- sim_SBM_list(num_seq=15, n=100, rho=0.0)

sim_result1 <- CPD_STERGM_list(SBM_list, directed=TRUE, network_stats, d_weight=FALSE, list_of_lambda=10^c(0:4))
sim_result2 <- CPD_STERGM_list(SBM_list, directed=TRUE, network_stats, update_alpha=FALSE, list_of_lambda=10^c(0:4))

final_results <- list(sim_result1, sim_result2)
save(final_results, file = 'SBMrho00_ablation.Rdata')


#########################
# Simulation 2 (STERGM) #
#########################


y1_target <- 1400
network_stats=c("edges", "mutual", "triangles")
coefs_pos <- matrix(c(-2, -1.5, -2, -1.5,
                      2,  1,    2,  1,
                      -2, -1,   -2, -1),
                    nrow=3, ncol=4, byrow = T)

coefs_neg <- matrix(c( -1,  2,    -1,  2,
                       2,  1,     2,  1,
                       1,  1.5,   1,  1.5),
                    nrow=3, ncol=4, byrow = T)



set.seed(1)
STERGM_list <- sim_STERGM_list(num_seq=15, n=100, network_stats, coefs_pos, coefs_neg, y1_stats=y1_target)

sim_result1 <- CPD_STERGM_list(STERGM_list, directed=TRUE, network_stats, d_weight=FALSE, list_of_lambda=10^c(0:4))
sim_result2 <- CPD_STERGM_list(STERGM_list, directed=TRUE, network_stats, update_alpha=FALSE, list_of_lambda=10^c(0:4))

final_results <- list(sim_result1, sim_result2)
save(final_results, file = 'STERGMp6_ablation.Rdata')



########################
# Simulation 3 (RDPGM) #
########################



network_stats=c("edges", "mutual")


set.seed(1)
RDPG_list <- sim_RDPG_list(num_seq=15, n=100, rho=0.9, d=10)

sim_result1 <- CPD_STERGM_list(RDPG_list, directed=TRUE, network_stats, d_weight=FALSE, list_of_lambda=10^c(0:4))
sim_result2 <- CPD_STERGM_list(RDPG_list, directed=TRUE, network_stats, update_alpha=FALSE, list_of_lambda=10^c(0:4))

final_results <- list(sim_result1, sim_result2)
save(final_results, file = 'RDPGd10_ablation.Rdata')


###############
# Print Table #
###############
load(file.choose()) #final_results
iter <- 2

means <- round(apply(final_results[[iter]], 2, mean), 2)
sds   <- round(apply(final_results[[iter]], 2, sd), 2)

formatted <- mapply(function(m, s) sprintf("$%.1f$ $(%.1f)$", m, s), means, sds)
formatted[4] <- sprintf("$%.2f\\%%$", mean(final_results[[iter]][, 4]) * 100)

latex_row <- paste0(paste(formatted, collapse = " & "), " \\\\")

cat(latex_row)
