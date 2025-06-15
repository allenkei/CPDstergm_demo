
library(devtools)
#install_github("allenkei/CPDstergm", force = TRUE)
library(CPDstergm)
source("EVAL.R")




######################
# Simulation 1 (SBM) #
######################


num_seq = 15
num_node <- c(50, 100, 200)
network_stats=c("edges", "mutual")


##################
# SBM, rho = 0.0 #
##################

final_results <- list()

for(i in 1:length(num_node)){
  set.seed(1)
  SBM_list <- sim_SBM_list(num_seq, n=num_node[i], rho=0.0)

  sim_result1 <- CPD_STERGM_list(SBM_list, directed=TRUE, network_stats, list_of_lambda=10^c(0:4))
  sim_result2 <- Evaluation_gSeg(SBM_list, p_threshold=0.05)
  sim_result3 <- Evaluation_kerSeg(SBM_list, p_threshold=0.001)
  sim_result4 <- Evaluation_gSeg_on_stats(SBM_list, p_threshold=0.05, num_stats=length(network_stats))
  sim_result5 <- Evaluation_kerSeg_on_stats(SBM_list, p_threshold=0.001, num_stats=length(network_stats))
  sim_result6 <- Evaluation_RDPG(SBM_list, M=50, d=5, delta=5)
  sim_result7 <- Evaluation_NBS(SBM_list, M=15, delta=5)

  final_results <- c(final_results, list(
    sim_result1, sim_result2, sim_result3, sim_result4, sim_result5, sim_result6, sim_result7
  ))
  
}


#save(final_results, file = 'SBMrho00.Rdata')


##################
# SBM, rho = 0.5 #
##################


final_results <- list()


for(i in 1:length(num_node)){
  set.seed(1)
  SBM_list <- sim_SBM_list(num_seq, n=num_node[i], rho=0.5)

  sim_result1 <- CPD_STERGM_list(SBM_list, directed=TRUE, network_stats, list_of_lambda=10^c(0:4))
  sim_result2 <- Evaluation_gSeg(SBM_list, p_threshold=0.05)
  sim_result3 <- Evaluation_kerSeg(SBM_list, p_threshold=0.001)
  sim_result4 <- Evaluation_gSeg_on_stats(SBM_list, p_threshold=0.05, num_stats=length(network_stats))
  sim_result5 <- Evaluation_kerSeg_on_stats(SBM_list, p_threshold=0.001, num_stats=length(network_stats))
  sim_result6 <- Evaluation_RDPG(SBM_list, M=50, d=5, delta=5)
  sim_result7 <- Evaluation_NBS(SBM_list, M=15, delta=5)

  final_results <- c(final_results, list(
    sim_result1, sim_result2, sim_result3, sim_result4, sim_result5, sim_result6, sim_result7
  ))
  
}


#save(final_results, file = 'SBMrho05.Rdata')


##################
# SBM, rho = 0.9 #
##################

final_results <- list()

for(i in 1:length(num_node)){
  set.seed(1)
  SBM_list <- sim_SBM_list(num_seq, n=num_node[i], rho=0.9)

  sim_result1 <- CPD_STERGM_list(SBM_list, directed=TRUE, network_stats, list_of_lambda=10^c(0:4))
  sim_result2 <- Evaluation_gSeg(SBM_list, p_threshold=0.05)
  sim_result3 <- Evaluation_kerSeg(SBM_list, p_threshold=0.001)
  sim_result4 <- Evaluation_gSeg_on_stats(SBM_list, p_threshold=0.05, num_stats=length(network_stats))
  sim_result5 <- Evaluation_kerSeg_on_stats(SBM_list, p_threshold=0.001, num_stats=length(network_stats))
  sim_result6 <- Evaluation_RDPG(SBM_list, M=50, d=5, delta=5)
  sim_result7 <- Evaluation_NBS(SBM_list, M=15, delta=5)
  
  final_results <- c(final_results, list(
    sim_result1, sim_result2, sim_result3, sim_result4, sim_result5, sim_result6, sim_result7
  ))
  
}


#save(final_results, file = 'SBMrho09.Rdata')


#########################
# Simulation 2 (STERGM) #
#########################


#################
# STERGM, p = 4 #
#################


num_seq = 15
num_node <- c(50, 100, 200)
y1_target <- c(250, 500, 1000)
network_stats=c("edges", "mutual")
coefs_pos <- matrix(c(-1, -1, -1, -1, -2, 1, -2, 1), nrow=2, ncol=4, byrow = T)
coefs_neg <- matrix(c( -1, -1, -1, -1, -2, -1, -2, -1), nrow=2, ncol=4, byrow = T)

final_results <- list()

for(i in 1:length(num_node)){
  set.seed(1)
  STERGM_list <- sim_STERGM_list(num_seq, n=num_node[i], network_stats, coefs_pos, coefs_neg, y1_stats=y1_target[i])

  sim_result1 <- CPD_STERGM_list(STERGM_list, directed=TRUE, network_stats, list_of_lambda=10^c(0:4))
  sim_result2 <- Evaluation_gSeg(STERGM_list, p_threshold=0.05)
  sim_result3 <- Evaluation_kerSeg(STERGM_list, p_threshold=0.001)
  sim_result4 <- Evaluation_gSeg_on_stats(STERGM_list, p_threshold=0.05, num_stats=length(network_stats))
  sim_result5 <- Evaluation_kerSeg_on_stats(STERGM_list, p_threshold=0.001, num_stats=length(network_stats))
  sim_result6 <- Evaluation_RDPG(STERGM_list, M=50, d=5, delta=5)
  sim_result7 <- Evaluation_NBS(STERGM_list, M=15, delta=5)
  
  final_results <- c(final_results, list(
    sim_result1, sim_result2, sim_result3, sim_result4, sim_result5, sim_result6, sim_result7
  ))
  
}


#save(final_results, file = 'STERGMp4.Rdata')


#################
# STERGM, p = 6 #
#################


num_seq = 15
num_node <- c(50, 100, 200)
y1_target <- c(500, 1400, 2500)
network_stats=c("edges", "mutual", "triangles")
coefs_pos <- matrix(c(-2, -1.5, -2, -1.5,
                       2,  1,    2,  1,
                      -2, -1,   -2, -1),
                    nrow=3, ncol=4, byrow = T)

coefs_neg <- matrix(c( -1,  2,    -1,  2,
                        2,  1,     2,  1,
                        1,  1.5,   1,  1.5),
                    nrow=3, ncol=4, byrow = T)


final_results <- list()

for(i in 1:length(num_node)){
  set.seed(1)
  STERGM_list <- sim_STERGM_list(num_seq, n=num_node[i], network_stats, coefs_pos, coefs_neg, y1_stats=y1_target[i])

  sim_result1 <- CPD_STERGM_list(STERGM_list, directed=TRUE, network_stats, list_of_lambda=10^c(0:4))
  sim_result2 <- Evaluation_gSeg(STERGM_list, p_threshold=0.05)
  sim_result3 <- Evaluation_kerSeg(STERGM_list, p_threshold=0.001)
  sim_result4 <- Evaluation_gSeg_on_stats(STERGM_list, p_threshold=0.05, num_stats=length(network_stats))
  sim_result5 <- Evaluation_kerSeg_on_stats(STERGM_list, p_threshold=0.001, num_stats=length(network_stats))
  sim_result6 <- Evaluation_RDPG(STERGM_list, M=50, d=5, delta=5)
  sim_result7 <- Evaluation_NBS(STERGM_list, M=15, delta=5)

  
  final_results <- c(final_results, list(
    sim_result1, sim_result2, sim_result3, sim_result4, sim_result5, sim_result6, sim_result7
  ))
  
}


#save(final_results, file = 'STERGMp6.Rdata')


#################
# STERGM, p = 8 #
#################


num_seq = 15
num_node <- c(50, 100, 200)
y1_target <- c(500, 1400, 2500)
network_stats <- c("edges", "mutual", "triangles", "nodematch(\"node_attr\")")
coefs_pos <- matrix(c(-2, -1.5, -2, -1.5,
                       2,  1,    2,  1,
                      -2, -1,   -2, -1,
                      -1,  1,   -1,  1),
                    nrow=4, ncol=4, byrow = T)

coefs_neg <- matrix(c( -1,  2,    -1,  2,
                        2,  1,     2,  1,
                        1,  1.5,   1,  1.5,
                        1,  2,     1,  2),
                    nrow=4, ncol=4, byrow = T)


final_results <- list()

for(i in 1:length(num_node)){
  set.seed(1)
  gender <- c(c("M", "F")[rbinom(num_node[i], 1, 0.5) + 1]) # fixed
  STERGM_list <- sim_STERGM_list(num_seq, n=num_node[i], network_stats, coefs_pos, coefs_neg, y1_stats=y1_target[i], node_attr=gender)

  sim_result1 <- CPD_STERGM_list(STERGM_list, directed=TRUE, network_stats, node_attr=gender, list_of_lambda=10^c(0:4))
  sim_result2 <- Evaluation_gSeg(STERGM_list, p_threshold=0.05)
  sim_result3 <- Evaluation_kerSeg(STERGM_list, p_threshold=0.001)
  sim_result4 <- Evaluation_gSeg_on_stats(STERGM_list, p_threshold=0.05, num_stats=length(network_stats))
  sim_result5 <- Evaluation_kerSeg_on_stats(STERGM_list, p_threshold=0.001, num_stats=length(network_stats))
  sim_result6 <- Evaluation_RDPG(STERGM_list, M=50, d=5, delta=5)
  sim_result7 <- Evaluation_NBS(STERGM_list, M=15, delta=5)

  
  final_results <- c(final_results, list(
    sim_result1, sim_result2, sim_result3, sim_result4, sim_result5, sim_result6, sim_result7
  ))
  
}


#save(final_results, file = 'STERGMp8.Rdata')


########################
# Simulation 3 (RDPGM) #
########################


source("EVAL.R")
num_seq = 15
num_node <- c(50, 100, 200)
network_stats=c("edges", "mutual")


#################
# RDPGM, d = 10 #
#################


final_results <- list()

for(i in 1:length(num_node)){
  set.seed(1)
  RDPG_list <- sim_RDPG_list(num_seq, n = num_node[i], rho = 0.9, d = 10)
  
  sim_result1 <- CPD_STERGM_list(RDPG_list, directed=TRUE, network_stats, list_of_lambda=10^c(0:4))
  sim_result2 <- Evaluation_gSeg(RDPG_list, p_threshold=0.05)
  sim_result3 <- Evaluation_kerSeg(RDPG_list, p_threshold=0.001)
  sim_result4 <- Evaluation_gSeg_on_stats(RDPG_list, p_threshold=0.05, num_stats=length(network_stats))
  sim_result5 <- Evaluation_kerSeg_on_stats(RDPG_list, p_threshold=0.001, num_stats=length(network_stats))
  sim_result6 <- Evaluation_RDPG(RDPG_list, M=50, d=5, delta=5)
  sim_result7 <- Evaluation_NBS(RDPG_list, M=15, delta=5)

  
  final_results <- c(final_results, list(
    sim_result1, sim_result2, sim_result3, sim_result4, sim_result5, sim_result6, sim_result7
  ))
  
}


#save(final_results, file = 'RDPGd10.Rdata')


#################
# RDPGM, d = 15 #
#################

final_results <- list()

for(i in 1:length(num_node)){
  set.seed(1)
  RDPG_list <- sim_RDPG_list(num_seq, n = num_node[i], rho = 0.9, d = 15)
  
  sim_result1 <- CPD_STERGM_list(RDPG_list, directed=TRUE, network_stats, list_of_lambda=10^c(0:4))
  sim_result2 <- Evaluation_gSeg(RDPG_list, p_threshold=0.05)
  sim_result3 <- Evaluation_kerSeg(RDPG_list, p_threshold=0.001)
  sim_result4 <- Evaluation_gSeg_on_stats(RDPG_list, p_threshold=0.05, num_stats=length(network_stats))
  sim_result5 <- Evaluation_kerSeg_on_stats(RDPG_list, p_threshold=0.001, num_stats=length(network_stats))
  sim_result6 <- Evaluation_RDPG(RDPG_list, M=50, d=5, delta=5)
  sim_result7 <- Evaluation_NBS(RDPG_list, M=15, delta=5)

  
  final_results <- c(final_results, list(
    sim_result1, sim_result2, sim_result3, sim_result4, sim_result5, sim_result6, sim_result7
  ))
  
}


#save(final_results, file = 'RDPGd15.Rdata')


#################
# RDPGM, d = 20 #
#################


final_results <- list()

for(i in 1:length(num_node)){
  set.seed(1)
  RDPG_list <- sim_RDPG_list(num_seq, n = num_node[i], rho = 0.9, d = 20)
  
  sim_result1 <- CPD_STERGM_list(RDPG_list, directed=TRUE, network_stats, list_of_lambda=10^c(0:4))
  sim_result2 <- Evaluation_gSeg(RDPG_list, p_threshold=0.05)
  sim_result3 <- Evaluation_kerSeg(RDPG_list, p_threshold=0.001)
  sim_result4 <- Evaluation_gSeg_on_stats(RDPG_list, p_threshold=0.05, num_stats=length(network_stats))
  sim_result5 <- Evaluation_kerSeg_on_stats(RDPG_list, p_threshold=0.001, num_stats=length(network_stats))
  sim_result6 <- Evaluation_RDPG(RDPG_list, M=50, d=5, delta=5)
  sim_result7 <- Evaluation_NBS(RDPG_list, M=15, delta=5)

  
  final_results <- c(final_results, list(
    sim_result1, sim_result2, sim_result3, sim_result4, sim_result5, sim_result6, sim_result7
  ))
  
}


#save(final_results, file = 'RDPGd20.Rdata')

