#library(devtools)
#install_github("allenkei/CPDstergm")
library(CPDstergm) 
library(ggplot2)


# SBM
set.seed(1)
SBM_rho05 <- sim_SBM_list(num_seq=1, n=100, rho=0.5)

par(mfrow=c(2,3)) # Figures: 9 by 6 inches
for(i in c(25,50,75, 26,51,76)){
  par(mar = c(2.2, 2.2, 2.2, 2.2)); 
  image(SBM_rho05[[1]][[i]], xaxt = "n", yaxt = "n")
  axis(side=1,at=seq(0,1,length.out = 4),labels=round(seq(1,100,length.out = 4)),xpd=NA,cex.axis=1.4)
  axis(side=2,at=seq(0,1,length.out = 4),labels=round(seq(1,100,length.out = 4)),xpd=NA,cex.axis=1.4)
}

# STERGM
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
STERGM3 <- sim_STERGM_list(num_seq=1, n=100, network_stats, coefs_pos, coefs_neg, y1_stats=1400)

par(mfrow=c(2,3)) # Figures: 9 by 6 inches
for(i in c(25,50,75, 26,51,76)){
  par(mar = c(2.2, 2.2, 2.2, 2.2)); 
  image(STERGM3[[1]][[i]], xaxt = "n", yaxt = "n")
  axis(side=1,at=seq(0,1,length.out = 3),labels=round(seq(1,100,length.out = 3)),xpd=NA,cex.axis=1.4)
  axis(side=2,at=seq(0,1,length.out = 3),labels=round(seq(1,100,length.out = 3)),xpd=NA,cex.axis=1.4)  
}

# RDPGM
set.seed(1)
RDPG_d15 <- sim_RDPG_list(num_seq=1, n=100, rho=0.9, d=15)

par(mfrow=c(2,3)) # Figures: 9 by 6 inches
for(i in c(25,50,75, 26,51,76)){
  par(mar = c(2.2, 2.2, 2.2, 2.2)); 
  image(RDPG_d15[[1]][[i]], xaxt = "n", yaxt = "n")
  axis(side=1,at=seq(0,1,length.out = 4),labels=round(seq(1,100,length.out = 4)),xpd=NA,cex.axis=1.4)
  axis(side=2,at=seq(0,1,length.out = 4),labels=round(seq(1,100,length.out = 4)),xpd=NA,cex.axis=1.4)
}



#################
# MITphone Data #
#################

# 8 by 5
data(MITphone)
network_list <- MITphone

n_time <- length(network_list) 
n_nodes <- nrow(network_list[[1]])

total_edges <- n_time * n_nodes * n_nodes
i_vec <- integer(total_edges)
j_vec <- integer(total_edges)
time_vec <- integer(total_edges)
value_vec <- integer(total_edges)

counter <- 1

for (t in seq_along(network_list)) {
  mat <- network_list[[t]]
  for (i in 1:n_nodes) {
    for (j in 1:n_nodes) {
      i_vec[counter] <- i
      j_vec[counter] <- j
      time_vec[counter] <- t
      value_vec[counter] <- mat[i, j]
      counter <- counter + 1
    }
  }
}


raster_data <- data.frame(i = i_vec, j = j_vec,
                          time = time_vec, value = value_vec)


raster_data <- raster_data[raster_data$i < raster_data$j, ]
raster_data$edge_id <- paste(raster_data$i, raster_data$j, sep = "-")
change_points <- c(29,40,63,94,191) 

library(ggplot2)


highlight_boxes <- data.frame(
  xmin = c(95,   1    , 192),
  xmax = c(115,  93   , 232),
  ymin = c(1,    500  , 1000),
  ymax = c(2459, 1600 , 2200)  # 2459 is the number of active edges
)

ggplot(raster_data[raster_data$value == 1, ], aes(x = time, y = as.numeric(factor(edge_id)))) +
  geom_tile(fill = "black") +
  geom_rect(data = highlight_boxes,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "blue", alpha = 0.15, color = "blue", inherit.aes = FALSE) +
  annotate("text",
           x = c(99,   5,    197),      
           y = c(2350, 1500, 2100),      
           label = c("I", "III","II"),
           size = 5, fontface = "bold", color = "blue") +
  
  labs(x = "Time points", y = "Active edges") +
  theme_minimal() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  geom_vline(xintercept = change_points, color = "red", linewidth = 0.5)







#####################
# Stock Market Data #
#####################
# 8 by 5

library(ecp)
data(DJIA)
market <- DJIA$market
date_vec <- DJIA$dates[1:1138]
rownames(market) <- date_vec

start <- which(date_vec == '2007-01-01')
end <- which(date_vec == '2010-01-04')
date_range <- start:end
mydata <- array(NA,c(29, 29, length(date_range)))
df <- list()
for(i in 1:length(date_range)){
  temp <- market[date_range[i]:(date_range[i]+3),]
  temp <- ifelse(cor(temp)< 0, 1, 0)
  diag(temp) <- 0
  df[[i]] <- temp
}


network_list <- df; rm(df, temp, start, end, mydata, market, DJIA)
n_time <- length(network_list) 
n_nodes <- nrow(network_list[[1]])

total_edges <- n_time * n_nodes * n_nodes
i_vec <- integer(total_edges)
j_vec <- integer(total_edges)
time_vec <- integer(total_edges)
value_vec <- integer(total_edges)

counter <- 1

for (t in seq_along(network_list)) {
  mat <- network_list[[t]]
  for (i in 1:n_nodes) {
    for (j in 1:n_nodes) {
      i_vec[counter] <- i
      j_vec[counter] <- j
      time_vec[counter] <- t
      value_vec[counter] <- mat[i, j]
      counter <- counter + 1
    }
  }
}


raster_data <- data.frame(i = i_vec, j = j_vec,
                          time = time_vec, value = value_vec)

raster_data <- raster_data[raster_data$i < raster_data$j, ]
raster_data$edge_id <- paste(raster_data$i, raster_data$j, sep = "-")
change_points <- c(17, 93, 121) 

highlight_boxes <- data.frame(
  xmin = c(18,  94,  122),
  xmax = c(30,  106, 134),
  ymin = c(1,   1,   1),
  ymax = c(440, 440, 440)  
)

ggplot(raster_data[raster_data$value == 1, ], aes(x = time, y = as.numeric(factor(edge_id)))) +
  geom_tile(fill = "black") +
  geom_rect(data = highlight_boxes,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "blue", alpha = 0.15, color = "blue", inherit.aes = FALSE) +
  annotate("text",
           x = c(20,  97, 126),      
           y = c(420, 420, 420),      
           label = c("I", "II", "III"),
           size = 5, fontface = "bold", color = "blue") +
  
  labs(x = "Time points", y = "Active edges") +
  theme_minimal() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  geom_vline(xintercept = change_points, color = "red", linewidth = 0.5)



###################
# Time Comparison #
###################

library(CPDstergm)
source("EVAL.R")

num_node <- c(50, 100)
network_stats=c("edges", "mutual")
SBM_time_by_node <- matrix(NA, nrow=7, ncol=3)

for(i in 1:length(num_node)){
  set.seed(1)
  SBM_list <- sim_SBM_list(num_seq=3, n=num_node[i], rho=0.0)
  
  time1 <- system.time( CPD_STERGM_list(SBM_list, directed=TRUE, network_stats, list_of_lambda=100) )
  time2 <- system.time( Evaluation_gSeg(SBM_list, p_threshold=0.05) )
  time3 <- system.time( Evaluation_kerSeg(SBM_list, p_threshold=0.001) )
  time4 <- system.time( Evaluation_gSeg_on_stats(SBM_list, p_threshold=0.05, num_stats=length(network_stats)) )
  time5 <- system.time( Evaluation_kerSeg_on_stats(SBM_list, p_threshold=0.001, num_stats=length(network_stats)) )
  time6 <- system.time( Evaluation_RDPG(SBM_list, M=50, d=5, delta=5) )
  time7 <- system.time( Evaluation_NBS(SBM_list, M=15, delta=5) )
  time8 <- system.time( Evaluation_CPDker(SBM_list, M=30, num_stats=length(network_stats)) )
  
  SBM_time_by_node[,i] <- c(time1['elapsed'],time2['elapsed'],time3['elapsed'],time4['elapsed'],
                            time5['elapsed'],time6['elapsed'],time7['elapsed'],time8['elapsed'])
  
}

#save(SBM_time_by_node, file = 'time_SBMrho0.Rdata')





num_node <- c(50, 100)
y1_target <- c(250, 500)
network_stats=c("edges", "mutual")
coefs_pos <- matrix(c(-1, -1, -1, -1, -2, 1, -2, 1), nrow=2, ncol=4, byrow = T)
coefs_neg <- matrix(c( -1, -1, -1, -1, -2, -1, -2, -1), nrow=2, ncol=4, byrow = T)


STERGM_time_by_node <- matrix(NA, nrow=7, ncol=3)

for(i in 1:length(num_node)){
  set.seed(1)
  STERGM_list <- sim_STERGM_list(num_seq=3, n=num_node[i], network_stats, coefs_pos, coefs_neg, y1_stats=y1_target[i])
  
  time1 <- system.time( CPD_STERGM_list(STERGM_list, directed=TRUE, network_stats, list_of_lambda=100) )
  time2 <- system.time( Evaluation_gSeg(STERGM_list, p_threshold=0.05) )
  time3 <- system.time( Evaluation_kerSeg(STERGM_list, p_threshold=0.001) )
  time4 <- system.time( Evaluation_gSeg_on_stats(STERGM_list, p_threshold=0.05, num_stats=length(network_stats)) )
  time5 <- system.time( Evaluation_kerSeg_on_stats(STERGM_list, p_threshold=0.001, num_stats=length(network_stats)) )
  time6 <- system.time( Evaluation_RDPG(STERGM_list, M=50, d=5, delta=5) )
  time7 <- system.time( Evaluation_NBS(STERGM_list, M=15, delta=5) )
  time8 <- system.time( Evaluation_CPDker(STERGM_list, M=30, num_stats=length(network_stats)) )
  
  STERGM_time_by_node[,i] <- c(time1['elapsed'],time2['elapsed'],time3['elapsed'],time4['elapsed'],
                               time5['elapsed'],time6['elapsed'],time7['elapsed'],time8['elapsed'])
  
}

#save(STERGM_time_by_node, file = 'time_STERGMp4.Rdata')




num_node <- c(50, 100)
network_stats=c("edges", "mutual")
RDPGM_time_by_node <- matrix(NA, nrow=7, ncol=3)

for(i in 1:length(num_node)){
  set.seed(1)
  RDPG_list <- sim_RDPG_list(num_seq = 3, n = num_node[i], rho = 0.9, d = 20)
  
  time1 <- system.time( CPD_STERGM_list(RDPG_list, directed=TRUE, network_stats, list_of_lambda=10) )
  time2 <- system.time( Evaluation_gSeg(RDPG_list, p_threshold=0.05) )
  time3 <- system.time( Evaluation_kerSeg(RDPG_list, p_threshold=0.001) )
  time4 <- system.time( Evaluation_gSeg_on_stats(RDPG_list, p_threshold=0.05, num_stats=length(network_stats)) )
  time5 <- system.time( Evaluation_kerSeg_on_stats(RDPG_list, p_threshold=0.001, num_stats=length(network_stats)) )
  time6 <- system.time( Evaluation_RDPG(RDPG_list, M=50, d=5, delta=5) )
  time7 <- system.time( Evaluation_NBS(RDPG_list, M=15, delta=5) )
  time8 <- system.time( Evaluation_CPDker(RDPG_list, M=30, num_stats=length(network_stats)) )
  
  RDPGM_time_by_node[,i] <- c(time1['elapsed'],time2['elapsed'],time3['elapsed'],time4['elapsed'],
                              time5['elapsed'],time6['elapsed'],time7['elapsed'],time8['elapsed'])
  
}

#save(RDPGM_time_by_node, file = 'time_RDPGMd20.Rdata')


