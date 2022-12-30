
library(devtools)
install_github("allenkei/CPDstergm")
library(CPDstergm)

source("EVAL.R")

################
# Simulation 1 #
################

num_node <- c(50, 100, 500)
network_stats=c("edges", "mutual")


result1 <- g_result1 <- ker_result1 <- g_stats_result1 <- ker_stats_result1 <- matrix(NA, nrow=3, ncol=4)
for(i in 1:length(num_node)){
  set.seed(1)
  SBM_list <- sim_SBM_list(num_seq=10, n=num_node[i], rho=0.0)
  
  sim_result <- CPD_STERGM_list(SBM_list, directed=TRUE, network_stats)
  result1[i,] <- colMeans(sim_result)
  
  sim_result <- Evaluation_gSeg(SBM_list, p_threshold=0.05)
  g_result1[i,] <- colMeans(sim_result)
  
  sim_result <- Evaluation_kerSeg(SBM_list, p_threshold=0.001)
  ker_result1[i,] <- colMeans(sim_result)
  
  sim_result <- Evaluation_gSeg_on_stats(SBM_list, p_threshold=0.05, num_stats=length(network_stats))
  g_stats_result1[i,] <- colMeans(sim_result)
  
  sim_result <- Evaluation_kerSeg_on_stats(SBM_list, p_threshold=0.001, num_stats=length(network_stats))
  ker_stats_result1[i,] <- colMeans(sim_result)
}



result2 <- g_result2 <- ker_result2 <- g_stats_result2 <- ker_stats_result2 <- matrix(NA, nrow=3, ncol=4)
for(i in 1:length(num_node)){
  set.seed(1)
  SBM_list <- sim_SBM_list(num_seq=10, n=num_node[i], rho=0.5)
  
  sim_result <- CPD_STERGM_list(SBM_list, directed=TRUE, network_stats)
  result2[i,] <- colMeans(sim_result)
  
  sim_result <- Evaluation_gSeg(SBM_list, p_threshold=0.05)
  g_result2[i,] <- colMeans(sim_result)
  
  sim_result <- Evaluation_kerSeg(SBM_list, p_threshold=0.001)
  ker_result2[i,] <- colMeans(sim_result)
  
  sim_result <- Evaluation_gSeg_on_stats(SBM_list, p_threshold=0.05, num_stats=length(network_stats))
  g_stats_result2[i,] <- colMeans(sim_result)
  
  sim_result <- Evaluation_kerSeg_on_stats(SBM_list, p_threshold=0.001, num_stats=length(network_stats))
  ker_stats_result2[i,] <- colMeans(sim_result)
}



result3 <- g_result3 <- ker_result3 <- g_stats_result3 <- ker_stats_result3 <- matrix(NA, nrow=3, ncol=4)
for(i in 1:length(num_node)){
  set.seed(1)
  SBM_list <- sim_SBM_list(num_seq=10, n=num_node[i], rho=0.9)
  
  sim_result <- CPD_STERGM_list(SBM_list, directed=TRUE, network_stats)
  result3[i,] <- colMeans(sim_result)
  
  sim_result <- Evaluation_gSeg(SBM_list, p_threshold=0.05)
  g_result3[i,] <- colMeans(sim_result)
  
  sim_result <- Evaluation_kerSeg(SBM_list, p_threshold=0.001)
  ker_result3[i,] <- colMeans(sim_result)
  
  sim_result <- Evaluation_gSeg_on_stats(SBM_list, p_threshold=0.05, num_stats=length(network_stats))
  g_stats_result3[i,] <- colMeans(sim_result)
  
  sim_result <- Evaluation_kerSeg_on_stats(SBM_list, p_threshold=0.001, num_stats=length(network_stats))
  ker_stats_result3[i,] <- colMeans(sim_result)
}



################
# Simulation 2 #
################

num_node <- c(50,100,500)
y1_target <- c(250,500,2500)
network_stats=c("edges", "mutual")
coefs_pos <- matrix(c(-1, -1, -1, -1, -2, 1, -2, 1), nrow=2, ncol=4, byrow = T)
coefs_neg <- matrix(c( -1, -1, -1, -1, -2, -1, -2, -1), nrow=2, ncol=4, byrow = T)


result1 <- g_result1 <- ker_result1 <- g_stats_result1 <- ker_stats_result1 <- matrix(NA, nrow=3, ncol=4)
for(i in 1:length(num_node)){
  set.seed(1)
  STERGM_list <- sim_STERGM_list(num_seq=10, n=num_node[i], network_stats, coefs_pos, coefs_neg, y1_stats=y1_target[i])
  
  sim_result <- CPD_STERGM_list(STERGM_list, directed=TRUE, network_stats)
  result1[i,] <- colMeans(sim_result)
  
  sim_result <- Evaluation_gSeg(STERGM_list, p_threshold=0.05)
  g_result1[i,] <- colMeans(sim_result)
  
  sim_result <- Evaluation_kerSeg(STERGM_list, p_threshold=0.001)
  ker_result1[i,] <- colMeans(sim_result)
  
  sim_result <- Evaluation_gSeg_on_stats(STERGM_list, p_threshold=0.05, num_stats=length(network_stats))
  g_stats_result1[i,] <- colMeans(sim_result)
  
  sim_result <- Evaluation_kerSeg_on_stats(STERGM_list, p_threshold=0.001, num_stats=length(network_stats))
  ker_stats_result1[i,] <- colMeans(sim_result)
}




num_node <- c(50,100,500)
y1_target <- c(500,1400,2500)
network_stats=c("edges", "mutual", "triangles")
coefs_pos <- matrix(c(-2, -1.5, -2, -1.5,
                       2,  1,    2,  1,
                      -2, -1,   -2, -1),
                    nrow=3, ncol=4, byrow = T)

coefs_neg <- matrix(c( -1,  2,    -1,  2,
                        2,  1,     2,  1,
                        1,  1.5,   1,  1.5),
                    nrow=3, ncol=4, byrow = T)



result2 <- g_result2 <- ker_result2 <- g_stats_result2 <- ker_stats_result2 <- matrix(NA, nrow=3, ncol=4)
for(i in 1:length(num_node)){
  set.seed(1)
  STERGM_list <- sim_STERGM_list(num_seq=10, n=num_node[i], network_stats, coefs_pos, coefs_neg, y1_stats=y1_target[i])
  
  sim_result <- CPD_STERGM_list(STERGM_list, directed=TRUE, network_stats)
  result2[i,] <- colMeans(sim_result)
  
  sim_result <- Evaluation_gSeg(STERGM_list, p_threshold=0.05)
  g_result2[i,] <- colMeans(sim_result)
  
  sim_result <- Evaluation_kerSeg(STERGM_list, p_threshold=0.001)
  ker_result2[i,] <- colMeans(sim_result)
  
  sim_result <- Evaluation_gSeg_on_stats(STERGM_list, p_threshold=0.05, num_stats=length(network_stats))
  g_stats_result2[i,] <- colMeans(sim_result)
  
  sim_result <- Evaluation_kerSeg_on_stats(STERGM_list, p_threshold=0.001, num_stats=length(network_stats))
  ker_stats_result2[i,] <- colMeans(sim_result)
}



num_node <- c(50,100,500)
y1_target <- c(500,1400,2500)
network_stats <- c("edges", "mutual", "triangles", "nodematch(\"Gender\")")
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



result3 <- g_result3 <- ker_result3 <- g_stats_result3 <- ker_stats_result3 <- matrix(NA, nrow=3, ncol=4)
for(i in 1:length(num_node)){
  set.seed(1)
  gender <- c(c("M", "F")[rbinom(num_node[i], 1, 0.5) + 1]) # fixed
  STERGM_list <- sim_STERGM_list(num_seq=10, n=num_node[i], network_stats, coefs_pos, coefs_neg, y1_stats=y1_target[i], node_attr=gender)
  
  sim_result <- CPD_STERGM_list(STERGM_list, directed=TRUE, network_stats, node_attr=gender)
  result3[i,] <- colMeans(sim_result)
  
  sim_result <- Evaluation_gSeg(STERGM_list, p_threshold=0.05)
  g_result3[i,] <- colMeans(sim_result)
  
  sim_result <- Evaluation_kerSeg(STERGM_list, p_threshold=0.001)
  ker_result3[i,] <- colMeans(sim_result)
  
  sim_result <- Evaluation_gSeg_on_stats(STERGM_list, p_threshold=0.05, num_stats=length(network_stats))
  g_stats_result3[i,] <- colMeans(sim_result)
  
  sim_result <- Evaluation_kerSeg_on_stats(STERGM_list, p_threshold=0.001, num_stats=length(network_stats))
  ker_stats_result3[i,] <- colMeans(sim_result)
}



######################
# MIT Cellphone Data #
######################

data("MITphone")
result <- CPD_STERGM(MITphone, directed=FALSE, network_stats=c("edges", "isolates", "triangles"))


gSeg_result <- Evaluation_gSeg_on_stats(MITphone, p_threshold=0.05, num_stats=3, is_experiment=TRUE)
kerSeg_result <- Evaluation_kerSeg_on_stats(MITphone, p_threshold=0.001, num_stats=3, is_experiment=TRUE)


theta_change <- result$theta_change; threshold <- result$threshold; xtick <- result$est_CP
seq_date <- seq(as.Date("2004-09-15"), as.Date("2005-05-04"), by="days"); tau <- length(seq_date)-1


par(mar=c(4, 4, 2, 1), fig=c(0,1,0,0.74))
plot(1:length(theta_change), theta_change, type='l',ylab="", xlab="", xaxt="n", yaxt="n")
abline(h = threshold, col='red',lwd=2)

# winter and spring vacation (duration): seq_date[95];seq_date[110];seq_date[188];seq_date[192]
win_spr_break <- c(95,110,188,192)

for(i in c(1,3)){
  rect(xleft = win_spr_break[i]-2, xright = win_spr_break[i+1]-2, ybottom = par("usr")[3], ytop = par("usr")[4], 
       border = NA, col = adjustcolor("grey", alpha = 0.4))
}


axis(side=1, at=xtick-2, labels = F, lwd = 0, lwd.ticks = 1)
text(x=xtick-8,  par("usr")[3]-1.8, labels = seq_date[xtick], srt=45, cex=0.7, xpd=TRUE)

ytick <- c(0,2,4,6,8,10)
axis(side=2, at=ytick, labels = FALSE)
text(par("usr")[1]-1.7, ytick, labels=ytick, pos=2, xpd=TRUE, cex=0.8)


par(mar=c(0, 4, 0, 1), fig=c(0,1,0.67,0.74), new=T)
plot(NULL, ylim=c(0,1), xlim=c(1,tau-1), ylab="", xlab="", xaxt="n", yaxt="n")
for(i in result$est_CP){abline(v=i-2, col='red', lwd=2)}
text(par("usr")[1]+1, 0.45, labels='CPDstergm', pos=2, xpd=TRUE, cex=0.8)

par(mar=c(0, 4, 0, 1), fig=c(0,1,0.75,0.82), new=T)
plot(NULL, ylim=c(0,1), xlim=c(1,tau-1), ylab="", xlab="", xaxt="n", yaxt="n")
for(i in kerSeg_result){abline(v=i-2, col='blue', lwd=2)}
text(par("usr")[1]+1, 0.45, labels='kerSeg', pos=2, xpd=TRUE, cex=0.8)

par(mar=c(0, 4, 0, 1), fig=c(0,1,0.83,0.90), new=T)
plot(NULL, ylim=c(0,1), xlim=c(1,tau-1), ylab="", xlab="", xaxt="n", yaxt="n")
for(i in gSeg_result){abline(v=i-2, col='blue', lwd=2)}
text(par("usr")[1]+1, 0.45, labels='gSeg', pos=2, xpd=TRUE, cex=0.8)



#####################
# Stock Market Data #
#####################



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


result <- CPD_STERGM(df, directed=FALSE, network_stats=c("edges", "triangles"))

gSeg_result <- Evaluation_gSeg(df, p_threshold=0.001, is_experiment=TRUE)
kerSeg_result <- Evaluation_kerSeg(df, p_threshold=0.001, is_experiment=TRUE)



theta_change <- result$theta_change; threshold <- result$threshold; xtick <- result$est_CP
seq_date <- rownames(market[start:end,]); tau <- length(seq_date)-1

par(mar=c(4, 4, 2, 1), fig=c(0,1,0,0.74))
plot(1:length(theta_change), theta_change, type='l',ylab="", xlab="", xaxt="n", yaxt="n")
abline(h = threshold, col='red',lwd=2)
axis(side=1, at=xtick-2, labels = F, lwd = 0, lwd.ticks = 1)
text(x=xtick-6,  par("usr")[3]-1.3, labels = seq_date[xtick], srt=45, cex=0.7, xpd=TRUE)

ytick <- c(0,2,4,6)
axis(side=2, at=ytick, labels = FALSE)
text(par("usr")[1]-1.7, ytick, labels=ytick, pos=2, xpd=TRUE, cex=0.8)


par(mar=c(0, 4, 0, 1), fig=c(0,1,0.67,0.74), new=T)
plot(NULL, ylim=c(0,1), xlim=c(1,tau-1), ylab="", xlab="", xaxt="n", yaxt="n")
for(i in result$est_CP){abline(v=i-2, col='red', lwd=2)}
text(par("usr")[1]+1, 0.45, labels='CPDstergm', pos=2, xpd=TRUE, cex=0.8)

par(mar=c(0, 4, 0, 1), fig=c(0,1,0.75,0.82), new=T)
plot(NULL, ylim=c(0,1), xlim=c(1,tau-1), ylab="", xlab="", xaxt="n", yaxt="n")
for(i in kerSeg_result){abline(v=i-2, col='blue', lwd=2)}
text(par("usr")[1]+1, 0.45, labels='kerSeg', pos=2, xpd=TRUE, cex=0.8)

par(mar=c(0, 4, 0, 1), fig=c(0,1,0.83,0.90), new=T)
plot(NULL, ylim=c(0,1), xlim=c(1,tau-1), ylab="", xlab="", xaxt="n", yaxt="n")
for(i in gSeg_result){abline(v=i-2, col='blue', lwd=2)}
text(par("usr")[1]+1, 0.45, labels='gSeg', pos=2, xpd=TRUE, cex=0.8)




library(tergm)
y_data <- df

y_before93 <- y_after93 <- y_after121 <- list()
holder_after93 <- holder_before93 <- holder_after121 <- c()
duration <- 24


for(i in 1:duration){
  holder_before93 <- c(holder_before93, 93-(duration+1-i))
  y_before93[[i]] <- as.network(y_data[[93-(duration+1-i)]], directed = F)
  
  holder_after93 <- c(holder_after93,93-1+i)
  y_after93[[i]] <- as.network(y_data[[93-1+i]], directed = F)
  
  holder_after121 <- c(holder_after121, 121-1+i)
  y_after121[[i]] <- as.network(y_data[[121-1+i]], directed = F)
}


before93.dyn <- NetSeries(y_before93)
after93.dyn <- NetSeries(y_after93)
after121.dyn <- NetSeries(y_after121)

before93.fit <- tergm(
  before93.dyn ~
    Form(~edges + triangles) +
    Persist(~edges + triangles),
  estimate = "CMLE",
  times = c(1:duration)
)


after93.fit <- tergm(
  after93.dyn ~
    Form(~edges + triangles) +
    Persist(~edges + triangles),
  estimate = "CMLE",
  times = c(1:duration)
)


after121.fit <- tergm(
  after121.dyn ~
    Form(~edges + triangles) +
    Persist(~edges + triangles),
  estimate = "CMLE",
  times = c(1:duration)
)



summary(before93.fit)
summary(after93.fit)
summary(after121.fit)

