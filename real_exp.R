
#library(devtools)
#install_github("allenkei/CPDstergm", force = TRUE)
library(CPDstergm)

source("EVAL.R")


######################
# MIT Cellphone Data #
######################

data("MITphone")



node_degree <- rep(0, dim(MITphone[[1]])[1])
for (t in seq_len(length(MITphone))) {
  mat <- MITphone[[t]]
  node_degree <- node_degree + rowSums(mat)
}

degree_median <- median(node_degree)
node_status <- ifelse(node_degree > degree_median, "A", "N")


result <- CPD_STERGM(MITphone, directed=FALSE, network_stats=c("edges", "isolates", "triangles", "nodematch(\"Gender\")"), 
                     node_attr=node_status, list_of_lambda=10^c(0:4))


gSeg_result <- Evaluation_gSeg_on_stats(MITphone, p_threshold=0.05, num_stats=3, is_experiment=TRUE)
kerSeg_result <- Evaluation_kerSeg_on_stats(MITphone, p_threshold=0.001, num_stats=3, is_experiment=TRUE)
rdpg_result <- Evaluation_RDPG(MITphone, M=100, d=5, delta=5, is_experiment=TRUE)
nbs_result <- Evaluation_NBS(MITphone, M=15, delta=7, is_experiment=TRUE)


theta_change <- result$theta_change; threshold <- result$threshold; xtick <- result$est_CP
seq_date <- seq(as.Date("2004-09-15"), as.Date("2005-05-04"), by="days"); tau <- length(seq_date)-1


par(mar=c(4, 4, 2, 1), fig=c(0,1,0,0.66))
plot(1:length(theta_change), theta_change, type='l',ylab="", xlab="", xaxt="n", yaxt="n")
abline(h = threshold, col='red',lwd=2)

# winter and spring vacation (duration): seq_date[95];seq_date[110];seq_date[188];seq_date[192]
win_spr_break <- c(95,110,188,192)

for(i in c(1,3)){
  rect(xleft = win_spr_break[i]-2, xright = win_spr_break[i+1]-2, ybottom = par("usr")[3], ytop = par("usr")[4], 
       border = NA, col = adjustcolor("grey", alpha = 0.4))
}


axis(side=1, at=xtick-2, labels = F, lwd = 0, lwd.ticks = 1)
text(x=xtick-8,  par("usr")[3]-2, labels = seq_date[xtick], srt=45, cex=0.7, xpd=TRUE)

ytick <- c(0,2,4,6,8,10)
axis(side=2, at=ytick, labels = FALSE)
text(par("usr")[1]-1.7, ytick, labels=ytick, pos=2, xpd=TRUE, cex=0.8)



par(mar=c(0, 4, 0, 1), fig=c(0,1,0.59,0.66), new=T)
plot(NULL, ylim=c(0,1), xlim=c(1,tau), ylab="", xlab="", xaxt="n", yaxt="n")
for(i in nbs_result){abline(v=i-1, col='blue', lwd=2)}
text(par("usr")[1]+1, 0.45, labels='CPDnbs', pos=2, xpd=TRUE, cex=0.8)

par(mar=c(0, 4, 0, 1), fig=c(0,1,0.67,0.74), new=T)
plot(NULL, ylim=c(0,1), xlim=c(1,tau), ylab="", xlab="", xaxt="n", yaxt="n")
for(i in rdpg_result){abline(v=i-1, col='blue', lwd=2)}
text(par("usr")[1]+1, 0.45, labels='CPDrdpg', pos=2, xpd=TRUE, cex=0.8)

par(mar=c(0, 4, 0, 1), fig=c(0,1,0.75,0.82), new=T)
plot(NULL, ylim=c(0,1), xlim=c(1,tau), ylab="", xlab="", xaxt="n", yaxt="n")
for(i in kerSeg_result){abline(v=i-1, col='blue', lwd=2)}
text(par("usr")[1]+1, 0.45, labels='kerSeg', pos=2, xpd=TRUE, cex=0.8)

par(mar=c(0, 4, 0, 1), fig=c(0,1,0.83,0.9), new=T)
plot(NULL, ylim=c(0,1), xlim=c(1,tau), ylab="", xlab="", xaxt="n", yaxt="n")
for(i in gSeg_result){abline(v=i-1, col='blue', lwd=2)}
text(par("usr")[1]+1, 0.45, labels='gSeg', pos=2, xpd=TRUE, cex=0.8)



df <- list(result$est_CP, gSeg_result, kerSeg_result, rdpg_result, nbs_result)
#save(df, file = 'MITphone_result.Rdata')



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


result <- CPD_STERGM(df, directed=FALSE, network_stats=c("edges", "triangles"), list_of_lambda=10^c(0:4))

# Use network stat has no CP for gSeg and kerSeg
gSeg_result <- Evaluation_gSeg(df, p_threshold=0.001, is_experiment=TRUE)
kerSeg_result <- Evaluation_kerSeg(df, p_threshold=0.001, is_experiment=TRUE)
rdpg_result <- Evaluation_RDPG(df, M=100, d=5, delta=5, is_experiment=TRUE)
nbs_result <- Evaluation_NBS(df, M=15, delta=7, is_experiment=TRUE)



theta_change <- result$theta_change; threshold <- result$threshold; xtick <- result$est_CP
seq_date <- rownames(market[start:end,]); tau <- length(seq_date)-1

par(mar=c(4, 4, 2, 1), fig=c(0,1,0,0.66))
plot(1:length(theta_change), theta_change, type='l',ylab="", xlab="", xaxt="n", yaxt="n")
abline(h = threshold, col='red',lwd=2)
axis(side=1, at=xtick-2, labels = F, lwd = 0, lwd.ticks = 1)
text(x=xtick-6,  par("usr")[3]-1.3, labels = seq_date[xtick], srt=45, cex=0.7, xpd=TRUE)

ytick <- c(0,2,4,6)
axis(side=2, at=ytick, labels = FALSE)
text(par("usr")[1]-1.7, ytick, labels=ytick, pos=2, xpd=TRUE, cex=0.8)


par(mar=c(0, 4, 0, 1), fig=c(0,1,0.59,0.66), new=T)
plot(NULL, ylim=c(0,1), xlim=c(1,tau), ylab="", xlab="", xaxt="n", yaxt="n")
for(i in nbs_result){abline(v=i-1, col='blue', lwd=2)}
text(par("usr")[1]+1, 0.45, labels='CPDnbs', pos=2, xpd=TRUE, cex=0.8)

par(mar=c(0, 4, 0, 1), fig=c(0,1,0.67,0.74), new=T)
plot(NULL, ylim=c(0,1), xlim=c(1,tau), ylab="", xlab="", xaxt="n", yaxt="n")
for(i in rdpg_result){abline(v=i-1, col='blue', lwd=2)}
text(par("usr")[1]+1, 0.45, labels='CPDrdpg', pos=2, xpd=TRUE, cex=0.8)

par(mar=c(0, 4, 0, 1), fig=c(0,1,0.75,0.82), new=T)
plot(NULL, ylim=c(0,1), xlim=c(1,tau), ylab="", xlab="", xaxt="n", yaxt="n")
for(i in kerSeg_result){abline(v=i-1, col='blue', lwd=2)}
text(par("usr")[1]+1, 0.45, labels='kerSeg', pos=2, xpd=TRUE, cex=0.8)

par(mar=c(0, 4, 0, 1), fig=c(0,1,0.83,0.9), new=T)
plot(NULL, ylim=c(0,1), xlim=c(1,tau), ylab="", xlab="", xaxt="n", yaxt="n")
for(i in gSeg_result){abline(v=i-1, col='blue', lwd=2)}
text(par("usr")[1]+1, 0.45, labels='gSeg', pos=2, xpd=TRUE, cex=0.8)


df <- list(result$est_CP, gSeg_result, kerSeg_result, rdpg_result, nbs_result)
#save(df, file = 'Stock_result.Rdata')




