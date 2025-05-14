
#library(devtools)
#install_github("allenkei/CPDstergm", force = TRUE)
library(CPDstergm)

source("EVAL.R")


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

