library(devtools)
install_github("allenkei/CPDstergm")
library(CPDstergm)

# 9 by 6 inches

set.seed(1)
SBM_rho05 <- sim_SBM_list(num_seq=1, n=100, rho=0.5)

par(mfrow=c(2,3))

for(i in c(25,50,75, 26,51,76)){
  par(mar = c(2.2, 2.2, 2.2, 2.2)); 
  image(SBM_rho05[[1]][[i]], xaxt = "n", yaxt = "n")
  axis(side=1,at=seq(0,1,length.out = 4),labels=round(seq(1,100,length.out = 4)),xpd=NA,cex.axis=1.4)
  axis(side=2,at=seq(0,1,length.out = 4),labels=round(seq(1,100,length.out = 4)),xpd=NA,cex.axis=1.4)
}





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

par(mfrow=c(2,3))

for(i in c(25,50,75, 26,51,76)){
  par(mar = c(2.2, 2.2, 2.2, 2.2)); 
  image(STERGM3[[1]][[i]], xaxt = "n", yaxt = "n")
  axis(side=1,at=seq(0,1,length.out = 3),labels=round(seq(1,100,length.out = 3)),xpd=NA,cex.axis=1.4)
  axis(side=2,at=seq(0,1,length.out = 3),labels=round(seq(1,100,length.out = 3)),xpd=NA,cex.axis=1.4)  
}










