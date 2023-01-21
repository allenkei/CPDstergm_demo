library(ecp)
library(ggplot2)
library(ggpubr)
library(reshape)
library(dplyr)

data(DJIA)
market <- DJIA$market
date_vec <- DJIA$dates[1:1138]
rownames(market) <- date_vec
selection <- market[266:109,] # 2007-01-01 to 2010-01-04

# index 93 is 2008-10-06
holder_after93 <- holder_before93 <- c()
return_before <- return_after <- return4_20 <- return8_16 <- return12_12 <- return16_8 <- return20_4 <- c()
duration <- 24
ones <- rep(1,10)
number_of_portfolio <- 100
number_of_stock <- 10


for(i in 1:duration){
  holder_before93 <- c(holder_before93, 93-(duration+1-i))
  holder_after93 <- c(holder_after93,93-1+i)
}


set.seed(1)
for(iter in 1:number_of_portfolio ){
  
  companies <- sample(29,number_of_stock) # randomly choose k companies
  
  
  full_log_return <- selection[c(holder_before93[5:24],holder_after93[1:4]), companies]
  full_return <- exp(full_log_return)-1
  full_cov <- cov(full_return)
  full_portfolio <- solve(full_cov) %*% ones / as.numeric(t(ones) %*% solve(full_cov) %*% ones)
  r_f <- t(full_portfolio) %*% colMeans(full_return) # Expected return 
  return20_4 <- c(r_f, return20_4)
  
  
  full_log_return <- selection[c(holder_before93[9:24],holder_after93[1:8]), companies]
  full_return <- exp(full_log_return)-1
  full_cov <- cov(full_return)
  full_portfolio <- solve(full_cov) %*% ones / as.numeric(t(ones) %*% solve(full_cov) %*% ones)
  r_f <- t(full_portfolio) %*% colMeans(full_return) # Expected return 
  return16_8 <- c(r_f, return16_8)
  
  
  full_log_return <- selection[c(holder_before93[13:24],holder_after93[1:12]), companies]
  full_return <- exp(full_log_return)-1
  full_cov <- cov(full_return)
  full_portfolio <- solve(full_cov) %*% ones / as.numeric(t(ones) %*% solve(full_cov) %*% ones)
  r_f <- t(full_portfolio) %*% colMeans(full_return) # Expected return 
  return12_12 <- c(r_f, return12_12)
  
  
  full_log_return <- selection[c(holder_before93[17:24],holder_after93[1:16]), companies]
  full_return <- exp(full_log_return)-1
  full_cov <- cov(full_return)
  full_portfolio <- solve(full_cov) %*% ones / as.numeric(t(ones) %*% solve(full_cov) %*% ones)
  r_f <- t(full_portfolio) %*% colMeans(full_return) # Expected return 
  return8_16 <- c(r_f, return8_16)
  
  
  full_log_return <- selection[c(holder_before93[21:24],holder_after93[1:20]), companies]
  full_return <- exp(full_log_return)-1
  full_cov <- cov(full_return)
  full_portfolio <- solve(full_cov) %*% ones / as.numeric(t(ones) %*% solve(full_cov) %*% ones)
  r_f <- t(full_portfolio) %*% colMeans(full_return) # Expected return 
  return4_20 <- c(r_f, return4_20)
  
  
  ####################
  # Before and After #
  ####################
  
  before93_log_return <- selection[holder_before93, companies]
  before93_return <- exp(before93_log_return)-1
  before93_cov <- cov(before93_return)
  before93_portfolio <- solve(before93_cov) %*% ones / as.numeric(t(ones) %*% solve(before93_cov) %*% ones)
  r_b <- t(before93_portfolio) %*% colMeans(before93_return) # Expected return 
  return_before <- c(r_b, return_before)
  
  
  after93_log_return <- selection[holder_after93, companies]
  after93_return <- exp(after93_log_return)-1
  after93_cov <- cov(after93_return)
  after93_portfolio <- solve(after93_cov) %*% ones / as.numeric(t(ones) %*% solve(after93_cov) %*% ones)
  r_a <- t(after93_portfolio) %*% colMeans(after93_return) # Expected return 
  return_after <- c(r_a, return_after)
  
}


data <- data.frame(returns = 10*c(mean(return_before), mean(return20_4), mean(return16_8), mean(return12_12), 
                               mean(return8_16), mean(return4_20), mean(return_after)), 
                   type = c('24:0','20:4','16:8','12:12','8:16','4:20','0:24'),
                   cols = c('deepskyblue','dimgray','dimgray','dimgray','dimgray','dimgray','orange'))


#barplot(data$returns)


p3 <- ggplot(data=data, aes(x=type, y=returns, fill=cols)) + 
  geom_bar(stat="identity", width=0.5) + scale_fill_identity() +
  scale_y_continuous(limits = c(-0.05, 0.05), breaks = c(-0.05,-0.025,0, 0.025, 0.05), labels = scales::percent(c(-0.05,-0.025,0, 0.025, 0.05))) + 
  ylab("return (x10)") + xlab("combination of weeks") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) + 
  ggtitle("mean of 100 portfolio returns") +
  scale_x_discrete(limits=c('24:0','20:4','16:8','12:12','8:16','4:20','0:24')) 



df1 <- selection[holder_before93,]
df <- data.frame(week = seq_along(df1[, 1]), df1)
df <- melt(df, id.vars = "week")
p1 <- ggplot(df, aes(x = week, y = value, color = variable)) + geom_line() +
  scale_color_manual(values = alpha(rep('deepskyblue',29),0.6)) +
  scale_y_continuous(limits = c(-0.6, 0.6)) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) + ylab("log return") +
  ggtitle("before 2008-10-06")


df1 <- selection[holder_after93,]
df <- data.frame(week = seq_along(df1[, 1]), df1)
df <- melt(df, id.vars = "week")
p2 <- ggplot(df, aes(x = week, y = value, color = variable)) + geom_line() +
  scale_color_manual(values = alpha(rep('orange',29),0.6)) +
  scale_y_continuous(limits = c(-0.6, 0.6)) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) + ylab("log return") +
  ggtitle("after 2008-10-06")



ggarrange(p1, p2, p3, ncol = 3, nrow = 1, widths = c(0.95, 0.95, 1.05))



