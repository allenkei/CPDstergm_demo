
################################
# Evaluation Function (kerSeg) #
################################


library(kerSeg)

Evaluation_kerSeg <- function(y_list, p_threshold, is_experiment=FALSE, true_CP=c(26, 51, 76)){
  
  binary_search <- function(temp, holder, lower_index, upper_index, p_threshold){
    
    if(dim(temp)[1] > 5){
      
      K = gaussiankernel(temp) # Gaussian kernel matrix
      a = kerseg1(dim(temp)[1], K, pval.perm=TRUE, B=1000)
      cur_index <- a$tauhat 
      cat('lower: ', lower_index, '\n')
      cat('upper: ', upper_index, '\n')
      cat('dim: ', dim(temp)[1], '\n')
      cat('est_CP: ', cur_index, '\n')
      cat('actual index: ', lower_index + cur_index - 1, '\n')
      cat('p-value: ',a$appr$fGKCP1_bon, '\n')
      
      if(a$appr$fGKCP1_bon < p_threshold){
        holder <- c(holder, lower_index+cur_index-1)
        cat('holder: ', holder, '\n')
        
        cat('left\n')
        temp1 <- temp[1:cur_index,]
        holder <- binary_search(temp1, holder, lower_index, lower_index+cur_index-1, p_threshold)
        cat('right\n')
        temp2 <- temp[(cur_index+1):dim(temp)[1],]
        holder <- binary_search(temp2, holder, lower_index+cur_index, upper_index, p_threshold)
      }
      
    }else{cat('less than 5\n')}
    
    return(holder)
    
  }
  
  if(is_experiment){
    
    y_data <- y_list
    num_time <- length(y_data)
    n <- dim(y_data[[1]])[1]
    tau <- num_time-1 
    
    temp <- matrix(NA, nrow=num_time, ncol=n*n)
    for(iter in 1:length(y_data)){temp[iter,] <- c(y_data[[iter]])}
    
    est_CP <- binary_search(temp, c(), 1, num_time, p_threshold); rm(temp)
    est_CP <- sort(est_CP)
    est_CP <- est_CP + 1
    return(est_CP)
    
  }else{
    
    repetition_iteration <- length(y_list)
    result <- matrix(NA, nrow=repetition_iteration, ncol = 4)
    
    for(rep_iter in 1:repetition_iteration){
      
      y_data <- y_list[[rep_iter]]
      num_time <- length(y_data)
      n <- dim(y_data[[1]])[1]
      tau <- num_time-1 # i.e. T=150, then num_theta=149 (tau), and num_beta=148
      
      temp <- matrix(NA, nrow=num_time, ncol=n*n)
      for(iter in 1:length(y_data)){temp[iter,] <- c(y_data[[iter]])}
      
      est_CP <- binary_search(temp, c(), 1, num_time, p_threshold); rm(temp)
      est_CP <- sort(est_CP)
      
      est_CP <- est_CP + 1
      num_CP <- length(est_CP)
      
      gt_CP_corrected <- c(1, true_CP, tau+2) # tau+2 = T+1
      est_CP_corrected <- c(1, est_CP, tau+2)
      
      gt_list <- est_list <- list();
      for(i in 2:length(gt_CP_corrected)){
        gt_list[[i-1]] <- gt_CP_corrected[i-1]:(gt_CP_corrected[i]-1)
      }
      for(i in 2:length(est_CP_corrected)){
        est_list[[i-1]] <- est_CP_corrected[i-1]:(est_CP_corrected[i]-1)
      }
      
      if(num_CP == 0){
        dist_est_gt <- Inf
        dist_gt_est <- -Inf
        covering_metric <- 0
      }else{
        
        holder <- c()
        for(i in true_CP){
          dist_diff <- c()
          for(j in est_CP){dist_diff <- c(dist_diff, abs(j-i))}
          holder <- c(holder, min(dist_diff))
        }
        dist_est_gt <- max(holder)
        
        holder <- c()
        for(i in est_CP){
          dist_diff <- c()
          for(j in true_CP){dist_diff <- c(dist_diff, abs(j-i))}
          holder <- c(holder, min(dist_diff))
        }
        dist_gt_est <- max(holder)
        
        covering_metric <- 0
        for(i in 1:length(gt_list)){
          A <- gt_list[[i]]
          jaccard <- c()
          for(j in 1:length(est_list)){
            A_prime <- est_list[[j]]
            jaccard <- c(jaccard,length(intersect(A,A_prime))/length(union(A,A_prime)))
          }
          covering_metric <- covering_metric + length(A)*max(jaccard)
        }
        covering_metric <- covering_metric/(tau+1) # tau+1 = T
      }
      
      abs_error <- abs(num_CP - length(true_CP))
      
      result[rep_iter, 1] <- abs_error
      result[rep_iter, 2] <- dist_est_gt
      result[rep_iter, 3] <- dist_gt_est
      result[rep_iter, 4] <- covering_metric
    }
    
    return(result)
    
  }
  
}


##############################
# Evaluation Function (gSeg) #
##############################


library(gSeg)
library(ade4)

Evaluation_gSeg <- function(y_list, p_threshold, is_experiment=FALSE, true_CP=c(26, 51, 76)){
  
  binary_search <- function(temp, holder, lower_index, upper_index, p_threshold){
    
    if(dim(temp)[1] > 5){
      
      temp_dist <- dist(temp)
      temp_mst <- mstree(temp_dist)
      temp_result <- gseg1(dim(temp)[1], temp_mst, statistics="o") # "m" "w" "g"
      cur_index <- temp_result$scanZ$ori$tauhat # temp_result$scanZ$max$tauhat # temp_result$scanZ$weighted$tauhat # temp_result$scanZ$generalized$tauhat
      p_value <- temp_result$pval.appr$ori # temp_result$pval.appr$max # temp_result$pval.appr$weighted # temp_result$pval.appr$generalized
      
      cat('lower: ', lower_index, '\n')
      cat('upper: ', upper_index, '\n')
      cat('dim: ', dim(temp)[1], '\n')
      cat('est_CP: ', cur_index, '\n')
      cat('actual index: ', lower_index + cur_index - 1, '\n')
      cat('p-value: ',p_value, '\n')
      
      if(p_value < p_threshold){
        holder <- c(holder, lower_index+cur_index-1)
        cat('holder: ', holder, '\n')
        
        cat('left\n')
        temp1 <- temp[1:cur_index,]
        holder <- binary_search(temp1, holder, lower_index, lower_index+cur_index-1, p_threshold)
        cat('right\n')
        temp2 <- temp[(cur_index+1):dim(temp)[1],]
        holder <- binary_search(temp2, holder, lower_index+cur_index, upper_index, p_threshold)
      }
      
    }else{cat('less than 5\n')}
    
    return(holder)
    
  }
  
  if(is_experiment){
    
    y_data <- y_list
    num_time <- length(y_data)
    n <- dim(y_data[[1]])[1]
    tau <- num_time-1 
    
    temp <- matrix(NA, nrow=num_time, ncol=n*n)
    for(iter in 1:length(y_data)){temp[iter,] <- c(y_data[[iter]])}
    
    est_CP <- binary_search(temp, c(), 1, num_time, p_threshold); rm(temp)
    est_CP <- sort(est_CP)
    est_CP <- est_CP + 1
    return(est_CP)
    
  }else{
    repetition_iteration <- length(y_list)
    result <- matrix(NA, nrow=repetition_iteration, ncol = 4)
    
    for(rep_iter in 1:repetition_iteration){
      
      y_data <- y_list[[rep_iter]]
      num_time <- length(y_data)
      n <- dim(y_data[[1]])[1]
      tau <- num_time-1 # i.e. T=150, then num_theta=149 (tau), and num_beta=148
      
      temp <- matrix(NA, nrow=num_time, ncol=n*n)
      for(iter in 1:length(y_data)){temp[iter,] <- c(y_data[[iter]])}
      
      est_CP <- binary_search(temp, c(), 1, num_time, p_threshold); rm(temp)
      est_CP <- sort(est_CP)
      
      est_CP <- est_CP + 1
      num_CP <- length(est_CP)
      
      gt_CP_corrected <- c(1, true_CP, tau+2) # tau+2 = T+1
      est_CP_corrected <- c(1, est_CP, tau+2)
      
      gt_list <- est_list <- list();
      for(i in 2:length(gt_CP_corrected)){
        gt_list[[i-1]] <- gt_CP_corrected[i-1]:(gt_CP_corrected[i]-1)
      }
      for(i in 2:length(est_CP_corrected)){
        est_list[[i-1]] <- est_CP_corrected[i-1]:(est_CP_corrected[i]-1)
      }
      
      if(num_CP == 0){
        dist_est_gt <- Inf
        dist_gt_est <- -Inf
        covering_metric <- 0
      }else{
        
        holder <- c()
        for(i in true_CP){
          dist_diff <- c()
          for(j in est_CP){dist_diff <- c(dist_diff, abs(j-i))}
          holder <- c(holder, min(dist_diff))
        }
        dist_est_gt <- max(holder)
        
        holder <- c()
        for(i in est_CP){
          dist_diff <- c()
          for(j in true_CP){dist_diff <- c(dist_diff, abs(j-i))}
          holder <- c(holder, min(dist_diff))
        }
        dist_gt_est <- max(holder)
        
        covering_metric <- 0
        for(i in 1:length(gt_list)){
          A <- gt_list[[i]]
          jaccard <- c()
          for(j in 1:length(est_list)){
            A_prime <- est_list[[j]]
            jaccard <- c(jaccard,length(intersect(A,A_prime))/length(union(A,A_prime)))
          }
          covering_metric <- covering_metric + length(A)*max(jaccard)
        }
        covering_metric <- covering_metric/(tau+1) # tau+1 = T
        
      }
      
      abs_error <- abs(num_CP - length(true_CP))
      
      result[rep_iter, 1] <- abs_error
      result[rep_iter, 2] <- dist_est_gt
      result[rep_iter, 3] <- dist_gt_est
      result[rep_iter, 4] <- covering_metric
      
    }
    
    return(result)
  }
  
}


######################################################
# Evaluation Function on Network Statistics (kerSeg) #
######################################################


Evaluation_kerSeg_on_stats <- function(y_list, p_threshold, num_stats, is_experiment=FALSE, true_CP=c(26, 51, 76)){
  
  library(ergm)
  binary_search <- function(temp, holder, lower_index, upper_index, p_threshold){
    
    if(dim(temp)[1] > 5){
      
      K = gaussiankernel(temp) # Gaussian kernel matrix
      a = kerseg1(dim(temp)[1], K, pval.perm=TRUE, B=1000)
      cur_index <- a$tauhat
      cat('lower: ', lower_index, '\n')
      cat('upper: ', upper_index, '\n')
      cat('dim: ', dim(temp)[1], '\n')
      cat('est_CP: ', cur_index, '\n')
      cat('actual index: ', lower_index + cur_index - 1, '\n')
      cat('p-value: ',a$appr$fGKCP1_bon, '\n')
      
      if(a$appr$fGKCP1_bon < p_threshold){
        holder <- c(holder, lower_index+cur_index-1)
        cat('holder: ', holder, '\n')
        
        cat('left\n')
        temp1 <- temp[1:cur_index,]
        holder <- binary_search(temp1, holder, lower_index, lower_index+cur_index-1, p_threshold)
        cat('right\n')
        temp2 <- temp[(cur_index+1):dim(temp)[1],]
        holder <- binary_search(temp2, holder, lower_index+cur_index, upper_index, p_threshold)
      }
      
    }else{cat('less than 5\n')}
    
    return(holder)
    
  }
  
  if(is_experiment){
    
    y_data <- y_list
    num_time <- length(y_data)
    n <- dim(y_data[[1]])[1]
    tau <- num_time-1 
    
    temp <- matrix(NA, nrow=num_time, ncol=num_stats)
    for(t in 1:num_time){
      if(num_stats == 2){
        yt <- network(y_data[[t]], directed = F)
        temp[t,] = summary(yt ~ edges + triangles)
      }else if(num_stats == 3){
        yt <- network(y_data[[t]], directed = F)
        temp[t,] = summary(yt ~ edges + isolates + triangles)
      }
    }
    
    est_CP <- binary_search(temp, c(), 1, num_time, p_threshold); rm(temp)
    est_CP <- sort(est_CP)
    est_CP <- est_CP + 1
    return(est_CP)
    
  }else{
    
    repetition_iteration <- length(y_list)
    result <- matrix(NA, nrow=repetition_iteration, ncol = 4)
    
    for(rep_iter in 1:repetition_iteration){
      
      y_data <- y_list[[rep_iter]]
      num_time <- length(y_data)
      n <- dim(y_data[[1]])[1]
      tau <- num_time-1 # i.e. T=150, then num_theta=149 (tau), and num_beta=148
      
      temp <- matrix(NA, nrow=num_time, ncol=num_stats)
      for(t in 1:num_time){
        if(num_stats == 2){
          temp[t,] = summary(y_data[[t]] ~ edges + mutual)
        }else if(num_stats == 3){
          temp[t,] = summary(y_data[[t]] ~ edges + mutual + triangles)
        }else if(num_stats == 4){
          yt <- network(y_data[[t]])
          network::set.vertex.attribute(yt, "Gender", gender)
          temp[t,] = summary(yt ~ edges + mutual + triangles + nodematch("Gender"))
        }
      }
      
      est_CP <- binary_search(temp, c(), 1, num_time, p_threshold); rm(temp)
      est_CP <- sort(est_CP)
      
      est_CP <- est_CP + 1
      num_CP <- length(est_CP)
      
      gt_CP_corrected <- c(1, true_CP, tau+2) # tau+2 = T+1
      est_CP_corrected <- c(1, est_CP, tau+2)
      
      gt_list <- est_list <- list();
      for(i in 2:length(gt_CP_corrected)){
        gt_list[[i-1]] <- gt_CP_corrected[i-1]:(gt_CP_corrected[i]-1)
      }
      for(i in 2:length(est_CP_corrected)){
        est_list[[i-1]] <- est_CP_corrected[i-1]:(est_CP_corrected[i]-1)
      }
      
      if(num_CP == 0){
        dist_est_gt <- Inf
        dist_gt_est <- -Inf
        covering_metric <- 0
      }else{
        
        holder <- c()
        for(i in true_CP){
          dist_diff <- c()
          for(j in est_CP){dist_diff <- c(dist_diff, abs(j-i))}
          holder <- c(holder, min(dist_diff))
        }
        dist_est_gt <- max(holder)
        
        holder <- c()
        for(i in est_CP){
          dist_diff <- c()
          for(j in true_CP){dist_diff <- c(dist_diff, abs(j-i))}
          holder <- c(holder, min(dist_diff))
        }
        dist_gt_est <- max(holder)
        
        covering_metric <- 0
        for(i in 1:length(gt_list)){
          A <- gt_list[[i]]
          jaccard <- c()
          for(j in 1:length(est_list)){
            A_prime <- est_list[[j]]
            jaccard <- c(jaccard,length(intersect(A,A_prime))/length(union(A,A_prime)))
          }
          covering_metric <- covering_metric + length(A)*max(jaccard)
        }
        covering_metric <- covering_metric/(tau+1) # tau+1 = T
        
      }
      
      abs_error <- abs(num_CP - length(true_CP))
      
      result[rep_iter, 1] <- abs_error
      result[rep_iter, 2] <- dist_est_gt
      result[rep_iter, 3] <- dist_gt_est
      result[rep_iter, 4] <- covering_metric
      
    }
    
    return(result)
  }
  
}


####################################################
# Evaluation Function on Network Statistics (gSeg) #
####################################################


Evaluation_gSeg_on_stats <- function(y_list, p_threshold, num_stats, is_experiment=FALSE, true_CP=c(26, 51, 76)){
  
  library(ergm)
  
  binary_search <- function(temp, holder, lower_index, upper_index, p_threshold){
    
    if(dim(temp)[1] > 5){
      
      temp_dist <- dist(temp)
      temp_mst <- mstree(temp_dist)
      temp_result <- gseg1(dim(temp)[1], temp_mst, statistics="o") # "o" "m" "w" "g"
      # temp_result$scanZ$ori$tauhat # temp_result$scanZ$max$tauhat # temp_result$scanZ$weighted$tauhat # temp_result$scanZ$generalized$tauhat
      # temp_result$pval.appr$ori    # temp_result$pval.appr$max    # temp_result$pval.appr$weighted    # temp_result$pval.appr$generalized
      cur_index <- temp_result$scanZ$ori$tauhat
      p_value <- temp_result$pval.appr$ori
      
      cat('lower: ', lower_index, '\n')
      cat('upper: ', upper_index, '\n')
      cat('dim: ', dim(temp)[1], '\n')
      cat('est_CP: ', cur_index, '\n')
      cat('actual index: ', lower_index + cur_index - 1, '\n')
      cat('p-value: ',p_value, '\n')
      
      if(p_value < p_threshold){
        holder <- c(holder, lower_index+cur_index-1)
        cat('holder: ', holder, '\n')
        
        cat('left\n')
        temp1 <- temp[1:cur_index,]
        holder <- binary_search(temp1, holder, lower_index, lower_index+cur_index-1, p_threshold)
        cat('right\n')
        temp2 <- temp[(cur_index+1):dim(temp)[1],]
        holder <- binary_search(temp2, holder, lower_index+cur_index, upper_index, p_threshold)
      }
      
    }else{cat('less than 5\n')}
    
    return(holder)
    
  }
  
  if(is_experiment){
    
    y_data <- y_list
    num_time <- length(y_data)
    n <- dim(y_data[[1]])[1]
    tau <- num_time-1
    
    temp <- matrix(NA, nrow=num_time, ncol=num_stats)
    for(t in 1:num_time){
      if(num_stats == 2){
        yt <- network(y_data[[t]], directed = F)
        temp[t,] = summary(yt ~ edges + triangles)
      }else if(num_stats == 3){
        yt <- network(y_data[[t]], directed = F)
        temp[t,] = summary(yt ~ edges + isolates + triangles)
      }
    }
    
    est_CP <- binary_search(temp, c(), 1, num_time, p_threshold); rm(temp)
    est_CP <- sort(est_CP)
    est_CP <- est_CP + 1
    return(est_CP)
    
  }else{
    
    repetition_iteration <- length(y_list)
    result <- matrix(NA, nrow=repetition_iteration, ncol = 4)
    
    for(rep_iter in 1:repetition_iteration){
      
      y_data <- y_list[[rep_iter]]
      num_time <- length(y_data)
      n <- dim(y_data[[1]])[1]
      tau <- num_time-1 # i.e. T=150, then num_theta=149 (tau), and num_beta=148
      
      temp <- matrix(NA, nrow=num_time, ncol=num_stats)
      for(t in 1:num_time){
        if(num_stats == 2){
          temp[t,] = summary(y_data[[t]] ~ edges + mutual)
        }else if(num_stats == 3){
          temp[t,] = summary(y_data[[t]] ~ edges + mutual + triangles)
        }else if(num_stats == 4){
          yt <- network(y_data[[t]])
          network::set.vertex.attribute(yt, "Gender", gender)
          temp[t,] = summary(yt ~ edges + mutual + triangles + nodematch("Gender"))
        }
      }
      
      est_CP <- binary_search(temp, c(), 1, num_time, p_threshold); rm(temp)
      est_CP <- sort(est_CP)
      
      est_CP <- est_CP + 1
      num_CP <- length(est_CP)
      
      gt_CP_corrected <- c(1, true_CP, tau+2) # tau+2 = T+1
      est_CP_corrected <- c(1, est_CP, tau+2)
      
      gt_list <- est_list <- list();
      for(i in 2:length(gt_CP_corrected)){
        gt_list[[i-1]] <- gt_CP_corrected[i-1]:(gt_CP_corrected[i]-1)
      }
      for(i in 2:length(est_CP_corrected)){
        est_list[[i-1]] <- est_CP_corrected[i-1]:(est_CP_corrected[i]-1)
      }
      
      if(num_CP == 0){
        dist_est_gt <- Inf
        dist_gt_est <- -Inf
        covering_metric <- 0
      }else{
        
        holder <- c()
        for(i in true_CP){
          dist_diff <- c()
          for(j in est_CP){dist_diff <- c(dist_diff, abs(j-i))}
          holder <- c(holder, min(dist_diff))
        }
        dist_est_gt <- max(holder)
        
        holder <- c()
        for(i in est_CP){
          dist_diff <- c()
          for(j in true_CP){dist_diff <- c(dist_diff, abs(j-i))}
          holder <- c(holder, min(dist_diff))
        }
        dist_gt_est <- max(holder)
        
        covering_metric <- 0
        for(i in 1:length(gt_list)){
          A <- gt_list[[i]]
          jaccard <- c()
          for(j in 1:length(est_list)){
            A_prime <- est_list[[j]]
            jaccard <- c(jaccard,length(intersect(A,A_prime))/length(union(A,A_prime)))
          }
          covering_metric <- covering_metric + length(A)*max(jaccard)
        }
        covering_metric <- covering_metric/(tau+1) # tau+1 = T
        
      }
      
      abs_error <- abs(num_CP - length(true_CP))
      
      result[rep_iter, 1] <- abs_error
      result[rep_iter, 2] <- dist_est_gt
      result[rep_iter, 3] <- dist_gt_est
      result[rep_iter, 4] <- covering_metric
      
    }
    
    return(result)
  }
  
}





##############################
# Evaluation Function (RDPG) #
##############################


Evaluation_RDPG <- function(y_list, M, d, delta, is_experiment=FALSE, true_CP=c(26, 51, 76)){
  
  library(changepoints)
  
  RDPG_CPD <- function(one_seq, M.=M, d.=d, delta.=delta){
    # one_seq: list of 100, with matrix n by n
    # M: number of random intervals for WBS
    # d: parameter for scaled PCA algorithm
    set.seed(1)
    
    one_seq <- do.call(cbind, lapply(one_seq, as.vector)) # matrix of n*n by 100
    intervals = WBS.intervals(M=M., lower = 1, upper = ncol(one_seq))
    WBS_result = WBS.nonpar.RDPG(one_seq, lowerdiag = FALSE, d.,
                                 Alpha = intervals$Alpha, Beta = intervals$Beta, delta)
    cpt_hat = tuneBSnonparRDPG(WBS_result, one_seq, lowerdiag = FALSE, d.)
    return(cpt_hat)
  }
  
  
  if(is_experiment){
    
    est_CP <- RDPG_CPD(y_list, M.=M, d.=d, delta.=delta)
    return(est_CP)
    
  }else{
    # y_list: 10 by 100 by n by n
    repetition_iteration <- length(y_list)
    result <- matrix(NA, nrow=repetition_iteration, ncol = 4)
    
    for(rep_iter in 1:repetition_iteration){
      
      y_data <- y_list[[rep_iter]]
      
      est_CP <- RDPG_CPD(y_data)
      num_CP <- length(est_CP)
      
      gt_CP_corrected <- c(1, true_CP, 100) 
      est_CP_corrected <- c(1, est_CP, 100)
      
      gt_list <- est_list <- list();
      for(i in 2:length(gt_CP_corrected)){
        gt_list[[i-1]] <- gt_CP_corrected[i-1]:(gt_CP_corrected[i]-1)
      }
      for(i in 2:length(est_CP_corrected)){
        est_list[[i-1]] <- est_CP_corrected[i-1]:(est_CP_corrected[i]-1)
      }
      
      if(num_CP == 0){
        dist_est_gt <- Inf
        dist_gt_est <- -Inf
        covering_metric <- 0
      }else{
        
        holder <- c()
        for(i in true_CP){
          dist_diff <- c()
          for(j in est_CP){dist_diff <- c(dist_diff, abs(j-i))}
          holder <- c(holder, min(dist_diff))
        }
        dist_est_gt <- max(holder)
        
        holder <- c()
        for(i in est_CP){
          dist_diff <- c()
          for(j in true_CP){dist_diff <- c(dist_diff, abs(j-i))}
          holder <- c(holder, min(dist_diff))
        }
        dist_gt_est <- max(holder)
        
        covering_metric <- 0
        for(i in 1:length(gt_list)){
          A <- gt_list[[i]]
          jaccard <- c()
          for(j in 1:length(est_list)){
            A_prime <- est_list[[j]]
            jaccard <- c(jaccard,length(intersect(A,A_prime))/length(union(A,A_prime)))
          }
          covering_metric <- covering_metric + length(A)*max(jaccard)
        }
        covering_metric <- covering_metric/100
      }
      
      abs_error <- abs(num_CP - length(true_CP))
      
      result[rep_iter, 1] <- abs_error
      result[rep_iter, 2] <- dist_est_gt
      result[rep_iter, 3] <- dist_gt_est
      result[rep_iter, 4] <- covering_metric
    }
    
    return(result)
    
  }
  
}



#############################
# Evaluation Function (NBS) #
#############################


Evaluation_NBS <- function(y_list, M, delta, is_experiment=FALSE, true_CP=c(26, 51, 76)){
  
  library(changepoints)
  set.seed(1)
  
  NBS_CPD <- function(one_seq, M.=M, delta.=delta){
    # one_seq: list of 100, with matrix n by n
    # M: number of random intervals for WBS
    
    n <- length(one_seq) # T = 100
    p <- dim(one_seq[[1]])[1] # num_node = 50 or 100
    
    data_mat <- do.call(cbind, lapply(one_seq, as.vector)) # matrix of n*n by 100
    data_mat1 = data_mat[,seq(1,ncol(data_mat),2)]
    data_mat2 = data_mat[,seq(2,ncol(data_mat),2)]
    
    intervals = WBS.intervals(M = M., lower = 1, upper = ncol(data_mat1))
    temp = WBS.network(data_mat1, data_mat2, 1, ncol(data_mat1), intervals$Alpha, intervals$Beta, delta = delta.)
    rho_hat = quantile(rowMeans(data_mat), 0.95)
    tau = p*rho_hat*(log(n))^2/20 # default threshold given in the paper
    cpt_init = unlist(thresholdBS(temp, tau)$cpt_hat[,1])
    cpt_WBS = 2*cpt_init
    return(sort(cpt_WBS))
  }
  
  
  if(is_experiment){
    
    est_CP <- NBS_CPD(y_list, M.=M, delta.=delta)
    return(est_CP)
    
  }else{
    # y_list: 10 by 100 by n by n
    
    repetition_iteration <- length(y_list)
    result <- matrix(NA, nrow=repetition_iteration, ncol = 4)
    
    for(rep_iter in 1:repetition_iteration){
      
      y_data <- y_list[[rep_iter]]
      
      est_CP <- NBS_CPD(y_data)
      num_CP <- length(est_CP)
      
      gt_CP_corrected <- c(1, true_CP, 100) 
      est_CP_corrected <- c(1, est_CP, 100)
      
      gt_list <- est_list <- list();
      for(i in 2:length(gt_CP_corrected)){
        gt_list[[i-1]] <- gt_CP_corrected[i-1]:(gt_CP_corrected[i]-1)
      }
      for(i in 2:length(est_CP_corrected)){
        est_list[[i-1]] <- est_CP_corrected[i-1]:(est_CP_corrected[i]-1)
      }
      
      if(num_CP == 0){
        dist_est_gt <- Inf
        dist_gt_est <- -Inf
        covering_metric <- 0
      }else{
        
        holder <- c()
        for(i in true_CP){
          dist_diff <- c()
          for(j in est_CP){dist_diff <- c(dist_diff, abs(j-i))}
          holder <- c(holder, min(dist_diff))
        }
        dist_est_gt <- max(holder)
        
        holder <- c()
        for(i in est_CP){
          dist_diff <- c()
          for(j in true_CP){dist_diff <- c(dist_diff, abs(j-i))}
          holder <- c(holder, min(dist_diff))
        }
        dist_gt_est <- max(holder)
        
        covering_metric <- 0
        for(i in 1:length(gt_list)){
          A <- gt_list[[i]]
          jaccard <- c()
          for(j in 1:length(est_list)){
            A_prime <- est_list[[j]]
            jaccard <- c(jaccard,length(intersect(A,A_prime))/length(union(A,A_prime)))
          }
          covering_metric <- covering_metric + length(A)*max(jaccard)
        }
        covering_metric <- covering_metric/100
      }
      
      abs_error <- abs(num_CP - length(true_CP))
      
      result[rep_iter, 1] <- abs_error
      result[rep_iter, 2] <- dist_est_gt
      result[rep_iter, 3] <- dist_gt_est
      result[rep_iter, 4] <- covering_metric
    }
    
    return(result)
    
  }
  
}




################################
# Evaluation Function (CPDker) #
################################


Evaluation_CPDker <- function(y_list, num_stats, M, is_experiment=FALSE, true_CP=c(26, 51, 76)){
  
  library(ergm)
  
  if(is_experiment){
    
    y_data <- y_list
    num_time <- length(y_data)
    n <- dim(y_data[[1]])[1]
    tau <- num_time-1
    
    temp <- matrix(NA, nrow=num_time, ncol=num_stats)
    for(t in 1:num_time){
      if(num_stats == 2){
        yt <- network(y_data[[t]], directed = F)
        temp[t,] = summary(yt ~ edges + triangles)
      }else if(num_stats == 3){
        yt <- network(y_data[[t]], directed = F)
        temp[t,] = summary(yt ~ edges + isolates + triangles)
      }
    }
    Y <- as.matrix(t(temp)) # p by T
    
    p = num_stats
    intervals = WBS.intervals(M = M, lower = 1, upper = ncol(Y)) #Random intervals
    K_max = 5
    h = 5*(K_max*log(num_time)/num_time)^{1/p} # bandwith
    temp = WBS.multi.nonpar(Y, Y, 1, ncol(Y), intervals$Alpha, intervals$Beta, h, delta = 10)
    temp = thresholdBS(temp, median(temp$Dval))
    est_CP <- c(temp$cpt_hat$location)
    est_CP <- sort(est_CP)
    est_CP <- est_CP + 1
    return(est_CP)
    
  }else{
    
    repetition_iteration <- length(y_list)
    result <- matrix(NA, nrow=repetition_iteration, ncol = 4)
    
    for(rep_iter in 1:repetition_iteration){
      
      y_data <- y_list[[rep_iter]]
      num_time <- length(y_data)
      n <- dim(y_data[[1]])[1]
      tau <- num_time-1 # i.e. T=150, then num_theta=149 (tau), and num_beta=148
      
      temp <- matrix(NA, nrow=num_time, ncol=num_stats)
      for(t in 1:num_time){
        yt <- network(y_data[[t]], directed = T)
        if(num_stats == 2){
          temp[t,] = summary(yt ~ edges + mutual)
        }else if(num_stats == 3){
          temp[t,] = summary(yt ~ edges + mutual + triangles)
        }else if(num_stats == 4){
          yt <- network(yt)
          network::set.vertex.attribute(yt, "Gender", gender)
          temp[t,] = summary(yt ~ edges + mutual + triangles + nodematch("Gender"))
        }
      }
      Y <- as.matrix(t(temp)) # p by T
      
      p = num_stats
      intervals = WBS.intervals(M = M, lower = 1, upper = ncol(Y)) #Random intervals
      K_max = 5
      h = 5*(K_max*log(num_time)/num_time)^{1/p} # bandwith
      temp = WBS.multi.nonpar(Y, Y, 1, ncol(Y), intervals$Alpha, intervals$Beta, h, delta = 10)
      temp = thresholdBS(temp, median(temp$Dval))
      est_CP <- c(temp$cpt_hat$location)
      est_CP <- sort(est_CP)
      
      est_CP <- est_CP + 1
      num_CP <- length(est_CP)
      
      gt_CP_corrected <- c(1, true_CP, tau+2) # tau+2 = T+1
      est_CP_corrected <- c(1, est_CP, tau+2)
      
      gt_list <- est_list <- list();
      for(i in 2:length(gt_CP_corrected)){
        gt_list[[i-1]] <- gt_CP_corrected[i-1]:(gt_CP_corrected[i]-1)
      }
      for(i in 2:length(est_CP_corrected)){
        est_list[[i-1]] <- est_CP_corrected[i-1]:(est_CP_corrected[i]-1)
      }
      
      if(num_CP == 0){
        dist_est_gt <- Inf
        dist_gt_est <- -Inf
        covering_metric <- 0
      }else{
        
        holder <- c()
        for(i in true_CP){
          dist_diff <- c()
          for(j in est_CP){dist_diff <- c(dist_diff, abs(j-i))}
          holder <- c(holder, min(dist_diff))
        }
        dist_est_gt <- max(holder)
        
        holder <- c()
        for(i in est_CP){
          dist_diff <- c()
          for(j in true_CP){dist_diff <- c(dist_diff, abs(j-i))}
          holder <- c(holder, min(dist_diff))
        }
        dist_gt_est <- max(holder)
        
        covering_metric <- 0
        for(i in 1:length(gt_list)){
          A <- gt_list[[i]]
          jaccard <- c()
          for(j in 1:length(est_list)){
            A_prime <- est_list[[j]]
            jaccard <- c(jaccard,length(intersect(A,A_prime))/length(union(A,A_prime)))
          }
          covering_metric <- covering_metric + length(A)*max(jaccard)
        }
        covering_metric <- covering_metric/(tau+1) # tau+1 = T
        
      }
      
      abs_error <- abs(num_CP - length(true_CP))
      
      result[rep_iter, 1] <- abs_error
      result[rep_iter, 2] <- dist_est_gt
      result[rep_iter, 3] <- dist_gt_est
      result[rep_iter, 4] <- covering_metric
      
    }
    
    return(result)
  }
  
}



