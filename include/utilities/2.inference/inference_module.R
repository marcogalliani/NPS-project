# (0) Dependecies --------------------------------------------------
library(pbapply)
library(purrr)
library(parallel)


# (1) Permutational Inference --------------------------------------------------
perm_t.test <- function(data1, data2, test_stat, B, seed, cl){
  
  set.seed(seed)
  
  #'test_stat -> function taking as input the two permuted vectors and returning
  #'              a double
  #'data1 and data2 are the data about the two groups that we want to compare
  #'cl is the result of a makeCluster operation
  
  test_results <- list()
  
  #(1) Observed Test statistic
  test_results$obs.T_stat <- test_stat(data1, data2)
  
  n1 <- nrow(data1)
  n2 <- nrow(data2)
  
  n <- n1 + n2

  clusterExport(cl=cl, list('data1', 'data2', "n1"), envir = environment())
  
  test_stat_eval <- function(permutation){
    data.pooled <- rbind(data1, data2)
    
    perm.pooled.data <- data.pooled[permutation, , drop=F]
    
    data1.perm <- perm.pooled.data[1:n1, , drop=F]
    data2.perm <- perm.pooled.data[(n1+1):n, , drop=F]
    
    return(test_stat(data1.perm, data2.perm))
  }
  
  #(2) Permutational distribution
  test_results$perm_distr.T_stat <- 
    pbreplicate(B, test_stat_eval(sample(n)), cl=cl)
  
  #(3) P-value of the test
  test_results$p_val <- sum(test_results$perm_distr.T_stat >= test_results$obs.T_stat)/B
  
  return(test_results)
}

perm_ctr_sym.test <- function(data, mu0, test_stat, B, seed, cl){
  
  set.seed(seed)
  
  #'test_stat: function taking as input the two permuted vectors and returning a double
  #'data: data that we want to test
  #'mu0: assumption on the center of symmetry
  #'cl: the result of a makeCluster operation
  
  test_results <- list()
  
  #(1) Observed Test statistic
  test_results$obs.T_stat <- test_stat(data, mu0)
  
  n <- nrow(data)
  p <- ncol(data)

  clusterExport(cl=cl, list("data","n","p","mu0"), envir = environment())
  
  test_stat_eval <- function(signs.perm){
    #permuted sample
    data.perm <- matrix(mu0,nrow=n,ncol=p,byrow=T) + (data - mu0) * matrix(signs.perm,nrow=n,ncol=p,byrow=F)
    
    return(test_stat(data.perm, mu0))
  }
  
  #(2) Permutational distribution
  test_results$perm_distr.T_stat <- 
    pbreplicate(B, test_stat_eval(rbinom(n, 1, 0.5)*2 - 1), cl=cl)
  
  #(3) P-value of the test
  test_results$p_val <- sum(test_results$perm_distr.T_stat >= test_results$obs.T_stat)/B
  
  return(test_results)
}

perm_paired_t.test <-function(data1, data2, test_stat, B, seed, cl){
  
  set.seed(seed)
  
  #'test_stat -> function taking as input the two permuted vectors and returning
  #'              a double
  #'data1 and data2 are the data about the two groups that we want to compare
  #'cl is the result of a makeCluster operation
  
  test_results <- list()
  
  #(1) Observed Test statistic
  test_results$obs.T_stat <- test_stat(data1, data2)
  
  # here naturally I can use data1 or data2 (paired samples)
  n <- dim(data1)[1]
  p <- dim(data2)[2] 
  
  data.pooled <- rbind(data1,data2)
  
  clusterExport(cl=cl, list('data.pooled', "n","p"), envir = environment())
  
  test_stat_eval <- function(perm.indices.data1){
    data1.perm <- data.pooled[perm.indices.data1, , drop=F]
    data2.perm <- data.pooled[-perm.indices.data1, , drop=F]
    
    return(test_stat(data1.perm, data2.perm))
  }
  
  #(2) Permutational distribution
  test_results$perm_distr.T_stat <- 
    pbreplicate(B, test_stat_eval(seq(1, n) + n * rbinom(n,1, 0.5)), cl=cl)
  
  #(3) P-value of the test
  test_results$p_val <- sum(test_results$perm_distr.T_stat >= test_results$obs.T_stat)/B
  
  return(test_results)
}


perm_anova.test <- function(fit.H0, fit.H1,data,test_stat, B, seed, cl){
  set.seed(seed)
  
  #'fitH0: result of fitting the model under H0 assumptions
  #'fitH1: result of fitting the model under H1 assumptions
  #'data: data on which the two models were fitted
  #'cl is the result of a makeCluster operation
  
  test_results <- list()
  
  #(1) Observed Test statistic
  test_results$obs.T_stat <- test_stat(fit.H0,fit.H1)
  
  #(2) Permutational distribution
  n <- nrow(data)
  target_name <- all.vars(terms(fit.H0))[1]

  clusterExport(cl=cl, list("data","fit.H0","fit.H1","n","target_name"), envir = environment())
  
  test_stat_eval <- function(permutation){
    # permuting the residuals
    perm.target <- fitted(fit.H0) + residuals(fit.H0)[permutation]
    data[[target_name]] <- perm.target
    
    # fitting on the permuted data
    perm.fit.H1 <- update(fit.H1, data=data)
    
    return(test_stat(fit.H0, perm.fit.H1))
  }
  
  test_results$perm_distr.T_stat <- 
    pbreplicate(B, test_stat_eval(sample(n)), cl=cl)
  
  #(3) P-value of the test
  test_results$p_val <- sum(test_results$perm_distr.T_stat >= test_results$obs.T_stat)/B
  
  return(test_results)
  
}

plot.perm_results <- function(perm_results){
  #p-value
  cat("P-value: ", perm_results$p_val,"\n")
  
  #histogram of the permutational distribution
  hist(perm_results$perm_distr.T_stat, 
       xlim=range(c(perm_results$perm_distr.T_stat,perm_results$obs.T_stat)))
  abline(v=perm_results$obs.T_stat,col=3,lwd=4)
  
  #empirical distribution function
  plot(ecdf(perm_results$perm_distr.T_stat))
  abline(v=perm_results$obs.T_stat,col=3,lwd=4)
}

# (2) Bootstrap inference --------------------------------------------------

boot.test.adapter <- function(boot_distr.T_stat,obs.T_stat,mu0){
  alpha_grid <- seq(0.001, 1, by=0.001)
  
  #computation of the confidence intervals
  CI.computation <- function(alpha_level){
    right.quantile <- quantile(boot_distr.T_stat, 1 - alpha_level/2)
    left.quantile  <- quantile(boot_distr.T_stat, alpha_level/2)
    
    out <- c(obs.T_stat - (right.quantile - obs.T_stat), 
             obs.T_stat - (left.quantile - obs.T_stat)
    )
    names(out) <- c('lwr','upr')
    return(out)
  }
  
  CI.list <- pblapply(alpha_grid, CI.computation)
  CI.mat <- dplyr::bind_rows(CI.list)
  
  #computing the p-value
  p_val <- (alpha_grid[CI.mat[,1]>mu0 | CI.mat[,2]<mu0])[1]
  
  return(p_val)
}

boot.estimation <- function(data, t_stat, alpha, B, seed, cl){
  
  #'test_stat: statistic of which we want to compute the distribution
  #'data: data on which we want to estimate the distribution of the test statistic
  #'alpha: confidence level to build confidence intervals
  #'cl: the result of a makeCluster operation
  
  set.seed(seed)
  estimation_results <- list()
  
  data <- as.data.frame(data)
  n <- dim(data)[1]
  
  #(1) Observed test statistic
  estimation_results$obs.T_stat <- t_stat(data)
  estimation_results$sd <- apply(data, MARGIN=2, FUN=sd)
  
  #(2) Bootstrap distribution
  cl <- makeCluster(parallel::detectCores())
  clusterExport(cl=cl, list('data','n','t_stat'), envir = environment())
  
  test_stat_eval <- function(boot.indices){
    return(t_stat(data[boot.indices,,drop=F]))
  }
  
  estimation_results$boot_distr.T_stat <- 
    pbreplicate(B,
                test_stat_eval(sample(1:n, replace=T)),
                cl=cl, simplify = "matrix")
  
  #(3) Bias-Variance
  estimation_results$bias <- mean(estimation_results$boot_distr.T_stat) - estimation_results$obs.T_stat
  estimation_results$variance <- var(estimation_results$boot_distr.T_stat)
  
  estimation_results$MSE <- estimation_results$bias^2 + estimation_results$variance
  
  #(4) Confidence intervals (-> reverse percentile intervals)
  left.quantile <- apply(estimation_results$boot_distr.T_stat, 
                         MARGIN=1, FUN=quantile, probs=alpha/2) 
  right.quantile  <- apply(estimation_results$boot_distr.T_stat, 
                           MARGIN=1, FUN=quantile, probs=1-alpha/2) 
  
  estimation_results$CI.RP <- 
    rbind(estimation_results$obs.T_stat - (right.quantile - estimation_results$obs.T_stat), 
          estimation_results$obs.T_stat - (left.quantile - estimation_results$obs.T_stat))
  
  return(estimation_results)
}


boot.t.test <- function(data1, data2, test_stat, alpha, B, seed, cl){
  
  test_results <- list()
  
  n1 <- nrow(data1)
  n2 <- nrow(data2)
  n <- n1 + n2
  
  #(1) Observed Test statistic
  test_results$obs.T_stat <- test_stat(data1, data2)
  
  #(2) Bootstrap distribution
  clusterExport(cl=cl, list('data1', 'data2', "n1"), envir = environment())
  
  test_stat_eval <- function(ind1.b, ind2.b){
    data1.boot <- data1[ind1.b, , drop=F]
    data2.boot <- data2[ind2.b, , drop=F]
    
    return(test_stat(data1.boot, data2.boot))
  }
  
  test_results$boot_distr.T_stat <- 
    pbreplicate(B, test_stat_eval(sample(1:n1, replace = T),
                                  sample(1:n2, replace = T)),
                cl=cl)
  
  #(3) Confidence intervals -> reverse percentile intervals
  right.quantile <- quantile(test_results$boot_distr.T_stat, 1 - alpha/2)
  left.quantile  <- quantile(test_results$boot_distr.T_stat, alpha/2)
  test_results$CI.RP <- c(test_results$obs.T_stat - (right.quantile - test_results$obs.T_stat), 
                          test_results$obs.T_stat - (left.quantile - test_results$obs.T_stat))

  
  #(4) p-value of the test
  test_results$p_val <- boot.test.adapter(test_results$boot_distr.T_stat,
                                          test_results$obs.T_stat,
                                          0)
  
  return(test_results)
}

boot.regression.inference <- function(fit,data,test_stat,alpha,B,seed,cl){
  test_results <- list()
  
  #(1) Observed Test statistic
  test_results$obs.T_stat <- test_stat(fit)
  
  n <- nrow(data)
  target_name <- all.vars(terms(fit))[1]
  
  #(2) Bootstrap distribution
  clusterExport(cl=cl, list("data","fit","n","target_name"), envir = environment())
  
  test_stat_eval <- function(boot.ind){
    # permuting the residuals
    boot.target <- fitted(fit) + residuals(fit)[boot.ind]
    data[[target_name]] <- boot.target
    
    # fitting on the permuted data
    boot.fit <- update(fit, data=data)
    
    return(test_stat(boot.fit))
  }
  test_results$boot_distr.T_stat <- 
    pbreplicate(B, test_stat_eval(sample(n,replace = T)), cl=cl)
  
  #(3) Bias-Variance
  test_results$bias <- mean(test_results$boot_distr.T_stat) - test_results$obs.T_stat
  test_results$variance <- var(test_results$boot_distr.T_stat)
  
  test_results$MSE <- test_results$bias^2 + test_results$variance
  
  #(4) Confidence intervals -> reverse percentile
  right.quantile <- quantile(test_results$boot_distr.T_stat, 1 - alpha/2)
  left.quantile  <- quantile(test_results$boot_distr.T_stat, alpha/2)
  test_results$CI.RP <- c(test_results$obs.T_stat - (right.quantile - test_results$obs.T_stat), 
                          test_results$obs.T_stat - (left.quantile - test_results$obs.T_stat))
  
  return(test_results)
}

boot.anova <- function(fit.H0, fit.H1,data,test_stat,alpha,B,seed,cl){
  test_results <- list()
  
  #(1) Observed Test statistic
  test_results$obs.T_stat <- test_stat(fit.H0,fit.H1)
  
  #(2) Bootstrap distribution
  n <- nrow(data)
  target_name <- all.vars(terms(fit.H0))[1]
  
  clusterExport(cl=cl, list("data","fit.H0","fit.H1","n","target_name"), envir = environment())
  
  test_stat_eval <- function(boot.ind){
    # permuting the residuals
    boot.target <- fitted(fit.H0) + residuals(fit.H0)[boot.ind]
    data[[target_name]] <- boot.target
    
    # fitting on the permuted data
    boot.fit.H1 <- update(fit.H1, data=data)
    
    return(test_stat(fit.H0, boot.fit.H1))
  }
  test_results$boot_distr.T_stat <- 
    pbreplicate(B, test_stat_eval(sample(n,replace = T)), cl=cl)
  
  #(3) Bias-Variance
  test_results$bias <- mean(test_results$boot_distr.T_stat) - test_results$obs.T_stat
  test_results$variance <- var(test_results$boot_distr.T_stat)
  
  test_results$MSE <- test_results$bias^2 + test_results$variance
  
  #(4) Confidence intervals -> reverse percentile
  right.quantile <- quantile(test_results$boot_distr.T_stat, 1 - alpha/2)
  left.quantile  <- quantile(test_results$boot_distr.T_stat, alpha/2)
  test_results$CI.RP <- c(test_results$obs.T_stat - (right.quantile - test_results$obs.T_stat), 
                          test_results$obs.T_stat - (left.quantile - test_results$obs.T_stat))
  
  return(test_results)
}

boot.coxph.inference <- function(fit, data, test_stat,alpha,B,seed,cl){
  test_results <- list()
  
  #(1) Observed Test statistic
  test_results$obs.T_stat <- test_stat(fit)
  
  n <- nrow(data)
  target_name <- all.vars(terms(fit))[1]
  
  #(2) Bootstrap distribution
  clusterExport(cl=cl, list("data","fit","n"), envir = environment())
  
  test_stat_eval <- function(boot.ind){
    library(survival)
    
    # fitting on the permuted data
    boot.fit <- update(fit, data=data[boot.ind,])
    
    return(test_stat(boot.fit))
  }
  test_results$boot_distr.T_stat <- 
    pbreplicate(B, test_stat_eval(sample(n,replace = T)), cl=cl)
  
  #(3) Bias-Variance
  test_results$bias <- mean(test_results$boot_distr.T_stat) - test_results$obs.T_stat
  test_results$variance <- var(test_results$boot_distr.T_stat)
  
  test_results$MSE <- test_results$bias^2 + test_results$variance
  
  #(4) Confidence intervals -> reverse percentile
  right.quantile <- quantile(test_results$boot_distr.T_stat, 1 - alpha/2)
  left.quantile  <- quantile(test_results$boot_distr.T_stat, alpha/2)
  test_results$CI.RP <- c(test_results$obs.T_stat - (right.quantile - test_results$obs.T_stat), 
                          test_results$obs.T_stat - (left.quantile - test_results$obs.T_stat))
  
  return(test_results)
}
