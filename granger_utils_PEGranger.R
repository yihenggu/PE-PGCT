# Dec 27, 2025 # 
library(plm)
library(lmtest)
library(harmonicmeanp)

C <- function(T_adj){
  return(max(1.5, log(log(T_adj))))
}


PE_PGCT  <-  function(data, order = 1L, index = NULL){
  
  ## panel information
  pdim <- pdim(data)
  balanced <- pdim$balanced
  N  <- pdim$nT$n
  T_adj <- pdim$nT$T - order
  T_obs <- pdim$nT$T
  Ti <- pdim$Tint$Ti
  indi <- unclass(index(data))[[1L]]
  
  ## DH
  {
    DH <- pgrangertest(y ~ x, data = data, order = order)
    pval_DH <- DH$p.value
    stat_DH <- DH$statistic
  }
  
  ## HPJ
  {
    X <- list()
    Z <- list()
    Mz <- list()
    y <- NULL
    bhat_left <- matrix(0, nrow = order, ncol = order)
    bhat_right <- matrix(0, nrow = order, ncol = 1)
    bhat12_left <- matrix(0, nrow = order, ncol = order)
    bhat12_right <- matrix(0, nrow = order, ncol = 1)
    bhat21_left <- matrix(0, nrow = order, ncol = order)
    bhat21_right <- matrix(0, nrow = order, ncol = 1)
    
    for (i in 1:N){
      Xi <- NULL
      Zi <- NULL
      
      xi <- data$x[data$id == i]
      yi <- data$y[data$id == i]
      for (t in order:(T_obs-1)){
        xit <- xi[t:(t-order+1)]
        Xi <- cbind(Xi,xit)
        
        zit <- c(1,yi[t:(t-order+1)])
        Zi <-cbind(Zi,zit)
      }
      
      #for bhat
      X[[i]] <- t(Xi)
      Z[[i]] <- t(Zi)
      y <-cbind(y,yi[-(1:order)])
      Mzi <- diag(T_adj) - Z[[i]] %*% solve(t(Z[[i]]) %*% Z[[i]]) %*% t(Z[[i]])
      Mz[[i]] <- Mzi
      bhat_left <- bhat_left + t(X[[i]]) %*% Mz[[i]] %*% X[[i]]
      bhat_right <- bhat_right + t(X[[i]]) %*% Mz[[i]] %*% y[,i]
      
      
      #for bhat12
      X12i <- NULL
      Z12i <- NULL
      for (t1 in order: (order+T_adj/2-1)){
        x12it <- xi[t1:(t1-order+1)]
        X12i <- cbind(X12i,x12it)
        z12it <- c(1,yi[t1:(t1-order+1)])
        Z12i <- cbind(Z12i, z12it)
        
      }
      X12i <- t(X12i)
      Z12i <- t(Z12i)
      y12i <- yi[(order+1):(T_adj/2 + order)]
      Mz12i <- diag(T_adj/2) - Z12i %*% solve(t(Z12i) %*% Z12i) %*% t(Z12i)
      bhat12_left <- bhat12_left + t(X12i) %*% Mz12i %*% X12i
      bhat12_right <- bhat12_right + t(X12i) %*% Mz12i %*% y12i
      
      
      #calculate bhat21
      X21i <- NULL
      Z21i <- NULL
      for (t2 in (T_adj/2 + order) : (T_obs-1)){
        x21it <- xi[t2:(t2-order+1)]
        X21i <- cbind(X21i,x21it)
        z21it <- c(1,yi[t2:(t2-order+1)])
        Z21i <- cbind(Z21i, z21it)
      }
      X21i <- t(X21i)
      Z21i <- t(Z21i)
      y21i <- yi[(order+1+T_adj/2):(T_adj+order)]
      Mz21i <- diag(T_adj- T_adj/2) - Z21i %*% solve(t(Z21i) %*% Z21i) %*% t(Z21i)
      bhat21_left <- bhat21_left + t(X21i) %*% Mz21i %*% X21i
      bhat21_right <- bhat21_right + t(X21i) %*% Mz21i %*% y21i
    }
    
    #calculate bhat
    bhat <- solve(bhat_left) %*% bhat_right
    #calculate bhat12, bhat21
    bhat12 <- solve(bhat12_left) %*% bhat12_right
    bhat21 <- solve(bhat21_left) %*% bhat21_right
    
    #calculate btilde
    btilde <- 2*bhat - (bhat12+bhat21)/2
    
    Jhat <- bhat_left / (N * T_adj)
    
    sigma2hat <- 0
    Vhatcs <- 0
    
    for (i in 1:N){
      sigma2hat <- sigma2hat + t(y[,i]-X[[i]] %*% bhat) %*% Mz[[i]] %*% (y[,i]-X[[i]] %*% bhat)
      Vhatcs <- Vhatcs + t(X[[i]]) %*% Mz[[i]] %*% (y[,i]-X[[i]]%*% bhat) %*% t(y[,i]-X[[i]]%*% bhat) %*% Mz[[i]] %*% X[[i]]
    }
    sigma2hat <- sigma2hat / (N*(T_adj-1-order)-order)
    Vhatcs <- Vhatcs / (N*(T_adj-1-order)-order)
    stat_HPJ <- N * T_adj * t(btilde) %*% solve(solve(Jhat) %*% Vhatcs %*% solve(Jhat)) %*% btilde
    pval_HPJ <- 1 - pchisq(stat_HPJ, order)
  }
  
  
  {
    listdata <- split(data[, c("time","y","x")], data$id) # split data per individual
    #listdata <- collapse::rsplit(data, indi, use.names = FALSE)
    ## use lmtest::grangertest for the individual Granger tests
    
    # for this, if necessary, expand order argument for lmtest::grangertest to full length (N)
    # [but leave variable 'order' in its current length for later decision making]
    order_grangertest <- if(length(order) == 1L) rep(order, N) else order
    
    # Dumitrescu/Hurlin (2012), p. 1453 use the Chisq definition of the Granger test
    grangertests_i <- mapply(function(data, order)
      lmtest::grangertest(y~x, data = data,
                          order = order, test = "F"),
      listdata, order_grangertest, SIMPLIFY = FALSE)
    
    # extract Wald/Chisq-statistics and p-values of individual Granger tests
    Wi   <- vapply(grangertests_i, function(g) g[["F"]][2L],        FUN.VALUE = 0.0, USE.NAMES = FALSE)
    pWi  <- vapply(grangertests_i, function(g) g[["Pr(>F)"]][[2L]], FUN.VALUE = 0.0, USE.NAMES = FALSE)
    # mean(Wi) == pgrangertest(y ~ x, data = data, order = order, test="Wbar")
    
    stat_Cauchy <- sum(tan((0.5-pWi)*pi))/length(pWi)
    pval_Cauchy <- 1-pcauchy(stat_Cauchy)
    
  }
  
  
  ## Power Enhancement Component
  {
    max_Wi <- max(Wi)
    PE_delta <- (order + 2 * sqrt(order * log(N)) + 2 * log(N)) * C(T_adj)
    
    J1 <- ifelse(max_Wi > PE_delta, sqrt(N) * max_Wi, 0)
    
    stat_PE_DH <- abs(stat_DH) + J1
    pval_PE_DH <- 2*pnorm(stat_PE_DH, lower.tail=FALSE)
    
    stat_PE_HPJ <- stat_HPJ + J1
    pval_PE_HPJ <- 1-pchisq(stat_PE_HPJ, order)
    
    if(T_adj > N)
    {
      stat_PE <- stat_PE_DH
      pval_PE <- pval_PE_DH
    }else{
      stat_PE <- stat_PE_HPJ
      pval_PE <- pval_PE_HPJ
    }
  }
  
  return(list(pval_PE = pval_PE, 
              stat_PE = stat_PE,
              J_PE = J1,
              PE_delta = PE_delta,
              max_Wi = max_Wi,
              Wi = Wi,
              pWi = pWi
  ))
}
