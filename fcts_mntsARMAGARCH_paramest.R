library(quadprog)
library(quantmod)
library(functional)
library(nloptr)
library(pracma)
library(spatstat)
library(Matrix)
library(foreach)
library(doParallel)
library(xtable)
library(rugarch)
library(temStaR)
source("readStockData.R")
source("miscellaneous_tools.R")


paramest_ind_armagarch_parallel <- function(ret, numofelements, yTickers=NULL, parallelSocketCluster = NULL ){
  flag_sock  <- FALSE
  if (is.null(parallelSocketCluster)){
    numofcluster  <- detectCores()
    parallelSocketCluster <- makePSOCKcluster(numofcluster)
    registerDoParallel(parallelSocketCluster)
    flag_sock <- TRUE
  }
  
  muvec <- matrix(nrow = numofelements, ncol = 1)
  sigvec <- matrix(nrow = numofelements, ncol = 1)
  ntsparamMtx <- matrix(nrow = numofelements, ncol = 5)
  colnames(ntsparamMtx) <- c("mu", "sigma", "alpha", "theta", "beta")
  rownames(ntsparamMtx) <- yTickers
  stdret <- matrix(nrow = dim(ret)[1], ncol = numofelements)
  tgarchparam <- matrix(nrow = numofelements, ncol = 7)
  rownames(tgarchparam) <- yTickers
  
  res <- foreach( n = 1:numofelements) %dopar% {
    library(rugarch)
    library(temStaR)
    tspec <- ugarchspec(mean.model = list(armaOrder = c(1,1)),
                        variance.model = list(garchOrder = c(1,1), 
                                              model = "sGARCH"), 
                        distribution.model = "std")
    tf <- ugarchfit(spec = tspec, data = ret[,n])
    tgaparam <- coef(tf)
    v1 <- tf@fit$sigma
    resi <- tf@fit$residuals/v1 #standard residual
    tforc <- ugarchforecast(tf, n.ahead=1)
    mu <- tforc@forecast$seriesFor[1]
    sig <- tforc@forecast$sigmaFor[1]
    stdntsparam <- fitstdnts(resi)
    return(list(mu = mu, sig = sig, tgaparam = tgaparam, 
                al = stdntsparam["alpha"], th = stdntsparam["theta"], beta = stdntsparam["beta"], 
                resi = resi))
  }
  
  for( n in 1:numofelements){
    tgarchparam[n,] <- res[[n]]$tgaparam
    stdret[,n] <- res[[n]]$resi
    ntsparamMtx[n,"mu"] <- res[[n]]$mu
    ntsparamMtx[n,"sigma"] <- res[[n]]$sig
    ntsparamMtx[n,"alpha"] <- res[[n]]$al 
    ntsparamMtx[n,"theta"] <- res[[n]]$th
    ntsparamMtx[n,"beta"] <- res[[n]]$beta
  }
  corrMtx <- cor(stdret)
  colnames(tgarchparam) <- names(res[[1]]$tgaparam)
  colnames(corrMtx) <- yTickers
  rownames(corrMtx) <- yTickers

  if(flag_sock){
    stopCluster(parallelSocketCluster)
  }
  return(list(tagparam=tgarchparam, ntsparam=ntsparamMtx, corrMtx = corrMtx))
}


paramest_mntsarmagarch <- function(ret, numofelements, yTickers=NULL){
  muvec <- matrix(nrow = numofelements, ncol = 1)
  sigvec <- matrix(nrow = numofelements, ncol = 1)
  stdret <- matrix(nrow = dim(ret)[1], ncol = numofelements)
  tgarchparam <- matrix(nrow = numofelements, ncol = 7)
  rownames(tgarchparam) <- yTickers
  
  tspec <- ugarchspec(mean.model = list(armaOrder = c(1,1)),
                      variance.model = list(garchOrder = c(1,1), 
                                            model = "sGARCH"), 
                      distribution.model = "std")
  
  for( n in 1:numofelements) {
    tf <- ugarchfit(spec = tspec, data = ret[,n])
    a<-coef(tf)
    if (n==1){
      colnames(tgarchparam) <- names(a)
    }
    tgarchparam[n,] <- a
    v1 <- tf@fit$sigma
    stdret[,n] <- tf@fit$residuals/v1 #standard residual
    tforc <- ugarchforecast(tf, n.ahead=1)
    muvec[n] <- tforc@forecast$seriesFor[1]
    sigvec[n] <- tforc@forecast$sigmaFor[1]
  }
  
  options(warn=-1)
  st <- fitmnts(returndata = stdret, n = numofelements, stdflag = TRUE)
  st$mu <- as.numeric(muvec)
  st$sigma <- as.numeric(sigvec)
  st$CovMtx <- diag(as.numeric(sigvec))%*%st$CovMtx%*%diag(as.numeric(sigvec))
  options(warn=0)  
  
  colnames(tgarchparam) <- names(res[[1]]$tgaparam)
  colnames(st$Rho) <- yTickers
  rownames(st$Rho) <- yTickers
  colnames(st$CovMtx) <- yTickers
  rownames(st$CovMtx) <- yTickers
  names(st$mu) <- yTickers
  names(st$sigma) <- yTickers
  names(st$beta) <- yTickers
  
  return(list(st=st, tagparam=tgarchparam))
}




paramest_mntsarmagarch_parallel <- function(ret, numofelements, yTickers=NULL, parallelSocketCluster = NULL ){
  flag_sock  <- FALSE
  if (is.null(parallelSocketCluster)){
    numofcluster  <- detectCores()
    parallelSocketCluster <- makePSOCKcluster(numofcluster)
    registerDoParallel(parallelSocketCluster)
    flag_sock <- TRUE
  }
  
  muvec <- matrix(nrow = numofelements, ncol = 1)
  sigvec <- matrix(nrow = numofelements, ncol = 1)
  stdret <- matrix(nrow = dim(ret)[1], ncol = numofelements)
  tgarchparam <- matrix(nrow = numofelements, ncol = 7)
  rownames(tgarchparam) <- yTickers
  
  res <- foreach( n = 1:numofelements) %dopar% {
    library(rugarch)
    tspec <- ugarchspec(mean.model = list(armaOrder = c(1,1)),
                        variance.model = list(garchOrder = c(1,1), 
                                              model = "sGARCH"), 
                        distribution.model = "std")
    tf <- ugarchfit(spec = tspec, data = ret[,n])
    tgaparam <- coef(tf)
    v1 <- tf@fit$sigma
    resi <- tf@fit$residuals/v1 #standard residual
    tforc <- ugarchforecast(tf, n.ahead=1)
    mu <- tforc@forecast$seriesFor[1]
    sig <- tforc@forecast$sigmaFor[1]
    return(list(mu = mu, sig = sig, tgaparam = tgaparam, resi = resi))
  }
  
  for( n in 1:numofelements){
    tgarchparam[n,] <- res[[n]]$tgaparam
    stdret[,n] <- res[[n]]$resi
    muvec[n] <- res[[n]]$mu
    sigvec[n] <- res[[n]]$sig
  }
  options(warn=-1)
  st <- fitmnts_par(returndata = stdret, n = numofelements, stdflag = TRUE)
  st$mu <- as.numeric(muvec)
  st$sigma <- as.numeric(sigvec)
  st$CovMtx <- diag(as.numeric(sigvec))%*%st$CovMtx%*%diag(as.numeric(sigvec))
  options(warn=0)  
  
  colnames(tgarchparam) <- names(res[[1]]$tgaparam)
  colnames(st$Rho) <- yTickers
  rownames(st$Rho) <- yTickers
  colnames(st$CovMtx) <- yTickers
  rownames(st$CovMtx) <- yTickers
  names(st$mu) <- yTickers
  names(st$sigma) <- yTickers
  names(st$beta) <- yTickers
  
  if(flag_sock){
    stopCluster(parallelSocketCluster)
  }
  return(list(st=st, tagparam=tgarchparam))
}

errfuncSkewMatch <- function(x, alth, indstdparam){
  s1 <- moments_stdNTS(indstdparam)
  s2 <- moments_stdNTS(c(alth[1], alth[2], x ))
  e <- sqrt((s1[3]-s2[3])^2)
  return(e)
}


indStdNTS2StdMNTS <-function(alphavec, thetavec, betavec, muvec, sigmavec, corrMtx, yTickers = NULL, alphaNtheta = NULL){
  if( is.null(alphaNtheta) ){
    al <- mean(alphavec)
    th <- mean(thetavec)
  } else {
    al <- alphaNtheta[1]
    th <- alphaNtheta[2]
  }
  
  N <- length(betavec)
  newbeta <- array(dim = N)
  for( n in 1:N ){
    opt <- optimize(f = errfuncSkewMatch,
                    interval = c( -sqrt(2*th/(2-al)), 
                                  sqrt(2*th/(2-al)) ),
                    alth = c(al, th), 
                    indstdparam = c(alphavec[n], thetavec[n], betavec[n])
    )
    newbeta[n] <- opt$minimum
  }
  
  gamma <- sqrt(1-newbeta^2*(2-al)/(2*th))
  bmat <- matrix(newbeta, nrow = 1, ncol = N)
  bmat <- t(bmat)%*%bmat
  bmat = ((2-al)/(2*th))*bmat
  Rho <- diag(1/gamma)%*%(corrMtx-bmat)%*%diag(1/gamma)
  #x <- nearPD(Rho)$mat
  Rho <- as.matrix(nearPD(Rho)$mat)
  covmatrix <- diag(as.numeric(sigmavec))%*%corrMtx%*%diag(as.numeric(sigmavec))
  stMNTS = list(mu = muvec, sigma = sigmavec, 
                alpha = al, theta = th,
                beta = newbeta, 
                Rho = Rho,
                CovMtx =  covmatrix )
  if( !is.null(yTickers) ){
    names(stMNTS$mu) <- yTickers
    names(stMNTS$sigma) <- yTickers
    names(stMNTS$beta) <- yTickers 
    colnames(stMNTS$Rho) <- yTickers
    rownames(stMNTS$Rho) <- yTickers
    colnames(stMNTS$CovMtx) <- yTickers
    rownames(stMNTS$CovMtx) <- yTickers
  }
  return(stMNTS)
}
