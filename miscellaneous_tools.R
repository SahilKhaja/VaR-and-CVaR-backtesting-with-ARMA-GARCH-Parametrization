library(nloptr)
library(rugarch)
library(pracma)
library(functional)
library(spatstat)
library(temStaR)
source("readHistoryStockDataFromYahoo.R")


get_dailyreturn <- function(symb, begindate, enddate){
  p <- readHistoryStockDataFromYahoo(symb, as.Date(begindate), as.Date(enddate))
  datevec <- p$Date
  n <- length(p$Adj.Close)
  ret <- c(0, as.numeric(diff(log(p$Adj.Close))))
  return( list( datearray = datevec, pricearray = p$Adj.Close, dailyretarray = ret ) )
}

CVaR_normal <- function( eps, mu=0, sigma=1, dt = 1){
  var <- -qnorm(eps, mean = mu*dt, sd = sigma*sqrt(dt));
  stdK <- qnorm(eps);
  avar <- sigma*sqrt(dt)/(eps*sqrt(2*pi))*exp(-stdK^2/2)-mu*dt;
  
  return(avar)
}

CVaR_T <- function( ep, nu = 5 ){
  var = -qt( p = ep, df = nu)
  if (nu>1)
    avar = gamma((nu+1)/2)/gamma(nu/2)*sqrt(nu)/((nu-1)*ep*sqrt(pi))*(1+(var)^2/nu)^((1-nu)/2)
  else
    avar = 0;
  
  return(avar)
}


emp_CVaR_ret <- function(alpha,ret){
  eps <- 1-alpha
  lossdata <- -ret
  n <- length(lossdata)
  sortlossdata <- sort(lossdata)
  CVaR = array(dim = length(eps), dimnames = list(alpha))
  for (j in 1:length(eps)){
    cce <- ceiling(n*eps[j])
    CVaR[j] <- (sum(sortlossdata[cce:n])/n+((cce-1)/n-eps[j])*sortlossdata[cce-1])/(1-eps[j]);
  }
  return(CVaR)
}