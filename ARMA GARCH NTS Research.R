# AMS 518 - Final Project


# Topic: Application of Tempered Stable Distribution in Risk Management. Comparative analysis of risk measures across different distributions


library(foreach)
library(doParallel)
library(xtable)
library(rugarch)
library(temStaR)
library(quadprog)
library(quantmod)
library(functional)
library(nloptr)
library(pracma)
library(NlcOptim)
library(evmix)
library(spatstat)
library(Matrix)
library(ggplot2)
library(MASS)
library(dplyr)
library(forecast)


# Set my working directory 
setwd("C:/Users/sahil/Downloads/MNTS_Optimization/MNTS_Optimization")
source("C:/Users/sahil/Downloads/MNTS_Optimization/MNTS_Optimization/readStockData.R")
source("C:/Users/sahil/Downloads/MNTS_Optimization/MNTS_Optimization/func_findOptport.R")
source("C:/Users/sahil/Downloads/MNTS_Optimization/MNTS_Optimization/miscellaneous_tools.R")
source("C:/Users/sahil/Downloads/MNTS_Optimization/MNTS_Optimization/fcts_mntsARMAGARCH_paramest.R")





# Pulling S&P500 data 
begindate <- as.Date("2000-03-01")
enddate<-as.Date("2014-12-31")
intervalDateNumber <- 1

yTickers<-c("SPY")
numofelements<-length(yTickers)
ret <- readdata(yTickers, numofelements, begindate, enddate, intervalDateNumber) # This is the log returns of the SPY



# Defining Helper functions
length(ret)
roll_var_norm <-function(r,alpha=0.01,window=500){
  n<-length(r); VaR <-rep(NA_real_,n)
  for(t in (window+1):n){
    x<-r[(t-window):(t-1)]
    m<-mean(x);s<-sd(x)
    VaR[t]<-m+s*qnorm(alpha)
  }
  return(VaR)
}

# roll_var_t <- function(r, alpha = 0.01, window = 500) {
#   n <- length(r)
#   VaR <- rep(NA_real_, n)
# 
#   for (t in (window + 1):n) {
#     x <- r[(t - window):(t - 1)]
#     fit <- try(
#       suppressWarnings(MASS::fitdistr(x, densfun = "t")),
#       silent = TRUE
#     )
# 
#     if (inherits(fit, "try-error")) {
#       m  <- mean(x)
#       s  <- sd(x)
#       nu <- 6
#     } else {
#       m  <- fit$estimate["m"]
#       s  <- fit$estimate["s"]
#       nu <- fit$estimate["df"]
#     }
# 
#     # Compute 1-step-ahead parametric VaR
#     VaR[t] <- m + s * qt(alpha, df = nu)
#   }
# 
#   return(VaR)
# }
  
# Trying an alternative code for fitting t distribution with 6 degrees of freedom

roll_var_t <- function(r, alpha = 0.01, window = 500, nu = nu) {
  n <- length(r)
  VaR <- rep(NA_real_, n)

  pb <- txtProgressBar(min = window + 1, max = n, style = 3)

  for (t in (window + 1):n) {
    
    # We are considering data from t-1 to create parameters for fitting VaR at time t
    x <- r[(t - window):(t - 1)]

    # Fit only location and scale, ignore estimated df
    fit <- MASS::fitdistr(x, densfun = "t")
    m <- fit$estimate["m"]
    s <- fit$estimate["s"]
    nu <- nu

    VaR[t] <- m + s * qt(alpha, df = nu)

    setTxtProgressBar(pb, t)
  }
  close(pb)
  return(VaR)
}

roll_var_arma_garch_norm <- function(r,alpha=0.01,window=500){
  n<-length(r)
  VaR<-rep(NA_real_,n)
  
  tspec<-rugarch::ugarchspec(
    variance.model=list(model="sGARCH",garchOrder=c(1,1)),
    mean.model = list(armaOrder = c(1,1)),
    distribution.model = "norm"
  )
  
  pb<-txtProgressBar(min=window+1,max=n,style=3)
  
  for (t in (window+1):n){
    
    x <- r[(t - window):(t - 1)]
    fit <- rugarch::ugarchfit(spec = tspec, data = x,solver="hybrid")
    fcast<-rugarch::ugarchforecast(fit, n.ahead = 1)
    mu <- fcast@forecast$seriesFor[1]
    sigma <- fcast@forecast$sigmaFor[1]
    VaR[t] <- mu + sigma*qnorm(alpha)
    
    # Periodic cleanup for very long runs
    if (t %% 100 == 0) {
      gc()
      try(detach("package:rugarch", unload = TRUE), silent = TRUE)
      suppressMessages(library(rugarch))
    }
    
    setTxtProgressBar(pb,t)
  }
  
  return(VaR)
}



roll_var_nts <- function(r, alpha = 0.01, window = 500) {
  n <- length(r)
  VaR <- rep(NA_real_, n)
  
  pb <- txtProgressBar(min = window + 1, max = n, style = 3)   # initialize progress bar
  
  for (t in (window + 1):n) {
    x <- r[(t - window):(t - 1)]
    fit <- fitnts(x)
    ntsparam <- as.numeric(fit[c("alpha", "theta", "beta", "gamma", "mu")])
    VaR[t] <- qnts(alpha, ntsparam)
    
    setTxtProgressBar(pb, t)   # update progress bar
  }
  
  close(pb)
  VaR
}

roll_CVaR_nts <- function(r,alpha=0.01,window=500){
  n<-length(r)
  CVaR<-rep(NA_real_,n)
  
  pb <- txtProgressBar(min=window+1,max=n,style=3) # initialize progress bar
  
  for(t in (window+1):n){
    x<-r[(t-window):(t-1)]
    fit<-fitnts(x)
    ntsparam <- as.numeric(fit[c("alpha", "theta", "beta", "gamma", "mu")])
    CVaR[t] <- cvarnts(alpha, ntsparam)
    
    setTxtProgressBar(pb,t)
  }
  
  close(pb)
  CVaR
}



kupiec_test<-function(hit,alpha){
  N<-length(hit)
  x<-sum(hit)
  pi_hat<-x/N
  
  # Calculating the likelihood ratio
  LR<- -2 * (log((1 - alpha)^(N - x) * alpha^x) - log((1 - pi_hat)^(N - x) * pi_hat^x))
  p_value <- 1- pchisq(LR,df=1)
  
  return(list(LR=LR, p_value = p_value, exceedances = x, ratio = pi_hat))
}

christoffersen_test <- function(hit){
  n00 <- sum(hit[-1] == 0 & hit[-length(hit)] == 0)  # 0 → 0 
  n01 <- sum(hit[-1] == 1 & hit[-length(hit)] == 0)  # 0 → 1
  n10 <- sum(hit[-1] == 0 & hit[-length(hit)] == 1)  # 1 → 0
  n11 <- sum(hit[-1] == 1 & hit[-length(hit)] == 1)  # 1 → 1
  
  pi0 <- n01/(n00+n01)
  pi1 <- n11/(n10+n11)
  pi <-  (n01 + n11) / (n00 + n01 + n10 + n11)
  
  L1 <- ((1-pi)^(n00+n10))*pi^(n01+n11)
  L2 <- ((1-pi0)^(n00))*(pi0^(n01))*((1-pi1)^(n10))*((pi1)^(n11))
  
  LR_ind <- -2 * log(L1/L2)
  p_val <- 1 - pchisq(LR_ind,df=1)
  return(list(LR_ind = LR_ind, p_value = p_val,counts = c(n00=n00,n01=n01,n10=n10,n11=n11),
              probs  = c(pi0=pi0, pi1=pi1, pi=pi)))
}

christoffersen_combined <- function(hit, alpha) {
  kup <- kupiec_test(hit, alpha)
  chr <- christoffersen_test(hit)
  LR_cc <- kup$LR + chr$LR_ind
  p_value <- 1 - pchisq(LR_cc, df = 2)
  list(LR_cc = LR_cc, p_value = p_value)
}

# Extending parametrization to ARMA-GARCH NTS to capture time varying mean and time varying volatility 

source("C:/Users/sahil/Downloads/MNTS_Optimization/MNTS_Optimization/readStockData.R")
source("C:/Users/sahil/Downloads/MNTS_Optimization/MNTS_Optimization/func_findOptport.R")
source("C:/Users/sahil/Downloads/MNTS_Optimization/MNTS_Optimization/miscellaneous_tools.R")
source("C:/Users/sahil/Downloads/MNTS_Optimization/MNTS_Optimization/fcts_mntsARMAGARCH_paramest.R")

# Modifying the function before using it

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
    library(spatstat)
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

# Creating the non parallel version for parameter estimation
paramest_ind_armagarch <- function(ret, numofelements, yTickers=NULL){
  muvec <- matrix(nrow = numofelements, ncol = 1)
  sigvec <- matrix(nrow = numofelements, ncol = 1)
  ntsparamMtx <- matrix(nrow = numofelements, ncol = 5)
  colnames(ntsparamMtx) <- c("mu", "sigma", "alpha", "theta", "beta")
  rownames(ntsparamMtx) <- yTickers
  
  for(n in 1:numofelements){
    tspec <- rugarch::ugarchspec(mean.model = list(armaOrder = c(1,1)),
                                 variance.model = list(garchOrder = c(1,1), model = "sGARCH"),
                                 distribution.model = "std")
    
    tf <- rugarch::ugarchfit(spec = tspec, data = ret[,n],solver="hybrid")
    
    v1 <- tf@fit$sigma
    resi <- tf@fit$residuals / v1
    tforc <- rugarch::ugarchforecast(tf, n.ahead=1)
    mu <- tforc@forecast$seriesFor[1]
    sig <- tforc@forecast$sigmaFor[1]
    
    
    stdntsparam <- temStaR::fitstdnts(resi)
    ntsparamMtx[n,] <- c(mu, sig, stdntsparam["alpha"], stdntsparam["theta"], stdntsparam["beta"])
  }
  return(list(ntsparam=ntsparamMtx))
}


roll_var_arma_garch_nts <- function(r, alpha = 0.01, window = 500,yTicker = "^SPY") {
  n <- length(r)
  VaR <- rep(NA_real_, n)
  numofelements <- 1
  
  pb <- txtProgressBar(min = window + 1, max = n, style = 3)   # initialize progress bar
  
  for (t in (window + 1):n) {
    x <- as.matrix(r[(t - window):(t - 1)])
    
    marketparam <- paramest_ind_armagarch(
      ret = x,
      numofelements = 1
    )
    mu    <- marketparam$ntsparam[1, "mu"]
    sigma <- marketparam$ntsparam[1, "sigma"]
    alpha_hat <- marketparam$ntsparam[1, "alpha"]
    theta_hat <- marketparam$ntsparam[1, "theta"]
    beta_hat  <- marketparam$ntsparam[1, "beta"]
    VaR[t] <-  mu + sigma * temStaR::qnts(alpha, c(alpha_hat, theta_hat, beta_hat))
    
    
    # Doing a memory cleanup every 50 iterations
    
    if (t %% 50 == 0) {
      message("🔄 Memory cleanup at t = ", t)
      gc()
      try(detach("package:rugarch", unload = TRUE), silent = TRUE)
      suppressMessages(library(rugarch))
    }
    
    setTxtProgressBar(pb, t)   # update progress bar
  }
  
  close(pb)
  VaR
}

roll_cvar_arma_garch_nts <- function(r, alpha = 0.01, window = 500,yTicker = "^SPY") {
  n <- length(r)
  CVaR <- rep(NA_real_, n)
  numofelements <- 1
  
  pb <- txtProgressBar(min = window + 1, max = n, style = 3)   # initialize progress bar
  
  for (t in (window + 1):n) {
    x <- as.matrix(r[(t - window):(t - 1)])
    
    marketparam <- paramest_ind_armagarch(
      ret = x,
      numofelements = 1
    )
    mu    <- marketparam$ntsparam[1, "mu"]
    sigma <- marketparam$ntsparam[1, "sigma"]
    alpha_hat <- marketparam$ntsparam[1, "alpha"]
    theta_hat <- marketparam$ntsparam[1, "theta"]
    beta_hat  <- marketparam$ntsparam[1, "beta"]
    
    CVaR[t] <- mu - sigma * cvarnts(alpha, c(alpha_hat, theta_hat, beta_hat)) # I made a change to this line of code 
    
    
    # Doing a memory cleanup every 50 iterations
    
    if (t %% 50 == 0) {
      message("🔄 Memory cleanup at t = ", t)
      gc()
      try(detach("package:rugarch", unload = TRUE), silent = TRUE)
      suppressMessages(library(rugarch))
    }
    
    setTxtProgressBar(pb, t)   # update progress bar
  }
  
  close(pb)
  CVaR
}



roll_cvar_norm <- function(r,alpha=0.01,window=500){
  n <- length(r)
  CVaR <- rep(NA_real_, n)
  numofelements <- 1
  z_alpha <-qnorm(alpha)
  phi_alpha <- dnorm(z_alpha)
  
  for (t in (window + 1):n){
    x<-r[(t-window):(t-1)]
    m<-mean(x)
    s<-sd(x)
    CVaR[t]<- m-s*(phi_alpha/alpha)
  }
 return(CVaR) 
}

VaR_arma_garch_nts_500 <-roll_var_arma_garch_nts(ret,window=500)
Var_nts_500 <- roll_var_nts(ret,window=500)
Var_norm_500 <- roll_var_norm(ret,window = 500)
VaR_t_500 <-roll_var_t(ret,window=500,nu = 5) # Our function uses 5 degrees of freedom
VaR_arma_garch_normal_500<-roll_var_arma_garch_norm(ret,window=500)



CVaR_nts_500 <- roll_CVaR_nts(ret,alpha=0.01,window=500)
CVaR_arma_garch_nts_500 <- roll_cvar_arma_garch_nts(ret,window=500)
CVaR_norm_500 <- roll_cvar_norm(ret,window = 500)



# Loading the data
VaR_arma_garch_nts_500<-readRDS("ARMA_GARCH_NTS_500_2000-2014.rds")
Var_nts_500<-readRDS("VAR_NTS_500_2000-2014.rds")

CVaR_nts_500<-readRDS("CVAR_NTS_500_2000-2014.rds")
CVaR_arma_garch_nts_500<-readRDS("CVAR_ARMA_GARCH_NTS_500_2000-2014.rds")



# Getting the realized return of the days past 500
window<-500
realized <- ret[(window + 1):length(ret)]


hit_norm <- ifelse(realized<Var_norm_500[(window + 1):length(Var_norm_500)],1,0)
hit_t <- ifelse(realized<VaR_t_500[(window + 1):length(VaR_t_500)],1,0)
hit_arma_garch_nts <- ifelse(realized<VaR_arma_garch_nts_500[(window + 1):length(VaR_arma_garch_nts_500)],1,0)
hit_nts <- ifelse(realized<Var_nts_500[(window + 1):length(Var_nts_500)],1,0)
hit_arma_garch_normal <- ifelse(realized<VaR_arma_garch_normal_500[(window + 1):length(VaR_arma_garch_normal_500)],1,0)

hit_nts


kupiec_test(hit_norm,0.01)
kupiec_test(hit_t,0.01)
kupiec_test(hit_nts,0.01)
kupiec_test(hit_arma_garch_nts,0.01)
kupiec_test(hit_arma_garch_normal,0.01)


christoffersen_test(hit_norm)
christoffersen_test(hit_t)
christoffersen_test(hit_arma_garch_nts)
christoffersen_test(hit_nts)
christoffersen_test(hit_arma_garch_normal)


christoffersen_combined(hit_norm,0.01)
christoffersen_combined(hit_t,0.01)
christoffersen_combined(hit_nts,0.01)
christoffersen_combined(hit_arma_garch_nts,0.01)
christoffersen_combined(hit_arma_garch_normal,0.01)


results <- data.frame(
  Model = c("Normal", "Student-t (5 df)", "NTS","ARMA GARCH NTS","ARMA GARCH Normal"),
  Kupiec_p = c(kupiec_test(hit_norm, 0.01)$p_value,
               kupiec_test(hit_t, 0.01)$p_value,
               kupiec_test(hit_nts, 0.01)$p_value,
               kupiec_test(hit_arma_garch_nts,0.01)$p_value,
               kupiec_test(hit_arma_garch_normal,0.01)$p_value),
  
  Kupiec_Ratio_result = c(kupiec_test(hit_norm, 0.01)$ratio,
                          kupiec_test(hit_t, 0.01)$ratio,
                          kupiec_test(hit_nts, 0.01)$ratio,
                          kupiec_test(hit_arma_garch_nts,0.01)$ratio,
                          kupiec_test(hit_arma_garch_normal,0.01)$ratio),
  Christoffersen_p = c(christoffersen_test(hit_norm)$p_value,
                       christoffersen_test(hit_t)$p_value,
                       christoffersen_test(hit_nts)$p_value,
                       christoffersen_test(hit_arma_garch_nts)$p_value,
                       christoffersen_test(hit_arma_garch_normal)$p_value)
)
print(results)

# Analysis of the result:
# We can see that Normal distribution and Student t distribution with 5 degrees of freedom have extremely low p values with the Kupiec Ratio results that are far from VaR Alpha %.
# NTS and ARMA GARCH NTS have higher P values with hit ratio of 1.45% and 1.6% respectively that's closer to alpha of 1% that was selected. The standard NTS fitted returns performed slighly better than ARMA GARCH NTS in the Kupiec test.
# For the Christoffersen test, normal and distribution show low p values with clustered volatility results, where as NTS and ARMA GARCH NTS show much better results. ARMA GARCH shows statistical significance in the christoffersen test, showing that it captures temporal dependency of in exceedance, capturing volatility clustering better. 

# Plotting the exceedance

realized <- ret[(window + 1):length(ret)]
VaR_norm_valid <- Var_norm_500[(window + 1):length(ret)]
VaR_t_valid    <- VaR_t_500[(window + 1):length(ret)]
VaR_nts_valid  <- Var_nts_500[(window + 1):length(ret)]
VaR_arma_garch_nts_valid  <- VaR_arma_garch_nts_500[(window + 1):length(ret)]
VaR_arma_garch_norm_valid <- VaR_arma_garch_normal_500[(window + 1):length(ret)]

CVaR_ARMA_GARCH_NTS_valid <- CVaR_arma_garch_nts_500[(window + 1):length(ret)]
CVaR_NTS_valid <- CVaR_nts_500[(window + 1):length(ret)]
CVaR_norm_500_valid <- CVaR_norm_500[(window+1):length(ret)]



exceed_norm <- realized < VaR_norm_valid
exceed_t    <- realized < VaR_t_valid
exceed_nts  <- realized < VaR_nts_valid
exceed_arma_garch_nts <- realized < VaR_arma_garch_nts_valid
exceed_arma_garch_norm <- realized < VaR_arma_garch_norm_valid

Date <- seq(from = begindate, to = enddate, length.out = length(ret))

plot_df <- data.frame(
  Date = Date[(window+1):length(ret)],
  Return = realized,
  VaR_Normal = VaR_norm_valid,
  VaR_t = VaR_t_valid,
  VaR_NTS = VaR_nts_valid,
  VaR_ARMA_GARCH_NTS = VaR_arma_garch_nts_valid,
  VaR_ARMA_GARCH_Normal = VaR_arma_garch_norm_valid,
  Exceed_Normal = exceed_norm,
  Exceed_t = exceed_t,
  Exceed_NTS = exceed_nts,
  Exceed_ARMA_GARCH_NTS = exceed_arma_garch_nts,
  Exceed_ARMA_GARCH_Norm = exceed_arma_garch_norm,
  CVaR_ARMA_GARCH_NTS = CVaR_ARMA_GARCH_NTS_valid,
  CVaR_NTS = CVaR_NTS_valid
)
plot_df



ggplot(plot_df, aes(x = Date)) +
  geom_line(aes(y = Return, color = "Return"), size = 0.4) +
  geom_line(aes(y = VaR_Normal, color = "VaR_Normal"), size = 0.7, alpha = 0.7) +
  geom_line(aes(y = VaR_t, color = "VaR_t"), size = 0.7, alpha = 0.7) +
  geom_line(aes(y = VaR_NTS, color = "VaR_NTS"), size = 0.7, alpha = 0.7) +
  geom_line(aes(y = VaR_ARMA_GARCH_NTS, color = "VaR_ARMA_GARCH_NTS"), size = 0.7, alpha = 0.7) +
  # geom_line(aes(y = VaR_ARMA_GARCH_Normal, color = "VaR_ARMA_GARCH_Normal"), size = 0.7, alpha = 0.7) +
  
  # labs(title = "Rolling 1% VaR Backtest (Normal vs t vs NTS vs ARMA GARCH NTS vs ARMA GARCH)",
  #      y = "Log Return", x = "Time") +
  scale_color_manual(values = c("Return" = "black", 
                                "VaR_Normal" = "red", 
                                "VaR_t" = "blue", 
                                "VaR_NTS" = "darkgreen", 
                                "VaR_ARMA_GARCH_NTS" = "orange",
                                "VaR_ARMA_GARCH_Normal"="purple")) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  theme_minimal() +
  theme(legend.title = element_blank())  

plot_df <- data.frame(
  Date = Date[(window+1):length(ret)],
  Return = realized,
  VaR_Normal = VaR_norm_valid,
  VaR_t = VaR_t_valid,
  VaR_NTS = VaR_nts_valid,
  VaR_ARMA_GARCH_NTS = VaR_arma_garch_nts_valid,
  CVaR_norm = CVaR_norm_500_valid,
  CVaR_ARMA_GARCH_NTS = CVaR_ARMA_GARCH_NTS_valid,
  CVaR_NTS = -CVaR_NTS_valid
)
plot_df

ggplot(plot_df, aes(x = Date)) +
  geom_line(aes(y = Return, color = "Return"), size = 0.4) +
  geom_line(aes(y = CVaR_ARMA_GARCH_NTS, color = "CVaR_ARMA_GARCH_NTS"), size = 0.7, alpha = 0.7) +
  # geom_line(aes(y = CVaR_NTS, color = "CVaR_NTS"), size = 0.7, alpha = 0.7) +
  # geom_line(aes(y= VaR_NTS,color = "VaR_NTS"), size = 0.7,alpha = 0.7)+
  # geom_line(aes(y= VaR_Normal,color = "VaR_Normal"), size = 0.7,alpha = 0.7)+
  # geom_line(aes(y= CVaR_norm,color = "CVaR_norm"), size = 0.7,alpha = 0.7)+
  
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  
  scale_color_manual(values = c("Return" = "black", 
                                "VaR_NTS" = "red", 
                                "CVaR_NTS" = "blue",
                                "CVaR_ARMA_GARCH_NTS" = "purple",
                                "VaR_Normal"="yellow",
                                "CVaR_norm"="green")) +
  theme_minimal() +
  theme(legend.title = element_blank())  




# Backtesting Expected Shortfall function

es_test1 <- function(ret,VaR,CVar){
  
  ret <- ret[501:length(ret)]
  VaR <- VaR[501:length(VaR)]
  CVar <- CVar[501:length(CVar)]
  
  # Defining the indicator function
  exceedances = as.numeric(ret< VaR)
  sum_exceedance = sum(exceedances)
  
  # Computing the Z1 statistic
  
  
  Z1 = (sum((ret/CVar)*exceedances)/sum_exceedance) + 1
  return(Z1)
}

es_test2 <- function(ret,VaR,CVar,alpha=0.01){
  ret <- ret[501:length(ret)]
  VaR <- VaR[501:length(VaR)]
  CVar <- CVar[501:length(CVar)]
  T_ <- length(ret)
  
  exceedances = as.numeric(ret<VaR)
  sum_exceedance = sum(exceedances)
  
  Z2 = sum((ret*exceedances)/(T_*alpha*CVar))+1
  
  return(Z2)
}


es_test1(ret,Var_nts_500,-CVaR_nts_500)
es_test1(ret,Var_norm_500,CVaR_norm_500)

es_test2(ret,Var_nts_500,-CVaR_nts_500)
es_test2(ret,Var_norm_500,CVaR_norm_500)







