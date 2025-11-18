
portRisk <- function(w, covMtx){
  w <- matrix(w, nrow = 1, ncol = length(w))
  return(w%*%covMtx%*%t(w))
}

findOptimalPortfolio <- function(meanvector, covMtx, benchmarkReturn){
  numofelements <- length(meanvector)
  w0 <- rep(1/numofelements, numofelements)
  optport <- fmincon(x0 = w0, 
                     fn = functional::Curry(portRisk, 
                                            covMtx = covMtx),
                     A = -matrix(meanvector, nrow = 1, ncol = numofelements), 
                     b = -benchmarkReturn,
                     Aeq = matrix(data = 1, nrow = 1, ncol = numofelements), 
                     beq = 1,
                     lb = rep(0, numofelements), 
                     ub = rep(1, numofelements))
  
  ExpectedReturn <- sum(optport$par*meanvector)
  returnvalue <- list(weight=optport$par, 
                      reward = ExpectedReturn, 
                      risk = optport$value)
  return( returnvalue )
}