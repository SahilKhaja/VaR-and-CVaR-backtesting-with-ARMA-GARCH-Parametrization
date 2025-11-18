source("readHistoryStockDataFromYahoo.R")

readLogReturn_interval <- function(yTickers, begindate, enddate, intervalDateNumber) {
  p = readHistoryStockDataFromYahoo(yTickers, begindate, enddate)
  sub_price = p$Adj.Close[seq(1,length(p$Adj.Close), by = intervalDateNumber)]
}

readdata <- function(yTickers, numofelements, begindate, enddate, intervalDateNumber ){
  p <- readLogReturn_interval(yTickers[1], begindate, enddate, intervalDateNumber)
  screensize <- length(p)-1
  ret <- matrix(0, screensize, numofelements)
  for( i in 1:numofelements)
  {
    p <- readLogReturn_interval(yTickers[i], begindate, enddate,intervalDateNumber)
    tempr <- matrix(diff(log(p)), 1, length(p)-1)
    ret[,i] <- tempr
  }
  return(ret)
}