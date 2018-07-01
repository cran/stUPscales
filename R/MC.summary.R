# Monte-Carlo summary
# author: J.A. Torres-Matallana
# organization: Luxembourg Institute of Science and Technology (LIST), Luxembourg
#               Wagenigen University and Research Centre (WUR), Wageningen, The Netherlands   
# date: 17.10.2015 - 10.10.2016

MC.summary <- function(p1, data){
  # library(parallel)
  # library(foreach)
  # library(doParallel)
  # 
  # data <- data.VChamber
  # cores <- 3
  # 
  # system.time({
  #   cl <- makeCluster(cores, outfile="")
  #   registerDoParallel(cl, cores=cores)
  #   numCores <- detectCores()
  #   # equivalent to lapply(1:3, sqrt)
  #   summary <- foreach(i=1:ncol(data), .verbose=FALSE, .combine = "rbind") %dopar% {
  #     data.frame(Mean = mean(data[,i]),
  #                Sd = sd(data[,i]),
  #                Variance = var(data[,i]),
  #                q05 = quantile(data[,i], probs = 0.05),
  #                q50 = quantile(data[,i], probs = 0.50),
  #                q95 = quantile(data[,i], probs = 0.95),
  #                q995 = quantile(data[,i], probs = 0.995),
  #                q999 = quantile(data[,i], probs = 0.999),
  #                Max = max(data[,i]),
  #                Sum = sum(data[,i])
  #     )
  #   }
  #   
  #   stopCluster(cl)
  #   closeAllConnections()
  # }
  # )

  # system.time({
  #   summ <- data.frame(
  #   Mean  = apply(data, 2, FUN = mean),
  #   Sd    = apply(data, 2, FUN = sd),
  #   Variance = apply(data, 2, FUN = var),
  #   #Min  = apply(data, 2, FUN = min),
  #   Max   = apply(data, 2, FUN = max),
  #   q05   = apply(data, 2, FUN=function(x)quantile(x, probs = 0.05)),
  #   #q25  = apply(data, 2, FUN=function(x)quantile(x, probs = 0.25)),
  #   q50   = apply(data, 2, FUN=function(x)quantile(x, probs = 0.5)),
  #   q95   = apply(data, 2, FUN=function(x)quantile(x, probs = 0.95)),
  #   q995  = apply(data, 2, FUN=function(x)quantile(x, probs = 0.995)),
  #   q999  = apply(data, 2, FUN=function(x)quantile(x, probs = 0.999)),
  #   Sum   = apply(data, 2, FUN = sum)
  #   )
  # })
  
  # head(summ); tail(summ)

## data.table implementation  
# system.time({
#   library(data.table)
#     data1 <- data.table(cbind(id=c(1:nrow(data)), data))
#     setkey(data1, id)
# 
#     Mean <- as.vector(data1[, lapply(.SD, mean), by = .I])
#     Sd   <- as.vector(data1[, lapply(.SD, sd), by = .I])
#     Variance <- as.vector(data1[, lapply(.SD, var), by = .I])
#     Max  <- as.vector(data1[, lapply(.SD, max), by = .I])
#     q05  <- as.vector(data1[, lapply(.SD, function(x)quantile(x, probs = 0.05)), by = .I])
#     q50  <- as.vector(data1[, lapply(.SD, function(x)quantile(x, probs = 0.50)), by = .I])
#     q95  <- as.vector(data1[, lapply(.SD, function(x)quantile(x, probs = 0.95)), by = .I])
#     q995  <- as.vector(data1[, lapply(.SD, function(x)quantile(x, probs = 0.995)), by = .I])
#     q999  <- as.vector(data1[, lapply(.SD, function(x)quantile(x, probs = 0.999)), by = .I])
#     Sum  <- as.vector(data1[, lapply(.SD, sum), by = .I])
# 
#     summ <- cbind(t(Mean), t(Sd), t(Variance), t(Max), t(q05), t(q50), t(q95), t(q995), t(q999), t(Sum))
#     summ <- summ[2:nrow(summ),]
#     colnames(summ) <- c("Mean", "Sd", "Variance", "Max", "q05", "q50", "q95", "q995", "q999","Sum")
# })  

## data.table implementation: way 2 
# system.time({
  # data <- data.VChamber
  requireNamespace("data.table")
  idx   <- rep(1,nrow(data))
  data1 <- data.table(cbind(idx, (data)))
  setkey(data1, idx)
  
  # summ1 <- data.table(Mean = numeric(0), Sd = numeric(0), Variance = numeric(0),
  #                     q05 = numeric(0), q50 = numeric(0), q95 = numeric(0),
  #                     q995 = numeric(0), q999 = numeric(0), Max = numeric(0),
  #                     Sum = numeric(0))
  
  data2 <- (t(data1))
  dim(data2)
  class(data2)
  x <- data2[3,]
  length(x)
  class(x)
  min(x); max(x)
  mean(x)
  
  Summ.mean <- function(x){
    mean(x[x > 1e-5]) 
  }
  
  summ2     <- data1[, list(
    Mean     = lapply(.SD, Summ.mean),
    Sum      = lapply(.SD, sum)),
    by       = idx]
  
  summ1     <- data1[, list(
    Mean     = lapply(.SD, mean),
    Sd       = lapply(.SD, sd),
    Variance = lapply(.SD, var),
    q05      = lapply(.SD, function(x)quantile(x, probs = 0.05)),
    q25      = lapply(.SD, function(x)quantile(x, probs = 0.25)),
    q50      = lapply(.SD, function(x)quantile(x, probs = 0.50)),
    q75      = lapply(.SD, function(x)quantile(x, probs = 0.75)),
    q95      = lapply(.SD, function(x)quantile(x, probs = 0.95)),
    q995     = lapply(.SD, function(x)quantile(x, probs = 0.995)),
    q999     = lapply(.SD, function(x)quantile(x, probs = 0.999)),
    Max      = lapply(.SD, max),
    Sum      = lapply(.SD, sum)),
    by       = idx]
# })
  col.names <- names(summ1)
  summ1 <- data.frame(matrix(unlist(summ1), nrow=nrow(summ1), byrow=F))
  colnames(summ1) <- col.names
# data1 <- data.table(cbind(id=c(1:ncol(data)), t(data)))
# setkey(data1, id)
# 
# Q95 <- function(x) quantile(x, probs = 0.95)
# 
# Mean <- data1[, list(Mean=sum(.SD)/(ncol(data1)-1), Sd=sd(.SD), Mean1=mean(.I)), by = "id"]

# Variance <- summ[,"Sd"]^2
# summ <- cbind(summ, Variance)
# max(summ1[,"Variance"] - Variance)

time <- p1$time
summ1 <- cbind(summ1, time)

p1 <- p1[,2]
# all.equal(mc$time, p1$time)
summ1 <- cbind(summ1, p1)

#   sum(summ[,"Sd"])
#   mean(summ[,"Variance"])
#   write.csv(summ, file="summ.csv")


#   head(summ)
#   Mean <- mean(summ[,"Mean"])
#   Sd <- mean(summ[,"Sd"])
#   q05 <- mean(summ[,"q05"])
#   q95 <- mean(summ[,"q95"])
#   q50 <- mean(summ[,"q50"])

return(summ1)
}

MC.summary.agg <- function(summ, det, delta, func.agg, func.agg.p){
  # summ <- summ.BNH4
  # summ <- summ.VChamber
  
  if(length(det) > 1) summ <- cbind(summ, det)
  
  tt <- as.POSIXct(summ$time, tz="UTC")
  
  # delta min
  dt <- 60/1*delta # 60_s/1_min * delta_min = dt_s
  tt.agg = (tt) - as.numeric(tt) %% dt
  
  summ1 <- as.data.frame(summ[, which(names(summ) != "time")])
  summ.agg <- aggregate(x=summ1[, names(summ1) != "p1"], by=list(tt.agg), FUN = func.agg)
  summP.agg <- aggregate(x=summ1$p1, by=list(tt.agg), FUN = func.agg.p)
  colnames(summ.agg)[1] <- "time"
  colnames(summP.agg)[2] <- "p1"
  summ.agg <- cbind(summ.agg, p1=summP.agg[,"p1"])
  summ.agg$month <- as.Date(cut(summ.agg$time, breaks = "month"))
  
  return(summ.agg)
}