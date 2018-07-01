# aggregation of data.frame function
# author: J.A. Torres-Matallana
# organization: Luxembourg Institute of Science and Technology (LIST), Luxembourg
#               Wagenigen University and Research Centre (WUR), Wageningen, The Netherlands   
# date: 15.07.2015 - 15.07.2015

# data     <- P1
# nameData <- deparse(substitute(P1))
# var      <- 1
# delta    <- 1
# func     <- "sum"

# data     <- wlt_obs
# nameData <- "wlt_obs"
# var      <- 1
# delta    <- 10
# func     <- "mean"

Agg.t <- function(data, nameData, delta, func, namePlot){
  # data <- var; nameData <- var.name; delta <- 60; func <- "mean"; namePlot <- "hourly"
  # data <- wlt_obs; nameData <- "wlt_obs"; delta <- 1; func <- "mean"; namePlot <- "test"
  #---------------------------------------------------------------------------------------------------------
  # aggregating to 10, 30, 60, etc. min resolution
  #---------------------------------------------------------------------------------------------------------
  tt <- as.POSIXct(data[,1], tz="UTC")
  
  # delta min
  dt <- 60/1*delta # 60_s/1_min * delta_min = dt_s
  bucket = (tt) - as.numeric(tt) %% dt
  namePlot <- paste(namePlot, "(res =", delta, "min)", sep=" ")
  
  if(nameData == "P1"){
    P1 <- data
    ts <- aggregate(P1[,2], list(bucket), func)
    # ts[,3] <- NA
    colnames(ts) <- c("time", "Rainfall")
    # length(P1$Rainfall)
    # length(ts$Rainfall)
    # head(ts)
    
    par(mfrow = c(2, 1))
    par(mar = rep(2, 4)) #------------------------------------------ added after MC set-up
    plot(P1[,1],P1[,2], type="l", main=namePlot) #------------------------------------------ commented after MC set-up
    plot(ts$time,ts$Rainfall, type="l")
      
    P1 <- ts
    # head(P1)
    # save(P1, file="P1.RData")
    
    return(P1)
  }else{
    obs <- data
    ts <- aggregate(data[-1], list(bucket), func)
    
    pdf(paste(namePlot, ".pdf", sep=""), pointsize=10)
    par(mfrow = c(2,1))
    par(cex.lab=1, cex.axis=1., cex.main = 1.5)
    # length(obs$time); length(obs$value)
    # plot(obs$tt,obs$value, type="l", main="Original time series", xlab = "Time", ylab = nameData)#------ commented after MC set-up
    plot(obs[,1],obs[,2], type="l", main="Original time series", xlab = "Time", ylab = nameData)#------ commented after MC set-up
    plot(ts[,1],ts[,2], type="l", main=namePlot, xlab = "Time", ylab = nameData)
    dev.off()
    
    colnames(ts) <- c("time", colnames(obs)[-1])
    aggregated <- ts
    #save(obs, file="obs.RData")
    
    return(aggregated)
  }
}