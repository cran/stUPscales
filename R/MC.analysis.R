# Monte-Carlo simulations analysis
# author: J.A. Torres-Matallana
# organization: Luxembourg Institute of Science and Technology (LIST), Luxembourg
#               Wagenigen University and Research Centre (WUR), Wageningen, The Netherlands
# date: 30.04.2016 - 22.10.2016


# setGeneric("MC.analysis", function(x, delta, qUpper, p1.det, sim.det, event.ini, event.end, ntick, summ) standardGeneric("MC.analysis"))
# setMethod("MC.analysis", signature = c("list", "numeric", "character", "data.frame", "list", "POSIXct", "POSIXct", "numeric", "list"), 

MC.analysis <- function(x, delta, qUpper, p1.det, sim.det, 
                        event.ini, event.end, ntick, summ.data = NULL){
  
  sim1 <- x[["sim1"]]
  
  P1 <- p1.det
  sum(P1[,2])
  
  sim <- sim.det
  
  #-----------------------------------------------------------------------------------------------------------
  #  define precipitation input data
  #-----------------------------------------------------------------------------------------------------------
  
  a <- IsReg.ts(P1, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  a[[1]]
  ts.P1 <- a[[2]]
  colnames(ts.P1) <- c("p")
  
  p1 <- P1
  colnames(p1) <- c("time", "p")
  
  #-----------------------------------------------------------------------------------------------------------
  #  deterministic sim
  #-----------------------------------------------------------------------------------------------------------
  det  <- data.frame(as.POSIXct(sim[[1]][[1]][[2]]))
  det  <- cbind(det, sim[[1]][[1]][[9]])
  det  <- cbind(det, sim[[1]][[1]][[10]])
  # det  <- cbind(det, sim[[1]][[1]][[22]])
  det  <- cbind(det, sim[[1]][[1]][[11]])
  det  <- cbind(det, sim[[1]][[1]][[12]])
  det  <- cbind(det, sim[[1]][[1]][[13]])
  det  <- cbind(det, sim[[1]][[1]][[14]])
  colnames(det) <- c("time", "VChamber", "Vsv", "BCODsv", "BNH4sv", "CCODsv", "CNH4sv")
  det.var <- c("VChamber", "Vsv", "", "BCODsv", "BNH4sv", "CCODsv", "CNH4sv")
  head(det)
  
  det.event <- det[((det$time >= event.ini) & (det$time <= event.end)),]
  
  ##============================================================================================================================
  ##  uncertainty confidence intervals
  ##============================================================================================================================
  
  ##-------------------------------------------------------------------------------------------
  ## VChamber, Vsv and Qsv
  ##-------------------------------------------------------------------------------------------
  ## VChamber
  ifelse(length(summ.data) == 0,{
    data.VChamber <- t(as.data.frame(sim1[,"mcVChamber"]))
    summ.VChamber <- MC.summary(p1, data.VChamber)
    summ <- list(summ.VChamber = summ.VChamber)
    # save(summ, file="summ.RData") # uncomment to save file
    summ.VChamber
  },{
    summ.VChamber <- summ[["summ.VChamber"]]
    summ.Vol <- summ[["summ.Vol"]]
    # summ.Q <- summ[["summ.Q"]]
    summ.BCOD <- summ[["summ.BCOD"]]
    summ.CCOD <- summ[["summ.CCOD"]]
    summ.BNH4 <- summ[["summ.BNH4"]]
    summ.CNH4 <- summ[["summ.CNH4"]]
  })
  
  summ.VChamber.agg <- MC.summary.agg(summ.VChamber, det, delta, mean, sum)
  
  namePlot.VChamber <- "VChamber_MC-uncertainty"; ylab="V chamber [m3]"
  PlotMC.season(summ.VChamber.agg, paste(namePlot.VChamber, "_", round(delta, 0), "min", sep=""), ylab, qUpper)
  
  VChamber.event <- summ.VChamber[((as.POSIXct(summ.VChamber$time) >= event.ini) & (as.POSIXct(summ.VChamber$time) <= event.end)),]
  VChamber.event.agg <- MC.summary.agg(VChamber.event, det.event, delta, mean, sum)
  
  ## Vsv
  if(length(summ.data) == 0){
    data.Vol <- t(as.data.frame(sim1[,"mcVsv"]))
    summ.Vol <- MC.summary(p1, data.Vol)
    summ <- append(summ, list(summ.Vol = summ.Vol))
    # save(summ, file="summ.RData") # uncomment to save file
    summ.Vol
  }
  
  summ.Vol.agg <- MC.summary.agg(summ.Vol, det, delta, mean, sum)
  
  namePlot.Vol <- paste("Vsv_MC-uncertainty", sep=""); ylab="V overflow [m3]"
  PlotMC.season(summ.Vol.agg, paste(namePlot.Vol, "_", round(delta, 0), "min", sep=""), ylab, qUpper)
  
  V.event <- summ.Vol[((as.POSIXct(summ.Vol$time) >= event.ini) & (as.POSIXct(summ.Vol$time) <= event.end)),]
  V.event.agg <- MC.summary.agg(V.event, det.event, delta, mean, sum)
  
  namePlot.VChamberV <- "VChamber_Vsv_MC-uncertainty"; ylab="V chamber [m3]"; ylab1="V overflow [m3]"
  PlotMC.event(summ = VChamber.event.agg, summ1 = V.event.agg, obs = 0, det.var = det.var[1], det.var1 = det.var[2], 
               namePlot = paste(namePlot.VChamberV, "_", round(delta, 0), "min", sep=""), 
               ylab = ylab, ylab1 = ylab1, ntick = ntick, qUpper = qUpper)
  
  rm(data.VChamber, data.Vol)
  
  ##-------------------------------------------------------------------------------------------
  ## CODsv
  ##-------------------------------------------------------------------------------------------
  ## BCODsv
  if(length(summ.data) == 0){
    data.BCOD <- t(as.data.frame(sim1[,"mcBCODsv"]))
    summ.BCOD <- MC.summary(p1, data.BCOD)
    summ <- append(summ, list(summ.BCOD=summ.BCOD))
    # save(summ, file="summ.RData") # uncomment to save file
    summ.BCOD
  }
  
  summ.BCOD.agg <- MC.summary.agg(summ.BCOD, det, delta, mean, sum)
  
  namePlot.BCOD <- "BCODsv_MC-uncertainty"; ylab="BCOD overflow [kg]"
  PlotMC.season(summ.BCOD.agg, paste(namePlot.BCOD, "_", round(delta, 0), "min", sep=""), ylab, qUpper)
  
  BCOD.event <- summ.BCOD[((as.POSIXct(summ.BCOD$time) >= event.ini) & (as.POSIXct(summ.BCOD$time) <= event.end)),]
  BCOD.event.agg <- MC.summary.agg(BCOD.event, det.event, delta, mean, sum)
  
  ## CCODsv
  if(length(summ.data) == 0){
    data.CCOD <- t(as.data.frame(sim1[,"mcCCODsv"]))
    summ.CCOD <- MC.summary(p1, data.CCOD)
    summ <- append(summ, list(summ.CCOD=summ.CCOD))
    # save(summ, file="summ.RData") # uncomment to save file
    summ.CCOD
  }
  
  summ.CCOD.agg <- MC.summary.agg(summ.CCOD, det, delta, mean, sum)
  
  namePlot.CCOD <- "CCODsv_MC-uncertainty"; ylab1="CCOD overflow [mg/l]"
  PlotMC.season(summ.CCOD.agg, paste(namePlot.CCOD, "_", round(delta, 0), "min", sep=""), ylab1, qUpper)
  
  namePlot.BCCOD <- "BCCODsv_MC-uncertainty";
  CCOD.event <- summ.CCOD[((as.POSIXct(summ.CCOD$time) >= event.ini) & (as.POSIXct(summ.CCOD$time) <= event.end)),]
  CCOD.event.agg <- MC.summary.agg(CCOD.event, det.event, delta, mean, sum)
  PlotMC.event(BCOD.event.agg, CCOD.event.agg, 0, det.var[4], det.var[6], paste(namePlot.BCCOD, "_", round(delta, 0), "min", sep=""), 
               ylab, ylab1, ntick, qUpper)
  
  ## removing no needed objects
  rm(data.BCOD, data.CCOD)
  
  ##-------------------------------------------------------------------------------------------
  ## NH4sv
  ##-------------------------------------------------------------------------------------------
  ## BNH4sv
  if(length(summ.data) == 0){
    data.BNH4 <- t(as.data.frame(sim1[,"mcBNH4sv"]))
    summ.BNH4 <- MC.summary(p1, data.BNH4)
    summ <- append(summ, list(summ.BNH4=summ.BNH4))
    # save(summ, file="summ.RData") # uncomment to save file
    summ.BNH4
  }
  
  summ.BNH4.agg <- MC.summary.agg(summ.BNH4, det, delta, mean, sum)
  
  namePlot.BNH4 <- "BNH4sv_MC-uncertainty"; ylab="BNH4 overflow [kg]"
  PlotMC.season(summ.BNH4.agg, paste(namePlot.BNH4, "_", round(delta, 0), "min", sep=""),ylab, qUpper)
  
  BNH4.event <- summ.BNH4[((as.POSIXct(summ.BNH4$time) >= event.ini) & (as.POSIXct(summ.BNH4$time) <= event.end)),]
  BNH4.event.agg <- MC.summary.agg(BNH4.event, det.event, delta, mean, sum)
  
  ## CNH4sv
  if(length(summ.data) == 0){
    data.CNH4 <- t(as.data.frame(sim1[,"mcCNH4sv"]))
    summ.CNH4 <- MC.summary(p1, data.CNH4)
    summ <- append(summ, list(summ.CNH4=summ.CNH4))
    # save(summ, file="summ.RData") # uncomment to save file
    summ.CNH4
  }
  
  summ.CNH4.agg <- MC.summary.agg(summ.CNH4, det, delta, mean, sum)
  
  namePlot.CNH4 <- "CNH4sv_MC-uncertainty"; ylab1="CNH4 overflow [mg/l]"
  PlotMC.season(summ.CNH4.agg, paste(namePlot.CNH4, "_", round(delta, 0), "min", sep=""),  ylab1, qUpper)
  
  namePlot.BCNH4 <- "BCNH4sv_MC-uncertainty";
  CNH4.event <- summ.CNH4[((as.POSIXct(summ.CNH4$time) >= event.ini) & (as.POSIXct(summ.CNH4$time) <= event.end)),]
  CNH4.event.agg <- MC.summary.agg(CNH4.event, det.event, delta, mean, sum)
  PlotMC.event(BNH4.event.agg, CNH4.event.agg, 0, det.var[5], det.var[7], paste(namePlot.BCNH4, "_", round(delta, 0), "min", sep=""), 
               ylab, ylab1, ntick, qUpper)
  
  namePlot.CCODCNH4 <- "CCODsv-CNH4sv_MC-uncertainty"; ylab="CCOD overflow [mg/l]"; ylab1="CNH4 overflow [mg/l]"
  PlotMC.event(CCOD.event.agg, CNH4.event.agg, 0, det.var[6], det.var[7], paste(namePlot.CCODCNH4, "_", round(delta, 0), "min", sep=""), 
               ylab, ylab1, ntick, qUpper)
  
  ## removing no needed objects
  rm(data.BNH4, data.CNH4)
  
  ##============================================================================================================================
  ##  variances
  ##============================================================================================================================
  
  ##-------------------------------------------------------------------------------------------
  ## computing variances
  ##-------------------------------------------------------------------------------------------
  
  ## computing variance matrix
  a <- as.data.frame(rapply(summ, f = mean, how="unlist", classes="numeric"))
  a.names <- rownames(a)
  
  b <- as.data.frame(rapply(summ,f = function(x)quantile(x, probs = 0.05), how="unlist", classes="numeric"))
  b.names <- rownames(b)
  
  c <- as.data.frame(rapply(summ,f = function(x)quantile(x, probs = 0.25), how="unlist", classes="numeric"))
  c.names <- rownames(c)
  
  d <- as.data.frame(rapply(summ,f = function(x)quantile(x, probs = 0.50), how="unlist", classes="numeric"))
  d.names <- rownames(d)
  
  e <- as.data.frame(rapply(summ,f = function(x)quantile(x, probs = 0.75), how="unlist", classes="numeric"))
  e.names <- rownames(e)
  
  f <- as.data.frame(rapply(summ,f = function(x)quantile(x, probs = 0.95), how="unlist", classes="numeric"))
  f.names <- rownames(f)
  
  g <- as.data.frame(rapply(summ,f = function(x)quantile(x, probs = 0.995), how="unlist", classes="numeric"))
  g.names <- rownames(g)
  
  h <- as.data.frame(rapply(summ,f = function(x)quantile(x, probs = 0.999), how="unlist", classes="numeric"))
  h.names <- rownames(h)
  
  ii<- as.data.frame(rapply(summ,f = max, how="unlist", classes="numeric"))
  i.names <- rownames(ii)
  
  jj <- as.data.frame(rapply(summ,f = sum, how="unlist", classes="numeric"))
  j.names <- rownames(jj)
  
  length(a.names)
  names <- cbind(a.names, b.names, c.names, d.names, e.names, f.names, g.names, h.names, i.names, j.names)
  
  getwd()
  
  all.equal(a.names, b.names, c.names, d.names, e.names, f.names, g.names, h.names, i.names, j.names)
  
  variable <- rownames(a)
  ab   <- cbind(variable, a, b, c, d, e, f, g, h, ii, jj)
  dim(ab)
  id1  <- grep("Variance", variable)
  id2  <- grep("Mean", variable)
  id3  <- grep("q05", variable)
  id4  <- grep("q25", variable)
  id5  <- grep("q50",  variable)
  id6  <- grep("q75",  variable)
  id7  <- grep("q95",  variable)
  id8  <- grep("q995", variable)
  id9  <- grep("q999", variable)
  id10 <- grep("Max",  variable)
  id11 <- grep("Sum",  variable)
  
  id  <- c(id1, id2, id3, id4, id5, id6, id7, id8, id9, id10, id11)
  variance <- ab[id,]
  dim(variance)
  colnames(variance)[2:11] <- c("Mean", "q05", "q25","q50", "q75", "q95", "q995", "q999", "Max", "Sum")
  rownames(variance) <- NULL
  
  ## saving variance matrix
  # save(variance, file="variance.RData") # uncomment to save file
  
  ## ending
  print("End MC analysis. Please check your output folder.")
  
  return(list(summ=summ, variance=variance))
}
#)