# Monte-Carlo simulations analysis
# author: J.A. Torres-Matallana
# organization: Luxembourg Institute of Science and Technology (LIST), Luxembourg
#               Wagenigen University and Research Centre (WUR), Wageningen, The Netherlands
# date: 30.04.2016 - 22.10.2016


# setGeneric("MC.analysis", function(x, delta, qUpper, p1.det, sim.det, event.ini, event.end, ntick, summ) standardGeneric("MC.analysis"))
# setMethod("MC.analysis", signature = c("list", "numeric", "character", "data.frame", "list", "POSIXct", "POSIXct", "numeric", "list"), 

MC.analysis <- function(x, delta, qUpper, p1.det, sim.det, 
                        event.ini, event.end, ntick, summ.data = NULL){
  # x <- sims

  # delta <- 2*60    # 2 hours 
  # delta <- 1*60/6 # 10 minutes 
  # delta <- 1      # 1 minutes 
  
  # qUpper <- "q995" #  "q999", q995", "q99", "q95" "q50", "Mean", "q05"
  
  # redo = "FALSE"  # to re-draw plots
  # redo = "TRUE"  # to re-draw plots
  
  #-----------------------------------------------------------------------------------------------------------
  #  load variables in R 
  #-----------------------------------------------------------------------------------------------------------
  ## MC matrices and data
  # load(paste(folder, "sim1.RData", sep=""))
  # load(paste(folder, "mc.RData", sep=""))
  
  sim1 <- x[["sim1"]]
  
  # if(redo == "TRUE") load(paste(folder, "summ.RData", sep=""))
  
  ## precipitation data
  # load("ts.winter.RData")
  
  # data("Esch_Sure2010")
  # P1 <- Esch_Sure2010
  
  # ts-s-EschSure2010 time series number 526
  # load(paste(folder, "MC_ts_s_EschSure2010.RData", sep=""))
  # P1 <- ts.s.EschSure2010[, c(1,527)]
  P1 <- p1.det
  sum(P1[,2])
  # plot(P1[,1], P1[,2], typ="l")
  
  
  # Event7-8
  # load(paste(folder, "P1_Event7-8.RData", sep=""))
  # load(paste(folder, "ts.s.Event7-8.1min.RData", sep=""))
  # time <- seq(as.POSIXct(P1[1,1]), by = "1 min", length.out = nrow(ts.s.Event7_8))
  # p1 <- time
  # p <- rep(x = 0, length(p1))
  # p1 <- cbind.data.frame(p1, p)
  # p1[1:nrow(P1),2] <- P1[,2]
  # P1 <- p1
  # save(P1, file="P1_Event7-8.RData")
  # plot.ts(P1[,2])
  
  ## deterministic simulation
  # load("sim_Goesdorf_winter2009-2010.RData")
  # load(paste(folder, "sim_Goesdorf_EschSure2010_(peTS_qsTS).RData", sep=""))
  # load(paste(folder, "sim_Goesdorf_ts-s-EschSure2010-526_(peTS_qsTS).RData", sep=""))  # simulation precipitation time series 526
  # load(paste(folder, "sim_Goesdorf_Event7-8.RData", sep=""))
  sim <- sim.det
  # plot.ts(sim[[1]][["out1"]]$C_COD_InTank)
  # plot.ts(sim[[1]][["out1"]]$C_NH4_InTank)
  # plot.ts(sim[[1]][["out1"]]$C_COD_sv)
  # plot.ts(sim[[1]][["out1"]]$V_Tank)
  # 
  # str(mc)
  # mc[["id"]] <- "MC-WQ_2pdf-acorr(NH4s-CODr)-Aimp-P_1500sims_Goesdorf_winter2009-2010"
  # save(mc,file="mc.RData")
  
  #-----------------------------------------------------------------------------------------------------------
  #  define precipitation input data
  #-----------------------------------------------------------------------------------------------------------
  
  # library(xts)
  # P1 <- as.data.frame(index(ts.winter))
  # P1 <- cbind(P1, coredata(ts.winter[,1]))
  # p1 <- P1
  # colnames(p1) <- c("time", "p")
  # head(p1)
  
  a <- IsReg.ts(P1, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  a[[1]]
  ts.P1 <- a[[2]]
  colnames(ts.P1) <- c("p")
  
  p1 <- P1
  colnames(p1) <- c("time", "p")
  
  # event.ini <- "2010-08-14"; event.end <- "2010-08-18"
  # event.ini <- "2010-11-06"; event.end <- "2010-11-14"
  
  # event.ini <- index(ts.winter)[1]; event.end <- index(ts.winter)[nrow(ts.winter)]; ntick <- 70
  # event.ini <- "2009-12-01 "; event.end <- "2009-12-8"; ntick <- 10
  # event.ini <- "2009-12-01 "; event.end <- "2009-12-12"; ntick <- 10
  # event.ini <- index(ts.P1)[1]; event.end <- index(ts.P1)[nrow(ts.P1)]; ntick <- 5000
  # event.ini <- "2009-12-01 "; event.end <- "2009-12-8"; ntick <- 10
  # event.ini <- "2009-12-01 "; event.end <- "2009-12-12"; ntick <- 10
  # event.ini <- "2010-06-16"; event.end <- "2010-06-18"; ntick <- 35
  
  #-----------------------------------------------------------------------------------------------------------
  #  deterministic sim
  #-----------------------------------------------------------------------------------------------------------
  # str(sim)
  # var <- names(sim[[1]][[1]])
  det  <- data.frame(as.POSIXct(sim[[1]][[1]][[2]]))
  det  <- cbind(det, sim[[1]][[1]][[9]])
  det  <- cbind(det, sim[[1]][[1]][[10]])
  # det  <- cbind(det, sim[[1]][[1]][[22]])
  det  <- cbind(det, sim[[1]][[1]][[11]])
  det  <- cbind(det, sim[[1]][[1]][[12]])
  det  <- cbind(det, sim[[1]][[1]][[13]])
  det  <- cbind(det, sim[[1]][[1]][[14]])
  # det.var <- colnames(det) <- c(var[2], var[9], var[10],var[22],var[11],var[12],var[13],var[14])
  # det.var <- colnames(det) <- c("time", "VChamber", "Vsv", "Qsv", "BCODsv", "BNH4sv", "CCODsv", "CNH4sv")
  colnames(det) <- c("time", "VChamber", "Vsv", "BCODsv", "BNH4sv", "CCODsv", "CNH4sv")
  det.var <- c("VChamber", "Vsv", "", "BCODsv", "BNH4sv", "CCODsv", "CNH4sv")
  head(det)
  # plot.ts(det[,"time"], det[,"CCODsv"], typ="l", xlab = "Time", ylab= "CCODsv")
  
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
    save(summ, file="summ.RData")
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
    save(summ, file="summ.RData")
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
  
  # ## Qsv
  # if(summ == NULL){
  #   data.Q <- t(as.data.frame(sim1[,"mcQsv"]))
  #   summ.Q <- MC.summary(p1, data.Q)
  #   summ <- append(summ, list(summ.Q = summ.Q))
  #   save(summ, file="summ.RData")
  #   summ.Q
  # }
  # 
  # 
  # summ.Q.agg <- MC.summary.agg(summ.Q, det, delta, mean, sum)
  # 
  # namePlot.Q <- "Qsv_MC-uncertainty"; ylab1="Q overflow [l/s]"
  # PlotMC.season(summ.Q.agg, paste(namePlot.Q, "_", round(delta/60, 2), "h", sep=""), ylab1, qUpper)
  # 
  # namePlot.QV <- "Vsv_Qsv_MC-uncertainty"; ylab="V overflow [m3]"
  # Q.event <- summ.Q[((as.POSIXct(summ.Q$time) >= event.ini) & (as.POSIXct(summ.Q$time) <= event.end)),]
  # Q.event.agg <- MC.summary.agg(Q.event, det.event, delta, mean, sum)
  # PlotMC.event(V.event.agg, Q.event.agg, 0, det.var[2], det.var[3], paste(namePlot.QV, "_", round(delta/60, 2), "h", sep=""), 
  #              ylab, ylab1, ntick, qUpper)
  # 
  ## removing no needed objects
  # rm(data.VChamber, data.Vol, data.Q)
  rm(data.VChamber, data.Vol)
  
  ##-------------------------------------------------------------------------------------------
  ## CODsv
  ##-------------------------------------------------------------------------------------------
  ## BCODsv
  if(length(summ.data) == 0){
    data.BCOD <- t(as.data.frame(sim1[,"mcBCODsv"]))
    summ.BCOD <- MC.summary(p1, data.BCOD)
    summ <- append(summ, list(summ.BCOD=summ.BCOD))
    save(summ, file="summ.RData")
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
    save(summ, file="summ.RData")
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
    save(summ, file="summ.RData")
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
    save(summ, file="summ.RData")
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
  
  # CNH4.aug <- summ1[summ1$time>="2010-08-09"&summ1$time<="2010-08-19",] 
  # head(CNH4.aug)
  # plot(CNH4.aug$time, CNH4.aug$q95, type="l", col="gray")
  # lines(CNH4.aug$time, CNH4.aug$q50, type="l", col="blue")
  # lines(CNH4.aug$time, CNH4.aug$q05, type="l", col="gray")
  
  ## removing no needed objects
  rm(data.BNH4, data.CNH4)
  
  ##-------------------------------------------------------------------------------------------
  ## saving summary data
  ##-------------------------------------------------------------------------------------------
  # summ <- list(summ.VChamber=summ.VChamber, summ.V = summ.V, summ.Q = summ.Q, 
  #              summ.BCOD=summ.BCOD, summ.BNH4=summ.BNH4, summ.CCOD=summ.CCOD, summ.CNH4=summ.CNH4)
  # save(summ, file="summ.RData")
  # # str(summ)
  
  # ## extracting values
  # load("summ.RData")
  # summ.names <- names(summ)
  # head(summ[[summ.names[1]]])
  # 
  # summ.VChamber <- summ[["summ.VChamber"]]
  # summ.V <- summ[["summ.V"]]
  # summ.Q <- summ[["summ.Q"]]
  # summ.BCOD <- summ[["summ.BCOD"]]
  # summ.CCOD <- summ[["summ.CCOD"]]
  # summ.BNH4 <- summ[["summ.BNH4"]]
  # summ.CNH4 <- summ[["summ.CNH4"]]
  
  # ##============================================================================================================================
  # ##  autocorrelated ts
  # ##============================================================================================================================
  # str(mc)
  # library(xts)
  # head(ts.winter)
  # nrow(ts.winter)
  # CODs <- mc[["par"]]$CODs[1,]
  # length(CODs)
  # ts.winter$CODs <- CODs
  # head(ts.winter)
  # 
  # plot(ts.winter$CODs)
  # acf(ts.winter$CODs)
  # 
  # summ.CODs <- MC.summary(p1, ts.wfunction(x)quantile(x, probs = 0.95)inter$CODs)
  # head(summ.CODs)
  # 
  # summ.CODs.agg <- MC.summary.agg(summ.CODs, det, det.var, delta)
  # 
  # namePlot.CODs <- "CODs_Goesdorf_winter2009-2010"; ylab1="CODs [g/(PE d)]"
  # PlotMC.season(summ.CODs.agg, paste(namePlot.CODs, "_", delta/60, "h", sep=""),  ylab1)
  # 
  # CODs1 <- index(ts.winter)
  # CODs1 <- cbind.data.frame(CODs1, CODs)
  # colnames(CODs1) <- c("tt", "value")
  # head(CODs1)
  # plot(CODs1)
  # CODs.agg <- Agg(CODs1, colnames(CODs1)[2], 24*60, mean, "CODs_Goesdorf_winter2009-2010_2")
  # head(CODs.agg)
  # plot(CODs.agg, typ="l", ylab="CODs [g/(PE d)]")
  # 
  # ## No autocorrelation
  # CODs <- mc[["par"]]$CODs[1,]
  # CODs1 <- mc[["time"]]
  # CODs1 <- cbind.data.frame(CODs1,function(x)quantile(x, probs = 0.95) CODs)
  # class(CODs1[,1][1])
  # CODs1 <- CODs1[CODs1[,1] < as.POSIXct("2010-04-01 00:00:00"),]
  # colnames(CODs1) <- c("tt", "value")
  # head(CODs1)
  # plot(CODs1, typ="l", ylab="CODs [g/(PE d)]")
  # CODs.agg <- Agg(CODs1, colnames(CODs1)[2], 24*60, mean, "CODs_Goesdorf_winter2009-2010_plusMarch")
  # head(CODs.agg)
  # plot(CODs.agg, typ="l", ylab="CODs [g/(PE d)]")
  # acf(CODs.agg[,2])
  ##============================================================================================================================
  ##  variances
  
  ##============================================================================================================================
  
  ##-------------------------------------------------------------------------------------------
  ## computing variances
  ##-------------------------------------------------------------------------------------------
  # c(mean(summ.BCOD$Variance)) 
  # c(mean(summ.CCOD$Variance))
  # c(mean(summ.BNH4$Variance)) 
  # c(mean(summ.CNH4$Variance)) 
  # 
  # c(sum(summ.BCOD$Variance)) function(x)quantile(x, probs = 0.95)
  # c(sum(summ.CCOD$Variance)) 
  # c(sum(summ.BNH4$Variance))
  # c(sum(summ.CNH4$Variance)) 
  # 
  # c(mean(summ.BCOD$Mean)) 
  # c(mean(summ.BNH4$Mean)) 
  # c(mean(summ.CCOD$Mean)) 
  # c(mean(summ.CNH4$Mean)) 
  # 
  # c(sum(summ.BCOD$Mean)) 
  # c(sum(summ.BNH4$Mean)) 
  # c(sum(summ.CCOD$Mean)) 
  # c(sum(summ.CNH4$Mean))
  # 
  # c(mean(summ.CCOD$Max)) 
  # c(mean(summ.CCOD$q95)) 
  
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
  # write.csv(names, file="names.csv")
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
  # variance.CODs.NH4s.CODr <- variance
  # save(variance.CODs.NH4s.CODr, file="variance_CODs-NH4s-CODr.RData")
  # 
  # variance.CODs.NH4s <- variance
  # save(variance.CODs.NH4s, file="variance_CODs-NH4s.RData")
  # 
  # variance.CODs.CODr <- variance
  # save(variance.CODs.CODr, file="variance_CODs-CODr.RData")
  # 
  # variance.NH4s.CODr <- variance
  # save(variance.NH4s.CODr, file="variance_NH4s-CODr.RData")
  # 
  # det.var
  # mean(det[,"CCODsv"])
  # max(det[,"CCODsv"])
  
  save(variance, file="variance.RData")
  
  # # ##-------------------------------------------------------------------------------------------
  # # ## NH4s
  # # ##-------------------------------------------------------------------------------------------
  # # class(mc)
  # # names(mc)
  # # data.NH4s <- mc$par$NH4s
  # # dim(dataNH4s)
  # # 
  # # summ.NH4s <- MC.summary(p1, data.NH4s)
  # # head(summ.NH4s)
  # # summ.NH4s.agg <- MC.summary.agg(summ.NH4s, 6*60)
  # # head(summ.NH4s.agg)
  # # 
  # # namePlot.NH4s <- "MC_NH4s_Confidence-intervals_Goesdorf_Esch-Sure2010"; ylab="NH4s [g/PE/d]"
  # # PlotMC.season(summ.NH4s.agg, namePlot.NH4s)
  # # 
  # # NH4s.event <- summ.NH4s[((summ.NH4s$time >= "2010-08-09") & (summ.NH4s$time <= "2010-08-19")),]
  # # head(NH4s.event)
  # # PlotMC.event(summ= NH4s.event, obs=0, namePlot=namePlot, ylab=ylab)
  # # 
  # # NH4s.event.agg <- MC.summary.agg(NH4s.event, 6*60)
  # # PlotMC.event(summ= NH4s.event.agg, obs=0, namePlot=namePlot, ylab=ylab)
  # # 
  # # ##-------------------------------------------------------------------------------------------
  # # ## volume analysis
  # # ##-------------------------------------------------------------------------------------------
  # # head(mctV)
  # # dim(mctV)
  # # meanV <- apply(mctV,1,FUN = mean)
  # # length(meanV)
  # # hist(meanV, breaks=20, freq=F)
  # # lines(density(meanV))
  # # plot(mctV[1,], type = "l")
  
  print("End MC analysis. Please check your output folder.")
  
  return(list(summ=summ, variance=variance))
}
#)