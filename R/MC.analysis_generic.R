# Monte-Carlo simulations analysis (generic case)
# author: J.A. Torres-Matallana
# organization: Wagenigen University and Research Centre (WUR), Wageningen, The Netherlands
# date: 13.02.2019 - 13.02.2019


# setGeneric("MC.analysis", function(x, delta, qUpper, p1.det, sim.det, event.ini, event.end, ntick, summ) standardGeneric("MC.analysis"))
# setMethod("MC.analysis", signature = c("list", "numeric", "character", "data.frame", "list", "POSIXct", "POSIXct", "numeric", "list"), 

MC.analysis_generic <- function(x, delta, qUpper, data.det, sim.det, 
                        event.ini, event.end, ntick, summ.data = NULL){
  
  sim1 <- x
  
  P1 <- data.det
  sum(P1[,2])
  
  #-----------------------------------------------------------------------------------------------------------
  #  define generic (precipitation) input data
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
  det <- sim.det[[1]]
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
    data.1 <- t(as.data.frame(sim1[[1]][,2:ncol(sim1[[1]])]))
    summ.1 <- MC.summary(p1, data.1)
    summ <- list(summ.1 = summ.1)
    # save(summ, file="summ.RData") # uncomment to save file
    summ.1
  },{
    summ.1 <- summ[["summ.1"]]
  })
  
  summ.1.agg <- MC.summary.agg(summ.1, det, delta, mean, sum)
  
  namePlot.1 <- "Output1_MC-uncertainty"; ylab="Output 1 [units]"
  
  ## not run
  ## creating the plot (uncomment to run)
  # PlotMC.season(summ.1.agg, paste(namePlot.1, "_", round(delta, 0), "min", sep=""), ylab, qUpper)
  
  output1.event <- summ.1[((as.POSIXct(summ.1$time) >= event.ini) & (as.POSIXct(summ.1$time) <= event.end)),]
  output1.event.agg <- MC.summary.agg(output1.event, det.event, delta, mean, sum)
  
  ## not run
  ## creating the plot (uncomment to run)
  # PlotMC.event(summ = output1.event.agg, summ1 = output1.event.agg, obs = 0, det.var = "output.det", det.var1 = "output.det", 
  #              namePlot = paste(namePlot.1, "_", round(delta, 0), "min", sep=""), 
  #              ylab = ylab, ylab1 = ylab, ntick = ntick, qUpper = qUpper)
   
  rm(data.1)
  
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