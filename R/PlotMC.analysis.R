# Monte-Carlo simulations analysis (plot)
# author: J.A. Torres-Matallana
# organization: Luxembourg Institute of Science and Technology (LIST), Luxembourg
#               Wagenigen University and Research Centre (WUR), Wageningen, The Netherlands   
# date: 17.10.2015 - 16.05.2016

PlotMC.event <- function(summ, summ1, obs, det.var, det.var1, namePlot, ylab, ylab1, ntick, qUpper){
  #summ <- CNH4.event
  #head(summ)
  # dev.new(width=10, height=10)
  pdf(paste0(namePlot, "-event.pdf"), width=10, height=10)
  par(mfrow = c(3, 1))
  par(mar = c(2, 7, 3, 6), cex.axis=1.5, cex.lab=1.5, cex.main=2) # bottom, left, top and right margins 
  ymax <- max(summ["p1"]); ymin <- min(summ["p1"]) 
  plot(summ[,"time"], summ[,"p1"], type="l", col="blue", ylab = "", xaxt="n", ylim=c(ymax, ymin), main=namePlot) # rainfall
  mtext("Precipitation [mm]", side=2, line=3.5, at=(ymax-ymin)/2, cex=1.)
  #axis(side=1, at=summ[,"time"], labels=as.Date(summ[,"time"], format = "%Y %M "))
  
  par(mar = c(4, 7, 1, 6), cex.axis=1.5, cex.lab=1.5, cex.main=2)
  x <- 1:dim(summ)[1]
  ylim <- c(0, max(max(summ[[qUpper]]), max(summ["Mean"], summ[det.var])))
  plot(x, t(summ["q05"]), typ="l", ylim=ylim, col="grey83", xlab="", ylab=ylab, xaxt="n")
  lines(x, t(summ[[qUpper]]), typ="l", col="grey83")
  polygon(c(x, rev(x)), c(t(summ[[qUpper]]), rev(t(summ["q05"]))), col = "grey83", border = NA)
  #polygon(c(x, rev(x)), c(t(summ["q75"]), rev(t(summ["q25"]))), col = "grey50", border = NA)
  lines((summ["Mean"]), typ = "l", col="black", lty="dashed", lwd=2)
  lines((summ[det.var]), typ = "l", col="blue", lty="solid", lwd=2)
  
  legend("topright", legend = c("Deterministic simulation",
                                "Mean MC simulation", "Confidence interval"),
         lty=c(1,2,1), col=c("blue", "black", "grey60"), lwd = c(1,1,10), cex=1.5, bty = "n") # "bottomright" "bottomleft"
  # legend("top", c(paste("sum mean =", floor(sum(summ["Mean"],na.rm = T)), sep=" "), 
  #        paste("sum det =", floor(sum(summ[det.var],na.rm = T)), sep=" "), 
  #        paste("diff =", floor(sum(summ["Mean"],na.rm = T) - sum(summ[det.var],na.rm = T)), sep=" ")), 
  #        lty=NULL, col=c("white"), cex=1.5, bty = "n") # "bottomright" "bottomleft"
  
  par(mar = c(5, 7, 0, 6), cex.axis=1.5, cex.lab=1.5, cex.main=2)
  x <- 1:dim(summ1)[1]
  ylim <- c(0, max(max(summ1[[qUpper]]), max(summ1["Mean"], summ1[det.var1])))
  plot(x, t(summ1["q05"]), typ="l", ylim=ylim, col="grey83", xlab="Time [month-day hour:min]", ylab=ylab1, xaxt="n")
  lines(x, t(summ1[[qUpper]]), typ="l", col="grey83")
  polygon(c(x, rev(x)), c(t(summ1[[qUpper]]), rev(t(summ1["q05"]))), col = "grey83", border = NA)
  #polygon(c(x, rev(x)), c(t(summ1["q75"]), rev(t(summ1["q25"]))), col = "grey50", border = NA)
  lines((summ1["Mean"]), typ = "l", col="black", lty="dashed", lwd=2)
  lines((summ[det.var1]), typ = "l", col="blue", lty="solid", lwd=2)
  legend("topright", legend = c("Deterministic simulation",
                                "Mean MC simulation", "Confidence interval"),
         lty=c(1,2,1), col=c("blue", "black", "grey60"), lwd = c(1,1,10), cex=1.5, bty = "n") # "bottomright" "bottomleft"
  # legend("top", c(paste("sum mean =", floor(sum(summ1["Mean"],na.rm = T)), sep=" "), 
  #                 paste("sum det =", floor(sum(summ1[det.var1],na.rm = T)), sep=" "), 
  #                 paste("diff =", floor(sum(summ1["Mean"],na.rm = T) - sum(summ1[det.var1],na.rm = T)), sep=" ")), 
  #        lty=NULL, col=c("white"), cex=1.5, bty = "n") # "bottomright" "bottomleft"
  
  #summ1 <- CCOD.event.agg
  #head(summ1)
  #min(summ1[,"time"]); max(summ1[,"time"])
  ix <- seq(1, length(x), ntick)
  labels <- strftime(summ1[,"time"], format = "%m-%d %H:%M")
  axis(side=1, at=x[ix], labels=labels[ix], cex=1.5, cex.axis=1.5)
  
  if(obs == 0){
    #  legend("topright", legend=c("q5-q95 confidence interval", "q25-q75 condifence interval", "Mean"), 
    #        fill=c("grey90", "grey70", "blue", "blue"), cex=2)
  }else{
    #       lines(obs, typ = "p", col="blue", pch=22, bg="blue", cex=1.0)
    #       legend("topright", legend=c("q5-q95 confidence interval", "q25-q75 condifence interval",
    #                                   "observed", "Mean simulated") , pch=c(NA,NA,22,NA), lty = c(1,1,NA, 2), lwd=c(20,20,1,2),
    #              col=c("grey90", "grey70", "blue", "red"), cex=2)  
  }
  
  # printing plot
  # dev.copy2pdf(file= paste(namePlot, "-event.pdf", sep=""))
  #dev.copy(png, filename=paste(namePlot, ".png", sep=""), res=600, height=10, width=10, units="in")
  dev.off()
}

PlotMC.season <- function(summ1, namePlot, ylab, qUpper){
  #summ1 <- summ.agg
  #summ1 <- summ.agg.NH4
  #summ1 <- summ.VChamber.agg1
  
  
  my.panel.bands <- function(x, y, upper, lower, fill, col,
                             subscripts, ..., font, fontface)
  {
    upper <- upper[subscripts]
    lower <- lower[subscripts]
    panel.polygon(c(x, rev(x)), c(upper, rev(lower)),
                  col = fill, border = FALSE,
                  ...)
  }
  
  lev <- (factor(summ1$month))
  nlev <- levels(lev)
  
  if(length(nlev==13)){
    summ1 <- summ1[2:nrow(summ1),]
    lev <- (factor(summ1$month))
    nlev <- levels(lev)
    # length(nlev)
  }
  
  requireNamespace("lattice")
  # dev.new(type="Xlib",width=16, height=10)
  pdf(file= paste0(namePlot, "-season.pdf"), width=10, height=5)
  par(mar = c(1, 1, 1, 1))
  #dev.new(width=16, height=10)
  #pdf(paste(namePlot, ".pdf", sep=""), pointsize=10)
  
  ylim <- c(min(min(summ1$q05, summ1$Mean), min(summ1[qUpper])), max(max(summ1$q05, summ1$Mean), max(summ1[qUpper]))*1.1)
  if(length(nlev) >= 12){lay1 <- 3; lay2 <- 4; index <- c(10,11,12,7,8,9,4,5,6,1,2,3) }else
    if(length(nlev) == 1){lay1 <- 1; lay2 <- 1; index <- c(1)}else
      if(length(nlev) == 2){lay1 <- 1; lay2 <- 2; index <- c(2,1)}else
      {lay1 <- 1; lay2 <- 3; index <- c(2,1,3)}
  
  print(xyplot(summ1[[qUpper]]~(summ1$time)|lev, data = summ1, scales=list(x="free", cex=1.1,
         y=list(cex=1.1)), ylim=ylim, groups = lev, xlab = list(label="Time", cex=1.9), ylab=list(label=ylab, cex=1.9),
         layout =c(lay1,lay2), index.cond=list(index),
         upper = summ1[[qUpper]], lower = summ1$q05,
         panel = function(x, y, ...){
           panel.superpose(x, y, panel.groups = my.panel.bands, type='l', fill='gray',col='gray', ...)
           panel.xyplot(x, y, type='l', cex=1.5, lty=1, col='gray',...)
           panel.xyplot(summ1$time, summ1$q05, type='l', cex=0.6, lty=1, col='gray',...)
           panel.xyplot(summ1$time, summ1$Mean, type='l', cex=0.6, lty=1, col='blue',...)
           #panel.xyplot(summ1$time, summ1$BNH4ov, type='l', cex=0.6, lty=1, col='black',...)
           }))   
  
  # dev.copy2pdf(out.type = "pdf", file= paste(namePlot, "-season.pdf", sep="")) 
  dev.off()
}

# PlotMC.variable <- function(summ, ylab, namePlot, q05, qUpper){
#   # summ: data.frame
#   # ylab: character
#   # namePlot: character
#   # summ <- p1.summ9
#   # requireNamespace("ggplot2")
#   # requireNamespace("lattice")
#   # dev.new(type="Xlib",width=16, height=10)
#   pdf(paste0(namePlot, ".pdf"),width=10, height=10)
#   par(mar = c(1, 1, 1, 1))
#   
#   myplot <- ggplot(summ, aes(time, summ[,"Mean"]))+
#     geom_line(data=summ[,c("time","Mean")], col="blue")+
#     geom_ribbon(data=summ,aes(ymin=q05,ymax=qUpper),alpha=0.3)+
#     ylab(ylab)
#   
#   print(myplot + theme(axis.text.x= element_text(size=22), axis.text.y= element_text(size=22)))
#   print(myplot + theme_bw())
#         
#   # dev.copy2pdf(out.type = "pdf", file= paste(namePlot, ".pdf", sep="")) 
#   dev.off()
# }

##-----------------------------------------------------------------------------------------------------------
## plus minus two standard errors against iterations for a single sequence of simulations
##-----------------------------------------------------------------------------------------------------------
# NN <- function(x) {
#   if(obs == 0 ){
#     0
#   }else{
#     of <- "ME"
#     of <- "RMSE"
#     of <- "NSE"
#     
#     idCrit <- 4
#     of     <- cstr[idCrit]
#     
#     
#     n <- dim(pars)[1]
#     
#     
#     # x <- critn[,of]
#     x <- crit[,of]
#     
#     min(x); mean(x); max(x)
#     head(x)
#     
#     estint=cumsum(x) / ( 1 : n)
#     esterr=sqrt ( cumsum( (x-estint) ^2) ) / ( 1 : n)
#     
#     dev.new(width=10, height=10)
#     par(mar = c(5, 5, 2, 2), cex.axis=2, cex.lab=2, cex.main=2)
#     plot ( estint , main="Mean and error range" , xlab="Number of simulations", ylab=of, type="l" , lwd=
#              + 2 , ylim=mean(x) +50*c ( -esterr [n] , esterr [n] ), xlim=c(0,1000))
#     lines ( estint+2*esterr, col="grey60" , lwd=2, lty="dashed")
#     lines ( estint-2*esterr, col="grey60" , lwd=2, lty="dashed")
#     legend("topright", legend=c("Cumulative sum", "Mean +- 2 standard error") , lty=c(1,2), 
#            col=c("black", "grey50"), cex=2)
#     
#     # printing plot
#     dev.copy2pdf(file= paste("MC_error-range_",of, ".pdf", sep=""))
#     dev.copy(png, filename=paste("MC_error-range_",of, ".png", sep=""), res=600, height=10, width=10, units="in"); dev.off()
#     
#   }
#   ##-----------------------------------------------------------------------------------------------------------
#   ## parameter dotty plots
#   ##-----------------------------------------------------------------------------------------------------------
#   if(obs == 0){
#     0
#   }else{
#     
#     dev.new(width=10, height=10)
#     par(mfrow = c(4, 3))
#     par(mar = c(5, 5, 2, 2), cex.axis=2, cex.lab=2, cex.main=2)
#     dim(pars)
#     pstr
#     for (i in 1:ncol(pars)){
#       plot(pars[,i], crit[,idCrit], xlab=pstr[i],ylab=cstr[idCrit])  
#     }
#     
#     # printing plot
#     dev.copy2pdf(file= paste("MC_dottyPlots_",of, ".pdf", sep=""))
#     dev.copy(png, filename=paste("MC_dottyPlots_",of, ".png", sep=""), res=600, height=10, width=10, units="in"); dev.off()
#     ##-----------------------------------------------------------------------------------------------------------
#     ## parameter histograms
#     ##-----------------------------------------------------------------------------------------------------------
#     dev.new(width=10, height=10)
#     par(mfrow = c(4, 3))
#     par(mar = c(5, 5, 2, 2), cex.axis=2, cex.lab=2, cex.main=2)
#     dim(pars)
#     pstr
#     for (i in 1:ncol(pars)){
#       hist(pars[,i], xlab=pstr[i], main="")  
#     }
#     
#     # printing plot
#     dev.copy2pdf(file= paste("MC_parHisto_",of, ".pdf", sep=""))
#     dev.copy(png, filename=paste("MC_parHisto_",of, ".png", sep=""), res=600, height=10, width=10, units="in"); dev.off()
#   }
#   return(list(summ, summ1))
# }
