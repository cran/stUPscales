# plot Eval data.frame
# author: J.A. Torres-Matallana
# organization: Luxembourg Institute of Science and Technology (LIST), Luxembourg
#               Wagenigen University and Research Centre (WUR), Wageningen, The Netherlands   
# date: 14.07.2015 - 14.07.2015

PlotEval <- function(eval, ts, gof1, namePlot, pos1, pos2, pos3){
            
  #=================================================================================
  # change languague to english for months in plot.zoo
  #=================================================================================
  Sys.setlocale("LC_TIME", "C") 
  
  #=================================================================================
  # plot.zoo
  #=================================================================================
  # dev.new(width=30, height=30)
  dev.new(paste0(namePlot, ".pdf"), width=30, height=30)
  par(mfrow = c(3, 1))
  # library(xts)
  
  roundUp <- function(x) round(x+5,-1)
  
  # ts rainfall
  par(mar = c(0, 6, 3, 6), cex.axis=2, cex.lab=2, cex.main=2) # bottom, left, top and right margins 
  ymax <- max(ts[,5]); ymin <- min(ts[,5]) 
  plot(eval[,"time"], ts[,"Rainfall"], type="l", col="blue", xlab = "time", ylab = "", xaxt="n", main=namePlot, ylim=c(ymax, ymin)) # rainfall
  mtext("Precipitation [mm]", side=2, line=3.5, at=(ymax-ymin)/2, cex=1.5)
  
  # ts sim vs obs
  gof_ME   <- round(gof1[1], digits=3)
  gof_RMSE <- round(gof1[4], digits=3)
  gof_NSE  <- round(gof1[9], digits=3)
  gof_R2   <- round(gof1[17], digits=3)
  gof_rNSE   <- round(gof1[11], digits=3)
  
  par(mar = c(0, 6, 0, 6), cex.axis=2, cex.lab=2, cex.main=2)
  head(ts)
  ymax <- roundUp(max(ts[,3], ts[,2])); ymin <- -roundUp(abs(min(ts[,3], ts[,2]))) 
  plot(eval[,"time"], ts[,3], type="b", lty="dashed", col="red", xlab = "time", ylab = "", xaxt="n", yaxt="n", 
           ylim= c(ymin, ymax), cex=0.5) # simulated
  lines(eval[,"time"], ts[,2], type="b", cex=.5) # observed
  axis(side=4, at=c(0,ymax/4, 2*ymax/4, 3*ymax/4, ymax))
  mtext("Volume, tank [m3]", side=4, line=3.5, at=(ymax-ymin)/2, cex=1.5)
  legend(pos2, c("simulated", "observed", paste("ME       =", gof_ME, sep=" "),
                       paste("RMSE =", gof_RMSE, sep=" "), paste("NSE    =", gof_NSE, sep=" ")),
         lty=c(2,1), pch=c(1), 
         col=c("red","black", "white", "white", "white"), cex=2)
  
  # ts difference
  par(mar = c(5, 6, 0, 6), cex.axis=2, cex.lab=2, cex.main=2)
  #ymax <- roundUp(max(ts[,4])); 
  # ymin <- -roundUp(abs(min(ts[,4]))) 
  plot(eval[,"time"], ts[,4],pch=0,type="b", col="black", xlab="Time", ylab="", ylim=c(-ymax, ymax), yaxt="n", cex=0.5)
  abline(a=0, b=0, lty=2, col="blue")
  axis(side=2, at=c(ymin, 0, ymax), labels=c(as.character(ymin), "0", as.character(ymax)))
  mtext("Diff [m3]", side=2, line=3.5, at=(ymax-abs(ymin))/2, cex=1.5)
  legend(pos3, c("zero line", "Diff = sim. - obs."),
         lty=c(2,1), pch=c(NA,0), col=c("blue", "black"), cex=2.0)
  
  # dev.copy2pdf(file= paste(namePlot, ".pdf", sep=""))
  # dev.copy(png, filename=paste(namePlot, ".png", sep=""), res=600, height=10, width=10, units="in"); dev.off()
  dev.off()
}