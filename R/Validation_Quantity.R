# model EmiStat-R, validation for water quantity
# author: J.A. Torres-Matallana
# organization: Luxembourg Institute of Science and Technology (LIST), Luxembourg
#               Wagenigen University and Research Centre (WUR), Wageningen, The Netherlands   
# date: 13.04.2017 - 13.04.2017

setGeneric("Validation_Quantity", function(x, y) standardGeneric("Validation_Quantity"))

setMethod("Validation_Quantity", signature = c("input", "inputObs"), 
          
          function(x, y){
            # x <- input.user
            # y <- inputObs(id = 1, plot = 1, delta = delta, observations = observations, lev2vol = lev2vol,
            #          namePlot = "Goesdorf Event 7-8", legendPosition = legendPosition)            
            #================================================================================
            # Rain data
            #================================================================================
            # for extracting windows of ts see: EmiStat-R_valida_extractRainfall.R
            # setwd("/home/atorres/Documents/02_working/06_LIST-Wageningen_PhD/18_models/01_EmiStat-R/R/EmiStat-R_input/input")
            
            # load("P1.RData")
            # str(x)
            # str(y)
            
            
            id <- slot(y, "id")
            
            P1 <- slot(x, "P1")
            # head(P1); tail(P1)
            # class(P1)
            # dim(P1)
            # plot(P1["time"], P1["rainfall"], type="l")
            # plot(P1["time"][1:1441,], P1["rainfall"][1:1441,], type="l")
            # plot(P1[,1:2])
            #---------------------------------------------------------------------------------
            # checking regularity of TS
            #---------------------------------------------------------------------------------
            a   <- IsReg.ts(data = P1, format="%Y-%m-%d %H:%M:%S", tz="UTC")
            # str(a)
            a[1]
            deltat(a[1])
            ts <- a[[2]]
            # class(ts)
            # head(ts)
            
            #---------------------------------------------------------------------------------
            # aggregation of rainfall
            #---------------------------------------------------------------------------------
            delta_P1 <- slot(y, "delta")$P1
            namePlot <- slot(y, "namePlot")
            # P1 <- Agg(P1, "P1", delta_P1, "sum", namePlot)  ## modify for spatial precipitation
            # slot(x, "P1") <- P1
            
            #=======================================================================================================
            # run EmiStat-R
            #=======================================================================================================
            sim <- EmiStatR(x)
            #str(sim)
            #summary(sim[[1]]$out1)
            
            #---------------------------------------------------------------------------------
            # plotting EmiStat-R output
            #---------------------------------------------------------------------------------
            # str(sim)
            #data <- read.csv(fileNameOut, header = TRUE, sep = ",", quote = "\"", dec = ".", fill = TRUE, comment.char = "")
            data <- sim[[id]]["out1"]$out1
            data1 <- data[,2:dim(data)[2]]
            # summary(data)
            # head(data1,10)
            # fileNameOut <-  slot(y, "id")
            # if(plots == 1){
            #   PlotCol(data=data1, fileName=paste("plot1_", fileNameOut, sep=""), format="%Y-%m-%d %H:%M:%S", tz="UTC", plot=c(1,7,8))
            #   PlotCol(data=data1, fileName=paste("plot2_",fileNameOut, sep=""), format="%Y-%m-%d %H:%M:%S", tz="UTC", plot=c(1,11,12))
            # }
            
            #=======================================================================================================
            # Goodness-of-fit volume tank
            #=======================================================================================================
            #--------------------------------
            # loading observed data
            #--------------------------------
            # for individual extraction see: EmiStat-R_valida_extractWaterLevelTank.R
            # setwd("/home/atorres/Documents/02_working/06_LIST-Wageningen_PhD/18_models/01_EmiStat-R/R/EmiStat-R_valida")
            # load(nameWLT_obs);
            # wlt_obs <- obs
            
            wlt_obs <- slot(y, "observations")$WLT[,c(1,4)]
            #--------------------------------
            # aggregation of volume
            #---------------------------------
            # # source("/home/atorres/Documents/02_working/06_LIST-Wageningen_PhD/18_models/01_EmiStat-R/R/agg.R")
            # delta_wlt_obs <- slot(y, "delta")$wlt_obs
            # wlt_obs <- Agg(wlt_obs, "wlt_obs", delta_wlt_obs, "mean", namePlot)

            colnames(wlt_obs) <- c("time", "wlt_obs")

            #-------------------------------------------------------------------------------------------------------
            # conversion level to volume
            #-------------------------------------------------------------------------------------------------------
            # source("/home/atorres/Documents/02_working/06_LIST-Wageningen_PhD/18_models/01_EmiStat-R/R/tanks.R")
            # head(wlt_obs)
            # levels <- wlt_obs[4]
            # # head(levels)
            # # tail(levels)
            # # class(levels)
            # # dim(levels)
            # 
            # volume <- apply(levels, 2, tankGOE)
            # # head(volume)
            # # dim(volume)
            # # class(volume)
            # # str(volume)
            # 
            # levels <- wlt_obs[1]
            # levels <- cbind.data.frame(levels, wlt_obs[4])
            # levels[,3] <- volume
            # # dim(levels)
            # colnames(levels) <- c("time", "levelT_obs")
            # colnames(levels) <- c("time", "levelT_obs", "volumeT_obs")
            # colnames(levels[,3]) <- c("volumeT_obs")
            # 
            # names(levels)[1] <- "time"
            # names(levels)[2] <- "levelT_obs"
            # names(levels)[3] <- "volumeT_obs"
            # 
            # plot(levels[,1],levels[,3], typ="l")
            # 
            # # head(levels); tail(levels)
            # # class(levels)
            # levels_obs <- levels
            
            #------------------------
            # 1, 10, 30, 60 minutes resolution
            #------------------------
            #source("/home/atorres/Documents/02_working/06_LIST-Wageningen_PhD/18_models/01_EmiStat-R/R/tanks.R")
            if(slot(y, "var") == "V_Chamber [m3]"){
              levels_obs <- wlt_obs[2]
              volume     <- apply(levels_obs/100, 2, Level2Volume, lev2vol = slot(y, "lev2vol"))
              levels_obs <- wlt_obs[1]
              levels_obs <- cbind.data.frame(levels_obs, wlt_obs[2])
              levels_obs[,3] <- volume
              
              # head(levels_obs)
              colnames(levels_obs) <- c("time", "level-tank_obs")
              colnames(levels_obs[,3]) <- c("volume-tank_obs")
              # names(levels_obs)[3] <-   "volume-tank_obs"
            }else{
              levels_obs <- cbind.data.frame(wlt_obs, wlt_obs[2])
              colnames(levels_obs) <- c("time", "obs", "obs1")
            }
            
            
            # head(levels_obs)
            # tail(levels_obs)
            
            # par(mfrow = c(2,1))
            # plot.zoo(levels_obs[,1],levels_obs[,2], typ="l")
            # plot.zoo(levels_obs[,1],levels_obs[,3], typ="l")
            # head(levels)
            # tail(levels)
            # class(levels)
            
            #-------------------------------------------------------------------------------------------------------
            # output tank filling-up volume
            #-------------------------------------------------------------------------------------------------------
            # setwd("/home/atorres/Documents/02_working/06_LIST-Wageningen_PhD/18_models/01_EmiStat-R/R/EmiStat-R_output")
            #data <- read.csv(fileNameOut, header = TRUE, sep = ",", quote = "\"", dec = ".", fill = TRUE, comment.char = "")
            
            # str(data)
            # head(data)
            # class(data)
            vol_sim <- as.data.frame(data[,2])
            # vol_sim <- cbind.data.frame(vol_sim, as.data.frame(data[,9])) # V_Chamber
            # vol_sim <- cbind.data.frame(vol_sim, as.data.frame(data[,17])) # V_InTank
            vol_sim <- cbind.data.frame(vol_sim, as.data.frame(data[,slot(y, "var")]))
            colnames(vol_sim) <- c("time", "value")
            # head(vol_sim); tail(vol_sim)
            # class(vol_sim)
            # dim(vol_sim)
            # plot.zoo(vol_sim[,1], vol_sim[,2], type="l")
            
            #---------------------------------
            # converting to 10, 30, 60 min resolution
            #---------------------------------
            # # source("/home/atorres/Documents/02_working/06_LIST-Wageningen_PhD/18_models/01_EmiStat-R/R/agg.R")
            # delta_vol_sim <- slot(y, "delta")$vol_sim
            # vol_sim_agg <- Agg(vol_sim, "vol_sim", delta_vol_sim, "mean", namePlot)  
            # colnames(vol_sim_agg) <- c("time", "volT_sim")
            
            colnames(vol_sim) <- c("time", "volT_sim")
            # vol_sim <- vol_sim_agg
          
            #-------------------------------------------------------------------------------------------------------
            # Goodness-of-fit
            #-------------------------------------------------------------------------------------------------------
            # head(levels_obs); tail(levels_obs)
            # head(vol_sim); tail(vol_sim)
            # dim(levels_obs)
            # dim(vol_sim)
            # class(levels_obs)
            # class(vol_sim)
            # head(vol_sim); tail(vol_sim)
            # ?merge
            # write.table(levels_obs, file = paste("test", "_levels_obs.csv", sep=""), sep = ",", qmethod = "double", row.names=FALSE)
            # write.table(vol_sim, file = paste("test", "_vol_sim.csv", sep=""), sep = ",", qmethod = "double", row.names=FALSE)
            # dim(levels_obs)
            eval <- merge(x=levels_obs, y=vol_sim, by="time")
            # write.table(eval, file = paste("test", "_eval-merge.csv", sep=""), sep = ",", qmethod = "double", row.names=FALSE)
            if(dim(eval)[1] > 1){}else{
              # for 1, 10,30, 60 min resolution
              eval <- cbind.data.frame(levels_obs, vol_sim[-1])
              print("done line 214")
            }
            
            # head(eval, 10); tail(eval,10)
            # dim(eval)  
            
            # logarithms
            # eval[,3] <- log1p(eval[,3])
            # eval[,4] <- log1p(eval[,4])
            
            # differences
            eval[,5] <- eval[,4]- eval[,3]
            # colnames(eval[,5]) <- "diff" ## <----------------- check if is it necessary this line
            # class(eval)
            # head(eval)
            
            
            #   # load previous result
            #   setwd("/home/atorres/Documents/02_working/06_LIST-Wageningen_PhD/18_models/01_EmiStat-R/R/EmiStat-R_output")
            #   load("eval.RData")
            #   load("gof.RData")
            
            #   merge with P1
            #   head(eval); tail(eval)
            #   head(P1); tail(P1)
            #   class(P1); class(eval)
            #   dim(eval); dim(P1)
            
            eval1 <- merge(x=eval, y=P1, by="time")
            if(dim(eval1)[1] > 1){}else{
              # for 1, 10,30, 60 min resolution
              eval1 <- cbind.data.frame(eval, P1[c(-1,-3)])
              print("done line 242")
            }
            # dim(P1)
            # dim(eval1)
            eval <- as.data.frame(eval1)
            
            
            eval[,3] <- eval[,3][c(1:dim(eval)[1])]
            eval[,5] <- eval[,5][c(1:dim(eval)[1])]
            colnames(eval)[3] <- "volT_obs"
            colnames(eval)[5] <- "diff"
            
            write.csv(eval, file = paste("evalVolT.csv", sep=""))
            
            # save(eval, file=paste("evalVolT.RData", sep=""))
            
            # eval to time series
            # source("/home/atorres/Documents/02_working/06_LIST-Wageningen_PhD/10_LIST-data/HauteSureData/R/isReg.R")
            a   <- IsReg.ts(data = eval, format="%Y-%m-%d %H:%M:%S", tz="UTC")
            ts <- a[[2]]
            # head(ts)
            
            # GoF
            # source("/home/atorres/Documents/02_working/06_LIST-Wageningen_PhD/18_models/01_EmiStat-R/R/gof1.R")
            gof1 <- GoF(eval, 4, 3, "VolT") # 4 = sim; 3 = obs
            # NSE(eval[,4], eval[,3])
            
            #-------------------------------------------------------------------------------------------------------
            # Goodness-of-fit (plotting)
            #-------------------------------------------------------------------------------------------------------
            # source("/home/atorres/Documents/02_working/06_LIST-Wageningen_PhD/18_models/01_EmiStat-R/R/plotEval.R")
            # namePlot <- paste(namePlot, "(res =", delta_P1, "-",delta_wlt_obs,"-",delta_vol_sim, "min)", sep=" ")
            namePlot <- paste(namePlot, "(res =", delta_P1, "_min)", sep=" ")
            
            # str(y)
            pos1  <- slot(y,"legendPosition")[[1]]$pos1
            pos2  <- slot(y,"legendPosition")[[1]]$pos2
            pos3  <- slot(y,"legendPosition")[[1]]$pos3
            plots <- slot(y,"plot")[[1]]
            if(plots == 1){
              colnames(eval)[6] <- "Rainfall"
              colnames(ts)[5]   <- "Rainfall"
              PlotEval(eval, ts, gof1, namePlot, pos1, pos2, pos3)  
            }
           
            #=======================================================================================================
            # Output
            #=======================================================================================================
            #-------------------------------------------------------------------------------------------------------
            # assembling output
            #-------------------------------------------------------------------------------------------------------
            # head(NH4_sim); dim(NH4_sim)
            # head(COD_sim); dim(COD_sim)
            # head(vol_sim); dim(vol_sim)
            #simu      <- c(vol_sim[,2], COD_sim[,2], NH4_sim[,2])
            simu       <- cbind.data.frame(VolT = vol_sim[,2])
            #sim_time  <- c(vol_sim[,1], COD_sim[,1], NH4_sim[,1])
            sim_time   <- vol_sim[,1]
            idd        <- 1:dim(simu)[1]
            # sim      <- cbind.data.frame(x=id, times=sim_time, sim=simu)
            sim        <- cbind.data.frame(time_step=idd, simu)
            # var1     <- rep("volT_[m3]", length(vol_sim[,2]))
            # var2     <- rep("COD_[mg/l]", length(COD_sim[,2]))
            # var3     <- rep("NH4_[mg/l]", length(NH4_sim[,2]))
            # var      <- c(var1, var2, var3)
            # #sim     <- cbind.data.frame(time=sim_time, sim=simu, var=var)
            output     <- list(sim=sim, sim_time=sim_time, evalVolT=eval, of.NSE=gof1[9])    
            
            # length(sim)/3
            
            # write.table(output, file = paste("output.csv", sep=""), 
            #             sep = ",", qmethod = "double", row.names=FALSE, col.names=TRUE)
            #save(output, file=paste("output.RData", sep=""))
            return(output)
          }
)
