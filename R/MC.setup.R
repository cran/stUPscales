# setup for Monte-Carlo runs with pdfs definition (formerly MC_WQ_pdf_pre.R)
# author: J.A. Torres-Matallana
# organization: Luxembourg Institute of Science and Technology (LIST), Luxembourg
#               Wagenigen University and Research Centre (WUR), Wageningen, The Netherlands   
# date: 22.04.2016 - 28.07.2017


setGeneric("MC.setup", function(x) standardGeneric("MC.setup"))

setMethod("MC.setup", signature = "setup", 
          
          function(x){
            # library(EmiStatR)
            # data("Esch_Sure2010")
            # id       = "MC_calibra1",
            # nsim     = 11,  # for discrete sampling 'dis' nsim is overrided
            # seed     = 123,
            # mcCores  = 4,
            # ts.input = RT,
            # rng      = rng <- list(
            #   qs   = 150,                # [l/PE/d]
            #   CODs = c(pdf = "nor", mu = 4.378, sigma = 0.751),    # log[g/PE/d]
            #   NH4s = c(pdf = "nor", mu = 1.473, sigma = 0.410),    # log[g/PE/d]
            #   qf   = 0.04,               # [l/s/ha]
            #   CODf = 0,                  # [g/PE/d]
            #   NH4f = 0,                  # [g/PE/d]
            #   CODr = 71,                 # log[mg/l]
            #   NH4r = 1,                  # [mg/l]
            #   nameCSO = "E1",            # [-]
            #   id      = 1,               # [-]
            #   ns      = "FBH Goesdorf",  # [-]
            #   nm      = "Goesdorf",      # [-]
            #   nc      = "Obersauer",     # [-]
            #   numc    = 1,               # [-]
            #   use     = "R/I",           # [-]
            #   Atotal  = 36,              # [ha]
            #   Aimp    = x[1], # c(pdf = "dis", min = 10, max = 25.2, nsample = 6),    # 25.2 [ha]
            #   Cimp    = c(pdf = "tnor", mu = 0.85, sigma = 0.01, lower = 0.65, upper = 0.95),  # 0.85 [-]
            #   Cper    = c(pdf = "tnor", mu = 0.60, sigma = 0.01, lower = 0.30, upper = 0.70),  # 0.60 [-]
            #   tfS     = 1,               # [time steps]
            #   pe      = 650,             # [PE]
            #   Qd      = 5,               # [l/s]
            #   Dd      = 0.150,           # orifice diameter [m]
            #   Cd      = 0.18,            # orifice discharge coefficient [-]
            #   V       = x[2]; # c(pdf = "dis", min = 100, max = 600, nsample = 6),    # 190 [m3]
            #   lev.ini = 1.8,             # initial water level in the chamber [m]
            #   lev2vol = lev2vol          # [m] - [m3]
            # ),
            # ar.model  = ar.model <- list(
            #   CODs    = 0.5,
            #   NH4s    = 0.5,
            #   CODr    = 0.7),
            # var.model = var.model <- list(
            #   inp     = c("", ""), # c("CODs", "NH4s"), # c("", ""),
            #   w       = c(0.04778205, 0.02079010),
            #   A       = matrix(c(9.916452e-01, -8.755558e-05,
            #                      -0.003189094, 0.994553910), nrow=2, ncol=2),
            #   C       = matrix(c(0.009126591, 0.002237936,
            #                      0.002237936, 0.001850941), nrow=2, ncol=2)),
            # folderOutput = paste(workingFolder,"/output",sep=""))
  
            id           <- slot(x, "id")
            nsim         <- slot(x, "nsim")
            seed         <- slot(x, "seed")
            mcCores      <- slot(x, "mcCores")
            ts.input     <- slot(x, "ts.input")
            rng          <- slot(x, "rng")
            ar.model     <- slot(x, "ar.model")
            var.model    <- slot(x, "var.model")
            # folderOutput <- slot(x, "folderOutput")
            
            #=============================================================================================
            # MC simulations and set-up
            #=============================================================================================
            par <- rng
            
            # generation of random numbers (uniform distribution)
            set.seed(seed)
            # i <-1
            npar <- 0
            
            # par[[18]]["pdf"]=="uni"
            # par[[2]][1]
            
            
            indexVAR <- 0
            l <- 1
            # k <- 1
            ifelse(length(var.model[[1]]) > 0,
                   {
                     if(var.model[["inp"]][1] != "") 
                     {for(k in 1:length(var.model[["inp"]])){
                       indexVAR[l] <- which(names(rng)==var.model[["inp"]][k])
                       l <- l+1
                     }
                       indexVAR
                     }else{indexVAR <- c(0,0)}
                   }, indexVAR <- c(0,0)
            )
            
            # ifelse(var.model[["inp"]][1] != "", 
            #        {for(k in 1:length(var.model[["inp"]])){
            #          indexVAR[l] <- which(names(rng)==var.model[["inp"]][k])
            #          l <- l+1}
            #        }, indexVAR <- c(0,0))
            
            # initial nsim.c6 and id.c6 (only case 6)
            nsim.c6 <- 0
            id.c6   <- 0
            
            for (i in 1:length(rng)){
              if(i == indexVAR[2] | is(rng[[i]][1], "list")){next}
              
              if(length(rng[[i]]) > 1 & rng[[i]][1] == "uni") {case <- 1}
              
              if(length(rng[[i]]) > 1 & rng[[i]][1] == "nor" & any(names(rng)[[i]] == names(ar.model)) & 
                 !any(names(rng)[[i]] == var.model[["inp"]])) {case <- 2}

              if(length(rng[[i]]) > 1 & rng[[i]][1] == "nor" & !any(names(rng)[[i]] == names(ar.model))) {case <- 3}
              
              if(length(rng[[i]]) == 1) {case <- 4}
              
              if(length(rng[[i]]) > 1 & rng[[i]][1] == "nor" & any(names(rng)[[i]] == names(ar.model)) & 
                 any(names(rng)[[i]] == var.model[["inp"]])) {case <- 5}
              
              # discrete sampling
              if(length(rng[[i]]) > 1 & rng[[i]][1] == "dis") {case <- 6}
              
              # truncated normal distribution
              # i <- 18 
              if(length(rng[[i]]) > 1 & rng[[i]][1] == "tnor") {case <- 7}
              
                
              # sampling cases definition
              switch(case,
                     # case 1: uniform sampling ---------------------------------------------------------------------------------------------------
                     {r <- matrix(NA, nrow=nsim, ncol = nrow(ts.input)) 
                     j <- 1
                     for(j in 1:nsim){
                       r[j,] <- runif(nrow(ts.input), min=as.numeric(rng[[i]]["min"]), max=as.numeric(rng[[i]]["max"]))    
                     }
                     par[[i]] <- r 
                     npar <- npar+1
                     }, 
                     
                     # case 2: normal autocorrelated time series (AR1 model) ----------------------------------------------------------------------
                     {requireNamespace("lmom")                                                                                      
                       r <- matrix(NA, nrow=nsim, ncol = nrow(ts.input))
                       j <- 1
                       for(j in 1:nsim){
                         # r[j,] <- quanor(runif(nrow(ts.input)), c(as.numeric(rng[[i]]["mu"]), as.numeric(rng[[i]]["sigma"])))
                         
                         y1 <- arima.sim(n = nrow(ts.input) , list(order=c(1,0,0), ar=ar.model[[names(rng)[[i]]]]))  # ar =.5
                         m1 <- mean(y1);    s1 <- sd(y1)
                         m2 <- as.numeric(rng[[i]]["mu"]); s2 <- as.numeric(rng[[i]]["sigma"])
                         y2 <- m2 +(y1-m1)*s2/s1
                         
                         r[j,] <- y2
                       }
                       par[[i]] <- r  
                       npar <- npar+1 
                     },
                     
                     # case 3: normal non-autocorrelated time series ----------------------------------------------------------------------------------
                     {requireNamespace("lmom") 
                       r <- matrix(NA, nrow=nsim, ncol = 1)
                       # i <- 7
                       # j <- 1
                       for(j in 1:nsim){
                         r[j,1] <- quanor(runif(1), c(as.numeric(rng[[i]]["mu"]), as.numeric(rng[[i]]["sigma"])))
                       }
                       par[[i]] <- r  
                       npar <- npar+1 
                     },
                     
                     # case 4: constant value ---------------------------------------------------------------------------------------------------------
                     {par[[i]] <- rng[[i]]
                     },
                     
                     # # case 5: normal auto- and cross-correlated time series (var.model) --------------------------------------------------------------
                     # {
                     #   library(lmom)
                     #   library(mAr)
                     #   r <- matrix(NA, nrow=nsim, ncol = nrow(ts.input))
                     #   
                     #   index <- 0
                     #   l <- 1
                     #   for(k in 1:length(var.model[["inp"]])){
                     #     
                     #     for(j in 1:nsim){
                     #       
                     #       y1 <- mAr.sim(w = var.model[["w"]], A = var.model[["A"]], C = var.model[["C"]], N = nrow(ts.input))
                     #       
                     #       index[l] <- which(names(rng)==var.model[["inp"]][k])
                     #       
                     #       r[j,] <- quanor(runif(nrow(ts.input)), c(as.numeric(rng[[index[l]]]["mu"]), as.numeric(rng[[index[l]]]["sigma"])))
                     #       
                     #       m1 <- mean(y1[,l]); s1 <- sd(y1[,l])
                     #       m2 <- mean(r[j,]); s2 <- sd(r[j,])
                     #       y2 <- m2 +(y1[,l]-m1)*s2/s1
                     #       
                     #       r[j,] <- y2
                     #     }
                     #     l <- l+1
                     #     
                     #     par[[index[l-1]]] <- r
                     #     npar <- npar+1
                     #   }
                     # }
                     
                     # case 5: normal auto- and cross-correlated time series (var.model), parallel code -------------------------------
                     {
                       var1 <- (matrix(NA, nrow=1, ncol = nrow(ts.input)) )
                       var2 <- (matrix(NA, nrow=1, ncol = nrow(ts.input)) )
                       
                       obj1 <- 1:nsim
                       
                       requireNamespace("parallel")
                       requireNamespace("doParallel")
                       requireNamespace("foreach")
                       cl <- makeCluster(mcCores, outfile="")
                       registerDoParallel(cl, cores=mcCores)
                       numCores <- detectCores()
                       
                       rp <- foreach(obj1 = obj1, .packages = c("lmom", "mAr", "parallel", "doParallel", "foreach"), .export=c("obj1"),
                                     .errorhandling = "pass", .verbose=TRUE, .combine = "rbind") %dopar% {
                                       
                                       # j <- 1
                                       j <- obj1
                                       print(paste(j, ", ", mcCores, " of ", numCores, " cores for creating VAR model", sep=""))
                                       
                                       y1 <- mAr.sim(w = var.model[["w"]], A = var.model[["A"]], C = var.model[["C"]], N = nrow(ts.input))
                                       
                                       # code for only bivariate case
                                       m1 <- mean(y1[,1]); s1 <- sd(y1[,1])
                                       m2 <- as.numeric(rng[[indexVAR[1]]]["mu"]); s2 <- as.numeric(rng[[indexVAR[1]]]["sigma"])
                                       y2 <- m2 +(y1[,1]-m1)*s2/s1
                                       var1[1,] <- y2
                                       
                                       m1 <- mean(y1[,2]); s1 <- sd(y1[,2])
                                       m2 <- as.numeric(rng[[indexVAR[2]]]["mu"]); s2 <- as.numeric(rng[[indexVAR[2]]]["sigma"])
                                       y2 <- t(as.data.frame(m2 +(y1[,2]-m1)*s2/s1))
                                       var2[1,] <- y2
                                       
                                       return(c(var1, var2))
                                     }
                       
                       stopCluster(cl)
                       # closeAllConnections()
                       
                       # asignation only for bivariate case:
                       par[[indexVAR[1]]] <- rp[,1:nrow(ts.input)]  
                       par[[indexVAR[2]]] <- rp[,(nrow(ts.input)+1):(nrow(ts.input)*2)]  
                       npar <- npar+2  
                     },
                    
                     # case 6: discrete sampling ------------------------------------------------------------------
                     {id.c6    <- c(id.c6, i)
                     
                     r <- seq(from = as.numeric(rng[[i]]["min"]), to = as.numeric(rng[[i]]["max"]), 
                              length.out=as.numeric(rng[[i]]["nsample"]))   
                     
                     nsim.c6  <- c(nsim.c6, length(r))
                     
                     par[[i]] <- r 
                     npar <- npar+1
                     },
                     
                     # case 7: truncated normal distribution ----------------------------------------------------------------------------------
                     {requireNamespace("msm") 
                       r <- matrix(NA, nrow=nsim, ncol = 1)
                       # i <- 18
                       # j <- 1
                       for(j in 1:nsim){
                         r[j,1] <- rtnorm(n=1, mean=c(as.numeric(rng[[i]]["mu"])),
                                          sd=as.numeric(rng[[i]]["sigma"]),
                                          lower=c(as.numeric(rng[[i]]["lower"])),
                                          upper=c(as.numeric(rng[[i]]["upper"])))
                       }
                       par[[i]] <- r  
                       npar <- npar+1 
                     } 
              )
            }
            
            npars <- length(rng)
            
            return(list(id = id, nsim=nsim, seed=seed, mcCores = mcCores, ts.input=ts.input, rng=rng, par=par, 
                        ar.model=ar.model, var.model=var.model, nsim.c6 = nsim.c6, 
                        id.c6 = id.c6))
          }
)
                         