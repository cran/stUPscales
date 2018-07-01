# Monte-Carlo simulation with pdfs definition or discrete sampling
# author: J.A. Torres-Matallana
# organization: Luxembourg Institute of Science and Technology (LIST), Luxembourg
#               Wagenigen University and Research Centre (WUR), Wageningen, The Netherlands   
# date: 13.04.2017 - 04.09.2017

setGeneric("MC.sim", function(x, EmiStatR.cores) standardGeneric("MC.sim"))

setMethod("MC.sim", signature = c("list", "numeric"), 
          
          function(x, EmiStatR.cores = 0){
            start    <- Sys.time()
            MC_setup = x
            
            #=============================================================================================
            # setting-up MC simulation 
            #=============================================================================================
            
            ## parallel loop for EmiStatR
            
            lap <- matrix(data = NA, nrow = 4, ncol = 2)
            
            cl <- makeCluster(MC_setup[["mcCores"]], outfile="")
            registerDoParallel(cl, cores=MC_setup[["mcCores"]])
            numCores <- detectCores()
            
            # registerDoParallel(cl, cores=i)
            
            # # create progress bar
            # library("tcltk2")
            # pb <- tkProgressBar(title = "progress bar", min = 0, max = nsim, width = 1000)
            nsim    <- MC_setup$nsim
            nsim.c6 <- MC_setup$nsim.c6
            
            if(length(nsim.c6) == 1) {case <- 1} # non-discrete sampling, Monte Carlo simulation
            if(length(nsim.c6) > 1)  {case <- 2} # discrete sampling
            
            switch(case,
                   # case 1: non-discrete sampling, Monte Carlo simulation
                   {obj1 <- 1:nsim
                   a <- system.time(
                     sim1 <- foreach(obj1 = obj1, .packages = c("EmiStatR", "stUPscales"), .export=c("obj1"), 
                                     .errorhandling = "pass", .verbose=TRUE, .combine = "rbind") %dopar% {
                                       # i <- 1
                                       i <- obj1
                                       print(paste(i, ", ", MC_setup[["mcCores"]], " cores", sep=""))
                                       #                     # update progress bar
                                       #                     setTkProgressBar(pb, j, label=paste( round(j/nsim*100, 0), "% done"))
                                       
                                       #for(i in 1:nsim){
                                       
                                       #-------------------------------------------------------------------------------------------------------
                                       # synthetic precipitation definition
                                       #-------------------------------------------------------------------------------------------------------
                                       P1 <- MC_setup$ts.input
                                       if(ncol(P1) > 3){
                                         P1 <- cbind.data.frame(P1[,1], P1[,i+1])
                                       }
                                       
                                       #-------------------------------------------------------------------------------------------------------
                                       # run model (to be modified according to model set-up)
                                       #-------------------------------------------------------------------------------------------------------
                                       #source("/home/atorres/Documents/02_working/06_LIST-Wageningen_PhD/18_models/01_EmiStat-R/R/genEmiStat-R_valida.R", echo=TRUE)
                                       #source(paste(folder,"genEmiStat-R_valida.R", sep=""), echo=TRUE)
                                       
                                       ## defining Input objet 
                                       
                                       # defining structure
                                       # i <- 1
                                       ##***********************************************************************************
                                       ## to be updated for general cases (levesl in not pars, exp CODs and NH4s)
                                       ##***********************************************************************************
                                       # csot <- list(id = setup[["par"]]$id, ns = levels(parsT[["ns"]][i]), nm = levels(parsT[["nm"]][i]),
                                       #              nc = levels(parsT[["nc"]][i]), numc = levels(parsT[["numc"]][i]),
                                       #              use = levels(parsT[["use"]][i]), Atotal = as.numeric(levels(parsT[["Atotal"]][i])),
                                       #              Aimp = as.numeric(parsT[["Aimp"]][i]), tfS = as.numeric(levels(parsT[["tfS"]][i])),
                                       #              pe = as.numeric(levels(parsT[["pe"]][i])), Qd = as.numeric(levels(parsT[["Qd"]][i])),
                                       #              V = as.numeric(levels(parsT[["V"]][i])))
                                       
                                       ifelse(length(MC_setup[["par"]]$Aimp) > 1, Aimpi <- sample(MC_setup[["par"]]$Aimp[i,],1), Aimpi <- MC_setup[["par"]]$Aimp) 
                                       ifelse(length(MC_setup[["par"]]$Atotal) > 1, Atotali <- sample(MC_setup[["par"]]$Atotal[i,],1), Atotali <- MC_setup[["par"]]$Atotal) 
                                       ifelse(length(MC_setup[["par"]]$Cimp) > 1, Cimpi <- sample(MC_setup[["par"]]$Cimp[i,],1), Cimpi <- MC_setup[["par"]]$Cimp)
                                       ifelse(length(MC_setup[["par"]]$Cper) > 1, Cperi <- sample(MC_setup[["par"]]$Cper[i,],1), Cperi <- MC_setup[["par"]]$Cper)
                                       ifelse(length(MC_setup[["par"]]$tfS)  > 1, tfSi  <- round(sample(MC_setup[["par"]]$tfS[i,],1)), tfSi <- round(MC_setup[["par"]]$tfS))
                                       ifelse(length(MC_setup[["par"]]$pe)   > 1, pei   <- sample(MC_setup[["par"]]$pe[i,],1), pei <- MC_setup[["par"]]$pe)
                                       ifelse(length(MC_setup[["par"]]$Qd)   > 1, Qdi   <- sample(MC_setup[["par"]]$Qd[i,],1), Qdi <- MC_setup[["par"]]$Qd)
                                       ifelse(length(MC_setup[["par"]]$V)    > 1, Vi    <- sample(MC_setup[["par"]]$V[i,],1), Vi <- MC_setup[["par"]]$V)
                                       ifelse(length(MC_setup[["par"]]$lev.ini) > 1, lev.inii   <- sample(MC_setup[["par"]]$lev.ini[i,],1), lev.inii <- MC_setup[["par"]]$lev.ini)
                                       ifelse(length(MC_setup[["par"]]$Dd)   > 1, Ddi   <- sample(MC_setup[["par"]]$Dd[i,],1), Ddi <- MC_setup[["par"]]$Dd)
                                       ifelse(length(MC_setup[["par"]]$Cd)   > 1, Cdi   <- sample(MC_setup[["par"]]$Cd[i,],1), Cdi <- MC_setup[["par"]]$Cd)
                                       
                                       if(length(MC_setup[["par"]][["lev2vol"]]) == 2){
                                         # updating Vmax
                                         lev2vol <- MC_setup[["par"]]["lev2vol"]
                                         lev2vol$lev2vol$vol[length(lev2vol$vol)] <- Vi
                                         MC_setup[["par"]]["lev2vol"] <- lev2vol
                                         
                                       }else if(length(MC_setup[["par"]][["lev2vol"]]) > 2){
                                         # # defining lev2vol from available list
                                         lev2vol.rng <- MC_setup[["rng"]]$lev2vol
                                         lev2vol.id <- lapply(lev2vol.rng, function(x){which(x$vol[length(x$vol)] >= Vi)})
                                         lev2vol.id <- which(sapply(lev2vol.id, function(x){match(1, x)}) == 1)
                                         lev2vol.id <- lev2vol.id[1]
                                         
                                         MC_setup[["par"]]$lev2vol <- lev2vol.rng[[lev2vol.id]]
                                         
                                         print(paste0("Aimp = ", Aimpi))
                                         print(paste0("Vi = ", Vi, " - lev2vol = ", MC_setup[["par"]]$lev2vol))
                                       }
                                       
                                       
                                       # ## Estructure 4 - Goesdorf (deterministic setup)
                                       E1 <- list(id = 1, ns = MC_setup[["par"]]$ns, nm = MC_setup[["par"]]$nm,
                                                  nc = MC_setup[["par"]]$nc, numc = MC_setup[["par"]]$numc,
                                                  use = MC_setup[["par"]]$use, Atotal = Atotali,
                                                  Aimp = Aimpi, 
                                                  tfS = tfSi,
                                                  pe = pei, Qd = Qdi,
                                                  V = Vi, lev2vol = MC_setup[["par"]]$lev2vol, 
                                                  lev.ini = lev.inii,
                                                  Dd = Ddi,
                                                  Cd = Cdi,
                                                  Cimp = Cimpi, Cper = Cperi)
                                       
                                       ## defining input objet 
                                       ifelse(length(MC_setup[["par"]][["qs"]])   > 1, qsi   <- sample(MC_setup[["par"]][["qs"]][i,],1),   qsi   <- MC_setup[["par"]][["qs"]][1]) 
                                       ifelse(length(MC_setup[["par"]][["CODs"]]) > 1, CODsi <- exp(sample(MC_setup[["par"]][["CODs"]][i,],1)), CODsi <- MC_setup[["par"]][["CODs"]][1])
                                       ifelse(length(MC_setup[["par"]][["NH4s"]]) > 1, NH4si <- exp(sample(MC_setup[["par"]][["NH4s"]][i,],1)), NH4si <- MC_setup[["par"]][["NH4s"]][1])
                                       ifelse(length(MC_setup[["par"]][["qf"]])   > 1, qfi   <- sample(MC_setup[["par"]][["qf"]][i,],1),   qfi   <- MC_setup[["par"]][["qf"]][1])
                                       ifelse(length(MC_setup[["par"]][["CODf"]]) > 1, CODfi <- sample(MC_setup[["par"]][["CODf"]][i,],1), CODfi <- MC_setup[["par"]][["CODf"]][1])
                                       ifelse(length(MC_setup[["par"]][["NH4f"]]) > 1, NH4fi <- sample(MC_setup[["par"]][["NH4f"]][i,],1), NH4fi <- MC_setup[["par"]][["NH4f"]][1])
                                       ifelse(length(MC_setup[["par"]][["CODr"]]) > 1, CODri <- exp(sample(MC_setup[["par"]][["CODr"]][i,],1)), CODri <- MC_setup[["par"]][["CODr"]][1])
                                       ifelse(length(MC_setup[["par"]][["NH4r"]]) > 1, NH4ri <- sample(MC_setup[["par"]][["NH4r"]][i,],1), NH4ri <- MC_setup[["par"]][["NH4r"]][1])
                                       
                                       input.user <- input(spatial = 0, zero = 1e-5, folder = system.file("shiny", package = "EmiStatR"),
                                                           cores = EmiStatR.cores,
                                                           ww = list(qs = qsi, CODs = CODsi, NH4s = NH4si), 
                                                           inf = list(qf= qfi, CODf = CODfi, NH4f = NH4fi),
                                                           rw = list(CODr = CODri, NH4r = NH4ri, stat = "Goesdorf"), P1 = P1, 
                                                           st = list(E1=E1), 
                                                           export = 0)
                                       
                                       ##***********************************************************************************
                                       ## end
                                       ##***********************************************************************************
                                       
                                       ## executing simulation 
                                       #library(EmiStatR)  
                                       sim <- EmiStatR(input.user)
                                       
                                       #-------------------------------------------------------------------------------------------------------
                                       # Monte-Carlo output time-series matrix  (to be modified according to model output)
                                       #-------------------------------------------------------------------------------------------------------
                                       data1       <- sim[[1]]$out1
                                       data1       <- data1[,2:ncol(data1)]
                                       # head(data1)
                                       mcVChamber  <- data1[,8]      # V_Chamber (V_Tank)
                                       mcVsv       <- data1[,9]      # V_Spill-volume (V_Ov)
                                       # mcQsv     <- data1[,16]     # Q_Sv (Q_Ov)
                                       mcBCODsv    <- data1[,10]     # B_COD_Sv (B_COD_Ov)
                                       mcBNH4sv    <- data1[,11]     # B_NH4_Sv (B_NH4_Ov)
                                       mcCCODsv    <- data1[,12]     # C_COD_Sv (C_COD_Ov)
                                       mcCNH4sv    <- data1[,13]     # C_NH4_Sv (C_NH4_Ov)
                                       
                                       mc1 <- list(mcVChamber=mcVChamber, mcVsv=mcVsv, #mcQov=mcQov, 
                                                   mcBCODsv=mcBCODsv, mcBNH4sv=mcBNH4sv, mcCCODsv=mcCCODsv, mcCNH4sv=mcCNH4sv,
                                                   time = P1[,1], st=E1, inputUser=input.user)
                                       
                                       return(mc1)
                                     } # end of foreach loop for Monte-Carlo simulations
                   )  # end system.time
                   }, # end case 1
                   
                   # case 2: discrete sampling, orifice
                   {
                     id.c6 <- MC_setup$id.c6
                     discrete.samples <- matrix(NA, nrow = 1, ncol = (length(id.c6)-1))
                     
                     if(length(id.c6) == 3){
                       ie <- 1
                       id <- 1
                       # ia <- 1
                       # ib <- 1
                       # ic <- 1
                       # id <- 1
                       for(ia in 1:length(MC_setup$par[id.c6[2]][[1]])){
                         for(ib in 1:length(MC_setup$par[id.c6[3]][[1]])){
                           if(ie == 1){ 
                             discrete.samples[id,] <- c(MC_setup$par[id.c6[2]][[1]][1], MC_setup$par[id.c6[3]][[1]][1])
                             ie <- ie+1
                           }else{discrete.samples <- rbind(discrete.samples, c(MC_setup$par[id.c6[2]][[1]][ia], 
                                                                               MC_setup$par[id.c6[3]][[1]][ib]))}
                         }
                       }
                       
                       colnames(discrete.samples) <- names(MC_setup$par[id.c6[2:3]])
                     } # end if(length(id.c6) == 3)
                     
                     if(length(id.c6) == 4){
                       ie <- 1
                       # ia <- 1
                       # ib <- 1
                       # ic <- 1
                       # id <- 1
                       for(ia in 1:length(MC_setup$par[id.c6[2]][[1]])){
                         for(ib in 1:length(MC_setup$par[id.c6[3]][[1]])){
                           for(ic in 1:length(MC_setup$par[id.c6[4]][[1]])){
                             # for(id in 1:length(MC_setup$par[id.c6[5]][[1]])){
                             if(ie == 1){ 
                               discrete.samples[id,] <- c(MC_setup$par[id.c6[2]][[1]][1], MC_setup$par[id.c6[3]][[1]][1], 
                                                          MC_setup$par[id.c6[4]][[1]][1])
                               ie <- ie+1
                             }else{discrete.samples <- rbind(discrete.samples, c(MC_setup$par[id.c6[2]][[1]][ia], 
                                                                                 MC_setup$par[id.c6[3]][[1]][ib], 
                                                                                 MC_setup$par[id.c6[4]][[1]][ic]))}
                             # }
                           }
                         }
                       }
                       
                       colnames(discrete.samples) <- names(MC_setup$par[id.c6[2:4]])
                     } # end if(length(id.c6) == 4)  
                     
                     if(length(id.c6) == 5){
                       ie <- 1
                       # ia <- 1
                       # ib <- 1
                       # ic <- 1
                       # id <- 1
                       for(ia in 1:length(MC_setup$par[id.c6[2]][[1]])){
                         for(ib in 1:length(MC_setup$par[id.c6[3]][[1]])){
                           for(ic in 1:length(MC_setup$par[id.c6[4]][[1]])){
                             for(id in 1:length(MC_setup$par[id.c6[5]][[1]])){
                               if(ie == 1){ 
                                 discrete.samples[id,] <- c(MC_setup$par[id.c6[2]][[1]][1], 
                                                            MC_setup$par[id.c6[3]][[1]][1], 
                                                            MC_setup$par[id.c6[4]][[1]][1],
                                                            MC_setup$par[id.c6[5]][[1]][1])
                                 ie <- ie+1
                               }else{discrete.samples <- rbind(discrete.samples, c(MC_setup$par[id.c6[2]][[1]][ia], 
                                                                                   MC_setup$par[id.c6[3]][[1]][ib], 
                                                                                   MC_setup$par[id.c6[4]][[1]][ic],
                                                                                   MC_setup$par[id.c6[5]][[1]][ic]))}
                             }
                           }
                         }
                       }
                       
                       colnames(discrete.samples) <- names(MC_setup$par[id.c6[2:5]])
                     } # end if(length(id.c6) == 5) 
                     
                     if(length(id.c6) == 7){
                       ih <- 1
                       # ia <- 1
                       # ib <- 1
                       # ic <- 1
                       # id <- 1
                       for(ia in 1:length(MC_setup$par[id.c6[2]][[1]])){
                         for(ib in 1:length(MC_setup$par[id.c6[3]][[1]])){
                           for(ic in 1:length(MC_setup$par[id.c6[4]][[1]])){
                             for(id in 1:length(MC_setup$par[id.c6[5]][[1]])){
                               for(ie in 1:length(MC_setup$par[id.c6[6]][[1]])){
                                 for(iff in 1:length(MC_setup$par[id.c6[7]][[1]])){
                                   # for(ig in 1:length(MC_setup$par[id.c6[8]][[1]])){
                                   if(ih == 1){ 
                                     discrete.samples[1,] <- c(MC_setup$par[id.c6[2]][[1]][1], MC_setup$par[id.c6[3]][[1]][1], 
                                                               MC_setup$par[id.c6[4]][[1]][1], MC_setup$par[id.c6[5]][[1]][1],
                                                               MC_setup$par[id.c6[6]][[1]][1], MC_setup$par[id.c6[7]][[1]][1])#,
                                     # MC_setup$par[id.c6[8]][[1]][1])
                                     ih <- ih+1
                                   }else{discrete.samples <- rbind(discrete.samples, c(MC_setup$par[id.c6[2]][[1]][ia], 
                                                                                       MC_setup$par[id.c6[3]][[1]][ib], 
                                                                                       MC_setup$par[id.c6[4]][[1]][ic],
                                                                                       MC_setup$par[id.c6[5]][[1]][id],
                                                                                       MC_setup$par[id.c6[6]][[1]][ie],
                                                                                       MC_setup$par[id.c6[7]][[1]][iff])#,
                                                                   # MC_setup$par[id.c6[8]][[1]][ig])
                                   )}
                                 }#}
                               }}}}}
                       
                       colnames(discrete.samples) <- names(MC_setup$par[id.c6[2:7]])
                     } # end if(length(id.c6) == 7)
                     
                     nsamples <- nrow(discrete.samples)
                     
                     if(ncol(MC_setup[["ts.input"]]) == 2){
                       obj1 <- 1:nsamples
                       a <- system.time(
                         sim1 <- foreach(obj1 = obj1, .packages = c("EmiStatR", "stUPscales"), .export=c("obj1"), # previously "AggProp" instead of "stUPscales"
                                         .errorhandling = "pass", .verbose=TRUE, .combine = "rbind") %dopar% {
                                           # i <- 1
                                           i <- obj1
                                           print(paste("sim-discrete ", i, "/", nsamples, ",", MC_setup[["mcCores"]], " cores, running in parallel mode...", sep=""))
                                           #                     # update progress bar
                                           #                     setTkProgressBar(pb, j, label=paste( round(j/nsim*100, 0), "% done"))
                                           
                                           #for(i in 1:nsim){
                                           
                                           #-------------------------------------------------------------------------------------------------------
                                           # synthetic precipitation definition
                                           #-------------------------------------------------------------------------------------------------------
                                           P1 <- MC_setup$ts.input
                                           if(ncol(P1) > 3){
                                             P1 <- cbind.data.frame(P1[,1], P1[,i+1])
                                           }
                                           
                                           #-------------------------------------------------------------------------------------------------------
                                           # run model (to be modified according to model set-up)
                                           #-------------------------------------------------------------------------------------------------------
                                           #source("/home/atorres/Documents/02_working/06_LIST-Wageningen_PhD/18_models/01_EmiStat-R/R/genEmiStat-R_valida.R", echo=TRUE)
                                           #source(paste(folder,"genEmiStat-R_valida.R", sep=""), echo=TRUE)
                                           
                                           ## defining Input objet 
                                           
                                           # defining structure
                                           # i <- 1
                                           ##***********************************************************************************
                                           ## to be updated for general cases (levesl in not pars, exp CODs and NH4s)
                                           ##***********************************************************************************
                                           # csot <- list(id = setup[["par"]]$id, ns = levels(parsT[["ns"]][i]), nm = levels(parsT[["nm"]][i]),
                                           #              nc = levels(parsT[["nc"]][i]), numc = levels(parsT[["numc"]][i]),
                                           #              use = levels(parsT[["use"]][i]), Atotal = as.numeric(levels(parsT[["Atotal"]][i])),
                                           #              Aimp = as.numeric(parsT[["Aimp"]][i]), tfS = as.numeric(levels(parsT[["tfS"]][i])),
                                           #              pe = as.numeric(levels(parsT[["pe"]][i])), Qd = as.numeric(levels(parsT[["Qd"]][i])),
                                           #              V = as.numeric(levels(parsT[["V"]][i])))
                                           
                                           # ## Estructure 4 - Goesdorf (deterministic setup)
                                           ifelse(length(MC_setup[["par"]]$Aimp) > 1, Aimpi <- discrete.samples[i,"Aimp"], Aimpi <- MC_setup[["par"]]$Aimp) 
                                           ifelse(length(MC_setup[["par"]]$Atotal) > 1, Atotali <- discrete.samples[i,"Atotal"], Atotali <- MC_setup[["par"]]$Atotal) 
                                           ifelse(length(MC_setup[["par"]]$Cimp) > 1, Cimpi <- discrete.samples[i,"Cimp"], Cimpi <- MC_setup[["par"]]$Cimp) 
                                           ifelse(length(MC_setup[["par"]]$Cper) > 1, Cperi <- discrete.samples[i,"Cper"], Cperi <- MC_setup[["par"]]$Cper) 
                                           ifelse(length(MC_setup[["par"]]$tfS)  > 1, tfSi  <- round(discrete.samples[i,"tfS"]), tfSi <- round(MC_setup[["par"]]$tfS))
                                           ifelse(length(MC_setup[["par"]]$pe)   > 1, pei   <- discrete.samples[i,"pe"], pei <- MC_setup[["par"]]$pe)
                                           ifelse(length(MC_setup[["par"]]$Qd)   > 1, Qdi   <- discrete.samples[i,"Qd"], Qdi <- MC_setup[["par"]]$Qd)
                                           ifelse(length(MC_setup[["par"]]$V)    > 1, Vi    <- discrete.samples[i, "V"], Vi <- MC_setup[["par"]]$V)
                                           ifelse(length(MC_setup[["par"]]$lev.ini) > 1, lev.inii   <- discrete.samples[i,"lev.ini"], lev.inii <- MC_setup[["par"]]$lev.ini)
                                           ifelse(length(MC_setup[["par"]]$Dd)   > 1, Ddi   <- discrete.samples[i, "Dd"], Ddi <- MC_setup[["par"]]$Dd)
                                           ifelse(length(MC_setup[["par"]]$Cd)   > 1, Cdi   <- discrete.samples[i, "Cd"], Cdi <- MC_setup[["par"]]$Cd)
                                           
                                           # updating Vmax
                                           if(length(MC_setup[["par"]][["lev2vol"]]) == 2){
                                             lev2vol <- MC_setup[["par"]]["lev2vol"]
                                             lev2vol$lev2vol$vol[length(lev2vol$vol)] <- Vi
                                             MC_setup[["par"]]["lev2vol"] <- lev2vol
                                           }
                                           # ## Estructure 4 - Goesdorf (deterministic setup)
                                           E1 <- list(id = 1, ns = MC_setup[["par"]]$ns, nm = MC_setup[["par"]]$nm,
                                                      nc = MC_setup[["par"]]$nc, numc = MC_setup[["par"]]$numc,
                                                      use = MC_setup[["par"]]$use, Atotal = Atotali,
                                                      Aimp = Aimpi, Cimp = Cimpi, Cper = Cperi, 
                                                      tfS = tfSi,
                                                      pe = pei, Qd = Qdi,
                                                      V = Vi,
                                                      lev2vol = MC_setup[["par"]]$lev2vol,
                                                      lev.ini = lev.inii,
                                                      Dd = Ddi,
                                                      Cd = Cdi)
                                           
                                           ## defining input objet 
                                           ifelse(length(MC_setup[["par"]][["qs"]])   > 1, qsi   <- discrete.samples[i,"qs"],   qsi   <- MC_setup[["par"]][["qs"]][1]) 
                                           ifelse(length(MC_setup[["par"]][["CODs"]]) > 1, CODsi <- discrete.samples[i,"CODs"], CODsi <- MC_setup[["par"]][["CODs"]][1])
                                           ifelse(length(MC_setup[["par"]][["NH4s"]]) > 1, NH4si <- discrete.samples[i,"NH4s"], NH4si <- MC_setup[["par"]][["NH4s"]][1])
                                           ifelse(length(MC_setup[["par"]][["qf"]])   > 1, qfi   <- discrete.samples[i,"qf"],   qfi   <- MC_setup[["par"]][["qf"]][1])
                                           ifelse(length(MC_setup[["par"]][["CODf"]]) > 1, CODfi <- discrete.samples[i,"CODf"], CODfi <- MC_setup[["par"]][["CODf"]][1])
                                           ifelse(length(MC_setup[["par"]][["NH4f"]]) > 1, NH4fi <- discrete.samples[i,"NH4f"], NH4fi <- MC_setup[["par"]][["NH4f"]][1])
                                           ifelse(length(MC_setup[["par"]][["CODr"]]) > 1, CODri <- discrete.samples[i,"CODr"], CODri <- MC_setup[["par"]][["CODr"]][1])
                                           ifelse(length(MC_setup[["par"]][["NH4r"]]) > 1, NH4ri <- discrete.samples[i,"NH4r"], NH4ri <- MC_setup[["par"]][["NH4r"]][1])
                                           
                                           input.user <- input(spatial = 0, zero = 1e-5, folder = system.file("shiny", package = "EmiStatR"),
                                                               cores = EmiStatR.cores,
                                                               ww = list(qs = qsi, CODs = CODsi, NH4s = NH4si), 
                                                               inf = list(qf= qfi, CODf = CODfi, NH4f = NH4fi),
                                                               rw = list(CODr = CODri, NH4r = NH4ri, stat = "Goesdorf"), P1 = P1, 
                                                               st = list(E1=E1), 
                                                               export = 0)
                                           
                                           ##***********************************************************************************
                                           ## end
                                           ##***********************************************************************************
                                           
                                           
                                           ## executing simulation 
                                           #library(EmiStatR)  
                                           sim <- EmiStatR(input.user)
                                           
                                           #-------------------------------------------------------------------------------------------------------
                                           # Monte-Carlo output time-series matrix  (to be modified according to model output)
                                           #-------------------------------------------------------------------------------------------------------
                                           data1       <- sim[[1]]$out1
                                           data1       <- data1[,2:ncol(data1)]
                                           # head(data1)
                                           mcVChamber  <- data1[,8]      # V_Chamber (V_Tank)
                                           mcVsv       <- data1[,9]      # V_Spill-volume (V_Ov)
                                           # mcQsv     <- data1[,16]     # Q_Sv (Q_Ov)
                                           mcBCODsv    <- data1[,10]     # B_COD_Sv (B_COD_Ov)
                                           mcBNH4sv    <- data1[,11]     # B_NH4_Sv (B_NH4_Ov)
                                           mcCCODsv    <- data1[,12]     # C_COD_Sv (C_COD_Ov)
                                           mcCNH4sv    <- data1[,13]     # C_NH4_Sv (C_NH4_Ov)
                                           
                                           mc1 <- list(mcVChamber=mcVChamber, mcVsv=mcVsv, #mcQov=mcQov, 
                                                       mcBCODsv=mcBCODsv, mcBNH4sv=mcBNH4sv, mcCCODsv=mcCCODsv, mcCNH4sv=mcCNH4sv,
                                                       time = P1[,1], st=E1, inputUser=input.user)
                                           
                                           return(mc1)
                                         } # end of foreach loop for Monte-Carlo simulations
                       ) # end system.time
                     }else{ # end if(ncol(MC_setup[["ts.input"]]) == 2)
                       
                       obj1 <- 1:MC_setup[["nsim"]]
                       a <- system.time(
                         sim1 <- foreach(obj1 = obj1, .packages = c("EmiStatR", "stUPscales"), .export=c("obj1"), # previously "AggProp" instead of "stUPscales"
                                         .errorhandling = "pass", .verbose=TRUE, .combine = "rbind") %dopar% {
                                           # j <- 1
                                           j <- obj1
                                           print(paste("simulation ", j, "/", MC_setup[["nsim"]], ",", MC_setup[["mcCores"]], " cores, running in parallel mode...", sep=""))
                                           #                     # update progress bar
                                           #                     setTkProgressBar(pb, j, label=paste( round(j/nsim*100, 0), "% done"))
                                           
                                           #for(i in 1:nsim){
                                           for(i in 1:nsamples){
                                             print(paste("sub-simulation ", i, "/", nsamples, sep=""))
                                             #-------------------------------------------------------------------------------------------------------
                                             # synthetic precipitation definition
                                             #-------------------------------------------------------------------------------------------------------
                                             P1 <- MC_setup$ts.input
                                             if(ncol(P1) > 3){
                                               P1 <- cbind.data.frame(P1[,1], P1[,i+1])
                                             }
                                             
                                             #-------------------------------------------------------------------------------------------------------
                                             # run model (to be modified according to model set-up)
                                             #-------------------------------------------------------------------------------------------------------
                                             #source("/home/atorres/Documents/02_working/06_LIST-Wageningen_PhD/18_models/01_EmiStat-R/R/genEmiStat-R_valida.R", echo=TRUE)
                                             #source(paste(folder,"genEmiStat-R_valida.R", sep=""), echo=TRUE)
                                             
                                             ## defining Input objet 
                                             
                                             # defining structure
                                             # i <- 1
                                             ##***********************************************************************************
                                             ## to be updated for general cases (levesl in not pars, exp CODs and NH4s)
                                             ##***********************************************************************************
                                             # csot <- list(id = setup[["par"]]$id, ns = levels(parsT[["ns"]][i]), nm = levels(parsT[["nm"]][i]),
                                             #              nc = levels(parsT[["nc"]][i]), numc = levels(parsT[["numc"]][i]),
                                             #              use = levels(parsT[["use"]][i]), Atotal = as.numeric(levels(parsT[["Atotal"]][i])),
                                             #              Aimp = as.numeric(parsT[["Aimp"]][i]), tfS = as.numeric(levels(parsT[["tfS"]][i])),
                                             #              pe = as.numeric(levels(parsT[["pe"]][i])), Qd = as.numeric(levels(parsT[["Qd"]][i])),
                                             #              V = as.numeric(levels(parsT[["V"]][i])))
                                             
                                             # ## Estructure 4 - Goesdorf (deterministic setup)
                                             ifelse(length(MC_setup[["par"]]$Aimp) > 1, Aimpi <- discrete.samples[i,"Aimp"], Aimpi <- MC_setup[["par"]]$Aimp) 
                                             ifelse(length(MC_setup[["par"]]$Atotal) > 1, Atotali <- discrete.samples[i,"Atotal"], Atotali <- MC_setup[["par"]]$Atotal) 
                                             ifelse(length(MC_setup[["par"]]$Cimp) > 1, Cimpi <- discrete.samples[i,"Cimp"], Cimpi <- MC_setup[["par"]]$Cimp) 
                                             ifelse(length(MC_setup[["par"]]$Cper) > 1, Cperi <- discrete.samples[i,"Cper"], Cperi <- MC_setup[["par"]]$Cper) 
                                             ifelse(length(MC_setup[["par"]]$tfS)  > 1, tfSi  <- round(discrete.samples[i,"tfS"]), tfSi <- round(MC_setup[["par"]]$tfS))
                                             ifelse(length(MC_setup[["par"]]$pe)   > 1, pei   <- discrete.samples[i,"pe"], pei <- MC_setup[["par"]]$pe)
                                             ifelse(length(MC_setup[["par"]]$Qd)   > 1, Qdi   <- discrete.samples[i,"Qd"], Qdi <- MC_setup[["par"]]$Qd)
                                             ifelse(length(MC_setup[["par"]]$V)    > 1, Vi    <- discrete.samples[i, "V"], Vi <- MC_setup[["par"]]$V)
                                             ifelse(length(MC_setup[["par"]]$lev.ini) > 1, lev.inii   <- discrete.samples[i,"lev.ini"], lev.inii <- MC_setup[["par"]]$lev.ini)
                                             ifelse(length(MC_setup[["par"]]$Dd)   > 1, Ddi   <- discrete.samples[i, "Dd"], Ddi <- MC_setup[["par"]]$Dd)
                                             ifelse(length(MC_setup[["par"]]$Cd)   > 1, Cdi   <- discrete.samples[i, "Cd"], Cdi <- MC_setup[["par"]]$Cd)
                                             
                                             ## defining lev2vol from available list
                                             if(length(MC_setup[["par"]][["lev2vol"]]) > 2){
                                               lev2vol.rng <- MC_setup[["rng"]]$lev2vol
                                               lev2vol.id <- lapply(lev2vol.rng, function(x){which(x$vol[length(x$vol)] >= Vi)})
                                               lev2vol.id <- which(sapply(lev2vol.id, function(x){match(1, x)}) == 1)
                                               lev2vol.id <- lev2vol.id[1]
                                               
                                               MC_setup[["par"]]$lev2vol <- lev2vol.rng[[lev2vol.id]]
                                               print(paste0("Vi = ", Vi, " - lev2vol = ", MC_setup[["par"]]$lev2vol))
                                             }
                                             
                                             E1 <- list(id = 1, ns = MC_setup[["par"]]$ns, nm = MC_setup[["par"]]$nm,
                                                        nc = MC_setup[["par"]]$nc, numc = MC_setup[["par"]]$numc,
                                                        use = MC_setup[["par"]]$use, Atotal = Atotali,
                                                        Aimp = Aimpi, Cimp = Cimpi, Cper = Cperi,
                                                        tfS = tfSi,
                                                        pe = pei, Qd = Qdi,
                                                        V = Vi,
                                                        lev2vol = MC_setup[["par"]]$lev2vol,
                                                        lev.ini = lev.inii,
                                                        Dd = Ddi,
                                                        Cd = Cdi)
                                             
                                             ## defining input objet 
                                             ifelse(length(MC_setup[["par"]][["qs"]])   > 1, qsi   <- discrete.samples[i,"qs"],   qsi   <- MC_setup[["par"]][["qs"]][1]) 
                                             ifelse(length(MC_setup[["par"]][["CODs"]]) > 1, CODsi <- discrete.samples[i,"CODs"], CODsi <- MC_setup[["par"]][["CODs"]][1])
                                             ifelse(length(MC_setup[["par"]][["NH4s"]]) > 1, NH4si <- discrete.samples[i,"NH4s"], NH4si <- MC_setup[["par"]][["NH4s"]][1])
                                             ifelse(length(MC_setup[["par"]][["qf"]])   > 1, qfi   <- discrete.samples[i,"qf"],   qfi   <- MC_setup[["par"]][["qf"]][1])
                                             ifelse(length(MC_setup[["par"]][["CODf"]]) > 1, CODfi <- discrete.samples[i,"CODf"], CODfi <- MC_setup[["par"]][["CODf"]][1])
                                             ifelse(length(MC_setup[["par"]][["NH4f"]]) > 1, NH4fi <- discrete.samples[i,"NH4f"], NH4fi <- MC_setup[["par"]][["NH4f"]][1])
                                             ifelse(length(MC_setup[["par"]][["CODr"]]) > 1, CODri <- discrete.samples[i,"CODr"], CODri <- MC_setup[["par"]][["CODr"]][1])
                                             ifelse(length(MC_setup[["par"]][["NH4r"]]) > 1, NH4ri <- discrete.samples[i,"NH4r"], NH4ri <- MC_setup[["par"]][["NH4r"]][1])
                                             
                                             input.user <- input(spatial = 0, zero = 1e-5, folder = system.file("shiny", package = "EmiStatR"),
                                                                 cores = EmiStatR.cores,
                                                                 ww = list(qs = qsi, CODs = CODsi, NH4s = NH4si), 
                                                                 inf = list(qf= qfi, CODf = CODfi, NH4f = NH4fi),
                                                                 rw = list(CODr = CODri, NH4r = NH4ri, stat = "Goesdorf"), P1 = P1, 
                                                                 st = list(E1=E1), 
                                                                 export = 0)
                                             
                                             ##***********************************************************************************
                                             ## end
                                             ##***********************************************************************************
                                             
                                             
                                             ## executing simulation 
                                             #library(EmiStatR)  
                                             sim <- EmiStatR(input.user)
                                             
                                             #-------------------------------------------------------------------------------------------------------
                                             # Monte-Carlo output time-series matrix  (to be modified according to model output)
                                             #-------------------------------------------------------------------------------------------------------
                                             data1       <- sim[[1]]$out1
                                             data1       <- data1[,2:ncol(data1)]
                                             # head(data1)
                                             # mcVChamber  <- data1[,8]      # V_Chamber (V_Tank)
                                             mcVsv       <- data1[,9]      # V_Spill-volume (V_Ov)
                                             # mcQsv     <- data1[,16]     # Q_Sv (Q_Ov)
                                             # mcBCODsv    <- data1[,10]     # B_COD_Sv (B_COD_Ov)
                                             mcBNH4sv    <- data1[,11]     # B_NH4_Sv (B_NH4_Ov)
                                             # mcCCODsv    <- data1[,12]     # C_COD_Sv (C_COD_Ov)
                                             mcCNH4sv    <- data1[,13]     # C_NH4_Sv (C_NH4_Ov)
                                             
                                             # assign(paste0("sample", i), list(mcVChamber=mcVChamber, mcVsv=mcVsv, #mcQov=mcQov, 
                                             #                mcBCODsv=mcBCODsv, mcBNH4sv=mcBNH4sv, mcCCODsv=mcCCODsv, mcCNH4sv=mcCNH4sv,
                                             #                time = P1[,1], inputUser=input.user))
                                             
                                             assign(paste0("sample", i), list(mcVsv=mcVsv, 
                                                                              mcBNH4sv=mcBNH4sv,mcCNH4sv=mcCNH4sv,
                                                                              time = P1[,1], inputUser=input.user))
                                             
                                           } # end for(i in 1:nsamples)
                                           samples <- sapply(ls(pattern = "sample"), function(x) list(get(x)))
                                           # return(list(samples = samples, discrete.samples = discrete.samples))
                                           return(samples)
                                         } # end of foreach loop for Monte-Carlo simulations
                       ) # end system.time
                     } # end else from if(ncol(MC_setup[["ts.input"]]) == 2)
                   } # end case 2
            )
            
            stopCluster(cl)
            # closeAllConnections()
            
            ## end timing
            lap[1,1] <- nsim
            lap[2,1] <- a[3]/60
            lap[3,1] <- EmiStatR.cores
            lap[4,1] <- MC_setup[["mcCores"]]
            lap[1,2] <- "sims"
            lap[2,2] <- "min"
            lap[3,2] <- "cores"
            lap[4,2] <- "MC cores"
            
            end     <- Sys.time()
            elapsed <- difftime(end, start, units="hours"); elapsed
            timing  <- list(start=start, end=end, elapsed_hours=elapsed)
            
            mc <- list(MC_setup = MC_setup)
            mc <- c(mc, list(timing = timing))
            mc <- c(mc, list(lap = lap))
            
            #str(mc)
            #mc["par"]
            #mc["mcVTank"]
            
            #-------------------------------------------------------------------------------------------------------
            # saving outputs
            #-------------------------------------------------------------------------------------------------------
            currentDir <- getwd()
            # setwd(MC_setup$folderOutput)
            save(sim1, file=paste("sim1.RData", sep=""))
            save(mc, file=paste("mc.RData", sep=""))
            setwd(currentDir)
            
            ## closing progress bar
            # close(pb)
            print(paste("End of", nsim[length(nsim)], "Monte-Carlo simulations.",sep=" "))
            
            return(list(mc=mc, sim1=sim1))
          }
)
