# Monte-Carlo calibration: water quantity (based on MC.sim.R)
# author: J.A. Torres-Matallana
# organization: Luxembourg Institute of Science and Technology (LIST), Luxembourg
#               Wagenigen University and Research Centre (WUR), Wageningen, The Netherlands   
# date: 13.04.2017 - 13.04.2017


setGeneric("MC.calibra", function(x, obs, EmiStatR.cores) standardGeneric("MC.calibra"))

setMethod("MC.calibra", signature = c("list", "inputObs", "numeric"), 
          
          function(x, obs, EmiStatR.cores = 0){
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
            if(length(nsim.c6) > 1)  {case <- 2} # discrete sampling, orifice
            
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
                                       
                                       ifelse(length(MC_setup[["par"]]$Aimp)    > 1, Aimpi    <- sample(MC_setup[["par"]]$Aimp[i,],1), Aimpi <- MC_setup[["par"]]$Aimp) 
                                       ifelse(length(MC_setup[["par"]]$Atotal)  > 1, Atotali  <- sample(MC_setup[["par"]]$Atotal[i,],1), Atotali <- MC_setup[["par"]]$Atotal) 
                                       ifelse(length(MC_setup[["par"]]$tfS)     > 1, tfSi     <- round(sample(MC_setup[["par"]]$tfS[i,],1)), tfSi <- round(MC_setup[["par"]]$tfS))
                                       ifelse(length(MC_setup[["par"]]$pe)      > 1, pei      <- sample(MC_setup[["par"]]$pe[i,],1), pei <- MC_setup[["par"]]$pe)
                                       ifelse(length(MC_setup[["par"]]$Qd)      > 1, Qdi      <- sample(MC_setup[["par"]]$Qd[i,],1), Qdi <- MC_setup[["par"]]$Qd)
                                       ifelse(length(MC_setup[["par"]]$V)       > 1, Vi       <- sample(MC_setup[["par"]]$V[i,],1), Vi <- MC_setup[["par"]]$V)
                                       ifelse(length(MC_setup[["par"]]$lev.ini) > 1, lev.inii <- sample(MC_setup[["par"]]$lev.ini[i,],1), lev.inii <- MC_setup[["par"]]$lev.ini)
                                       ifelse(length(MC_setup[["par"]]$Dd)      > 1, Ddi      <- sample(MC_setup[["par"]]$Dd[i,],1), Ddi <- MC_setup[["par"]]$Dd)
                                       ifelse(length(MC_setup[["par"]]$Cd)      > 1, Cdi      <- sample(MC_setup[["par"]]$Cd[i,],1), Cdi <- MC_setup[["par"]]$Cd)
                                       ifelse(length(MC_setup[["par"]]$Cimp)    > 1, Cimpi    <- sample(MC_setup[["par"]]$Cimp[i,],1), Cimpi <- MC_setup[["par"]]$Cimp)
                                       ifelse(length(MC_setup[["par"]]$Cper)    > 1, Cperi    <- sample(MC_setup[["par"]]$Cper[i,],1), Cperi <- MC_setup[["par"]]$Cper)
                                       
                                       # ## Estructure 4 - Goesdorf (deterministic setup)
                                       E1 <- list(id = 1, ns = MC_setup[["par"]]$ns, nm = MC_setup[["par"]]$nm,
                                                  nc = MC_setup[["par"]]$nc, numc = MC_setup[["par"]]$numc,
                                                  use = MC_setup[["par"]]$use, Atotal = Atotali,
                                                  Aimp = Aimpi, tfS = tfSi,
                                                  pe = pei, Qd = Qdi,
                                                  V = Vi, lev2vol = MC_setup[["par"]]$lev2vol, 
                                                  lev.ini = lev.inii,
                                                  Dd = Ddi,
                                                  Cd = Cdi,
                                                  Cimp = Cimpi, Cper = Cperi)
                                       
                                       ## defining input objet 
                                       ifelse(length(MC_setup[["par"]][["qs"]])   > 1, qsi   <- sample(MC_setup[["par"]][["qs"]][i,],1),   qsi   <- MC_setup[["par"]][["qs"]][1]) 
                                       ifelse(length(MC_setup[["par"]][["CODs"]]) > 1, CODsi <- sample(MC_setup[["par"]][["CODs"]][i,],1), CODsi <- MC_setup[["par"]][["CODs"]][1])
                                       ifelse(length(MC_setup[["par"]][["NH4s"]]) > 1, NH4si <- sample(MC_setup[["par"]][["NH4s"]][i,],1), NH4si <- MC_setup[["par"]][["NH4s"]][1])
                                       ifelse(length(MC_setup[["par"]][["qf"]])   > 1, qfi   <- sample(MC_setup[["par"]][["qf"]][i,],1),   qfi   <- MC_setup[["par"]][["qf"]][1])
                                       ifelse(length(MC_setup[["par"]][["CODf"]]) > 1, CODfi <- sample(MC_setup[["par"]][["CODf"]][i,],1), CODfi <- MC_setup[["par"]][["CODf"]][1])
                                       ifelse(length(MC_setup[["par"]][["NH4f"]]) > 1, NH4fi <- sample(MC_setup[["par"]][["NH4f"]][i,],1), NH4fi <- MC_setup[["par"]][["NH4f"]][1])
                                       ifelse(length(MC_setup[["par"]][["CODr"]]) > 1, CODri <- sample(MC_setup[["par"]][["CODr"]][i,],1), CODri <- MC_setup[["par"]][["CODr"]][1])
                                       ifelse(length(MC_setup[["par"]][["NH4r"]]) > 1, NH4ri <- sample(MC_setup[["par"]][["NH4r"]][i,],1), NH4ri <- MC_setup[["par"]][["NH4r"]][1])
                                       
                                       input.user <- input(spatial = 0, zero = 1e-5, folder = system.file("shiny", package = "EmiStatR"),
                                                           cores = EmiStatR.cores,
                                                           ww = list(qs = qsi, CODs = CODsi, NH4s = NH4si), 
                                                           inf = list(qf= qfi, CODf = CODfi, NH4f = NH4fi),
                                                           rw = list(CODr = CODri, NH4r = NH4ri, stat = "Goesdorf"), P1 = P1, 
                                                           st = list(E1=E1), 
                                                           export = 0)
                                       
                                       lev2vol <- slot(obs, "lev2vol")
                                       lev2vol$vol[length(lev2vol$vol)] <- E1$V
                                       slot(obs, "lev2vol") <- lev2vol

                                       ##***********************************************************************************
                                       ## end
                                       ##***********************************************************************************
                                       
                                       ## executing simulation 
                                       # library(EmiStatR)  
                                       # sim  <- EmiStatR(x = input.user)
                                       
                                       eval  <- Validation_Quantity(x = input.user, y = obs)
                                       
                                       #-------------------------------------------------------------------------------------------------------
                                       # Monte-Carlo output time-series matrix  (to be modified according to model output)
                                       #-------------------------------------------------------------------------------------------------------
                                       # data1 <- sim[[1]]$out1
                                       # data1 <- data1[,2:ncol(data1)]
                                       # head(data1)
                                       mcVTank  <- eval[["evalVolT"]]   # V_Tank
                                       
                                       mc1 <- list(mcVTank=mcVTank, of = eval[["of.NSE"]], st=E1, inputUser=input.user)
                                       
                                       return(mc1)
                                     } # end of foreach loop for Monte-Carlo simulations
                   ) # end system.time
                   },
                   
                   # case 2: discrete sampling, orifice
                   {
                     id.c6 <- MC_setup$id.c6
                     discrete.samples <- matrix(NA, nrow = 1, ncol = (length(id.c6)-1))
                     
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
                               discrete.samples[id,] <- c(MC_setup$par[id.c6[2]][[1]][1], 
                                                          MC_setup$par[id.c6[3]][[1]][1], 
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
                     }  
                     
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
                     }
                     
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
                     }  
                     
                     if(length(id.c6) == 8){
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
                                   for(ig in 1:length(MC_setup$par[id.c6[8]][[1]])){
                                   if(ih == 1){ 
                                     discrete.samples[1,] <- c(MC_setup$par[id.c6[2]][[1]][1], MC_setup$par[id.c6[3]][[1]][1], 
                                                               MC_setup$par[id.c6[4]][[1]][1], MC_setup$par[id.c6[5]][[1]][1],
                                                               MC_setup$par[id.c6[6]][[1]][1], MC_setup$par[id.c6[7]][[1]][1],
                                                               MC_setup$par[id.c6[8]][[1]][1])
                                     ih <- ih+1
                                   }else{discrete.samples <- rbind(discrete.samples, c(MC_setup$par[id.c6[2]][[1]][ia], 
                                                                                       MC_setup$par[id.c6[3]][[1]][ib], 
                                                                                       MC_setup$par[id.c6[4]][[1]][ic],
                                                                                       MC_setup$par[id.c6[5]][[1]][id],
                                                                                       MC_setup$par[id.c6[6]][[1]][ie],
                                                                                       MC_setup$par[id.c6[7]][[1]][iff],
                                                                                       MC_setup$par[id.c6[8]][[1]][ig])
                                   )}
                                 }}
                               }}}}}
                       
                       colnames(discrete.samples) <- names(MC_setup$par[id.c6[2:8]])
                     } 
                     
                     nsamples <- nrow(discrete.samples)
                     print(paste0("starting discrete sampling calibration, nsamples = ", nsamples))
                     
                     obj1 <- 1:nsamples
                     a <- system.time(
                       sim1 <- foreach(obj1 = obj1, .packages = c("EmiStatR", "stUPscales"), .export=c("obj1"), # previously "AggProp" instead of "stUPscales"
                                       .errorhandling = "pass", .verbose=TRUE, .combine = "rbind") %dopar% {
                                         # i <- 1
                                         i <- obj1
                                         print(paste("simulation ", i, ", ", MC_setup[["mcCores"]], " cores, running in parallel mode...", sep=""))
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
                                         ifelse(length(MC_setup[["par"]]$tfS)  > 1, tfSi  <- round(discrete.samples[i,"tfS"]), tfSi <- round(MC_setup[["par"]]$tfS))
                                         ifelse(length(MC_setup[["par"]]$pe)   > 1, pei   <- discrete.samples[i,"pe"], pei <- MC_setup[["par"]]$pe)
                                         ifelse(length(MC_setup[["par"]]$Qd)   > 1, Qdi   <- discrete.samples[i,"Qd"], Qdi <- MC_setup[["par"]]$Qd)
                                         ifelse(length(MC_setup[["par"]]$V)    > 1, Vi    <- discrete.samples[i, "V"], Vi <- MC_setup[["par"]]$V)
                                         ifelse(length(MC_setup[["par"]]$lev.ini) > 1, lev.inii   <- discrete.samples[i,"lev.ini"], lev.inii <- MC_setup[["par"]]$lev.ini)
                                         ifelse(length(MC_setup[["par"]]$Dd)   > 1, ADi   <- discrete.samples[i, "Dd"], Ddi <- MC_setup[["par"]]$Dd)
                                         ifelse(length(MC_setup[["par"]]$Cd)   > 1, Cdi   <- discrete.samples[i, "Cd"], Cdi <- MC_setup[["par"]]$Cd)
                                         ifelse(length(MC_setup[["par"]]$Cimp)   > 1, Cimpi   <- discrete.samples[i, "Cimp"], Cimpi <- MC_setup[["par"]]$Cimp)
                                         ifelse(length(MC_setup[["par"]]$Cper)   > 1, Cperi   <- discrete.samples[i, "Cper"], Cperi <- MC_setup[["par"]]$Cper)
                                         
                                         E1 <- list(id = 1, ns = MC_setup[["par"]]$ns, nm = MC_setup[["par"]]$nm,
                                                    nc = MC_setup[["par"]]$nc, numc = MC_setup[["par"]]$numc,
                                                    use = MC_setup[["par"]]$use, Atotal = Atotali,
                                                    Aimp = Aimpi, tfS = tfSi,
                                                    pe = pei, Qd = Qdi,
                                                    V = Vi,
                                                    lev2vol = MC_setup[["par"]]$lev2vol,
                                                    lev.ini = lev.inii,
                                                    Dd = Ddi,
                                                    Cd = Cdi,
                                                    Cimp = Cimpi, Cper = Cperi)
                                         
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
                                                             folderOutput = MC_setup[["folderOutput"]], cores = EmiStatR.cores,
                                                             ww = list(qs = qsi, CODs = CODsi, NH4s = NH4si), 
                                                             inf = list(qf= qfi, CODf = CODfi, NH4f = NH4fi),
                                                             rw = list(CODr = CODri, NH4r = NH4ri, stat = "Goesdorf"), P1 = P1, 
                                                             st = list(E1=E1), 
                                                             export = 0)
                                         
                                         lev2vol <- slot(obs, "lev2vol")
                                         lev2vol$vol[length(lev2vol$vol)] <- E1$V
                                         slot(obs, "lev2vol") <- lev2vol
                                         
                                         ##***********************************************************************************
                                         ## end
                                         ##***********************************************************************************
                                         
                                         ## executing simulation 
                                         # library(EmiStatR)  
                                         # sim  <- EmiStatR(x = input.user)
                                         
                                         eval  <- Validation_Quantity(x = input.user, y = obs)
                                         
                                         #-------------------------------------------------------------------------------------------------------
                                         # Monte-Carlo output time-series matrix  (to be modified according to model output)
                                         #-------------------------------------------------------------------------------------------------------
                                         # data1 <- sim[[1]]$out1
                                         # data1 <- data1[,2:ncol(data1)]
                                         # head(data1)
                                         mcVTank  <- eval[["evalVolT"]]   # V_Tank
                                         
                                         mc1 <- list(mcVTank=mcVTank, of = eval[["of.NSE"]], st=E1, inputUser=input.user)
                                         
                                         return(mc1)
                                       } # end of foreach loop for Monte-Carlo simulations
                     ) # end system.time
                   }
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
            setwd(MC_setup$folderOutput)
            # save(sim1, file=paste("sim1.RData", sep="")) # uncomment to save file
            # save(mc, file=paste("mc.RData", sep="")) # uncomment to save file
            setwd(currentDir)
            
            ## closing progress bar
            # close(pb)
            print(paste("End of", nsim[length(nsim)], "Monte-Carlo simulations.",sep=" "))
            
            return(list(mc=mc, sim1=sim1))
          }
)
