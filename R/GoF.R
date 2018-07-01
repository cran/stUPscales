# goodness-of-fit from data.frame
# author: J.A. Torres-Matallana
# organization: Luxembourg Institute of Science and Technology (LIST), Luxembourg
#               Wagenigen University and Research Centre (WUR), Wageningen, The Netherlands   
# date: 14.07.2015 - 14.07.2015

GoF <- function(eval, col_sim, col_obs, name){
  requireNamespace("hydroGOF")
  gof1 <- gof(sim = as.numeric(eval[,col_sim]), obs = as.numeric(eval[,col_obs]), digits=9) # should be gof(sim, obs)
  
  colnames(gof1) <- c("")
  
  # if(name!=""){
  if(!missing(name)){
    write.table(gof1, file = paste("gof",name,".csv", sep=""), 
                sep = ",", qmethod = "double", row.names=TRUE, col.names=FALSE)
    save(gof1, file=paste("gof",name,".RData", sep=""))  
  }
  return(gof1)
}
