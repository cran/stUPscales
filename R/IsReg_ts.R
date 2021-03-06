# Plot routine for Harmonised data
# author: J.A. Torres-Matallana
# organization: Luxembourg Institute of Science and Technology (LIST), Luxembourg
#               Wagenigen University and Research Centre (WUR), Wageningen, The Netherlands   
# date: 21.05.2015 

IsReg.ts <- function(data, format, tz){
  #=======================================================================================================
  # creating a xts object from data
  #=======================================================================================================
  requireNamespace("xts")
  
  dat  <- data[,1]
  val  <- data[,-1] 
  #date <- strptime(dat, format=format, tz=tz) # use if duration is less than 1 day
  date <- as.POSIXct(dat, format=format, tz = tz)
  # save(date, file="date.RData")
  ts   <- as.xts(val, order.by=date)
  # head(ts, 20)
  #=======================================================================================================
  # determining regularity of time series
  #=======================================================================================================
  ifelse(is.regular(ts, strict=TRUE)=="TRUE", reg <- "_TSregular", reg <- "_TSirregular")
  return(list(reg, ts))
}