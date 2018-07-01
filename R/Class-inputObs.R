# class: input for observations
# author: J.A. Torres-Matallana
# organization: Luxembourg Institute of Science and Technology (LIST), Luxembourg
#               Wagenigen University and Research Centre (WUR), Wageningen, The Netherlands   
# date: 22.04.2016 - 28.07.2017

inputObs <- setClass("inputObs", 
         slots = c(
           id = "numeric",
           plot = "numeric",
           delta = "list",
           observations = "list",
           lev2vol = "list",
           namePlot = "character",
           legendPosition = "list",
           var = "character"),
         prototype = list(
           id = 1,
           plot = 1,
           delta = NULL,
           observations = NULL,
           lev2vol = NULL,
           namePlot = "Validation plot",
           legendPosition = NULL,
           var = NULL)
         )
