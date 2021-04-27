#get arguments from command line
args <- commandArgs(trailingOnly = TRUE)

#load functions and required packages for TimeCycle
source("/home/emn6548/TimeCycleV3/Scripts/TimeCycleV3.R")
reqPkgs = c('plyr', 'TDAstats', 'parallel', 'imputeTS')
for(pkg in reqPkgs){
  if(!require(pkg, character.only=TRUE)){
    install.packages(pkg)
  }
}
loadedPkgs <-lapply(reqPkgs,library, character.only = TRUE)

#select inputs file
data <- read.delim(args[2], row.names = 1)

#remove point from each row depending on number
data <- as.data.frame(t(apply(data, 1, function(expr){
  toNA <- sample(x = 1:length(expr), size = args[7],replace = F)
  expr[toNA] <- NA
  expr
})))

#set rep labels
repLabel <- rep(args[5],args[4])

#Run TimeCycle
#set maxLag
TimeCycleResults <- TimeCycle(data = data, repLabel = repLabel, cores = 6, minLag = 2, maxLag = args[6], resamplings = 10000)

#set data to preprocessedData for use in analysis App
preprocessedData <- dataAvg
colnames(preprocessedData) <- paste0("ZT_",colnames(preprocessedData))

#remove all other data storage except results and initial data
rm(list=setdiff(ls(), c("TimeCycleResults","preprocessedData", "args")))
save.image(file = args[3])