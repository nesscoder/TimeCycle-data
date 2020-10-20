source("~/Desktop/TimeCycle-data/Scripts/TimeCycleV3.R")
pkg = c('TDA', 'parallel','pROC', 'imputeTS')
loadedPkgs <- lapply(pkg,library, character.only = TRUE)

#set inputs and output file names
data <- read.delim("~/Desktop/TimeCycle-data/Data/Processed/2_48_NoiseLV_0.1_BioRep_1.txt", row.names = 1)

#set rep labels
repLabel <- rep(1,dim(data)[2])

#Run TimeCycle
TimeCycleResults <- TimeCycle(data = data, repLabel = repLabel, period = 24, cores = 6, minLag = 2, maxLag = 5, resamplings = 5000)

#check output
expected <- c(rep(1,7000),rep(0,4000))
pred <- as.numeric(as.vector(unlist(TimeCycleResults$pVals.adj)))
ROC <- roc(expected,pred)
plot(ROC)
ROC$auc



