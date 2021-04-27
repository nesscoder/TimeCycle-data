#get arguments from command line
args <- commandArgs(trailingOnly = TRUE)

#load functions and required packages for GeneCycle
reqPkgs = c('GeneCycle', 'plyr')
for(pkg in reqPkgs){
  if(!require(pkg, character.only=TRUE)){
    install.packages(pkg)
  }
}
loadedPkgs <- lapply(reqPkgs,library, character.only = TRUE)

###############
## Process Data
###############
data <- read.delim(args[2], row.names = 1)

#set timepoints
timePoints <- as.numeric(gsub(pattern = "ZT",replacement = "",x =  colnames(data)))

#Run GeneCycle
GeneCyclePvals <- robust.spectrum(x = t(data),
                                  t = timePoints,
                                  periodicity.time = 24,
                                  noOfPermutations = 5000,
                                  algorithm = "regression")

#return raw and fdr corrected p-values
GeneCycleResults <- data.frame(Gene = rownames(data),
                               pVals = GeneCyclePvals,
                               pVals.adj = p.adjust(GeneCyclePvals, method = 'fdr'))

#remove all other data storage except results and initial data
rm(list=setdiff(ls(), c("GeneCycleResults","args")))
save.image(file = args[3])
