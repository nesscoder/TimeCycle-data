#get arguments from command line
args <- commandArgs(trailingOnly = TRUE)

#load functions and required packages for GeneCycle
reqPkgs = c('GeneCycle', 'plyr')
for(pkg in reqPkgs){
  if(!require(pkg, character.only=TRUE)){
    install.packages(pkg)
  }
}
loadedPkgs <-lapply(reqPkgs,library, character.only = TRUE)

#set input file
data <- read.delim(args[2], row.names = 1)

#set timepoints
interval <- 2
len <- 48
reps <- 1
timePoints <- seq(from = 0, to = len, by = interval)

#remove point from each row depending on number
data <- as.data.frame(t(apply(data, 1, function(expr){
  toNA <- sample(x = 1:length(expr), size = args[4], replace = F)
  expr[toNA] <- NA
  expr
})))

#Run GeneCycle
GeneCyclePvals <- robust.spectrum(x = t(data),
                                  t = timePoints,
                                  periodicity.time = 24,
                                  algorithm = "regression")

#return raw and fdr corrected p-values
GeneCycleResults <- data.frame(Gene = rownames(data),
                               pVals = GeneCyclePvals,
                               pVals.adj = p.adjust(GeneCyclePvals, method = 'fdr'))

#remove all other data storage except results and initial data
rm(list=setdiff(ls(), c("GeneCycleResults","args")))
save.image(file = args[3])
