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
## Functions
###############
getRepAvgedDataFrame <- function(data, repLabel) {

  # get Column Names
  colnames <- colnames(data)
  # remove replicate label from colnames if they exist
  colnames <- gsub(pattern = "_rep.", replacement = "", colnames)
  # get unique ZT time for each point
  colnames <- unique(as.numeric(unlist(regmatches(colnames, gregexpr("[[:digit:]]+\\.*[[:digit:]]*", colnames)))))
  # get average of Replicate Labels
  output <- t(apply(data, MARGIN = 1, FUN = function(geneExprRow) {
    averageReps(geneExprRow, repLabel)
  }))
  output <- as.data.frame(output)
  colnames(output) <- colnames
  return(output)
}

averageReps <- function(geneExpr, Reps) {
  splitBy <- unlist(sapply(1:length(Reps), function(seq) {
    rep(seq, Reps[seq])
  }))
  as.vector(unlist(lapply(split(geneExpr, splitBy), mean)))
}

###############
## Process Data
###############
data <- read.delim(args[2], row.names = 1)

#set timepoints
interval <- as.numeric(args[4])
len <- as.numeric(args[5])
reps <- as.numeric(args[6])

#set timepoints
timePoints <- seq(from = 0, to = len, by = interval)

#GeneCycle does not have a built in function to deal with replicates
#if their are replicates, average together by ZT time point
if(reps > 1){
  data <- getRepAvgedDataFrame(data = data , repLabel = rep(reps, length(timePoints)))
}

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
