
args <- commandArgs(trailingOnly = TRUE)

source("JTK_CYCLEv3.1.R")

project <- args[1]

options(stringsAsFactors=FALSE)
data <- read.delim(args[2])
annot <- data[,1]

rownames(data) <- data[,1]
data <- data[,-1]
jtkdist(as.numeric(args[4]), as.numeric(args[5]))       # 13 total time points, 2 replicates per time point

periods <- as.numeric(args[6]):as.numeric(args[7])      # looking for rhythms between 20-28 hours (i.e. between 5 and 7 time points per cycle).
jtk.init(periods,as.numeric(args[8]))  #4 is the number of hours between time points

cat("JTK analysis started on",date(),"\n")
flush.console()

st <- system.time({
  res <- apply(data,1,function(z) {
    jtkx(z)
    c(JTK.ADJP,JTK.PERIOD,JTK.LAG,JTK.AMP)
  })
  res <- as.data.frame(t(res))
  bhq <- p.adjust(unlist(res[,1]),"BH")
  res <- cbind(bhq,res)
  colnames(res) <- c("BH.Q","ADJ.P","PER","LAG","AMP")
  results <- cbind(annot,res)
  #results <- results[order(res$ADJ.P,-res$AMP),]
})
print(st)
results
sum(results$ADJ.P<0.5)

#remove all other data storage except results and initial data
rm(list=setdiff(ls(), c("st","results","args")))
save.image(file = args[3])
