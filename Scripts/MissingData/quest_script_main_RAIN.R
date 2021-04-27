
args <- commandArgs(trailingOnly = TRUE)
library("rain")
print(c(args[1],args[2],args[3],args[4],args[5],args[6],args[7],args[8]))
project <- args[1]

options(stringsAsFactors=FALSE)
data <- read.delim(args[2],row.names = 1)

#remove point from each row depending on number
data <- as.data.frame(t(apply(data, 1, function(expr){
  toNA <- sample(x = 1:length(expr), size = args[8],replace = F)
  expr[toNA] <- NA
  expr
})))

st <- system.time({
results <- rain(t(data), deltat = as.numeric(args[4]), period = as.numeric(args[5]), peak.border = c(as.numeric(args[6]),as.numeric(args[7])), verbose = F)
})
#remove all other data storage except results and initial data
rm(list=setdiff(ls(), c("st","results","args")))
save.image(file = args[3])