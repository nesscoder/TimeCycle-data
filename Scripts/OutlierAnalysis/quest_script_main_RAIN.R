
args <- commandArgs(trailingOnly = TRUE)
library("rain")
print(c(args[1],args[2],args[3],args[4],args[5],args[6],args[7]))
project <- args[1]

options(stringsAsFactors=FALSE)
data <- read.delim(args[2],row.names = 1)

st <- system.time({
results <- rain(t(data), deltat = as.numeric(args[4]), period = as.numeric(args[5]), nr.series = as.numeric(args[6]), peak.border = c(as.numeric(args[7]),as.numeric(args[8])), verbose = F)
})
#remove all other data storage except results and initial data
rm(list=setdiff(ls(), c("st","results","args")))
save.image(file = args[3])
