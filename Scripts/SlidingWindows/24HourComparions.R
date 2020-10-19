#################################################
# Install packages and Dependencies
#################################################
reqPkgs <- c("tidyverse", "lubridate", "ggplot2", "gplots", "TDA", "parallel", "pROC", "signal", "pracma", "TSA", "imputeTS", "rain", "UpSetR")
for (pkg in reqPkgs) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
  }
}
loadedPkgs <- lapply(reqPkgs, library, character.only = TRUE)
# load functions and required packages for TimeCycle
source("~/Desktop/TimeCycleV3/Scripts/TimeCycleV3.R")
source("~/Desktop/TimeCycleV3/Scripts/Synthetic/JTK_CYCLEv3.1.R")

#################################################
# Load Data
#################################################
# Load the Data - Hogenesch 1H sampling
data_hogenesch <- read.delim("~/Desktop/TimeCycleV3/Data/RealData/Hogenesch2009.all.data.ann.txt", row.names = 1)
# Load the Data - Hughes 2H sampling
data_hughes <- read.delim("~/Desktop/TimeCycleV3/Data/RealData/Hughes2012.all.data.ann.txt", row.names = 1)
# Load the Data - Zhang 2H sampling
data_zhang <- read.delim("~/Desktop/TimeCycleV3/Data/RealData/Zhang2014.all.data.ann.txt", row.names = 1)


#################################################
# Functions for 24 Hour Comparison Analysis
#################################################
# function to select the correct ZT times from the dataframe, 13 points over 2 hours for the 24 hour cycle
selData <- function(data, ZT_start) {
  sel <- seq(from = ZT_start, by = 2, length.out = 13)
  return(data[, sel])
}

# Processing Data with TimeCycle
process_TC <- function(data) {
  repLabel <- rep(1, ncol(data))
  results <- TimeCycleV5(
    data = data,
    repLabel = repLabel,
    period = 24,
    laplacian = T,
    cores = 6,
    minLag = 2,
    maxLag = 3,
    experiments = 10000
  )
  return(results$pVals.adj)
}

# Processing Data withJTK_CYCLE
process_JTK <- function(data) {
  annot <- rownames(data)
  jtkdist(13, 1) # 13 total time points, 1 replicates per time point

  periods <- 12 # looking for rhythms between 20-28 hours (i.e. between 5 and 7 time points per cycle).
  jtk.init(periods, 2) # 4 is the number of hours between time points

  cat("JTK analysis started on", date(), "\n")
  flush.console()

  st <- system.time({
    res <- apply(data, 1, function(z) {
      jtkx(z)
      c(JTK.ADJP, JTK.PERIOD, JTK.LAG, JTK.AMP)
    })
    res <- as.data.frame(t(res))
    bhq <- p.adjust(unlist(res[, 1]), "BH")
    res <- cbind(bhq, res)
    colnames(res) <- c("BH.Q", "ADJ.P", "PER", "LAG", "AMP")
    results <- cbind(annot, res)
  })
  return(results$BH.Q)
}

# Processing Data with RAIN
process_RAIN <- function(data) {
  results <- rain(t(data),
    deltat = 2,
    period = 24,
    nr.series = 1,
    peak.border = c(0.2, 0.8),
    verbose = F
  )
  return(results$pVal)
}


##########################################
# Sliding Window 24 Hour Time Comparison
##########################################
# ZT times to cycle through starting at 1-24
ZT_startTimes <- 1:24

# Run Analysis on 24 hour windows using Timecycle and Save the Adjusted p-values
TC_results <- sapply(ZT_startTimes, function(startTime) {
  dataToProcess <- selData(data = data_hogenesch, ZT_start = startTime)
  process_TC(dataToProcess)
})
colnames(TC_results) <- paste0("TC_ZT_", ZT_startTimes)
rownames(TC_results) <- rownames(data_hogenesch)

# Run Analysis on 24 hour windows using JTK_CYCLE and Save the Adjusted p-values
JTK_results <- sapply(ZT_startTimes, function(startTime) {
  dataToProcess <- selData(data = data_hogenesch, ZT_start = startTime)
  process_JTK(dataToProcess)
})
colnames(JTK_results) <- paste0("JTK_ZT_", ZT_startTimes)
rownames(JTK_results) <- rownames(data_hogenesch)

# Run Analysis on 24 hour windows using RAIN and Save the Adjusted p-values
RAIN_results <- sapply(ZT_startTimes, function(startTime) {
  dataToProcess <- selData(data = data_hogenesch, ZT_start = startTime)
  process_RAIN(dataToProcess)
})
colnames(RAIN_results) <- paste0("RAIN_ZT_", ZT_startTimes)
rownames(RAIN_results) <- rownames(data_hogenesch)


##########################################
# 24 Hour Same ZT Time Comparison
##########################################
ZT_Times_Overlap <- colnames(selData(data = data_hogenesch, ZT_start = 1))

Zhang <- data_zhang[, ZT_Times_Overlap]
Hughes <- data_hughes[, ZT_Times_Overlap]
Hogenesch <- data_hogenesch[, ZT_Times_Overlap]

RAIN_results_Dataset_Comp <- data.frame(
  Zhang = process_RAIN(Zhang),
  Hughes = process_RAIN(Hughes),
  Hogenesch = process_RAIN(Hogenesch)
)

JTK_results_Dataset_Comp <- data.frame(
  Zhang = process_JTK(Zhang),
  Hughes = process_JTK(Hughes),
  Hogenesch = process_JTK(Hogenesch)
)

TC_results_Dataset_Comp <- data.frame(
  Zhang = process_TC(Zhang),
  Hughes = process_TC(Hughes),
  Hogenesch = process_TC(Hogenesch)
)


rm(list = setdiff(ls(), c("TC_results", "JTK_results", "RAIN_results", "RAIN_results_Dataset_Comp", "JTK_results_Dataset_Comp", "TC_results_Dataset_Comp", "data_hogenesch", "Hughes", "Hogenesch", "Zhang")))
#save.image("~/Desktop/TimeCycleV3/Results/Downsampling/24HourComparison.Rdata")
