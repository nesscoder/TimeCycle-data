## --------------------------------------- Generating false Positive Rate Plot --------------------------------
library("pROC")
library("ggplot2")
library("gridExtra")
library("tidyverse")

## --------------------------------------- Preprocessing TimeCycle Data --------------------------------

#################################
# LOAD THE DATA INTO R TIME CYCLE
#################################
setwd("~/Desktop/TimeCycle-data/Results/Synthetic/TimeCycle_comp/")
all.files <- as.list(dir())
# get Timecycle Results
count <- 0
data <- lapply(all.files, FUN = function(x) {
  count <<- count + 1
  load(all.files[[count]])
  attach(all.files[[count]], warn = FALSE)
  df <- get("TimeCycleResults", pos = paste0("file:", all.files[[count]]))
  assign(paste(all.files[[count]], sep = ""), df)
})
names(data) <- all.files

###############################
# EXTRACT DATA FRAME OF INTEREST
###############################
# get the adjusted pvalues and make them into a single dataframe
TC_results <- lapply(data, function(df) {
  df$pVals.adj
})
TC_results <- do.call(cbind, TC_results)
rownames(TC_results) <- rownames(data[[1]])
TC_results <- as.data.frame(TC_results)

## --------------------------------------- Preprocessing JTK Data --------------------------------

#################################
# LOAD THE DATA INTO R TIME CYCLE
#################################
setwd("~/Desktop/TimeCycle-data/Results/Synthetic/JTK/")
all.files <- as.list(dir())
# get JTK Results
count <- 0
data <- lapply(all.files, FUN = function(x) {
  count <<- count + 1
  load(all.files[[count]])
  attach(all.files[[count]], warn = FALSE)
  df <- get("results", pos = paste0("file:", all.files[[count]]))
  assign(paste(all.files[[count]], sep = ""), df)
})
names(data) <- all.files

###############################
# EXTRACT DATA FRAME OF INTEREST
###############################
# get the adjusted pvalues and make them into a single dataframe
JTK_results <- lapply(data, function(df) {
  df$ADJ.P
})
JTK_results <- do.call(cbind, JTK_results)
rownames(JTK_results) <- rownames(data[[1]])
JTK_results <- as.data.frame(JTK_results)

## --------------------------------------- Preprocessing RAIN Data --------------------------------

#################################
# LOAD THE DATA INTO R TIME CYCLE
#################################
setwd("~/Desktop/TimeCycle-data/Results/Synthetic/RAIN/")
all.files <- as.list(dir())
# get RAIN Results
count <- 0
data <- lapply(all.files, FUN = function(x) {
  count <<- count + 1
  load(all.files[[count]])
  attach(all.files[[count]], warn = FALSE)
  df <- get("results", pos = paste0("file:", all.files[[count]]))
  assign(paste(all.files[[count]], sep = ""), df)
})
names(data) <- all.files

###############################
# EXTRACT DATA FRAME OF INTEREST
###############################
# get the adjusted pvalues and make them into a single dataframe
RAIN_results <- lapply(data, function(df) {
  p.adjust(df$pVal, 'fdr')
})
RAIN_results <- do.call(cbind, RAIN_results)
rownames(RAIN_results) <- rownames(data[[1]])
RAIN_results <- as.data.frame(RAIN_results)

## ------------------------------------------- Preprocessing GeneCycle --------------------------------

#################################
# LOAD THE DATA INTO R RAIN
#################################
setwd("~/Desktop/TimeCycle-data/Results/Synthetic/GeneCycle/")
all.files <- as.list(dir())
# get GeneCycle Results
count <- 0
data <- lapply(all.files, FUN = function(x) {
  count <<- count + 1
  load(all.files[[count]])
  attach(all.files[[count]], warn = FALSE)
  df <- get("GeneCycleResults", pos = paste0("file:", all.files[[count]]))
  assign(paste(all.files[[count]], sep = ""), df)
})
names(data) <- all.files

###############################
# EXTRACT DATA FRAME OF INTEREST
###############################
# get the adjusted pvalues and make them into a single dataframe
GeneCycle_results <- lapply(data, function(df){
  df$pVals.adj
})
GeneCycle_results <- do.call(cbind, GeneCycle_results)
rownames(GeneCycle_results) <- rownames(data[[1]])
GeneCycle_results <- as.data.frame(GeneCycle_results)

## ------------------------------------------- Processing FDR Threshold --------------------------------

# FDR Threshold
typeErrorsPlot <- function(df) {
  fdrThreshold <- seq(0, 1, 0.005)
  nl <- seq(0.1, 0.4, 0.1)

  
  for (i in 1:length(nl)) {
    falsePositives <- sapply(fdrThreshold, function(fdr) {
      sum(df[7001:11000, i] <= fdr) / 4000
    })
    falseNegativeRate <- sapply(fdrThreshold, function(fdr) {
      sum(df[1:7000, i] > fdr) / 7000
    })
    plot(fdrThreshold,
      falsePositives,
      main = paste0("Noise Level: ", nl[i]),
      xlab = "FDR Threshold",
      ylab = "Rate (%)",
      type = "l",
      lwd = 3,
      ylim = c(0, 1),
      col = "#ff7f00"
    )
    points(fdrThreshold,
      falseNegativeRate,
      lwd = 3,
      type = "l",
      col = "#377eb8"
    )
    abline(0, 1)
    abline(v = 0.05, col = "gray", lwd = 2)
    box(lwd = 3)
  }
}

#remove noise level 0 from the analysis
TC_results <- TC_results %>%
  dplyr::select(-"2_48_NoiseLV_0_BioRep_1.Rdata")

pdf(file = "~/Desktop/TimeCycle-data/Results/Figures/supplemental_Type1and2ErrorAnalysis.pdf", width = 9, height = 9)
par(mfrow = c(4, 4))
typeErrorsPlot(TC_results)
typeErrorsPlot(JTK_results[,colnames(TC_results)])
typeErrorsPlot(RAIN_results[,colnames(TC_results)])
typeErrorsPlot(GeneCycle_results[,colnames(TC_results)])
dev.off()

pdf(file = "~/Desktop/TimeCycle-data/Results/Figures/supplemental_Type1and2ErrorAnalysis_legend.pdf")
par(mfrow = c(1, 1))
plot.new()
legend("center",
       legend = c("False Postives", "False Negatives"),
       col = c(
         "#ff7f00",
         "#377eb8"
       ),
       pch = 15,
       bty = "n",
       pt.cex = 2,
       cex = 1.2,
       text.col = "black",
       horiz = F,
       inset = c(0.1, 0.1)
)
dev.off()
