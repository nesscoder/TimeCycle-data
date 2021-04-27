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
  df$pVals
})
TC_results <- do.call(cbind, TC_results)
rownames(TC_results) <- rownames(data[[1]])
TC_results <- as.data.frame(TC_results)

## --------------------------------------- Preprocessing JTK_Cycle Data --------------------------------

##########################
# LOAD THE DATA INTO R JTK
##########################
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

## ------------------------------------------- Preprocessing RAIN Data --------------------------------

#################################
# LOAD THE DATA INTO R RAIN
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
  df$pVal
})
RAIN_results <- do.call(cbind, RAIN_results)
rownames(RAIN_results) <- rownames(data[[1]])
RAIN_results <- as.data.frame(RAIN_results)

## ------------------------------------------- Preprocessing GeneCycle --------------------------------

#################################
# LOAD THE DATA INTO R GeneCycle
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
GeneCycle_results <- lapply(data, function(df) {
  df$pVals
})
GeneCycle_results <- do.call(cbind, GeneCycle_results)
rownames(GeneCycle_results) <- rownames(data[[1]])
GeneCycle_results <- as.data.frame(GeneCycle_results)


## ------------------------------------------- Generate Plots --------------------------------

pdf("~/Desktop/TimeCycle-data/Results/Figures/supplement_pvalDistributions.pdf",
    width = 10,
    height = 8
)

par(font = 2, font.axis = 2, font.lab = 2)

layout.matrix <- matrix(c(
  1, 2, 3, 4,
  5, 6, 7, 8,
  9, 10, 11, 12,
  13, 14, 15, 16
), byrow = T, nrow = 4, ncol = 4)

layout(mat = layout.matrix)

histPlots <- function(data, col, title ){
  hist(data, xlab = 'p-value',col = col, main = title, right = T)
  abline(h = 400, lty =2)
}

# For plotting purposes any pval >= 1 set to .999999 to prevent left sided binning
TC_results[(TC_results >= 1)] <- TC_results[(TC_results >= 1)] - .000001

histPlots(TC_results$`2_48_NoiseLV_0.1_BioRep_1.Rdata`[7001:11000], "#ff7f00", title = "Noise Level: 0.1")
histPlots(TC_results$`2_48_NoiseLV_0.2_BioRep_1.Rdata`[7001:11000], "#ff7f00", title = "Noise Level: 0.2")
histPlots(TC_results$`2_48_NoiseLV_0.3_BioRep_1.Rdata`[7001:11000], "#ff7f00", title = "Noise Level: 0.3")
histPlots(TC_results$`2_48_NoiseLV_0.4_BioRep_1.Rdata`[7001:11000], "#ff7f00", title = "Noise Level: 0.4")

histPlots(JTK_results$`2_48_NoiseLV_0.1_BioRep_1.Rdata`[7001:11000], "#377eb8", title = "Noise Level: 0.1")
histPlots(JTK_results$`2_48_NoiseLV_0.2_BioRep_1.Rdata`[7001:11000], "#377eb8", title = "Noise Level: 0.2")
histPlots(JTK_results$`2_48_NoiseLV_0.3_BioRep_1.Rdata`[7001:11000], "#377eb8", title = "Noise Level: 0.3")
histPlots(JTK_results$`2_48_NoiseLV_0.4_BioRep_1.Rdata`[7001:11000], "#377eb8", title = "Noise Level: 0.4")

histPlots(RAIN_results$`2_48_NoiseLV_0.1_BioRep_1.Rdata`[7001:11000], "Black", title = "Noise Level: 0.1")
histPlots(RAIN_results$`2_48_NoiseLV_0.2_BioRep_1.Rdata`[7001:11000], "Black", title = "Noise Level: 0.2")
histPlots(RAIN_results$`2_48_NoiseLV_0.3_BioRep_1.Rdata`[7001:11000], "Black", title = "Noise Level: 0.3")
histPlots(RAIN_results$`2_48_NoiseLV_0.4_BioRep_1.Rdata`[7001:11000], "Black", title = "Noise Level: 0.4")

histPlots(GeneCycle_results$`2_48_NoiseLV_0.1_BioRep_1.Rdata`[7001:11000], "green3", title = "Noise Level: 0.1")
histPlots(GeneCycle_results$`2_48_NoiseLV_0.2_BioRep_1.Rdata`[7001:11000], "green3", title = "Noise Level: 0.2")
histPlots(GeneCycle_results$`2_48_NoiseLV_0.3_BioRep_1.Rdata`[7001:11000], "green3", title = "Noise Level: 0.3")
histPlots(GeneCycle_results$`2_48_NoiseLV_0.4_BioRep_1.Rdata`[7001:11000], "green3", title = "Noise Level: 0.4")

dev.off()
