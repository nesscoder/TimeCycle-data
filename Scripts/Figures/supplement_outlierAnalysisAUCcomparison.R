
## --------------------------------------- Generating AUC Comparison plot --------------------------------
library("pROC")
library("ggplot2")
library("tidyverse")
library("Rmisc")

## --------------------------------------- Preprocessing TimeCycle Data --------------------------------

#################################
# LOAD THE DATA INTO R TIME CYCLE
#################################
setwd("~/Desktop/TimeCycle-data/Results/OutlierAnalysis/TimeCycle/")
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

###############################
# GET AUC Value OF DATA
###############################
AUC_TimeCycle <- apply(TC_results, 2, FUN = function(exprResult) {
  expected <- c(rep(1, 7000), rep(0, 4000))
  pred <- as.numeric(as.vector(unlist(exprResult)))
  ROC <- roc(expected, pred)
  as.vector(ROC$auc)
})

## --------------------------------------- Preprocessing JTK_Cycle Data --------------------------------

##########################
# LOAD THE DATA INTO R JTK
##########################
setwd("~/Desktop/TimeCycle-data/Results/OutlierAnalysis/JTK/")
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

###############################
# GET AUC Value OF JTK_Results
###############################
AUC_JTK <- apply(JTK_results, 2, FUN = function(exprResult) {
  expected <- c(rep(1, 7000), rep(0, 4000))
  pred <- as.numeric(as.vector(unlist(exprResult)))
  ROC <- roc(expected, pred)
  as.vector(ROC$auc)
})

## ------------------------------------------- Preprocessing Sw1pers Data --------------------------------

##########################
# LOAD THE DATA INTO R Sw1pers
##########################
setwd("~/Desktop/TimeCycle-data/Results/OutlierAnalysis/Sw1pers/")
all.files <- as.list(dir())
# get Sw1pers Results
data <- lapply(all.files, function(x) {
  read.delim(x, quote = "")
})
names(data) <- all.files

###############################
# EXTRACT DATA FRAME OF INTEREST
###############################
# get the adjusted pvalues and make them into a single dataframe
sw1per_results <- lapply(data, function(df) {
  df$score
})
sw1per_results <- do.call(cbind, sw1per_results)
rownames(sw1per_results) <- data[[1]]$id
sw1per_results <- as.data.frame(sw1per_results)

###############################
# GET AUC Value OF Sw1pers
###############################
AUC_sw1pers <- apply(sw1per_results, 2, FUN = function(exprResult) {
  expected <- c(rep(1, 7000), rep(0, 4000))
  pred <- as.numeric(as.vector(unlist(exprResult)))
  ROC <- roc(expected, pred)
  as.vector(ROC$auc)
})

## ------------------------------------------- Preprocessing RAIN Data --------------------------------

#################################
# LOAD THE DATA INTO R RAIN
#################################
setwd("~/Desktop/TimeCycle-data/Results/OutlierAnalysis/RAIN/")
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

###############################
# GET AUC Value OF DATA
###############################
AUC_RAIN <- apply(RAIN_results, 2, FUN = function(exprResult) {
  expected <- c(rep(1, 7000), rep(0, 4000))
  pred <- as.numeric(as.vector(unlist(exprResult)))
  ROC <- roc(expected, pred)
  as.vector(ROC$auc)
})

## ------------------------------------------- Preprocessing GeneCycle --------------------------------

#################################
# LOAD THE DATA INTO R RAIN
#################################
setwd("~/Desktop/TimeCycle-data/Results/OutlierAnalysis/GeneCycle/")
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

###############################
# GET AUC Value OF DATA
###############################
AUC_GeneCycle <- apply(GeneCycle_results, 2, FUN = function(exprResult) {
  expected <- c(rep(1, 7000), rep(0, 4000))
  pred <- as.numeric(as.vector(unlist(exprResult)))
  ROC <- roc(expected, pred)
  as.vector(ROC$auc)
})


## ------------------------------------------- META DATA --------------------------------

##################################################
# LOAD META DATA INTO R SAME FOR JTK, TIMECYCLE, GeneCycle and Sw1pers
##################################################

meltDataForPlotting <- function(AUCscores, method) {
  exprLabels <- names(AUCscores)
  metaInfo <- regmatches(exprLabels, gregexpr("[[:digit:]]+\\.*[[:digit:]]*", exprLabels))
  metaInfo <- as.data.frame(do.call(rbind, metaInfo))
  metaInfo <- as.data.frame(apply(metaInfo, 2, function(x) {
    as.numeric(as.character(x))
  }))
  colnames(metaInfo) <- c("Interval", "Length", "NoiseLevel", "Replicates")
  metaInfo$Method <- rep(method, length(AUCscores))
  metaInfo$AUC <- as.numeric(unname(unlist(AUCscores)))
  metaInfo$numSamples <- (metaInfo$Length / metaInfo$Interval + 1) * metaInfo$Replicates
  return(metaInfo)
}



## --------------------------------------- Plot Results --------------------------------
# select points to be plotted
JTK_CYCLE <- meltDataForPlotting(AUC_JTK, "JTK_CYCLE")
RAIN <- meltDataForPlotting(AUC_RAIN, "RAIN")
TimeCycle <- meltDataForPlotting(AUC_TimeCycle, "TimeCycle")
SW1PERs <- meltDataForPlotting(AUC_sw1pers, "SW1PERs")
GeneCycle <- meltDataForPlotting(AUC_GeneCycle, "GeneCycle")


toPlot <- rbind(JTK_CYCLE, TimeCycle, RAIN, SW1PERs, GeneCycle)
dataToPlot <- as.data.frame(toPlot) %>%
  dplyr::filter(NoiseLevel != 0)

dataToPlot$Method <- factor(dataToPlot$Method, levels = c("JTK_CYCLE", "RAIN", "SW1PERs", "GeneCycle", "TimeCycle"))


pdf("~/Desktop/TimeCycle-data/Results/Figures/supplement_outlierAnalysisAUCcomp.pdf", width = 7 * 1.5 / 3, height = 3.5 * 2)

theme_set(new = theme_light())
theme_replace(
  panel.border = element_rect(colour = "black", fill = NA, size = 2),
  text = element_text(size = 16, face = "bold")
)

ggplot(
  data = dataToPlot,
  aes(
    x = as.numeric(numSamples),
    y = as.numeric(AUC),
    fill = Method,
    colour = Method
  )
) +
  geom_point() +
  geom_point(shape = 21, colour = "black", show.legend = F) +
  geom_smooth(method = lm, se = FALSE, size = 1) +
  scale_fill_manual(
    labels = c("JTK_CYCLE", "RAIN", "SW1PERS", "GeneCycle", "TimeCycle"),
    values = c("#377eb8", "black", "grey", "green3", "#ff7f00")
  ) +
  scale_colour_manual(
    labels = c("JTK_CYCLE", "RAIN", "SW1PERS", "GeneCycle", "TimeCycle"),
    values = c("#377eb8", "black", "grey", "green3", "#ff7f00")
  ) +
  ylim(0.5, 1) + 
  facet_grid(Interval ~ Length) +
  xlab("Number of Samples") +
  ylab("AUC") +
  guides(colour = guide_legend(override.aes = list(linetype = 0, shape = 22, size = 5)), fill = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom", legend.title = element_blank())
dev.off()

AUCscoreSummary<- dataToPlot %>%
  summarySE(measurevar = "AUC", groupvars = c("Method", "Length", "Interval", "Replicates")) %>%
  mutate(numSamples = (Length / Interval + 1) * Replicates)

# FDR Threshold
typeErrorsPlot <- function(df) {
  fdrThreshold <- seq(0, 1, 0.005)
  nl <- seq(0.1, 0.4, 0.1)
  
  df <- apply(df,2,function(pvals){p.adjust(pvals,method = "fdr")})
  
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
  dplyr::select(c("2_48_NoiseLV_0.1_BioRep_1.Rdata","2_48_NoiseLV_0.2_BioRep_1.Rdata","2_48_NoiseLV_0.3_BioRep_1.Rdata","2_48_NoiseLV_0.4_BioRep_1.Rdata"))

pdf(file = "~/Desktop/TimeCycle-data/Results/Figures/supplemental_OutlierAnalysisWithType1and2ErrorAnalysis.pdf", width = 9, height = 9)
par(mfrow = c(4, 4))
typeErrorsPlot(TC_results)
typeErrorsPlot(JTK_results[,colnames(TC_results)])
typeErrorsPlot(RAIN_results[,colnames(TC_results)])
typeErrorsPlot(GeneCycle_results[,colnames(TC_results)])
dev.off()

pdf(file = "~/Desktop/TimeCycle-data/Results/Figures/supplemental_OutlierAnalysisWithType1and2ErrorAnalysis_legend.pdf")
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
