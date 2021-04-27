## --------------------------------------- Generating MissingData plot --------------------------------
library("pROC")
library("ggplot2")
library("gridExtra")
library("tidyverse")
library("grid")
library("gridBase")

## --------------------------------------- Preprocessing TimeCycle Data --------------------------------

#################################
# LOAD THE DATA INTO R TIME CYCLE
#################################
setwd("~/Desktop/TimeCycle-data/Results/MissingData/TimeCycle/")
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

## ------------------------------------------- META DATA

##############################################################
# LOAD META DATA INTO R SAME FOR TIMECYCLE
##############################################################

metaInfo <- regmatches(names(AUC_TimeCycle), gregexpr("[[:digit:]]+\\.*[[:digit:]]*", names(AUC_TimeCycle)))
metaInfo <- do.call(rbind, metaInfo)
noiseLevels <- as.numeric(metaInfo[, 3])
missingness <- as.numeric(metaInfo[, 5])

plotData <- data.frame(AUC_TimeCycle, noiseLevels, missingness)

## --------------------------------------- Plot Results --------------------------------

pdf(file = "~/Desktop/TimeCycle-data/Results/Figures/missingness.pdf", width = 10, height = 8)
# plot AUC CURVES and ROC curves layout
layout(matrix(
  data = c(
    1, 1, 1, 1,
    1, 1, 1, 1,
    2, 3, 4, 5,
    6, 6, 6, 6
  ),
  ncol = 4,
  byrow = T
))

# colors for Methods
col_pal_fill <- c("grey", "#fff800", "#f37200", "#cf0000", "#000000")
col_pal_outline <- rep("black", 5)
trnplvl <- 1 # transparency Level

# set up plotting area
par(font = 2)
plot(
  x = NULL,
  y = NULL,
  pch = 16,
  cex = 1,
  xlim = c(0, max(missingness) + 1),
  ylim = c(0, 1),
  main = "AUC Score as a Measure of Percent Missingness",
  xlab = "Percent Missingness Per Gene",
  ylab = "AUC",
  font.axis = 2,
  font.lab = 2,
  axes = F
)
abline(v = 12, col = "orange", lwd = 5)
axis(side = 1, lwd = 2, tick = 2, at = seq(0, 24, 6), labels = seq(0, 24, 6) / 24)
axis(side = 2, lwd = 2)
box(lwd = 3)

# plot TC results
for (i in unique(plotData$noiseLevels)[-1]) {
  x <- dplyr::filter(plotData, noiseLevels == i) %>%
    dplyr::arrange(missingness) %>%
    dplyr::select(missingness)

  y <- dplyr::filter(plotData, noiseLevels == i) %>%
    dplyr::arrange(missingness) %>%
    dplyr::select(AUC_TimeCycle)


  x <- as.vector(unlist(x))
  y <- as.vector(unlist(y))


  points(
    x = x,
    y = y,
    bg = alpha(col_pal_fill[which(unique(noiseLevels) == i)], trnplvl),
    col = col_pal_outline[1],
    type = "o",
    pch = 21,
    cex = 2
  )
}

# add Legends
legend("bottomleft",
  legend = unique(noiseLevels)[-1],
  fil = col_pal_fill[-1],
  title = "Noise Level",
  bg = NA,
  box.lwd = 0
)


# ROCplot of pVals
# extract uncorrected pvalues
pVal <- lapply(data, function(df) {
  df$pVals
})
pVal <- do.call(cbind, pVal)
rownames(pVal) <- rownames(data)

# extract metaData
metaData <- regmatches(colnames(pVal), gregexpr("[[:digit:]]+\\.*[[:digit:]]*", colnames(pVal)))
metaData <- lapply(metaData, function(x) {
  as.numeric(as.character(x))
})
metaData <- do.call(rbind, metaData)

rocPlots <- sapply(unique(plotData$noiseLevels)[-1], function(NLtoCheck) {

  # select pVals of at specific Noise Level
  keep <- which(plotData[, 2] == NLtoCheck)
  sortKeep <- order(plotData[keep, 3])
  keep <- keep[sortKeep]

  # set coloring
  count <- 0
  pal <- colorRampPalette(c("#377eb8", "grey"))
  color_pal <- pal(length(sortKeep))
  p50 <- length(color_pal) / 2 # mark 50 percent missingness
  color_pal[p50] <- "#ff7f00" # special color of 50 percent

  # generate roc Plot
  expected <- c(rep(1, 7000), rep(0, 4000))
  apply(pVal[, keep], 2, function(pvals) {
    count <<- count + 1
    if (count == 1) {
      plot(roc(expected, pvals), col = color_pal[count], main = paste0("Noise Level: ", NLtoCheck))
    } else {
      plot(roc(expected, pvals), add = T, col = color_pal[count])
    }
  })
  rocCheck <<- roc(expected, pVal[, keep[p50]])
  plot(roc(expected, pVal[, keep[p50]]), add = T, lwd = 4, col = color_pal[p50])
})



# histograms of pVals
p <- list()
count <- 1
keep <- which(as.numeric(missingness) == 12)[-1]

for (i in keep) {
  TCResults <- data[[i]]
  ordering <- c("sin", "peak", "saw", "lTrend", "damp", "amped", "contract", "flat", "linear", "sigmoid", "exp")
  waveFormCols <- c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#FFFF33", "#a65628", "#cccccc", "#969696", "#636363", "#252525")
  wave <- gsub(pattern = "_.*", replacement = "", rownames(TCResults))
  TCResults$wave <- wave
  TCResults$wave <- fct_relevel(wave, ordering)
  TCResults$pVal_log <- -log(TCResults$pVals)

  p[[count]] <- ggplot(TCResults, aes(x = pVal_log, fill = wave, color = wave)) +
    geom_histogram() +
    scale_color_manual(values = waveFormCols) +
    scale_fill_manual(values = waveFormCols) +
    xlab(label = "-log(pVal)") +
    theme_light(base_rect_size = 2) +
    theme(
      legend.position = "none",
      panel.border = element_rect(colour = "black", fill = NA, size = 2)
    )
  count <- count + 1
}


lay <- matrix(c(
  7, 7, 7, 7,
  7, 7, 7, 7,
  7, 7, 7, 7,
  2, 3, 4, 5
),
nrow = 4, byrow = T
)
histPlts <- grid.arrange(grobs = p, layout_matrix = lay, newpage = FALSE)

# plot the histograms of each ggplot
plot.new()
vps <- baseViewports()
pushViewport(vps$figure)
vp1 <- plotViewport(c(1.8, 1, 0, 1))
print(histPlts, vp = vp1)
dev.off()
