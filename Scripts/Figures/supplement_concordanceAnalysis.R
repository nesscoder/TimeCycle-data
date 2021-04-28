## --------------------------------------- Getting Data ready for Plotting ---------------------------------------
library(plotrix)

# get RAIN and JTK FDR adjusted Results
load("~/Desktop/TimeCycle-data/Results/timeTrialAdjRealDataComplete.RData")
JTK_results <- resultsAdj[[3]]

RAIN_results <- resultsAdj[[4]]
RAIN_results <- apply(RAIN_results, 2, function(pVal) {
  p.adjust(pVal, "fdr")
})

#################################
# LOAD THE DATA INTO R GeneCycle
#################################
setwd("~/Desktop/TimeCycle-data/Results/Biological/GeneCycle/")
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
  df$pVals.adj
})
GeneCycle_results <- do.call(cbind, GeneCycle_results)
rownames(GeneCycle_results) <- data[[1]]$Gene
GeneCycle_results <- as.data.frame(GeneCycle_results)

# if the pvalue is zero, set to R machine percision (2.2e-16)
GeneCycle_results <- GeneCycle_results + (GeneCycle_results == 0) * .Machine$double.eps

setwd("~/Desktop/TimeCycle-data/Results/Biological/TimeCycle/")
all.files <- as.list(dir())

##########################
# LOAD THE DATA INTO R TimeCycle
##########################

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
TimeCycle_results <- lapply(data, function(df) {
  df$pVals.adj
})
TimeCycle_results <- do.call(cbind, TimeCycle_results)
rownames(TimeCycle_results) <- rownames(data[[1]])

###############################
# functions
###############################
getConcordantDiscontantRatio <- function(s1, s2, thres) {
  discordant <- sum(xor((s1 <= thres), (s2 <= thres)))
  concordant <- length(s1) - discordant
  percentConcordant <- concordant / (discordant + concordant)
  return(percentConcordant)
}

concordancePlot <- function(dataset1, dataset2, title) {
  # TimeCycle
  thresholds <- seq(0, 0.1, 0.001)
  output <- vector()
  lineWidth <- 2.5

  for (i in 1:length(thresholds)) {
    output[i] <- getConcordantDiscontantRatio(TimeCycle_results[, dataset1], TimeCycle_results[, dataset2], thresholds[i])
  }
  plot(
    x = thresholds,
    y = output,
    xlab = "FDR",
    ylab = "Percent Concordance",
    main = title,
    axes = FALSE,
    ylim = c(0.5, 1),
    col = "#ff7f00",
    type = "l",
    lwd = lineWidth
  )
  box(lwd = 1)
  axis(1)
  axis(2, at = c(0.5, 0.6, 0.7, 0.8, 0.9, 1.0), labels = c(0, 0.6, 0.7, 0.8, 0.9, 1.0))
  axis.break(axis = 2, breakpos = 0.55)

  # JTK_CYCLE
  for (i in 1:length(thresholds)) {
    output[i] <- getConcordantDiscontantRatio(JTK_results[, dataset1], JTK_results[, dataset2], thresholds[i])
  }
  lines(
    x = thresholds,
    y = output,
    col = "#377eb8",
    lwd = lineWidth
  )

  # RAIN
  for (i in 1:length(thresholds)) {
    output[i] <- getConcordantDiscontantRatio(RAIN_results[, dataset1], RAIN_results[, dataset2], thresholds[i])
  }
  lines(
    x = thresholds,
    y = output,
    col = "Black",
    lwd = lineWidth
  )

  # GENECycle
  for (i in 1:length(thresholds)) {
    output[i] <- getConcordantDiscontantRatio(GeneCycle_results[, dataset1], GeneCycle_results[, dataset2], thresholds[i])
  }
  lines(
    x = thresholds,
    y = output,
    col = "green3",
    lwd = lineWidth
  )

  legend("bottomleft",
    legend = c("TimeCycle", "JTK_CYCLE", "RAIN", "GeneCycle"),
    col = c("#ff7f00", "#377eb8", "Black", "green3"), lty = 1, lwd = 2, cex = 0.8, box.lwd = NA
  )
}

###############################
# Create Plot
###############################

pdf("~/Desktop/TimeCycle-data/Results/Figures/supplement_concordanceAnalysis.pdf",
  width = 10,
  height = 11
)

par(font = 2)

layout.matrix <- matrix(c(
  0, 0, 0, 0,
  1, 0, 0, 0,
  2, 4, 0, 0,
  3, 5, 6, 0
), byrow = T, nrow = 4, ncol = 4)

layout(mat = layout.matrix)

# pairwise comparison of concordance plots Hogenesch 2A, Hogenesch 2B, Hughes, Zhang
concordancePlot(2, 3, "Hogenesch 2A vs Hogenesch 2B")
concordancePlot(2, 8, "Hogenesch 2A vs Hughes")
concordancePlot(2, 11, "Hogenesch 2A vs Zhang")
concordancePlot(3, 8, "Hogenesch 2B vs Hughes")
concordancePlot(3, 11, "Hogenesch 2B vs Zhang")
concordancePlot(8, 11, "Hughes vs Zhang")

dev.off()
