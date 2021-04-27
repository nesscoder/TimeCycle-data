## --------------------------------------- Getting Data ready for Plotting ---------------------------------------

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

concordancePlot <- function(dataset1, dataset2) {
  # TimeCycle
  thresholds <- seq(0, 0.2, 0.001)
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
    col = "#ff7f00",
    type = "l",
    lwd = lineWidth,
    ylim = c(0, 1)
  )
  box(lwd = 1.75)

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
}

###############################
# Create Plot
###############################

pdf("~/Desktop/TimeCycle-data/Results/Figures/supplement_conordanceAnalysis.pdf",
  width = 10,
  height = 11
)

par(font = 2, font.axis = 2, font.lab = 2)

layout.matrix <- matrix(c(
  7, 0, 11, 11,
  1, 8, 11, 11,
  2, 4, 9, 0,
  3, 5, 6, 10
), byrow = T, nrow = 4, ncol = 4)

layout(mat = layout.matrix)

#pairwise comparison of concordance plots Hogenesch 2A, Hogenesch 2B, Hughes, Zhang
concordancePlot(2, 3)
concordancePlot(2, 8)
concordancePlot(2, 11)
concordancePlot(3, 8)
concordancePlot(3, 11)
concordancePlot(8, 11)

plot.new()
text(0.5, 0.5, "Hogenesch 2A", font = 2, cex = 1.5)
plot.new()
text(0.5, 0.5, "Hogenesch 2B", font = 2, cex = 1.5)
plot.new()
text(0.5, 0.5, "Hughes", font = 2, cex = 1.5)
plot.new()
text(0.5, 0.5, "Zhang", font = 2, cex = 1.5)
plot.new()
legend(0, 0.5,
  legend = c("TimeCycle", "JTK_CYCLE", "RAIN", "GeneCycle"),
  col = c("#ff7f00", "#377eb8", "Black", "green3"), lty = 1, lwd = 4, cex = 2, box.lwd = NA
)
dev.off()
