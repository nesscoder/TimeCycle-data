## --------------------------------------- Generating pVal Comparison plot --------------------------------
library("pROC")
library("ggplot2")
library("gridExtra")
library("tidyverse")
library("rlang")
library("Cairo")

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
# get the adjusted pvalues and make them into a singel dataframe
TC_results <- lapply(data, function(df) {
  df$pVals
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
# get the adjusted pvalues and make them into a singel dataframe
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
# get the adjusted pvalues and make them into a singel dataframe
RAIN_results <- lapply(data, function(df) {
  df$pVal
})
RAIN_results <- do.call(cbind, RAIN_results)
rownames(RAIN_results) <- rownames(data[[1]])
RAIN_results <- as.data.frame(RAIN_results)

## --------------------------------------- Preprocessing Sw1pers Data --------------------------------

##########################
# LOAD THE DATA INTO R Sw1pers
##########################
setwd("~/Desktop/TimeCycle-data/Results/Synthetic/Sw1pers/")
all.files <- as.list(dir())
# get Sw1pers Results
data <- lapply(all.files, function(x) {
  read.delim(x, quote = "")
})
names(data) <- gsub(pattern = "_sw_scores.txt", replacement = ".Rdata", x = all.files)

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
#if the pvalue is zero, set to R machine percision (2.2e-16)
GeneCycle_results <- GeneCycle_results + (GeneCycle_results == 0)*.Machine$double.eps


# keep only what is needed for plotting
rm(list = setdiff(ls(), c("JTK_results", "RAIN_results", "TC_results", "sw1per_results","GeneCycle_results")))


## --------------------------------------- Plot Functions --------------------------------

best_ROC <- function(data) {
  expected <- c(rep(1, 7000), rep(0, 4000))
  pred <- data
  ROC <- roc(expected, pred)
  out <- coords(ROC, "best", ret = "threshold", transpose = F, best.method = "youden")
  as.vector(unlist(out))
}

plot_scatter <- function(data, xAxis = "Method 1", yAxis = "Method 2", coloring = "wave", xAxisLab = NULL, yAxisLab = NULL, waveFormCols = waveFormCols) {

  # log Transform the data
  data <- -log(data)
  #remove p-values that are zero
  data <- data*is.finite(as.matrix(data))

  wave <- gsub(pattern = "_.*", replacement = "", rownames(data))
  data$wave <- wave
  data$wave <- fct_relevel(wave, ordering)


  x <- enquo(xAxis)
  y <- enquo(yAxis)

  data_mean <- data %>%
    group_by(wave) %>%
    summarise(x = mean(!!x, na.rm = T),
              y = mean(!!y, na.rm = T))


  if (is.null(xAxisLab)) {
    xAxisLab <- xAxis
  }
  if (is.null(yAxisLab)) {
    yAxisLab <- yAxis
  }

  ggplot(data, aes_string(x = xAxis, y = yAxis, fill = coloring, color = coloring)) +
    geom_point(alpha = .4, size = 0.1) +
    geom_point(data = data_mean, aes(x = x, y = y), size = 4, shape = 21, color = "black", alpha = .9) +
    scale_color_manual(values = waveFormCols) +
    scale_fill_manual(values = waveFormCols) +
    xlab(xAxisLab) +
    ylab(yAxisLab) +
    theme(legend.position = "none")
}

plot_hist <- function(data, xAxis = NULL, coloring = "wave", xAxisLab = NULL, yAxisLab = NULL, waveFormCols = waveFormCols) {

  # log Transform the data
  data <- -log(data)
  #remove p-values that are zero
  data <- data*is.finite(as.matrix(data))

  wave <- gsub(pattern = "_.*", replacement = "", rownames(data))
  data$wave <- wave
  data$wave <- fct_relevel(wave, ordering)

  x <- enquo(xAxis)

  if (is.null(xAxisLab)) {
    xAxisLab <- xAxis
  }

  p <- ggplot(data, aes_string(x = xAxis, fill = coloring, color = coloring)) +
    geom_histogram() +
    xlab(xAxisLab) +
    scale_color_manual(values = waveFormCols) +
    scale_fill_manual(values = waveFormCols) +
    theme(legend.position = "none")

  hcenter <- max(data[, xAxis],na.rm = T) / 2
  vcenter <- layer_scales(p)$y$range$range[2] / 2
  p + annotate(geom = "text", x = hcenter, y = vcenter, label = xAxisLab, size = 10, fontface = "bold")
}


plot_percentCorrect <- function(data, bottomAxis = "Method 1", topAxis = "Method 2", threshold_xAxis = NULL, threshold_yAxis = NULL, coloring = "wave", xAxisLab = NULL, yAxisLab = NULL, waveFormCols = waveFormCols) {

  # log Transform the data
  wave <- gsub(pattern = "_.*", replacement = "", rownames(data))
  data$waves <- wave
  data$waves <- fct_relevel(wave, ordering)

  x <- enquo(bottomAxis)
  y <- enquo(topAxis)

  if (is.null(threshold_xAxis)) {
    threshold_xAxis <- best_ROC(data[, bottomAxis])
  }
  if (is.null(threshold_yAxis)) {
    threshold_yAxis <- best_ROC(data[, topAxis])
  }

  per_correct_cycling <- data[1:7000, ] %>%
    group_by(waves) %>%
    summarise(
      method1 = sum(!!x <= threshold_xAxis) / 1000 * 100,
      method2 = sum(!!y <= threshold_yAxis) / 1000 * 100
    )

  per_correct_noncycling <- data[7001:11000, ] %>%
    group_by(waves) %>%
    summarise(
      method1 = sum(!!x > threshold_xAxis) / 1000 * 100,
      method2 = sum(!!y > threshold_yAxis) / 1000 * 100
    )

  correctClass <- ungroup(bind_rows(per_correct_cycling, per_correct_noncycling))
  data <- correctClass %>% gather(method, perCorrect, c(method1, method2))

  if (is.null(xAxisLab)) {
    xAxisLab <- bottomAxis
  }
  if (is.null(yAxisLab)) {
    yAxisLab <- topAxis
  }


  ggplot(data, aes(y = perCorrect, x = waves, color = interaction(waves, method), fill = interaction(waves, method)), group = method) +
    geom_histogram(position = "dodge", stat = "identity") +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = "Waveform Shape", y = "Correct Classification (%)") +
    scale_y_continuous(breaks = c(0, 25, 50, 75, 100)) +
    scale_fill_manual(values = c(rep("white", length(waveFormCols)), waveFormCols)) +
    scale_fill_manual(values = c(alpha(waveFormCols, 0.4), waveFormCols)) +
    scale_color_manual(values = c(waveFormCols, waveFormCols)) +
    coord_flip() +
    scale_x_discrete(limits = rev(levels(data$waves)), labels = rev(c("Sin", "Peak", "Saw", "Linear Trend", "Damped", "Amped", "Contractile", "Flat", "Linear", "Sigmoid", "Exponential"))) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.position = "none"
    )
}


## --------------------------------------- Plot Results --------------------------------

#Data to plot
toPlot <- "2_48_NoiseLV_0.1_BioRep_1.Rdata"

ordering <- c("sin", "peak", "saw", "lTrend", "damp", "amped", "contract", "flat", "linear", "sigmoid", "exp")
waveFormCols <- c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#FFFF33", "#a65628", "#cccccc", "#969696", "#636363", "#252525")


theme_set(new = theme_light())
theme_replace(
  panel.border = element_rect(colour = "black", fill = NA, size = 2),
  text = element_text(size = 16, face = "bold")
)

plottingResults <- data.frame(list(
  TimeCycle = TC_results[, toPlot],
  JTK = JTK_results[, toPlot],
  RAIN = RAIN_results[, toPlot],
  GeneCycle = GeneCycle_results[, toPlot],
  SW1PERS = sw1per_results[, toPlot]
))
rownames(plottingResults) <- rownames(TC_results)

TimeCycle <- "TimeCycle"
JTK <- "JTK"
RAIN <- "RAIN"
GeneCycle <- "GeneCycle"
SW1PERS <- "SW1PERS"


p1 <- plot_hist(
  data = plottingResults,
  xAxis = TimeCycle,
  waveFormCols = waveFormCols
)

p2 <- plot_hist(
  data = plottingResults,
  xAxis = JTK,
  waveFormCols = waveFormCols
)


p3 <- plot_hist(
  data = plottingResults,
  xAxis = RAIN,
  waveFormCols = waveFormCols
)

p4 <- plot_hist(
  data = plottingResults,
  xAxis = SW1PERS,
  waveFormCols = waveFormCols,
  xAxisLab = "SW1PERS*"
)

p5 <- plot_hist(
  data = plottingResults,
  xAxis = GeneCycle,
  waveFormCols = waveFormCols,
  xAxisLab = expression(paste("GeneCycle",""^"\u2020"))
)

p6 <- plot_scatter(
  data = plottingResults,
  xAxis = TimeCycle,
  yAxis = JTK,
  waveFormCols = waveFormCols
)

p7 <- plot_scatter(
  data = plottingResults,
  xAxis = TimeCycle,
  yAxis = RAIN,
  waveFormCols = waveFormCols
)

p8 <- plot_scatter(
  data = plottingResults,
  xAxis = JTK,
  yAxis = RAIN,
  waveFormCols = waveFormCols
)

p9 <- plot_scatter(
  data = plottingResults,
  xAxis = TimeCycle,
  yAxis = SW1PERS,
  waveFormCols = waveFormCols,
  yAxisLab = "SW1PERS*"
)

p10 <- plot_scatter(
  data = plottingResults,
  xAxis = JTK,
  yAxis = SW1PERS,
  waveFormCols = waveFormCols,
  yAxisLab = "SW1PERS*"
)

p11 <- plot_scatter(
  data = plottingResults,
  xAxis = RAIN,
  yAxis = SW1PERS,
  waveFormCols = waveFormCols,
  yAxisLab = "SW1PERS*"
)

p12 <- plot_scatter(
  data = plottingResults,
  xAxis = TimeCycle,
  yAxis = GeneCycle,
  waveFormCols = waveFormCols,
  yAxisLab = expression("GeneCycle"^"\u2020")
)

p13 <- plot_scatter(
  data = plottingResults,
  xAxis = JTK,
  yAxis = GeneCycle,
  waveFormCols = waveFormCols,
  yAxisLab = expression("GeneCycle"^"\u2020")
)

p14 <- plot_scatter(
  data = plottingResults,
  xAxis = RAIN,
  yAxis = GeneCycle,
  waveFormCols = waveFormCols,
  yAxisLab = expression("GeneCycle"^"\u2020")
)

p15 <- plot_scatter(
  data = plottingResults,
  xAxis = SW1PERS,
  yAxis = GeneCycle,
  waveFormCols = waveFormCols,
  xAxisLab = "SW1PERS*",
  yAxisLab = expression("GeneCycle"^"\u2020")
)

p16 <- plot_percentCorrect(data = plottingResults, topAxis = TimeCycle, bottomAxis = JTK, waveFormCols = waveFormCols)
p17 <- plot_percentCorrect(data = plottingResults, topAxis = TimeCycle, bottomAxis = RAIN, waveFormCols = waveFormCols)
p18 <- plot_percentCorrect(data = plottingResults, topAxis = TimeCycle, bottomAxis = SW1PERS, waveFormCols = waveFormCols)
p19 <- plot_percentCorrect(data = plottingResults, topAxis = TimeCycle, bottomAxis = GeneCycle, waveFormCols = waveFormCols)

p20 <- plot_percentCorrect(data = plottingResults, topAxis = JTK, bottomAxis = RAIN, waveFormCols = waveFormCols)
p21 <- plot_percentCorrect(data = plottingResults, topAxis = JTK, bottomAxis = SW1PERS, waveFormCols = waveFormCols)
p22 <- plot_percentCorrect(data = plottingResults, topAxis = JTK, bottomAxis = GeneCycle, waveFormCols = waveFormCols)

p23 <- plot_percentCorrect(data = plottingResults, topAxis = RAIN, bottomAxis = SW1PERS, waveFormCols = waveFormCols)
p24 <- plot_percentCorrect(data = plottingResults, topAxis = RAIN, bottomAxis = GeneCycle, waveFormCols = waveFormCols)

p25 <- plot_percentCorrect(data = plottingResults, topAxis = SW1PERS, bottomAxis = GeneCycle, waveFormCols = waveFormCols)

# plot the layout
lay <- rbind(
  c(1, 16, 17, 18, 19),
  c(6, 2, 20, 21, 22),
  c(7, 8, 3, 23, 24),
  c(9, 10, 11, 4, 25),
  c(12, 13, 14, 15,5)
)

cairo_pdf(paste0("~/Desktop/TimeCycle-data/Results/Figures/", gsub(pattern = ".Rdata", replacement = ".pdf", toPlot)),
  width = 16, height = 15
)

grid.arrange(grobs = list(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16,p17,p18,p19,p20,p21,p22,p23,p24,p25), layout_matrix = lay)
dev.off()

