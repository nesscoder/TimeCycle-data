library(UpSetR)
library(ggplot2)
library(gplots)
library(tidyverse)
library(gridExtra)
library(ggpubr)

setwd("~/Desktop/TimeCycleV3/Results/Biological/")

## --------------------------------------- Function for generating upset plot --------------------------------

##########################
# Input: Df of pvalues
# Description:
#   gets list of genes that are statistically significant and are above a Amp threshold and within a period threshold
# Output: List of Signficantly Cycling Genes in each of the 13 datasets
##########################

getSigGeneList <- function(df, pval) {
  colnames <- c("Hogenesch_1", "Hogenesch_2A", "Hogenesch_2B", "Hogenesch_4A", "Hogenesch_4B", "Hogenesch_4C", "Hogenesch_4D", "Hughes_2", "Hughes_4A", "Hughes_4B", "Zhang_2", "Zhang_4A", "Zhang_4B")

  names <- rownames(df)
  # genes That have NA for circadaian Period or Amp are removed from analysis
  geneLists <- apply(df < pval, 2, function(x) names[x])
  names(geneLists) <- colnames
  return(geneLists)
}

##########################
# Input: Df of genes
# Description:
#   Creates an upset plot that marks the circadian cycling genes by an orange triangle
# Output: Upset Plot of Gene Interactions
##########################
plotUpset <- function(df, showPlot = T) {
  methodNames <- c("TimeCycle")
  # identify the Cycling Genes in the Dataset
  Cycle <- c("Per1", "Per2", "Per3", "Cry1", "Cry2", "Npas2", "Clock", "Arntl", "Arntl2")
  Genes <- unique(unlist(df))
  collected <- as.numeric(Genes %in% Cycle)
  Cyclers <- c(0, 2)[collected + 1]

  # Bind Data together to Use for Plotting
  data <- cbind(fromList(df), Cyclers)

  dataQuery <- list(list(
    query = elements,
    params = list("Cyclers", 2),
    color = "orange",
    active = F,
    query.name = "Circadian Clock Genes"
  ))

  upset(data,
    nsets = 4,
    order.by = "freq",
    mb.ratio = c(0.6, 0.4),
    sets = c("Hogenesch_2A", "Hogenesch_2B", "Hughes_2", "Zhang_2"),
    keep.order = T,
    number.angles = 0,
    text.scale = 1.1,
    point.size = 2.8,
    set_size.show = T,
    line.size = 1,
    set_size.scale_max = 8000,
    query.legend = "bottom",
    empty.intersections = F,
    group.by = "degree",
    queries = dataQuery
  )
}

##########################
# Input: dataframe names, subset name to plot
# Description:
#   Creates an heatmap of significantly cycling genes
# Output: Upset Plot of Gene Interactions
##########################
plotHeatmap <- function(dataset, datasetName) {

  # select Representative heatmap of cycling genes
  # with the  significantly cycling genes as plotted in the upset Plot
  sigGenesToPlot <- timeSeries[[dataset]][sigGenes[[datasetName]], ]
  metaDatasigGenesToPlot <- data[[dataset]][sigGenes[[datasetName]], ]
  xVals <<- as.numeric(unlist(regmatches(colnames(sigGenesToPlot), gregexpr("[[:digit:]]+\\.*[[:digit:]]*", colnames(sigGenesToPlot)))))

  # creates a own color palette from lightgrey to red
  my_palette <- colorRampPalette(c("lightgrey", "#ff7f00"))(n = 3)

  # defines the color breaks manually for color transition
  col_breaks <- c(
    seq(-2, 0.0, length = 2), # for grey
    seq(0.0, 2, length = 2)
  ) # for red

  # detrend, mean center, and scale the data
  plotData <- apply(sigGenesToPlot, 1, function(y) {
    fit <- lm(y ~ xVals) # define linear trend
    yDetrend <- zapsmall(y - (xVals * fit$coefficients[2] + fit$coefficients[1]), 10) # subract trend from data
    yNormed <- yDetrend / max(yDetrend) # normalize data
    yNormed - mean(yNormed) # mean center
  })

  # plot genes in Phase order
  plotData <- t(plotData) # transpose matrix
  plotData <- plotData[order(metaDatasigGenesToPlot["Phase.in.Hours"]), ]
  colnames(plotData) <- gsub(pattern = "_AM", replacement = "", colnames(plotData))

  # plot sign heatmap of Data
  heatmap(sign(plotData),
    Colv = NA,
    Rowv = NA,
    col = my_palette,
    main = paste0("Significantly Cycling Genes: ", datasetName),
    breaks = col_breaks
  )
}
##########################
# Input: dataframe of  pvalues
# Description:
#   two functions used to generate the correlation and lower panel of pairwise scatter plot
# Output: Pairwise Scatter Plot of Data
##########################
# Correlation panel
panel.cor <- function(x, y) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- round(cor(x, y, method = "spearman"), digits = 2)
  text(0.5, 0.5, r, cex = 2)
}
# Customize lower panel
lower.panel <- function(x, y) {
  points(x, y, pch = 19, col = alpha("black", 0.1))
}


# Check to see Correlation between genes
getCorrelation <- function(controlData, comparisonData) {
  output <- sapply(1:nrow(controlData), function(i) {
    set1 <- unlist(unname(controlData[i, ]))
    set2 <- unlist(unname(comparisonData[i, ]))
    cor(set1, set2, method = "spearman")
  })
  return(output)
}

## --------------------------------------- Getting Data ready for Plotting ---------------------------------------

all.files <- as.list(dir())

##########################
# LOAD THE DATA INTO R
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

# get Timecycle Results
count <- 0
timeSeries <- lapply(all.files, FUN = function(x) {
  count <<- count + 1
  load(all.files[[count]])
  attach(all.files[[count]], warn = FALSE)
  df <- get("preprocessedData", pos = paste0("file:", all.files[[count]]))
  assign(paste(all.files[[count]], sep = ""), df)
})
names(timeSeries) <- all.files

###############################
# EXTRACT DATA FRAME OF INTEREST
###############################

# get the adjusted pvalues and make them into a singel dataframe
adjustpVal <- lapply(data, function(df) {
  df$pVals.adj
})
adjustpVal <- do.call(cbind, adjustpVal)
rownames(adjustpVal) <- rownames(data[[1]])


## --------------------------------------- Plot the data ---------------------------------------
# create the upset Plot for TimeCycle
sigGenes <- getSigGeneList(df = adjustpVal, pval = 0.05)
sigGenes <- lapply(sigGenes, na.omit) # remove NA genes, see getSigGeneList Function
# only plot list if it is significant
isSig <- as.vector(unlist(lapply(sigGenes, function(x) {
  length(x) > 1
})))
isSig[1] <- F


pdf(file = "~/Desktop/TimeCycleV3/Results/Figures/upsetPlot.pdf", width = 12, height = 8)
plotUpset(sigGenes[isSig])
dev.off()

# create plots of cycling genes in each Dataset
pdf(file = "~/Desktop/TimeCycleV3/Results/Figures/heatMaps.pdf", width = 4, height = 4)
dataNames <- c("Hogenesch_1", "Hogenesch_2A", "Hogenesch_2B", "Hogenesch_4A", "Hogenesch_4B", "Hogenesch_4C", "Hogenesch_4D", "Hughes_2", "Hughes_4A", "Hughes_4B", "Zhang_2", "Zhang_4A", "Zhang_4B")
p <- sapply(which(isSig), FUN = function(ds) {
  plotHeatmap(all.files[[ds]], datasetName = dataNames[ds])
})
dev.off()

# Cycle <- c("Per1", "Per2", "Per3", "Cry1", "Cry2", "Npas2", "Clock", "Arntl", "Arntl2")
# sigGenes[isSig]$Hogenesch_2B[which(sigGenes[isSig]$Hogenesch_2B %in% Cycle)]


# plotTsExpression <- function(geneName) {
#   par(mfrow = c(2, 2))
#   plot(x = seq(18, 64, 2), unname(unlist(timeSeries$Hogenesch2009.sub2A.data.ann_Results.Rdata[geneName, ])), pch = 21, bg = "orange", type = "o", main = paste0("Hogenesch_2A\n", geneName), xlab = "ZT Time", ylab = "Expression")
#   plot(x = seq(19, 65, 2), unname(unlist(timeSeries$Hogenesch2009.sub2B.data.ann_Results.Rdata[geneName, ])), pch = 21, bg = "orange", type = "o", main = paste0("Hogenesch_2B\n", geneName), xlab = "ZT Time", ylab = "Expression")
#   plot(x = seq(0, 46, 2), y = unname(unlist(timeSeries$Hughes2012.all.data.ann_Results.Rdata[geneName, ])), pch = 21, bg = "orange", type = "o", main = paste0("Hughes_2\n", geneName), xlab = "ZT Time", ylab = "Expression")
#   plot(x = seq(18, 64, 2), y = unname(unlist(timeSeries$Zhang2014.all.data.ann_Results.Rdata[geneName, ])), pch = 21, bg = "orange", type = "o", main = paste0("Zhang_2\n", geneName), xlab = "ZT Time", ylab = "Expression")
#   box(lwd = 3)
#   box(lwd = 1, col = "green")
# }
# 
# plotTsExpression("Cry1")
# plotTsExpression("Per1")
# plotTsExpression("Per3")
# plotTsExpression("Clock")
# plotTsExpression("Npas2")

# pairwise plots of Data with Correlation

load("~/Desktop/TimeCycleV3/Results/timeTrialAdjRealDataComplete.RData")
JTK_results <- resultsAdj[[3]]
genes <- rownames(JTK_results)
JTK_sigGene <- apply(JTK_results, 2, function(pVal) {
  genes[pVal < 0.05]
})
names(JTK_sigGene) <- dataNames
#plotUpset(df = JTK_sigGene[isSig])


RAIN_results <- resultsAdj[[4]]
genes <- rownames(RAIN_results)
RAIN_sigGene <- apply(RAIN_results, 2, function(pVal) {
  genes[p.adjust(pVal,"fdr") < 0.05]
})
names(RAIN_sigGene) <- dataNames
#plotUpset(df = RAIN_sigGene[isSig])


load("~/Desktop/TimeCycleV3/Results/timeTrialRealDataComplete.RData")


# Cycle <- c("Clock", "Cry1", "Npas2", "Per1", "Per3")
# data$Hogenesch2009.all.data.ann_Results.Rdata[Cycle, ]
# data$Hogenesch2009.sub2A.data.ann_Results.Rdata[Cycle, ]
# data$Hogenesch2009.sub2B.data.ann_Results.Rdata[Cycle, ]
# data$Hughes2012.all.data.ann_Results.Rdata[Cycle, ]
# data$Zhang2014.all.data.ann_Results.Rdata[Cycle, ]

# Get Overlaps between SigGene Lists
overlap <- intersect(sigGenes$Hughes_2, sigGenes$Zhang_2)
overlap_JTK <- intersect(JTK_sigGene$Hughes_2, JTK_sigGene$Zhang_2)
overlap_RAIN <- intersect(RAIN_sigGene$Hughes_2, RAIN_sigGene$Zhang_2)

dataForPlotting <- as.data.frame(data)
dataForPlotting <- dataForPlotting[complete.cases(dataForPlotting), ]

## Plot Amp, period, Phase for TimeCycle Overlap

par(mfcol = c(1, 3))
p1 <- ggscatterhist(
  data = dataForPlotting[overlap, ],
  x = "Hughes2012.all.data.ann_Results.Rdata.Amp",
  y = "Zhang2014.all.data.ann_Results.Rdata.Amp",
  xlab = "Hughes 2012",
  add = "reg.line", 
  add.params = list(color = "red", alpha = 0.5),
  conf.inf = T,
  alpha  = 0.5,
  cor.coef = T, 
  
  cor.method = "pearson",
  margin.plot = c("histogram"),
  ylab = "Zhang 2014", 
  asp = 1, 
  ggtheme = theme_bw(base_rect_size = 2)
) +
  ggtitle("Log(Amplitude)") +
  theme(plot.title = element_text(hjust = 0.5))

p2 <- ggscatterhist(
  data = dataForPlotting[overlap, ],
  x = "Hughes2012.all.data.ann_Results.Rdata.Period.in.Hours",
  y = "Zhang2014.all.data.ann_Results.Rdata.Period.in.Hours",
  xlab = "Hughes 2012",
  margin.plot = c("histogram"),
  alpha  = 0.5,
  ylab = "Zhang 2014", asp = 1, ggtheme = theme_bw(base_rect_size = 2)
) +
  ggtitle("Period (h)") +
  theme(plot.title = element_text(hjust = 0.5))


p3 <- ggscatterhist(
  data = dataForPlotting[overlap, ],
  x = "Hughes2012.all.data.ann_Results.Rdata.Phase.in.Hours",
  y = "Zhang2014.all.data.ann_Results.Rdata.Phase.in.Hours",
  xlab = "Hughes 2012",
  margin.plot = c("histogram"),
  alpha  = 0.5,
  ylab = "Zhang 2014", asp = 1, ggtheme = theme_bw(base_rect_size = 2)
) +
  ggtitle("Phase (h)") +
  theme(plot.title = element_text(hjust = 0.5)) 

pdf(file = "~/Desktop/TimeCycleV3/Results/Figures/TC_scatterPlots.pdf", width = 12, height = 4.3)
grid.arrange(grobs = list(p1,p2,p3), nrow = 1)
dev.off()


pdf(file = "~/Desktop/TimeCycleV3/Results/Figures/TC_corrPlots.pdf", width = 4, height = 4.5)
for(dataSet in c("Hogenesch_2A","Hogenesch_2B","Hughes_2", "Zhang_2") ){
overlap <- sigGenes[[dataSet]]
overlap_JTK <- JTK_sigGene[[dataSet]]
overlap_RAIN <- RAIN_sigGene[[dataSet]]

# Null Distribution - No ordering for gene sample. randomly drawn with replacement.
Null_Cor_Zhang_Hughes <- getCorrelation(controlData = sample_n(timeSeries$Zhang2014.all.data.ann_Results.Rdata, size = 12868, replace = T), comparisonData = sample_n(timeSeries$Hughes2012.all.data.ann_Results.Rdata, size = 12868, replace = T)[, c(22:24, 1:21)])
Null_Cor_Zhang_Hog_A <- getCorrelation(controlData = sample_n(timeSeries$Zhang2014.all.data.ann_Results.Rdata, size = 12868, replace = T), comparisonData = sample_n(timeSeries$Hogenesch2009.sub2A.data.ann_Results.Rdata, size = 12868, replace = T))
Null_Cor_Zhang_Hog_B <- getCorrelation(controlData = sample_n(timeSeries$Zhang2014.all.data.ann_Results.Rdata, size = 12868, replace = T), comparisonData = sample_n(timeSeries$Hogenesch2009.sub2B.data.ann_Results.Rdata, size = 12868, replace = T))
Null_Cor_Hughes_Hog_A <- getCorrelation(controlData = sample_n(timeSeries$Hughes2012.all.data.ann_Results.Rdata, size = 12868, replace = T)[, c(22:24, 1:21)], comparisonData = sample_n(timeSeries$Hogenesch2009.sub2A.data.ann_Results.Rdata, size = 12868, replace = T))
Null_Cor_Hughes_Hog_B <- getCorrelation(controlData = sample_n(timeSeries$Hughes2012.all.data.ann_Results.Rdata, size = 12868, replace = T)[, c(22:24, 1:21)], comparisonData = sample_n(timeSeries$Hogenesch2009.sub2B.data.ann_Results.Rdata, size = 12868, replace = T))
Null_Cor_Hog_A_Hog_B <- getCorrelation(controlData = sample_n(timeSeries$Hogenesch2009.sub2A.data.ann_Results.Rdata, size = 12868, replace = T), comparisonData = sample_n(timeSeries$Hogenesch2009.sub2B.data.ann_Results.Rdata, size = 12868, replace = T))

# Distribution given gene ordered samples, but no filter for sig genes
Ordered_Cor_Zhang_Hughes <- getCorrelation(controlData = timeSeries$Zhang2014.all.data.ann_Results.Rdata, comparisonData = timeSeries$Hughes2012.all.data.ann_Results.Rdata[, c(22:24, 1:21)])
Ordered_Cor_Zhang_Hog_A <- getCorrelation(controlData = timeSeries$Zhang2014.all.data.ann_Results.Rdata, comparisonData = timeSeries$Hogenesch2009.sub2A.data.ann_Results.Rdata)
Ordered_Cor_Zhang_Hog_B <- getCorrelation(controlData = timeSeries$Zhang2014.all.data.ann_Results.Rdata, comparisonData = timeSeries$Hogenesch2009.sub2B.data.ann_Results.Rdata)
Ordered_Cor_Hughes_Hog_A <- getCorrelation(controlData = timeSeries$Hughes2012.all.data.ann_Results.Rdata[, c(22:24, 1:21)], comparisonData = timeSeries$Hogenesch2009.sub2A.data.ann_Results.Rdata)
Ordered_Cor_Hughes_Hog_B <- getCorrelation(controlData = timeSeries$Hughes2012.all.data.ann_Results.Rdata[, c(22:24, 1:21)], comparisonData = timeSeries$Hogenesch2009.sub2B.data.ann_Results.Rdata)
Ordered_Cor_Hog_A_Hog_B <- getCorrelation(controlData = timeSeries$Hogenesch2009.sub2A.data.ann_Results.Rdata, comparisonData = timeSeries$Hogenesch2009.sub2B.data.ann_Results.Rdata)

# TimeCycle - Distribution given gene ordered samples, and filtered for sig genes
Sig_Cor_Zhang_Hughes <- getCorrelation(controlData = timeSeries$Zhang2014.all.data.ann_Results.Rdata[overlap, ], comparisonData = timeSeries$Hughes2012.all.data.ann_Results.Rdata[overlap, c(22:24, 1:21)])
Sig_Cor_Zhang_Hog_A <- getCorrelation(controlData = timeSeries$Zhang2014.all.data.ann_Results.Rdata[overlap, ], comparisonData = timeSeries$Hogenesch2009.sub2A.data.ann_Results.Rdata[overlap, ])
Sig_Cor_Zhang_Hog_B <- getCorrelation(controlData = timeSeries$Zhang2014.all.data.ann_Results.Rdata[overlap, ], comparisonData = timeSeries$Hogenesch2009.sub2B.data.ann_Results.Rdata[overlap, ])
Sig_Cor_Hughes_Hog_A <- getCorrelation(controlData = timeSeries$Hughes2012.all.data.ann_Results.Rdata[overlap, c(22:24, 1:21)], comparisonData = timeSeries$Hogenesch2009.sub2A.data.ann_Results.Rdata[overlap, ])
Sig_Cor_Hughes_Hog_B <- getCorrelation(controlData = timeSeries$Hughes2012.all.data.ann_Results.Rdata[overlap, c(22:24, 1:21)], comparisonData = timeSeries$Hogenesch2009.sub2B.data.ann_Results.Rdata[overlap, ])
Sig_Cor_Hog_A_Hog_B <- getCorrelation(controlData = timeSeries$Hogenesch2009.sub2A.data.ann_Results.Rdata[overlap, ], comparisonData = timeSeries$Hogenesch2009.sub2B.data.ann_Results.Rdata[overlap, ])

# RAIN - Distribution given gene ordered samples and filtered for sig genes
Sig_Cor_Zhang_Hughes_RAIN <- getCorrelation(controlData = timeSeries$Zhang2014.all.data.ann_Results.Rdata[overlap_RAIN, ], comparisonData = timeSeries$Hughes2012.all.data.ann_Results.Rdata[overlap_RAIN, c(22:24, 1:21)])
Sig_Cor_Zhang_Hog_A_RAIN <- getCorrelation(controlData = timeSeries$Zhang2014.all.data.ann_Results.Rdata[overlap_RAIN, ], comparisonData = timeSeries$Hogenesch2009.sub2A.data.ann_Results.Rdata[overlap_RAIN, ])
Sig_Cor_Zhang_Hog_B_RAIN <- getCorrelation(controlData = timeSeries$Zhang2014.all.data.ann_Results.Rdata[overlap_RAIN, ], comparisonData = timeSeries$Hogenesch2009.sub2B.data.ann_Results.Rdata[overlap_RAIN, ])
Sig_Cor_Hughes_Hog_A_RAIN <- getCorrelation(controlData = timeSeries$Hughes2012.all.data.ann_Results.Rdata[overlap_RAIN, c(22:24, 1:21)], comparisonData = timeSeries$Hogenesch2009.sub2A.data.ann_Results.Rdata[overlap_RAIN, ])
Sig_Cor_Hughes_Hog_B_RAIN <- getCorrelation(controlData = timeSeries$Hughes2012.all.data.ann_Results.Rdata[overlap_RAIN, c(22:24, 1:21)], comparisonData = timeSeries$Hogenesch2009.sub2B.data.ann_Results.Rdata[overlap_RAIN, ])
Sig_Cor_Hog_A_Hog_B_RAIN <- getCorrelation(controlData = timeSeries$Hogenesch2009.sub2A.data.ann_Results.Rdata[overlap_RAIN, ], comparisonData = timeSeries$Hogenesch2009.sub2B.data.ann_Results.Rdata[overlap_RAIN, ])

# JTK - Distribution given gene ordered samples and filtered for sig genes
Sig_Cor_Zhang_Hughes_JTK <- getCorrelation(controlData = timeSeries$Zhang2014.all.data.ann_Results.Rdata[overlap_JTK, ], comparisonData = timeSeries$Hughes2012.all.data.ann_Results.Rdata[overlap_JTK, c(22:24, 1:21)])
Sig_Cor_Zhang_Hog_A_JTK <- getCorrelation(controlData = timeSeries$Zhang2014.all.data.ann_Results.Rdata[overlap_JTK, ], comparisonData = timeSeries$Hogenesch2009.sub2A.data.ann_Results.Rdata[overlap_JTK, ])
Sig_Cor_Zhang_Hog_B_JTK <- getCorrelation(controlData = timeSeries$Zhang2014.all.data.ann_Results.Rdata[overlap_JTK, ], comparisonData = timeSeries$Hogenesch2009.sub2B.data.ann_Results.Rdata[overlap_JTK, ])
Sig_Cor_Hughes_Hog_A_JTK <- getCorrelation(controlData = timeSeries$Hughes2012.all.data.ann_Results.Rdata[overlap_JTK, c(22:24, 1:21)], comparisonData = timeSeries$Hogenesch2009.sub2A.data.ann_Results.Rdata[overlap_JTK, ])
Sig_Cor_Hughes_Hog_B_JTK <- getCorrelation(controlData = timeSeries$Hughes2012.all.data.ann_Results.Rdata[overlap_JTK, c(22:24, 1:21)], comparisonData = timeSeries$Hogenesch2009.sub2B.data.ann_Results.Rdata[overlap_JTK, ])
Sig_Cor_Hog_A_Hog_B_JTK <- getCorrelation(controlData = timeSeries$Hogenesch2009.sub2A.data.ann_Results.Rdata[overlap_JTK, ], comparisonData = timeSeries$Hogenesch2009.sub2B.data.ann_Results.Rdata[overlap_JTK, ])

NullDist <- c(Null_Cor_Zhang_Hughes, Null_Cor_Zhang_Hog_A, Null_Cor_Zhang_Hog_B, Null_Cor_Hughes_Hog_A, Null_Cor_Hughes_Hog_B, Null_Cor_Hog_A_Hog_B)
OrderedDist <- c(Ordered_Cor_Zhang_Hughes, Ordered_Cor_Zhang_Hog_A, Ordered_Cor_Zhang_Hog_B, Ordered_Cor_Hughes_Hog_A, Ordered_Cor_Hughes_Hog_B, Ordered_Cor_Hog_A_Hog_B)
SigDist <- c(Sig_Cor_Zhang_Hughes, Sig_Cor_Zhang_Hog_A, Sig_Cor_Zhang_Hog_B, Sig_Cor_Hughes_Hog_A, Sig_Cor_Hughes_Hog_B, Sig_Cor_Hog_A_Hog_B)
SigDist_RAIN <- c(Sig_Cor_Zhang_Hughes_RAIN, Sig_Cor_Zhang_Hog_A_RAIN, Sig_Cor_Zhang_Hog_B_RAIN, Sig_Cor_Hughes_Hog_A_RAIN, Sig_Cor_Hughes_Hog_B_RAIN, Sig_Cor_Hog_A_Hog_B_RAIN)
SigDist_JTK <- c(Sig_Cor_Zhang_Hughes_JTK, Sig_Cor_Zhang_Hog_A_JTK, Sig_Cor_Zhang_Hog_B_JTK, Sig_Cor_Hughes_Hog_A_JTK, Sig_Cor_Hughes_Hog_B_JTK, Sig_Cor_Hog_A_Hog_B_JTK)

par(mfrow = c(1, 1), font = 2)
lineWidth <- 6
plot(ecdf(SigDist_JTK), pch = NA ,xlim = range(c(SigDist, SigDist_JTK, SigDist_RAIN, OrderedDist, NullDist)), lwd = lineWidth, col = "#377eb8", main = paste0(dataSet," Rank Corr eCDF"), ylab = "eCDF", xlab = expression(paste("Rank Correlation (", rho, ")")))
plot(ecdf(SigDist), pch = NA, add = TRUE, lwd = lineWidth, col = "#ff7f00")
plot(ecdf(SigDist_RAIN), add = TRUE, lwd = lineWidth, col = "black")
plot(ecdf(OrderedDist), add = TRUE, lwd = lineWidth, col = "green3")
plot(ecdf(NullDist), add = TRUE, lwd = lineWidth, col = "grey")
legend(
  x = "topleft", inset = c(0.01, 0.05), legend = c("TimeCycle", "JTK_CYCLE", "RAIN", "All Gene Correlation", "Null Distribution"),
  col = c("#ff7f00", "#377eb8", "Black", "green3", "grey"), lty = 1, lwd = 4, cex = 0.75, box.lwd = 0
)
box(lwd = 2)
}
dev.off()

ks.test(x = SigDist, y = SigDist_JTK)
ks.test(x = SigDist_RAIN, y = SigDist)
ks.test(OrderedDist, SigDist_RAIN)

pdf("~/Desktop/blank2.pdf", width = 24, height = 16)
plot.new()
dev.off()
