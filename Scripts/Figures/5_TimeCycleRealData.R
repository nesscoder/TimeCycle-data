library(ggplot2)
library(UpSetR)
library(gplots)
library(tidyverse)
library(gridExtra)
library(ggpubr)

# packageurl <- "https://cran.r-project.org/src/contrib/Archive/ggpubr/ggpubr_0.3.0.tar.gz"
# install.packages(packageurl, repos=NULL, type="source")

require(AnnotationDbi)
require(org.Mm.eg.db)

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
plotUpset <- function(df, showPlot = T, onlyKnownCircadianGenes = F) {
  methodNames <- c("TimeCycle")
  # identify the Cycling Genes in the Dataset
  Cycle <- getKnownCircadainGenes()

  Genes <- unique(unlist(df))
  collected <- as.numeric(Genes %in% Cycle)
  print(length(Cycle %in% unique(Genes)))
  Cyclers <- c(0, 2)[collected + 1]

  # Bind Data together to Use for Plotting
  data <- cbind(fromList(df), Cyclers)
  if (onlyKnownCircadianGenes) {
    data <- data[which(data$Cyclers == 2), ]
  }

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
    group.by = "degree"
  )
}

##########################
# Description:
#   Take in the CGBD database and filers for Mus Genes in the Liver that have been experimentally Validated
# Output: Returns list of genes
##########################

getKnownCircadainGenes <- function() {
  cgbd_database <- read.delim("~/Desktop/TimeCycle-data/Data/RealData/cgdb.txt", comment.char = "#")

  # only select genes that are from the mouse liver and validated by small-scale experiments
  df <- cgbd_database %>%
    filter(Organism == "Mus musculus (Mouse)") %>%
    filter(str_detect(str_to_lower(Tissue.Cell), "\\liver") & !str_detect(str_to_lower(Tissue.Cell), "\\liver, heart, kidney")) %>%
    filter(str_detect(str_to_lower(Data.type), "\\Small-scale"))

  # convert genes to Gene Names
  geneSymbols <- mapIds(org.Mm.eg.db, keys = unique(df$Uniprot.Ensembl.ID), column = "SYMBOL", keytype = "UNIPROT", multiVals = "first")
  inds <- which(!is.na(geneSymbols))
  found_genes <- geneSymbols[inds]
  return(unname(found_genes))
}

##########################
#Input: Data set of pvals for a given method and dataset name (Hogenesch_2a, Hogenesch_2B, Hughes_2, Zhang_2)
# Description:
#   Creates probability density histogram of the known circadian genes vs all genes
# Output: Histogram plot
##########################

plotKnownCircadainHist <- function(data, dataset, KS = F) {
  breaks <- seq(0, 1, length.out = 100)
  data1 <- data[, dataset]
  data2 <- data[intersect(rownames(data), getKnownCircadainGenes()), dataset]

  h1 <- hist(data1, breaks = breaks, plot = F)
  h2 <- hist(data2, breaks = breaks, plot = F)

  h1$counts <- h1$counts / sum(h1$counts)
  h2$counts <- h2$counts / sum(h2$counts)
  if(KS){
  ksTestResult <- ks.test(x = h1$counts, y = h2$counts)
  }

  ylim <- range(c(h1$counts, h2$counts))
  plot(h2, main = gsub(pattern = "_", replacement = " ", dataset), xlab = "FDR", ylab = "Density", col = "#ff7f00", xlim = c(0, 1), ylim = ylim, lty="blank", axes = F, xaxs = "i", yaxs="i")
  plot(h1, col = ggplot2::alpha("grey45", 0.8), add = T, lty="blank")
  if(KS){
  adjPval <- ksTestResult$p.value*4
  text(x = 0.5, y = ylim[2] / 2, paste0("KS Test: p = ", signif(adjPval,1)), font = 2)
  }
  axis(side=1, lwd=2, xpd=TRUE)
  axis(side=2, lwd=2, xpd=TRUE)
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

  # creates a own color palette from lightgrey to orange
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

# get the adjusted pvalues and make them into a single dataframe
adjustpVal <- lapply(data, function(df) {
  df$pVals.adj
})
adjustpVal <- do.call(cbind, adjustpVal)
rownames(adjustpVal) <- rownames(data[[1]])

## --------------------------------------- Plot the data ---------------------------------------

########################################################
## TimeCycle Upset Plots
########################################################
sigGenes <- getSigGeneList(df = adjustpVal, pval = 0.05)
sigGenes <- lapply(sigGenes, na.omit) # remove NA genes, see getSigGeneList Function
# only plot list if it is significant
isSig <- as.vector(unlist(lapply(sigGenes, function(x) {
  length(x) > 1
})))
isSig[1] <- F

dataNames <- c("Hogenesch_1", "Hogenesch_2A", "Hogenesch_2B", "Hogenesch_4A", "Hogenesch_4B", "Hogenesch_4C", "Hogenesch_4D", "Hughes_2", "Hughes_4A", "Hughes_4B", "Zhang_2", "Zhang_4A", "Zhang_4B")
colnames(adjustpVal) <- dataNames

pdf(file = "~/Desktop/TimeCycle-data/Results/Figures/upsetPlot.pdf", width = 12, height = 8)
plotUpset(sigGenes[isSig])
dev.off()

pdf(file = "~/Desktop/TimeCycle-data/Results/Figures/hist_knownCyclers.pdf", width = 12, height = 8)
par(mfrow = c(2, 2), font.lab = 2)
plotKnownCircadainHist(data = adjustpVal, dataset = "Hogenesch_2A", KS = T)
plotKnownCircadainHist(data = adjustpVal, dataset = "Hogenesch_2B", KS = T)
plotKnownCircadainHist(data = adjustpVal, dataset = "Hughes_2", KS = T)
plotKnownCircadainHist(data = adjustpVal, dataset = "Zhang_2", KS = T)
dev.off()

########################################################
# create plots of cycling genes in each Dataset
########################################################
pdf(file = "~/Desktop/TimeCycle-data/Results/Figures/heatMaps.pdf", width = 4, height = 4)
p <- sapply(which(isSig), FUN = function(ds) {
  plotHeatmap(all.files[[ds]], datasetName = dataNames[ds])
})
dev.off()

########################################################
## JTK Upset Plots
########################################################
load("~/Desktop/TimeCycle-data/Results/timeTrialAdjRealDataComplete.RData")
JTK_results <- resultsAdj[[3]]
genes <- rownames(JTK_results)
JTK_sigGene <- apply(JTK_results, 2, function(pVal) {
  genes[pVal < 0.05]
})
names(JTK_sigGene) <- dataNames
colnames(JTK_results) <- dataNames
  
pdf(file = "~/Desktop/TimeCycle-data/Results/Figures/supplement_upsetPlot_JTK.pdf", width = 12, height = 8)
plotUpset(df = JTK_sigGene[isSig])
dev.off()


pdf(file = "~/Desktop/TimeCycle-data/Results/Figures/hist_knownCyclers_JTK.pdf", width = 12, height = 8)
par(mfrow = c(2, 2), font.lab = 2)
plotKnownCircadainHist(data = JTK_results, dataset = "Hogenesch_2A")
plotKnownCircadainHist(data = JTK_results, dataset = "Hogenesch_2B")
plotKnownCircadainHist(data = JTK_results, dataset = "Hughes_2")
plotKnownCircadainHist(data = JTK_results, dataset = "Zhang_2")
dev.off()

########################################################
## RAIN Upset Plots
########################################################
RAIN_results <- resultsAdj[[4]]
genes <- rownames(RAIN_results)
RAIN_sigGene <- apply(RAIN_results, 2, function(pVal) {
  genes[p.adjust(pVal, "fdr") < 0.05]
})
names(RAIN_sigGene) <- dataNames

RAIN_results <- apply(RAIN_results, 2, function(pVal) {
  p.adjust(pVal, "fdr")
})
colnames(RAIN_results) <- dataNames

pdf(file = "~/Desktop/TimeCycle-data/Results/Figures/supplement_upsetPlot_RAIN.pdf", width = 12, height = 8)
plotUpset(df = RAIN_sigGene[isSig])
dev.off()

pdf(file = "~/Desktop/TimeCycle-data/Results/Figures/hist_knownCyclers_RAIN.pdf", width = 12, height = 8)
par(mfrow = c(2, 2), font.lab = 2)
plotKnownCircadainHist(data = RAIN_results, dataset = "Hogenesch_2A")
plotKnownCircadainHist(data = RAIN_results, dataset = "Hogenesch_2B")
plotKnownCircadainHist(data = RAIN_results, dataset = "Hughes_2")
plotKnownCircadainHist(data = RAIN_results, dataset = "Zhang_2")
dev.off()

##########################################
## GeneCycle Upset Plots
##########################################
genes <- rownames(GeneCycle_results)
GeneCycle_sigGene <- apply(GeneCycle_results, 2, function(pVal) {
  genes[pVal < 0.05]
})
names(GeneCycle_sigGene) <- dataNames
colnames(GeneCycle_results) <- dataNames

pdf(file = "~/Desktop/TimeCycle-data/Results/Figures/supplement_upsetPlot_GeneCycle.pdf", width = 12, height = 8)
plotUpset(df = GeneCycle_sigGene[isSig])
dev.off()

pdf(file = "~/Desktop/TimeCycle-data/Results/Figures/hist_knownCyclers_GeneCycle.pdf", width = 12, height = 8)
par(mfrow = c(2, 2), font.lab = 2)
plotKnownCircadainHist(data = GeneCycle_results, dataset = "Hogenesch_2A")
plotKnownCircadainHist(data = GeneCycle_results, dataset = "Hogenesch_2B")
plotKnownCircadainHist(data = GeneCycle_results, dataset = "Hughes_2")
plotKnownCircadainHist(data = GeneCycle_results, dataset = "Zhang_2")
dev.off()


########################################################
# Get Overlaps between SigGene Lists
########################################################
load("~/Desktop/TimeCycle-data/Results/timeTrialRealDataComplete.RData")
overlap <- intersect(sigGenes$Hughes_2, sigGenes$Zhang_2)
overlap_JTK <- intersect(JTK_sigGene$Hughes_2, JTK_sigGene$Zhang_2)
overlap_RAIN <- intersect(RAIN_sigGene$Hughes_2, RAIN_sigGene$Zhang_2)
overlap_GeneCycle <- intersect(GeneCycle_sigGene$Hughes_2, GeneCycle_sigGene$Zhang_2)

dataForPlotting <- as.data.frame(data)
dataForPlotting <- dataForPlotting[complete.cases(dataForPlotting), ]

######################################################## 
## Plot Amp, period, Phase for TimeCycle Overlap
########################################################
par(mfcol = c(1, 3))
p1 <- ggscatterhist(
  data = dataForPlotting[overlap, ],
  x = "Hughes2012.all.data.ann_Results.Rdata.Amp",
  y = "Zhang2014.all.data.ann_Results.Rdata.Amp",
  xlab = "Hughes 2012",
  add = "reg.line",
  add.params = list(color = "red", alpha = 0.5),
  conf.inf = T,
  alpha = 0.5,
  cor.coef = T,
  title = "LogFC Amplitude",
  cor.method = "pearson",
  margin.plot = c("histogram"),
  ylab = "Zhang 2014",
  asp = 1,
  ggtheme = theme_bw(base_rect_size = 2) + theme(plot.title = element_text(hjust = 0.5))
)

p2 <- ggscatterhist(
  data = dataForPlotting[overlap, ],
  x = "Hughes2012.all.data.ann_Results.Rdata.Period.in.Hours",
  y = "Zhang2014.all.data.ann_Results.Rdata.Period.in.Hours",
  xlab = "Hughes 2012",
  margin.plot = c("histogram"),
  alpha = 0.5,
  title = "Period (h)",
  ylab = "Zhang 2014",
  asp = 1,
  ggtheme = theme_bw(base_rect_size = 2) + theme(plot.title = element_text(hjust = 0.5))
)

p3 <- ggscatterhist(
  data = dataForPlotting[overlap, ],
  x = "Hughes2012.all.data.ann_Results.Rdata.Phase.in.Hours",
  y = "Zhang2014.all.data.ann_Results.Rdata.Phase.in.Hours",
  xlab = "Hughes 2012",
  margin.plot = c("histogram"),
  alpha = 0.5,
  title = "Phase (h)",
  ylab = "Zhang 2014", asp = 1, ggtheme = theme_bw(base_rect_size = 2) + theme(plot.title = element_text(hjust = 0.5))
)

pdf(file = "~/Desktop/TimeCycle-data/Results/Figures/TC_scatterPlots.pdf", width = 12, height = 4.3)
grid.arrange(grobs = list(p1, p2, p3), nrow = 1)
dev.off()




pdf(file = "~/Desktop/TimeCycle-data/Results/Figures/TC_corrPlots.pdf", width = 4, height = 4.5)
names(timeSeries) <- dataNames
for (dataSet in c("Hogenesch_2A", "Hogenesch_2B", "Hughes_2", "Zhang_2")) {
  amplitude <- apply(timeSeries[[dataSet]], 1, function(i) {
    ampChange <- range(i)
    ampChange[2] - ampChange[1]
  })
  overlap <- intersect(names(which(amplitude > 2)), sigGenes[[dataSet]])
  overlap_JTK <- intersect(names(which(amplitude > 2)), JTK_sigGene[[dataSet]])
  overlap_RAIN <- intersect(names(which(amplitude > 2)), RAIN_sigGene[[dataSet]])
  overlap_GeneCycle <- intersect(names(which(amplitude > 2)), GeneCycle_sigGene[[dataSet]])

  # Null Distribution - No ordering for gene sample. randomly drawn with replacement.
  Null_Cor_Zhang_Hughes <- getCorrelation(controlData = sample_n(timeSeries$Zhang_2, size = 12868, replace = T), comparisonData = sample_n(timeSeries$Hughes_2, size = 12868, replace = T)[, c(22:24, 1:21)])
  Null_Cor_Zhang_Hog_A <- getCorrelation(controlData = sample_n(timeSeries$Zhang_2, size = 12868, replace = T), comparisonData = sample_n(timeSeries$Hogenesch_2A, size = 12868, replace = T))
  Null_Cor_Zhang_Hog_B <- getCorrelation(controlData = sample_n(timeSeries$Zhang_2, size = 12868, replace = T), comparisonData = sample_n(timeSeries$Hogenesch_2B, size = 12868, replace = T))
  Null_Cor_Hughes_Hog_A <- getCorrelation(controlData = sample_n(timeSeries$Hughes_2, size = 12868, replace = T)[, c(22:24, 1:21)], comparisonData = sample_n(timeSeries$Hogenesch_2A, size = 12868, replace = T))
  Null_Cor_Hughes_Hog_B <- getCorrelation(controlData = sample_n(timeSeries$Hughes_2, size = 12868, replace = T)[, c(22:24, 1:21)], comparisonData = sample_n(timeSeries$Hogenesch_2B, size = 12868, replace = T))
  Null_Cor_Hog_A_Hog_B <- getCorrelation(controlData = sample_n(timeSeries$Hogenesch_2A, size = 12868, replace = T), comparisonData = sample_n(timeSeries$Hogenesch_2B, size = 12868, replace = T))

  # Distribution given gene ordered samples, but no filter for sig genes
  Ordered_Cor_Zhang_Hughes <- getCorrelation(controlData = timeSeries$Zhang_2, comparisonData = timeSeries$Hughes_2[, c(22:24, 1:21)])
  Ordered_Cor_Zhang_Hog_A <- getCorrelation(controlData = timeSeries$Zhang_2, comparisonData = timeSeries$Hogenesch_2A)
  Ordered_Cor_Zhang_Hog_B <- getCorrelation(controlData = timeSeries$Zhang_2, comparisonData = timeSeries$Hogenesch_2B)
  Ordered_Cor_Hughes_Hog_A <- getCorrelation(controlData = timeSeries$Hughes_2[, c(22:24, 1:21)], comparisonData = timeSeries$Hogenesch_2A)
  Ordered_Cor_Hughes_Hog_B <- getCorrelation(controlData = timeSeries$Hughes_2[, c(22:24, 1:21)], comparisonData = timeSeries$Hogenesch_2B)
  Ordered_Cor_Hog_A_Hog_B <- getCorrelation(controlData = timeSeries$Hogenesch_2A, comparisonData = timeSeries$Hogenesch_2B)

  # TimeCycle - Distribution given gene ordered samples, and filtered for sig genes
  Sig_Cor_Zhang_Hughes <- getCorrelation(controlData = timeSeries$Zhang_2[overlap, ], comparisonData = timeSeries$Hughes_2[overlap, c(22:24, 1:21)])
  Sig_Cor_Zhang_Hog_A <- getCorrelation(controlData = timeSeries$Zhang_2[overlap, ], comparisonData = timeSeries$Hogenesch_2A[overlap, ])
  Sig_Cor_Zhang_Hog_B <- getCorrelation(controlData = timeSeries$Zhang_2[overlap, ], comparisonData = timeSeries$Hogenesch_2B[overlap, ])
  Sig_Cor_Hughes_Hog_A <- getCorrelation(controlData = timeSeries$Hughes_2[overlap, c(22:24, 1:21)], comparisonData = timeSeries$Hogenesch_2A[overlap, ])
  Sig_Cor_Hughes_Hog_B <- getCorrelation(controlData = timeSeries$Hughes_2[overlap, c(22:24, 1:21)], comparisonData = timeSeries$Hogenesch_2B[overlap, ])
  Sig_Cor_Hog_A_Hog_B <- getCorrelation(controlData = timeSeries$Hogenesch_2A[overlap, ], comparisonData = timeSeries$Hogenesch_2B[overlap, ])

  # RAIN - Distribution given gene ordered samples and filtered for sig genes
  Sig_Cor_Zhang_Hughes_RAIN <- getCorrelation(controlData = timeSeries$Zhang_2[overlap_RAIN, ], comparisonData = timeSeries$Hughes_2[overlap_RAIN, c(22:24, 1:21)])
  Sig_Cor_Zhang_Hog_A_RAIN <- getCorrelation(controlData = timeSeries$Zhang_2[overlap_RAIN, ], comparisonData = timeSeries$Hogenesch_2A[overlap_RAIN, ])
  Sig_Cor_Zhang_Hog_B_RAIN <- getCorrelation(controlData = timeSeries$Zhang_2[overlap_RAIN, ], comparisonData = timeSeries$Hogenesch_2B[overlap_RAIN, ])
  Sig_Cor_Hughes_Hog_A_RAIN <- getCorrelation(controlData = timeSeries$Hughes_2[overlap_RAIN, c(22:24, 1:21)], comparisonData = timeSeries$Hogenesch_2A[overlap_RAIN, ])
  Sig_Cor_Hughes_Hog_B_RAIN <- getCorrelation(controlData = timeSeries$Hughes_2[overlap_RAIN, c(22:24, 1:21)], comparisonData = timeSeries$Hogenesch_2B[overlap_RAIN, ])
  Sig_Cor_Hog_A_Hog_B_RAIN <- getCorrelation(controlData = timeSeries$Hogenesch_2A[overlap_RAIN, ], comparisonData = timeSeries$Hogenesch_2B[overlap_RAIN, ])

  # JTK - Distribution given gene ordered samples and filtered for sig genes
  Sig_Cor_Zhang_Hughes_JTK <- getCorrelation(controlData = timeSeries$Zhang_2[overlap_JTK, ], comparisonData = timeSeries$Hughes_2[overlap_JTK, c(22:24, 1:21)])
  Sig_Cor_Zhang_Hog_A_JTK <- getCorrelation(controlData = timeSeries$Zhang_2[overlap_JTK, ], comparisonData = timeSeries$Hogenesch_2A[overlap_JTK, ])
  Sig_Cor_Zhang_Hog_B_JTK <- getCorrelation(controlData = timeSeries$Zhang_2[overlap_JTK, ], comparisonData = timeSeries$Hogenesch_2B[overlap_JTK, ])
  Sig_Cor_Hughes_Hog_A_JTK <- getCorrelation(controlData = timeSeries$Hughes_2[overlap_JTK, c(22:24, 1:21)], comparisonData = timeSeries$Hogenesch_2A[overlap_JTK, ])
  Sig_Cor_Hughes_Hog_B_JTK <- getCorrelation(controlData = timeSeries$Hughes_2[overlap_JTK, c(22:24, 1:21)], comparisonData = timeSeries$Hogenesch_2B[overlap_JTK, ])
  Sig_Cor_Hog_A_Hog_B_JTK <- getCorrelation(controlData = timeSeries$Hogenesch_2A[overlap_JTK, ], comparisonData = timeSeries$Hogenesch_2B[overlap_JTK, ])

  # Genecycle - Distribution given gene ordered samples and filtered for sig genes
  Sig_Cor_Zhang_Hughes_GeneCycle <- getCorrelation(controlData = timeSeries$Zhang_2[overlap_GeneCycle, ], comparisonData = timeSeries$Hughes_2[overlap_GeneCycle, c(22:24, 1:21)])
  Sig_Cor_Zhang_Hog_A_GeneCycle <- getCorrelation(controlData = timeSeries$Zhang_2[overlap_GeneCycle, ], comparisonData = timeSeries$Hogenesch_2A[overlap_GeneCycle, ])
  Sig_Cor_Zhang_Hog_B_GeneCycle <- getCorrelation(controlData = timeSeries$Zhang_2[overlap_GeneCycle, ], comparisonData = timeSeries$Hogenesch_2B[overlap_GeneCycle, ])
  Sig_Cor_Hughes_Hog_A_GeneCycle <- getCorrelation(controlData = timeSeries$Hughes_2[overlap_GeneCycle, c(22:24, 1:21)], comparisonData = timeSeries$Hogenesch_2A[overlap_GeneCycle, ])
  Sig_Cor_Hughes_Hog_B_GeneCycle <- getCorrelation(controlData = timeSeries$Hughes_2[overlap_GeneCycle, c(22:24, 1:21)], comparisonData = timeSeries$Hogenesch_2B[overlap_GeneCycle, ])
  Sig_Cor_Hog_A_Hog_B_GeneCycle <- getCorrelation(controlData = timeSeries$Hogenesch_2A[overlap_GeneCycle, ], comparisonData = timeSeries$Hogenesch_2B[overlap_GeneCycle, ])


  NullDist <- c(Null_Cor_Zhang_Hughes, Null_Cor_Zhang_Hog_A, Null_Cor_Zhang_Hog_B, Null_Cor_Hughes_Hog_A, Null_Cor_Hughes_Hog_B, Null_Cor_Hog_A_Hog_B)
  OrderedDist <- c(Ordered_Cor_Zhang_Hughes, Ordered_Cor_Zhang_Hog_A, Ordered_Cor_Zhang_Hog_B, Ordered_Cor_Hughes_Hog_A, Ordered_Cor_Hughes_Hog_B, Ordered_Cor_Hog_A_Hog_B)
  SigDist <- c(Sig_Cor_Zhang_Hughes, Sig_Cor_Zhang_Hog_A, Sig_Cor_Zhang_Hog_B, Sig_Cor_Hughes_Hog_A, Sig_Cor_Hughes_Hog_B, Sig_Cor_Hog_A_Hog_B)
  SigDist_RAIN <- c(Sig_Cor_Zhang_Hughes_RAIN, Sig_Cor_Zhang_Hog_A_RAIN, Sig_Cor_Zhang_Hog_B_RAIN, Sig_Cor_Hughes_Hog_A_RAIN, Sig_Cor_Hughes_Hog_B_RAIN, Sig_Cor_Hog_A_Hog_B_RAIN)
  SigDist_JTK <- c(Sig_Cor_Zhang_Hughes_JTK, Sig_Cor_Zhang_Hog_A_JTK, Sig_Cor_Zhang_Hog_B_JTK, Sig_Cor_Hughes_Hog_A_JTK, Sig_Cor_Hughes_Hog_B_JTK, Sig_Cor_Hog_A_Hog_B_JTK)
  SigDist_GeneCycle <- c(Sig_Cor_Zhang_Hughes_GeneCycle, Sig_Cor_Zhang_Hog_A_GeneCycle, Sig_Cor_Zhang_Hog_B_GeneCycle, Sig_Cor_Hughes_Hog_A_GeneCycle, Sig_Cor_Hughes_Hog_B_GeneCycle, Sig_Cor_Hog_A_Hog_B_GeneCycle)

  par(mfrow = c(1, 1), font = 2)
  lineWidth <- 3
  plot(ecdf(SigDist_JTK), pch = NA, xlim = range(c(SigDist, SigDist_JTK, SigDist_RAIN, OrderedDist, NullDist)), lwd = lineWidth, col = "#377eb8", main = paste0(dataSet, " Rank Corr eCDF"), ylab = "eCDF", xlab = expression(paste("Rank Correlation (", rho, ")")))
  plot(ecdf(SigDist_RAIN), pch = NA, add = TRUE, lwd = lineWidth, col = "black")
  plot(ecdf(SigDist_GeneCycle), pch = NA, add = TRUE, lwd = lineWidth, col = "green3")
  plot(ecdf(SigDist), pch = NA, add = TRUE, lwd = lineWidth, col = "#ff7f00")
  plot(ecdf(OrderedDist), add = TRUE, lwd = lineWidth, col = "grey50")
  plot(ecdf(NullDist), add = TRUE, lwd = lineWidth, col = "grey")
  legend(
    x = "topleft", inset = c(0.01, 0.05), legend = c("TimeCycle", "JTK_CYCLE", "RAIN", "GeneCycle", "All Gene Correlation", "Null Distribution"),
    col = c("#ff7f00", "#377eb8", "Black", "green3", "grey50", "grey"), lty = 1, lwd = 4, cex = 0.75, box.lwd = 0
  )
  box(lwd = 2)

  print(ks.test(x = SigDist, y = SigDist_RAIN))
  print(ks.test(x = SigDist_JTK, y = SigDist_RAIN))
  print(ks.test(x = SigDist_GeneCycle, y = SigDist_RAIN))

  print(ks.test(x = SigDist, y = NullDist))
  print(ks.test(x = SigDist_JTK, y = NullDist))
  print(ks.test(x = SigDist_GeneCycle, y = NullDist))

  print(ks.test(x = SigDist, y = OrderedDist))
  print(ks.test(x = SigDist_JTK, y = OrderedDist))
  print(ks.test(x = SigDist_GeneCycle, y = OrderedDist))
}
dev.off()
