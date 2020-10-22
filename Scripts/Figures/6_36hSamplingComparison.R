##########################################
# Plot Results
##########################################
load("~/Desktop/TimeCycle-data/Results/SlidingWindows/36HourComparison.Rdata")
library(UpSetR)
library(tidyverse)
library(ggalluvial)
library(gridExtra)

# Check to see Correlation between genes
getCorrelation <- function(controlData, comparisonData) {
  output <- sapply(1:nrow(controlData), function(i) {
    set1 <- unlist(unname(controlData[i, ]))
    set2 <- unlist(unname(comparisonData[i, ]))
    cor(set1, set2, method = "spearman")
  })
  return(output)
}

getSigGenesCount <- function(data, alpha1, alpha2) {
  geneNames <- rownames(data)
  sigGenes_0.05 <- colSums(t(data[, START:END] < alpha1))
  sigGenes_0.1 <- colSums(t(data[, START:END] < alpha2))

  df_sigGenes_0.05 <- data.frame(gene = geneNames, overlaps = sigGenes_0.05, threshold = "alpha1")
  df_sigGenes_0.1 <- data.frame(gene = geneNames, overlaps = sigGenes_0.1, threshold = "alpha2")

  return(rbind(df_sigGenes_0.05, df_sigGenes_0.1))
}

#############################
# TC Results
#############################
START <- 3
END <- 6

TC_sigGenes <- getSigGenesCount(data = TC_results, alpha1 = 0.05, alpha2 = 0.1)
JTK_sigGenes <- getSigGenesCount(data = JTK_results, alpha1 = 0.05, alpha2 = 0.1)
RAIN_sigGenes <- getSigGenesCount(data = RAIN_results, alpha1 = 0.05, alpha2 = 0.1)



plotAlluvial <- function(data, plotTitle, col) {
  theme_set(new = theme_light())
  theme_replace(
    panel.border = element_rect(colour = "black", fill = NA, size = 2),
    text = element_text(size = 16, face = "bold"))
  
  data$overlaps <- as.factor(as.character(data$overlaps))
  data$gene <- as.factor(data$gene)
  data$threshold <- as.factor(data$threshold)
  
  ggplot(data = data, mapping = aes(
           x = threshold, 
           stratum = overlaps, 
           alluvium = gene,
           fill = overlaps, 
           label = overlaps)) +
    # scale_fill_brewer(type = "qual", palette = "Set1") +
    scale_fill_manual(values = col) +
    geom_flow() +
    geom_stratum() +
    scale_x_discrete(labels = c("FDR < 0.05", "FDR < 0.1")) +
    theme(legend.position = "bottom", panel.grid.major.x = element_blank()) +
    guides(fill = guide_legend(override.aes = list(linetype = 0, shape = 22, size = 5), title = "Number of Overlaps", title.position = "top", ncol = 5, byrow = TRUE)) +
    xlab(label = "") +
    ggrepel::geom_text_repel(data = subset(data, subset = as.numeric(threshold) == 1), stat = "stratum", aes(label = 100 * round(after_stat(prop), 3)),
                             size = 4, direction = "y", nudge_x = -0.5) +
    ggrepel::geom_text_repel(data = subset(data, subset = as.numeric(threshold) == 2), stat = "stratum", aes(label = 100 * round(after_stat(prop), 3)),
                             size = 4, direction = "y", nudge_x = 0.5) +
    ggtitle(plotTitle)
}
# customColors <- c("#7F3C8D","#11A579","#3969AC","#F2B701","#E73F74","#80BA5A","#E68310","#008695","#CF1C90","#f97b72","#4b4b8f","#A5AA99")
customColors <- c("#cccccc", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#FFFF33", "#a65628", "#cccccc", "#969696", "#636363", "#252525")

pdf("~/Desktop/TimeCycle-data/Results/Figures/AluvOut36.pdf", width = 14, height = 5)
p1 <- plotAlluvial(data = TC_sigGenes, plotTitle = "TimeCycle", col = customColors[1:5])
p2 <- plotAlluvial(data = JTK_sigGenes, plotTitle = "JTK_CYCLE", col = customColors[1:5])
p3 <- plotAlluvial(data = RAIN_sigGenes, plotTitle = "RAIN", col = customColors[1:5])
grid.arrange(grobs = list(p1, p2, p3), nrow = 1)
dev.off()

# Correlation Across Data Sets
alpha <- 0.1
geneNames <- rownames(TC_results)
tc_sg <- apply(TC_results_Dataset_Comp, 2, function(pVals) {
  geneNames[pVals < alpha]
})
jtk_sg <- apply(JTK_results_Dataset_Comp, 2, function(pVals) {
  geneNames[pVals < alpha]
})
rain_sg <- apply(RAIN_results_Dataset_Comp, 2, function(pVals) {
  geneNames[p.adjust(pVals, "fdr") < alpha]
})


pdf("~/Desktop/TimeCycle-data/Results/Figures/36hCorr.pdf", width = 14, height = 5)
par(mfrow = c(1, 3), font = 2)
for (ds in c("Hogenesch", "Hughes", "Zhang")) {
  overlap_JTK <- jtk_sg[[ds]]
  overlap_RAIN <- rain_sg[[ds]]
  overlap_TC <- tc_sg[[ds]]

  # Null Distribution - No ordering for gene sample. randomly drawn with replacement.
  Null_Cor_Zhang_Hughes <- getCorrelation(controlData = sample_n(Zhang, size = 12868, replace = T), comparisonData = sample_n(Hughes, size = 12868, replace = T))
  Null_Cor_Zhang_Hogenesch <- getCorrelation(controlData = sample_n(Zhang, size = 12868, replace = T), comparisonData = sample_n(Hogenesch, size = 12868, replace = T))
  Null_Cor_Hughes_Hogenesch <- getCorrelation(controlData = sample_n(Hughes, size = 12868, replace = T), comparisonData = sample_n(Hogenesch, size = 12868, replace = T))

  # Distribution given gene ordered samples, but no filter for sig genes
  Ordered_Cor_Zhang_Hughes <- getCorrelation(controlData = Zhang, comparisonData = Hughes)
  Ordered_Cor_Zhang_Hogenesch <- getCorrelation(controlData = Zhang, comparisonData = Hogenesch)
  Ordered_Cor_Hughes_Hogenesch <- getCorrelation(controlData = Hughes, comparisonData = Hogenesch)

  # RAIN - Distribution given gene ordered samples and filtered for sig genes
  Sig_Cor_Zhang_Hughes_TC <- getCorrelation(controlData = Zhang[overlap_TC, ], comparisonData = Hughes[overlap_TC, ])
  Sig_Cor_Zhang_Hogenesch_TC <- getCorrelation(controlData = Zhang[overlap_TC, ], comparisonData = Hogenesch[overlap_TC, ])
  Sig_Cor_Hughes_Hogenesch_TC <- getCorrelation(controlData = Hughes[overlap_TC, ], comparisonData = Hogenesch[overlap_TC, ])

  # RAIN - Distribution given gene ordered samples and filtered for sig genes
  Sig_Cor_Zhang_Hughes_RAIN <- getCorrelation(controlData = Zhang[overlap_RAIN, ], comparisonData = Hughes[overlap_RAIN, ])
  Sig_Cor_Zhang_Hogenesch_RAIN <- getCorrelation(controlData = Zhang[overlap_RAIN, ], comparisonData = Hogenesch[overlap_RAIN, ])
  Sig_Cor_Hughes_Hogenesch_RAIN <- getCorrelation(controlData = Hughes[overlap_RAIN, ], comparisonData = Hogenesch[overlap_RAIN, ])

  # JTK - Distribution given gene ordered samples and filtered for sig genes
  Sig_Cor_Zhang_Hughes_JTK <- getCorrelation(controlData = Zhang[overlap_JTK, ], comparisonData = Hughes[overlap_JTK, ])
  Sig_Cor_Zhang_Hogenesch_JTK <- getCorrelation(controlData = Zhang[overlap_JTK, ], comparisonData = Hogenesch[overlap_JTK, ])
  Sig_Cor_Hughes_Hogenesch_JTK <- getCorrelation(controlData = Hughes[overlap_JTK, ], comparisonData = Hogenesch[overlap_JTK, ])

  # merge correlations into single vector
  NullDist <- c(Null_Cor_Zhang_Hughes, Null_Cor_Zhang_Hogenesch, Null_Cor_Hughes_Hogenesch)
  OrderedDist <- c(Ordered_Cor_Zhang_Hughes, Ordered_Cor_Zhang_Hogenesch, Ordered_Cor_Hughes_Hogenesch)
  SigDist_TC <- c(Sig_Cor_Zhang_Hughes_TC, Sig_Cor_Zhang_Hogenesch_TC, Sig_Cor_Hughes_Hogenesch_TC)
  SigDist_RAIN <- c(Sig_Cor_Zhang_Hughes_RAIN, Sig_Cor_Zhang_Hogenesch_RAIN, Sig_Cor_Hughes_Hogenesch_RAIN)
  SigDist_JTK <- c(Sig_Cor_Zhang_Hughes_JTK, Sig_Cor_Zhang_Hogenesch_JTK, Sig_Cor_Hughes_Hogenesch_JTK)

  # par(mfrow = c(1, 5))
  # hist(NullDist)
  # hist(OrderedDist)
  # hist(SigDist_TC)
  # hist(SigDist_RAIN)
  # hist(SigDist_JTK)

  # par(mfrow = c(1, 1), font = 2)
  lineWidth <- 6
  plot(ecdf(SigDist_JTK), pch = NA, xlim = range(c(SigDist_TC, SigDist_JTK, SigDist_RAIN, OrderedDist, NullDist)), lwd = lineWidth, col = "#377eb8", main = paste0(ds, " Rank Corr eCDF"), ylab = "eCDF", xlab = expression(paste("Rank Correlation (", rho, ")")))
  plot(ecdf(SigDist_TC), add = TRUE, lwd = lineWidth, col = "#ff7f00")
  plot(ecdf(SigDist_RAIN), add = TRUE, lwd = lineWidth, col = "black")
  plot(ecdf(OrderedDist), add = TRUE, lwd = lineWidth, col = "green3")
  plot(ecdf(NullDist), add = TRUE, lwd = lineWidth, col = "grey")
  legend(
    x = "topleft", inset = c(0.01, 0.05), legend = c("TimeCycle", "JTK_CYCLE", "RAIN", "All Gene Correlation", "Null Distribution"),
    col = c("#ff7f00", "#377eb8", "Black", "green3", "grey"), lty = 1, lwd = 4, cex = 1, box.lwd = 0
  )
  box(lwd = 2)
}
dev.off()



# pdf("~/Desktop/reproducibilityAnalysis.pdf", width = 32, height = 20)
# plot.new()
# dev.off()
