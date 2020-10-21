##########################################
# Plot Results
##########################################
load("~/Desktop/TimeCycle-data/Results/Downsampling/24HourComparison.Rdata")
library(UpSetR)
library(dplyr)

# Check to see Correlation between genes
getCorrelation <- function(controlData, comparisonData) {
  output <- sapply(1:nrow(controlData), function(i) {
    set1 <- unlist(unname(controlData[i, ]))
    set2 <- unlist(unname(comparisonData[i, ]))
    cor(set1, set2, method = "spearman")
  })
  return(output)
}

#############################
# JTK Results
#############################

# upset Plots
geneNames <- rownames(JTK_results)
JTK_sigGenes <- apply(JTK_results, 2, function(pVals) {
  geneNames[pVals < 0.05]
})

upset(
  data = fromList(JTK_sigGenes),
  order.by = "freq",
  mb.ratio = c(0.6, 0.4),
  nsets = 24,
  sets = paste0("JTK_ZT_", seq(from = 1, by = 1, length.out = 24)),
  keep.order = T,
  number.angles = 0,
  text.scale = 1.1,
  set_size.scale_max = 2000,
  point.size = 2.8,
  set_size.show = T,
  line.size = 1,
  group.by = "degree"
)

#############################
# Rain Results
#############################


RAIN_sigGenes <- apply(RAIN_results, 2, function(pVals) {
  geneNames[p.adjust(p = pVals, method = "fdr") < 0.05]
})

upset(
  data = fromList(RAIN_sigGenes),
  order.by = "freq",
  mb.ratio = c(0.6, 0.4),
  nsets = 24,
  sets = paste0("RAIN_ZT_", seq(from = 1, by = 1, length.out = 24)),
  keep.order = T,
  number.angles = 0,
  set_size.scale_max = 2000,
  text.scale = 1.1,
  point.size = 2.8,
  set_size.show = T,
  line.size = 1,
  group.by = "degree"
)


# ecdf Plots
#overlap analysis counts of sig genes across datasets
RAIN_Count <- table(colSums(t(RAIN_results < 0.05)))
JTK_Count <- table(colSums(t(JTK_results < 0.05)))
JTK_Per <- JTK_Count / 12868
RAIN_Per <- RAIN_Count / 12868

# ecdf Plots

par(mfrow = c(1, 3), font.lab = 2, font.axis = 2, mar = c(5.1,5.1, 4.1, 2.1))
plot(JTK_Count,
     type = "o",
     pch = 16,
     xlab = "Number of Overlaps Between Datasets",
     ylab = "Number of Detected Cycling Genes",
     col = "#377eb8",
     ylim = range(c(JTK_Count, RAIN_Count)),
     main = "Detected Genes By Overlap",
     axes = F
)
lines(RAIN_Count,
      type = "o",
      pch = 16,
      lwd = 2
)
axis(side = 1, at = seq(0,24,2),labels = seq(0,24,2))
axis(side = 2)
box(lwd = 3)
plot(JTK_Per,
     ylim = c(0,1),
     type = "o",
     pch = 16,
     xlab = "Number of Overlaps Between Datasets",
     ylab = "Percentage Genes Detected Cycling\n Across all Cycling Genes Detected",
     col = "#377eb8",
     main = "Percent Detected By Overlap",
     axes = F
)
lines(RAIN_Per,
      type = "o",
      pch = 16,
      lwd = 2
)
axis(side = 1, at = seq(0,24,2),labels = seq(0,24,2))
axis(side = 2)
box(lwd = 3)
plot(cumsum(rev(JTK_Per)),
     cex = 0, 
     xlab = "Number of Overlaps Between Datasets",
     ylab = "Percentage Genes Detected Cycling\n Across all Cycling Genes Detected",
     main = "CDF of Percent Detected",
     ylim = c(0, 1),
     axes = F
)
abline(v = 13, col = "grey", lwd = 1.5)
lines(cumsum(rev(JTK_Per)),
      col = "#377eb8",
      type = "o",
      pch = 16,
      lwd = 2,
)
lines(cumsum(rev(RAIN_Per)),
      type = "o",
      pch = 16,
      lwd = 2
)
axis(side = 1, at = seq(-1,26,2),labels = rev(seq(0,26,2)))
axis(side = 2)
box(lwd = 3)



#ecdf Plots
#Correlation Across Data Sets
jtk_sg <- apply(JTK_results_Dataset_Comp, 2, function(pVals) {
  geneNames[pVals < 0.05]
})
rain_sg <- apply(RAIN_results_Dataset_Comp, 2, function(pVals) {
  geneNames[p.adjust(pVals, "fdr") < 0.05]
})

upset(fromList(jtk_sg),
      order.by = "freq",
      mb.ratio = c(0.6, 0.4),
      nsets = 24,
      keep.order = T,
      number.angles = 0,
      set_size.scale_max = 2000,
      text.scale = 1.1,
      point.size = 2.8,
      set_size.show = T,
      line.size = 1,
      group.by = "degree")
upset(fromList(rain_sg),
      order.by = "freq",
      mb.ratio = c(0.6, 0.4),
      nsets = 24,
      keep.order = T,
      number.angles = 0,
      set_size.scale_max = 2000,
      text.scale = 1.1,
      point.size = 2.8,
      set_size.show = T,
      line.size = 1,
      group.by = "degree")

for (ds in c("Hughes", "Hogenesch", "Zhang")){
  
overlap_JTK <- jtk_sg[[ds]]
overlap_RAIN <- rain_sg[[ds]]

# Null Distribution - No ordering for gene sample. randomly drawn with replacement.
Null_Cor_Zhang_Hughes <- getCorrelation(controlData = sample_n(Zhang, size = 12868, replace = T), comparisonData = sample_n(Hughes, size = 12868, replace = T))
Null_Cor_Zhang_Hogenesch <- getCorrelation(controlData = sample_n(Zhang, size = 12868, replace = T), comparisonData = sample_n(Hogenesch, size = 12868, replace = T))
Null_Cor_Hughes_Hogenesch <- getCorrelation(controlData = sample_n(Hughes, size = 12868, replace = T), comparisonData = sample_n(Hogenesch, size = 12868, replace = T))

# Distribution given gene ordered samples, but no filter for sig genes
Ordered_Cor_Zhang_Hughes <- getCorrelation(controlData = Zhang, comparisonData = Hughes)
Ordered_Cor_Zhang_Hogenesch <- getCorrelation(controlData = Zhang, comparisonData = Hogenesch)
Ordered_Cor_Hughes_Hogenesch <- getCorrelation(controlData = Hughes, comparisonData = Hogenesch)

# RAIN - Distribution given gene ordered samples and filtered for sig genes
Sig_Cor_Zhang_Hughes_RAIN <- getCorrelation(controlData = Zhang[overlap_RAIN,], comparisonData = Hughes[overlap_RAIN,])
Sig_Cor_Zhang_Hogenesch_RAIN <- getCorrelation(controlData = Zhang[overlap_RAIN,], comparisonData = Hogenesch[overlap_RAIN,])
Sig_Cor_Hughes_Hogenesch_RAIN <- getCorrelation(controlData = Hughes[overlap_RAIN,], comparisonData = Hogenesch[overlap_RAIN,])

# JTK - Distribution given gene ordered samples and filtered for sig genes
Sig_Cor_Zhang_Hughes_JTK <- getCorrelation(controlData = Zhang[overlap_JTK,], comparisonData = Hughes[overlap_JTK,])
Sig_Cor_Zhang_Hogenesch_JTK <- getCorrelation(controlData = Zhang[overlap_JTK,], comparisonData = Hogenesch[overlap_JTK,])
Sig_Cor_Hughes_Hogenesch_JTK <- getCorrelation(controlData = Hughes[overlap_JTK,], comparisonData = Hogenesch[overlap_JTK,])

#merge correlations into single vector
NullDist <- c(Null_Cor_Zhang_Hughes, Null_Cor_Zhang_Hogenesch, Null_Cor_Hughes_Hogenesch)
OrderedDist <- c(Ordered_Cor_Zhang_Hughes, Ordered_Cor_Zhang_Hogenesch, Ordered_Cor_Hughes_Hogenesch)
SigDist_RAIN <- c(Sig_Cor_Zhang_Hughes_RAIN, Sig_Cor_Zhang_Hogenesch_RAIN, Sig_Cor_Hughes_Hogenesch_RAIN)
SigDist_JTK <- c(Sig_Cor_Zhang_Hughes_JTK, Sig_Cor_Zhang_Hogenesch_JTK, Sig_Cor_Hughes_Hogenesch_JTK)

par(mfrow = c(1, 4))
hist(NullDist)
hist(OrderedDist)
hist(SigDist_RAIN)
hist(SigDist_JTK)

par(mfrow = c(1, 1), font = 2)
lineWidth <- 6
plot(ecdf(SigDist_JTK), pch = NA, xlim = range(c(SigDist_JTK, SigDist_RAIN, OrderedDist, NullDist)), lwd = lineWidth, col = "#377eb8", main = "eCDF")
plot(ecdf(SigDist_RAIN), add = TRUE, lwd = lineWidth, col = "black")
plot(ecdf(OrderedDist), add = TRUE, lwd = lineWidth, col = "green3")
plot(ecdf(NullDist), add = TRUE, lwd = lineWidth, col = "grey")
legend(
  x = "topleft", inset = c(0.01, 0.05), legend = c("JTK_CYCLE", "RAIN", "All Gene Correlation", "Null Distribution"),
  col = c("#377eb8", "Black", "green3", "grey"), lty = 1, lwd = 4, cex = 1, box.lwd = 0
)
box(lwd = 2)
}

