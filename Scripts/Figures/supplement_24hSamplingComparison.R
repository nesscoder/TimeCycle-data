##########################################
# Plot Results
##########################################
source("~/Desktop/TimeCycle-data/Scripts/Figures/2_AUCcomparisonOverview.R")
rm(list=setdiff(ls(), c("toPlot")))
load("~/Desktop/TimeCycle-data/Results/SlidingWindows/24HourComparison.Rdata")
library(dplyr)
library(ggalluvial)
library(gridExtra)
library(pals)

getSigGenesCount <- function(data, alpha1, alpha2) {
  geneNames <- rownames(data)
  sigGenes_0.05 <- colSums(t(data[, START:END] < alpha1))
  sigGenes_0.1 <- colSums(t(data[, START:END] < alpha2))
  
  df_sigGenes_0.05 <- data.frame(gene = geneNames, overlaps = sigGenes_0.05, threshold = "alpha1")
  df_sigGenes_0.1 <- data.frame(gene = geneNames, overlaps = sigGenes_0.1, threshold = "alpha2")
  
  return(rbind(df_sigGenes_0.05, df_sigGenes_0.1))
}


plotAlluvial <- function(data, plotTitle, col) {
  theme_set(new = theme_light())
  theme_replace(
    panel.border = element_rect(colour = "black", fill = NA, size = 2),
    text = element_text(size = 16, face = "bold"))
  
  data$overlaps <- factor(as.character(data$overlaps), levels = 0:24)
  data$gene <- as.factor(data$gene)
  data$threshold <- as.factor(data$threshold)
  
  ggplot(data = data, mapping = aes(
    x = threshold, 
    stratum = overlaps, 
    alluvium = gene,
    fill = overlaps, 
    label = overlaps)) +
    scale_fill_manual(values = col) +
    geom_flow() +
    geom_stratum() +
    scale_x_discrete(labels = c("FDR < 0.05", "FDR < 0.1")) +
    theme(legend.position = "none", panel.grid.major.x = element_blank()) +
    guides(fill = guide_legend(override.aes = list(linetype = 0, shape = 22, size = 5), title = "Number of Overlaps", title.position = "top", ncol = 5, byrow = TRUE)) +
    xlab(label = "") +
    ggrepel::geom_text_repel(data = subset(data, subset = as.numeric(threshold) == 1), stat = "stratum", aes(label = 100 * round(after_stat(prop), 3)),
                             size = 4, direction = "y", nudge_x = -0.5) +
    ggrepel::geom_text_repel(data = subset(data, subset = as.numeric(threshold) == 2), stat = "stratum", aes(label = 100 * round(after_stat(prop), 3)),
                             size = 4, direction = "y", nudge_x = 0.5) +
    ggtitle(plotTitle)
}
START <- 1
END <- 24

dataToPlot <- as.data.frame(toPlot) %>%
  filter(Length == 24)

theme_set(new = theme_light())
theme_replace(
  panel.border = element_rect(colour = "black", fill = NA, size = 2),
  text = element_text(size = 16, face = "bold")
)

pdf("~/Desktop/TimeCycle-data/Results/Figures/AUC24Plot.pdf", width = 7 * 1.5 / 4.5, height = 3.5 * 1.5)
ggplot(data = dataToPlot,
  aes( x = as.numeric(numSamples),
    y = as.numeric(AUC),
    fill = Method,
    colour = Method)
) +
  geom_point() +
  geom_point(shape = 21, colour = "black", show.legend = F) +
  geom_smooth(method = lm, se = FALSE, size = 1) +
  scale_fill_manual(
    labels = c("JTK_CYCLE", "RAIN", "SW1PERS", "TimeCycle"),
    values = c("#377eb8", "black", "grey", "#ff7f00")
  ) +
  scale_colour_manual(
    labels = c("JTK_CYCLE", "RAIN", "SW1PERS", "TimeCycle"),
    values = c("#377eb8", "black", "grey", "#ff7f00")
  ) +
  facet_grid(Interval ~ Length) +
  xlab("Number of Samples") +
  ylab("AUC") +
  guides(colour = guide_legend(override.aes = list(linetype = 0, shape = 22, size = 5)), fill = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom", legend.title = element_blank())
dev.off()

TC_sigGenes <- getSigGenesCount(data = TC_results, alpha1 = 0.05, alpha2 = 0.1)
JTK_sigGenes <- getSigGenesCount(data = JTK_results, alpha1 = 0.05, alpha2 = 0.1)
RAIN_sigGenes <- getSigGenesCount(data = RAIN_results, alpha1 = 0.05, alpha2 = 0.1)

customColors <- cols25(n = 25)
pdf("~/Desktop/TimeCycle-data/Results/Figures/AluvOut24.pdf", width = 14, height = 5)
p1 <- plotAlluvial(data = TC_sigGenes, plotTitle = "TimeCycle", col = customColors[1])
p2 <- plotAlluvial(data = JTK_sigGenes, plotTitle = "JTK_CYCLE", col = customColors)
p3 <- plotAlluvial(data = RAIN_sigGenes, plotTitle = "RAIN", col = customColors)
grid.arrange(grobs = list(p1, p2, p3), nrow = 1)
dev.off()