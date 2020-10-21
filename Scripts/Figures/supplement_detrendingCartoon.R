
# initialize dataset
tp <- seq(from = 0, to = 24, by = .1)
amp <- 1
growth <- -1
shift <- 12
sigmoid <- amp / (1 + exp(growth * (tp - shift)))

data <- data.frame(
  x = tp,
  y = sigmoid
)
linearFit <- lm(y ~ x, data = data)

# detrend sigmoid with
data <- cbind(data, data.frame(detrend = sigmoid - predict(linearFit)))
checkDetrended <- lm(detrend ~ x, data = data)

# plotting parameters
colLT <- "#ff7f00"
colSig <- "#377eb8"
lwdLT <- 3
lwdSig <- 6
ltyLT <- 5

# generate the plot
pdf(file = "~/Desktop/TimeCycle-data/Results/Figures/detrend.pdf", width = 8, height = 4)

par(mfrow = c(1, 2), font.lab = 2, font.axis = 2)
plot(
  x = data$x,
  y = data$y,
  type = "l",
  lwd = lwdSig,
  col = colSig,
  main = "Sigmoid",
  xlab = "ZT Time",
  ylab = "Expression",
  axes = F
)
abline(linearFit, lwd = lwdLT, lty = ltyLT, col = colLT)
axis(side = 1, at = seq(0, 24, 4), labels = seq(0, 24, 4))
axis(side = 2)
box(lwd = 3)

plot(
  x = data$x,
  y = data$detrend,
  type = "l",
  lwd = lwdSig,
  col = colSig,
  main = "Detrended Sigmoid",
  xlab = "ZT Time",
  ylab = "Expression",
  axes = F
)
abline(checkDetrended, lwd = lwdLT, lty = ltyLT, col = colLT)
axis(side = 1, at = seq(0, 24, 4), labels = seq(0, 24, 4))
axis(side = 2)
box(lwd = 3)


dev.off()
