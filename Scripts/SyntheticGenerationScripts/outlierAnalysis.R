# select inputs file

reps <- list(1, 2, 3)
lengths <- list(48)
splRate <- list(1, 2, 4)
noiseLevels <- list(seq(0, 0.4, .1))
set.seed(123)

fileNames <- sapply(reps, function(rep) {
  sapply(lengths, function(lng) {
    sapply(noiseLevels, function(nl) {
      sapply(splRate, function(smps) {
        runName <- paste(smps, lng, "NoiseLV", nl, "BioRep", rep, sep = "_") # create a Name for the sample Run
      })
    })
  })
})

outlierAnalysisFiles <- sapply(as.vector(unlist(fileNames)), function(fileName) {
  fileLocation <- paste0("../../Data/Processed/", fileName, ".txt")
  data <- read.delim(file = fileLocation, row.names = 1)
  outlierData <- data

  # randomly select 1 percent of points to be outliers
  nPoints <- nrow(data) * ncol(data)
  randomPoints <- sample(x = 1:nPoints, size = (nPoints * .01), replace = F)
  randomMatrixOfPoints <- matrix(data = 1:nPoints, nrow = nrow(data), ncol = ncol(data), byrow = T) %in% randomPoints

  # create a matrix of randomly selected points to be outliers
  outliers <- matrix(data = randomMatrixOfPoints, nrow = nrow(data), ncol = ncol(data), byrow = T)

  # outlier points are between 3 to 4 standard deviation above or below the time-series diurnal mean
  for (i in 1:nrow(outlierData)) {
    for (j in 1:ncol(outlierData)) {
      if (outliers[i, j]) {
        tsMean <- mean(unname(unlist(outlierData[i, ])))
        tsSd <- sd(unname(unlist(outlierData[i, ])))
        lowerBound <- runif(n = 1, min = tsMean - (4 * tsSd), max = tsMean - (3 * tsSd))
        upperBound <- runif(n = 1, min = tsMean + (3 * tsSd), max = tsMean + (4 * tsSd))
        outlierData[i, j] <- sample(x = c(lowerBound, upperBound), size = 1, replace = F)
      }
    }
  }


  # Create Dataframe of Samples
  outlierData <- as.data.frame(outlierData)
  outlierData$"probes" <- rownames(data)
  outlierData <- outlierData[c(ncol(outlierData), 1:(ncol(outlierData) - 1))]

  write.table(outlierData, file = paste("../../Data/Outliers/", fileName, ".txt", sep = ""), sep = "\t", row.names = F, col.names = T, quote = FALSE)
})

source("../TimeCycleV3.R")
rm(list = setdiff(ls(), c("fileNames", "averageReps", "getRepAvgedDataFrame")))

outlierAnalysisFilesSw1pers <- sapply(as.vector(unlist(fileNames)), function(fileName) {
  fileLocation <- paste0("../../Data/Outliers/", fileName, ".txt")
  data <- read.delim(file = fileLocation, row.names = 1)

  # Get Meta Data From File Name
  metaData <- as.numeric(unlist(regmatches(fileName, gregexpr("[[:digit:]]+\\.*[[:digit:]]*", fileName))))
  timePoints <- seq(from = 0, to = metaData[2], by = metaData[1])

  data <- getRepAvgedDataFrame(data = data, repLabel = rep(metaData[4], length(timePoints)))

  # Add Column for Names of Samples
  data$"#" <- rownames(data)
  data <- data[c(ncol(data), 1:(ncol(data) - 1))]

  write.table(data, file = paste("../../Data/sw1per_Outliers/", fileName, ".txt", sep = ""), sep = "\t", row.names = F, col.names = T, quote = FALSE)
})
