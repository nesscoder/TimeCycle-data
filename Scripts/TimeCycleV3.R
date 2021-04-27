#========================================================#
# Author: Elan Ness-Cohn                                 #
# Email: elanness-cohn2017@u.northwestern.edu            #
# Title: Time Cycle V3                                   #
# Code Produced: February 19, 2019                       #
# Code Updated: October  19, 2020                        #
#========================================================#


#' Main TimeCycle Function
#'
#' Main function of the \pkg{TimeCycle} package used for detecting rhythmic signals in time-series gene expression sets.
#' For additional help with parameter selection, see TimeCycle's vignette:
#' \code{vignette("TimeCycle")}.
#'
#' @param data a \code{data.frame} of numeric gene expression over time (row = genes \emph{x} col = ZT times).
#' @param repLabel a \code{vector} defining the number of replicates at each time points.
#' @param resamplings a \code{numeric} specifying the number of resamplings to use in the null-distribution. Default is \code{10000}.
#' @param minLag a \code{numeric} specifying the min lag to check in the 3-D embedding. Default is \code{2}.
#' @param maxLag a \code{numeric} specifying the max lag to check in the 3-D embedding. Default is  \code{5}.
#' @param cores a \code{numeric} specifying the number of parallel cores to use. Default number of cores is \code{parallel::detectedCores() - 2}.
#' @param period a \code{numeric} specifying the period of interest in hours for rhythm detection. Default is \code{24}.
#' @param laplacian  a \code{logical} scalar. Should the Laplacian Eigenmaps be used for dimensionality reduction? Default \code{TRUE}
#' @param linearTrend  a \code{logical} scalar. Should TimeCycle Prioritize detecting linear trending signals? Default \code{FALSE}. Not recommended to change from default \code{FALSE} - will increases false positives rate. See vignette("TimeCycle") for more details.
#'
#' @references{
#' \itemize{
#'    \item{A pre-print version of TimeCycle is available on BioRxiv at \url{https://doi.org/10.1101/2020.11.19.389981}}
#'    }
#'    \subsection{TDA Package References}{
#'    \itemize{
#'      \item Wadhwa RR, Williamson DFK, Dhawan A, Scott JG. (2018). "TDAstats: R pipeline for computing persistent homology in topological data analysis." \emph{Journal of Open Source Software}. 2018; 3(28): 860. doi:\href{https://doi.org/10.21105/joss.00860}{[10.21105/joss.00860]}
#'      \item Bauer U. (2019). "Ripser: Efficient computation of Vietoris-Rips persistence barcodes." \emph{arXiv}: 1908.02518.
#'    }
#'
#'
#'    }
#'}
#'
#'
#'
#' @return a tidy \code{data.frame} of processed results for each gene:
#' \tabular{lcccccc}{
#'  \strong{sampleName} \tab  \strong{perScore} \tab \strong{pVals} \tab \strong{pVals.adj} \tab \strong{Period.in.Hours} \tab \strong{Amp} \tab \strong{Phase.In.Hours} \cr
#'  the gene name \tab  the median persistence score across all lags (min to max)  \tab raw empirical \emph{p}-value \tab FDR adjusted \emph{p}-value \tab period (h) \tab amplitude \tab phase (h) \cr
#' }
#'
#' @examples
#' # use built in zhang2014 data set sampled every
#' # 2 hours for 48 hours (i.e. 24 time points with 1 replicate each).
#' # Search for genes with period of 24 hours.
#'
#' #set seed for reproducibility with random variables in example usage
#' set.seed(1234)
#'
#' TimeCycleResults <- TimeCycle(data = zhang2014[1:100,],
#'                               repLabel = rep(1,24),
#'                               period = 24,
#'                               cores = 2,
#'                               resamplings = 10)
#'
#' # Check number of genes with FDR < 0.05 and period between 22 to 26 hours.
#' library(tidyverse)
#'
#' TimeCycleResults %>%
#'    filter(22 < Period.in.Hours & Period.in.Hours < 26) %>%
#'    filter(pVals.adj < 0.05) %>%
#'    glimpse()
#'
#'@export
TimeCycle <- function(data,
                      repLabel,
                      resamplings = 10000,
                      minLag = 2,
                      maxLag = 5,
                      period = 24,
                      cores = parallel::detectCores()-2,
                      laplacian = T,
                      linearTrend = F
){

  ## -----------------------------------pre-process Data ----------------------------

  cat("
      ########################################################################################
      ###      ████████ ██ ███    ███ ███████  ██████ ██    ██  ██████ ██      ███████     ###
      ###         ██    ██ ████  ████ ██      ██       ██  ██  ██      ██      ██          ###
      ###         ██    ██ ██ ████ ██ █████   ██        ████   ██      ██      █████       ###
      ###         ██    ██ ██  ██  ██ ██      ██         ██    ██      ██      ██          ###
      ###         ██    ██ ██      ██ ███████  ██████    ██     ██████ ███████ ███████     ###
      ########################################################################################\n")

  #remove before launching
  # set.seed(123)

  startTime <- Sys.time()
  print("Starting TimeCycle")
  print("Pre-Processing Data")

  #reorder the Data By Replicate
  dataAvg <<- getRepAvgedDataFrame(data = data, repLabel = repLabel)

  #impute Missing Values
  dataAvg <- imputeMissingData(dataAvg)

  #center the data and merge data back together
  centeredData <- meanCenter(dataAvg)

  #detrend the data and ONLY USE FOR Compute the period
  dataDetrend <- detrend(centeredData)

  print("Computing Periods")
  #Compute Periods of Replicate 1 Time Series from the detrended data - does not break FFT assumptions
  periods <- periodFinder(movingAverageDF(dataDetrend))

  #scale Data Between 0 and 1 across all genes
  #allow us to compare all genes at once to the null distribution rather than just one at a time
  dataScaled  <- as.data.frame(t(apply(centeredData, 1, scaleTimeSeries)))
  colnames(dataScaled) <- colnames(centeredData)

  #Smooth data with Autocorrelation
  preProcessedData <- preprocess_acf(dataScaled, period, linearTrend = linearTrend)

  ##----------------------- pre-process NUll Distribution Data -----------------------
  print("Pre-Processing Null Distribution")

  #create the Null resampling of the Data
  resampledDataForNull <- nullResampling(dataScaled, numExperiments = resamplings)

  #Smooth -> Mean Center
  resampledprocessedData <- preprocess_acf(resampledDataForNull, period, linearTrend = linearTrend)

  ##------------------------ Compute the NUll Distribution  ------------------------
  print("Computing Null Distribution")

  nullDist <- computePersistence(resampledprocessedData, minLag = minLag, maxLag = maxLag, cores = cores, laplacian = laplacian)

  #Order the Null Distribtion to Speed up Ranking
  nullDistOrder <- as.numeric(nullDist[order(as.vector(nullDist))])

  ##--------------------Calculate the Peristence Score for Each Sample ------------------------
  print("Computing Persistence Scores")
  dataPS <- computePersistence(preProcessedData, minLag = minLag, maxLag = maxLag, cores = cores, laplacian = laplacian)

  ##-----------------------    Get pVals for Data ------------------------
  print("Calculating p-values")
  perScore <- as.numeric(dataPS)
  sampleNames <- rownames(preProcessedData)
  vect <- 1:length(perScore)

  #Compute Emperical pVal
  pVals <- parallel::mclapply(X = vect, mc.cores = cores, function(x){
    return(rank(-c(perScore[x],nullDistOrder),ties.method = "max")[1])
  })

  #get FDR Adjusted pVal
  pVals <- unlist(pVals)/as.numeric(resamplings)
  pVals.adj <- stats::p.adjust(pVals,method = 'fdr')
  names(pVals.adj) <- sampleNames

  endTime <- Sys.time()
  runTime <- getTimeSpan(startTime,endTime)

  results <- data.frame(sampleNames,perScore,pVals,pVals.adj,t(periods[,]))

  print(paste0("TimeCycle Completed"))
  print(paste0("Analysis Time: ", runTime))

  return(results)
}




#' Handles Sample Replicates in Single Time-Series
#'
#' Averages expression values for a single time-series across replicate groups.
#'
#' @param geneExpr a \code{vector} of \code{numeric} time-series expression values with replicates.
#' @param Reps a \code{vector} defining the number of replicates at each time point.
#'
#' @return a \code{vector} of average  \code{numeric} time-series expression values by replicate time-points.
#' @seealso \code{\link{getRepAvgedDataFrame}}
#'
#' @export
#'
#' @examples
#'
#' geneExpr <- c(1, 5, 3, 4, 5, 8)
#' reps <- c(2, 3, 1)
#' averageReps(geneExpr = geneExpr, Reps = reps)
averageReps <- function(geneExpr, Reps) {
  splitBy <- unlist(sapply(1:length(Reps), function(seq) {
    rep(seq, Reps[seq])
  }))
  as.vector(unlist(lapply(split(geneExpr, splitBy), mean)))
}



#' Generates Takens' Embedding For a Time-Series
#'
#' Generates the Takens' embedding for a time-series given a specified delay (i.e. lag) and dimension (i.e. 2-D or 3-D).
#'
#' @param x a \code{vector} of \code{numeric} time-series expression values.
#' @param dim a \code{numeric} specifying the dimension to use for in the time-delayed embedding (i.e. 2-D or 3-D).
#' @param delay a \code{numeric} specifying the lag to use for in the n-dimensional time delayed embedding specified by \code{dim}.
#'
#' @seealso \code{\link{getPersistence}}.
#' @return a \code{data.frame} of the n-dimensional Takens' embedding. Columns defined from (1-D to n-D).
#' @export
#'
#'
buildTakens_ndim <- function(x, dim, delay = 1) {
  n <- length(x) - (dim-1)*delay
  X <- seq_along(x)
  if(n <= 0)
    stop("Insufficient observations for the requested embedding")
  out <- matrix(rep(X[seq_len(n)], dim), ncol = dim)
  out[,-1] <- out[,-1, drop = FALSE] +
    rep(seq_len(dim - 1) * delay, each = nrow(out))

  out <- matrix(x[out], ncol = dim)

  return(out)
}




#' Computes the Laplacian Eigenmap For Dimension Reduction
#'
#' Takes a  n-D Takens' Embedding and returns the 2-D Laplacian Eigenmap.
#'
#' @param takens a \code{data.frame} of the n-dimensional Takens' embedding. Columns defined from (1-D to n-D). TimeCycle defaults to \code{n = 3}.
#'
#' @return a \code{data.frame} of the 2-dimensional Takens' embedding. Columns defined by eigen vectors with last two non-trivial eigen values.
#'
#' @seealso
#' \itemize{
#'      \item \code{\link{getPersistence}} for usage
#'      \item \code{\link{matrixLaplacian}} for defining laplacian eigenmaps.
#'}
#'
#'
#' @export
#'
computeLaplacianEmbedding <- function(takens) {
  A <- as.matrix(stats::dist(x = takens, method = "euclidian"))
  A <- 1 - A / max(A) # normalize

  Z <- matrixLaplacian(A)
  # select last two non trivial eigen values
  eigen1 <- dim(Z$eigenvector)[2] - 2
  eigen2 <- dim(Z$eigenvector)[2] - 1
  return(Z$eigenvector[, c(eigen1, eigen2)])
}




#' Computes Persistence Scores For a Data.Frame of Time-Series Across Multiple Lags
#'
#' Takes a \code{data.frame} of numeric gene expression over time (genes X ZT times) and computes the persistence score using \code{\link{getPersistence}}.
#' For a given gene, each lag (min to max) is used to transform the expression into a 3-D embedded space via time-delay embedding.
#' A non-linear dimension reduction technique (laplacian eigenmaps) is used to transfrom the 3-D embedding to a 2-D embedding.
#' Finally, the persistence score of the 2-D embedding is calculated via persistence homology.
#' The median persistence score across all lags (min to max) for each gene is returned as a numeric vector.
#' For more details see TimeCycle's vignette:
#' \code{vignette("TimeCycle")}.
#'
#' @param data a \code{data.frame} of \code{numeric} gene expression over time (row = genes \emph{x} col = ZT times).
#' @param minLag a \code{numeric} specifying the min lag to check in the 3-D embedding. Default is \code{2}.
#' @param maxLag a \code{numeric} specifying the max lag to check in the 3-D embedding. Default is \code{5}.
#' @param cores a \code{numeric} specifying the number of parallel cores to use. Default number of cores is \code{parallel::detectedCores() - 2}.
#' @param laplacian a \code{logical} scalar. Should the Laplacian Eigenmaps be used for dimensionality reduction? Default \code{TRUE}.
#'
#' @references{
#'    \itemize{
#'      \item Wadhwa RR, Williamson DFK, Dhawan A, Scott JG. (2018). "TDAstats: R pipeline for computing persistent homology in topological data analysis." \emph{Journal of Open Source Software}. 2018; 3(28): 860. doi:\href{https://doi.org/10.21105/joss.00860}{[10.21105/joss.00860]}
#'      \item Bauer U. (2019). "Ripser: Efficient computation of Vietoris-Rips persistence barcodes." \emph{arXiv}: 1908.02518.
#'    }
#' }
#' @seealso
#' \itemize{
#'      \item \code{\link[TDAstats]{calculate_homology}} for Persistence Homology calculation.
#'      \item \code{\link{buildTakens_ndim}} for for generating time-delay embedding.
#'      \item \code{\link{computeLaplacianEmbedding}} for 3-D to 2-D laplacian eigenmaps dimension reduction.
#'      \item \code{\link{getPersistence}} for use on a single gene expression time-series.
#'}
#'
#' @return a \code{vector} of the median persistence score across lags (minLag to maxLag) for each gene in data
#' @export
#'
computePersistence <- function(data, minLag = 2, maxLag = 5, cores = parallel::detectCores() - 2, laplacian = T){
  vect <- as.list(minLag:maxLag)
  #compute PS at each lag
  output <- parallel::mclapply(vect, mc.cores = cores, function(lag){
    perTSoutput <- apply(data,1,function(TS){
      return(getPersistence(t(as.matrix(TS)), lag = lag, laplacian = laplacian))
    })
    perTSoutput <- as.data.frame(perTSoutput)
  })

  #save persistence score at each lag
  PSatEachLag <- do.call(cbind,output)

  return(apply(PSatEachLag,1, stats::median))
}




#' Example data set of normalized microarray time-series expression set from Zhang et. al. of WT mouse liver sampled every 2-h for 48-h
#'
#' Data Preprocessing: In the accompanying paper, the CEL files from three mouse liver Affymetrix microarray time-series expression sets
#'  \eqn{-} Hogenesch 2009 - GSE11923 (Hughes et al., 2009), Hughes 2012 - GSE30411 (Hughes et al., 2012),
#' Zhang 2014 - GSE54650 (Zhang et al., 2014) \eqn{-} were downloaded from the Gene Expression Omnibus database (GEO).
#' In each experiment, wild-type C57BL/6J mice were entrained to a 12-h light, 12-h dark environment before
#' being released into constant darkness. Mouse age, length of entrainment, time of sampling,
#' and sampling resolution vary by experiment. The data were subsequently normalized by robust
#' multi-array average (rma) using the Affy R Package (Gautier et al., 2004) and
#' checked for quality control using the Oligo R Package (Carvalho and Irizarry, 2010),
#' following each package’s vignette, respectively. Since each GEO data set used a different microarray
#' platform (affy_mouse430_2, affy_moex_1_0_st_v1, affy_mogene_1_0_st_v1), each had a different set of probes.
#' A common set of features needed to be identified to compare across microarrays.
#' Probes for each data set were mapped to genes based on prealigned databases specific to each
#' microarray (mouse4302.db, moex10sttranscriptcluster.db, mogene10sttranscriptcluster.db).
#' Multiple probes corresponding to a single gene were aggregated by taking the mean expression.
#' A final 12,868 common set of genes across all three microarray platforms were used for subsequent analysis.
#' See the supplement for code.
#'
#' @docType data
#'
#'
#' @format A data frame with 12868 rows (genes) and 24 variables (ZT time)
#'
#' @keywords datasets
#'
#' @references Zhang, R, Lahens, NF, Ballance, HI, Hughes, ME, Hogenesch, JB (2014) A circadian gene expression atlas in mammals: implications for biology and medicine. Proc Natl Acad Sci USA 111:16219-16224.
#' \href{https://doi.org/10.1073/pnas.1408886111}{DOI: 10.1073/pnas.1408886111}
#'
#' @source \href{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE54650}{GSE54650}
#'
"zhang2014"
#' Linearly Detrends a \code{data.frame} of Gene Expression Time-Series
#'
#' Detrends a \code{data.frame} of gene expression over time (row = genes \emph{x} col = ZT times)
#' by fiting a linear model to each gene and removing the linear trend.
#'
#' @param data a \code{data.frame} of \code{numeric} gene expression over time (row = genes \emph{x} col = ZT times).
#'
#' @return a detrended \code{data.frame} of \code{numeric} gene expression over time (row = genes \emph{x} col = ZT times).
#'
#' @seealso \code{\link{periodFinder}}
#'
#' @export
#'
detrend <- function(data){

  sampleNames <- rownames(data)
  colNames <- colnames(data)
  #get numeric timepoints values from column Names
  xVals <- as.numeric(unlist(regmatches(colNames, gregexpr("[[:digit:]]+\\.*[[:digit:]]*",colNames))))

  #Use timeseries to compute the trend fit
  dataDetrend <- apply(data, 1, function(TS) {
    fit <- stats::lm(TS~xVals)
    y <- zapsmall(TS - (xVals*fit$coefficients[2] + fit$coefficients[1]),10)
  })

  #return the detrended data as a data.frame
  dataDetrend <- as.data.frame(t(as.data.frame(dataDetrend)))
  rownames(dataDetrend) <- sampleNames
  colnames(dataDetrend) <- colNames

  return(dataDetrend)
}



#' Computes the Fast Fourier Transform of a Time-Series
#'
#' Takes a time-series and sample frequency and computes the fast Fourier transform of the time-series.
#'
#' @param y a \code{vector} of numeric time-series expression values.
#' @param sampFreq a \code{vector} of numeric frequencies.
#'
#' @return a \code{data.frame} harmonic fits by frequency.
#'
#' @seealso \code{\link[stats]{fft}}
#' @export
#'
getFFT <- function(y, sampFreq) {
  N <- length(y)
  fk <- stats::fft(y) / N # normalize Data
  fk <- 2 * fk[1:((length(fk) / 2) + 1)] # DC comp + half of positives
  freq <- (0:(N - 1)) * sampFreq / N
  freq <- freq[(1:(length(fk)))]
  return(data.frame(fur = fk, freq = freq))
}




#' Computes Persistence Scores For a Single Time-Series Across a Single Lag
#'
#' Takes a \code{vector} of numeric gene expression over time and computes the persistence score.
#' The specified lag is used to transform the expression into a 3-D embedded space via time-delay embedding.
#' A non-linear dimension reduction technique (laplacian eigenmaps) is used to transfrom the 3-D embedding to a 2-D embedding.
#' Finally, the persistence score of the 2-D embedding is calculated via persistence homology.
#' Returns the Max persistence score, returns 0 if no persistence score exists.
#' For more details see TimeCycle's vignette:
#' \code{vignette("TimeCycle")}.
#'
#' @param timeSeries a \code{vector} of \code{numeric} time-series expression values.
#' @param lag a \code{numeric} specifying the Lag to use for in the 3-D time delayed embedding.
#' @param laplacian a \code{logical} scalar. Should the Laplacian Eigenmaps be used for dimensionality reduction? Default \code{TRUE}.
#'
#' @return the max persistence score at the specified lag, returns 0 if no persistence score exists.
#'
#' @references{
#'    \itemize{
#'      \item Wadhwa RR, Williamson DFK, Dhawan A, Scott JG. (2018). "TDAstats: R pipeline for computing persistent homology in topological data analysis." \emph{Journal of Open Source Software}. 2018; 3(28): 860. doi:\href{https://doi.org/10.21105/joss.00860}{[10.21105/joss.00860]}
#'      \item Bauer U. (2019). "Ripser: Efficient computation of Vietoris-Rips persistence barcodes." \emph{arXiv}: 1908.02518.
#'    }
#' }
#' @seealso
#' \itemize{
#'      \item \code{\link[TDAstats]{calculate_homology}} for Persistence Homology calculation.
#'      \item \code{\link{buildTakens_ndim}} for for generating time-delay embedding.
#'      \item \code{\link{computeLaplacianEmbedding}} for 3-D to 2-D laplacian eigenmaps dimension reduction.
#'      \item \code{\link{computePersistence}} for use parallelized function for a \code{data.frame} of gene expression.
#'    }
#' @export
#'
getPersistence <- function(timeSeries, lag, laplacian = T) {
  
  
  # embed the Points
  embedding <- tryCatch(
    {
      if (laplacian) {
        buildTakens_ndim(x = as.vector(unlist(timeSeries)), dim = 3, delay = lag)
      } else {
        buildTakens_ndim(x = as.vector(unlist(timeSeries)), dim = 2, delay = lag)
      }
    },
    warning = function(warning_condition) {
    },
    error = function(error_condition) {
      # print("Insufficient observations for the requested embedding")
      return(NULL)
    },
    finally = {
    }
  )
  
  
  if (is.null(embedding)) {
    # print("Insufficient observations for the requested embedding, setting Persistence Score to 0")
    return(0)
  }
  
  # compute the Laplacian Embedding about the merged embedding
  embeddingLap <- tryCatch(
    {
      if (laplacian) {
        embeddingLap <- computeLaplacianEmbedding(embedding)
      } else {
        embeddingLap <- embedding
      }
    },
    warning = function(warning_condition) {
      "warning"
    },
    error = function(error_condition) {
      # print("Unable to Compute Eigen Vector, embedding colapsed")
      return(NULL)
    },
    finally = {
    }
  )
  
  if (is.null(embeddingLap)) {
    return(0) # setting Persistence Score to 0
  }
  
  # compute Rips Complex on the laplacian matrix
  # scale of Persistence to Check up to in persistence
  maxScale <- ceiling(2 * max(abs(as.vector(unlist(timeSeries)))))
  maxScale <- max(maxScale, 1, na.rm = T)
  persistence <- TDAstats::calculate_homology(embeddingLap, dim = 1, threshold = maxScale)
  # create a vectors loops formed by complex
  loops <- which(persistence[, 1] == 1)
  
  if (length(loops) > 0) {
    
    # select persistence loops with max diff between Birth and Death
    diff <- persistence[loops, 3] - persistence[loops, 2]
    maxdiff <- which.max(diff)
    maxPersit <- diff[maxdiff]
    return(maxPersit)
  } else {
    return(0)
  }
}




#' Handles Sample Replicates in a \code{data.frame} of Time-Series
#'
#' Averages expression values for a single time-series across replicate groups in a \code{data.frame}.
#'
#'
#' @param data a \code{data.frame} of \code{numeric} gene expression over time (row = genes \emph{x} col = ZT times).
#' @param repLabel a \code{vector} defining the number of replicates at each time point.
#'
#' @return a \code{data.frame} of average  \code{numeric} time-series expression values by replicate time-points for each gene.
#'
#' @seealso \code{\link{averageReps}}
#'
#' @export
#'
getRepAvgedDataFrame <- function(data, repLabel) {

  # get Column Names
  colnames <- colnames(data)
  # remove replicate label from colnames if they exist
  colnames <- gsub(pattern = "_rep.", replacement = "", colnames)
  # get unique ZT time for each point
  colnames <- unique(as.numeric(unlist(regmatches(colnames, gregexpr("[[:digit:]]+\\.*[[:digit:]]*", colnames)))))

  # get average of Replicate Labels
  output <- t(apply(data, MARGIN = 1, FUN = function(geneExprRow) {
    averageReps(geneExprRow, repLabel)
  }))
  output <- as.data.frame(output)
  colnames(output) <- colnames

  return(output)
}



#' Computes Algorithm Run Time
#'
#' Converts system time difference between TimeCycle's start and end time to
#' hh:mm:ss format.
#'
#' @param start a \code{numeric} defining the algorithm start time.
#' @param end a \code{numeric} defining the algorithm end time.
#'
#' @return a \code{numeric} defining the algorithm run time.
#' @export
#'
#' @examples
#' getTimeSpan(Sys.time(), Sys.time() + 600)
getTimeSpan <- function(start, end) {
  dsec <- round(abs(as.numeric(difftime(end, start, units = "secs"))))
  hours <- floor(dsec / 3600)
  minutes <- floor((dsec - 3600 * hours) / 60)
  seconds <- dsec - 3600 * hours - 60 * minutes
  paste0(
    sapply(c(hours, minutes, seconds), function(x) {
      formatC(x, width = 2, format = "d", flag = "0")
    }),
    collapse = ":"
  )
}



#' Imputates Missing Time-Points in Data
#'
#' Imputes \code{numeric} values for time-points with an \code{NA} by computing the linear path between missing points
#'
#' @param data a \code{data.frame} of \code{numeric} gene expression over time (row = genes \emph{x} col = ZT times) with missing values.
#'
#' @return a imputed \code{data.frame} of \code{numeric} gene expression over time (row = genes \emph{x} col = ZT times).
#' @export
#'
#' @seealso \code{\link[imputeTS]{na_interpolation}} for imputation procedure.
#'
#' @examples
#' a <- c(10, 12, 14, NA, 18)
#' b <- c(1, 2, NA, NA, 5)
#' data <- t(data.frame(a, b))
#' imputeMissingData(data)
imputeMissingData <- function(data) {
  colNames <- colnames(data)
  rowNames <- rownames(data)

  out <- t(apply(data, 1, function(geneExpr) {
    x <- as.numeric(geneExpr)
    imputed <- imputeTS::na_interpolation(x)
  }))

  colnames(out) <- colNames
  rownames(out) <- rowNames

  return(out)
}




#' Computes Laplacian Eigenmaping
#'
#' Computes the Laplacian matrix and eigenvectors of the n-D Takens' embedding.
#'
#' @param A a distance \code{matrix} of the n-dimensional Takens' embedding.
#'
#' @return a \code{list} of the laplacian matrix and eigenvectors.
#'
#' @seealso \code{\link{computeLaplacianEmbedding}}
#'
#' @export
#'
matrixLaplacian <- function(A) {
  B <- A
  D <- matrix(0, nrow = dim(A)[1], ncol = dim(A)[1])
  diag(D) <- B %*% rep(1, dim(A)[1])
  diag(D) <- 1 / sqrt(diag(D))
  Q <- D %*% B %*% D
  N <- diag(1, dim(Q)[1]) - Q
  Eigen <- eigen(N)$vectors

  object <- list(LaplacianMatrix = N, eigenvector = Eigen)
  return(object)
}



#' Mean Centers a \code{data.frame}
#'
#' Mean centers a \code{data.frame} by row.
#'
#' @param  df a \code{data.frame} of \code{numerics}.
#'
#' @return a mean center \code{data.frame} of \code{numerics} by row.
#'
#' @seealso \code{\link{preprocess_acf}}, \code{\link{scaleTimeSeries}}
#'
#' @export
#'
meanCenter <- function(df){
  as.data.frame(t(apply(df, 1, function(y) y-mean(y))))
}



#' Computes the Moving Average of a Single Time-Series
#'
#' Computes the moving average about a time-series defined by a specified number of points.
#'
#' @param x a \code{vector} of \code{numeric} time-series expression values.
#' @param n a \code{numeric} specifying the number of points to use in the moving average. Default \code{n = 3}.
#' @param centered a \code{logical} scalar. Should the moving average be centered about the current points? Default \code{TRUE} (i.e. average of current point (\code{p}) with  \code{p - n/2} and \code{p + n/2}).
#'
#' @return a \code{vector} containing the smoothed \code{numeric} moving average time-series expression values.
#'
#' @seealso \code{\link{movingAverageDF}}
#'
#' @export
#'
movingAverage <- function(x, n = 3, centered = TRUE) {
  if (centered) {
    before <- floor((n - 1) / 2)
    after <- ceiling((n - 1) / 2)
  } else {
    before <- n - 1
    after <- 0
  }

  # Track the sum and count of number of non-NA items
  s <- rep(0, length(x))
  count <- rep(0, length(x))

  # Add the centered data
  new <- x
  # Add to count list wherever there isn't a
  count <- count + !is.na(new)
  # Now replace NA_s with 0_s and add to total
  new[is.na(new)] <- 0
  s <- s + new

  # Add the data from before
  i <- 1
  while (i <= before) {
    # This is the vector with offset values to add
    new <- c(rep(NA, i), x[1:(length(x) - i)])

    count <- count + !is.na(new)
    new[is.na(new)] <- 0
    s <- s + new

    i <- i + 1
  }

  # Add the data from after
  i <- 1
  while (i <= after) {
    # This is the vector with offset values to add
    new <- c(x[(i + 1):length(x)], rep(NA, i))

    count <- count + !is.na(new)
    new[is.na(new)] <- 0
    s <- s + new

    i <- i + 1
  }

  # return sum divided by count
  return(s / count)
}




#' Computes the Moving Average By Row in a \code{data.frame} of Time-Series
#'
#' Computes the moving average about each time-series in a \code{data.frame}.
#'
#' @param data a \code{data.frame} of \code{numeric} time-series expression values.
#'
#' @return a \code{data.frame} containing the smoothed \code{numeric} moving average time-series expression values by row.
#'
#' @seealso \code{\link{movingAverage}} for parameter definitions
#'
#' @export
#'
movingAverageDF <- function(data){

  sampleNames <- rownames(data)
  colNames <- colnames(data)

  #smooth the data with an moving average filter
  dataDenoise <- apply(data, 1, function(TS) {
    movingAverage(x = as.numeric(TS), n = 3, centered = T)
  })

  output <- as.data.frame(t(as.data.frame(dataDenoise)))
  rownames(output) <- sampleNames
  colnames(output) <- colNames

  return(output)
}




#' Generates a \code{data.frame} of Resampled Time-Series for Computing the Null Distribution
#'
#' Generates a \code{data.frame} null resampling time-series as defined by the difference between consecutive points across all time-series.
#' For additional information about the null distribution used by TimeCycle, see TimeCycle's vignette:
#' \code{vignette("TimeCycle")}.
#'
#' @param data a \code{data.frame} of \code{numeric} gene expression over time (row = genes \emph{x} col = ZT times).
#' @param numExperiments a \code{numeric} specifying the number of resampling to use in the null-distribution. Default is \code{10000}.
#'
#' @return a \code{data.frame} of \code{numeric} resampled time-series (row = numExperiments \emph{x} col = ZT times).
#'
#' @seealso \code{\link{resampleTimeSeries}}
#'
#' @export
#'
nullResampling <- function(data, numExperiments = 10000){

  #Get Average Value of Data across Replicate
  dataForResampling <- data
  colNames <- colnames(data)

  #get Resampled Data for Use in generating Null Distribution
  diffBetweenTP <- t(apply(dataForResampling, 1, function(geneName){
    #Get Difference Between Time Points
    diff(geneName)
  }))

  #Resample Differences
  resampledDataNull <- resampleTimeSeries(diffBetweenTP, numTP = dim(data)[2], numRsmps = numExperiments)
  colnames(resampledDataNull) <- colNames

  return(resampledDataNull)
}




#' Estimates the period, amp, and phase of a \code{data.frame} of time-series
#'
#' Estimates the period, amp, and phase of a \code{data.frame} of time-series via fast Fourier transform (FFT).
#' A model fit for the time-series is generated using the first 3 harmonics. The period, amp, and phase are
#' computed based on the aggregate fit.
#'
#' @param data a \code{data.frame} of \code{numeric} gene expression over time (row = genes \emph{x} col = ZT times).
#'
#' @return  a \code{data.frame} of \code{numeric} period, amplitude, and phase estimates for each gene.
#' @seealso
#' \itemize{
#'      \item \code{\link{getFFT}} for FFT calculation.
#'      \item \code{\link{detrend}} for linear detrending of time-series prior to periodFinder processing.
#'      \item \code{\link{movingAverageDF}} for smoothing of time-series prior to periodFinder processing.
#'}
#'
#'
#' @export
#'
periodFinder <- function(data){

  colNames <- colnames(data)

  #get Time Points from Colnames
  xVals <- as.numeric(unlist(regmatches(colNames, gregexpr("[[:digit:]]+\\.*[[:digit:]]*",colNames))))

  #get difference between Time Points
  avgDiff <- mean(diff(xVals))

  dataPeriod <- apply(data, 1, function(TS) {

    #compute the FFT on the Time series
    fft <- getFFT(TS,1/avgDiff)

    #create variable to hold FFT fits with first 3 harmonics
    fit <- list()
    for(i in 1:4){

      #Get period
      per <- (1/fft$freq[i])

      #Get Amplitude
      amp <- Mod(fft$fur)[i]

      #Get Phase
      phase <- Arg(fft$fur)[i]

      #get Cos Wave
      tp <- seq(0,diff(range(xVals)),0.1)
      harmonicFit <- amp*cos(2*pi*(tp/per)+phase)

      #save the harmonic cos fit
      fit[[i]] <- as.vector(harmonicFit)
    }

    #compute Signal using first 3 FFT harmonics
    fftFullfit <- rowSums(do.call(cbind,fit))

    #get local Max and local Min of signal
    localMax <- (which(diff(sign(diff(fftFullfit)))==-2)+1)
    localMin <- (which(diff(sign(diff(fftFullfit)))== 2)+1)
    minMax <- c(localMax,localMin)
    minMaxOrdered <- minMax[order(minMax)]
    deriv <- diff(fftFullfit)

    #search for outlier points identified as min and maxes
    pos <- c(1,minMaxOrdered)
    pat <- rep(seq_along(pos), times=diff(c(pos, length(fftFullfit))))
    z <- split(abs(deriv)/max(abs(deriv)), pat)

    #for each point compute the area under the derivative curve
    areas <- as.vector(unlist(lapply(z,sum)))/sum(abs(deriv)/max(abs(deriv)))

    #if normalized area is less than 0.05 of total, consider outlier
    outliers <- which(areas[2:length(areas)] < 0.05)
    #if normalized area is more than 0.05 of total, keep point
    keep <- which(areas[2:length(areas)] > 0.05)

    for(i in 1:length(keep)){
      if(length(outliers) > 0){
        toAvg <- which(outliers < keep[i]) #which outlier points are less than points to keep
        if(length(toAvg) > 0 ){ #check if there are points to average
          start <- outliers[toAvg][1] #start from outlier point less than keep point
          stop <- keep[i] #upto the first keep point
          selPoints <- start:stop
          minMaxOrdered[keep[i]] <- mean(minMaxOrdered[selPoints]) #average points together
          outliers <- outliers[-toAvg] #remove the points from outlier
        }
      }

    }
    minMaxOrdered <- minMaxOrdered[keep]

    #Get overall period
    if(length(localMax)+length(localMin) < 2){
      #if only one max and min, the period is the sampling length
      per <- abs(diff(range(xVals)))
    } else{
      #else the period is the avg diff between all local maxima and minima
      per <- round(mean(2*abs(diff(minMaxOrdered)))*0.1,2)
    }

    #Get overall Amplitude
    #mean of the absolute value of the mean center signal. taking all local min and max into account
    amp <- mean(abs(fftFullfit[minMaxOrdered] - mean(fftFullfit)))

    #Get overall Phase
    #time to first localMax. ie. how far is the cos wave shifted from zero
    if(length(localMax) < 1){
      phase <-  0
    }else{
      phase <- localMax[1]
    }
    phaseInHours <- (phase*0.1)%%per

    return(round(c(per, amp, phaseInHours),2))
  })

  dataPeriod <- as.data.frame(dataPeriod)
  rownames(dataPeriod) <- c("Period.in.Hours", "Amp", "Phase.in.Hours")
  return(dataPeriod)
}




#' Computes Time-Series Autocovariance Function
#'
#' Computing the autocovariance function from a \code{data.frame} of time-series. The resulting autocovariance
#' function is smoothed using a moving average defined by the period and sampling scheme.
#' Returns the smoothed time-series as a \code{data.frame}.
#'
#' @param data a \code{data.frame} of \code{numeric} gene expression over time (row = genes \emph{x} col = ZT times).
#' @param period a \code{numeric} specifying the period of interest in hours for rhythm detection. Default is \code{24}.
#' @param linearTrend  a \code{logical} scalar. Should TimeCycle Prioritize detecting linear trending signals? Default \code{FALSE}. Not recommended to change from default \code{FALSE} - will increases false positives rate. See vignette("TimeCycle") for more details.
#'
#' @return a smoothed \code{data.frame} of \code{numeric} gene expression covariance over time (row = genes \emph{x} col = ZT times).
#'
#' @seealso \code{\link{meanCenter}}, \code{\link{scaleTimeSeries}}
#'
#' @export
#'
#'
preprocess_acf <- function(data, period = 24, linearTrend = F){
  data <- as.data.frame(data)
  xVals <- as.numeric(colnames(data))
  len <- max(xVals)
  interval <- mean(diff(xVals))
  maxAcfLag <- min(dim(data)[2]-1, 2*period/interval)
  output <- apply(data, 1, function(ts){
    ts <- unlist(unname(ts))
    corr <- stats::acf(as.vector(ts), lag = maxAcfLag, plot = F, type = "covariance")
    corr <- as.vector(corr$acf)
    corr <- scale(corr)

    if(linearTrend){
      # better for Linear Trends, increases Sigmoid false positive
      toCheck <- diff(as.vector(movingAverage(corr, n = period/ceiling(len/period)/interval+1, centered =  T)))
    } else{
      toCheck <- as.vector(movingAverage(corr, n = period/ceiling(len/period)/interval+1, centered =  T))
    }

  })

  output <- as.data.frame(t(output))
  #return the smoothed, mean centered data
  return(output)
}




#' Generates a \code{data.frame} of Resampled Time-Series
#'
#' Converts a \code{data.frame} of \code{numeric} values into a single vector and generates a random resampling with dimension
#' \code{numRsmps} by \code{numTP}.
#'
#' @param distTP a \code{data.frame} of \code{numeric} values.
#' @param numTP a \code{numeric} specifying the number of columns in the outputted \code{data.frame}.
#' @param numRsmps a \code{numeric} specifying the number of rows in the outputted \code{data.frame}.
#'
#' @return a \code{data.frame} of randomly sampled \code{numerics} time-series (row = \code{numRsamps} \emph{x} col = \code{numTP}).
#'
#' @seealso \code{\link{nullResampling}}
#'
#' @export
#'
resampleTimeSeries <- function(distTP, numTP, numRsmps) {
  dist <- as.numeric(as.vector(unlist(distTP)))
  return(as.data.frame(matrix(data = sample(dist, size = numTP * numRsmps, replace = T), ncol = numTP, nrow = numRsmps)))
}




#' Scales a Time-Series
#'
#' Scales a time-series between 0 and 1.
#'
#' @param TimeSeries a \code{vector} of \code{numeric} time-series expression values.
#'
#' @return a scaled \code{vector} of \code{numeric} time-series expression values between 0  and 1.
#'
#' @seealso \code{\link{preprocess_acf}}, \code{\link{meanCenter}}
#'
#' @export
#'
scaleTimeSeries <- function(TimeSeries){
  minVal <- min(TimeSeries)

  if(sum(TimeSeries-mean(TimeSeries) == 0) == length(TimeSeries)){
    output <- TimeSeries
  } else if(minVal < 0){
    output <- TimeSeries + abs(minVal)
  } else{
    output <- TimeSeries - minVal
  }

  maxVal <- max(output)

  if(maxVal == 0){
    output <- rep(0, length(output))
  } else{
    output <- output/maxVal
  }

  return(output)
}
