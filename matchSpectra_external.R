# load required libraries
library(MsBackendMgf)
library(MetaboAnnotation)

# helper function to normalize intensities
norm_int <- function(x, ...) {
  
  maxint <- max(x[, "intensity"], na.rm = TRUE)
  x[, "intensity"] <- 100 * x[, "intensity"] / maxint
  x
  
}

# load example query data and normalize intensities
query_specs <- Spectra("exampleData/fecal_ms2data.mgf",
                       source = MsBackendMgf(),
                       backend = MsBackendDataFrame())

query_specs <- addProcessing(query_specs, norm_int)
query_specs <- applyProcessing(query_specs)

# load external database and normalize intensities
target_external <- Spectra("exampleData/BILELIB19.mgf",
                           source = MsBackendMgf(),
                           backedn = MsBackendDataFrame())

target_external <- addProcessing(target_external, norm_int)
target_external <- applyProcessing(target_external)

# create parameter object
param <- MatchForwardReverseParam(tolerance = 0.005,
                                  ppm = 0.0,
                                  requirePrecursor = TRUE,
                                  requirePrecursorPeak = FALSE,
                                  TRESHFUN = function(x) which(x >= 0.7),
                                  TRESHFUN_REVERSE = function(x) which(x >= 0.7),
                                  toleranceRt = Inf)
  
# perform matching
match_results <- matchSpectra(query_specs,
                              target_external,
                              param)

# retrieve results
matchedData(match_results)

match_results_filtered <- match_results[whichQuery(match_results)]

matchedData(match_results_filtered)

# create mirrorplot
plotSpectraMirror(match_results[whichQuery(match_results)][2])

# browse through results
validateMatchedSpectra(match_results_filtered)
