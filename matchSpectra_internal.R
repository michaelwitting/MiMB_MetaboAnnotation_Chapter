# load required libraries
library(MsBackendMgf)
library(MsBackendMassbank)
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

# load internal database and normalize intensities
target_internal <- Spectra(list.files("exampleData/in-house_MassBank/", full.names = TRUE),
                           source = MsBackendMassbank(),
                           backedn = MsBackendDataFrame())

target_internal <- addProcessing(target_internal, norm_int)
target_internal <- applyProcessing(target_internal)

# create parameter object
param <- MatchForwardReverseParam(tolerance = 0.005,
                                  ppm = 0.0,
                                  requirePrecursor = TRUE,
                                  requirePrecursorPeak = FALSE,
                                  TRESHFUN = function(x) which(x >= 0.7),
                                  TRESHFUN_REVERSE = function(x) which(x >= 0.7),
                                  toleranceRt = 0.5 * 60)

# perform matching
match_results <- matchSpectra(query_specs,
                              target_internal,
                              param)

# retrieve results
matchedData(match_results)

match_results_filtered <- match_results[whichQuery(match_results)]

matchedData(match_results_filtered)

# create mirrorplot
plotSpectraMirror(match_results[whichQuery(match_results)][1])

# browse through results
validateMatchedSpectra(match_results_filtered)
