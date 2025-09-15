# load required libraries
library(MsBackendMgf)
library(SpectraQL)

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

# search for LPC species
query_results_spec <- query(query_specs, "QUERY * WHERE MS2PROD=184.07:TOLERANCEMZ=0.005")
query_results_df <- query(query_specs, "QUERY scaninfo(MS2DATA) WHERE MS2PROD=184.07:TOLERANCEMZ=0.005")

# plot spectra
plotSpectra(query_results_spec)

# plot LPC features
plot(query_results_df$rtime,
     query_results_df$precursorMz)

