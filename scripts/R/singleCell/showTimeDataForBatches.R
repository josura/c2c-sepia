library(dplyr)
library(readr)

timeSeriesDirectory <- "/home/josura/Projects/ccc/c2c-sepia/outputsTimeSeries/datiIdoResults_corrected" 

readTimeSeries <- function(timeSeriesDirectory) {
  timeSeries <- list.files(timeSeriesDirectory, full.names = TRUE) %>%
    lapply(read_tsv)
  return(timeSeries)
}

timeSeries <- readTimeSeries(timeSeriesDirectory)
