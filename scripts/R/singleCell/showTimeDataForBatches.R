library(dplyr)
library(readr)

timeSeriesDirectory <- "/home/josura/Projects/ccc/c2c-sepia/outputsTimeSeries/risultatiIdo_noSat" 

readTimeSeries <- function(timeSeriesDirectory) {
  files <- list.files(timeSeriesDirectory, full.names = TRUE)
  #remove the fullgraph_output.tsv file
  files <- files[!grepl("fullGraph_output.tsv", files)]
  timeSeries <- files  %>%
    lapply(read_tsv)
  # remove .x from the column names, when present
  for (i in 1:length(timeSeries)) {
    colnames(timeSeries[[i]]) <- gsub("\\.x", "", colnames(timeSeries[[i]]))
        }
  # get the common columns of all the data frames in the list
  commonCols <- Reduce(intersect, lapply(timeSeries, colnames))
  # select only the common columns for eevery data frame in the list
  for (i in 1:length(timeSeries)) {
    timeSeries[[i]] <- timeSeries[[i]][, commonCols]
  }
  # add the name of the file without the path and the extension to every data frame in the list
  typeNames <- sapply(files, function(x) {
    basename(x) %>%
      gsub("\\_outputAll.tsv", "", .)
  })
  # remove names from the typeNames vector
  names(typeNames) <- NULL
  for (i in 1:length(timeSeries)) {
    timeSeries[[i]]$type <- typeNames[i]
  } 
  return(timeSeries)
}

timeSeries <- readTimeSeries(timeSeriesDirectory)



timeSeries_all <- timeSeries[[1]]
for (i in 2:length(timeSeries)) {
  timeSeries_all <- rbind(timeSeries_all, timeSeries[[i]])
}

features_to_visualize <- c("Kras","Trp53","Egfr","Rb1","Nf1","Myc","Braf")

install.packages("gridExtra")
library(gridExtra)

plotList <- list()
for (feature in features_to_visualize) {
  p <- ggplot(timeSeries_all, aes_string(x = "time", y = feature, color = "type")) +
    geom_point() +
    geom_line() +
    labs(title = feature) +
    theme_minimal()
  plotList[[feature]] <- p
}   

do.call(grid.arrange, c(plotList, ncol = 3))
