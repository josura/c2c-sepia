library(dplyr)
library(readr)


getperturbationTimeSeries <- function(path,output){
  # get a list of all files in the folder
  file_list <- list.files(path)
  
  # initialize an empty list to store the dataframes
  df_list <- list()
  type_list <- list()
  firstFileName <- ""
  # iterate over the files
  for (file_name in file_list) {
    # check if the file ends with ".tsv"
    print(file_name)
    if (endsWith(file_name, ".tsv")) {
      # get the prefix and suffix of the file (e.g. "file1" and "1")
      file_parts <- strsplit(file_name, "--")[[1]]
      file_prefix <- file_parts[1]
      file_suffix <- sub("\\.tsv", "", file_parts[2])
      
      type_list <- c(type_list,file_prefix)
      
    }
    print(file_name)
  }
  
  
  type_list <- unique(type_list)
  for (typeName in type_list){
    node_list <- list()
    firstFileName <- paste(path,paste(typeName,"--0.tsv",sep=""),sep = "")
    node_list <- read.csv(firstFileName,header=TRUE,sep = "\t")$nodeName
    outputFileName <- paste(paste(output,typeName,sep = ""),"_outputAll.tsv",sep = "")
    print(outputFileName)
    # create the file
    #file.create(outputFileName)
    write(paste(c("iteration",node_list),collapse = "\t"), file = outputFileName, append = FALSE, sep = "\t")
  }
  
  for (file_name in file_list) {
    # check if the file ends with ".tsv"
    print(paste("writing values for ",file_name))
    if (endsWith(file_name, ".tsv")) {
      # get the prefix and suffix of the file (e.g. "file1" and "1")
      file_parts <- strsplit(file_name, "--")[[1]]
      file_prefix <- file_parts[1]
      file_suffix <- sub("\\.tsv", "", file_parts[2])
      
      # read the file into a dataframe
      file_path <- paste0(path, file_name)
      outputFileName <- paste(paste(output,file_prefix,sep = ""),"_outputAll.tsv",sep = "")
      df <- read.delim(file_path)
      colValues <- df$nodeValue
      write(paste(c(file_suffix,colValues),collapse = "\t"), file = outputFileName, append = TRUE, sep = "\t")
      # write the results in the corresponding file
      
    }
    print(file_name)
  }
  
}

getAllGraph <- function(pathToSingleFiles,pathToOutput){
  # read the files written in the previous function in a dataframe
  # get a list of all files in the folder
  file_list <- list.files(pathToSingleFiles)
  
  # initialize an empty list to store the dataframes
  df_list <- list()
  
  # iterate over the files
  for (file_name in file_list) {
    # check if the file ends with ".tsv"
    if (endsWith(file_name, "_outputAll.tsv")) {
      # read the file into a dataframe
      file_path <- paste0(pathToSingleFiles, file_name)
      df <- read.delim(file_path, check.names = FALSE)
      # add the dataframe to the list
      df_list[[file_name]] <- df
    }
  }
  
  # filter all the columns in the dataframes that contain in their names v.in or v.out
  # remove the columns from all the dataframes
  for (i in 1:length(df_list)) {
    # get the names of the columns that contain v.in or v.out
    columns_to_remove <- grep("v.in|v.out", names(df_list[[i]]))
    df_list[[i]] <- df_list[[i]][-columns_to_remove]
  }
  
  # join the dataframe to create a single dataframe, join is done on the iteration column
  df <- Reduce(function(x, y) left_join(x, y, by = "iteration"), df_list)
  
  # write the dataframe to the output file
  outputFileNameAll <- paste0(pathToOutput,"fullGraph_output.tsv")
  write_tsv(df, outputFileNameAll, quote = "none", append = FALSE)
}

graphModels <- c("barabasiAlbert","erdosRenyi","barabasiAlbertAging")

nodesNumber <- c(100,1000,10000)

for(graphModel in graphModels){
  pathToPerturbationFolder <- paste0("~/Projects/ccc/c2c-sepia/outputs/",graphModel,"/")
  pathToTimeSeriesFolder <- paste0("~/Projects/ccc/c2c-sepia/outputsTimeSeries/",graphModel,"/")
  for(nodes in nodesNumber){
    for(graphName in range(1,30)){
      # path to output perturbation folder (results of the simulations)
      pathToPerturbationFolder <- paste0(pathToPerturbationFolder,nodes,"Nodes/",graphName,"/")
      #path to output time series folder
      pathToTimeSeriesFolder <- paste0(pathToTimeSeriesFolder,nodes,"Nodes/",graphName,"/")

      # get the subfolders in the folder
      subfolders <- list.dirs(pathToPerturbationFolder, recursive = FALSE)

      # iterate over the subfolders
      for (subfolder in subfolders) {
        # get the name of the subfolder
        subfolder_name <- basename(subfolder)

        # set the path to the folder containing the files
        pathToPerturbation <- paste0(pathToPerturbationFolder, subfolder_name, "/")
        # set the path to the output folder 
        pathToOutput <- paste0(pathToTimeSeriesFolder, subfolder_name, "/")
        pathToSingleFiles <- paste0(pathToOutput, "singleFiles/")

        # create the output folder if it does not exist
        dir.create(pathToOutput, showWarnings = FALSE, recursive = TRUE)
        dir.create(pathToSingleFiles, showWarnings = FALSE, recursive = TRUE)

        # get the perturbation time series
        getperturbationTimeSeries(path = pathToPerturbation,output = pathToSingleFiles)

        # get a single file with all the values like before, but without the virtual nodes
        pathToOutputAll <- paste0(pathToOutput, "allFiles/")
        dir.create(pathToOutputAll, showWarnings = FALSE, recursive = TRUE)

        getAllGraph(pathToSingleFiles,pathToOutputAll)
      } 
    }
  }
}
