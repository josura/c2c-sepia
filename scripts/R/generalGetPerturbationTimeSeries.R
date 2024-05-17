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
    write(paste(c("time",node_list),collapse = "\t"), file = outputFileName, append = FALSE, sep = "\t")
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
    if(length(columns_to_remove) != 0){
      # remove the columns from the dataframe
      df_list[[i]] <- df_list[[i]][-columns_to_remove]
    }
  }
  
  # join the dataframe to create a single dataframe, join is done on the iteration column
  df <- Reduce(function(x, y) left_join(x, y, by = "time"), df_list)
  
  # write the dataframe to the output file
  outputFileNameAll <- paste0(pathToOutput,"fullGraph_output.tsv")
  write_tsv(df, outputFileNameAll, quote = "none", append = FALSE)
}

# control if the command line arguments are provided
if (length(commandArgs(trailingOnly = TRUE)) != 2) {
    stop("Please provide the path to the output folder and the path to the output time series folder, usage is: Rscript generalGetPerturbationTimeSeries.R <path_to_outputs_folder> <path_to_outputs_time_series_folder>")
}

# read output folder and the output time series folder from the command line
outputsFolder <- commandArgs(trailingOnly = TRUE)[1]
outputsTimeSeriesFolder <- commandArgs(trailingOnly = TRUE)[2]

# control if the output folder is empty
if (outputsFolder == "" | outputsTimeSeriesFolder == "") {
    stop("Please provide the path to the output folder and the path to the output time series folder")
}



pathToPerturbationFolder <- paste0(outputsFolder,"/")
pathToTimeSeriesFolder <- paste0(outputsTimeSeriesFolder,"/")
# get the subfolders in the folder
subfolders <- list.dirs(pathToPerturbationFolder, recursive = FALSE)
# get only the files in the folder not the subfolders
subFiles <- list.files(pathToPerturbationFolder, pattern=".tsv" ,full.names = TRUE, recursive = FALSE)

# control if the folder is empty or if the folder contains the time series files (no subfolders)
if (length(subfolders) == 0) {
    print("The folder does not contain any subfolders")
    if (length(subFiles) == 0) {
        stop("The folder does not contain any time series files, aborting")
    }
    print("The folder contains time series files, using them")
}

# control if the output folder for the time series exists
if (!dir.exists(pathToTimeSeriesFolder)) {
    stop("The output folder for the time series does not exist, aborting")
}

if(length(subFiles) != 0){
  print("generating time series for the files")
  getperturbationTimeSeries(path = pathToPerturbationFolder,output = pathToTimeSeriesFolder)
  getAllGraph(pathToTimeSeriesFolder,pathToTimeSeriesFolder)
}else {
  if (length(subfolders) != 0){
    print("generating time series for the subfolders")
    # iterate over the subfolders
    for (subfolder in subfolders) {
        # get the name of the subfolder
        subfolder_name <- basename(subfolder)

        ## TESTING
        print(subfolder_name)

        # set the path to the folder containing the files
        pathToPerturbation <- paste0(pathToPerturbationFolder, subfolder_name, "/")
        # set the path to the output folder 
        pathToOutput <- paste0(pathToTimeSeriesFolder, subfolder_name, "/")
        pathToSingleFiles <- paste0(pathToOutput, "singleFiles/")

        # create the output folder if it does not exist
        dir.create(pathToOutput, showWarnings = TRUE, recursive = TRUE)
        dir.create(pathToSingleFiles, showWarnings = TRUE, recursive = TRUE)

        # get the perturbation time series
        getperturbationTimeSeries(path = pathToPerturbation,output = pathToSingleFiles)

        # get a single file with all the values like before, but without the virtual nodes
        pathToOutputAll <- paste0(pathToOutput, "allFiles/")
        dir.create(pathToOutputAll, showWarnings = TRUE, recursive = TRUE)

        getAllGraph(pathToSingleFiles,pathToOutputAll)
    }
  } 

}

    
