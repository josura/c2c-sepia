library(dplyr)
library(readr)


getperturbationTimeSeries <- function(path,output){
  # get a list of all files in the folder
  file_list <- list.files(path)
  
  # initialize an empty list to store the dataframes
  df_list <- list()
  cell_list <- list()
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
      
      cell_list <- c(cell_list,file_prefix)
      
    }
    print(file_name)
  }
  
  
  cell_list <- unique(cell_list)
  for (cellName in cell_list){
    gene_list <- list()
    firstFileName <- paste(path,paste(cellName,"--0.tsv",sep=""),sep = "")
    gene_list <- read.csv(firstFileName,header=TRUE,sep = "\t")$nodeName
    outputFileName <- paste(paste(output,cellName,sep = ""),"_outputAll.tsv",sep = "")
    print(outputFileName)
    # create the file
    #file.create(outputFileName)
    write(paste(c("iteration",gene_list),collapse = "\t"), file = outputFileName, append = FALSE, sep = "\t")
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


#100 nodes

# set the path to the folder containing the files
pathToPerturbationFolder.100Nodes <- "~/Projects/ccc/c2c-sepia/outputs/100NodesEpidemics/"
# set the path to the output folder 
pathToOutputFolder.100Nodes <- "~/Projects/ccc/c2c-sepia/outputsTimeSeries/epidemics100Nodes/"

# get the subfolders in the folder
subfolders <- list.dirs(pathToPerturbationFolder.100Nodes, recursive = FALSE)

# iterate over the subfolders
for (subfolder in subfolders) {
  # get the name of the subfolder
  subfolder_name <- basename(subfolder)
  
  # set the path to the folder containing the files
  pathToPerturbation.100Nodes <- paste0(pathToPerturbationFolder.100Nodes, subfolder_name, "/")
  # set the path to the output folder 
  pathToOutput.100Nodes <- paste0(pathToOutputFolder.100Nodes, subfolder_name, "/")
  pathToSingleFiles.100Nodes <- paste0(pathToOutput.100Nodes, "singleFiles/")
  
  # create the output folder if it does not exist
  dir.create(pathToOutput.100Nodes, showWarnings = FALSE, recursive = TRUE)
  dir.create(pathToSingleFiles.100Nodes, showWarnings = FALSE, recursive = TRUE)
  
  # get the perturbation time series
  getperturbationTimeSeries(path = pathToPerturbation.100Nodes,output = pathToSingleFiles.100Nodes)
  
  # get a single file with all the values like before, but without the virtual nodes
  pathToOutputAll.100Nodes <- paste0(pathToOutput.100Nodes, "allFiles/")
  dir.create(pathToOutputAll.100Nodes, showWarnings = FALSE, recursive = TRUE)
  
  getAllGraph(pathToSingleFiles.100Nodes,pathToOutputAll.100Nodes)
}

#1000 nodes

# set the path to the folder containing the files
pathToPerturbationFolder.1000Nodes <- "~/Projects/ccc/c2c-sepia/outputs/1000NodesEpidemics/"
# set the path to the output folder 
pathToOutputFolder.1000Nodes <- "~/Projects/ccc/c2c-sepia/outputsTimeSeries/epidemics1000Nodes/"

# get the subfolders in the folder
subfolders <- list.dirs(pathToPerturbationFolder.1000Nodes, recursive = FALSE)

# iterate over the subfolders
for (subfolder in subfolders) {
  # get the name of the subfolder
  subfolder_name <- basename(subfolder)
  
  # set the path to the folder containing the files
  pathToPerturbation.1000Nodes <- paste0(pathToPerturbationFolder.1000Nodes, subfolder_name, "/")
  # set the path to the output folder 
  pathToOutput.1000Nodes <- paste0(pathToOutputFolder.1000Nodes, subfolder_name, "/")
  pathToSingleFiles.1000Nodes <- paste0(pathToOutput.1000Nodes, "singleFiles/")
  
  # create the output folder if it does not exist
  dir.create(pathToOutput.1000Nodes, showWarnings = FALSE, recursive = TRUE)
  dir.create(pathToSingleFiles.1000Nodes, showWarnings = FALSE, recursive = TRUE)
  
  # get the perturbation time series
  getperturbationTimeSeries(path = pathToPerturbation.1000Nodes,output = pathToSingleFiles.1000Nodes)
  
  # get a single file with all the values like before, but without the virtual nodes
  pathToOutputAll.1000Nodes <- paste0(pathToOutput.1000Nodes, "allFiles/")
  dir.create(pathToOutputAll.1000Nodes, showWarnings = FALSE, recursive = TRUE)
  
  getAllGraph(pathToSingleFiles.1000Nodes,pathToOutputAll.1000Nodes)
}

#10000 nodes

# set the path to the folder containing the files
pathToPerturbationFolder.10000Nodes <- "~/Projects/ccc/c2c-sepia/outputs/10000NodesEpidemics/"
# set the path to the output folder 
pathToOutputFolder.10000Nodes <- "~/Projects/ccc/c2c-sepia/outputsTimeSeries/epidemics10000Nodes/"

# get the subfolders in the folder
subfolders <- list.dirs(pathToPerturbationFolder.10000Nodes, recursive = FALSE)

# iterate over the subfolders
for (subfolder in subfolders) {
  # get the name of the subfolder
  subfolder_name <- basename(subfolder)
  
  # set the path to the folder containing the files
  pathToPerturbation.10000Nodes <- paste0(pathToPerturbationFolder.10000Nodes, subfolder_name, "/")
  # set the path to the output folder 
  pathToOutput.10000Nodes <- paste0(pathToOutputFolder.10000Nodes, subfolder_name, "/")
  pathToSingleFiles.10000Nodes <- paste0(pathToOutput.10000Nodes, "singleFiles/")
  
  # create the output folder if it does not exist
  dir.create(pathToOutput.10000Nodes, showWarnings = FALSE, recursive = TRUE)
  dir.create(pathToSingleFiles.10000Nodes, showWarnings = FALSE, recursive = TRUE)
  
  # get the perturbation time series
  getperturbationTimeSeries(path = pathToPerturbation.10000Nodes,output = pathToSingleFiles.10000Nodes)
  
  # get a single file with all the values like before, but without the virtual nodes
  pathToOutputAll.10000Nodes <- paste0(pathToOutput.10000Nodes, "allFiles/")
  dir.create(pathToOutputAll.10000Nodes, showWarnings = FALSE, recursive = TRUE)
  
  getAllGraph(pathToSingleFiles.10000Nodes,pathToOutputAll.10000Nodes)
}
