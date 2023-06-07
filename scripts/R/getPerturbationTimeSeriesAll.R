

# set the path to the folder containing the files
pathToPerturbation <- "~/Projects/ccc/c2c-sepia/outputs/"
# set the path to the output folder
pathToOutput <- "~/Projects/ccc/c2c-sepia/outputsTimeSeries/"
library(dplyr)

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
getperturbationTimeSeries(path = pathToPerturbation,output = pathToOutput)
