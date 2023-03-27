

# set the path to the folder containing the files
pathToPerturbation <- "path/to/folder/"


getperturbationTimeSeries <- function(path){
  # get a list of all files in the folder
  file_list <- list.files(path)
  
  # initialize an empty list to store the dataframes
  df_list <- list()
  cell_list <- list()
  
  # iterate over the files
  for (file_name in file_list) {
    # check if the file ends with ".tsv"
    if (endsWith(file_name, ".tsv")) {
      # get the prefix and suffix of the file (e.g. "file1" and "1")
      file_parts <- strsplit(file_name, "--")[[1]]
      file_prefix <- file_parts[1]
      file_suffix <- sub("\\.tsv", "", file_parts[2])
      
      cell_list <- c(file_prefix)
      
      # read the file into a dataframe
      file_path <- paste0(path, file_name)
      df <- read.delim(file_path)
      
      # write the results in the corresponding file
      
    }
  }  
}

