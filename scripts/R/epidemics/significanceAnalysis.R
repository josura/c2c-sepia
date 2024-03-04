library(ggplot2)

path_to_data_folder <- "path"

data <- read.csv(paste0(path_to_data_folder, "outputAll.tsv"), sep = "\t")

# Create a plot of the number of infected people over time
ggplot(data, aes(x = time, y = infected)) +
  geom_line() +
  xlab("Time") +
  ylab("Number of Infected People") +
  ggtitle("Temporal Evolution of Infected People")
