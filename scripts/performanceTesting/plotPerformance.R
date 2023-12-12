library(ggplot2)
library(grid)
library(dplyr)

# Read in the data from tsv files
data.differentNodes <- read.table("timesDifferentNodes.tsv", header=TRUE, sep="\t")
data.differentCommunities <- read.table("timesDifferentCommunities.tsv", header=TRUE, sep="\t")
data.differentProcessors <- read.table("timesDifferentProcessors.tsv", header=TRUE, sep="\t")

# Plot the data
# the data have the following format:
# graphNodesNumber	numberProcesses	numberTypes	numberIterations	time
# every point in the plot is a graph represented by the number of nodes and for every graph we have the average as the center point and the standard deviation as the error bar
# different colors of lines represent different number of iterations
# x axis is the number of nodes
plotTimesDifferentNodes <- function(data, title){
    # prepare the data
    data.prepared <- data %>% 
        group_by(graphNodesNumber, numberIterations) %>% 
        summarise(value = mean(time), sd = sd(time))
    # plot the data
    p <- ggplot(data.prepared, aes(x=graphNodesNumber, y=value, color=numberIterations, group=numberIterations)) + 
        #geom_line() + 
        geom_line() +
        geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.1) +
        ggtitle(title) +
        xlab("number of nodes") +
        ylab("time (milliseconds)") +
        theme(plot.title = element_text(hjust = 0.5))

    return(p)

}


# Plot the data
# the data have the following format:
# graphNodesNumber	numberProcesses	numberTypes	numberIterations	time
# every point in the plot is a graph represented by the number of nodes and for every graph we have the average as the center point and the standard deviation as the error bar
# different colors of lines represent different number of iterations
# x axis name is the number of types (communities)
plotTimesDifferentCommunities <- function(data, title){
    # prepare the data
    data.prepared <- data %>% 
        group_by(numberTypes, numberIterations) %>% 
        summarise(value = mean(time), sd = sd(time))
    # plot the data
    p <- ggplot(data.prepared, aes(x=numberTypes, y=value)) + 
        #geom_line() + 
        geom_line() +
        geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.1) +
        ggtitle(title) +
        xlab("number of communities") +
        ylab("time (milliseconds)") +
        theme(plot.title = element_text(hjust = 0.5))

    return(p)

}


# Plot the data
# the data have the following format:
# graphNodesNumber	numberProcesses	numberTypes	numberIterations	time
# every point in the plot is a graph represented by the number of nodes and for every graph we have the average as the center point and the standard deviation as the error bar
# different colors of lines represent different number of iterations
# x axis name is the number of types (communities)
plotTimesDifferentProcessors <- function(data, title){
  # prepare the data
  data.prepared <- data %>% 
    group_by(numberProcesses, numberIterations) %>% 
    summarise(value = mean(time), sd = sd(time))
  # plot the data
  p <- ggplot(data.prepared, aes(x=numberProcesses, y=value)) + 
    #geom_line() + 
    geom_line() +
    geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.1) +
    ggtitle(title) +
    xlab("number of processors") +
    ylab("time (milliseconds)") +
    theme(plot.title = element_text(hjust = 0.5))
  
  return(p)
  
}

plotTimesDifferentNodes(data.differentNodes, "weak-scaling: different number of nodes")
plotTimesDifferentCommunities(data.differentCommunities, "weak-scaling: different number of communities for the same number of nodes")
plotTimesDifferentProcessors(data.differentProcessors, "strong-scaling: different number of processors for the same graph")
