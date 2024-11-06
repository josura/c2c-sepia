library(ggplot2)
library(grid)
library(dplyr)

# Read in the data from tsv files
data.differentNodes <- read.table("timesDifferentNodes.tsv", header=TRUE, sep="\t")
data.differentCommunities <- read.table("timesDifferentCommunities.tsv", header=TRUE, sep="\t")
data.differentProcessors <- read.table("timesDifferentProcessors.tsv", header=TRUE, sep="\t")
data.spaceOccupied <- read.table("spaceOccupied.tsv", header=TRUE, sep="\t")
data.spaceOccupied.finer <- read.table("spaceOccupiedFiner.tsv", header=TRUE, sep="\t")

data.plottedForBMC <- read.table("/home/josura/Projects/ccc/MASFENON/scripts/bash/cluster/performanceAnalysis/localPerformanceTimesMultipleExperiments.tsv", header=TRUE, sep="\t")

data.plottedForBMC.10000Nodes <- data.plottedForBMC[data.plottedForBMC$graphNodesNumber == 10000,]

# Plot the data
# the data have the following format:
# graphNodesNumber	numberProcesses	numberTypes	numberIterations	time
# every point in the plot is a graph represented by the number of nodes and for every graph we have the average as the center point and the standard deviation as the error bar
# different colors of lines represent different number of iterations
# x axis is the number of nodes
plotTimesDifferentNodesDifferentNumberOfIterations <- function(data, title){
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

plotTimesDifferentNodesDifferentProcessors <- function(data, title){
  # prepare the data
  data.prepared <- data %>% 
    group_by(graphNodesNumber, numberProcesses) %>% 
    summarise(value = mean(time), sd = sd(time))
  # plot the data
  p <- ggplot(data.prepared, aes(x=graphNodesNumber, y=value, color=numberProcesses, group=numberProcesses)) + 
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



plotTimesDifferentNodesDifferentNumberOfIterations(data.differentNodes, "weak-scaling: different number of nodes")
plotTimesDifferentCommunities(data.differentCommunities, "weak-scaling: different number of communities for the same number of nodes(10000)")
plotTimesDifferentProcessors(data.differentProcessors, "strong-scaling: different number of processors for the same graph(30000 nodes, 113 communities)")

# BMC performance
plotTimesDifferentNodesDifferentProcessors(data.plottedForBMC,"weak-scaling: different number of nodes")
plotTimesDifferentProcessors(data.plottedForBMC.10000Nodes,"strong-scaling: different number of processors for one of the 10000 nodes Erdos-Renyi graph")

# plot speedup for strong scaling
# the data is the same as for the plotTimesDifferentProcessors, but we need to calculate the speedup
# speedup = avg_time(1 processor) / avg_time(n processors)
plotSpeedup <- function(data, title){
  # prepare the data
  data.prepared <- data %>% 
    group_by(numberProcesses, numberIterations) %>% 
    summarise(value = mean(time), sd = sd(time))
  # calculate the speedup
  data.prepared$speedup <- data.prepared$value[1] / data.prepared$value
  # plot the data
  p <- ggplot(data.prepared, aes(x=numberProcesses, y=speedup)) + 
    #geom_line() + 
    geom_line() +
    ggtitle(title) +
    xlab("number of processors") +
    ylab("speedup") +
    theme(plot.title = element_text(hjust = 0.5))
  
  return(p)
}

plotSpeedup(data.differentProcessors, "strong-scaling: average speedup for different number of processors for the same graph")

# plot space occupied
# the data have the following format:
# nodes max_memory max_order	number_of_communities	interactions_inter_community

plotSpaceOccupied <- function(data, title){
  # plot the data
  p <- ggplot(data, aes(x=nodes, y=max_memory)) + 
    geom_line() +
    ggtitle(title) +
    xlab("number of nodes") +
    ylab("space occupied (kB)") +
    theme(plot.title = element_text(hjust = 0.5))
  
  # fit the data to a sublinear function, because the space occupied should not grow linearly with the number of nodes
  p <- p + geom_smooth(method="lm", formula=y~x*log(x), se=FALSE)
  # p <- p + geom_smooth(method="lm", formula=y~(x*x)*log(x), color = "red")

  # add another line based on the formula (max_order + interactions_inter_community/number_of_communities)^2 *number_of_communities
  # p <- p + geom_line(aes(x=nodes, y=(nodes*nodes) *log(nodes)), color="yellow")

  
  return(p)
}

plotSpaceOccupied(data.spaceOccupied.finer, "space occupied by the graph representation")

# estimate the non linear regression function that estimates the space occupied by the finer data

# fit the data to a sublinear function, because the space occupied should not grow linearly with the number of nodes
fit <- lm(max_memory ~ nodes*log(nodes), data=data.spaceOccupied.finer)
summary(fit)

fit.2 <- lm(max_memory ~ (max_order + interactions_inter_community/number_of_communities)^2 *number_of_communities, data = data.spaceOccupied.finer)
summary(fit.2)

# plot the function estimated, add a legend for every line
funcEstimate <- function(x){
  return(100000+0.00005*(x^2)*log(25*x))
}

funcEstimate.2 <- function(x){
  return(100000+0.00000005*(x^3))
}

spaceEstimate <- function(max_order, interactions_inter_community,number_of_communities){
  return(83000 + 0.011*(max_order + interactions_inter_community/number_of_communities)^2 *number_of_communities)
}

x <- seq(1000, 30000, 1000)
y <- funcEstimate(x)

y.spaceEstimate <- spaceEstimate(data.spaceOccupied.finer$max_order, data.spaceOccupied.finer$interactions_inter_community, data.spaceOccupied.finer$number_of_communities)

p <- ggplot(data.spaceOccupied.finer, aes(x=nodes, y=max_memory)) + 
  geom_line(aes(color = "Actual Memory")) + # Assigning a color aesthetic for the legend
  ggtitle("space occupied by the graph representation") +
  xlab("number of nodes in the network of networks") +
  ylab("space occupied (kB)") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_line(aes(x=x, y=y, color = "nodes^2*log(nodes)")) + # Assigning a color aesthetic for the legend
  geom_line(aes(x=nodes, y=y.spaceEstimate, color = "Estimated Memory")) + # Assigning a color aesthetic for the legend
  scale_color_manual(name = "Legend", values = c("Actual Memory" = "black", "nodes^2*log(nodes)" = "red", "Estimated Memory" = "blue")) # Manual color and label assignment for the legend

p

