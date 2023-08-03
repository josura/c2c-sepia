library(ggplot2)

timeSeries.t0 <- read.csv("/home/josura/Projects/ccc/c2c-sepia/outputsTimeSeries/t0_outputAll.tsv",sep = "\t",header = TRUE)
timeSeries.t1 <- read.csv("/home/josura/Projects/ccc/c2c-sepia/outputsTimeSeries/t1_outputAll.tsv",sep = "\t",header = TRUE)

library(dplyr)


x.t0 <- timeSeries.t0[, 1]
y.t0 <- timeSeries.t0[, 2:ncol(timeSeries.t0)]

x.t1 <- timeSeries.t1[, 1]
y.t1 <- timeSeries.t1[, 2:ncol(timeSeries.t1)]


# Select only the interesting columns
#interesting.list <- c("v.in.t3","v.in.t1","v.in.t2","v.out.t3","v.out.t1","v.out.t2")
#interesting.list <- c("v.in.t3","v.in.t1","v.out.t3","v.out.t1","v.out.t2")
#y.int <- timeSeries.t2[, interesting.list]



# my reshaping since everything else doesn't work
df_new_corrected.t0 <- data.frame(x = c(0),variable = c("0"),value = c(0.0))
colnames(df_new_corrected.t0) = c("x","variable","value")


for (colname in colnames(y.t0)) {
  tmp_dataframe.t0 <- data.frame(x = x.t0,variable = c(colname),value = timeSeries.t0[,colnames(timeSeries.t0) == colname])
  df_new_corrected.t0 <- rbind(df_new_corrected.t0, tmp_dataframe.t0)
}

df_new_corrected.t0 <- df_new_corrected.t0[-1,]

#t1
df_new_corrected.t1 <- data.frame(x = c(0),variable = c("0"),value = c(0.0))
colnames(df_new_corrected.t1) = c("x","variable","value")


for (colname in colnames(y.t1)) {
  tmp_dataframe.t1 <- data.frame(x = x.t1,variable = c(colname),value = timeSeries.t1[,colnames(timeSeries.t1) == colname])
  df_new_corrected.t1 <- rbind(df_new_corrected.t1, tmp_dataframe.t1)
}

df_new_corrected.t1 <- df_new_corrected.t1[-1,]


# Plot the matrix of scatter plots
ggplot(df_new_corrected.t0, aes(x = x, y = value, color = variable)) +
  geom_point() +
  labs(x = "X Axis", y = "Y Axis") +
  scale_color_discrete(name = "Y Variables") +
  facet_wrap(~ variable, ncol = 3)  # Adjust the number of columns as needed

# Plot the matrix of scatter plots
ggplot(df_new_corrected.t1, aes(x = x, y = value, color = variable)) +
  geom_point() +
  labs(x = "X Axis", y = "Y Axis") +
  scale_color_discrete(name = "Y Variables") +
  facet_wrap(~ variable, ncol = 3)  # Adjust the number of columns as needed


plot(df_new_corrected[df_new_corrected$variable == "a",]$x,df_new_corrected[df_new_corrected$variable == "a",]$value)
plot( timeSeries.t2$iteration, timeSeries.t2$v.in.t1)
