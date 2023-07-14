library(ggplot2)

timeSeries.t2 <- read.csv("/home/josura/Projects/ccc/c2c-sepia/outputsTimeSeries/t2_outputAll.tsv",sep = "\t",header = TRUE)

library(dplyr)


x <- timeSeries.t2[, 1]
y <- timeSeries.t2[, (ncol(timeSeries.t2)-11):ncol(timeSeries.t2)]

# Select only the interesting columns
#interesting.list <- c("v.in.t3","v.in.t1","v.in.t2","v.out.t3","v.out.t1","v.out.t2")
interesting.list <- c("v.in.t3","v.in.t1","v.out.t3","v.out.t1","v.out.t2")
y.int <- timeSeries.t2[, interesting.list]

# Combine x and y into a new dataframe
df_new <- data.frame(x = rep(x, ncol(y)), y)

# Convert y columns into a single column
df_new <- df_new %>%
  pivot_longer(cols = starts_with("v."), names_to = "variable", values_to = "value")

# Plot the matrix of scatter plots
ggplot(df_new, aes(x = x, y = value, color = variable)) +
  geom_point() +
  labs(x = "X Axis", y = "Y Axis") +
  scale_color_discrete(name = "Y Variables") +
  facet_wrap(~ variable, ncol = 3)  # Adjust the number of columns as needed


# Combine x and y into a new dataframe
df_new <- data.frame(x = rep(x, ncol(y)))
df_new <- cbind(df_new, y)

# Reshape the data by converting columns into a single column
df_new <- data.frame(
  x = rep(df_new$x, times = ncol(y)),
  variable = rep(colnames(y), each = nrow(y)),
  value = c(t(y))
)


# my reshaping since everything else doesn't work
df_new_corrected <- data.frame(x = c(0),variable = c("0"),value = c(0.0))
colnames(df_new_corrected) = c("x","variable","value")


for (colname in colnames(y)) {
  tmp_dataframe <- data.frame(x = x,variable = c(colname),value = timeSeries.t2[,colnames(timeSeries.t2) == colname])
  df_new_corrected <- rbind(df_new_corrected, tmp_dataframe)
}

df_new_corrected <- df_new_corrected[-1,]

# Plot the matrix of scatter plots
ggplot(df_new_corrected, aes(x = x, y = value, color = variable)) +
  geom_point() +
  labs(x = "X Axis", y = "Y Axis") +
  scale_color_discrete(name = "Y Variables") +
  facet_wrap(~ variable, ncol = 3)  # Adjust the number of columns as needed


plot(df_new[df_new$variable == "v.in.t1",]$x,df_new[df_new$variable == "v.in.t1",]$value)
plot( timeSeries.t2$iteration, timeSeries.t2$v.in.t1)
