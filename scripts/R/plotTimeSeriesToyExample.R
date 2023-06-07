library(ggplot2)

timeSeries.t0 <- read.csv("/home/josura/Projects/ccc/c2c-sepia/outputsTimeSeries/t3_outputAll.tsv",sep = "\t",header = TRUE)

library(reshape2)
x <- timeSeries.t0[, 1]
y <- timeSeries.t0[, (ncol(timeSeries.t0)-11):ncol(timeSeries.t0)]

# Select only the interesting columns
#interesting.list <- c("v.in.t3","v.in.t1","v.in.t2","v.out.t3","v.out.t1","v.out.t2")
interesting.list <- c("v.in.t3","v.in.t1","v.out.t3","v.out.t1","v.out.t2")
y.int <- timeSeries.t0[, interesting.list]

# Combine x and y into a new dataframe
df_new <- data.frame(x = rep(x, ncol(y)), y.int = c(y))

# Convert y columns into a single column
df_new <- melt(df_new, id.vars = "x")

# Plot the matrix of scatter plots
ggplot(df_new, aes(x = x, y = value, color = variable)) +
  geom_point() +
  labs(x = "X Axis", y = "Y Axis") +
  scale_color_discrete(name = "Y Variables") +
  facet_wrap(~ variable, ncol = 3)  # Adjust the number of columns as needed
