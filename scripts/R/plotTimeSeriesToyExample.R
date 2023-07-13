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


# Plot the matrix of scatter plots
ggplot(df_new, aes(x = x, y = value, color = variable)) +
  geom_point() +
  labs(x = "X Axis", y = "Y Axis") +
  scale_color_discrete(name = "Y Variables") +
  facet_wrap(~ variable, ncol = 3)  # Adjust the number of columns as needed


plot(df_new[df_new$variable == "a",]$x,df_new[df_new$variable == "a",]$value)
