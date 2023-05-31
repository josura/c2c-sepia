library(ggplot2)

timeSeries.CD4Tconv <- read.csv("/home/josura/Projects/ccc/c2c-sepia/outputsTimeSeries/periodicSatNotfull/CD4Tconv_outputAll.tsv",sep = "\t",header = TRUE)
### visualization all iterations virtual inputs
ggplot() +
  ggtitle("Virtual input in the CD4Tconv augmented metapathway for Fibroblast (log10)") +
  geom_point(data=timeSeries.CD4Tconv, aes(x=iteration, y=v.in.Fibroblasts),size=0.5)+ 
  scale_y_log10()
##without scale 
ggplot() +
  ggtitle("Virtual input in the CD4Tconv augmented metapathway for Fibroblast") +
  geom_point(data=timeSeries.CD4Tconv, aes(x=iteration, y=v.in.Fibroblasts),size=0.5)
ggplot() +
  ggtitle("Virtual input in the CD4Tconv augmented metapathway for Malignant") +
  geom_point(data=timeSeries.CD4Tconv, aes(x=iteration, y=v.in.Malignant),size=0.5)
ggplot() +
  ggtitle("Virtual input in the CD4Tconv augmented metapathway for DC") +
  geom_point(data=timeSeries.CD4Tconv, aes(x=iteration, y=v.in.DC),size=0.5)
### visualization all iterations
ggplot() +
  ggtitle("Virtual input in the CD4Tconv augmented metapathway for Fibroblast") +
  geom_point(data=timeSeries.CD4Tconv[timeSeries.CD4Tconv$iteration < 50,], aes(x=iteration, y=v.in.Fibroblasts),size=0.5)
ggplot() +
  ggtitle("Virtual input in the CD4Tconv augmented metapathway for Malignant") +
  geom_point(data=timeSeries.CD4Tconv[timeSeries.CD4Tconv$iteration < 50,], aes(x=iteration, y=v.in.Malignant),size=0.5)
ggplot() +
  ggtitle("Virtual input in the CD4Tconv augmented metapathway for DC") +
  geom_point(data=timeSeries.CD4Tconv[timeSeries.CD4Tconv$iteration < 50,], aes(x=iteration, y=v.in.DC),size=0.5)
### visualization all iterations
ggplot() +
  ggtitle("Virtual input in the CD4Tconv augmented metapathway for Fibroblast") +
  ### visualization all iterations
ggplot() +
  ggtitle("Virtual input in the CD4Tconv augmented metapathway for Fibroblast") +
  geom_point(data=timeSeries.CD4Tconv[timeSeries.CD4Tconv$iteration < 50,], aes(x=iteration, y=v.in.Fibroblasts),size=0.5)
ggplot() +
  ggtitle("Virtual input in the CD4Tconv augmented metapathway for Malignant") +
  geom_point(data=timeSeries.CD4Tconv[timeSeries.CD4Tconv$iteration < 50,], aes(x=iteration, y=v.in.Malignant),size=0.5)
ggplot() +
  ggtitle("Virtual input in the CD4Tconv augmented metapathway for DC") +
  geom_point(data=timeSeries.CD4Tconv[timeSeries.CD4Tconv$iteration < 50,], aes(x=iteration, y=v.in.DC),size=0.5)

  geom_point(data=timeSeries.CD4Tconv[timeSeries.CD4Tconv$iteration < 30,], aes(x=iteration, y=v.in.Fibroblasts),size=0.5)
ggplot() +
  ggtitle("Virtual input in the CD4Tconv augmented metapathway for Malignant") +
  geom_point(data=timeSeries.CD4Tconv[timeSeries.CD4Tconv$iteration < 30,], aes(x=iteration, y=v.in.Malignant),size=0.5)
ggplot() +
  ggtitle("Virtual input in the CD4Tconv augmented metapathway for DC") +
  geom_point(data=timeSeries.CD4Tconv[timeSeries.CD4Tconv$iteration < 30,], aes(x=iteration, y=v.in.DC),size=0.5)


#virtual outputs

ggplot() +
  ggtitle("Virtual output in the CD4Tconv augmented metapathway for Fibroblast") +
  geom_point(data=timeSeries.CD4Tconv, aes(x=iteration, y=v.out.Fibroblasts),size=0.5)
ggplot() +
  ggtitle("Virtual output in the CD4Tconv augmented metapathway for Malignant") +
  geom_point(data=timeSeries.CD4Tconv, aes(x=iteration, y=v.out.Malignant),size=0.5)
ggplot() +
  ggtitle("Virtual output in the CD4Tconv augmented metapathway for DC") +
  geom_point(data=timeSeries.CD4Tconv, aes(x=iteration, y=v.out.DC),size=0.5)


# Select the first column as x-axis and 20 columns at the end
library(reshape2)
x <- timeSeries.CD4Tconv[, 1]
y <- timeSeries.CD4Tconv[, (ncol(timeSeries.CD4Tconv)-19):ncol(timeSeries.CD4Tconv)]

# Select only the interesting columns
interesting.list <- c("v.in.Fibroblasts","v.in.Malignant","v.in.DC","v.out.Fibroblasts","v.out.Malignant","v.out.DC")
y.int <- timeSeries.CD4Tconv[, interesting.list]

# Combine x and y into a new dataframe
df_new <- data.frame(x = rep(x, ncol(y.int)), y.int = c(y.int))

# Convert y columns into a single column
df_new <- melt(df_new, id.vars = "x")

# Plot the matrix of scatter plots
ggplot(df_new, aes(x = x, y = value, color = variable)) +
  geom_point() +
  labs(x = "X Axis", y = "Y Axis") +
  scale_color_discrete(name = "Y Variables") +
  facet_wrap(~ variable, ncol = 3)  # Adjust the number of columns as needed
