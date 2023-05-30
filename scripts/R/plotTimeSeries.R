library(ggplot2)

timeSeries.CD4Tconv <- read.csv("/home/josura/Projects/ccc/c2c-sepia/outputsTimeSeries/CD4Tconv_outputAll.tsv",sep = "\t",header = TRUE)
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
