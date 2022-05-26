#Load packages
library(Rtsne)
library(ggplot2)

#Set seed
set.seed(2)

#load data
data<- read.table('https://raw.githubusercontent.com/PineBiotech/bioinformatics/master/15gene_transposed.txt', sep='\t', header=TRUE, check.names=FALSE, stringsAsFactors=TRUE, row.names=1)

#Take data without groups
df <- data[2:16]

#Create t-SNE object
tsne_realData <- Rtsne(df, perplexity=15, verbose=TRUE, max_iter = 550, check_duplicates = FALSE)

#Add colors
colors = rainbow(length(unique(data$Group)))

#Create colored plot
plot(tsne_realData$Y,  col=colors)

#Check dimensions of t-SNE input data
dim(df)


Code for colored t-SNE plot with Labels

#Load packages
library(Rtsne)
library(ggplot2)

#Set seed
set.seed(2)

#load data
data<- read.table('https://raw.githubusercontent.com/PineBiotech/bioinformatics/master/15gene_transposed.txt', sep='\t', header=TRUE, check.names=FALSE, stringsAsFactors=TRUE, row.names=1)

#Take data without groups
df <- data[2:16]

#Create t-SNE object
tsne_realData <- Rtsne(df, perplexity=15, verbose=TRUE, max_iter = 550, check_duplicates = FALSE)

#Add colors
colors = rainbow(length(unique(data$Group)))

#Make tSNE object as data frame
tsne_df <- as.data.frame(tsne_realData$Y)

#Extract labels of samples
lab <- data$Group

#Combine samples with t-SNE plot
tsne_df1 <- cbind(tsne_df,lab)

#Add column names
colnames(tsne_df1) <- c("X","Y","Labels")
#head(tsne_df1)

#Check dimensions of t-SNE input data
dim(tsne_df1)

#Create colored plot with labels of sample groups
ggplot(tsne_df1, aes(X, Y, colour = Labels)) +geom_point()




