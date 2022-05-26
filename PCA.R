#input file
input_file = 'https://raw.githubusercontent.com/PineBiotech/omicslogic/master/CellLines_15Genes_1.txt'

#define the data
full_table <- read.table (input_file, header = TRUE, row.names=1, stringsAsFactors=FALSE, check.names=FALSE)
expressions <- as.matrix(full_table[2:nrow(full_table), 2:ncol(full_table)]) 

#transpose the matrix
expressionst <- t(expressions)

#prepare the 'pca' object
pca <- prcomp(expressionst)

#print the pca object:
print(pca)

#plot the pca using native scatterplot function plot()
plot(pca$x, pca$y)


#More "rich" plotting can be done using the autoplot() function. It does not need any object because it creates objects of the ggplot() .

pca <- prcomp(expressionst)

autoplot(pca)

autoplot(pca, label = TRUE, label.size = 3 , colour = "red" )

#The following code will add loading arrows (or correlations between variables and PCs). The loadings (or correlations) allow you to get a sense of the relationships between variables, as well as their associations with the extracted PCs.

autoplot(pca, label = TRUE, label.size = 3 , loadings = TRUE, loadings.label = TRUE, loadings.colour = "blue" )

#Next, to view the plot in 3D, we can use the following code leveraging the scatterplot3d library:
  
library(scatterplot3d)
PCA_table <- as .matrix(pca$x)
scatterplot3d(PCA_table[, 1 : 3 ], angle= 55 )





#now try tthe same with autoplot
library (ggfortify)

#Create plot
plot2 <- autoplot(pca)

#Show plot
plot(plot2)

#print SD of PCA object
print(pca$sd)



#Also you can try scaling and centering data, this is especially useful when your input data is non-normalized before analysis:
  
pca <- prcomp(expressionst, scale. = TRUE, center= TRUE)


#Now that you see visual outputs, let's explore some additional functionality PCA offers to understand how this method works and what information can it display. PCA results are essentially a new table that shows the new position of each sample in a multi-dimensional space created by the analysis. However, the output gives us PC values for genes ('ENSG...'). So ,before we can get the full output, we need to assign an object with sample names and pass it to a new table. Let's produce this results table by adding this code and saving the output table as a txt file:
  
#assign sample names
Samples <- colnames(expressions)

#run PCA and view results
pca <- prcomp(expressionst)
pca_res <- data.frame(Samples, pca$x)


#Save PCA table
write.table(pca_res, file="PCA_table1.txt", row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE, append=TRUE)



#Try running this code by adding this part of the script to what is already in the code block:
  
library (ggfortify)

#input file
input_file = 'https://raw.githubusercontent.com/PineBiotech/omicslogic/master/CellLines_15Genes_1.txt'

#define the data
full_table <- read.table (input_file, header = TRUE, row.names=1, stringsAsFactors=FALSE, check.names=FALSE)
expressions <- as.matrix(full_table[2:nrow(full_table), 2:ncol(full_table)]) 

#assign sample names
Samples <- colnames(expressions)

#run PCA and test length (number of rows)
pca <- prcomp((t(expressions)), scale. = TRUE, center= TRUE)
pca_res <- data.frame (Samples, pca$x)

#plot the pca using native scatterplot function plot()
plot(pca$x, pca$y)

#now try tthe same with autoplot
library (ggfortify)

#Create plot
plot2 <- autoplot(pca)

#Show plot
plot(plot2)

#print SD of PCA object
print(pca$sd)

#Now that you see visual outputs, let's explore some additional functionality PCA offers to understand how this method works and what information can it display. PCA results are essentially a new table that shows the new position of each sample in a multi-dimensional space created by the analysis. However, the output gives us PC values for genes ('ENSG...'). So ,before we can get the full output, we need to assign an object with sample names and pass it to a new table. Let's produce this results table by adding this code and saving the output table as a txt file:
  
#assign sample names
Samples <- colnames(expressions)
#run PCA and view results
pca <- prcomp(expressionst)
pca_res <- data.frame(Samples, pca$x)


#Save PCA table
write.table(pca_res, file="PCA_table1.txt", row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE, append=TRUE)


#print the length of pca_res object to make sure you have 15 rows as expected
print(length())

#add here the script to print the results table:
write.table(....
            
#make the plot from the previous step
plot(....
                 
                 
#These are the principal components that explain variance between samples. Now, let's study these components to understand what information they carry. Add the following to your script:

barplot barplot(pca$x)


#Make better pca plots
autoplot is a generic function to visualize various data object, it tries to give better default graphics and customized choices for each data type, quick and convenient to explore your genomic data compare to low level ggplot method, it is much simpler and easy to produce fairly complicated graphics, though you may loose some flexibility for each layer.

#perform PCA: pca <- prcomp(expressionst) #basic PCA scatter plot: autoplot(pca) #adding labels: autoplot(pca, label = TRUE , label.size = 3 )
But PCA only shows us how to best position the samples in a 2 or 3-dimensional space. Another useful option includes adding clusters - these are boundaries around samples that are most similar. Clustering can be done with the same approach and visualized using a scatterplot:

#clustering
pca <- kmeans(expressionst, 4)
autoplot(pca, data=expressionst, label = TRUE , label.size = 3 )

#When you run this analysis, you can also add a frame around each cluster, for example:
autoplot(pca, data=expressionst, frame= TRUE , label = TRUE , label.size = 3 , frame.type= 'norm' )


#As a result, you will see a border around each cluster, as shown here on the left. The cluster assignment is shown in the legend.

#However, the border for Normal-like and Basal samples here combines them into a single cluster, which is incorrect. That is because K-Means Clustering is relying on the cluster's center point as the 'mean' of that cluster and the others points are the observations that are nearest to the mean value.

#However, in Clustering Large Applications (CLARA) clustering, medoids are used as their center points for cluster, and rest of the observations in a cluster are near to that center point. As a result, CLARA clustering will work on large datasets as compared to K-Medoids and K-Means clustering algorithms, as it will select the random observations from the dataset and perform Partitioning Around Medioids (PAM) algorithm on it.

#Let's see how these methods work and show us which features (genes) are best for separation between clusters. First, we can change clustering from k-means to "clara" using the following method:
                   
scatterplotobject < - clara(expressionst , 4)
autoplot(scatterplotobject, label = TRUE, label.size = 3 , frame = TRUE, frame.type = 'norm' )

#Then, we can also add "loadings" to see how genes are positioned in the PCA space:
                   
scatterplotobject < - clara(expressionst , 4)
autoplot(scatterplotobject,label = TRUE, label.size = 3 , frame = TRUE, frame.type = 'norm' , loadings = TRUE, loadings.label = TRUE,loadings.colour = 'blue' , loadings.label.size = 4 )