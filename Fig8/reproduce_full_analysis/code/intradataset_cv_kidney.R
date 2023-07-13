library(tidyverse)

kidney_embedding <- read_tsv("/project2/gilad/jpopp/count-splitting-nb/data/fca_kidney_full_pcs.tsv")
kidney_clusters <- read_tsv("/project2/gilad/jpopp/count-splitting-nb/data/fca_kidney_full_clusters_1e-6.tsv")
stopifnot(all.equal(kidney_clusters$cell, kidney_embedding$cell))

X <- as.matrix(kidney_embedding[,2:ncol(kidney_embedding)])
clusterLabs <- kidney_clusters$value
folds <- 5

n <- NROW(X)
numClust <- length(unique(clusterLabs))
confusionMatrix <- matrix(0, nrow=numClust, ncol=numClust)

crossValIndices <- sample(1:n, size=n, replace=FALSE)
counter <- 1

for (f in 1:folds) {
  testIndices <-  crossValIndices[counter:(counter+n/folds-1)]
  counter <- counter+n/folds
  Xtrain <- X[-testIndices,]
  Xtest <- X[testIndices,]
  clustersTrain <- clusterLabs[-testIndices]
  clustersTest <- clusterLabs[testIndices]
  
  #### Multiclass SVM.
  classifier <- e1071::svm(x=Xtrain, y=as.factor(clustersTrain), kernel="linear")
  preds <- predict(classifier, newdata=Xtest)
  confusionMatrix <- confusionMatrix + table(preds, clustersTest)
}

saveRDS(confusionMatrix, "/project2/gilad/jpopp/count-splitting-nb/data/fca_kidney_full_confusionmatrix.rds")
