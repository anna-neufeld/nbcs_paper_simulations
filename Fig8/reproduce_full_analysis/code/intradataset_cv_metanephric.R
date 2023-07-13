library(tidyverse)

metanephric_embedding <- read_tsv("/project2/gilad/jpopp/count-splitting-nb/data/fca_metanephric_full_pcs.tsv")
metanephric_clusters <- read_tsv("/project2/gilad/jpopp/count-splitting-nb/data/fca_metanephric_full_clusters_2e-5.tsv")
stopifnot(all.equal(metanephric_clusters$cell, metanephric_embedding$cell))

X <- as.matrix(metanephric_embedding[,2:ncol(metanephric_embedding)])
clusterLabs <- metanephric_clusters$value
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

saveRDS(confusionMatrix, "/project2/gilad/jpopp/count-splitting-nb/data/fca_metanephric_full_confusionmatrix.rds")
