constant_b_function <- function(Lambda_bar_j, constant=3) {
  return(Lambda_bar_j/constant)
}

get_sctransform_bs <- function(X) {
  rownames(X) <- 1:NROW(X)
  colnames(X) <- 1:NCOL(X)
  fit <- sctransform::vst(t(X), verbosity=0, min_cells=0)
  b_j1 <- fit$model_pars_fit[,1]
  if (length(b_j1) != NCOL(X)) {
    temp <- b_j1
    b_j1 <- rep(Inf, NCOL(X))
    b_j1[rownames(temp)] <- temp
  }
  return(b_j1)
}

mynbglm <- function(u, lhat, sizeFacs) {
  if (length(unique(u))==1) {
    return(c(NA, 0, NA))
  }
  
  if (is.null(sizeFacs)) {
    try1 <- try(suppressWarnings(MASS::glm.nb(u~lhat)))
  } else {
    try1 <- try(suppressWarnings(MASS::glm.nb(u~lhat+offset(log(sizeFacs)))))
  }
  if (class(try1)[1]=="try-error") {
    return(c(NA, NA, NA))
  } else {
    return(c(summary(try1)$coefficients[2,4], summary(try1)$coefficients[2,1], try1$theta))
  }
}

betaBinSample <- function(x, b,eps) {
  p <- rbeta(length(x),eps*b,(1-eps)*b)
  return(rbinom(length(x),x,p))
}

constant_b_function <- function(Lambda_bar_j, constant=3) {
  return(Lambda_bar_j/constant)
}

# Feeling better about this one now. 
cluster.sse.log <- function(trainDat, testDat, clusters.train, clusters.test,
                            eps.train, eps.test) {
  totSS <- 0
  for (lab in unique(clusters.train)) {
    clustdat.test <- testDat[clusters.test==lab,, drop='F']
    clustdat.train <- trainDat[clusters.train==lab,, drop='F']
    
    #### This is done on the scale of the original data X. 
    pred.means <- 1/eps.train*colMeans(clustdat.train)
    ss <- apply(1/eps.test*clustdat.test, 1,  function(u) sum((log(u+1)-log(pred.means+1))^2))
    totSS <- totSS+sum(ss)
  }
  return(totSS)
}


#### Am I 100% sure that this works??? I should really test it haha.
#### But I think I did test!
dirMulSample <- function(x, folds,b, epsilon=1/folds) {
  gammas <- rgamma(folds,epsilon*b,1)
  if (sum(gammas)==0) {
    gammas[sample(1:length(gammas),1)] <- 1
  }
  if (is.infinite(sum(gammas))) {
    gammas <- rep(1, length(gammas))
  }
  ps <- gammas/sum(gammas)
  x_partition <- rmultinom(1, x, ps)
  return(x_partition)
}

# Always two clusters for simplicty. 
diffExp_naive <- function(X, clusters) {
  
  #clusters <- sample(1:2, size=n, replace=TRUE)
  #X <- genData(n,p,clusters,B0s, B1, overdisp)
  
  
  clusters.full<- kmeans(log(X+1), centers=2, nstart=10)$cluster
  rands.full <- mclust::adjustedRandIndex(clusters.full, clusters)
  res <- t(apply(X, 2, function(u) mynbglm(u, clusters.full, sizeFacs=NULL)))

  return(cbind(res, rands.full))
}

diffExp_sampsplit <- function(X, clusters, eps.train=0.5) {
  train <- sample(1:NROW(X), size=round(eps.train*n))
  Xtrain <- X[train,]
  Xtest <- X[-train,]
  
  clusters.train <- kmeans(log(Xtrain+1), centers=2, nstart=10)$cluster
  rands.train <- mclust::adjustedRandIndex(clusters.train, clusters[train])
  
  clusters.test <- class::knn(log(Xtrain+1), log(Xtest+1), clusters.train, 3)
  res <- t(apply(Xtest, 2, function(u) mynbglm(u, clusters.test, sizeFacs=NULL)))
  
  return(cbind(res, rands.train))
}
  

diffExp_datathin <- function(X, clusters, overdisps.hat, eps.train=0.5) {
  Xtrain <- t(apply(X, 1, function(u) betaBinSample(u, overdisps.hat, eps.train)))
  Xtest <- X-Xtrain
  clusters.train <- kmeans(log(Xtrain+1), centers=2, nstart=10)$cluster
  rands.train <- mclust::adjustedRandIndex(clusters.train, clusters)
  res <- t(apply(Xtest, 2, function(u) mynbglm(u, clusters.train, sizeFacs=NULL)))

  return(cbind(res, rands.train))
}

  
  
  
 

numClust_naive <- function(X, clustTry, clusters) {
  
  rands <- rep(NA, length(clustTry))
  res <- rep(NA, length(clustTry))
    
  for (j in 1:length(clustTry)) {
    ktry <- clustTry[j]
    clusters.full <- kmeans(log(X+1), centers=ktry)$cluster
    rands[j] <- mclust::adjustedRandIndex(clusters.full, clusters)
    res[j] <- cluster.sse.log(X, X, clusters.full, clusters.full, 1,1)
  }
  kGuess <- clustTry[which.min(res)]
  randBest <- which.max(rands)
  
  return(cbind(res/(n*p), rands, clustTry, randBest, kGuess))
} 

numClust_sampSplit <- function(X, clustTry, eps.train=0.5,clusters) {
  
  train <- sample(1:NROW(X), size=floor(eps.train*n))
  Xtrain <- X[train,]
  Xtest <- X[-train,]
  
  rands <- rep(NA, length(clustTry))
  res <- rep(NA, length(clustTry))
  
  for (j in 1:length(clustTry)) {
    ktry <- clustTry[j]
    clusters.train <- kmeans(log(Xtrain+1), centers=ktry)$cluster
    clusters.test <- class::knn(log(Xtrain+1), log(Xtest+1), clusters.train, 3)
    
    rands[j] <- mclust::adjustedRandIndex(clusters.train, clusters[train])
    res[j] <- cluster.sse.log(Xtrain, Xtest, clusters.train, clusters.test,
                              eps.train, 1-eps.train)
  }
  kGuess <- clustTry[which.min(res)]
  randBest <- which.max(rands)
  
  return(cbind(res/(n*p), rands, clustTry, randBest, kGuess))
}

numClust_crossVal <- function(X, folds, clustTry, overdisps.hat,clusters) {
  
  partition <- array(NA, dim=c(NROW(X), NCOL(X),folds))
  for (i in 1:NROW(X)) {
    for (j in 1:NCOL(X)) {
      partition[i,j,] <- dirMulSample(X[i,j], folds, overdisps.hat[j])
    }
  }
  
  fullRes <- matrix(NA, nrow=folds, ncol=length(clustTry))
  rands.trains <- matrix(NA, nrow=folds, ncol=length(clustTry))
  
  for (fold in 1:folds) {
  
    trainDat <- apply(partition[,,-c(fold)], 1:2, sum)
    testDat <- partition[,,fold]
    
    for (j in 1:length(clustTry)) {
      ktry <- clustTry[j]
      clusters.train <- kmeans(log(trainDat+1), centers=ktry, nstart=10)$cluster
      rands.trains[fold, j] <- mclust::adjustedRandIndex(clusters.train, clusters)
      fullRes[fold,j] <- cluster.sse.log(trainDat, testDat, clusters.train, clusters.train, 1-1/folds, 1/folds)
    }
  }

  #true.sse <- cluster.sse.log(trainDat, testDat, clusters, clusters, eps.train, 1-eps.train)
  kGuess <- clustTry[which.min(colSums(fullRes))]
  randBest <- which.max(colMeans(rands.trains))
  
  return(cbind(colSums(fullRes)/(n*p*folds), colMeans(rands.trains), clustTry, randBest, kGuess))
}

numClust_datathin <- function(X, eps.train, clustTry, overdisps.hat,clusters) {
  
  Xtrain <- t(apply(X, 1, function(u) betaBinSample(u, overdisps.hat, eps.train)))
  Xtest <- X-Xtrain
  
  res <- rep(NA, length(clustTry))
  rands <- rep(NA, length(clustTry))
  for (j in 1:length(clustTry)) {
      ktry <- clustTry[j]
      clusters.train <- kmeans(log(Xtrain+1), centers=ktry, nstart=10)$cluster
      rands[j] <- mclust::adjustedRandIndex(clusters.train, clusters)
      res[j] <- cluster.sse.log(Xtrain, Xtest, clusters.train, clusters.train, eps.train, 1-eps.train)
  }

  true.sse <- cluster.sse.log(Xtrain, Xtest, clusters, clusters, eps.train, 1-eps.train)
  kGuess <- clustTry[which.min(res)]
  randBest <- which.max(rands)
  
  return(cbind(res/(n*p), rands, clustTry, randBest, kGuess))
}



















