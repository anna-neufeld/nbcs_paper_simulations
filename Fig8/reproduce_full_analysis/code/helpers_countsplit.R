#' Takes one matrix of counts and returns a list containing two matrices of counts: a training set and a test set.
#'
#' The training set and the test set are independent under the assumption that the original data follow a Poisson distribution. 
#' Be sure to see the countsplitting tutorial vignette for details of how to use this correctly with existing
#' single cell RNA-seq pipelines.
#'
#' @export
#' @importFrom stats rbinom
#'
#' @param X A cell-by-gene matrix of integer counts
#' @param epsilon The thinning parameter for count splitting. Must be between 0 and 1.
countsplit <- function(X, epsilon=0.5) {
  if (epsilon <= 0 | epsilon >= 1) {
    stop('Parameter epsilon must be in (0,1)')
  }
  
  ## sparse formulation which samples only non-zero entries
  ## and replace sampled value inplace
  if(class(X) == 'dgCMatrix'){
    Xtrain <- X
    Xtrain@x <- as.numeric(rbinom(n=length(X@x), size=X@x, prob=epsilon))
    Xtest <- X-Xtrain
  } else {
    ## dense formulation
    Xtrain <-apply(X,2,function(u) rbinom(n=length(u), size=u, prob=epsilon))
    rownames(Xtrain) <- rownames(X)
    colnames(Xtrain) <- colnames(X)
    Xtest <- X-Xtrain
    rownames(Xtest) <- rownames(X)
    colnames(Xtest) <- colnames(X)
  }
  return(list(train=Xtrain, test=Xtest))
}


#### This is the only important function!!!! 
#### b is the overdispersion parameter. 
betaBinSample <- function(x, b,eps) {
  p <- rbeta(length(x),eps*b,(1-eps)*b)
  return(rbinom(length(x),x,p))
}

#' Takes one matrix of counts and returns a list containing two matrices of counts: a training set and a test set.
#'
#' The training set and the test set are independent under the assumption that the original data follow a Poisson distribution. 
#' Be sure to see the countsplitting tutorial vignette for details of how to use this correctly with existing
#' This function is still in testing mode. Tutorials for this function are coming soon.
#'
#' @importFrom stats rbinom
#'
#' @param X A cell-by-gene matrix of integer counts
#' @param epsilon The thinning parameter for count splitting. Must be between 0 and 1.
#' @param overdisps A vector of length p, where p is the number of columns in X. This vector stores the 
#' gene specific overdispersion parameters. 
nb.countsplit <- function(X, epsilon=0.5, overdisps=NULL) {
  if (epsilon <= 0 | epsilon >= 1) {
    stop('Parameter epsilon must be in (0,1)')
  }
  
  ### If no overdispersion parameters are provided, just call the Poisson count split function.
  if (is.null(overdisps)) {
    return(countsplit(X, epsilon))
  }
  
  if(class(X) == 'dgCMatrix'){
    ##### Overdisp_ij = b_j if X_ij is nonzero, and is 0 otherwise. 
    ## sparse formulation
    Xtrain <- X
    # for only the non-zero entries in X, map the associated (gene-specific) overdispersion param
    mapped_overdisps <- overdisps[which(Xtrain != 0, arr.ind=T)[,"row"]] 
    Xtrain@x <- as.numeric(betaBinSample(Xtrain@x, mapped_overdisps, epsilon))
    Xtrain <- drop0(Xtrain)
    Xtest <- drop0(X-Xtrain)
  } else {
    ## dense formulation
    Xtrain <- t(apply(X, 1, function(u) betaBinSample(u, overdisps, epsilon)))
    Xtest <- X-Xtrain
    rownames(Xtrain) <- rownames(X)
    colnames(Xtrain) <- colnames(X)
    rownames(Xtest) <- rownames(X)
    colnames(Xtest) <- colnames(X)
  }
  return(list(train=Xtrain, test=Xtest))
}

# ##### Generate dirichlet as normalized sum of gammas. 
# dirMulSample <- function(x, folds,b) {
#   epsilon = 1/folds
#   gammas <- rgamma(folds,epsilon*b,1)
#   if (sum(gammas)==0) {
#     gammas[sample(1:length(gammas),1)] <- 1
#   }
#   if (is.infinite(sum(gammas))) {
#     gammas <- rep(1, length(gammas))
#   }
#   ps <- gammas/sum(gammas)
#   x_partition <- rmultinom(1, x, ps)
#   return(x_partition)
# }



# ### Returns a 3D matrix storing the folds of data. 
# multifold.nb.countsplit <- function(X, folds=2, overdisps=NULL) {
#   ### If no overdispersion parameters are provided, just call the Poisson count split function.
#   if (is.null(overdisps)) {
#     stop("Need to write an error function or have this call Poisson multifold. ")
#   }
#   
#   if(class(X) == 'dgCMatrix'){
#     ### ASK JOSH how to make this work for sparse data. 
#   } else {
#     ## dense formulation.
#     ## This is going to be SUPER SLOW. Who can help with this? 
#     partition <- array(NA, dim=c(NROW(X), NCOL(X),folds))
#     for (i in 1:NROW(X)) {
#       for (j in 1:NCOL(X)) {
#         partition[i,j,] <- dirMulSample(X[i,j], folds, overdisps[j])
#       }
#     }
#     return(partition)
#   }
# }
