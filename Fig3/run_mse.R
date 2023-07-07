#!/usr/local/bin/Rscript

## Run the simulation!

## -----------------------------------------
## load user-defined functions, packages
## -----------------------------------------

## get command line arguments
library("argparse")
## tidy stuff
library("tidyr")
## Needed functions.
source("simulation_functions.R")

## -----------------------------------------
## load any command line arguments
## -----------------------------------------
parser <- ArgumentParser()
parser$add_argument("--simname", default = "test",
                    help = "name of simulation")
args <- parser$parse_args()

overdisps <- c(1,5)
Ks <- c(1,3,5)


clustTry = 1:10

param_grid <- expand.grid(
              K=Ks,
              overdisp=overdisps)
  
 
jobid <- as.numeric(Sys.getenv("SGE_TASK_ID"))
filename <- paste("res/",args$simname, jobid, ".txt", sep="")

for (trial in 1:NROW(param_grid)) {
    current_dynamic_args <- param_grid[trial, ]

    n=1000
    p=1000
    B1 = 1.5

    K <- current_dynamic_args$K
    overdisp <- current_dynamic_args$overdisp
  B0s <- rnorm(p)
  clusters <- sample(1:K, size=n, replace=TRUE)

  logLambda <- matrix(B0s, nrow=n, ncol=p, byrow = TRUE)
  K <- length(unique(clusters))

  c <- 1
  if (K > 1) {
    for (clust in 1:(K-1)) {
      logLambda[clusters==clust,c:(c+p/20-1)] <-  logLambda[clusters==clust,c:(c+p/20-1)]+B1
      c <- c+p/20
    }
  }

  Lambda <- exp(logLambda)
  Lambda_bar_js <- colMeans(Lambda)
  overdisps <- sapply(Lambda_bar_js, function(u) constant_b_function(u, overdisp))

  X <- sapply(1:p, function(u) rnbinom(length(Lambda[,u]), mu=Lambda[,u], size=overdisps[u]))



#### RES_NAIVE                                                                                                                                                                 
  try({
  res_naive <- numClust_naive(X, clustTry, clusters)
  write(t(cbind(res_naive, "naive", n,p,B1,K, overdisp, 1)), file=filename, ncol=12, append=TRUE)
  })

#### RES_SAMPSPLIT                                                                                                                                                             
  #try({                                                                                                                                                                         
  #res_ss <- numClust_sampSplit(X, clustTry, 0.5, clusters)                                                                                                                      
  #write(t(cbind(res_ss, "sampsplit", n,p,B1,K, overdisp, 0.5)), file=filename, ncol=12, append=TRUE)                                                                            
  #})                                                                                                                                                                            

  ### RES_COUNTSPLIT                                                                                                                                                             

  try({
  res_cs_known <- numClust_datathin(X,0.5,clustTry, overdisps, clusters)
  write(t(cbind(res_cs_known, "known.onefold", n,p,B1,K, overdisp, 0.5)), file=filename, ncol=12, append=TRUE)
  })

  try({
  overdisps.sct <- get_sctransform_bs(X)
  res_cs_sct <- numClust_datathin(X,0.5,clustTry,  overdisps.sct, clusters)
  write(t(cbind(res_cs_sct, "sct.onefold", n,p,B1,K, overdisp, 0.5)), file=filename, ncol=12, append=TRUE)
  })

  try({
  res_cs_inf <- numClust_datathin(X,0.5,clustTry, rep(Inf,p), clusters)
  write(t(cbind(res_cs_inf, "inf.onefold", n,p,B1,K, overdisp, 0.5)), file=filename, ncol=12, append=TRUE)
})
}
