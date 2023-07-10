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
#parser$add_argument("--nreps", type = "double", default = 10,
#                    help = "number of replicates for each set of params")
args <- parser$parse_args()


## -----------------------------------------
## set up a grid of parameters to cycle over
## -----------------------------------------
ns <- c(500)
ps <- c(40)
B1s <- c(2,3,4,5,6)
overdisps <- c(1)
Ks <- c(5)
### 

#nreps_per_combo <- args$nreps

param_grid <- expand.grid(n=ns,
              p=ps,
              K=Ks,
              B1=B1s,
              overdisp=overdisps)


 
jobid <- as.numeric(Sys.getenv("SGE_TASK_ID"))

#sizeGrid <- NROW(param_grid)
#current_dynamic_args <- param_grid[(jobid%%sizeGrid)+1, ]


## -----------------------------------------
## run the simulation nreps_per_job times
## -----------------------------------------
#set.seed(1)
#B0s <- rnorm(current_dynamic_args$p)
#set.seed(jobid)

#B0s <- rnorm(current_dynamic_args$p)     

filename <- paste("res/",args$simname, jobid, ".txt", sep="")

#filename <- "testfig1_2.text"


for (trial in 1:NROW(param_grid)) {
  current_dynamic_args <- param_grid[trial, ]                                                                                                                                 

  n=current_dynamic_args$n
  p=current_dynamic_args$p
                                        #B0s=B0s
B1 = current_dynamic_args$B1
K = current_dynamic_args$K
clustTry = 1:10
overdisp=current_dynamic_args$overdisp
 B0s <- rnorm(current_dynamic_args$p) 
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
  

  epsilons <- seq(0,1,length.out=41)[2:40]
  for (eps in epsilons) {  ### RES_COUNTSPLIT 50/50
  try({
  res_cs_known <- numClust_datathin(X,eps,clustTry, overdisps, clusters)
  write(t(cbind(res_cs_known, "known.onefold", n,p,B1,K, overdisp, eps)), file=filename, ncol=12, append=TRUE)
  })
  }
  
  #try({
  #overdisps.sct <- get_sctransform_bs(X)
  #res_cs_sct <- numClust_datathin(X,0.5,clustTry,  overdisps.sct, clusters)
  #write(t(cbind(res_cs_sct, "sct.onefold", n,p,B1,K, overdisp, 0.5)), file=filename, ncol=12, append=TRUE)
  #})
  
  ### RES_COUNTSPLIT 90/10
  #try({
  #res_cs_known <- numClust_datathin(X,0.9,clustTry, overdisps, clusters)
  #write(t(cbind(res_cs_known, "known.onefold", n,p,B1,K, overdisp, 0.9)), file=filename, ncol=12, append=TRUE)
  #})
  
  #try({
  #overdisps.sct <- get_sctransform_bs(X)
  #res_cs_sct <- numClust_datathin(X,0.9,clustTry,  overdisps.sct, clusters)
  #write(t(cbind(res_cs_sct, "sct.onefold", n,p,B1,K, overdisp, 0.9)), file=filename, ncol=12, append=TRUE)
  #})
  
  
  ### RES_CROSSVAL
  #try({
  #res_cross_known <- numClust_crossVal(X,10,clustTry, overdisps, clusters)
  #write(t(cbind(res_cross_known, "known.10fold", n,p,B1,K, overdisp, 0.9)), file=filename, ncol=12, append=TRUE)
  #})
  
  #try({
  #overdisps.sct <- get_sctransform_bs(X)
  #res_cross_sct <- numClust_crossVal(X,10,clustTry,  overdisps.sct, clusters)
  #write(t(cbind(res_cross_sct, "sct.10fold", n,p,B1,K, overdisp, 0.9)), file=filename, ncol=12, append=TRUE)
#})
  
}
