q#!/usr/local/bin/Rscript

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
parser$add_argument("--nreps", type = "double", default = 10,
                    help = "number of replicates for each set of params")
args <- parser$parse_args()


## -----------------------------------------
## set up a grid of parameters to cycle over
## -----------------------------------------
ns <- c(1000)
ps <- c(1000)
B1s <- seq(0,3,length.out=32)
overdisps <- c(1,5)
K <- 2

nreps_per_combo <- args$nreps

param_grid <- expand.grid(n=ns,
              p=ps,
              B1=B1s,
              overdisp=overdisps)


 
jobid <- as.numeric(Sys.getenv("SGE_TASK_ID"))


filename <- paste("res/",args$simname, jobid, ".txt", sep="")

for (i in 1:NROW(param_grid)) {


  current_dynamic_args <- param_grid[i, ]
  B0s <- rnorm(current_dynamic_args$p)

  
 n=current_dynamic_args$n
p=current_dynamic_args$p
B1 = current_dynamic_args$B1
K = 2
clustTry = 1:10
overdisp=current_dynamic_args$overdisp
  clusters <- sample(1:K, size=n, replace=TRUE)
  
  logLambda <- matrix(B0s, nrow=n, ncol=p, byrow = TRUE)
  
  
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
  

 
  #res_naive <- diffExp_naive(X,clusters)
  #write(t(cbind(1:p, res_naive, "naive", n,p,B1,K, B0s, overdisp, overdisps, NA, NA)), file=filename, ncol=15, append=TRUE)
  
  #res_ss <- diffExp_sampsplit(X,clusters, eps.train=0.5)
  #write(t(cbind(1:p, res_ss, "sampsplit", n,p,B1,K, B0s, overdisp, overdisps, NA, 0.5)), file=filename, ncol=15, append=TRUE)
  
  overdisps.sct <- try(get_sctransform_bs(X))
  
  for (eps in c(0.1, 0.25, 0.5, 0.75, 0.9)) {
    res_cs_known <- diffExp_datathin(X,clusters, overdisps, eps)
    write(t(cbind(1:(p/10), res_cs_known[1:(p/10),] , "known", n,p,B1,K, overdisp, eps)), file=filename, ncol=12, append=TRUE)

    if (class(overdisps.sct) != "try-error") {
  res_cs_sct <- diffExp_datathin(X,clusters, overdisps.sct, eps)
   write(t(cbind(1:(p/10), res_cs_sct[1:(p/10),] , "sct", n,p,B1,K, overdisp, eps)), file=filename, ncol=12, append=TRUE)
    }
  }
}
