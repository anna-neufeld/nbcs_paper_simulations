### Implements NBCS. Splits x, where the overdispersion parameter is b, and epsilon*b = b1. 
betaBinSample <- function(x, b,b1) {
  p <- rbeta(length(x),b1,b-b1)
  return(rbinom(length(x),x,p))
}

### Analytical expression for variance of Xtrain, when NB(lambda, b) is split assuming that the 
### overdispersion parameter is actually bprime. Using tuning parameter epsilon. 
exactvarXtrain <- function(ep, lambda, b, bprime) {
  varXtrain <- lambda*ep*(1-ep)*(bprime + 1 + lambda/b + lambda)/(bprime+1)+ep^2*(lambda+lambda^2/b)
  return(varXtrain)
}

### Same as above, but for Xtest
exactvarXtest <- function(ep, lambda, b, bprime) {
  varXtest <- lambda*ep*(1-ep)*(bprime + 1 + lambda/b + lambda)/(bprime+1)+(1-ep)^2*(lambda+lambda^2/b)
  return(varXtest)
}

### The exact correlation expression, derived in the paper. 
exactCorExpression <- function(ep, lambda, b, bprime) {
  varX <- lambda+lambda^2/b
  varXtrain <- lambda*ep*(1-ep)*(bprime + 1 + lambda/b + lambda)/(bprime+1)+ep^2*(lambda+lambda^2/b)
  varXtest <- lambda*ep*(1-ep)*(bprime + 1 + lambda/b + lambda)/(bprime+1)+(1-ep)^2*(lambda+lambda^2/b)
  covar <- (varX - varXtrain - varXtest)/2
  cor <- covar / sqrt(varXtrain*varXtest)
  return(cor)
}


#### Empirically verify. 
n=100000
lambda <- 25
b=8
##### Generating negative binomial data using the poisson-gamma formulation. 
overdisps <- rgamma(n, b,b)
X <- rpois(n, lambda*overdisps)
ep = 0.3

##### Different values of b-prime to test out. 
bps <- 10^(seq(-log10(100000), log10(1000000), length.out=50)[-1])
cors <- rep(NA, length(bps))
varXtrain <- rep(NA, length(bps))
varXtest <- rep(NA, length(bps))
c <- 1

for (bp in bps) {
  Xtrain <- sapply(X,function(u) betaBinSample(u, bp, ep*bp))
  Xtest <- X-Xtrain
  cors[c] <-cor(Xtrain, Xtest)
  varXtrain[c] <- var(Xtrain)
  varXtest[c] <- var(Xtest)
  c <- c+1
}

trueCors <- sapply(bps, function(u) exactCorExpression(ep, lambda, b, u))
trueVarXtrain <- sapply(bps, function(u) exactvarXtrain(ep, lambda, b, u))
trueVarXtest <- sapply(bps, function(u) exactvarXtest(ep, lambda, b, u))


p1 <- ggplot(data=NULL) +
  geom_line(aes(x=bps, y=trueCors, col="Closed Form"), lwd=1.5)+
  geom_point( aes(x=bps, y=cors, col="Empirical"),lwd=2)+
  theme_bw()+
  geom_vline(xintercept=b, col="red")+
  xlab("b' (log scale)") + 
  ylab("Correlation(Xtrain, Xtest)")+
  scale_x_log10()+
  #geom_hline(yintercept=0, col="green")+
  geom_hline(yintercept = sqrt(ep*(1-ep)/(ep*(1-ep)+b/lambda+b^2/lambda^2)), col="green")+
  #geom_hline(col="green", yintercept = (-lambda*sqrt(ep*(1-ep)))/sqrt((1+2*lambda/b+lambda^2/b^2+ep*(1-ep)*lambda^2+(lambda+lambda^2/b))), col="green")+
  labs(col="")

p1
