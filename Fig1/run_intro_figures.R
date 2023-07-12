library(tidyverse)
library(patchwork)


### This is the main function that will be used for NBCS.
### Also implemented in the forthcoming R package for NBCS. 
betaBinSample <- function(x, b,eps) {
  p <- rbeta(length(x),eps*b,(1-eps)*b)
  return(rbinom(length(x),x,p))
}


#### Make panel a of figure 1
set.seed(2)
X <- matrix(rnbinom(200, mu=5, size=5), ncol=2)
hX <- log(X+1)
predMean <- colMeans(hX)
totalMSE <- sum(apply(hX, 1, function(u) (u-predMean)^2))/200
p1 <- ggplot(data=NULL, aes(x=hX[,1], y=hX[,2]))+geom_point()+xlab("Gene 1 (log scale)")+ ylab("Gene 2 (log scale)")+
  theme_bw()+coord_fixed()+theme(panel.grid.major=element_blank(), 
                                 panel.grid.minor=element_blank())+
  ggtitle("(a) Example data realization")+
  xlim(0,3)+ylim(0,3)

### Make panel b of Figure 1
set.seed(2)
cluster.full1 <- kmeans(hX, centers=2, nstart=20)$cluster
k2_SSE <- 0
for (k in 1:2) {
  clusterDat <- hX[cluster.full1==k,]
  clusterMean <- colMeans(clusterDat)
  k2_SSE <- k2_SSE + sum(apply(clusterDat,1,function(u) (u-clusterMean)^2))
}
k2_MSE <- k2_SSE/200
pval <- summary(MASS::glm.nb(X[,1]~as.factor(cluster.full)))$coefficients[2,4]
p2 <- ggplot(data=NULL, aes(x=hX[,1], y=hX[,2], col=as.factor(cluster.full1)))+
  geom_point()+xlab("Gene 1 (log scale)")+ ylab("Gene 2 (log scale)")+
  theme_bw()+coord_fixed()+theme(panel.grid.major=element_blank(), 
                                 panel.grid.minor=element_blank())+
  guides(col="none")+
  xlim(0,3)+ylim(0,3)+
  ggtitle("(b) Example data realization,", "colored by estimated cluster")

#### Carries out the MSE simulation for Panel C of Figure 1
ktry <- 10
nTrials <- 1000
allMSEs <- matrix(NA, nrow=nTrials, ncol=ktry)
sampSplitMSEs <- matrix(NA, nrow=nTrials, ncol=ktry)
countSplitMSEs <- matrix(NA, nrow=nTrials, ncol=ktry)

for (t in 1:nTrials) {
  set.seed(t) 
  X <- matrix(rnbinom(200, mu=5, size=5), ncol=2)
  Xtrain <- X[1:50,]
  Xtest <- X[51:100,]
  XtrainCS <- apply(X, 1:2, function(u) betaBinSample(u, 5, 0.5))
  XtestCS <- X - XtrainCS
  
  hX <- log(X+1)
  hXtrain <- hX[1:50,]
  hXtest <- hX[51:100,]
  hXtrainCS <- log(XtrainCS+1)
  hXtestCS <- log(XtestCS+1)
  
  for (k in 1:ktry) {
    cluster.full.try <- kmeans(hX, centers=k)$cluster
    cluster.train.try <- kmeans(hXtrain, centers=k)$cluster
    cluster.test.try <- class::knn(hXtrain, hXtest, cluster.train.try,3)
    cluster.cs.try <- kmeans(hXtrainCS, centers=k)$cluster
    k_SSE <- 0
    k_SSE_train <- 0
    k_SSE_cs <- 0
    for (kp in 1:k) {
      clusterDat <- X[cluster.full.try==kp,,drop='F']
      clusterMean <- colMeans(clusterDat)
      
      clusterDatTrain <- Xtrain[cluster.train.try == kp,,drop='F']
      clusterDatTest <- Xtest[cluster.test.try== kp,,drop='F']
      clusterDatCS.Train <- XtrainCS[cluster.cs.try== kp,,drop='F']
      clusterDatCS.Test <- XtestCS[cluster.cs.try== kp,,drop='F']
      
      clusterMeanTrain <- colMeans(clusterDatTrain)
      clusterMeanTrainCS <- colMeans(clusterDatCS.Train)
      
      k_SSE <- k_SSE + sum(apply(clusterDat,1,function(u) (log(u+1)-log(clusterMean+1))^2))
      k_SSE_train <- k_SSE_train + sum(apply(clusterDatTest,1,function(u) (log(u+1)-log(clusterMeanTrain+1))^2))
      k_SSE_cs <- k_SSE_cs + sum(apply(clusterDatCS.Test,1,function(u) (log(u+0.5)-log( clusterMeanTrainCS+0.5))^2))
    }
    allMSEs[t, k] <- k_SSE/200
    sampSplitMSEs[t,k] <- k_SSE_train/100
    countSplitMSEs[t,k] <- k_SSE_cs/200
  }
}

naiveCol <- "purple"
sampSplitCol <- "grey"
csCol <- "blue"
pMSEthree <- ggplot(data=NULL, aes(x=1:ktry, y=colMeans(sampSplitMSEs, na.rm=TRUE), col="bSample Splitting"))+
  geom_line()+theme_bw()+geom_point(aes(col="bSample Splitting"))+
  geom_point(aes(x=1:ktry, y=colMeans(allMSEs), col="aUse data twice"))+
  geom_line(aes(x=1:ktry, y=colMeans(allMSEs), col="aUse data twice"))+
  geom_point(aes(x=1:ktry, y=colMeans(countSplitMSEs), col="cOur proposed method"))+
  geom_line(aes(x=1:ktry, y=colMeans(countSplitMSEs), col="cOur proposed method"))+
  theme(panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank())+
  xlab("Number of clusters")+ylab("Average MSE")+
  ggtitle("(c) Average within cluster MSE,", "over 1000 datasets")+
  scale_color_manual(values=c(naiveCol, sampSplitCol,csCol))+
  guides(col="none")

### Carries out the differential expression simulation for Panel D of Figure 1
nTrials <- 1000
pvalsNaive <- rep(NA, nrow=nTrials)
pvals_sampsplit <- rep(NA, nrow=nTrials)
pvals_cs <- rep(NA, nrow=nTrials)

for (t in 1:nTrials) {
  set.seed(t) 
  X <- matrix(rnbinom(200, mu=5, size=5), ncol=2)
  hX <- log(X+1)
  hXtrain <- hX[1:50,]
  hXtest <- hX[51:100,]
  Xtest <- X[51:100,]
  XtrainCS <- apply(X, 1:2, function(u) betaBinSample(u, 5, 0.5))
  XtestCS <- X - XtrainCS
  
  
  cluster.full.try <- kmeans(hX, centers=2)$cluster
  cluster.train.try <- kmeans(hXtrain, centers=2)$cluster
  cluster.test.try <- class::knn(hXtrain, hXtest, cluster.train.try,3)
  cluster.cs.try <- kmeans(log(XtrainCS+1), centers=2)$cluster
  
  
  pvalsNaive[t] <- summary(MASS::glm.nb(X[,1]~as.factor(cluster.full.try)))$coefficients[2,4]
  pvals_sampsplit[t] <- summary(MASS::glm.nb(Xtest[,1]~as.factor(cluster.test.try)))$coefficients[2,4]
  pvals_cs[t] <- summary(MASS::glm.nb(XtestCS[,1]~as.factor(cluster.cs.try)))$coefficients[2,4]
}

pQQthree <- ggplot(data=NULL)+
  geom_qq(aes(sample=as.numeric(pvals_sampsplit), col="bSample Splitting"), distribution=stats::qunif)+
  geom_qq(aes(sample=as.numeric(pvalsNaive), col="aUse data twice"), distribution=stats::qunif)+
  geom_qq(aes(sample=as.numeric(pvals_cs), col="cOur proposed method"), distribution=stats::qunif)+
  theme_bw()+
  scale_color_manual(values=c(naiveCol, sampSplitCol, csCol),
                     labels=c("Use full data twice", "Sample Splitting", "Our proposed method"))+
  geom_abline(intercept=0,slope=1)+
  ggtitle("(d) P-values for gene 1,", "over 1000 datasets")+
  xlab("Unif(0,1) Quantiles")+ylab("Sample Quantiles")+
  coord_fixed()+labs(col="Method")


#### Carries out the IDCV procedure for Panel E of Figure 1. 
set.seed(2)
n <- 100
p <- 2
Xnull <- matrix(rnbinom(n*p, mu=5, size=5), ncol=p)
clusters <- kmeans(log(Xnull+1), centers=5)$cluster
Y <- log(Xnull+1)
numClust <- 5
folds <- 5
clusterLabs <- as.factor(kmeans(Y, centers=numClust, nstart=50)$cluster)
confusionMatrix <- matrix(0, nrow=numClust, ncol=numClust)
crossValIndices <- sample(1:n, size=n, replace=FALSE)
counter <- 1
for (f in 1:folds) {
  testIndices <-  crossValIndices[counter:(counter+n/folds-1)]
  counter <- counter+n/folds
  Ytrain <- Y[-testIndices,]
  Ytest <- Y[testIndices,]
  clustersTrain <- clusterLabs[-testIndices]
  clustersTest <- clusterLabs[testIndices]
  classifier <- e1071::svm(x=Ytrain, y=as.factor(clustersTrain), kernel="linear")
  preds <- predict(classifier, newdata=Ytest)
  confusionMatrix <- confusionMatrix + table(preds, clustersTest)
}
p.IDCV.cao <- ggplot(data=data.frame(confusionMatrix/colSums(confusionMatrix)), aes(x= clustersTest, y= preds, fill=Freq))+geom_tile()+
  xlab("Cluster predicted by SVM")+
  ylab("Cluster estimated with kmeans")+
  ggtitle("(e) Intradataset cross validation", 
          "that uses data twice")+
  scale_fill_viridis_c(begin=0, end=1, limits=c(0,1))+
  labs(fill="Proportion") + coord_fixed()

sum(diag(confusionMatrix))/sum(confusionMatrix)

#### Carries out the simulation for Panel F of Figure 1
Xtrain <- t(apply(Xnull, 1, function(u) betaBinSample(u, 1, 0.5)))
Xtest <- Xnull-Xtrain
#### Cluster both datasets on log scale. 
Ytrain <- log(Xtrain+1)
Ytest <- log(Xtest+1)
clustersTrain <- kmeans(Ytrain, centers=numClust, nstart=50)$cluster
clustersTest <- kmeans(Ytest, centers=numClust, nstart=50)$cluster
myTab <- table(clustersTrain, clustersTest)
p.IDCV.cs <- ggplot(data=data.frame(myTab/colSums(myTab)), 
                    aes(x= clustersTrain, y= clustersTest, fill=Freq))+geom_tile()+
  xlab(expression("Cluster estimated from"~X^{(1)}))+
  ylab(expression("Cluster estimated from"~X^{(2)}))+
  ggtitle("(f) Intradataset cross validation", 
          "using our proposed method")+
  coord_fixed()+
  labs(fill="Proportion")+
  scale_fill_viridis_c(begin=0, end=1, limits=c(0,1))




#### MAKE THE FIGURE
library(viridis)
p1+p2+pMSEthree+pQQthree+p.IDCV.cao+p.IDCV.cs+
  plot_layout(guides="collect", nrow=2, widths=c(1,1,1)) & 
  theme(plot.title=element_text(size=12), 
        plot.subtitle = element_text(size=12),
        axis.title = element_text(size=10),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10)) &
  labs(fill = "Proportion of cells \nin column \nbelonging to row")
ggsave("~/Dropbox/nbCS_paper/v16-biometrics/figures/intro_MEGA.png",
       width=11, height=6)











