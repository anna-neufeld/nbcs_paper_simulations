library(tidyverse)
library(patchwork)
setwd("~/nbcs_paper/res")


### Read in the results which were printed to files. 
file_names <- dir("~/nbcs_paper/res", pattern="jan9_mse*") #where you have your files
res <- do.call(rbind,lapply(file_names,read.csv,sep="",header=FALSE))

names(res) <- c("mse", "rand", "k", "bestRand", "kGuess",
               "method", "n", "p", "B1", "K", "overdisp", "eps")

res$K2 <- paste("True K=", res$K)
res$K2 <- ordered(res$K2, levels=c("True K= 1", "True K= 3","True K= 5"))
res$method2 <- ordered(res$method, 
                       levels=c("naive","inf.onefold", "known.onefold", "sct.onefold"))
res$overdisp2 <- "Severe Overdispersion"
res$overdisp2[res$overdisp==1] <- "Mild Overdispersion"


consRes <- res %>% group_by(n,p,k, K, K2, B1 ,overdisp, overdisp2, method2) %>% 
  summarize(meanMSE=mean(mse))
consRes2 <- consRes %>% group_by(n,p,method2, K2, B1, overdisp, overdisp2) %>% 
  mutate(normalizeMeanMSE = (meanMSE - min(meanMSE))/(max(meanMSE)-min(meanMSE)))

naiveCol <- "#E69F00" 
InfCol <- "#009E73"
knownCol <- "#0072B2"
sctCol <- "#CC79A7"

ggplot(data=consRes2, aes(x=k, y=normalizeMeanMSE+0.001, col=method2,
                         lty=method2))+
  geom_line(lwd=1)+
  facet_grid(cols=vars(K2), vars(as.factor(overdisp2)))+
  theme_bw()+
  scale_color_manual(
    values=c(naiveCol, InfCol, knownCol, sctCol),
    labels=c("Naive", "PCS",
             expression("NBCS (known),"~epsilon~"= 0.5"),
             expression("NBCS (estimated),"~epsilon~"= 0.5")))+
  guides(lty="none")+
  geom_vline(aes(xintercept=K))+
  xlab("Number of estimated clusters")+
  ylab("Scaled within cluster MSE (log scale)")+
  scale_x_continuous(breaks=seq(1,10,1))+
  labs(col="Method")+
  scale_y_log10() +
  theme(legend.text.align = 0)
ggsave("~/Dropbox/nbCS_paper/v14/figures/main_MSE_highDim.png", width=10, height=6)

