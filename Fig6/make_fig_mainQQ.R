library(tidyverse)
library(patchwork)
setwd("~/Dropbox/nbCS_paper/sims/New_figures_for_paper/res")

file_names_genewise <- dir("~/Dropbox/nbCS_paper/sims/New_figures_for_paper/res", pattern="Jan9_QQ*") 
res_genewise <- do.call(rbind,lapply(file_names_genewise,read.csv,sep="",header=FALSE))
names(res_genewise) <- c("j", "pval", "coeff", "b.est.glm", 
                         "rand", "method", "n", "p" , "B1", "K", "B0", "overdisp" ,"overdisps",
                         "overdisps.hat", "eps")

res_genewise$method2 <- ordered(res_genewise$method, 
                       levels=c("naive","Inf", "known", "sct"))


res_genewise <- res_genewise %>% 
  mutate(null = (j > p/20 | B1==0),
  null2 = paste("Null= ", null),
  overdisp2 = ifelse(overdisp==1, "Mild overdispersion", "Severe overdispersion"),
  globalNull = ifelse(B1==0, "Global Null", "Not Global Null"))

nullRes_genewise <- (res_genewise %>% filter(null==TRUE))
nullRes_genewise %>% group_by(B1) %>% summari

#### LOW DIM FIGURE.

naiveCol <- "#E69F00" 
sampSplitCol <- "#56B4E9"
InfCol <- "#009E73"
knownCol <- "#0072B2"
sctCol <- "#CC79A7"



nullRes_genewise[nullRes_genewise$method=="known",]$pval <- 
  nullRes_genewise[nullRes_genewise$method=="known",]$pval+0.02

subset <- sample(1:NROW(nullRes_genewise), size= 50000)
nullRes_genewise_sub <- nullRes_genewise[subset,]

ggplot(data=nullRes_genewise_sub %>% filter(B1 < 3), aes(sample=as.numeric(pval), col=method2))+
  geom_qq(distribution="qunif")+
  facet_grid(col=vars(overdisp2))+
               geom_abline(intercept=0,slope=1)+theme_bw()+
  scale_color_manual(
    values=c(naiveCol, InfCol, knownCol, sctCol),
    labels=c("Naive", "PCS",
             expression("NBCS (known),"~epsilon~"= 0.5"),
             expression("NBCS (estimated),"~epsilon~"= 0.5")))+
  coord_fixed()+
  xlab("Unif(0,1) Quantiles")+
  ylab("Sample Quantiles")+
  labs(col="Method")+
  theme(legend.text.align = 0)

ggsave("~/Dropbox/nbCS_paper/v7/figures/low_dim_qq.png", width=8, height=5)
