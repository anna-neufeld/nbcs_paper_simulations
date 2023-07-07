library(tidyverse)
library(patchwork)
setwd("~/nbcs_paper/res")

file_names <- dir("~/nbcs_paper/res", pattern="jan9_propCor*") #where you have your files
res <- do.call(rbind,lapply(file_names,read.csv,sep="",header=FALSE))

names(res) <- c("mse", "rand", "k", "bestRand", "kGuess",
               "method", "n", "p", "B1", "K", "overdisp", "eps")

res$method2 <- paste(res$method, res$eps, sep="")
table(res$method2)


res$method2 <- ordered(res$method2, 
                       levels=c("known.onefold0.5",
                                "sct.onefold0.5",
                                "known.onefold0.9",
                                "sct.onefold0.9",
                                "known.10fold0.9",
                                "sct.10fold0.9"))

propCorrect_medium  <- res  %>% filter(k==K) %>% group_by(method2, K, eps, B1, overdisp) %>%
  summarize(correct = mean(K==kGuess), 
            rand = mean(rand), 
            correctRand = mean(K==bestRand)) %>%
  mutate(K2 = paste("True K* =", K), 
         overdisp2 = ifelse(overdisp==5, "Severe Overdispersion", "Mild Overdispersion"),
         method2 = ifelse(method2=="known.onefold0.9", "NBCS-known, \u03B5 = 0.9", "NBCV-known, 10 folds"))


p1 <- ggplot(data=propCorrect_medium %>% filter(overdisp !=10, eps > 0.5, K %in% c(1,3,5), 
                                                overdisp %in% c(1,5)), 
             aes(x=B1, y=correct, col=method2))+geom_line()+
  facet_grid( vars(overdisp2),vars(K2))+xlab(expression(beta))+ylab("Proportion of times we chose correct K")+
  theme_bw()+
  labs(col="Method")+ggtitle("Proportion of times we select correct K")+
  scale_color_manual(values=c("cyan", "navy"))


p1 + plot_layout(guides="collect") + theme(plot.title=element_text(size=12), 
                                         plot.subtitle = element_text(size=12),
                                         axis.title = element_text(size=10),
                                         legend.title = element_text(size=10),
                                         legend.text = element_text(size=10)) 
ggsave("~/Dropbox/nbCS_paper/v16-biometrics/figures/multifold_comp_low.png", width=8, height=5)




