library(tidyverse)
library(patchwork)
setwd("~/Dropbox/nbCS_paper/sims/New_figures_for_paper/res")

file_names <- dir("~/Dropbox/nbCS_paper/sims/New_figures_for_paper/res", pattern="jan9_roleEps*") #where you have your files
res <- do.call(rbind,lapply(file_names,read.csv,sep="",header=FALSE))

names(res) <- c("mse", "rand", "k", "bestRand", "kGuess", 
                "method", "n", "p", "B1", "K", "overdisp", "eps")

res$K2 <- paste("True K=", res$K)
res$K2 <- ordered(res$K2, levels=c("True K= 1", "True K= 3","True K= 5"))
res$overdisp2 <- "Severe Overdispersion"
res$overdisp2[res$overdisp==1] <- "Mild Overdispersion"

consRes <- res %>% filter(k==K, K ==5, n==500, p==40, overdisp==1) %>% group_by(K, K2, eps, B1, overdisp, overdisp2) %>%
  summarize(correctK = mean(K==kGuess), ### once per dataset. 
            meanRand = mean(rand), 
            correctRand = mean(K==bestRand)) 

consRes2 <- res %>% filter(k==K, rand > 0.8, K==5, n==500, p==40, overdisp==1) %>% group_by(K, K2, eps, B1, overdisp, overdisp2) %>%
  summarize(correctK = mean(K==kGuess), ### once per dataset. 
            meanRand = mean(rand), 
            correctRand = mean(K==bestRand), size=n()) 

p1 <- ggplot(data=consRes %>% filter(overdisp==1), aes(x=eps, y=meanRand, col=as.factor(B1), group=as.factor(B1)))+
  geom_point()+geom_line()+coord_fixed()+ylab("Adjusted Rand index")+
  ggtitle("", "Ability to estimate true clusters")+
  facet_grid(vars(K2), vars(overdisp2))

p2 <- ggplot(data=consRes2 %>% filter(overdisp==1, size > 2), aes(x=eps, y=correctK, col=as.factor(B1), group=as.factor(B1)))+
  geom_point()+geom_line()+coord_fixed()+ylab("Proportion of times we select correct K")+
  ggtitle("Ability to select K, among datasets where", "we estimate the clusters accurately")+
  facet_grid(vars(K2), vars(overdisp2))+ylim(0,1)

p3 <- ggplot(data=consRes %>% filter(overdisp==1), aes(x=eps, y=correctK, col=as.factor(B1), group=as.factor(B1)))+
  geom_point()+geom_line()+coord_fixed()+ylab("Proportion of times we select correct K")+ggtitle("", "Ability to select K, among all datasets")+
  facet_grid(vars(K2), vars(overdisp2))

p1+p2+p3+plot_layout(guides="collect", nrow=1) & theme_bw() &
  labs(col = expression(paste(beta, "*"))) &
  xlab(expression(epsilon)) &
  theme(
    axis.title.x = element_text(size=16),
    plot.title= element_text(size=11),
    plot.subtitle= element_text(size=12),
  )
ggsave("~/Dropbox/nbCS_paper/v14/figures/rollEps_cluster.png", width=12.5, height=4)
