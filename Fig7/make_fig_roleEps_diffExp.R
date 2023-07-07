library(tidyverse)
library(patchwork)
library(data.table)
setwd("~/Dropbox/nbCS_paper/sims/New_figures_for_paper/res")

file_names_genewise <- dir("~/Dropbox/nbCS_paper/sims/New_figures_for_paper/res", pattern="jan18_diffExpEps*") 
res_genewise <- do.call(rbind,lapply(file_names_genewise , read.csv ,sep="", header=FALSE))
names(res_genewise) <- c("j", "pval", "coeff", "b.est.glm", 
                         "rand", "method", "n", "p" , "B1", "K","overdisp", "eps")
res_genewise$eps2 <- paste("U+03B5", res_genewise$eps)

res_genewise$overdisp2 <- "Mild overdispersion"
res_genewise$overdisp2[res_genewise$overdisp==5] <- "Severe overdispersion"

table(res_genewise$B1)
table(res_genewise$method)

res_genewise$method2 <- "NBCS (known)"
res_genewise$method2[res_genewise$method == "sct"] <- "NBCS (estimated)"

res_genewise$method2 <- ordered(res_genewise$method2, 
                                levels=c("NBCS (known)","NBCS (estimated)"))

res_genewise <- res_genewise %>% 
  mutate(null = (j > 1000/20 | B1==0),
  null2 = paste("Null= ", null),
  overdisp2 = ifelse(overdisp==1, "Mild overdispersion", "Severe overdispersion"),
  globalNull = ifelse(B1==0, "Global Null", "Not Global Null"),
  eps2 = paste("\u03B5", "=", eps))



naiveCol <- "#E69F00" 
sampSplitCol <- "#56B4E9"
InfCol <- "#009E73"
knownCol <- "#0072B2"
sctCol <- "#CC79A7"


consRes_Rand <- res_genewise %>% filter(j==1) %>%
  group_by(j, B1, method2, overdisp2, eps) %>% summarize(meanRand = abs(mean(rand)))


#### Anna- remind yourself why you did this!!! Was it better than the altnerative? 
consRes_Power <- res_genewise %>% filter(j < p/20) %>%
  mutate(roundedSig = round(abs(coeff), 1)) %>%
  group_by(roundedSig, method2, overdisp2, eps) %>% summarize(power = mean(pval < 0.05))

Res_Power <- res_genewise %>% filter(j < p/20) 
random_indices <- sample(1:NROW(Res_Power), size=300000)

p_smooth <- ggplot(data=Res_Power[random_indices,], 
       aes(x=abs(coeff), y=as.numeric(pval < 0.05), col=as.factor(eps)))+
  geom_smooth(SE=F, method="glm", method.args=list(family="binomial"))+
  facet_grid(vars(overdisp2), row=vars(method2))+xlim(0,3)+
  ggtitle("Power")+theme_bw()+
  ylab("Power")+xlab(expression(widehat(beta[j]^`*`)))+
  guides(col="none")




p22 <- ggplot(data=consRes_Rand %>% filter(j==1), 
              aes(x=B1, y=meanRand, col=as.factor(eps)))+
  #geom_point()+
  geom_line(lwd=1.2)+
  #geom_smooth(method="glm", method.args=list(family="binomial"), se=F)+
  facet_grid(col=vars(overdisp2), row=vars(method2))+
  ggtitle("Gene-by-gene")+theme_bw()+
  ylab("Adjusted Rand Index")+
  xlab(expression(beta^`*`))+
  ggtitle("Ability to detect true clusters")

p33 <- ggplot(data=consRes_Power, 
              aes(x=roundedSig, y=power, col=as.factor(eps)))+
  geom_point()+geom_line()+
  facet_grid(vars(overdisp2), row=vars(method2))+xlim(0,3)+
  ggtitle("Power")+theme_bw()+
  ylab("Power")+xlab(expression(hat(beta)))

p22+p_smooth+plot_layout(guides="collect") & guides(lty="none") &
  labs(col=expression(epsilon))
ggsave("~/Dropbox/nbCS_paper/v14/figures/detect_power.png", width=10, height=7)
ggsave("~/Dropbox/Dissertation/nbcs_files/figures/detect_power.png", width=10, height=7)


