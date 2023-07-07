library(tidyverse)
library(patchwork)
library(viridis)
setwd("~/nbcs_paper/Fig8/real_data_results/")

#### PANEL A
main_idcv_table <- readRDS("3_fca_kidney_full_confusionmatrix.rds")
main_idcv_prop <- as.table(apply(main_idcv_table , 2, function(u) u/sum(u)))
main_idcv_df <- data.frame(main_idcv_prop)
names(main_idcv_df ) <- c("SVM_cluster", "clustering_cluster", "Freq")

fake_data <- c(NA, NA)
for (row in 1:NROW(main_idcv_table)){
  for (col in 1:NCOL(main_idcv_table)) {
    fake_data <- rbind(fake_data, matrix(
     rep(
       c(rownames(main_idcv_table)[row], 
         colnames(main_idcv_table)[col]
         ), main_idcv_table[row, col]), 
     ncol=2, byrow=TRUE))
  }
}
fake_data <- data.frame(fake_data[-1,])
names(fake_data) <- c("SVM", "cluster")

pA <- ggplot(data=main_idcv_df, 
       aes(x=clustering_cluster, y=SVM_cluster, fill=Freq))+
  geom_tile()+
  scale_fill_viridis_c(
    begin=0, end=1, limits=c(0,1))+
  xlab("Original cluster assignment")+
  ylab("Cluster predicted by SVM")+
  ggtitle(paste("(a) ARI = ", signif(mclust::adjustedRandIndex(fake_data$SVM, fake_data$cluster),4)))

#### PANEL E
sub_idcv_table <- readRDS("3_fca_metanephric_full_confusionmatrix.rds")
sub_idcv_prop <- as.table(apply(sub_idcv_table , 2, function(u) u/sum(u)))
sub_idcv_df <- data.frame(sub_idcv_prop)
names(sub_idcv_df ) <- c("SVM_cluster", "clustering_cluster", "Freq")

fake_data_sub <- c(NA, NA)
for (row in 1:NROW(sub_idcv_table)){
  for (col in 1:NCOL(sub_idcv_table)) {
    fake_data_sub <- rbind(fake_data_sub, matrix(
      rep(
        c(rownames(sub_idcv_table)[row], 
          colnames(sub_idcv_table)[col]
        ), sub_idcv_table[row, col]), 
      ncol=2, byrow=TRUE))
  }
}
fake_data_sub <- data.frame(fake_data_sub[-1,])
names(fake_data_sub) <- c("SVM", "cluster")
pE <- ggplot(data=sub_idcv_df, 
             aes(x=clustering_cluster, y=SVM_cluster, fill=Freq))+
  geom_tile()+
  scale_fill_viridis_c(
    begin=0, end=1, limits=c(0,1))+
  xlab("Original cluster assignment")+
  ylab("Cluster predicted by SVM")+
  ggtitle(paste("(e) ARI = ", signif(mclust::adjustedRandIndex(fake_data_sub$SVM, fake_data_sub$cluster),4)))



  
### PANEL B
main_pois_test <- read.csv("1_fca_kidney_poistest_clusters_1e-6_75.tsv", sep="\t")
names(main_pois_test) <- c("name", "test_cluster")
main_pois_train <- read.csv("1_fca_kidney_poistrain_clusters_1e-6_75.tsv", sep="\t")
names(main_pois_train) <- c("name", "train_cluster")
main_pois_confusion_tab <- table(main_pois_test$test_cluster,main_pois_train$train_cluster)
main_pois_confusion_prop <- as.table(apply(main_pois_confusion_tab, 2, function(u) u/sum(u)))
main_pois_confusion_df <- data.frame(main_pois_confusion_prop)
names(main_pois_confusion_df) <- c("test_cluster", "train_cluster", "Freq")
main_pois_confusion_df$test_cluster <- 
  factor(main_pois_confusion_df$test_cluster, ordered=TRUE,
         levels=c("1", "2", "3", "4", "5", "6",
                  "8", "7"))
main_pois_confusion_df$train_cluster <- 
  factor(main_pois_confusion_df$train_cluster, ordered=TRUE,    
     levels=c("1", "2", "3", "4", "5", "8", "6", "7"))
pB <- ggplot(data=main_pois_confusion_df, aes(x=train_cluster, y=test_cluster, fill=Freq))+
  geom_tile()+
  scale_fill_viridis_c(
    begin=0, end=1, limits=c(0,1))+
  xlab(expression("Cluster estimated from"~X^{(1)}))+
  ylab(expression("Cluster estimated from"~X^{(2)}))+
  ggtitle(paste("(b) ARI = ", 
                signif(mclust::adjustedRandIndex(
                  main_pois_test$test_cluster,main_pois_train$train_cluster),4)))

### PANEL C
main_nb_test <- read.csv("1_fca_kidney_nbtest_clusters_1e-6_75.tsv", sep="\t")
names(main_nb_test) <- c("name", "test_cluster")
main_nb_train <- read.csv("1_fca_kidney_nbtrain_clusters_1e-6_75.tsv", sep="\t")
names(main_nb_train) <- c("name", "train_cluster")
main_nb_confusion_tab <- table(main_nb_test$test_cluster,main_nb_train$train_cluster)
main_nb_confusion_prop <- as.table(apply(main_nb_confusion_tab, 2, function(u) u/sum(u)))
main_nb_confusion_df <- data.frame(main_nb_confusion_prop)
names(main_nb_confusion_df) <- c("test_cluster", "train_cluster", "Freq")
main_nb_confusion_df$test_cluster <- 
  factor(main_nb_confusion_df$test_cluster, ordered=TRUE,
  levels = c("1", "2", "3", "6","4", "5", "7"))
main_nb_confusion_df$train_cluster <- 
  factor(main_nb_confusion_df$train_cluster, ordered=TRUE,
         levels = c("1", "2", "7", "3", "4", "5", "9", "6", "8"))
pC <- ggplot(data=main_nb_confusion_df, aes(x=train_cluster, y=test_cluster, fill=Freq))+
  geom_tile()+
  scale_fill_viridis_c(
    begin=0, end=1, limits=c(0,1))+
  xlab(expression("Cluster estimated from"~X^{(1)}))+
  ylab(expression("Cluster estimated from"~X^{(2)}))+
  ggtitle(paste("(c) ARI = ", 
                signif(mclust::adjustedRandIndex(
                  main_nb_test$test_cluster,main_nb_train$train_cluster),4)))


### PANEL E
sub_pois_test <- read.csv("2_fca_metanephric_poistest_clusters_2e-5_75.tsv", sep="\t")
names(sub_pois_test) <- c("name", "test_cluster")
sub_pois_train <- read.csv("2_fca_metanephric_poistrain_clusters_2e-5_75.tsv", sep="\t")
names(sub_pois_train) <- c("name", "train_cluster")
sub_pois_confusion_tab <- table(sub_pois_test$test_cluster,sub_pois_train$train_cluster)
sub_pois_confusion_prop <- as.table(apply(sub_pois_confusion_tab, 2, function(u) u/sum(u)))
sub_pois_confusion_df <- data.frame(sub_pois_confusion_prop)
names(sub_pois_confusion_df) <- c("test_cluster", "train_cluster", "Freq")
sub_pois_confusion_df$test_cluster <- 
  factor(sub_pois_confusion_df$test_cluster, ordered=TRUE,
      levels=c("4", "1", "2", "3", "5", "6", 
               "7", "8", "9", "10", "11", "12"))
sub_pois_confusion_df$train_cluster <- 
  factor(sub_pois_confusion_df$train_cluster, ordered=TRUE,
         levels=c("1", "2", "7", "3", "4", "6", 
                  "5",  "8","9", "10"))

pF <- ggplot(data=sub_pois_confusion_df, aes(x=train_cluster, y=test_cluster, fill=Freq))+
  geom_tile()+
  scale_fill_viridis_c(
    begin=0, end=1, limits=c(0,1))+
  xlab(expression("Cluster estimated from"~X^{(1)}))+
  ylab(expression("Cluster estimated from"~X^{(2)}))+
  ggtitle(paste("(f) ARI = ", 
                signif(mclust::adjustedRandIndex(
                  sub_pois_test$test_cluster,sub_pois_train$train_cluster),4)))

### PANEL G
sub_nb_test <- read.csv("2_fca_metanephric_nbtest_clusters_2e-5_75.tsv", sep="\t")
names(sub_nb_test) <- c("name", "test_cluster")
sub_nb_train <- read.csv("2_fca_metanephric_nbtrain_clusters_2e-5_75.tsv", sep="\t")
names(sub_nb_train) <- c("name", "train_cluster")
sub_nb_confusion_tab <- table(sub_nb_test$test_cluster,sub_nb_train$train_cluster)
sub_nb_confusion_prop <- as.table(apply(sub_nb_confusion_tab, 2, function(u) u/sum(u)))
sub_nb_confusion_df <- data.frame(sub_nb_confusion_prop)
names(sub_nb_confusion_df) <- c("test_cluster", "train_cluster", "Freq")
sub_nb_confusion_df$test_cluster <- 
  factor(sub_nb_confusion_df$test_cluster, ordered=TRUE,
         levels = c("4", "1", "6", "7", "3", "5", "2", 
                "8", "9", "10", "11", "12"))
sub_nb_confusion_df$train_cluster <- 
  factor(sub_nb_confusion_df$train_cluster, ordered=TRUE,
         levels=c("1", "2", "3", "4", 
                  "5", "6", "8","7","9", "10"))


pG <- ggplot(data=sub_nb_confusion_df, aes(x=train_cluster, y=test_cluster, fill=Freq))+
  geom_tile()+
  scale_fill_viridis_c(
  begin=0, end=1, limits=c(0,1))+
  xlab(expression("Cluster estimated from"~X^{(1)}))+
  ylab(expression("Cluster estimated from"~X^{(2)}))+
  ggtitle(paste("(g) ARI = ", 
                signif(mclust::adjustedRandIndex(
                  sub_nb_test$test_cluster,sub_nb_train$train_cluster),4)))


load("4_fca_metanephric_poistrain_clusters_2e-5_ari_all.tsv")
ari_pois_met <- ari_list
load("4_fca_kidney_poistrain_clusters_1e-6_ari_all.tsv")
ari_pois_kidney <- ari_list
load("4_fca_metanephric_nbtrain_clusters_2e-5_ari_all.tsv")
ari_nb_met <- ari_list
load("4_fca_kidney_nbtrain_clusters_1e-6_ari_all.tsv")
ari_nb_kidney  <- ari_list

ari_met <- data.frame(
  "ari" = c(ari_pois_met, ari_nb_met), 
  "method" = c(rep("Assume Poisson", 10), rep("sctransform",10)),
  "dataset" = "Metenephric Cells"
)
ari_kidney <- data.frame(
  "ari" = c(ari_pois_kidney, ari_nb_kidney), 
  "method" = c(rep("Assume Poisson", 10), rep("sctransform",10)),
  "dataset" = "Kidney Cells"
)

pH <- ggplot(data=ari_met, aes(x=method, y=ari))+
  geom_boxplot()+
  ggtitle("(h) ARI over 10 random splits")+
  ylab("Adjusted Rand Index")+
  xlab("Method")+
  theme_bw()+
  ylim(0,1)

pD <- ggplot(data=ari_kidney, aes(x=method, y=ari))+
  geom_boxplot()+
  ggtitle("(d) ARI over 10 random splits")+
  ylab("Adjusted Rand Index")+
  xlab("Method")+
  theme_bw()+
  ylim(0,1)


ari_both <- rbind(ari_kidney, ari_met)
pBoth <- ggplot(data=ari_both %>% filter(method=="sctransform"), 
                aes(y=ari, x=dataset))+geom_boxplot()+
  ylim(0,1)+theme_bw()+ylab("Adjusted Rand Index")+xlab("")
pBoth 
ggsave("boxplot.png", width=3, height=3)

pC + ggtitle("")+ guides(fill="none")
ggsave("~/Dropbox/Spring 2023 talks/main.png", width=3, height=3)

pG+ ggtitle("") + guides(fill="none")
ggsave("~/Dropbox/Spring 2023 talks/sub.png", width=3, height=3)




pA+pB+pC+pD+pE+pF+pG+pH+
  plot_layout(guides="collect", nrow=2) &
  labs(fill = "Proportion of cells \nin column \nbelonging to row") &
  theme(legend.title = element_text(size=12))

setwd("~/Dropbox/nbCS_paper/v14/figures/")
ggsave("realData.png", height=6, width=14)


#### For slides
pC
pG


