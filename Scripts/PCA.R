library(ggplot2)
library(reshape2)
library(tidyverse)
library(cluster)
library(pvclust)


counts.merged.long <- readRDS('data/processed/counts.merged.long.rds')
bulk.meta.edit <- read.csv("data/raw/meta.bulk.edited.csv")

#Format
bulk <- counts.merged.long[,1:3] %>%  
  dcast(gene~sample, sum, value.var = 'bulk') 
row.names(bulk) <- bulk$gene
bulk <- bulk[,-1]
bulkized <- counts.merged.long[,-3:-4] %>%
  dcast(gene~sample, sum, value.var = 'bulkized.norm')
row.names(bulkized) <- bulkized$gene
bulkized <- bulkized[,-1]

bulk <- as.matrix(bulk)
bulkized <- as.matrix(bulkized)

#Log-transform 
bulkized.log <- log2(1+ bulkized)
bulk.log <- log2(1+ bulk)


#Histograms
png("results/figures/bulk.hist.expr.png")
hist(bulk.log, col = "gray", xlab="Log Expression", main="Distribution of expression in bulk")
dev.off()

png("results/figures/bulkized.hist.expr.png")
hist(bulkized.log, col = "gray", xlab="Log Expression", main="Distribution of expression in bulkized")
dev.off()

#PCs
pcs.bulk <- prcomp(t(bulk.log))
pcs.bulkized <- prcomp(t(bulkized.log))

#Plot PCs
png("results/figures/bulk.hist.PCs.png")
plot(pcs.bulk, main="PCs of Bulk")
dev.off()

png("results/figures/bulkized.hist.PCs.png")
plot(pcs.bulkized, main="PCs of Bulkized")
dev.off()

#Format metadata
meta <- bulk.meta.edit[,-1:-2]
rownames(meta) <- bulk.meta.edit$Sample.name

#Get PCs
x.bulk <- as.data.frame(pcs.bulk$x)
x.bulkized <- as.data.frame(pcs.bulkized$x)

#Combine 
prinComps.bulk <- cbind(meta, x.bulk)
prinComps.bulkized <- cbind(meta,x.bulkized)
prinComps.bulk <- prinComps.bulk[,1:13]
prinComps.bulkized <- prinComps.bulkized[,1:13]


#Format for plots
prinComps.bulk <- prinComps.bulk %>% mutate(PMI_group= case_when(PMI < 10 ~ "0-9",PMI < 15 ~ "10-14", PMI < 20 ~ "15-19", PMI < 25 ~ "20-24", PMI < 30 ~ "25-29", PMI < 30 ~ "25-29",
                                                                 PMI < 35 ~ "30-34",PMI < 40 ~ "35-49", PMI >= 40 ~ ">40"))
prinComps.bulkized <- prinComps.bulkized %>% mutate(PMI_group= case_when(PMI < 10 ~ "0-9",PMI < 15 ~ "10-14", PMI < 20 ~ "15-19", PMI < 25 ~ "20-24", PMI < 30 ~ "25-29", PMI < 30 ~ "25-29",
                                                                         PMI < 35 ~ "30-34",PMI < 40 ~ "35-49", PMI >= 40 ~ ">40"))

prinComps.bulk <- prinComps.bulk %>% mutate(RIN_group= case_when(RIN < 7 ~ "6.5-6.9",RIN < 7.5 ~ "7.0-7.4",RIN < 8 ~ "7.5-7.0",RIN < 8.5 ~ "8.0-8.4", RIN < 9 ~ "8.5-8.9", RIN == 9 ~ "9"))
prinComps.bulkized <- prinComps.bulkized %>% mutate(RIN_group= case_when(RIN < 7 ~ "6.5-6.9",RIN < 7.5 ~ "7.0-7.4",RIN < 8 ~ "7.5-7.0",RIN < 8.5 ~ "8.0-8.4", RIN < 9 ~ "8.5-8.9", RIN == 9 ~ "9"))
library(ggplot2)
library(reshape2)
library(tidyverse)
library(cluster)
library(pvclust)


counts.merged.long <- readRDS('data/processed/counts.merged.long.rds')
bulk.meta.edit <- read.csv("data/raw/meta.bulk.edited.csv")

#Format
bulk <- counts.merged.long[,1:3] %>%  
  dcast(gene~sample, sum, value.var = 'bulk') 
row.names(bulk) <- bulk$gene
bulk <- bulk[,-1]
bulkized <- counts.merged.long[,-3:-4] %>%
  dcast(gene~sample, sum, value.var = 'bulkized.norm')
row.names(bulkized) <- bulkized$gene
bulkized <- bulkized[,-1]

bulk <- as.matrix(bulk)
bulkized <- as.matrix(bulkized)

#Log-transform 
bulkized.log <- log2(1+ bulkized)
bulk.log <- log2(1+ bulk)


#Histograms
png("results/figures/bulk.hist.expr.png")
hist(bulk.log, col = "gray", xlab="Log Expression", main="Distribution of expression in bulk")
dev.off()

png("results/figures/bulkized.hist.expr.png")
hist(bulkized.log, col = "gray", xlab="Log Expression", main="Distribution of expression in bulkized")
dev.off()

#PCs
pcs.bulk <- prcomp(t(bulk.log))
pcs.bulkized <- prcomp(t(bulkized.log))

#Plot PCs
png("results/figures/bulk.hist.PCs.png")
plot(pcs.bulk, main="PCs of Bulk")
dev.off()

png("results/figures/bulkized.hist.PCs.png")
plot(pcs.bulkized, main="PCs of Bulkized")
dev.off()

#Format metadata
meta <- bulk.meta.edit[,-1:-2]
rownames(meta) <- bulk.meta.edit$Sample.name

#Get PCs
x.bulk <- as.data.frame(pcs.bulk$x)
x.bulkized <- as.data.frame(pcs.bulkized$x)

#Combine 
prinComps.bulk <- cbind(meta, x.bulk)
prinComps.bulkized <- cbind(meta,x.bulkized)
prinComps.bulk <- prinComps.bulk[,1:13]
prinComps.bulkized <- prinComps.bulkized[,1:13]


#Format for plots
prinComps.bulk <- prinComps.bulk %>% mutate(PMI_group= case_when(PMI < 10 ~ "0-9",PMI < 15 ~ "10-14", PMI < 20 ~ "15-19", PMI < 25 ~ "20-24", PMI < 30 ~ "25-29", PMI < 30 ~ "25-29",
                                                                 PMI < 35 ~ "30-34",PMI < 40 ~ "35-49", PMI >= 40 ~ ">40"))
prinComps.bulkized <- prinComps.bulkized %>% mutate(PMI_group= case_when(PMI < 10 ~ "0-9",PMI < 15 ~ "10-14", PMI < 20 ~ "15-19", PMI < 25 ~ "20-24", PMI < 30 ~ "25-29", PMI < 30 ~ "25-29",
                                                                         PMI < 35 ~ "30-34",PMI < 40 ~ "35-49", PMI >= 40 ~ ">40"))

prinComps.bulk <- prinComps.bulk %>% mutate(RIN_group= case_when(RIN < 7 ~ "6.5-6.9",RIN < 7.5 ~ "7.0-7.4",RIN < 8 ~ "7.5-7.0",RIN < 8.5 ~ "8.0-8.4", RIN < 9 ~ "8.5-8.9", RIN == 9 ~ "9"))
prinComps.bulkized <- prinComps.bulkized %>% mutate(RIN_group= case_when(RIN < 7 ~ "6.5-6.9",RIN < 7.5 ~ "7.0-7.4",RIN < 8 ~ "7.5-7.0",RIN < 8.5 ~ "8.0-8.4", RIN < 9 ~ "8.5-8.9", RIN == 9 ~ "9"))


prinComps.bulk <- prinComps.bulk %>% mutate(age_group= case_when(Age < 5 ~ "<5", Age < 10 ~ "5-9" ,Age < 15 ~ "10-14", Age < 20 ~ "15-19", Age < 25 ~ "20-22"))
prinComps.bulkized <- prinComps.bulkized %>% mutate(age_group= case_when(Age < 5 ~ "<5", Age < 10 ~ "5-9" ,Age < 15 ~ "10-14", Age < 20 ~ "15-19", Age < 25 ~ "20-22"))
                                                                         

prinComps.bulk <- prinComps.bulk %>% mutate(region_2= case_when(Region == "ACC" ~ "ACC", Region == "PFC" ~ "PFC", Region == "BA24" ~ "ACC"))
prinComps.bulkized <- prinComps.bulkized %>% mutate(region_2= case_when(Region == "ACC" ~ "ACC", Region == "PFC" ~ "PFC",Region == "BA24" ~ "ACC"))
                                                                        

#Plot 
ggplot(prinComps.bulk, aes(PC1,PC2, color=region_2)) +
  geom_point() +
  labs(title="", x= "PC1", y="PC2", color="Region") +
  theme_bw()
ggsave("results/figures/PCA_bulk_region.png")

ggplot(prinComps.bulkized, aes(PC1,PC2, color=region_2)) +
  geom_point() +
  labs(title="", x= "PC1", y="PC2",color="Region") +
  theme_bw()
ggsave("results/figures/PCA_bulkized_region.png")


prinComps.bulk <- prinComps.bulk %>% mutate(age_group= case_when(Age < 5 ~ "<5", Age < 10 ~ "5-9" ,Age < 15 ~ "10-14", Age < 20 ~ "15-19", Age < 25 ~ "20-22"))
prinComps.bulkized <- prinComps.bulkized %>% mutate(age_group= case_when(Age < 5 ~ "<5", Age < 10 ~ "5-9" ,Age < 15 ~ "10-14", Age < 20 ~ "15-19", Age < 25 ~ "20-22"))
                                                                         

prinComps.bulk <- prinComps.bulk %>% mutate(region_2= case_when(Region == "ACC" ~ "ACC", Region == "PFC" ~ "PFC", Region == "BA24" ~ "ACC"))
prinComps.bulkized <- prinComps.bulkized %>% mutate(region_2= case_when(Region == "ACC" ~ "ACC", Region == "PFC" ~ "PFC",Region == "BA24" ~ "ACC"))
                                                                        

#Plot 
ggplot(prinComps.bulk, aes(PC1,PC2, color=region_2)) +
  geom_point() +
  labs(title="", x= "PC1", y="PC2", color="Region") +
  theme_bw()
ggsave("results/figures/PCA_bulk_region.png")

ggplot(prinComps.bulkized, aes(PC1,PC2, color=region_2)) +
  geom_point() +
  labs(title="", x= "PC1", y="PC2",color="Region") +
  theme_bw()
ggsave("results/figures/PCA_bulkized_region.png")
