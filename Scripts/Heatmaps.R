
library(reshape2)
library(dplyr)
library(pheatmap)

#load-------------------------------------------------------
counts.merged.long <- readRDS('data/processed/counts.merged.long.rds')
#Format
counts.bulk <- counts.merged.long[,1:3] %>%  
  dcast(gene~sample, sum, value.var = 'bulk') 

counts.bulkized <- counts.merged.long[,-3:-4] %>%
  dcast(gene~sample, sum, value.var = 'bulkized.norm')

counts.bulk.log <- log2(1+ counts.bulk[,-1])
counts.bulkized.log <- log2(1+ counts.bulkized[,-1])

rownames(counts.bulk.log) <- counts.bulk$gene
rownames(counts.bulkized.log) <-counts.bulkized$gene

#Format metadata
bulk.meta <-read.csv("data/raw/meta.bulk.edited.csv")
sc.meta <- read.csv("data/raw/meta.singlecell.csv")

bulk.meta <- bulk.meta %>% mutate(region_2 = case_when(Region == "BA24" ~ "ACC",Region == "PFC" ~ "PFC",Region == "ACC" ~ "ACC"))
bulk.meta$Region <- bulk.meta$region_2
meta.combined <- inner_join(bulk.meta, sc.meta, by=c("Individual"="individual","Region"="region")) %>% 
  select(Sample.name, Region, Age, Diagnosis, PMI, Sex, Seqbatch, Capbatch) %>% unique()



#Make annotations---------------------------------------------
annotation_bulkized <- data.frame(row.names=meta.combined$Sample.name, Region=meta.combined$Region, Age=meta.combined$Age,Sex=meta.combined$Sex, Diagnosis=meta.combined$Diagnosis,
                                  PMI=meta.combined$PMI, Seqbatch= meta.combined$Seqbatch,Capbatch=meta.combined$Capbatch)

annotation_bulkized <- annotation_bulkized %>% mutate(age_group= case_when(Age < 5 ~ "0-4", Age < 10 ~ "5-9", Age < 15 ~ "10-14" ,Age < 20 ~ "15-19",  Age < 20 ~ "15-19",  Age >= 20 ~ "20-25"))
annotation_bulkized <- annotation_bulkized %>% mutate(PMI_group= case_when(PMI < 10 ~ "0-9", PMI < 15 ~ "10-14", PMI < 20 ~ "15-19", PMI < 25 ~ "20-24", PMI < 30 ~ "25-29", PMI < 35 ~ "30-34",
                                                                           PMI < 40 ~ "35-49", PMI >= 40 ~ ">40"))

annotation_bulkized <- annotation_bulkized %>% select(-Age, -PMI) %>% rename("Age" = "age_group",  "PMI" = "PMI_group")
row.names(annotation_bulkized) <- meta.combined$Sample.name

annotation_bulk <- annotation_bulkized %>% select(-Seqbatch, -Capbatch)

#Annotation colors
annotation_colors_bulkized <- list(
  Region=c(ACC="royalblue1", PFC="orchid3"), Age=c("0-4"="indianred2", '5-9'="deepskyblue2",'10-14'="plum3",'15-19'="springgreen3",'20-25'="green4"),
  PMI = c("10-14"="darkorchid4", "0-9"="azure1", "20-24"="palegreen3","30-34"="coral2","35-49"="dodgerblue2",">40"="khaki2", "15-19"="lightblue4","25-29"="darkorange3"),
  Diagnosis=c(ASD="lightpink", Control="darkslategray2"), Sex=c(Male="cyan3",Female="brown3"), Seqbatch=c(SB1="grey25",SB2="grey77", SB3="honeydew2"), 
  Capbatch=c(CB1="coral1",CB2="hotpink2",CB3="brown1",CB4="turquoise",CB5="darkseagreen",CB6="red",CB7="orchid3",CB8="palegreen3",CB9="honeydew2"))

annotation_colors_bulk <- list(
  Region=c(ACC="royalblue1", PFC="orchid3"), Age=c("0-4"="indianred2", '5-9'="deepskyblue2",'10-14'="plum3",'15-19'="springgreen3",'20-25'="green4"),
  PMI = c("10-14"="darkorchid4", "0-9"="azure1", "20-24"="palegreen3","30-34"="coral2","35-49"="dodgerblue2",">40"="khaki2", "15-19"="lightblue4","25-29"="darkorange3"),
  Diagnosis=c(ASD="lightpink", Control="darkslategray2"), Sex=c(Male="cyan3",Female="brown3"))

#Make heatmaps--------------------------------------------
#Correlations between samples on bulk:
png("results/figures/bulk_heatmap.png")
corr_bulk <- round(cor(counts.bulk[,-1]),2) %>% pheatmap(border_color=NA, cluster_rows = T, cluster_cols = T, fontsize=6, annotation_col = annotation_bulk, annotation_colors=annotation_colors_bulk) 
dev.off()

#Correlations between samples on bulkized:
png("results/figures/bulkized_heatmap.png")
corr_bulkized <- round(cor(counts.bulkized[,-1]),2) %>% pheatmap(border_color=NA, cluster_rows = T, cluster_cols = T, fontsize=6, annotation_col = annotation_bulkized, annotation_colors=annotation_colors_bulkized)
dev.off()

#LOGGED----------------------
#Correlations between samples on bulk:
png("results/figures/bulk_heatmap.log.png")
corr_bulk <- round(cor(counts.bulk.log),2) %>% pheatmap(border_color=NA, cluster_rows = T, cluster_cols = T, fontsize=6, annotation_col = annotation_bulk, annotation_colors=annotation_colors_bulk) 
dev.off()

#Correlations between samples on bulkized:
png("results/figures/bulkized_heatmap.log.png")
corr_bulkized <- round(cor(counts.bulkized.log),2) %>% pheatmap(border_color=NA, cluster_rows = T, cluster_cols = T, fontsize=6, annotation_col = annotation_bulkized, annotation_colors=annotation_colors_bulkized)
dev.off()

