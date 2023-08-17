# Library:
library(tidyverse)
library(limma)
library(dplyr)

# Data:
load("/home/liuwent/04-Full_Model/fitConfounders.RData")
load('/home/liuwent/04d-DM_analysis_after_adj/annEPIC_gene_top1.RData')
load('/home/liuwent/04d-DM_analysis_after_adj/annEPIC_gene_top2.RData')
load('/home/liuwent/04d-DM_analysis_after_adj/annEPIC_gene_top3.RData')
load('/home/liuwent/04d-DM_analysis_after_adj/annEPIC_gene_top4.RData')
load('/home/liuwent/04d-DM_analysis_after_adj/annEPIC_gene_top5.RData')
load('/home/liuwent/04d-DM_analysis_after_adj/annEPIC_gene_top6.RData')
load('/home/liuwent/04d-DM_analysis_after_adj/annEPIC_gene_top7.RData')
load('/home/liuwent/04d-DM_analysis_after_adj/annEPIC_gene_top8.RData')
load('/home/liuwent/04d-DM_analysis_after_adj/annEPIC_gene_top9.RData')
load('/home/liuwent/04d-DM_analysis_after_adj/annEPIC_gene_top10.RData')

# Get logFC values:
table_fit = topTable(fitConfounders, num=Inf, coef=2, adjust.method = "BH")
table_fit$Type = ifelse(table_fit$logFC>0, "Hyper", "Hypo")

################################################################################
top1_RPTOR <- merge.data.frame(table_fit, annEPIC_gene_top1, by = 'row.names')
top1_RPTOR <- top1_RPTOR[c('CpG', 'Gene', 'Type')]
dim(top1_RPTOR) # 9931    3

Hyper <- top1_RPTOR %>% 
  group_by(Gene) %>% 
  summarise(Hyper = sum(Type=='Hyper')) 
Hypo <- top1_RPTOR %>% 
  group_by(Gene) %>% 
  summarise(Hyper = sum(Type=='Hypo'))
temp <- merge(Hyper, Hypo, by='Gene')
colnames(temp) <- c('Gene', 'Hyper', 'Hypo')

for(i in 1:nrow(temp)){
  if(temp$Hyper[i] >= temp$Hypo[i]){
    temp$Type[i] = 'Hyper'
  }else{temp$Type[i] = 'Hypo'}
}
top1_RPTOR_new <- temp
write.csv(top1_RPTOR_new, file='top1_RPTOR_new.csv')
read.csv('/home/liuwent/04d-DM_analysis_after_adj/top1_RPTOR_new.csv', header=TRUE)
save(top1_RPTOR_new, file='top1_RPTOR_new.RData')
load('/home/liuwent/04d-DM_analysis_after_adj/top1_RPTOR_new.RData')

################################################################################
top2_PDS5A <- merge.data.frame(table_fit, annEPIC_gene_top2, by = 'row.names')
top2_PDS5A <- top2_PDS5A[c('CpG', 'Gene', 'Type')]
dim(top2_PDS5A) # 6436    3

Hyper <- top2_PDS5A %>% 
  group_by(Gene) %>% 
  summarise(Hyper = sum(Type=='Hyper')) 
Hypo <- top2_PDS5A %>% 
  group_by(Gene) %>% 
  summarise(Hyper = sum(Type=='Hypo'))
temp <- merge(Hyper, Hypo, by='Gene')
colnames(temp) <- c('Gene', 'Hyper', 'Hypo')

for(i in 1:nrow(temp)){
  if(temp$Hyper[i] >= temp$Hypo[i]){
    temp$Type[i] = 'Hyper'
  }else{temp$Type[i] = 'Hypo'}
}
top2_PDS5A_new <- temp
write.csv(top2_PDS5A_new, file='top2_PDS5A_new.csv')
read.csv('/home/liuwent/04d-DM_analysis_after_adj/top2_PDS5A_new.csv', header=TRUE)
save(top2_PDS5A_new, file='top2_PDS5A_new.RData')
load('/home/liuwent/04d-DM_analysis_after_adj/top2_PDS5A_new.RData')

################################################################################
top3_ADAR <- merge.data.frame(table_fit, annEPIC_gene_top3, by = 'row.names')
top3_ADAR <- top3_ADAR[c('CpG', 'Gene', 'Type')]
dim(top3_ADAR) # 5575    3

Hyper <- top3_ADAR %>% 
  group_by(Gene) %>% 
  summarise(Hyper = sum(Type=='Hyper')) 
Hypo <- top3_ADAR %>% 
  group_by(Gene) %>% 
  summarise(Hyper = sum(Type=='Hypo'))
temp <- merge(Hyper, Hypo, by='Gene')
colnames(temp) <- c('Gene', 'Hyper', 'Hypo')

for(i in 1:nrow(temp)){
  if(temp$Hyper[i] >= temp$Hypo[i]){
    temp$Type[i] = 'Hyper'
  }else{temp$Type[i] = 'Hypo'}
}
top3_ADAR_new <- temp
write.csv(top3_ADAR_new, file='top3_ADAR_new.csv')
read.csv('/home/liuwent/04d-DM_analysis_after_adj/top3_ADAR_new.csv', header=TRUE)
save(top3_ADAR_new, file='top3_ADAR_new.RData')
load('/home/liuwent/04d-DM_analysis_after_adj/top3_ADAR_new.RData')

################################################################################
top4_SAFB2 <- merge.data.frame(table_fit, annEPIC_gene_top4, by = 'row.names')
top4_SAFB2 <- top4_SAFB2[c('CpG', 'Gene', 'Type')]
dim(top4_SAFB2) # 686    3

Hyper <- top4_SAFB2 %>% 
  group_by(Gene) %>% 
  summarise(Hyper = sum(Type=='Hyper')) 
Hypo <- top4_SAFB2 %>% 
  group_by(Gene) %>% 
  summarise(Hyper = sum(Type=='Hypo'))
temp <- merge(Hyper, Hypo, by='Gene')
colnames(temp) <- c('Gene', 'Hyper', 'Hypo')

for(i in 1:nrow(temp)){
  if(temp$Hyper[i] >= temp$Hypo[i]){
    temp$Type[i] = 'Hyper'
  }else{temp$Type[i] = 'Hypo'}
}
top4_SAFB2_new <- temp
write.csv(top4_SAFB2_new, file='top4_SAFB2_new.csv')
read.csv('/home/liuwent/04d-DM_analysis_after_adj/top4_SAFB2_new.csv', header=TRUE)
save(top4_SAFB2_new, file='top4_SAFB2_new.RData')
load('/home/liuwent/04d-DM_analysis_after_adj/top4_SAFB2_new.RData')

################################################################################
top5_AMOTL2 <- merge.data.frame(table_fit, annEPIC_gene_top5, by = 'row.names')
top5_AMOTL2 <- top5_AMOTL2[c('CpG', 'Gene', 'Type')]
dim(top5_AMOTL2) # 2624    3

Hyper <- top5_AMOTL2 %>% 
  group_by(Gene) %>% 
  summarise(Hyper = sum(Type=='Hyper')) 
Hypo <- top5_AMOTL2 %>% 
  group_by(Gene) %>% 
  summarise(Hyper = sum(Type=='Hypo'))
temp <- merge(Hyper, Hypo, by='Gene')
colnames(temp) <- c('Gene', 'Hyper', 'Hypo')

for(i in 1:nrow(temp)){
  if(temp$Hyper[i] >= temp$Hypo[i]){
    temp$Type[i] = 'Hyper'
  }else{temp$Type[i] = 'Hypo'}
}
top5_AMOTL2_new <- temp
write.csv(top5_AMOTL2_new, file='top5_AMOTL2_new.csv')
read.csv('/home/liuwent/04d-DM_analysis_after_adj/top5_AMOTL2_new.csv', header=TRUE)
save(top5_AMOTL2_new, file='top5_AMOTL2_new.RData')
load('/home/liuwent/04d-DM_analysis_after_adj/top5_AMOTL2_new.RData')

################################################################################
top6_CRYAB <- merge.data.frame(table_fit, annEPIC_gene_top6, by = 'row.names')
top6_CRYAB <- top6_CRYAB[c('CpG', 'Gene', 'Type')]
dim(top6_CRYAB) # 2569    3

Hyper <- top6_CRYAB %>% 
  group_by(Gene) %>% 
  summarise(Hyper = sum(Type=='Hyper')) 
Hypo <- top6_CRYAB %>% 
  group_by(Gene) %>% 
  summarise(Hyper = sum(Type=='Hypo'))
temp <- merge(Hyper, Hypo, by='Gene')
colnames(temp) <- c('Gene', 'Hyper', 'Hypo')

for(i in 1:nrow(temp)){
  if(temp$Hyper[i] >= temp$Hypo[i]){
    temp$Type[i] = 'Hyper'
  }else{temp$Type[i] = 'Hypo'}
}
top6_CRYAB_new <- temp
write.csv(top6_CRYAB_new, file='top6_CRYAB_new.csv')
read.csv('/home/liuwent/04d-DM_analysis_after_adj/top6_CRYAB_new.csv', header=TRUE)
save(top6_CRYAB_new, file='top6_CRYAB_new.RData')
load('/home/liuwent/04d-DM_analysis_after_adj/top6_CRYAB_new.RData')

################################################################################
top7_CSTF2 <- merge.data.frame(table_fit, annEPIC_gene_top7, by = 'row.names')
top7_CSTF2 <- top7_CSTF2[c('CpG', 'Gene', 'Type')]
dim(top7_CSTF2) # 1945    3

Hyper <- top7_CSTF2 %>% 
  group_by(Gene) %>% 
  summarise(Hyper = sum(Type=='Hyper')) 
Hypo <- top7_CSTF2 %>% 
  group_by(Gene) %>% 
  summarise(Hyper = sum(Type=='Hypo'))
temp <- merge(Hyper, Hypo, by='Gene')
colnames(temp) <- c('Gene', 'Hyper', 'Hypo')

for(i in 1:nrow(temp)){
  if(temp$Hyper[i] >= temp$Hypo[i]){
    temp$Type[i] = 'Hyper'
  }else{temp$Type[i] = 'Hypo'}
}
top7_CSTF2_new <- temp
write.csv(top7_CSTF2_new, file='top7_CSTF2_new.csv')
read.csv('/home/liuwent/04d-DM_analysis_after_adj/top7_CSTF2_new.csv', header=TRUE)
save(top7_CSTF2_new, file='top7_CSTF2_new.RData')
load('/home/liuwent/04d-DM_analysis_after_adj/top7_CSTF2_new.RData')

################################################################################
top8_RPL37 <- merge.data.frame(table_fit, annEPIC_gene_top8, by = 'row.names')
top8_RPL37 <- top8_RPL37[c('CpG', 'Gene', 'Type')]
dim(top8_RPL37) # 1744    3

Hyper <- top8_RPL37 %>% 
  group_by(Gene) %>% 
  summarise(Hyper = sum(Type=='Hyper')) 
Hypo <- top8_RPL37 %>% 
  group_by(Gene) %>% 
  summarise(Hyper = sum(Type=='Hypo'))
temp <- merge(Hyper, Hypo, by='Gene')
colnames(temp) <- c('Gene', 'Hyper', 'Hypo')

for(i in 1:nrow(temp)){
  if(temp$Hyper[i] >= temp$Hypo[i]){
    temp$Type[i] = 'Hyper'
  }else{temp$Type[i] = 'Hypo'}
}
top8_RPL37_new <- temp
write.csv(top8_RPL37_new, file='top8_RPL37_new.csv')
read.csv('/home/liuwent/04d-DM_analysis_after_adj/top8_RPL37_new.csv', header=TRUE)
save(top8_RPL37_new, file='top8_RPL37_new.RData')
load('/home/liuwent/04d-DM_analysis_after_adj/top8_RPL37_new.RData')

################################################################################
top9_MAD1L1 <- merge.data.frame(table_fit, annEPIC_gene_top9, by = 'row.names')
top9_MAD1L1 <- top9_MAD1L1[c('CpG', 'Gene', 'Type')]
dim(top9_MAD1L1) # 2733    3

Hyper <- top9_MAD1L1 %>% 
  group_by(Gene) %>% 
  summarise(Hyper = sum(Type=='Hyper')) 
Hypo <- top9_MAD1L1 %>% 
  group_by(Gene) %>% 
  summarise(Hyper = sum(Type=='Hypo'))
temp <- merge(Hyper, Hypo, by='Gene')
colnames(temp) <- c('Gene', 'Hyper', 'Hypo')

for(i in 1:nrow(temp)){
  if(temp$Hyper[i] >= temp$Hypo[i]){
    temp$Type[i] = 'Hyper'
  }else{temp$Type[i] = 'Hypo'}
}
top9_MAD1L1_new <- temp
write.csv(top9_MAD1L1_new, file='top9_MAD1L1_new.csv')
read.csv('/home/liuwent/04d-DM_analysis_after_adj/top9_MAD1L1_new.csv', header=TRUE)
save(top9_MAD1L1_new, file='top9_MAD1L1_new.RData')
load('/home/liuwent/04d-DM_analysis_after_adj/top9_MAD1L1_new.RData')

################################################################################
top10_HLTF <- merge.data.frame(table_fit, annEPIC_gene_top10, by = 'row.names')
top10_HLTF <- top10_HLTF[c('CpG', 'Gene', 'Type')]
dim(top10_HLTF) # 513    3

Hyper <- top10_HLTF %>% 
  group_by(Gene) %>% 
  summarise(Hyper = sum(Type=='Hyper')) 
Hypo <- top10_HLTF %>% 
  group_by(Gene) %>% 
  summarise(Hyper = sum(Type=='Hypo'))
temp <- merge(Hyper, Hypo, by='Gene')
colnames(temp) <- c('Gene', 'Hyper', 'Hypo')

for(i in 1:nrow(temp)){
  if(temp$Hyper[i] >= temp$Hypo[i]){
    temp$Type[i] = 'Hyper'
  }else{temp$Type[i] = 'Hypo'}
}
top10_HLTF_new <- temp
write.csv(top10_HLTF_new, file='top10_HLTF_new.csv')
read.csv('/home/liuwent/04d-DM_analysis_after_adj/top10_HLTF_new.csv', header=TRUE)
save(top10_HLTF_new, file='top10_HLTF_new.RData')
load('/home/liuwent/04d-DM_analysis_after_adj/top10_HLTF_new.RData')

### Network plots
library(zoo)
library(igraph)
library(repr)

multi_full = read.csv('/home/liuwent/04d-DM_analysis_after_adj/top9_MAD1L1_new.csv', header=TRUE)
head(multi_full)
multi_full$value = 1
multi_full$hub = "MAD1L1"

rptor_df = with(multi_full, tapply(value, list(hub, Gene), FUN = sum ))
dim(rptor_df)
rptor_df <- replace(rptor_df, is.na(rptor_df), 0)
(rptor_df)

# build the graph object
network <- graph_from_incidence_matrix(rptor_df, weighted = TRUE)

# Change the color of th vertices according to index
V(network)$color <- c(rep("yellow",1), rep("orange", 61))
V(network)[1]$color = 'yellow'
for (i in 1:length(V(network))){
  # change the vertices according to hypo/hyper
  if ("Hypo" %in% multi_full$Type[match(names(V(network)[i]),multi_full$Gene)]){
    V(network)$color[i] =  "cadetblue1"
  }
  if ("Hyper" %in% multi_full$Type[match(names(V(network)[i]),multi_full$Gene)]){
    V(network)$color[i] =  "red"
  }
}
V(network)[1]$color = 'yellow'

jpeg("/home/liuwent/04d-DM_analysis_after_adj/top9_MAD1L1.jpg", width = 10000, height = 10000, res=300)
options(repr.plot.width=40,  repr.plot.height=40)

# change the size and font according to index
plot(network, mark.expand = 15, edge.width= 6, vertex.size= c(rep(20,1), rep(7, 296)), vertex.label.font = c(rep(2,1), rep(3, 296)), vertex.label.cex  = c(rep(2,182), rep(2, 17)),vertex.label.color  = c(rep("black",1), rep("dark blue", 296)),vertex.frame.width= 0, 
     edge.color= "gray40")
# add legend
legend("topright", legend = c("Gene(Hyper)","Gene(Hypo)"), pch=c(21,21),
       col=c("red","cadetblue1"), pt.bg=c("red","cadetblue1"), pt.cex=4, cex=4,bty="n", ncol=1)
dev.off()