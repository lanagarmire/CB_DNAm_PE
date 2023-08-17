# # R/4.1.2 on server:
# ## open:
# source /home/yhdu/anaconda3/bin/activate 
# conda init 
# conda activate r4
# R
# ## Quit:
# q() 
# conda deactivate

# Library:
library(limma)
library("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
library(dbplyr)
library(ggplot2)
library(gtools)
library(tidyverse)
library(stringi)
library("stringr")

# Shortcut
load("/home/liuwent/04-Full_Model/fitConfounders.RData")
load("/home/liuwent/04d-DM_analysis_after_adj/sig_DMP.RData")
load("/home/liuwent/04-Full_Model/Mvalues.RData")

# split the gene names
for(i in 1:length(sig_DMP$UCSC_RefGene_Name)){sig_DMP$UCSC_RefGene_Name_new[i] = strsplit(sig_DMP$UCSC_RefGene_Name[i], ";")[[1]][1]}
sig_DMP$UCSC_RefGene_Name_new

# split the gene features
for(i in 1:length(sig_DMP$UCSC_RefGene_Group)){sig_DMP$UCSC_RefGene_Group_new[i] = strsplit(sig_DMP$UCSC_RefGene_Group[i], ";")[[1]][1]}
sig_DMP$UCSC_RefGene_Group_new

# reorder chr
num = '0'
for(i in 1:length(sig_DMP$chr)){
  for(j in c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9')){
    if(sig_DMP$chr[i] == j){
      stri_sub(sig_DMP$chr[i], 4, 3) <- num
    }
  }
}

# write.csv(sig_DMP, file="sig_DMP.csv")
# save(sig_DMP, file ="sig_DMP.RData")
# load("/home/liuwent/04d-DM_analysis_after_adj/sig_DMP.RData")

## differential methylation on probe-level (DMP)
### Plot significant cpg distribution in chr, cgi(Relation_to_Island) and feature(UCSC_RefGene_Group)
annEPIC = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
table_fit = topTable(fitConfounders , num=Inf, coef=2, genelist=annEPIC, adjust.method = "BH")
sig_DMP <- subset(table_fit, adj.P.Val < 0.05)
sig_DMP$direction = ifelse(sig_DMP$logFC>0, "Hyper", "Hypo")

pdf(file = "/home/liuwent/04d-DM_analysis_after_adj/sig_DMP.pdf")
ggplot(sig_DMP, aes(x = chr, fill = direction)) + 
  geom_bar(stat = "count", position = "dodge") + 
  stat_count(geom = "text", size = 3.5, 
             aes(label = ..count..), position=position_dodge(width=0.9), vjust=-0.25) + 
  scale_x_discrete(guide = guide_axis(n.dodge=2))

ggplot(sig_DMP, aes(x = Relation_to_Island, fill = direction)) +
  geom_bar(stat = "count", position = "dodge") +
  stat_count(geom = "text", size = 3.5, 
             aes(label = ..count..), position=position_dodge(width=0.9), vjust=-0.25)

ggplot(sig_DMP, aes(x = UCSC_RefGene_Group_new, fill = direction)) +
  geom_bar(stat = "count", position = "dodge") +
  stat_count(geom = "text", size = 3.5, 
             aes(label = ..count..), position=position_dodge(width=0.9), vjust=-0.25)
dev.off()

## differential methylation on gene-level
### GSEA:
#### create Expression Datasets for GSEA:
exp_data_cpg <- Mvalues[rownames(Mvalues)%in%rownames(sig_DMP),]
exp_data_cpg <- exp_data_cpg[order(row.names(exp_data_cpg)), ]

#### extract gene names from sig_DMP:
gene_names <- sig_DMP['UCSC_RefGene_Name_new']
gene_names <- gene_names[order(row.names(gene_names)), ]

#### merge gene names into exp_data_cpg:
exp_data <- cbind(gene_names, exp_data_cpg)
write.csv(exp_data, file="exp_data.csv")

#### create Preranked dataset for GSEA:
logFC <- sig_DMP['logFC']
logFC <- as.numeric(logFC[order(row.names(logFC)), ])
GSEA_prernk <- cbind.data.frame(gene_names, logFC)
GSEA_prernk <- GSEA_prernk[order(logFC, decreasing=TRUE),]
write.csv(GSEA_prernk, file="GSEA_prernk.csv")

### missmethyl:
library(missMethyl)
sig = rownames(sig_DMP)
all = rownames(annEPIC)

pdf("/home/liuwent/04d-DM_analysis_after_adj/gometh.pdf")
par(mfrow=c(1,1))
gst = gometh(sig.cpg = sig, collection = "GO", array.type="EPIC", 
             plot.bias=TRUE, prior.prob=TRUE)
dev.off()
# Top 10 GO categories
topGSA(gst, number=10)

# make sure to use limma version >= 3.52.2, because KEGG API moving to HTTPS as of June 1, 2022.
pdf("/home/liuwent/04d-DM_analysis_after_adj/gometh_kegg.pdf")
par(mfrow=c(1,1))
kegg <- gometh(sig.cpg = sig, collection = "KEGG", array.type="EPIC", 
               plot.bias=TRUE, prior.prob=TRUE)
dev.off()
# Top 10 KEGG categories
topGSA(kegg, number=10)







