## Library:
library(dplyr)

## Data:
load("/home/liuwent/04d-DM_analysis_after_adj/sig_DMP.RData")
sig_DMP_sub <- sig_DMP[, c('chr','Relation_to_Island','UCSC_RefGene_Group_new','UCSC_RefGene_Name_new',
                       'logFC','AveExpr','t','P.Value','adj.P.Val','B','direction')]

PPI<-read.delim(file="/home/liuwent/04d-DM_analysis_after_adj/BIOGRID-ORGANISM-Homo_sapiens-4.4.217.tab3.txt", 
                sep = "\t", header = T)

## PPI Analysis:
### Both genes in the pair:
ppi_pair_both = PPI%>%filter(Official.Symbol.Interactor.A%in%sig_DMP_sub$UCSC_RefGene_Name_new & 
                          Official.Symbol.Interactor.B%in%sig_DMP_sub$UCSC_RefGene_Name_new)
ppi_pair_both = ppi_pair_both[which(ppi_pair_both$Official.Symbol.Interactor.A != 
                                      ppi_pair_both$Official.Symbol.Interactor.B),]
rownames(ppi_pair_both) = paste0(ppi_pair_both[,1],"_",ppi_pair_both[,2])
ppi_pair_both <- ppi_pair_both[c("Official.Symbol.Interactor.A","Official.Symbol.Interactor.B")]
ppi_pair_both <- unique(ppi_pair_both) # 2
write.csv(ppi_pair_both, file='ppi_pair_both.csv')
save(ppi_pair_both, file='ppi_pair_both.RData')
load("/home/liuwent/04d-DM_analysis_after_adj/ppi_pair_both.RData")

### One of the genes in the pair:
#### one of the genes in Official.Symbol.Interactor.A:
ppi_pair_A <- PPI%>%filter(Official.Symbol.Interactor.A%in%sig_DMP_sub$UCSC_RefGene_Name_new)
ppi_pair_A <- ppi_pair_A[which(ppi_pair_A$Official.Symbol.Interactor.A != 
                                 ppi_pair_A$Official.Symbol.Interactor.B),]
rownames(ppi_pair_A) = paste0(ppi_pair_A[,1],"_",ppi_pair_A[,2])
ppi_pair_A <- ppi_pair_A[c("Official.Symbol.Interactor.A","Official.Symbol.Interactor.B")]
ppi_pair_A <- unique(ppi_pair_A) # 957
write.csv(ppi_pair_A, file='ppi_pair_A.csv')
save(ppi_pair_A, file='ppi_pair_A.RData')
load("/home/liuwent/04d-DM_analysis_after_adj/ppi_pair_A.RData")

#### one of the genes in Official.Symbol.Interactor.B:
ppi_pair_B <- PPI%>%filter(Official.Symbol.Interactor.B%in%sig_DMP_sub$UCSC_RefGene_Name_new)
ppi_pair_B <- ppi_pair_B[which(ppi_pair_B$Official.Symbol.Interactor.A != 
                                 ppi_pair_B$Official.Symbol.Interactor.B),]
rownames(ppi_pair_B) = paste0(ppi_pair_B[,1],"_",ppi_pair_B[,2])
ppi_pair_B <- ppi_pair_B[c("Official.Symbol.Interactor.A","Official.Symbol.Interactor.B")]
ppi_pair_B <- unique(ppi_pair_B) # 1849
write.csv(ppi_pair_B, file='ppi_pair_B.csv')
save(ppi_pair_B, file='ppi_pair_B.RData')
load("/home/liuwent/04d-DM_analysis_after_adj/ppi_pair_B.RData")

## Overlab of significang genes and hub genes by using MCC:
ppi_pairs <- rbind(ppi_pair_both, ppi_pair_A, ppi_pair_B)
ppi_pairs <- unique(ppi_pairs)
write.csv(ppi_pairs, file='ppi_pairs.csv')
save(ppi_pairs, file='ppi_pairs.RData')
load("/home/liuwent/04d-DM_analysis_after_adj/ppi_pairs.RData")

load("/home/liuwent/04d-DM_analysis_after_adj/sig_DMP.RData")
diff_cpg_gene <- sig_DMP[, c(19, 55, 54, 47:53)]
colnames(diff_cpg_gene)[1:3] <- c('Island','Gene','Group')
colnames(diff_cpg_gene)[10] <- 'Type'
diff_gene = diff_cpg_gene[complete.cases(diff_cpg_gene),]
dim(diff_gene) # 41 10

ppi_pairs_rank <- read.csv(file="/home/liuwent/04d-DM_analysis_after_adj/ppi_pairs_MCC_top30.csv",header=TRUE)

sig <- diff_gene$Gene
write.csv(sig, file='sig.csv')

# # Edges (use Yuheng's panda_cobre code):
# ## Library:
# library(PANDA)
# 
# ## Data:
# load("/home/liuwent/04d-DM_analysis_after_adj/ppi_pair_both.RData")
# load("/home/liuwent/04d-DM_analysis_after_adj/ppi_pair_A.RData")
# load("/home/liuwent/04d-DM_analysis_after_adj/ppi_pair_B.RData")
# 
# ## Get the edge between two paired genes:
# ### For ppi_pair_both:
# ppis <- ppi_pair_both
# colnames(ppis) <- c('from_symbol', 'to_symbol')
# pandappi <- ppis
# Edges_both <- SignificantPairs(PPIdb = pandappi) #0
# 
# ### For ppi_pair_A:
# ppis <- ppi_pair_A
# colnames(ppis) <- c('from_symbol', 'to_symbol')
# pandappi <- ppis
# Edges_A <- SignificantPairs(PPIdb = pandappi) #0
# 
# ### For ppi_pair_B:
# ppis <- ppi_pair_B
# colnames(ppis) <- c('from_symbol', 'to_symbol')
# pandappi <- ppis
# Edges_B <- SignificantPairs(PPIdb = pandappi) #1
