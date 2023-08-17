# Library:
library(limma)
library("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
library(dbplyr)
library(lumi)

# Data:
load("/home/liuwent/04d-DM_analysis_after_adj/ppi_pairs.RData")
load("/home/liuwent/04-Full_Model/myCombat.RData")
load("/home/liuwent/04d-DM_analysis_after_adj/annEPIC_gene.RData")

# Gene to cpg:
annEPIC = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
## split the gene names
for(i in 1:length(annEPIC$UCSC_RefGene_Name)){annEPIC$UCSC_RefGene_Name_new[i] = 
  strsplit(annEPIC$UCSC_RefGene_Name[i], ";")[[1]][1]}
## split the gene features
for(i in 1:length(annEPIC$UCSC_RefGene_Group)){annEPIC$UCSC_RefGene_Group_new[i] = strsplit(annEPIC$UCSC_RefGene_Group[i], ";")[[1]][1]}
annEPIC$UCSC_RefGene_Group_new
save(annEPIC, file="annEPIC.RData")
# load("/home/liuwent/04d-DM_analysis_after_adj/annEPIC.RData")
annEPIC_gene <- annEPIC[, c('UCSC_RefGene_Name_new','UCSC_RefGene_Group_new'), drop=FALSE]

# Get top1 hub gene (RPTOR) module list:
top_in_A <- ppi_pairs[ppi_pairs$Official.Symbol.Interactor.A == 'RPTOR', ]
top_in_B <- ppi_pairs[ppi_pairs$Official.Symbol.Interactor.B == 'RPTOR', ]
top <- c('RPTOR', top_in_A[, 2], top_in_B[, 1])
top <- unique(top)
length(top) # 310
sum(top=='RPTOR')

# Get top2 hub gene (PDS5A) module list:
top_in_A <- ppi_pairs[ppi_pairs$Official.Symbol.Interactor.A == 'PDS5A', ]
to2_in_B <- ppi_pairs[ppi_pairs$Official.Symbol.Interactor.B == 'PDS5A', ]
top <- c('PDS5A', top_in_A[, 2], to2_in_B[, 1])
top <- unique(top)
length(top) # 240
sum(top=='PDS5A')

# Get top3 hub gene (ADAR) module list:
top_in_A <- ppi_pairs[ppi_pairs$Official.Symbol.Interactor.A == 'ADAR', ]
top_in_B <- ppi_pairs[ppi_pairs$Official.Symbol.Interactor.B == 'ADAR', ]
top <- c('ADAR', top_in_A[, 2], top_in_B[, 1])
top <- unique(top)
length(top) # 201
sum(top=='ADAR')

# Get top4 hub gene (SAFB2) module list:
top_in_A <- ppi_pairs[ppi_pairs$Official.Symbol.Interactor.A == 'SAFB2', ]
top_in_B <- ppi_pairs[ppi_pairs$Official.Symbol.Interactor.B == 'SAFB2', ]
top <- c(top_in_A[, 2], top_in_A[, 1])
top <- unique(top)
length(top) # 18
sum(top=='SAFB2')

# Get top5 hub gene (AMOTL2) module list:
top_in_A <- ppi_pairs[ppi_pairs$Official.Symbol.Interactor.A == 'AMOTL2', ]
top_in_B <- ppi_pairs[ppi_pairs$Official.Symbol.Interactor.B == 'AMOTL2', ]
top <- c(top_in_A[, 2], top_in_A[, 1])
top <- unique(top)
length(top) # 66
sum(top=='AMOTL2')

# Get top6 hub gene (AMOTL2) module list:
top_in_A <- ppi_pairs[ppi_pairs$Official.Symbol.Interactor.A == 'CRYAB', ]
top_in_B <- ppi_pairs[ppi_pairs$Official.Symbol.Interactor.B == 'CRYAB', ]
top <- c(top_in_A[, 2], top_in_A[, 1])
top <- unique(top)
length(top) # 79
sum(top=='CRYAB')

# Get top7 hub gene (CSTF2) module list:
top_in_A <- ppi_pairs[ppi_pairs$Official.Symbol.Interactor.A == 'CSTF2', ]
top_in_B <- ppi_pairs[ppi_pairs$Official.Symbol.Interactor.B == 'CSTF2', ]
top <- c(top_in_A[, 2], top_in_A[, 1])
top <- unique(top)
length(top) # 72
sum(top=='CSTF2')

# Get top8 hub gene (RPL37) module list:
top_in_A <- ppi_pairs[ppi_pairs$Official.Symbol.Interactor.A == 'RPL37', ]
top_in_B <- ppi_pairs[ppi_pairs$Official.Symbol.Interactor.B == 'RPL37', ]
top <- c(top_in_A[, 2], top_in_A[, 1])
top <- unique(top)
length(top) # 86
sum(top=='RPL37')

# Get top9 hub gene (MAD1L1) module list:
top_in_A <- ppi_pairs[ppi_pairs$Official.Symbol.Interactor.A == 'MAD1L1', ]
top_in_B <- ppi_pairs[ppi_pairs$Official.Symbol.Interactor.B == 'MAD1L1', ]
top <- c(top_in_A[, 2], top_in_A[, 1])
top <- unique(top)
length(top) # 61
sum(top=='MAD1L1')

# Get top10 hub gene (HLTF) module list:
top_in_A <- ppi_pairs[ppi_pairs$Official.Symbol.Interactor.A == 'HLTF', ]
top_in_B <- ppi_pairs[ppi_pairs$Official.Symbol.Interactor.B == 'HLTF', ]
top <- c(top_in_A[, 2], top_in_A[, 1])
top <- unique(top)
length(top) # 22
sum(top=='HLTF')

################################################################################
annEPIC_gene_top10 <- annEPIC_gene[which(annEPIC_gene$UCSC_RefGene_Name_new%in%top), , drop=FALSE]
colnames(annEPIC_gene_top10) <- c('CpG', 'Gene', 'Group')
dim(annEPIC_gene_top10)
save(annEPIC_gene_top10, file='annEPIC_gene_top10.RData')
load('/home/liuwent/04d-DM_analysis_after_adj/annEPIC_gene_top10.RData')
dim(annEPIC_gene_top10)

# Beta values:
myCombat_sub = myCombat[rownames(myCombat)[rownames(myCombat)%in%rownames(annEPIC_gene_top)],] 
dim(myCombat_sub)

# Take the geometric mean of all cpgs that related to the same gene
cpg_geomean_beta_top <- NULL
for(gene_name in unique(annEPIC_gene_top$Gene)){
  cpg = annEPIC_gene_top[annEPIC_gene_top$Gene==gene_name,]$CpG
  cpg_df = NULL
  if(sum(rownames(myCombat_sub)%in%cpg)!=0){
    cpg_sub = myCombat_sub[rownames(myCombat_sub)%in%cpg,]
    if(length(cpg_sub)==62){
      cpg_df = as.data.frame(t(cpg_sub))
      rownames(cpg_df) = c(gene_name)
    }else{
      cpg_df = as.data.frame(t(apply(cpg_sub, 2, function(x) exp(mean(log(x)))))) #geometric mean
      rownames(cpg_df) = c(gene_name)
    }
  }
  cpg_geomean_beta_top = rbind(cpg_geomean_beta_top, cpg_df)
}
cpg_geomean_mval_top = beta2m(cpg_geomean_beta_top) 
dim(cpg_geomean_mval_top)
sum(rownames(cpg_geomean_mval_top)=='RPTOR')
sum(rownames(cpg_geomean_mval_top)=='PDS5A')
sum(rownames(cpg_geomean_mval_top)=='ADAR')
sum(rownames(cpg_geomean_mval_top)=='SAFB2')
sum(rownames(cpg_geomean_mval_top)=='AMOTL2')
sum(rownames(cpg_geomean_mval_top)=='CRYAB')
sum(rownames(cpg_geomean_mval_top)=='CSTF2')
sum(rownames(cpg_geomean_mval_top)=='RPL37')
sum(rownames(cpg_geomean_mval_top)=='MAD1L1')
sum(rownames(cpg_geomean_mval_top)=='HLTF')

# GSEA:
## create Expression Datasets for GSEA:
exp_data_top10 <- cpg_geomean_mval_top
exp_data_top10 <- exp_data_top10[order(row.names(exp_data_top10)), ]
write.csv(exp_data_top10, file="exp_data_top10.csv")
save(exp_data_top10, file='exp_data_top10.RData')
load('/home/liuwent/04d-DM_analysis_after_adj/exp_data_top10.RData')