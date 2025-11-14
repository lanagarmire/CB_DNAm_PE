##################################################
# Preprocessing and analysis for Herzog et al. data
###################################################

#include code for figure 3H, 3I

#Herzog data (whole cord blood)
library(limma)
library(GEOquery)
library(dplyr)
library(ggplot2)
library(lumi)
library(EpiDISH)
library(tidyr)


# Load Herzog(EOLO) data from NCBI:
gset <- getGEO("GSE103253", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GSE103253", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
#
pd_EOLO_450K <- pData(gset)
dim(pd_EOLO_450K)   # 245  33
save(pd_EOLO_450K, file = "./pd_EOLO_450K.RData")
beta_EOLO_450K <- exprs(gset)
dim(beta_EOLO_450K)   # 482750    245
save(beta_EOLO_450K, file = "./beta_EOLO_450K.RData")

## only extract whole cord blood samples
pd_EOLO <- pd_EOLO_450K[pd_EOLO_450K$tissue == "UC-WBC", ]
dim(pd_EOLO) # 88 33
save(pd_EOLO, file = "./pd_EOLO.RData")
beta_EOLO <- beta_EOLO_450K[, colnames(beta_EOLO_450K) %in% rownames(pd_EOLO)]
dim(beta_EOLO) # 482750     88
save(beta_EOLO, file = "./beta_EOLO.RData")

### EOPE&LOPE -- Select preeclampsia related data only
pd_EOLO <- pd_EOLO[pd_EOLO$group %in% c("Early-onset preeclampsia", "Late-onset preeclampsia", "Uncomplicated control"), ]
dim(pd_EOLO) # 48 33
pd_EOLO$group <- ifelse(pd_EOLO$group %in% c("Early-onset preeclampsia", "Late-onset preeclampsia"), "Case", "Control")

pd_EOLO$Slide <- as.character(pd_EOLO$Slide)
save(pd_EOLO, file = "./pd_EOLO.RData")
beta_EOLO <- beta_EOLO_450K[, colnames(beta_EOLO_450K) %in% rownames(pd_EOLO)]
dim(beta_EOLO) # 482750     88
save(beta_EOLO, file = "./beta_EOLO.RData")


## EOLO: regress only on sample group-----------------------------
load("/nfs/dcmb-lgarmire/liuwent/15-Additional_Datasets/pd_EOLO.RData")
load("/nfs/dcmb-lgarmire/liuwent/15-Additional_Datasets/beta_EOLO.RData")
load("/nfs/dcmb-lgarmire/xtyang/CordBlood/multi_hit_new_450k.RData")

pd_EOLO <- pd_EOLO %>%
  separate(title, into = c("Slide", "Array"), sep = "_")

# use the new cross-hybridization reference by Zhou et al..2017
dim(beta_EOLO)#482750     48
beta_EOLO = beta_EOLO[!rownames(beta_EOLO)%in%multi_hit_new_450k$probeID,]
dim(beta_EOLO)#465472     48
champ.SVD(beta = as.data.frame(beta_EOLO), pd = pd_EOLO, 
          resultsDir = "/nfs/dcmb-lgarmire/xtyang/CordBlood/03-wenting_results/03-EOLO_data/CHAMP_SVDimages/")


pd_EOLO$group = relevel(factor(pd_EOLO$group), ref="Control")
design0 = model.matrix(~group, data = pd_EOLO)
fit0 = lmFit(beta2m(beta_EOLO), design0)
fit00 = eBayes(fit0)

## Get log_2 fold change and p value
log2fc = fit00$coefficients[,2]
pval0 = p.adjust(fit00$p.value[,2], "BH")
sum(pval0<0.05)
sig_list = pval0[pval0<0.05]
write.csv(names(pval<0.05))
## Make a volcano plot (label DMP with different colors)
col_list = sapply(seq(1, dim(fit00$coefficients)[1]), function(i) {
  if (pval0[i]<0.05) {return("red")}
  else {return("black")} })
table(col_list)

png(file = "/nfs/dcmb-lgarmire/xtyang/CordBlood/03-wenting_results/03-EOLO_data/EOLO_volcano3H.png")
volcanoplot(fit00, coef=2, col=col_list, highlight=0,
            names=rownames(fit00$coefficients))
# abline(h=-log10(0.05), lty=2, col="red")
abline(v=-1, lty=2, col="grey")
abline(v=1, lty=2, col="grey")
text(-1, 5, "x=-1", cex=1.2)
text(1, 5, "x=1", cex=1.2)
dev.off()

# Volcano plot with adjustment (clinical variables and cell types):
## PBMC: no significant cpg, probably because the data is unbalanced, 20 cases vs 90 controls.



## EOLO: no clinical variables, only can adjust for cell types.
### Cell type deconvolution:
library(EpiDISH)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(ggsignif)

load("/nfs/dcmb-lgarmire/xtyang/CordBlood/03-wenting_results/FlowSorted.CordBloodCombined.450k.compTable.rda")
ref = FlowSorted.CordBloodCombined.450k.compTable
out<-epidish(beta.m=beta_EOLO, ref.m=as.matrix(ref), method="CP")
pd_EOLO = data.frame(cbind(out$estF, pd_EOLO))

library(compositions)
library(zCompositions)
cell_type = out$estF
cell_type[cell_type < 0] <- 0          # guard against negatives
cell_type_zr <- cmultRepl(cell_type, method = "CZM", output = "prop")  # compositional zero replacement
cell_ilr <- ilr(acomp(cell_type))
colnames(cell_ilr) = c("ct1", "ct2", "ct3", "ct4", "ct5","ct6")

allConfounders = data.frame(pd_EOLO[,c("group")], cell_ilr)
colnames(allConfounders)[1] = "Sample_Group"
formula = paste0(names(allConfounders), collapse = ' + ')
formula = paste0("~", formula)
formula = formula(formula)
#design matrix
design2 = model.matrix(formula, data = allConfounders)
colnames(design2)[2] = "Cases"
design2[,2] = 1-design2[,2]
fitConfounders = lmFit(beta_EOLO, design2)
fitConfounders = eBayes(fitConfounders)
#save(fitConfounders, file = "fitConfounders_EOLO.RData")

t_stats <- fitConfounders$t[,"Cases"] #random shuffled group
df<- fitConfounders$df.total[1]
z_score2 = zscoreT(t_stats, df, approx=FALSE, method = "bailey")
bc2 <- bacon(z_score2)
limma_inflation_score2 = inflation(bc2)#0.957

#make a toptable on second column sample_group
allg.limma <- topTable(fitConfounders, coef="Cases", n=dim(fitConfounders)[1], adjust.method="fdr")
sigg.limma <- subset(allg.limma, adj.P.Val < 0.05)


#make a toptable on second column sample_group
col_list2 = sapply(seq(1, dim(fitConfounders$coefficients)[1]), function(i) {
  if (allg.limma$adj.P.Val[i]<0.05) {return("red")}
  else {return("black")} })

png("/nfs/dcmb-lgarmire/xtyang/CordBlood/03-wenting_results/03-EOLO_data/volcano_3I.png")
par(cex.axis = 1.5, cex.lab = 1.5)

volcanoplot(fitConfounders,coef=2,col=col_list2,highlight=0,
            names=rownames(fitConfounders$coefficients))
# abline(h=-log10(0.05), lty=2, col="red")
abline(v=-1, lty=2, col="grey")
abline(v=1, lty=2, col="grey")
text(-1, 5, "x=-1", cex=1.2)
text(1, 5, "x=1", cex=1.2)
dev.off()