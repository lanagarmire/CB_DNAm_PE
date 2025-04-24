Library:
library(GEOquery)
library(dplyr)
library(ggplot2)
library(lumi)
library(EpiDISH)

#Load PBMC data from NCBI:
gset <- getGEO("GSE110828", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL13534", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
#
pd_PBMC_450K <- pData(gset)
dim(pd_PBMC_450K)   # 157  82
save(pd_PBMC_450K, file = "./pd_PBMC_450K.RData")
beta_PBMC_450K <- exprs(gset)
dim(beta_PBMC_450K)   # 410735    157
save(beta_PBMC_450K, file = "./beta_PBMC_450K.RData")

# # Load EOPE&LOPE data from NCBI:
gset <- getGEO("GSE103253", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GSE103253", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
#
pd_EOLO_450K <- pData(gset)
dim(pd_EOLO_450K)   #245  33
save(pd_EOLO_450K, file = "./pd_EOLO_450K.RData")
beta_EOLO_450K <- exprs(gset)
dim(beta_EOLO_450K)   # 482750    245
save(beta_EOLO_450K, file = "./beta_EOLO_450K.RData")

#Subset cord blood related samples only for the two datasets:
# PBMC
pd_PBMC <- pd_PBMC_450K[pd_PBMC_450K$source_name_ch1 == "CBMC", ]
dim(pd_PBMC) # 110  82
pd_PBMC$preeclampsia <- ifelse(pd_PBMC$preeclampsia == "Yes", "Case", "Control")
save(pd_PBMC, file = "./pd_PBMC.RData")
beta_PBMC <- beta_PBMC_450K[, colnames(beta_PBMC_450K) %in% rownames(pd_PBMC)]
dim(beta_PBMC) # 410735    110
save(beta_PBMC, file = "./beta_PBMC.RData")

## EOPE&LOPE
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
save(pd_EOLO, file = "./pd_EOLO.RData")
beta_EOLO <- beta_EOLO_450K[, colnames(beta_EOLO_450K) %in% rownames(pd_EOLO)]
dim(beta_EOLO) # 482750     88
save(beta_EOLO, file = "./beta_EOLO.RData")


# Volcano plot without adjustment:
## PBMC:regress only on sample group-----------------------------
pd_PBMC$preeclampsia = relevel(factor(pd_PBMC$preeclampsia), ref="Control")
design0 = model.matrix(~pd_PBMC$preeclampsia)
fit0 = lmFit(beta2m(beta_PBMC), design0)#use Mvalues
fit00 = eBayes(fit0)
# save(fit00, file = "fit00.RData")

## Get log_2 fold change and p value
log2fc = fit00$coefficients[,2]
##pval = p.adjust(fit$p.value[,2], method="bonferroni")
pval0 = p.adjust(fit00$p.value[,2], "BH")
length(pval0[pval0<0.05]) # 0

cpg <- data.frame(logFC=log2fc, pval=pval0)
sig_hyper <- rownames(cpg[cpg$logFC>0 & cpg$pval<0.05,])
length(sig_hyper) # 0
# save(sig_hyper, file='sig_hyper.RData')
sig_hypo <- rownames(cpg[cpg$logFC<0 & cpg$pval<0.05,])
# save(sig_hypo, file='sig_hypo.RData')
length(sig_hypo) # 0

table_fit00 = topTable(fit00, num=Inf, coef=2, adjust.method = "BH")
table_fit00$Type = ifelse(table_fit00$logFC>0, "Hyper", "Hypo")
sum(table_fit00$adj.P.Val<0.05) # 0

## Make a volcano plot (label DMP with different colors)
col_list = sapply(seq(1, dim(fit00$coefficients)[1]), function(i) {
  if (pval0[i]<0.05) {return("red")}
  else {return("black")} })
table(col_list)

pdf(file = "./PBMC_fit00.pdf")
par(cex.axis = 1.5, cex.lab = 1.5) 

volcanoplot(fit00, coef=2, col=col_list, highlight=5,
            names=rownames(fit00$coefficients))
# abline(h=-log10(0.05), lty=2, col="red")
abline(v=-1, lty=2, col="grey")
abline(v=1, lty=2, col="grey")
text(-1, 5, "x=-1", cex=1.2)
text(1, 5, "x=1", cex=1.2)
dev.off()


## EOLO: regress only on sample group-----------------------------
pd_EOLO$group = relevel(factor(pd_EOLO$group), ref="Control")
design0 = model.matrix(~pd_EOLO$group)
fit0 = lmFit(beta2m(beta_EOLO), design0)
fit00 = eBayes(fit0)
# save(fit00, file = "fit00.RData")

## Get log_2 fold change and p value
log2fc = fit00$coefficients[,2]
##pval = p.adjust(fit$p.value[,2], method="bonferroni")
pval0 = p.adjust(fit00$p.value[,2], "BH")
length(pval0[pval0<0.05]) # 24597

cpg <- data.frame(logFC=log2fc, pval=pval0)
sig_hyper <- rownames(cpg[cpg$logFC>0 & cpg$pval<0.05,])
length(sig_hyper) # 12466
# save(sig_hyper, file='sig_hyper.RData')
sig_hypo <- rownames(cpg[cpg$logFC<0 & cpg$pval<0.05,])
# save(sig_hypo, file='sig_hypo.RData')
length(sig_hypo) # 12131

table_fit00 = topTable(fit00, num=Inf, coef=2, adjust.method = "BH")
table_fit00$Type = ifelse(table_fit00$logFC>0, "Hyper", "Hypo")
sum(table_fit00$adj.P.Val<0.05) # 24597

## Make a volcano plot (label DMP with different colors)
col_list = sapply(seq(1, dim(fit00$coefficients)[1]), function(i) {
  if (pval0[i]<0.05) {return("red")}
  else {return("black")} })
table(col_list)

pdf(file = "./EOLO_fit00.pdf")
par(cex.axis = 1.5, cex.lab = 1.5) 

volcanoplot(fit00, coef=2, col=col_list, highlight=5,
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

#load("/home/liuwent/04b-cell_type_deconvolution/newnewSetofMarkers_final.RData")
ref = FlowSorted.CordBloodCombined.450k.compTable
out<-epidish(beta.m=beta_EOLO, ref.m=as.matrix(ref), method="CP")
count <- out$estF
pd_EOLO = data.frame(cbind(out$estF, pd_EOLO))


### Volcano plot with adjustment (for cell types) using Gervin ref
# Library:

allConfounders = pd_EOLO[,c("group", colnames(pd_EOLO)[1:7])]
formula = paste0(names(allConfounders), collapse = ' + ')
formula = paste0("~", formula)
formula = formula(formula)
#design matrix
design2 = model.matrix(formula, data = allConfounders)
colnames(design2)[2] = "Cases"
design2[,2] = 1-design2[,2]
#fit linear model
# fitConfounders = lmFit(Mvalues[, which(!is.na(allConfounders$BMI))], design2)
fitConfounders = lmFit(beta_EOLO, design2)
fitConfounders = eBayes(fitConfounders)
save(fitConfounders, file = "fitConfounders_EOLO.RData")

#make a toptable on second column sample_group
allg.limma <- topTable(fitConfounders, coef="Cases", n=dim(fitConfounders)[1], adjust.method="fdr")
sigg.limma <- subset(allg.limma, adj.P.Val < 0.05)
save(sigg.limma, file='sigg.limma.RData')
nonsigg.limma <- subset(allg.limma, adj.P.Val >= 0.05)
dim(sigg.limma)[1]

# Make a volcano plot (label DMP with different colors)
pval2 = p.adjust(fitConfounders$p.value[,2], method = "fdr")
col_list2 = sapply(seq(1, dim(fitConfounders$coefficients)[1]), function(i) {
  if (pval2[i]<0.05) {return("red")}
  else {return("black")} })
cpg <- data.frame(logFC=fitConfounders$coefficients[,2],pval=pval2)
sig_hyper <- rownames(cpg[cpg$logFC>0 & cpg$pval<0.05,])
length(sig_hyper)
# save(sig_hyper, file='sig_hyper.RData')
sig_hypo <- rownames(cpg[cpg$logFC<0 & cpg$pval<0.05,])
# save(sig_hypo, file='sig_hypo.RData')
length(sig_hypo)

png(file = "./EOLO-fitConfounders_cvct_large.png")
par(cex.axis = 1.5, cex.lab = 1.5)

volcanoplot(fitConfounders,coef=2,col=col_list2,highlight=1,
            names=rownames(fitConfounders$coefficients), xlim = c(-5,5))
# abline(h=-log10(0.05), lty=2, col="red")
abline(v=-1, lty=2, col="grey")
abline(v=1, lty=2, col="grey")
text(-1, 5, "x=-1", cex=1.2)
text(1, 5, "x=1", cex=1.2)
dev.off()

length(pval2[pval2<0.05])
length(sig_hyper)
length(sig_hypo)







