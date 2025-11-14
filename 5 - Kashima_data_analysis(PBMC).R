######################################
# Data Analysis of Kashima data(PBMC)
######################################

# include code for Kashima data analysis, supplementary figure 5
# include PBMC cell proportion data for figure 5a, 5b


# Library:
library(GEOquery)
library(dplyr)
library(ggplot2)
library(lumi)
library(EpiDISH)

# Load PBMC data from NCBI:
gset <- getGEO("GSE110828", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL13534", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

pd_PBMC_450K <- pData(gset)
dim(pd_PBMC_450K)   # 157  82
save(pd_PBMC_450K, file = "./pd_PBMC_450K.RData")
beta_PBMC_450K <- exprs(gset)
dim(beta_PBMC_450K)   # 410735    157
save(beta_PBMC_450K, file = "./beta_PBMC_450K.RData")

# Subset cord blood related samples only for the two datasets:
## PBMC
pd_PBMC <- pd_PBMC_450K[pd_PBMC_450K$source_name_ch1 == "CBMC", ]
dim(pd_PBMC) # 110  82
pd_PBMC$preeclampsia <- ifelse(pd_PBMC$preeclampsia == "Yes", "Case", "Control")
save(pd_PBMC, file = "./pd_PBMC.RData")
beta_PBMC <- beta_PBMC_450K[, colnames(beta_PBMC_450K) %in% rownames(pd_PBMC)]
dim(beta_PBMC) # 410735    110

save(beta_PBMC, file = "./beta_PBMC.RData")


## PBMC:regress only on sample group-----------------------------
library(limma)
library(lumi)
load("/nfs/dcmb-lgarmire/liuwent/15-Additional_Datasets/beta_PBMC.RData")
load("/nfs/dcmb-lgarmire/liuwent/15-Additional_Datasets/pd_PBMC.RData")
#exclude multi-hit probes
load("/nfs/dcmb-lgarmire/xtyang/CordBlood/multi_hit_new_450k.RData")
beta_PBMC = beta_PBMC[!rownames(beta_PBMC)%in%multi_hit_new_450k$probeID,]


# compare case and control group 
pd_PBMC$preeclampsia = relevel(factor(pd_PBMC$preeclampsia), ref="Control")
design0 = model.matrix(~pd_PBMC$preeclampsia)
fit0 = lmFit(beta2m(beta_PBMC), design0)
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


#adjust for confounders
library(EpiDISH)
#get cell proportion
load("/nfs/dcmb-lgarmire/xtyang/CordBlood/03-wenting_results/FlowSorted.CordBloodCombined.450k.compTable.rda")
ref = FlowSorted.CordBloodCombined.450k.compTable
out<-epidish(beta.m=beta_PBMC, ref.m=as.matrix(ref), method="CP")
count <- out$estF
pd_PBMC= data.frame(cbind(out$estF, pd_PBMC))
pd_PBMC$GA = as.numeric(sub(".*: ", "", pd_PBMC$characteristics_ch1.10))
pd_PBMC$BMI = as.numeric(sub(".*: ", "", pd_PBMC$characteristics_ch1.15))
pd_PBMC$parity = sub(".*: ", "", pd_PBMC$characteristics_ch1.12)
pd_PBMC$Age = as.numeric(sub(".*: ", "", pd_PBMC$characteristics_ch1.14))
pd_PBMC$batch = sub(".*: ", "", pd_PBMC$characteristics_ch1.1)
pd_PBMC$babysex = sub(".*: ", "", pd_PBMC$characteristics_ch1)
pd_PBMC$paternal_age = as.numeric(sub(".*: ", "", pd_PBMC$characteristics_ch1.16))
pd_PBMC$paternal_BMI = as.numeric(sub(".*: ", "", pd_PBMC$characteristics_ch1.17))


pd_PBMC$preeclampsia = relevel(factor(pd_PBMC$preeclampsia), ref="Control")
design1 = model.matrix(~preeclampsia+GA+BMI+parity+Age+batch+babysex+Gran+CD8T+CD4T+NK+Bcell+Mono, data = pd_PBMC)
fit1 = lmFit(beta2m(beta_PBMC), design1)
fit11 = eBayes(fit1)

## Get log_2 fold change and p value
log2fc = fit11$coefficients[,2]
##pval = p.adjust(fit$p.value[,2], method="bonferroni")
pval0 = p.adjust(fit11$p.value[,2], "BH")
length(pval0[pval0<0.05]) # 0

col_list = sapply(seq(1, dim(fit11$coefficients)[1]), function(i) {
  if (pval0[i]<0.05) {return("red")}
  else {return("black")} })
table(col_list)

png(file = "/nfs/dcmb-lgarmire/xtyang/CordBlood/03-wenting_results/PBMC_volcano_full_model.png")
par(cex.axis = 1.5, cex.lab = 1.5) 

volcanoplot(fit11, coef=2, col=col_list, highlight=5,
            names=rownames(fit1$coefficients))
# abline(h=-log10(0.05), lty=2, col="red")
abline(v=-1, lty=2, col="grey")
abline(v=1, lty=2, col="grey")
text(-1, 5, "x=-1", cex=1.2)
text(1, 5, "x=1", cex=1.2)
dev.off()
