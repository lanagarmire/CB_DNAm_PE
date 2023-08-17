# Library:
library(limma)
library(dbplyr)
library(lumi)
library(EpiDISH)

# Data:
load("/home/liuwent/04-Full_Model/myLoad.Rdata")
load("/home/liuwent/04-Full_Model/pd.RData")
load("/home/liuwent/04-Full_Model/myCombat.RData")
load("/home/liuwent/04d-DM_analysis_after_adj/annEPIC_gene.RData")
load("/home/liuwent/04d-DM_analysis_after_adj/order_new.Rdata")


# Get annotations on promoter regions:
annEPIC_gene <- as.data.frame(annEPIC_gene)
colnames(annEPIC_gene) <- c('CpG', 'Gene', 'Group')
annEPIC_gene_promoter = annEPIC_gene[annEPIC_gene$Group%in%c('TSS200','TSS1500'),]
annEPIC_gene_promoter = annEPIC_gene_promoter[complete.cases(annEPIC_gene_promoter),]
dim(annEPIC_gene_promoter)   # 172,173      3

# Match gene annotations and beta values:
myCombat_sub = myCombat[rownames(myCombat)[rownames(myCombat)%in%rownames(annEPIC_gene_promoter)],] 
dim(myCombat_sub)   #157,985     62

# Take the geometric mean of all cpgs that related to the same gene
cpg_geomean_beta <- NULL
for(gene_name in unique(annEPIC_gene_promoter$Gene)){
  cpg = annEPIC_gene_promoter[annEPIC_gene_promoter$Gene==gene_name,]$CpG
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
  cpg_geomean_beta = rbind(cpg_geomean_beta, cpg_df)
}
cpg_geomean_beta[cpg_geomean_beta < 0.001]=0.001
cpg_geomean_beta[cpg_geomean_beta > 0.999]=0.999
cpg_geomean_mval = beta2m(cpg_geomean_beta) 
dim(cpg_geomean_mval)   #24353    62
# save(cpg_geomean_mval, file='cpg_geomean_mval.RData')
# load('/home/liuwent/04d-DM_analysis_after_adj/cpg_geomean_mval.RData')

## model1: regress only on sample group-----------------------------------------
myLoad$pd$Sample_Group = relevel(factor(myLoad$pd$Sample_Group), ref="Controls")
design0 = model.matrix(~myLoad$pd$Sample_Group)
fit0 = lmFit(cpg_geomean_mval, design0)
fit_genemean00 = eBayes(fit0)
save(fit_genemean00, file = "fit_genemean00.RData")
load('/home/liuwent/04d-DM_analysis_after_adj/fit_genemean00.RData')

table_fit_genemean00 = topTable(fit_genemean00, num=Inf, coef=2, adjust.method = "BH")
table_fit_genemean00$Type = ifelse(table_fit_genemean00$logFC>0, "Hyper", "Hypo")
sum(table_fit_genemean00$adj.P.Val<0.05) 

## Get log_2 fold change and p value
log2fc = fit_genemean00$coefficients[,2]
##pval = p.adjust(fit$p.value[,2], method="bonferroni")
pval0 = p.adjust(fit_genemean00$p.value[,2], "BH")
length(pval0[pval0<0.05])   # 4,767

## Make a volcano plot (label DMP with different colors)
col_list = sapply(seq(1, dim(fit_genemean00$coefficients)[1]), function(i) {
  if (pval0[i]<0.05) {return("red")}
  else {return("black")} })

pdf(file = "/home/liuwent/04d-DM_analysis_after_adj/fit_genemean00.pdf", width = 10, height = 10)
volcanoplot(fit_genemean00, coef=2, col=col_list, highlight=5,
            names=rownames(fit_genemean00$coefficients))
# abline(h=-log10(0.05), lty=2, col="red")
abline(v=-1, lty=2, col="grey")
abline(v=1, lty=2, col="grey")
text(-1, 5, "x=-1")
text(1, 5, "x=1")
dev.off()

### model2: six significant confounders according to SOV plot-------------------
allConfounders = pd[, c("Sample_Group", "GA", "BMI", "Eth2", "Parity", "Age")]

ind <- sapply(allConfounders, is.numeric)
#scale numerical variables
f = function(x){scale(x, center = FALSE)}
allConfounders[ind] <- lapply(allConfounders[ind],f)
formula = paste0(names(allConfounders), collapse = ' + ')
formula = paste0("~", formula)
formula = formula(formula)
#design matrix
design2 = model.matrix(formula, data = allConfounders)
#???#default setting control = 1, need to set case = 1???#
# design2[,2] = 1-design2[,2]
colnames(design2)[2] = "Cases"
#fit linear model
fitConfounders = lmFit(cpg_geomean_mval, design2)
fitConfounders_genemean = eBayes(fitConfounders)
save(fitConfounders_genemean, file = "fitConfounders_genemean.RData")
load('/home/liuwent/04d-DM_analysis_after_adj/fitConfounders_genemean.RData')

table_fitConfounders_genemean = topTable(fitConfounders_genemean, num=Inf, coef=2, adjust.method = "BH")
table_fitConfounders_genemean$Type = ifelse(table_fitConfounders_genemean$logFC>0, "Hyper", "Hypo")
sum(table_fitConfounders_genemean$adj.P.Val<0.05)

# Make a volcano plot (label DMP with different colors)
pval2 = p.adjust(fitConfounders_genemean$p.value[,2], "BH")
col_list2 = sapply(seq(1, dim(fitConfounders_genemean$coefficients)[1]), function(i) {
  if (pval2[i]<0.05) {return("red")}
  else {return("black")} })
cpg <- data.frame(logFC=fitConfounders_genemean$coefficients[,2],pval=pval2)
sig_hyper <- rownames(cpg[cpg$logFC>0 & cpg$pval<0.05,])
# save(sig_hyper, file='sig_hyper.RData')
sig_hypo <- rownames(cpg[cpg$logFC<0 & cpg$pval<0.05,])
# save(sig_hypo, file='sig_hypo.RData')

pdf("/home/liuwent/04d-DM_analysis_after_adj/fitConfounders_genemean.pdf", width = 10, height = 10)
volcanoplot(fitConfounders_genemean,coef=2,col=col_list2,highlight=5,names=rownames(fitConfounders_genemean$coefficients))
# abline(h=-log10(0.05), lty=2, col="red")
abline(v=-1, lty=2, col="grey")
abline(v=1, lty=2, col="grey")
text(-1, 5, "x=-1")
text(1, 5, "x=1")
dev.off()

length(pval2[pval2<0.05]) # 2
length(sig_hyper) # 1
length(sig_hypo) # 1

### model3: all confounders-----------------------------------------------------
pd = order_new
allConfounders = pd%>%dplyr::select("Sample_Group", "GA", "BMI", "Eth2", "Parity", "Age", 
                                    "CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran", "nRBC")

ind <- sapply(allConfounders, is.numeric)
#scale numerical variables
f = function(x){scale(x, center = FALSE)}
allConfounders[ind] <- lapply(allConfounders[ind],f)
formula = paste0(names(allConfounders), collapse = ' + ')
formula = paste0("~", formula)
formula = formula(formula)
#design matrix
design2 = model.matrix(formula, data = allConfounders)
colnames(design2)[2] = "Cases"
#fit linear model
fitConfounders = lmFit(cpg_geomean_mval, design2)
fitAllConfounders_genemean = eBayes(fitConfounders)
save(fitAllConfounders_genemean, file = "fitAllConfounders_genemean.RData")
load('/home/liuwent/04d-DM_analysis_after_adj/fitAllConfounders_genemean.RData')

table_fitAllConfounders_genemean = topTable(fitAllConfounders_genemean, num=Inf, coef=2, adjust.method = "BH")
table_fitAllConfounders_genemean$Type = ifelse(table_fitAllConfounders_genemean$logFC>0, "Hyper", "Hypo")
sum(table_fitAllConfounders_genemean$adj.P.Val<0.05)

# Make a volcano plot (label DMP with different colors)
pval2 = p.adjust(fitAllConfounders_genemean$p.value[,2], method = "fdr")
col_list2 = sapply(seq(1, dim(fitAllConfounders_genemean$coefficients)[1]), function(i) {
  if (pval2[i]<0.05) {return("red")}
  else {return("black")} })
cpg <- data.frame(logFC=fitAllConfounders_genemean$coefficients[,2],pval=pval2)
sig_hyper <- rownames(cpg[cpg$logFC>0 & cpg$pval<0.05,])
# save(sig_hyper, file='sig_hyper.RData')
sig_hypo <- rownames(cpg[cpg$logFC<0 & cpg$pval<0.05,])
# save(sig_hypo, file='sig_hypo.RData')

pdf("/home/liuwent/04d-DM_analysis_after_adj/fitAllConfounders_genemean.pdf", width = 10, height = 10)
volcanoplot(fitAllConfounders_genemean,coef=2,col=col_list2,highlight=5,names=rownames(fitAllConfounders_genemean$coefficients))
# abline(h=-log10(0.05), lty=2, col="red")
abline(v=-1, lty=2, col="grey")
abline(v=1, lty=2, col="grey")
text(-1, 5, "x=-1")
text(1, 5, "x=1")
dev.off()

length(pval2[pval2<0.05]) # 0
length(sig_hyper) # 0
length(sig_hypo) # 0