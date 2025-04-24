# Analysis using data from Ching et al.

# 2. regress only on sample group-----------------------------
# Data:
load("/nfs/dcmb-lgarmire/liuwent/13-Travers_data/Trav_pd_all_minfi.RData")
load("/nfs/dcmb-lgarmire/liuwent/13-Travers_data/m.RData")#m values
load("/nfs/dcmb-lgarmire/liuwent/13-Travers_data/beta.RData")

design0 = model.matrix(~Trav_pd_all_minfi$Sample_Group)
design0[,2] = 1-design0[,2]
colnames(design0)[2] = "Cases"
fit0 = lmFit(m, design0)
fit00 = eBayes(fit0)

## Get log_2 fold change and p value
log2fc = fit00$coefficients[,2]
##pval = p.adjust(fit$p.value[,2], method="bonferroni")
pval0 = p.adjust(fit00$p.value[,2], "BH")
length(pval0[pval0<0.05]) # 328,614
cpg <- data.frame(logFC=log2fc, pval=pval0)
sig_hyper <- rownames(cpg[cpg$logFC>0 & cpg$pval<0.05,])
length(sig_hyper) # 211,371
sig_hypo <- rownames(cpg[cpg$logFC<0 & cpg$pval<0.05,])
length(sig_hypo) # 117,243
table_fit00 = topTable(fit00, num=Inf, coef=2, adjust.method = "BH")
table_fit00$Type = ifelse(table_fit00$logFC>0, "Hyper", "Hypo")
sum(table_fit00$adj.P.Val<0.05) # 328,614

## Make a volcano plot (label DMP with different colors)
col_list = sapply(seq(1, dim(fit00$coefficients)[1]), function(i) {
  if (pval0[i]<0.05) {return("red")}
  else {return("black")} })
table(col_list)

png(file = "/nfs/dcmb-lgarmire/xtyang/CordBlood/03-wenting_results/02-Travis_data/Trav_fit00_minfi.png")
par(cex.axis = 1.5, cex.lab = 1.5) 
volcanoplot(fit00, coef=2, col=col_list, highlight=5,
            names=rownames(fit00$coefficients), xlim = c(-4,4))

abline(v=-1, lty=2, col="grey")
abline(v=1, lty=2, col="grey")
text(-1, 5, "x=-1", cex=1.2)
text(1, 5, "x=1", cex=1.2)

dev.off()

# 3. Cell type deconvolution---------------------------------------
## 3.1 Regular version:
library(EpiDISH)
library(ggplot2)
library(reshape2)
library(dplyr)
library(limma)
library(ggpubr)
library(ggsignif)

load("/nfs/dcmb-lgarmire/xtyang/CordBlood/03-wenting_results/FlowSorted.CordBloodCombined.450k.compTable.rda")#Gervin ref
ref = as.matrix(FlowSorted.CordBloodCombined.450k.compTable)
out<-epidish(beta.m=beta, ref.m=ref, method="CP")
Trav_estF <- out$estF
Trav_estF[Trav_estF<0.00001]=0

Trav_pd_all = data.frame(cbind(Trav_estF, Trav_pd_all_minfi[,c(8:25)]))

dim(Trav_pd_all) # 19 18
save(Trav_pd_all, file = "Trav_pd_all.RData")

## 3.2 Relationship between phenotypes and cell types-----------------------------
# Library
library(dplyr)
library(limma)
library(ggplot2)
library(tidyverse)


pd = Trav_pd_all%>%dplyr::select("Sample_Group", "GAWeek", "MomAge", "BMI", "EthnicMom1")
pdnames = names(pd)
formstr <- paste0(pdnames, collapse = ' + ')
formstr <- paste0('~', formstr)
formstr <- as.formula(formstr)
design = model.matrix(formstr, data=pd)
fit_estF_cv = lmFit(t(Trav_estF), design)
fit_estF_cv = eBayes(fit_estF_cv)

hm_data <- data.frame(t(-log10(fit_estF_cv$p.value[, -1])))

hm_data <- hm_data %>% 
  rownames_to_column() %>%
  gather(colname, value, -rowname) %>% 
  setNames(c("Clinical_Variable", "Cell_Type", "Value"))
hm_data$Stars <- cut(hm_data$Value, breaks=c(-Inf, -log10(0.05), -log10(0.01), -log10(0.001), Inf), label=c(" ", "*", "**", "***"))

# Define the order of the factor levels for the x-axis
hm_data$Cell_Type <- factor(hm_data$Cell_Type, 
                            levels=c("Gran", "NK", "Bcell", "nRBC", "CD4T", "CD8T", "Mono"))

pdf("/home/liuwent/13-Travers_data/03a-Trav_hm_estF_cv.pdf")
ggplot(data = hm_data, aes(x=Cell_Type, y=Clinical_Variable, fill=Value)) +
  geom_tile()+
  scale_fill_gradient(low = "white", high = "#FF3300")+
  geom_text(aes(label = Stars), col = "#000000")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))+
  guides(fill=guide_legend(title="-log10(pval)"))+
  theme(text = element_text(size = 15))
dev.off()

## 3.3 Residual version:
# fit for all significant confounders except PE by using Trav_estF:
pd_noPE = Trav_pd_all%>%dplyr::select("GAWeek", "MomAge", "BMI", "EthnicMom1")
pdnames = names(pd_noPE)
formstr <- paste0(pdnames, collapse = ' + ')
formstr <- paste0('~', formstr)
formstr <- as.formula(formstr)
design = model.matrix(formstr, data=pd_noPE)
fit_noPE = lmFit(t(Trav_estF), design)
fit_noPE = eBayes(fit_noPE)

Residuals_newnew_markers <- residuals(fit_noPE, t(Trav_estF))
# Add mean:
Residuals_newnew_markers <- Residuals_newnew_markers+matrix(apply(t(Trav_estF), 1, mean),
                                                            nrow=ncol(Trav_estF),
                                                            ncol=nrow(Trav_estF))
count <- t(Residuals_newnew_markers)

df = data.frame(cbind(count, Trav_pd))

Sample_Group = df$Sample_Group
df1 = data.frame(cbind(df[, 1:7], Sample_Group))

df1_long = melt(df1, id.vars = "Sample_Group")

colnames(df1_long)[1] = "Sample_Group"
colnames(df1_long)[2] = "Cell_Type"
colnames(df1_long)[3] = "Cell_Type_Proportion_Residuals"

# Define the order of the factor levels for the x-axis
df1_long$Cell_Type <- factor(df1_long$Cell_Type, 
                             levels=c("Gran", "NK", "Bcell", "nRBC", "CD4T", "CD8T", "Mono"))

pdf("/home/liuwent/13-Travers_data/03b-Trav_Residuals_newnewSetofMarkers_CP(151,794).pdf")
ggplot(df1_long, aes(Cell_Type, Cell_Type_Proportion_Residuals, fill=Sample_Group)) + 
  geom_boxplot() + 
  # labs(title="Cell Type Residuals by Sample Group") + 
  # stat_compare_means(aes(label = after_stat(p.signif))) +
  annotate("text", label = c("***","ns","ns","ns","*","***","ns"), size = 4, col = "black",  x=c(1.1,2.1,3.1,4.1,5.1,6.1,7.1), y = 0.45, vjust=1, hjust=1)+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))+
  theme(text = element_text(size = 15))
dev.off()


# 4. Adjust for all confounders---------------------------------------------------------------
allConfounders = Trav_pd_all%>%dplyr::select("Sample_Group", "EthnicMom1", "MomAge", "BMI", "GAWeek", "BabySex",  
                                    "Mono", "Bcell", "CD8T", "nRBC", "Gran", "CD4T", )

ind <- sapply(allConfounders, is.numeric)
#scale numerical variables
f = function(x){scale(x, center = FALSE)}
allConfounders[ind] <- lapply(allConfounders[ind],f)
formula = paste0(names(allConfounders), collapse = ' + ')
formula = paste0("~", formula)
formula = formula(formula)

#design matrix
design2 = model.matrix(formula, data = allConfounders)
design2[,2] = 1-design2[,2]
colnames(design2)[2] = "Cases"
#fit linear model
fitConfounders = lmFit(m, design2)
fitConfounders = eBayes(fitConfounders)

#make a toptable on second column sample_group
allg.limma <- topTable(fitConfounders, coef="Cases", n=dim(fitConfounders)[1], adjust.method="fdr")
sigg.limma <- subset(allg.limma, adj.P.Val < 0.05)
nonsigg.limma <- subset(allg.limma, adj.P.Val >= 0.05)
dim(sigg.limma) # 50,605 # 51,486 (reported by Travers paper)

# Make a volcano plot (label DMP with different colors)
pval2 = p.adjust(fitConfounders$p.value[,2], method = "fdr")
col_list2 = sapply(seq(1, dim(fitConfounders$coefficients)[1]), function(i) {
  if (pval2[i]<0.05) {return("red")}
  else {return("black")} })
cpg <- data.frame(logFC=fitConfounders$coefficients[,2],pval=pval2)
sig_hyper <- rownames(cpg[cpg$logFC>0 & cpg$pval<0.05,])
length(sig_hyper) # 31,201 # 12,563 (reported by Travers paper)
sig_hypo <- rownames(cpg[cpg$logFC<0 & cpg$pval<0.05,])
length(sig_hypo) # 19,404

pdf("/nfs/dcmb-lgarmire/xtyang/CordBlood/03-wenting_results/02-Travis_data/Travis_volcano_large_new.pdf")
par(cex.axis = 1.5, cex.lab = 1.5)

volcanoplot(fitConfounders,coef=2,col=col_list2,highlight=1,
            names=rownames(fitConfounders$coefficients))
# abline(h=-log10(0.05), lty=2, col="red")
abline(v=-1, lty=2, col="grey")
abline(v=1, lty=2, col="grey")
text(-1, 5, "x=-1", cex=1.2)
text(1, 5, "x=1", cex=1.2)
dev.off()
