# Library:
library(lumi)
library(limma)
library(dplyr)
library(ggplot2)
library(ggrepel)

#load from 2.0 - cell_type_deconvolution_Gervin.R
load("/nfs/dcmb-lgarmire/liuwent/04-Full_Model/pd_all.RData")
# Shortcut from 1 - Raw data pre-processing.Rmd
load("/nfs/dcmb-lgarmire/liuwent/04-Full_Model/Mvalues.RData")

#load clinical covariate and estimated cell proportions
allConfounders = pd%>%dplyr::select("Sample_Group", "GA", "Age", "Parity", "Eth2", 
                                    "CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran", "nRBC", "Smoker")

formula = paste0(names(allConfounders), collapse = ' + ')
formula = paste0("~", formula)
formula = formula(formula)
#design matrix
design2 = model.matrix(formula, data = allConfounders)
#???#default setting control = 1, need to set case = 1???#
colnames(design2)[2] = "Cases"
#fit linear model
fitConfounders = lmFit(Mvalues, design2)
fitConfounders = eBayes(fitConfounders)
save(fitConfounders, file = "fitConfounders_Garvin.RData")

#make a toptable on second column sample_group
allg.limma <- topTable(fitConfounders, coef="Cases", n=dim(fitConfounders)[1], adjust.method="fdr")
sigg.limma <- subset(allg.limma, adj.P.Val < 0.05)
save(sigg.limma, file='sigg.limma.RData')
nonsigg.limma <- subset(allg.limma, adj.P.Val >= 0.05)
dim(sigg.limma)[1]#0

# Make a volcano plot (label DMP with different colors)
col_of_interest = "Cases"
pval2 = p.adjust(fitConfounders$p.value[,col_of_interest], method = "fdr")
col_list2 = sapply(seq(1, dim(fitConfounders$coefficients)[1]), function(i) {
  if (pval2[i]<0.05) {return("red")}
  else {return("black")} })
cpg <- data.frame(logFC=fitConfounders$coefficients[,2],pval=pval2)
sig_hyper <- rownames(cpg[cpg$logFC>0 & cpg$pval<0.05,])
length(sig_hyper)#0
# save(sig_hyper, file='sig_hyper.RData')
sig_hypo <- rownames(cpg[cpg$logFC<0 & cpg$pval<0.05,])
# save(sig_hypo, file='sig_hypo.RData')
length(sig_hypo)

png("/nfs/dcmb-lgarmire/xtyang/CordBlood/03-wenting_results/01-new_result_Garvin_ref/volcano_Housedata_Gervin.png")
par(cex.axis = 1.5, cex.lab = 1.5)

volcanoplot(fitConfounders,coef=2,col=col_list2,highlight=1,
            names=rownames(fitConfounders$coefficients), xlim=c(-4,4))
# abline(h=-log10(0.05), lty=2, col="red")
abline(v=-1, lty=2, col="grey")
abline(v=1, lty=2, col="grey")
text(-1, 5, "x=-1", cex=1.2)
text(1, 5, "x=1", cex=1.2)
dev.off()
