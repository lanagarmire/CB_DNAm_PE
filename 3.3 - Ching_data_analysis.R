##########################################
# Data Analysis of Ching et al. Dataset
##########################################
# include code for figure 3e, 3f


library(lumi)
library(limma)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(EpiDISH)
library(reshape2)
library(ggpubr)

# 1. cell type deconvolution with Gervin et al.'s reference
load('/nfs/dcmb-lgarmire/xtyang/CordBlood/03-wenting_results/02-Travis_data/myNorm2.RData')
load('/nfs/dcmb-lgarmire/xtyang/CordBlood/03-wenting_results/02-Travis_data/pd2.RData')
load("/nfs/dcmb-lgarmire/xtyang/CordBlood/03-wenting_results/FlowSorted.CordBloodCombined.450k.compTable.rda")
ref = FlowSorted.CordBloodCombined.450k.compTable
out<-epidish(beta.m=myNorm, ref.m=as.matrix(ref), method="CP")
count <- out$estF
pd = data.frame(count, pd)

# 2. regress only on sample group-----------------------------
# Data:
m = beta2m(myNorm)#convert to m value

design0 = model.matrix(~pd$Sample_Group)
colnames(design0)[2] = "Cases"
design0[,2] = 1- design0[,2]#set control as 0, case as 1
fit0 = lmFit(m, design0)
fit00 = eBayes(fit0)

## Get log_2 fold change and p value
log2fc = fit00$coefficients[,2]
pval0 = p.adjust(fit00$p.value[,2], "BH")
allg.limma <- limma::topTable(fit00, coef="Cases", n=dim(fit00)[1], adjust.method="BH")
sigg.limma <- subset(allg.limma, adj.P.Val < 0.05)
sig_cpg = rownames(sigg.limma)
write.csv(sig_cpg, file = "/nfs/dcmb-lgarmire/xtyang/CordBlood/03-wenting_results/02-Travis_data/sig_cpg_before_adjust.csv")

## Make a volcano plot (label DMP with different colors)
col_list = sapply(seq(1, dim(fit00$coefficients)[1]), function(i) {
  if (pval0[i]<0.05) {return("red")}
  else {return("black")} })

png(file = "/nfs/dcmb-lgarmire/xtyang/CordBlood/03-wenting_results/02-Travis_data/Trav_figure3E.png")
par(cex.axis = 1.5, cex.lab = 1.5) 
volcanoplot(fit00, coef=2, col=col_list, highlight=0,
            names=rownames(fit00$coefficients))
# abline(h=-log10(0.05), lty=2, col="red")
abline(v=-1, lty=2, col="grey")
abline(v=1, lty=2, col="grey")
text(-1, 5, "x=-1", cex=1.2)
text(1, 5, "x=1", cex=1.2)
dev.off()


# 3. Adjust for all confounders---------------------------------------------------------------
library(compositions)
library(zCompositions)

cell_type = pd[,1:7]
cell_type[cell_type < 0] <- 0          # guard against negatives
cell_type_zr <- cmultRepl(cell_type, method = "CZM", output = "prop")  # compositional zero replacement
cell_ilr <- ilr(acomp(cell_type))
colnames(cell_ilr) = c("ct1", "ct2", "ct3", "ct4", "ct5","ct6")

#m = m[, rownames(pd)]
allConfounders = pd%>%dplyr::select("Sample_Group", "Eth2", "MomAge", "BMI", "GAWeek", "BabySex")%>%cbind(cell_ilr)

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

#calculate bacon inflation score
t_stats <- fitConfounders$t[,"Cases"] #random shuffled group
df<- fitConfounders$df.total[1]
z_score2 = zscoreT(t_stats, df, approx=FALSE, method = "bailey")
bc2 <- bacon(z_score2)
limma_inflation_score2 = inflation(bc2)#1.125

#make a toptable on second column sample_group
allg.limma <- limma::topTable(fitConfounders, coef="Cases", n=dim(fitConfounders)[1], adjust.method="BH")
sigg.limma <- subset(allg.limma, adj.P.Val < 0.05)


# Make a volcano plot (label DMP with different colors)
col_list2 = sapply(seq(1, dim(fitConfounders$coefficients)[1]), function(i) {
  if (allg.limma$adj.P.Val[i]<0.05) {return("red")}
  else {return("black")} })
cpg <- data.frame(logFC=fitConfounders$coefficients[,2],pval=allg.limma$adj.P.Val)
sig_hyper <- rownames(cpg[cpg$logFC>0 & cpg$pval<0.05,])
length(sig_hyper)#0
sig_hypo <- rownames(cpg[cpg$logFC<0 & cpg$pval<0.05,])
length(sig_hypo)

png("/nfs/dcmb-lgarmire/xtyang/CordBlood/03-wenting_results/02-Travis_data/Trav_volcano3F.png")
par(cex.axis = 1.5, cex.lab = 1.5)

volcanoplot(fitConfounders,coef=2,col=col_list2,highlight=0,xlim = c(-15,15),
            names=rownames(fitConfounders$coefficients))
# abline(h=-log10(0.05), lty=2, col="red")
abline(v=-1, lty=2, col="grey")
abline(v=1, lty=2, col="grey")
text(-1, 5, "x=-1", cex=1.2)
text(1, 5, "x=1", cex=1.2)
dev.off()

