library(GEOquery)
library(limma)
library(EpiDISH)
library(ChAMP)
library(lumi)
library(bacon)

# load series and platform data from GEO---------------------------------------------------
# Fernado et al.
gset <- getGEO("GSE66459", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL13534", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
pheno <- pData(gset)
beta <- exprs(gset)
save(beta, file = "beta_Fernado.RData")
save(pheno, file = "pheno_Fernado.RData")

#reformat pd
pd_Fernado = pheno[, c("gestational_age (days):ch1","geo_accession", "description", "characteristics_ch1")]
colnames(pd_Fernado) = c("GA", "Sample_Name", "Preterm", "baby_sex")
pd_Fernado$GA = floor(as.integer(pd_Fernado$GA)/7)
pd_Fernado$Sample_Group = rep("Controls", nrow(pd_Fernado))
#pd_Fernado = cbind(count, pd_Fernado)
pd_Fernado$dataset = rep("Fernado", nrow(pd_Fernado))
pd_Fernado$Preterm = ifelse(pd_Fernado$Preterm=="Preterm sample", "preterm", "fullterm")
#save(pd_Fernado, file = "pd_Fernando.RData")

#load house data
load("/nfs/dcmb-lgarmire/liuwent/04-Full_Model/myCombat.RData")
load("/nfs/dcmb-lgarmire/xtyang/CordBlood/03-wenting_results/pd_all.RData")

#merge house data and Fernado et al.
common_probes = intersect(rownames(beta), rownames(myCombat))
merged_data = cbind(myCombat[common_probes, ], beta[common_probes, ])
#normalize
merged_data_norm = champ.norm(merged_data, method = "BMIQ")
save(merged_data_norm, file = "merged_data_norm.RData")

df$Preterm = ifelse(df$GA<37, "preterm", "fullterm")
df$dataset = rep("House", nrow(df))
df$baby_sex = rep("gender: Female", nrow(df))

merged_pd = rbind(df[, c("GA", "Sample_Name", "Preterm","baby_sex", "Sample_Group", "dataset")], pd_Fernado)
champ.QC(beta = merged_data_norm,
         pheno=merged_pd$dataset,
         mdsPlot=TRUE,
         densityPlot=FALSE,
         dendrogram=FALSE,
         PDFplot=TRUE,
         Rplot=TRUE,
         Feature.sel="None",
         resultsDir="./CHAMP_QCimages/")


#batch correction------------------------------------------------------
library(sva)

mod = model.matrix(~as.factor(Sample_Group)+as.factor(Preterm), data=merged_pd)
corrected_data = ComBat(dat=merged_data_norm, batch=merged_pd$dataset, mod=mod, par.prior=TRUE, prior.plots=FALSE)
champ.QC(beta = corrected_data,
         pheno=merged_pd$Sample_Group,
         mdsPlot=TRUE,
         densityPlot=FALSE,
         dendrogram=FALSE,
         PDFplot=TRUE,
         Rplot=TRUE,
         Feature.sel="None",
         resultsDir="./CHAMP_QCimages/")
corrected_data[corrected_data<0]=0.0001
corrected_data[corrected_data>1]=0.9999

#cell type deconvolution using corrected data----------------------------------------
load("/nfs/dcmb-lgarmire/xtyang/CordBlood/03-wenting_results/FlowSorted.CordBloodCombined.450k.compTable.rda")
ref = FlowSorted.CordBloodCombined.450k.compTable
out<-epidish(beta.m=corrected_data, ref.m=as.matrix(ref), method="CP")
count <- out$estF
merged_pd = data.frame(count, merged_pd)
save(merged_pd, file = "merged_pd.RData")

#full model -------------------------------------------------------------
corrected_data[corrected_data>1] = 0.9999
corrected_data[corrected_data<0] = 0.0001
Mvalues = beta2m(corrected_data)
save(Mvalues, file = "merged_batch_corrected_Mvalues.RData")

design = model.matrix(formula("~Sample_Group+GA+baby_sex+CD8T+CD4T+NK+Bcell+Mono+Gran+nRBC"), data = merged_pd)
fitConfounders = lmFit(Mvalues, design)
fitConfounders = eBayes(fitConfounders)
allg.limma <- limma::topTable(fitConfounders, coef="Sample_GroupDisease", n=dim(fitConfounders)[1], adjust.method="fdr")
sigg.limma <- subset(allg.limma, adj.P.Val < 0.05)

#bacon
t_stats <- fitConfounders$t[,"Sample_GroupDisease"] #random shuffled group
df<- fitConfounders$df.total[1]
z_score2 = zscoreT(t_stats, df, approx=FALSE, method = "bailey")
bc2 <- bacon(z_score2)
limma_inflation_score2 = inflation(bc2)
plot(bc2, type="qq")

p_corrected <- pval(bc2)
p_adj_corrected = p.adjust(p_corrected, method = "fdr")
sum(p_adj_corrected<0.05)

# Make a volcano plot (label DMP with different colors)
col_list2 = sapply(seq(1, dim(fitConfounders$coefficients)[1]), function(i) {
  if (p_adj_corrected[i]<0.05) {return("red")}
  else {return("black")} })
cpg <- data.frame(logFC=fitConfounders$coefficients[,2],pval=p_adj_corrected)
sig_hyper <- rownames(cpg[cpg$logFC>0 & cpg$pval<0.05,])
length(sig_hyper)#0
sig_hypo <- rownames(cpg[cpg$logFC<0 & cpg$pval<0.05,])
length(sig_hypo)

png("/nfs/dcmb-lgarmire/xtyang/CordBlood/10-preterm_data/add_Fernado.png")
par(cex.axis = 1.5, cex.lab = 1.5)

volcanoplot(fitConfounders,coef=2,col=col_list2,highlight=1,
            names=rownames(fitConfounders$coefficients), xlim=c(-4,4))
# abline(h=-log10(0.05), lty=2, col="red")
abline(v=-1, lty=2, col="grey")
abline(v=1, lty=2, col="grey")
text(-1, 5, "x=-1", cex=1.2)
text(1, 5, "x=1", cex=1.2)
dev.off()


#calculate cell proportion-----------------------------------------
library(tidyverse)
library(ggpubr)
new_pd_long = merged_pd%>%pivot_longer(1:7, names_to = "CellType", values_to = "Proportion")


p1 = ggplot(new_pd_long, aes(x = CellType, y = Proportion, fill = Sample_Group))+
  geom_boxplot()+stat_compare_means(label = "p.signif")
ggsave("cell_prop_merged.png",p1,width = 5,height = 5,dpi = 1200)

fit_CD8T = lm(CD8T~Sample_Group*GA, data = merged_pd)
fit_CD4T = lm(CD4T~Sample_Group*GA, data = merged_pd)
fit_Bcell = lm(Bcell~Sample_Group*GA, data = merged_pd)
fit_NK = lm(NK~Sample_Group*GA, data = merged_pd)
fit_Gran= lm(Gran~Sample_Group*GA, data = merged_pd)
fit_nRBC = lm(nRBC~Sample_Group*GA, data = merged_pd)
fit_Mono = lm(Mono~Sample_Group*GA, data = merged_pd)

p2 = ggplot(new_pd_long, aes(x = GA, y = Proportion, color = Sample_Group)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE) +
  facet_wrap(~ CellType, scales = "free_y", nrow = 1) +
  theme_minimal()+
  theme(strip.text.y.left = element_text(angle = 0),  # Horizontal labels
        strip.placement = "outside") +
  labs(x = "GA", y = NULL)
ggsave("cell_trend_merged.png",p2,width = 25,height = 5,dpi = 1200)
