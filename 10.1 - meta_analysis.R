library(tidyr)
library(dplyr)
library(ChAMP)
library(limma)
library(lumi)
library(EpiDISH)
library(ggplot2)

setwd("/nfs/dcmb-lgarmire/xtyang/CordBlood/09-meta_analysis/")
#load data from Travis_data_analysis.R
load("/nfs/dcmb-lgarmire/liuwent/13-Travers_data/myLoad.RData")
travis_data_load = myLoad
travis_data_pd = myLoad$pd
rm(myLoad)
#load data from 1-Raw data preprocess.R
load("/nfs/dcmb-lgarmire/liuwent/04-Full_Model/myLoad.Rdata")
epic_load = myLoad
epic_pd = myLoad$pd
rm(myLoad)

#beta_eolo has only blood data beta_eolo_450k has both blood and tissue
#EOLO data was normalized average beta, epic and travis data are not normalized
# load Herzog et al. from 7.1-Herzog-data-analysis.R
beta = merge(epic_load$beta, travis_data_load$beta, by = "row.names", all = FALSE)
rownames(beta) = beta[,1]
beta = as.matrix(beta[,2:ncol(beta)])

champ.QC(beta = beta, pheno = c(rep("travis", 19), rep("epic", 62)), 
         dendrogram = FALSE, 
         resultsDir = "/nfs/dcmb-lgarmire/xtyang/CordBlood/09-meta_analysis/")



#Merge pd data------------------------------------------------------
load("/nfs/dcmb-lgarmire/liuwent/15-Additional_Datasets/pd_EOLO.RData")

pd_EOLO <- pd_EOLO %>%
  separate(title, into = c("Slide", "Array"), sep = "_") %>%
  select(Slide, Array, group, `group:ch1`) %>%
  mutate(cohort = "EOLO", 
         Sample_Name = rownames(.), 
         Sample_Group = group,
         detailed_group = `group:ch1`) %>%
  select(-`group:ch1`, -group)
pd_EOLO$Slide <- as.numeric(pd_EOLO$Slide)

epic_pd <- epic_pd %>%
  mutate(cohort = "epic", 
         Sample_Name = as.character(Sample_Name)) %>%
  select(Sample_Name, Slide, Array, Sample_Group, BMI, Eth2, Age, cohort)

travis_data_pd <- travis_data_pd %>%
  mutate(cohort = "travis",
         Eth2 = EthnicMom1,
         Age = MomAge,
         Sample_Name = as.character(Sample_Name)) %>%
  select(Sample_Name, Slide, Array, Sample_Group, BMI, Eth2, Age, cohort)

#save(pd_EOLO, file = "/nfs/dcmb-lgarmire/xtyang/CordBlood/09-meta_analysis/pd_EOLO_brief.RData")

merged_pd <- bind_rows(epic_pd, travis_data_pd, pd_EOLO)
for(i in 1:nrow(merged_pd)){
  if(merged_pd$Sample_Group[i] == "CASE"){
    merged_pd$Sample_Group[i] = "Case"
  }else if(merged_pd$Sample_Group[i] == "CONTROL"){
    merged_pd$Sample_Group[i] = "Control"
  }else if(merged_pd$Sample_Group[i] == "Controls"){
    merged_pd$Sample_Group[i] = "Control"
  }else if(merged_pd$Sample_Group[i] == "Disease"){
    merged_pd$Sample_Group[i] = "Case"
  }
}
save(merged_pd, file = "/nfs/dcmb-lgarmire/xtyang/CordBlood/09-meta_analysis/merged_pd.RData")


#normalize and plotting--------------------------------------------------------------------------------
load("/nfs/dcmb-lgarmire/liuwent/15-Additional_Datasets/beta_EOLO.RData")
load("/nfs/dcmb-lgarmire/xtyang/CordBlood/09-meta_analysis/merged_pd.RData")

#normalize epic and travis data together
my_norm = champ.norm(beta = beta, arraytype = "450k", cores = 8)
merged_norm = merge(my_norm, beta_EOLO,by = "row.names", all = FALSE)
rownames(merged_norm) = merged_norm$Row.names
merged_norm = merged_norm[,2:ncol(merged_norm)]
save(merged_norm, file = "/nfs/dcmb-lgarmire/xtyang/CordBlood/09-meta_analysis/merged_norm.RData")

#check for batch effects
champ.QC(beta = as.matrix(merged_norm),pheno = merged_pd$cohort, 
dendrogram = FALSE, 
resultsDir = "/nfs/dcmb-lgarmire/xtyang/CordBlood/09-meta_analysis/")

champ.SVD(beta = merged_norm, 
          pd = merged_pd,
          resultsDir = "/nfs/dcmb-lgarmire/xtyang/CordBlood/09-meta_analysis/")


#Harmonizing across datasets-------------------------------------------
load("/nfs/dcmb-lgarmire/xtyang/CordBlood/09-meta_analysis/merged_norm.RData")
load("/nfs/dcmb-lgarmire/xtyang/CordBlood/09-meta_analysis/merged_pd.RData")
harmonized_norm = champ.runCombat(beta = merged_norm, pd = merged_pd, 
                                  variablename = "Sample_Group",
                                  batchname = c("cohort", "Array"))#slide is correlated with cohort
champ.SVD(beta = as.data.frame(harmonized_norm), 
          pd = merged_pd,
          resultsDir = "/nfs/dcmb-lgarmire/xtyang/CordBlood/09-meta_analysis/")

champ.QC(beta = as.matrix(harmonized_norm),pheno = merged_pd$cohort, 
         dendrogram = FALSE, 
         resultsDir = "/nfs/dcmb-lgarmire/xtyang/CordBlood/09-meta_analysis/")

save(harmonized_norm, file = "/nfs/dcmb-lgarmire/xtyang/CordBlood/09-meta_analysis/harmonized_norm.RData")

# identify significant cpgs with null model-------------------------------------------
design0 = model.matrix(~merged_pd$Sample_Group)
design0[,2] = 1-design0[,2]
colnames(design0)[2] = "Cases"
fit0 = lmFit(beta2m(harmonized_norm), design0)#convert beta to m values
fit00 = eBayes(fit0)

pval0 = p.adjust(fit00$p.value[,2], "BH")
table_fit00 = topTable(fit00, num=Inf, coef=2, adjust.method = "BH")
table_fit00$Type = ifelse(table_fit00$logFC>0, "Hyper", "Hypo")
sum(table_fit00$adj.P.Val<0.05)#113375

col_list = sapply(seq(1, dim(fit00$coefficients)[1]), function(i) {
  if (pval0[i]<0.05) {return("red")}
  else {return("black")} })
table(col_list)

png(file = "/nfs/dcmb-lgarmire/xtyang/CordBlood/09-meta_analysis/fit00.png")
par(cex.axis = 1.5, cex.lab = 1.5) 
volcanoplot(fit00, coef=2, col=col_list, highlight=5,
            names=rownames(fit00$coefficients))
# abline(h=-log10(0.05), lty=2, col="red")
abline(v=-1, lty=2, col="grey")
abline(v=1, lty=2, col="grey")
text(-1, 5, "x=-1", cex=1.2)
text(1, 5, "x=1", cex=1.2)
dev.off()


# cell type deconvolution-----------------------------------------------------------

#load Gervin reference
load("/nfs/dcmb-lgarmire/xtyang/CordBlood/03-wenting_results/FlowSorted.CordBloodCombined.450k.compTable.rda")
ref = FlowSorted.CordBloodCombined.450k.compTable
out<-epidish(beta.m=harmonized_norm, ref.m=as.matrix(ref), method="CP")
count <- out$estF

pd = data.frame(cbind(merged_pd,out$estF))
save(pd, file = "/nfs/dcmb-lgarmire/xtyang/CordBlood/09-meta_analysis/pd_all.RData")

count_long = pd%>%pivot_longer(10:16, names_to = "cell_type", values_to = "proportion")
count_long$Sample_Group = ifelse(count_long$Sample_Group == "Case", "Disease", "Countrol")
count_long$cell_type = factor(count_long$cell_type, levels = c("Gran", "NK", "Bcell", "nRBC", "CD4T","CD8T","Mono"))

cell_type_pval <- c(
  t.test(pd$Gran ~ pd$Sample_Group)$p.value,
  t.test(pd$NK ~ pd$Sample_Group)$p.value,
  t.test(pd$Bcell ~ pd$Sample_Group)$p.value,
  t.test(pd$nRBC ~ pd$Sample_Group)$p.value,
  t.test(pd$CD4T ~ pd$Sample_Group)$p.value,
  t.test(pd$CD8T ~ pd$Sample_Group)$p.value,
  t.test(pd$Mono ~ pd$Sample_Group)$p.value
)
pval_star = c("***", "ns", "***", "***", "***", "***", "**")
pdf("/nfs/dcmb-lgarmire/xtyang/CordBlood/09-meta_analysis/cell_prop.pdf")
ggplot(count_long, aes(x = cell_type, y = proportion, fill = Sample_Group))+
  geom_boxplot()+
  annotate("text", label = pval_star, size = 4, col = "black",  x=c(1.1,2.1,3.1,4.1,5.1,6.1,7.1), y = 0.75, vjust=1, hjust=1)+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))+
  theme(text = element_text(size = 15))
dev.off()


# identify significant cpgs with cell type proportion model-------------------------------------------
allConfounders = pd%>%dplyr::select("Sample_Group", "CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran", "nRBC")
formula = paste0(names(allConfounders), collapse = ' + ')
formula = paste0("~", formula)
formula = formula(formula)
#design matrix
design2 = model.matrix(formula, data = allConfounders)
design2[,2] = 1-design2[,2]
colnames(design2)[2] = "Cases"
fit2 = lmFit(beta2m(harmonized_norm), design2)#convert beta to m values
fit22 = eBayes(fit2)

pval2 = p.adjust(fit22$p.value[,2], "BH")
table_fit22 = topTable(fit22, num=Inf, coef=2, adjust.method = "BH")
table_fit22$Type = ifelse(table_fit22$logFC>0, "Hyper", "Hypo")
sum(table_fit22$adj.P.Val<0.05)

col_list = sapply(seq(1, dim(fit22$coefficients)[1]), function(i) {
  if (pval2[i]<0.05) {return("red")}
  else {return("black")} })
table(col_list)

png(file = "/nfs/dcmb-lgarmire/xtyang/CordBlood/09-meta_analysis/fit22.png")
par(cex.axis = 1.5, cex.lab = 1.5) 
volcanoplot(fit22, coef=2, col=col_list, highlight=5,
            names=rownames(fit22$coefficients))
# abline(h=-log10(0.05), lty=2, col="red")
abline(v=-1, lty=2, col="grey")
abline(v=1, lty=2, col="grey")
text(-1, 5, "x=-1", cex=1.2)
text(1, 5, "x=1", cex=1.2)
dev.off()
