#####################################
# Pooled analysis of three datasets
#####################################

# include code for figure 3J, 3K and 3L

library(tidyr)
library(dplyr)
library(ChAMP)
library(limma)
library(lumi)
library(EpiDISH)
library(ggplot2)

setwd("/nfs/dcmb-lgarmire/xtyang/CordBlood/09-meta_analysis/")
load('/nfs/dcmb-lgarmire/xtyang/CordBlood/03-wenting_results/02-Travis_data/myLoad2.RData')
travis_data_load = myLoad
travis_data_pd = myLoad$pd
rm(myLoad)

load('/nfs/dcmb-lgarmire/xtyang/CordBlood/03-wenting_results/01-new_result_Garvin_ref/myLoad_new.RData')
epic_load = myLoad
epic_pd = myLoad$pd
rm(myLoad)

#beta_eolo has only blood data beta_eolo_450k has both blood and tissue
#EOLO data was normalized average beta, epic and travis data are not normalized
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
  dplyr::select(Slide, Array, group, `group:ch1`) %>%
  mutate(cohort = "EOLO", 
         Sample_Name = rownames(.), 
         Sample_Group = group,
         detailed_group = `group:ch1`) %>%
  dplyr::select(-`group:ch1`, -group)
pd_EOLO$Slide <- as.character(pd_EOLO$Slide)

epic_pd <- epic_pd %>%
  mutate(cohort = "epic", 
         Sample_Name = as.character(Sample_Name)) %>%
  dplyr::select(Sample_Name, Slide, Array, Sample_Group, BMI, Eth2, Age, cohort)
epic_pd$Slide = as.character(epic_pd$Slide)

travis_data_pd <- travis_data_pd %>%
  mutate(cohort = "travis",
         Eth2 = EthnicMom1,
         Age = MomAge,
         Sample_Name = as.character(Sample_Name)) %>%
  dplyr::select(Sample_Name, Slide, Array, Sample_Group, BMI, Eth2, Age, cohort)
travis_data_pd$Slide = as.character(travis_data_pd$Slide)

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
save(merged_pd, file = "/nfs/dcmb-lgarmire/xtyang/CordBlood/09-meta_analysis/merged_pd_new.RData")


#normalize and plotting--------------------------------------------------------------------------------
load("/nfs/dcmb-lgarmire/liuwent/15-Additional_Datasets/beta_EOLO.RData")
load("/nfs/dcmb-lgarmire/xtyang/CordBlood/09-meta_analysis/merged_pd.RData")
load("/nfs/dcmb-lgarmire/xtyang/CordBlood/multi_hit_new_450k.RData")
beta_EOLO = beta_EOLO[!rownames(beta_EOLO)%in%multi_hit_new_450k$probeID,]

#normalize epic and travis data together
my_norm = champ.norm(beta = beta, arraytype = "450k", cores = 8)
merged_norm = merge(my_norm, beta_EOLO,by = "row.names", all = FALSE)
rownames(merged_norm) = merged_norm$Row.names
merged_norm = merged_norm[,2:ncol(merged_norm)]
save(merged_norm, file = "/nfs/dcmb-lgarmire/xtyang/CordBlood/09-meta_analysis/merged_norm_new.RData")

#check for batch effects
champ.QC(beta = as.matrix(merged_norm),pheno = merged_pd$cohort, 
dendrogram = FALSE, 
resultsDir = "/nfs/dcmb-lgarmire/xtyang/CordBlood/09-meta_analysis/CHAMP_QC")

champ.SVD(beta = merged_norm, 
          pd = merged_pd,
          resultsDir = "/nfs/dcmb-lgarmire/xtyang/CordBlood/09-meta_analysis/CHAMP_SVD")


#Harmonizing across datasets-------------------------------------------
load("/nfs/dcmb-lgarmire/xtyang/CordBlood/09-meta_analysis/merged_norm_new.RData")
load("/nfs/dcmb-lgarmire/xtyang/CordBlood/09-meta_analysis/merged_pd.RData")
merged_pd = merged_pd[merged_pd$Sample_Name!="8", ]#exclude one low quality sample from Ching et al.data
harmonized_norm = champ.runCombat(beta = merged_norm, pd = merged_pd[, c("Sample_Group","cohort", "Slide", "Array")], 
                                  variablename = "Sample_Group",
                                  batchname = c("cohort"))#slide is correlated with cohort
champ.SVD(beta = as.data.frame(harmonized_norm), 
          pd = merged_pd,
          resultsDir = "/nfs/dcmb-lgarmire/xtyang/CordBlood/09-meta_analysis/")

champ.QC(beta = as.matrix(harmonized_norm),pheno = merged_pd$cohort, 
         dendrogram = FALSE, 
         resultsDir = "/nfs/dcmb-lgarmire/xtyang/CordBlood/09-meta_analysis/CHAMP_QC/")

save(harmonized_norm, file = "/nfs/dcmb-lgarmire/xtyang/CordBlood/09-meta_analysis/harmonized_norm_new.RData")

# identify significant cpgs with null model-------------------------------------------
load("/nfs/dcmb-lgarmire/xtyang/CordBlood/09-meta_analysis/harmonized_norm_new.RData")
design0 = model.matrix(~merged_pd$Sample_Group)
design0[,2] = 1-design0[,2]
colnames(design0)[2] = "Cases"
fit0 = lmFit(beta2m(harmonized_norm), design0)#convert beta to m values
fit00 = eBayes(fit0)

pval0 = p.adjust(fit00$p.value[,2], "BH")
table_fit00 = limma::topTable(fit00, num=Inf, coef=2, adjust.method = "BH")
table_fit00$Type = ifelse(table_fit00$logFC>0, "Hyper", "Hypo")
sum(table_fit00$adj.P.Val<0.05)#92386

col_list = sapply(seq(1, dim(fit00$coefficients)[1]), function(i) {
  if (pval0[i]<0.05) {return("red")}
  else {return("black")} })
table(col_list)

png(file = "/nfs/dcmb-lgarmire/xtyang/CordBlood/09-meta_analysis/volcano3K.png")
par(cex.axis = 1.5, cex.lab = 1.5) 
volcanoplot(fit00, coef=2, col=col_list, highlight=0,
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
save(pd, file = "/nfs/dcmb-lgarmire/xtyang/CordBlood/09-meta_analysis/pd_all_new.RData")

count_long = pd%>%pivot_longer(10:16, names_to = "cell_type", values_to = "proportion")
count_long$Sample_Group = ifelse(count_long$Sample_Group == "Case", "Disease", "Control")
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
pval_star = c("***", "ns", "***", "***", "*", "*", "***")
pdf("/nfs/dcmb-lgarmire/xtyang/CordBlood/09-meta_analysis/cell_prop_1005.pdf")
ggplot(count_long, aes(x = cell_type, y = proportion, fill = Sample_Group))+
  geom_boxplot()+
  annotate("text", label = pval_star, size = 4, col = "black",  
           x=c(1.1,2.1,3.1,4.1,5.1,6.1,7.1), y = 0.75, vjust=1, hjust=1)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5),
        axis.title.x=element_blank())+
  theme(text = element_text(size = 15))+
  ylab("Cell Proportion (unadjusted)")+
  labs(fill='Sample Group') 
dev.off()


# identify significant cpgs with cell type proportion model-------------------------------------------
library(compositions)
library(zCompositions)
cell_type = pd%>%dplyr::select("CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran", "nRBC")
cell_type[cell_type < 0] <- 0          # guard against negatives
cell_type_zr <- cmultRepl(cell_type, method = "CZM", output = "prop")  # compositional zero replacement
cell_ilr <- ilr(acomp(cell_type))
colnames(cell_ilr) = c("ct1", "ct2", "ct3", "ct4", "ct5","ct6")

allConfounders = data.frame("Sample_Group" = pd[,"Sample_Group"], cell_ilr)
formula = paste0(names(allConfounders), collapse = ' + ')
formula = paste0("~", formula)
formula = formula(formula)
#design matrix
design2 = model.matrix(formula, data = allConfounders)
design2[,2] = 1-design2[,2]
colnames(design2)[2] = "Cases"
fit2 = lmFit(beta2m(harmonized_norm), design2)#convert beta to m values
fit22 = eBayes(fit2)
pval22 = p.adjust(fit22$p.value[,2], "BH")
sum(pval22<0.05)

#calculate bacon inflation score
t_stats <- fit22$t[,"Cases"] #random shuffled group
df<- fit22$df.total[1]
bc2 <- bacon(t_stats)
limma_inflation_score2 = inflation(bc2)#0.971


#make a toptable on second column sample_group
col_list2 = sapply(seq(1, dim(fit22$coefficients)[1]), function(i) {
  if (pval22[i]<0.05) {return("red")}
  else {return("black")} })
cpg <- data.frame(logFC=fit22$coefficients[,2],pval=pval22)
sig_hyper <- rownames(cpg[cpg$logFC>0 & cpg$pval<0.05,])
length(sig_hyper)#0
sig_hypo <- rownames(cpg[cpg$logFC<0 & cpg$pval<0.05,])
length(sig_hypo)

png("/nfs/dcmb-lgarmire/xtyang/CordBlood/09-meta_analysis/volcano_3L")
par(cex.axis = 1.5, cex.lab = 1.5)

volcanoplot(fit22,coef=2,col=col_list2,highlight=0,
            names=rownames(fitConfounders$coefficients))
# abline(h=-log10(0.05), lty=2, col="red")
abline(v=-1, lty=2, col="grey")
abline(v=1, lty=2, col="grey")
text(-1, 5, "x=-1", cex=1.2)
text(1, 5, "x=1", cex=1.2)
dev.off()