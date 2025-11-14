# Library:
library(pathifier)
library(limma)
library(dbplyr)
library(lumi)
library(EpiDISH)
library(pheatmap)
library(KEGGREST)

# Data:
#load from 5.4 DM_on_gene_level.R
load('/home/liuwent/04d-DM_analysis_after_adj/cpg_geomean_mval.RData')
load('/home/liuwent/04d-DM_analysis_after_adj/fit_genemean00.RData')
load('/home/liuwent/04d-DM_analysis_after_adj/fitConfounders_genemean.RData')
load('/home/liuwent/04d-DM_analysis_after_adj/fitAllConfounders_genemean.RData')
# load from 2.0 - cell_type_deconvolution_Gervin.R
load("/home/liuwent/04-Full_Model/pd_all.RData")

# KEGG Pathway:
kegg_pathway_data = readRDS('/home/liuwent/04d-DM_analysis_after_adj/kegg_pathway_data_FULL_symbol.RDS')

# Take the residuals:
## For model1 - no adjustment:
model1_residuals <- cpg_geomean_mval


## For model2 - fit on clinical variables except Sample_Groups:
allConfounders = pd[, c("GA", "BMI", "Eth2", "Parity", "Age")]
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
model2 = lmFit(cpg_geomean_mval, design2)
model2 = eBayes(model2)
model2_residuals <- residuals(model2, cpg_geomean_mval)
model2_residuals <- model2_residuals+matrix(apply(cpg_geomean_mval, 1, mean), nrow=nrow(model2_residuals), ncol=ncol(model2_residuals))


## For model3 - fit on all (clinical variables and cell types):

allConfounders = pd%>%dplyr::select("GA", "BMI", "Eth2", "Parity", "Age", 
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
model3 = lmFit(cpg_geomean_mval, design2)
model3 = eBayes(model3)
model3_residuals <- residuals(model3, cpg_geomean_mval)
model3_residuals <- model3_residuals+matrix(apply(cpg_geomean_mval, 1, mean), nrow=nrow(model3_residuals), ncol=ncol(model3_residuals))


# pd file:
load("/home/liuwent/04-Full_Model/pd.RData")
for(i in 1:nrow(pd)){
  if(pd$Sample_Group[i] == 'Controls'){
    pd$Sample_Group[i] = 1
  }else{pd$Sample_Group[i] = 0}
}
pd2 <- as.logical(pd$Sample_Group%>%as.numeric%>%-1 +1)
pd2

# pathifier PDS score:
#CpG sites -> genes -> pathifier PDS score -> t-test/ heatmap -> FDR adjustment#
## For model1:
model1_PDS <- quantify_pathways_deregulation(as.matrix(model1_residuals), 
                                             row.names(model1_residuals), 
                                             kegg_pathway_data, 
                                             names(kegg_pathway_data), 
                                             pd2, 
                                             attempts = 5, min_exp = 0, min_std = 0)


## For model2:
model2_PDS <- quantify_pathways_deregulation(as.matrix(model2_residuals), 
                                             row.names(model2_residuals), 
                                             kegg_pathway_data, 
                                             names(kegg_pathway_data), 
                                             pd2, 
                                             attempts = 5, min_exp = 0, min_std = 0)


## For model3:
model3_PDS <- quantify_pathways_deregulation(as.matrix(model3_residuals), 
                                             row.names(model3_residuals), 
                                             kegg_pathway_data, 
                                             names(kegg_pathway_data), 
                                             pd2, 
                                             attempts = 5, min_exp = 0, min_std = 0)


# t-test:
## For model1:
model1_PDS_ttest <- matrix(as.data.frame(model1_PDS$scores),
                           nrow = length(names(model1_PDS$scores)),
                           byrow = TRUE)
colnames(model1_PDS_ttest) <- paste0(pd$Sample_Group,'_',colnames(cpg_geomean_mval))
rownames(model1_PDS_ttest) <- names(model1_PDS$scores)
mode(model1_PDS_ttest) <- "numeric"
dim(model1_PDS_ttest) # 206 62

### Use limma
load("/home/liuwent/04-Full_Model/pd.RData")
design1 = model.matrix(~pd$Sample_Group)
PDS_fit1 = lmFit(model1_PDS_ttest, design1)
PDS_fit1 = eBayes(PDS_fit1)

table1 = topTable(PDS_fit1, num=Inf, coef=2, adjust.method = "BH")
sum(table1$adj.P.Val<0.05) # 200

### Use for loop
controls = model1_PDS_ttest[, grepl('Controls', colnames(model1_PDS_ttest))]
cases = model1_PDS_ttest[, grepl('Disease', colnames(model1_PDS_ttest))]

pathway_ttest_result1 = data.frame()
for(i in 1:nrow(model1_PDS_ttest)){
  path_id = rownames(model1_PDS_ttest)[i]
  pathway_ttest = t.test(cases[i,],controls[i,])
  pathway_ttest_result1[path_id, 1] = pathway_ttest$p.value
  pathway_ttest_result1[path_id, 2] = pathway_ttest$estimate[1]-pathway_ttest$estimate[2]
}
colnames(pathway_ttest_result1) = c('p.val','diff')

pathway_ttest_result1[pathway_ttest_result1$p.val<0.05,]
sum(pathway_ttest_result1$p.val<0.05) # 193
pathway_ttest_result1$adj.p.val = p.adjust(pathway_ttest_result1$p.val, method = 'BH', 
                                           n = length(pathway_ttest_result1$p.val))
sum(pathway_ttest_result1$adj.p.val<0.05) # 190

## For model2:
model2_PDS_ttest <- matrix(as.data.frame(model2_PDS$scores),
                           nrow = length(names(model2_PDS$scores)),
                           byrow = TRUE)
colnames(model2_PDS_ttest) <- paste0(pd$Sample_Group,'_',colnames(cpg_geomean_mval))
rownames(model2_PDS_ttest) <- names(model2_PDS$scores)
mode(model2_PDS_ttest) <- "numeric"
dim(model2_PDS_ttest) # 255 62

### Use limma
design2 = model.matrix(~pd$Sample_Group)
PDS_fit2 = lmFit(model2_PDS_ttest, design2)
PDS_fit2 = eBayes(PDS_fit2)

table2 = topTable(PDS_fit2, num=Inf, coef=2, adjust.method = "BH")
sum(table2$adj.P.Val<0.05) # 0

### Use for loop
controls = model2_PDS_ttest[, grepl('Controls', colnames(model2_PDS_ttest))]
cases = model2_PDS_ttest[, grepl('Disease', colnames(model2_PDS_ttest))]

pathway_ttest_result2 = data.frame()
for(i in 1:nrow(model2_PDS_ttest)){
  path_id = rownames(model2_PDS_ttest)[i]
  pathway_ttest = t.test(cases[i,],controls[i,])
  pathway_ttest_result2[path_id, 1] = pathway_ttest$p.value
  pathway_ttest_result2[path_id, 2] = pathway_ttest$estimate[1]-pathway_ttest$estimate[2]
}
colnames(pathway_ttest_result2) = c('p.val','diff')

pathway_ttest_result2[pathway_ttest_result2$p.val<0.05,]
sum(pathway_ttest_result2$p.val<0.05) # 66
pathway_ttest_result2$adj.p.val = p.adjust(pathway_ttest_result2$p.val, method = 'BH', 
                                           n = length(pathway_ttest_result2$p.val))
sum(pathway_ttest_result2$adj.p.val<0.05) # 0


## For model3:
model3_PDS_ttest <- matrix(as.data.frame(model3_PDS$scores),
                           nrow = length(names(model3_PDS$scores)),
                           byrow = TRUE)
colnames(model3_PDS_ttest) <- paste0(pd$Sample_Group,'_',colnames(cpg_geomean_mval))
rownames(model3_PDS_ttest) <- names(model3_PDS$scores)
mode(model3_PDS_ttest) <- "numeric"
dim(model3_PDS_ttest) # 265 62

### Use limma
design3 = model.matrix(~pd$Sample_Group)
PDS_fit3 = lmFit(model3_PDS_ttest, design3)
PDS_fit3 = eBayes(PDS_fit3)

table3 = topTable(PDS_fit3, num=Inf, coef=2, adjust.method = "BH")
sum(table3$adj.P.Val<0.05) # 0

### Use for loop
controls = model3_PDS_ttest[, grepl('Controls', colnames(model3_PDS_ttest))]
cases = model3_PDS_ttest[, grepl('Disease', colnames(model3_PDS_ttest))]

pathway_ttest_result3 = data.frame()
for(i in 1:nrow(model3_PDS_ttest)){
  path_id = rownames(model3_PDS_ttest)[i]
  pathway_ttest = t.test(cases[i,],controls[i,])
  pathway_ttest_result3[path_id, 1] = pathway_ttest$p.value
  pathway_ttest_result3[path_id, 2] = pathway_ttest$estimate[1]-pathway_ttest$estimate[2]
}
colnames(pathway_ttest_result3) = c('p.val','diff')

pathway_ttest_result3[pathway_ttest_result3$p.val<0.05,]
sum(pathway_ttest_result3$p.val<0.05) # 89
pathway_ttest_result3$adj.p.val = p.adjust(pathway_ttest_result3$p.val, method = 'BH', 
                                           n = length(pathway_ttest_result3$p.val))
sum(pathway_ttest_result3$adj.p.val<0.05) # 0
