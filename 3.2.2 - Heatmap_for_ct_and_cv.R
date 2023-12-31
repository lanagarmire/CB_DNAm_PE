# Library
library(dplyr)
library(limma)
library(ggplot2)
library(tidyverse)

# Shortcut
load("/home/liuwent/04-Full_Model/pd.RData")
load("/home/liuwent/04b-cell_type_deconvolution/estF.RData")
load("/home/liuwent/04b-cell_type_deconvolution/Residuals.RData")

# fit for all clinical variables by using out$estF:
pd = pd%>%dplyr::select("Sample_Group", "GA", "BMI", "Eth2", "Parity", "Age")
pdnames = names(pd)
formstr <- paste0(pdnames, collapse = ' + ')
formstr <- paste0('~', formstr)
formstr <- as.formula(formstr)
design = model.matrix(formstr, data=pd)
fit_estF_cv = lmFit(t(estF), design)
fit_estF_cv = eBayes(fit_estF_cv)
#save(fit_estF_cv, file="fit_estF_cv.RData")

# Heatmap for cell types and clinical variables:
colnames(fit_estF_cv$p.value)[2] = "Preeclampsia"
colnames(fit_estF_cv$p.value)[5] = "Caucasian"
colnames(fit_estF_cv$p.value)[6] = "Pacific Islander"

hm_data <- data.frame(t(-log10(fit_estF_cv$p.value[, -1])))

hm_data <- hm_data %>% 
  rownames_to_column() %>%
  gather(colname, value, -rowname) %>% 
  setNames(c("Clinical_Variable", "Cell_Type", "Value"))
hm_data$Stars <- cut(hm_data$Value, breaks=c(-Inf, -log10(0.05), -log10(0.01), -log10(0.001), Inf), label=c(" ", "*", "**", "***"))

# Define the order of the factor levels for the x-axis
hm_data$Cell_Type <- factor(hm_data$Cell_Type, 
                                       levels=c("Gran", "NK", "Bcell", "nRBC", "CD4T", "CD8T", "Mono"))

# Define the order of the factor levels for the y-axis
hm_data$Clinical_Variable = factor(hm_data$Clinical_Variable, 
                                   levels=c("Age", "Parity", "Pacific Islander", "Caucasian", "BMI", "GA", "Preeclampsia"))

# for(i in 1:nrow(hm_data)){
#   if(hm_data$Stars[i] != "***"){
#     hm_data$Value[i] = NA
#   }else{hm_data$Value[i] = hm_data$Value[i]}
# }
# 
# hm_data$Value <- as.numeric(hm_data$Value)
# head(hm_data)

pdf("/home/liuwent/04b-cell_type_deconvolution/hm_estF_cv.pdf")
ggplot(data = hm_data, aes(x=Cell_Type, y=Clinical_Variable, fill=Value)) +
  geom_tile()+
  scale_fill_gradient(low = "white", high = "#FF3300")+
  geom_text(aes(label = Stars), col = "#000000")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))+
  guides(fill=guide_legend(title="-log10(pval)"))+
  theme(text = element_text(size = 15))
dev.off()

# # fit for all clinical variables by using residuals:
# pd = pd%>%dplyr::select("Sample_Group", "GA", "BMI", "Eth2", "Parity", "Age")
# pdnames = names(pd)
# formstr <- paste0(pdnames, collapse = ' + ')
# formstr <- paste0('~', formstr)
# formstr <- as.formula(formstr)
# design = model.matrix(formstr, data=pd)
# fit_resi_cv = lmFit(Residuals, design)
# fit_resi_cv = eBayes(fit_resi_cv)
# save(fit_resi_cv, file="fit_resi_cv.RData")
# # load("/home/liuwent/04b-cell_type_deconvolution/fit_resi_cv.RData")
# 
# # Heatmap for cell types and clinical variables:
# hm_data <- data.frame(t(-log10(fit_resi_cv$p.value[, -1])))
# 
# hm_data <- hm_data %>%
#   rownames_to_column() %>%
#   gather(colname, value, -rowname) %>%
#   setNames(c("Clinical_Variable_Residuals", "Cell_Type_Proportion", "Value"))
# hm_data$Stars <- cut(hm_data$Value, breaks=c(-Inf, -log10(0.05), Inf), label=c("", "***"))
# head(hm_data)
# 
# for(i in 1:nrow(hm_data)){
#   if(hm_data$Stars[i] != "***"){
#     hm_data$Value[i] = NA
#   }else{hm_data$Value[i] = hm_data$Value[i]}
# }
# 
# hm_data$Value <- as.numeric(hm_data$Value)
# head(hm_data)
# 
# pdf("/home/liuwent/04b-cell_type_deconvolution/hm_resi_cv.pdf")
# my_colors <- c("#f7fcf0", "#e0f3db", "#ccebc5", "#a8ddb5", 
#                "#7bccc4", "#4eb3d3", "#2b8cbe", "#0868ac", "#084081")
# ggplot(data = hm_data, aes(x=Cell_Type_Proportion, y=Clinical_Variable_Residuals, fill=Value)) +
#   geom_tile() + 
#   # scale_fill_gradient(low = "white", high = "#1b98e0", na.value='grey') +
#   scale_fill_gradient(low = 'white', high = '#1B98E0', na.value='grey') +
#   # geom_text(aes(label = Stars))
#   geom_text(aes(label = round(Value, 2)))
# dev.off()