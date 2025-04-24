#code for figure 5

# Library:
library(lumi)
library(EpiDISH)
library(dplyr)
library(ggplot2)
library(limma)
library(tidyverse)


#load in-house whole blood data
#load pd_all from 2.0 - cell_type_deconvolution_Gervin.R
load("/nfs/dcmb-lgarmire/xtyang/CordBlood/03-wenting_results/pd_all.RData")
names(data_wb)[names(data_wb) == 'pd.Sample_Group'] <- 'Sample_Group'
names(data_wb)[names(data_wb) == 'pd.GA'] <- 'GA'
data_wb$Group <- 'WB'
data_wb_PE <- data_wb[which(data_wb$Sample_Group=='Disease'), ]
data_wb_Controls <- data_wb[which(data_wb$Sample_Group=='Controls'), ]

# Load PBMC data from NCBI:
gset <- getGEO("GSE110828", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL13534", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

pd_450K <- pData(gset)
beta_450K <- exprs(gset)

# Extract the PE and GA:
pd_PBMC <- pd_450K[c('preeclampsia','gestational.age')]
colnames(pd_PBMC) <- c('Sample_Group','GA')
for(i in 1:nrow(pd_PBMC)){
  if(pd_PBMC$Sample_Group[i]==' Yes'){
    pd_PBMC$Sample_Group[i]='Disease'
  }else{pd_PBMC$Sample_Group[i]='Controls'}}
rownames(pd_PBMC) <- pd_450K_sub$geo_accession

# Manipulate cell type proportions and cell type deconvolution references:
## Extract the original PBMC data cell type proportions:
ctp_PBMC <- pd_450K_sub[,c(5:11)]
colnames(ctp_PBMC) <- c('Bcell','CD4T','CD8T','Gran','Mono','NK','nRBC')
rownames(ctp_PBMC) <- rownames(pd_PBMC)

## Omit Gran and re-calculate the cell type proportion in original PBMC data:
ctp_PBMC_nogran <- data.frame(matrix(nrow=nrow(ctp_PBMC), ncol=ncol(ctp_PBMC)-1))
colnames(ctp_PBMC_nogran) <- c('Bcell','CD4T','CD8T','Mono','NK','nRBC')
rownames(ctp_PBMC_nogran) <- rownames(ctp_PBMC)
for (i in 1:nrow(ctp_PBMC_nogran)){
  ctp_PBMC_nogran$Bcell[i] = round(ctp_PBMC$Bcell[i]/(1-ctp_PBMC$Gran[i]),6)
  ctp_PBMC_nogran$CD4T[i] = round(ctp_PBMC$CD4T[i]/(1-ctp_PBMC$Gran[i]),6)
  ctp_PBMC_nogran$CD8T[i] = round(ctp_PBMC$CD8T[i]/(1-ctp_PBMC$Gran[i]),6)
  ctp_PBMC_nogran$Mono[i] = round(ctp_PBMC$Mono[i]/(1-ctp_PBMC$Gran[i]),6)
  ctp_PBMC_nogran$NK[i] = round(ctp_PBMC$NK[i]/(1-ctp_PBMC$Gran[i]),6)
  ctp_PBMC_nogran$nRBC[i] = round(ctp_PBMC$nRBC[i]/(1-ctp_PBMC$Gran[i]),6)
}

data_pbmc <- data.frame(cbind(pd_PBMC$Sample_Group, pd_PBMC$GA, ctp_PBMC_nogran))
names(data_pbmc)[names(data_pbmc) == 'pd_PBMC.Sample_Group'] <- 'Sample_Group'
names(data_pbmc)[names(data_pbmc) == 'pd_PBMC.GA'] <- 'GA'
data_pbmc$Group <- 'PBMC'
data_pbmc_PE <- data_pbmc[which(data_pbmc$Sample_Group=='Disease'), ]
data_pbmc_Controls <- data_pbmc[which(data_pbmc$Sample_Group=='Controls'), ]

data <- rbind(data_wb, data_pbmc)
data_PE <- rbind(data_wb_PE, data_pbmc_PE)
data_Controls <- rbind(data_wb_Controls, data_pbmc_Controls)



##plot figure 5-----------------------------------------------------------
## PE:
pdf("/home/liuwent/04b-cell_type_deconvolution/PE_WB_PBMC_GA_Bcell_scatterplot.pdf")
ggplot(data = data_PE, aes(x = as.numeric(GA), y = as.numeric(Bcell), color=Group)) + 
  geom_point() +
  geom_smooth(aes(group=Group), method='lm') +
  scale_y_continuous(breaks = scales::breaks_extended(n = 5)) + 
  scale_color_manual(values = c("#d95f02", "#7570b3")) + 
  annotate("text", label = "p=0.1280", size = 8, col = "black",  x=Inf, y = Inf, vjust=1, hjust=1)+
  theme(text = element_text(size = 18))+
  labs(x = "GA (PE)", y = "Bcell")
dev.off()

pdf("/home/liuwent/04b-cell_type_deconvolution/PE_WB_PBMC_GA_CD4T_scatterplot.pdf")
ggplot(data = data_PE, aes(x = as.numeric(GA), y = as.numeric(CD4T), color=Group)) + 
  geom_point() +
  geom_smooth(aes(group=Group), method='lm') +
  scale_y_continuous(breaks = scales::breaks_extended(n = 5)) +
  scale_color_manual(values = c("#d95f02", "#7570b3")) +
  annotate("text", label = "p=0.2400", size = 8, col = "black",  x=Inf, y = Inf, vjust=1, hjust=1)+
  theme(text = element_text(size = 18))+
  labs(x = "GA (PE)", y = "CD4T")
dev.off()

pdf("/home/liuwent/04b-cell_type_deconvolution/PE_WB_PBMC_GA_CD8T_scatterplot.pdf")
ggplot(data = data_PE, aes(x = as.numeric(GA), y = as.numeric(CD8T), color=Group)) + 
  geom_point() +
  geom_smooth(aes(group=Group), method='lm') +
  scale_y_continuous(breaks = scales::breaks_extended(n = 5)) +
  scale_color_manual(values = c("#d95f02", "#7570b3")) +
  annotate("text", label = "p=0.4565", size = 8, col = "black",  x=Inf, y = Inf, vjust=1, hjust=1)+
  theme(text = element_text(size = 18))+
  labs(x = "GA (PE)", y = "CD8T")
dev.off()

pdf("/home/liuwent/04b-cell_type_deconvolution/PE_WB_PBMC_GA_Mono_scatterplot.pdf")
ggplot(data = data_PE, aes(x = as.numeric(GA), y = as.numeric(Mono), color=Group)) + 
  geom_point() +
  geom_smooth(aes(group=Group), method='lm') +
  scale_y_continuous(breaks = scales::breaks_extended(n = 5)) +
  scale_color_manual(values = c("#d95f02", "#7570b3")) +
  annotate("text", label = "p=0.8848", size = 8, col = "black",  x=Inf, y = Inf, vjust=1, hjust=1)+
  theme(text = element_text(size = 18))+
  labs(x = "GA (PE)", y = "Mono")
dev.off()

pdf("/home/liuwent/04b-cell_type_deconvolution/PE_WB_PBMC_GA_NK_scatterplot.pdf")
ggplot(data = data_PE, aes(x = as.numeric(GA), y = as.numeric(NK), color=Group)) + 
  geom_point() +
  geom_smooth(aes(group=Group), method='lm') +
  scale_y_continuous(breaks = scales::breaks_extended(n = 5)) +
  scale_color_manual(values = c("#d95f02", "#7570b3")) +
  annotate("text", label = "p=0.6899", size = 8, col = "black",  x=Inf, y = Inf, vjust=1, hjust=1)+
  theme(text = element_text(size = 18))+
  labs(x = "GA (PE)", y = "NK")
dev.off()

pdf("/home/liuwent/04b-cell_type_deconvolution/PE_WB_PBMC_GA_nRBC_scatterplot.pdf")
ggplot(data = data_PE, aes(x = as.numeric(GA), y = as.numeric(nRBC), color=Group)) + 
  geom_point() +
  geom_smooth(aes(group=Group), method='lm') +
  scale_y_continuous(breaks = scales::breaks_extended(n = 5)) +
  scale_color_manual(values = c("#d95f02", "#7570b3")) +
  annotate("text", label = "p=0.3620", size = 8, col = "black",  x=Inf, y = Inf, vjust=1, hjust=1)+
  theme(text = element_text(size = 18))+
  labs(x = "GA (PE)", y = "nRBC")
dev.off()

## Controls:
pdf("/home/liuwent/04b-cell_type_deconvolution/WB_PBMC_GA_Bcell_scatterplot.pdf")
ggplot(data = data_Controls, aes(x = as.numeric(GA), y = as.numeric(Bcell), color=Group)) + 
  geom_point() +
  geom_smooth(aes(group=Group), method='lm') +
  scale_y_continuous(breaks = scales::breaks_extended(n = 5)) + 
  scale_color_manual(values = c("#d95f02", "#7570b3")) +
  annotate("text", label = "p=0.5676", size = 8, col = "black",  x=Inf, y = Inf, vjust=1, hjust=1)+
  theme(text = element_text(size = 18))+
  labs(x = "GA (Control)", y = "Bcell")
dev.off()

pdf("/home/liuwent/04b-cell_type_deconvolution/WB_PBMC_GA_CD4T_scatterplot.pdf")
ggplot(data = data_Controls, aes(x = as.numeric(GA), y = as.numeric(CD4T), color=Group)) + 
  geom_point() +
  geom_smooth(aes(group=Group), method='lm') +
  scale_y_continuous(breaks = scales::breaks_extended(n = 5)) +
  scale_color_manual(values = c("#d95f02", "#7570b3")) +
  annotate("text", label = "p=0.7617", size = 8, col = "black",  x=Inf, y = Inf, vjust=1, hjust=1)+
  theme(text = element_text(size = 18))+
  labs(x = "GA (Control)", y = "CD4T")
dev.off()

pdf("/home/liuwent/04b-cell_type_deconvolution/WB_PBMC_GA_CD8T_scatterplot.pdf")
ggplot(data = data_Controls, aes(x = as.numeric(GA), y = as.numeric(CD8T), color=Group)) + 
  geom_point() +
  geom_smooth(aes(group=Group), method='lm') +
  scale_y_continuous(breaks = scales::breaks_extended(n = 5)) +
  scale_color_manual(values = c("#d95f02", "#7570b3")) +
  annotate("text", label = "p=0.7903", size = 8, col = "black",  x=Inf, y = Inf, vjust=1, hjust=1)+
  theme(text = element_text(size = 18))+
  labs(x = "GA (Control)", y = "CD8T")
dev.off()


pdf("/home/liuwent/04b-cell_type_deconvolution/WB_PBMC_GA_Mono_scatterplot.pdf")
ggplot(data = data_Controls, aes(x = as.numeric(GA), y = as.numeric(Mono), color=Group)) + 
  geom_point() +
  geom_smooth(aes(group=Group), method='lm') +
  scale_y_continuous(breaks = scales::breaks_extended(n = 5)) +
  scale_color_manual(values = c("#d95f02", "#7570b3")) +
  annotate("text", label = "p=0.3712", size = 8, col = "black",  x=Inf, y = Inf, vjust=1, hjust=1)+
  theme(text = element_text(size = 18))+
  labs(x = "GA (Control)", y = "Mono")
dev.off()

pdf("/home/liuwent/04b-cell_type_deconvolution/WB_PBMC_GA_NK_scatterplot.pdf")
ggplot(data = data_Controls, aes(x = as.numeric(GA), y = as.numeric(NK), color=Group)) + 
  geom_point() +
  geom_smooth(aes(group=Group), method='lm') +
  scale_y_continuous(breaks = scales::breaks_extended(n = 5)) +
  scale_color_manual(values = c("#d95f02", "#7570b3")) +
  annotate("text", label = "p=0.3082", size = 8, col = "black",  x=Inf, y = Inf, vjust=1, hjust=1)+
  theme(text = element_text(size = 18))+
  labs(x = "GA (Control)", y = "NK")
dev.off()

pdf("/home/liuwent/04b-cell_type_deconvolution/WB_PBMC_GA_nRBC_scatterplot.pdf")
ggplot(data = data_Controls, aes(x = as.numeric(GA), y = as.numeric(nRBC), color=Group)) + 
  geom_point() +
  geom_smooth(aes(group=Group), method='lm') +
  scale_y_continuous(breaks = scales::breaks_extended(n = 5)) +
  scale_color_manual(values = c("#d95f02", "#7570b3")) +
  annotate("text", label = "p=0.6630", size = 8, col = "black",  x=Inf, y = Inf, vjust=1, hjust=1)+
  theme(text = element_text(size = 18))+
  labs(x = "GA (Control)", y = "nRBC")
dev.off()

# Statistical test for slopes:
## For Bcell:
### Original:
model_Bcell <- lm(data = data, Bcell~as.factor(data$Group)+as.numeric(data$GA))
summary(model_Bcell)
### Disease group only:
model_PE_Bcell <- lm(data = data_PE, Bcell~as.factor(data_PE$Group)+as.numeric(data_PE$GA)+GA*Group)
summary(model_PE_Bcell)
### Control group only:
model_PE_Bcell <- lm(data = data_Controls, Bcell~as.factor(data_Controls$Group)+as.numeric(data_Controls$GA)+GA*Group)
summary(model_PE_Bcell)

## For CD4T:
model_CD4T <- lm(data = data, CD4T~as.factor(data$Group)+as.numeric(data$GA))
summary(model_CD4T)
### Disease group only:
model_PE_CD4T <- lm(data = data_PE, CD4T~as.factor(data_PE$Group)+as.numeric(data_PE$GA)+GA*Group)
summary(model_PE_CD4T)
### Control group only:
model_PE_CD4T <- lm(data = data_Controls, CD4T~as.factor(data_Controls$Group)+as.numeric(data_Controls$GA)+GA*Group)
summary(model_PE_CD4T)

## For CD8T:
model_CD8T <- lm(data = data, CD8T~as.factor(data$Group)+as.numeric(data$GA))
summary(model_CD8T)
### Disease group only:
model_PE_CD8T <- lm(data = data_PE, CD8T~as.factor(data_PE$Group)+as.numeric(data_PE$GA)+GA*Group)
summary(model_PE_CD8T)
### Control group only:
model_PE_CD8T <- lm(data = data_Controls, CD8T~as.factor(data_Controls$Group)+as.numeric(data_Controls$GA)+GA*Group)
summary(model_PE_CD8T)

## For Mono:
model_Mono <- lm(data = data, Mono~as.factor(data$Group)+as.numeric(data$GA))
summary(model_Mono)
### Disease group only:
model_PE_Mono <- lm(data = data_PE, Mono~as.factor(data_PE$Group)+as.numeric(data_PE$GA)+GA*Group)
summary(model_PE_Mono)
### Control group only:
model_PE_Mono <- lm(data = data_Controls, Mono~as.factor(data_Controls$Group)+as.numeric(data_Controls$GA)+GA*Group)
summary(model_PE_Mono)

## For NK:
model_NK <- lm(data = data, NK~as.factor(data$Group)+as.numeric(data$GA))
summary(model_NK)
### Disease group only:
model_PE_NK <- lm(data = data_PE, NK~as.factor(data_PE$Group)+as.numeric(data_PE$GA)+GA*Group)
summary(model_PE_NK)
### Control group only:
model_PE_NK <- lm(data = data_Controls, NK~as.factor(data_Controls$Group)+as.numeric(data_Controls$GA)+GA*Group)
summary(model_PE_NK)

## For nRBC:
model_nRBC <- lm(data = data, nRBC~as.factor(data$Group)+as.numeric(data$GA)+GA*Group)
summary(model_nRBC)
### Disease group only:
model_PE_nRBC <- lm(data = data_PE, nRBC~as.factor(data_PE$Group)+as.numeric(data_PE$GA)+GA*Group)
summary(model_PE_nRBC)
### Control group only:
model_PE_nRBC <- lm(data = data_Controls, nRBC~as.factor(data_Controls$Group)+as.numeric(data_Controls$GA)+GA*Group)
summary(model_PE_nRBC)
