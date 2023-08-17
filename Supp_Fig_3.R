# Library
library(dplyr)
library(limma)
library(ggplot2)
library(tidyverse)

# Shortcut
load("/home/liuwent/04-Full_Model/pd.RData")
load("/home/liuwent/04b-cell_type_deconvolution/estF.RData")
load("/home/liuwent/04b-cell_type_deconvolution/Residuals.RData")

pd2 = cbind(pd, estF)
pe_ga_long = pivot_longer(pd2, 15:21, names_to = "cell_Type", values_to = "cell_Type_Proportion")

# nRBC and BMI
pdf("/home/liuwent/04b-cell_type_deconvolution/BMI_nRBC.pdf")
ggplot(pe_ga_long%>%filter(cell_Type == "nRBC"), aes(x = BMI, y = cell_Type_Proportion, col = BMI))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_gradient(low = "blue", high = "#FF3300")+
  annotate("text", label = "p=0.0067", size = 5, x = Inf, y = Inf, col = "black", vjust=1, hjust=1)+
  theme(text = element_text(size = 18))+
  xlab("BMI") + 
  ylab("nRBC Proportion")
dev.off()

# CD8T and BMI
pdf("/home/liuwent/04b-cell_type_deconvolution/BMI_CD8T.pdf")
ggplot(pe_ga_long%>%filter(cell_Type == "CD8T"), aes(x = BMI, y = cell_Type_Proportion, col = BMI))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_color_gradient(low = "blue", high = "#FF3300")+
  annotate("text", label = "p=0.002", size = 5, x = Inf, y = Inf, col = "black", vjust=1, hjust=1)+
  theme(text = element_text(size = 18))+
  xlab("BMI") + 
  ylab("CD8T Proportion")
dev.off()

# Pacific Islander vs Cell Types:
library(EpiDISH)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(ggsignif)

load("/home/liuwent/04-Full_Model/pd.RData")
load("/home/liuwent/04-Full_Model/myCombat.RData")
load("/home/liuwent/04b-cell_type_deconvolution/newnewSetofMarkers_final.RData")

out<-epidish(beta.m=myCombat, ref.m=as.matrix(newnewSetofMarkers_final),method="CP")
count <- out$estF

df = data.frame(cbind(count, pd))

# Encoding Pacific Islander:
df$PI <- ifelse(df$Eth2=="PacificIslander", 'Y', 'N')

PI = df$PI
df1 = data.frame(cbind(df[, 1:7], PI))

df1_long = melt(df1, id.vars = "PI")

colnames(df1_long)[1] = "Pacific_Islander"
colnames(df1_long)[2] = "Cell_Type"
colnames(df1_long)[3] = "Cell_Type_Proportion"

# Define the order of the factor levels for the x-axis
df1_long$Cell_Type <- factor(df1_long$Cell_Type, 
                             levels=c("Gran", "NK", "Bcell", "nRBC", "CD4T", "CD8T", "Mono"))

pdf("/home/liuwent/04b-cell_type_deconvolution/PI_CT.pdf")
ggplot(df1_long, aes(Cell_Type, Cell_Type_Proportion, fill=Pacific_Islander)) + 
  geom_boxplot() +
  annotate("Text", label = "**", x = 4, y = 0.6)+
  annotate("Text", label = "*", x = 6, y = 0.6)+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))+
  theme(text = element_text(size = 15))+
  guides(fill=guide_legend(title="Pacific Islander"))
dev.off()