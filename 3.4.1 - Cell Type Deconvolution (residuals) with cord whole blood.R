# newnewSetofMarkers (151,794) residual version:
# Library:
library(EpiDISH)
library(ggplot2)
library(reshape2)
library(dplyr)
library(limma)
library(ggpubr)
library(ggsignif)

load("/home/liuwent/04-Full_Model/pd.RData")
load("/home/liuwent/04-Full_Model/myCombat.RData")
load("/home/liuwent/04b-cell_type_deconvolution/newnewSetofMarkers_final.RData")

# load("/home/liuwent/04b-cell_type_deconvolution/newnewSetofMarkers_pd.RData")
# load("/home/liuwent/04b-cell_type_deconvolution/newnewSetofMarkers_beta.RData")
# 
# Bcell = rowMeans(newnewSetofMarkers_beta[, which(newnewSetofMarkers_pd$CellType == "Bcell")])
# CD4T = rowMeans(newnewSetofMarkers_beta[, which(newnewSetofMarkers_pd$CellType == "CD4T")])
# CD8T = rowMeans(newnewSetofMarkers_beta[, which(newnewSetofMarkers_pd$CellType == "CD8T")])
# Gran = rowMeans(newnewSetofMarkers_beta[, which(newnewSetofMarkers_pd$CellType == "Gran")])
# Mono = rowMeans(newnewSetofMarkers_beta[, which(newnewSetofMarkers_pd$CellType == "Mono")])
# NK = rowMeans(newnewSetofMarkers_beta[, which(newnewSetofMarkers_pd$CellType == "NK")])
# nRBC = rowMeans(newnewSetofMarkers_beta[, which(newnewSetofMarkers_pd$CellType == "nRBC")])
# 
# # newnewSetofMarkers_final = cbind(Bcell, CD4T, CD8T, Gran, Mono, NK, nRBC)
# newnewSetofMarkers_final = cbind(CD8T, CD4T, NK, Bcell, Mono, Gran, nRBC)
# dim(newnewSetofMarkers_final)
# # 381,923 or 151,794   7
# # save(newnewSetofMarkers_final, file = "newnewSetofMarkers_final.RData")
# # load("/home/liuwent/04b-cell_type_deconvolution/newnewSetofMarkers_final.RData")

out<-epidish(beta.m=myCombat, ref.m=as.matrix(newnewSetofMarkers_final),method="CP")
estF <- out$estF

# fit for all significant confounders except PE by using out$estF:
pd_noPE = pd%>%dplyr::select("GA", "BMI", "Eth2", "Parity", "Age")
pdnames = names(pd_noPE)
formstr <- paste0(pdnames, collapse = ' + ')
formstr <- paste0('~', formstr)
formstr <- as.formula(formstr)
design = model.matrix(formstr, data=pd_noPE)
fit_noPE = lmFit(t(estF), design)
fit_noPE = eBayes(fit_noPE)
# save(fit_noPE, file="fit_noPE.RData")

Residuals_newnew_markers <- residuals(fit_noPE, t(estF))
# Add mean:
Residuals_newnew_markers <- Residuals_newnew_markers+matrix(apply(t(estF), 1, mean),
                                                            nrow=ncol(estF),
                                                            ncol=nrow(estF))
# # Add intercept:
# Residuals_newnew_markers <- Residuals_newnew_markers + fit_noPE$coefficients[,1]
# # save(Residuals_newnew_markers, file="Residuals_newnew_markers.RData")

count <- t(Residuals_newnew_markers)

df = data.frame(cbind(count, pd))

Sample_Group = df$Sample_Group
df1 = data.frame(cbind(df[, 1:7], Sample_Group))

df1_long = melt(df1, id.vars = "Sample_Group")

colnames(df1_long)[1] = "Sample_Group"
colnames(df1_long)[2] = "Cell_Type"
colnames(df1_long)[3] = "Cell_Type_Proportion_Residuals"

# Define the order of the factor levels for the x-axis
df1_long$Cell_Type <- factor(df1_long$Cell_Type, 
                             levels=c("Gran", "NK", "Bcell", "nRBC", "CD4T", "CD8T", "Mono"))

pdf("/home/liuwent/04b-cell_type_deconvolution/23-Residuals_newnewSetofMarkers_CP(151,794).pdf")
ggplot(df1_long, aes(Cell_Type, Cell_Type_Proportion_Residuals, fill=Sample_Group)) + 
  geom_boxplot() + 
  # labs(title="Cell Type Residuals by Sample Group") + 
  # stat_compare_means(aes(label = after_stat(p.signif))) +
  annotate("text", label = c("ns","ns","ns","ns","ns","**","*"), size = 4, col = "black",  x=c(1.1,2.1,3.1,4.1,5.1,6.1,7.1), y = 0.65, vjust=1, hjust=1)+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))+
  theme(text = element_text(size = 15))
dev.off()