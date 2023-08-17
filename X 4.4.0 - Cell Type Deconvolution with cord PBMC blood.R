# Library:
library(EpiDISH)
library(ggplot2)
library(reshape2)

# Shortcut
load("/home/liuwent/06-New_Model_AddOn/00-Dataset/pd_450K_sub2.RData")
load("/home/liuwent/06-New_Model_AddOn/00-Dataset/Combat_450K.RData")
load("/home/liuwent/06-New_Model_AddOn/00-Dataset/CBref_450K.RData")

# Cell type reference without gran and nrbc:
CBref_450K <- CBref_450K[, -c(6,7)]

# 4.4.0 - Cell Type Deconvolution with cord PBMC blood
out<-epidish(beta.m=Combat_450K, ref.m=as.matrix(CBref_450K),method="CP")
count <- out$estF

df = data.frame(cbind(count, pd_450K_sub2))

Sample_Group1 = df$Sample_Group
df1 = data.frame(cbind(df[, 1:5], Sample_Group1))

df1_long = melt(df1, id.vars = "Sample_Group1")

pdf("/home/liuwent/04b-cell_type_deconvolution/23-lin_without_gran_nrbc_CP.pdf")
ggplot(df1_long, aes(variable, value, fill=Sample_Group1)) + 
  geom_boxplot() + 
  labs(title="Cell Type Proportion by Sample Group")
dev.off()






