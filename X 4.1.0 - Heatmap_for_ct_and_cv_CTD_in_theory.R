# Cell Type Deconvolution Markers:

# Lin Ref.:
load("/home/liuwent/01-data/FlowSorted.CordBlood.EPIC.compTable.rda")
Lin_ref <- FlowSorted.CordBlood.EPIC.compTable

# identified 215,937 cell type specific CpGs at a Bonferroni threshold of 8×10−8
Lin_CT_markers <- Lin_ref[which(Lin_ref$p.value < 8e-8), ]
Lin_CT_markers <- Lin_CT_markers[, 3:8]


# Bakulski Ref.:
load("/home/liuwent/01-data/FlowSorted.CordBlood.450k.compTable.rda")
B_ref <- FlowSorted.CordBlood.450k.compTable

# Cord blood had 211,849 sites that distinguish at least one cell type (P < 10−8)
B_CT_markers <- B_ref[which(B_ref$p.value < 1e-8), ]
B_CT_markers <- B_CT_markers[, 3:9]


# The overlap of two refs markers:
Lin_CT_markers <- Lin_CT_markers[which(rownames(Lin_CT_markers)%in%rownames(B_CT_markers)), ]
Lin_CT_markers <- Lin_CT_markers[order(rownames(Lin_CT_markers)), ]

B_CT_markers <- B_CT_markers[which(rownames(B_CT_markers)%in%rownames(Lin_CT_markers)), ]
B_CT_markers <- B_CT_markers[order(rownames(B_CT_markers)), ]

all.equal(rownames(Lin_CT_markers),rownames(B_CT_markers))
# True

# Combine Bukuski's nRBC to Lin:
ct_ref_markers <- cbind(Lin_CT_markers, B_CT_markers[, 7])
colnames(ct_ref_markers)[7] <- 'nRBC'

save(ct_ref_markers, file="ct_ref_markers.RData")

load("/home/liuwent/04b-cell_type_deconvolution/ct_ref_markers.RData")
# 79049 CpGs, 7 cell types

# Library:
library(ggplot2)
library(EpiDISH)
library(tidyverse)

# Shortcut
load("/home/liuwent/04-Full_Model/myLoad.Rdata")
load("/home/liuwent/04-Full_Model/pd.RData")
load("/home/liuwent/04-Full_Model/myCombat.RData")
load("/home/liuwent/04-Full_Model/Mvalues.RData")
load("/home/liuwent/04-Full_Model/CBref_lin.RData")
load("/home/liuwent/04-Full_Model_New/order.Rdata")
load("/home/liuwent/04-Full_Model/estF.RData")
load("/home/liuwent/04b-cell_type_deconvolution/ct_ref_markers.RData")

#fit for all cell types and clinical variables by CTD in theory:
out<-epidish(beta.m=myCombat, ref.m=as.matrix(CBref),method="CP")
count <- out$estF

count = as.data.frame(count)
count$Sample_Name = rownames(count)
count$id <- 1:nrow(count)
merged <- merge(count, myLoad$pd, by="Sample_Name", sort=F)
order <- merged[order(merged$id),]
save(order, file="order.Rdata")

pd = order
pd = pd%>%dplyr::select("Sample_Group", "GA", "BMI", "Parity", "Age", "Eth2")
pdnames = names(pd)
formstr <- paste0(pdnames, collapse = ' + ')
formstr <- paste0('~', formstr)
formstr <- as.formula(formstr)
design = model.matrix(formstr, data=pd)
fit_CTD_in_theory = lmFit(t(estF), design)
fit_CTD_in_theory = eBayes(fit_CTD_in_theory)
save(fit_CTD_in_theory, file="fit_CTD_in_theory.Rdata")

# Heatmap for cell types and clinical variables:
hm_data <- data.frame(t(-log(fit_CTD_in_theory$p.value[, -1])))

hm_data <- hm_data %>%
           rownames_to_column() %>%
           gather(colname, value, -rowname)
head(hm_data)

pdf("/home/liuwent/04-Full_Model_New/fit_CTD_in_theory.pdf")
ggplot(data = hm_data, aes(x=colname, y=rowname, fill=value)) +
  geom_tile() + 
  scale_fill_gradient(low = "white", high = "#1b98e0") + 
  geom_text(aes(label = round(value, 2))) 
dev.off()


