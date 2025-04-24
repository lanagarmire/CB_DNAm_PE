# cell type deconvolution
# Xiaotong Yang

#cell type deconvolution is done using Houseman's CP from EpiDISH package, cord blood cell reference is obtained from Gervin et al.
#Gervin, K., Salas, L.A., Bakulski, K.M. et al. Systematic evaluation and validation of reference and library selection methods for 
#deconvolution of cord blood DNA methylation data. Clin Epigenet 11, 125 (2019). https://doi.org/10.1186/s13148-019-0717-y
library(EpiDISH)

#load phenotype data and batch-corrected beta matrix
# See 1 - Raw data pre-processing.Rmd
load("/nfs/dcmb-lgarmire/liuwent/04-Full_Model/pd.RData")
load("/nfs/dcmb-lgarmire/liuwent/04-Full_Model/myCombat.RData")

#load Gervins et al. cord blood reference
# see https://github.com/immunomethylomics/FlowSorted.CordBloodCombined.450k
load("/nfs/dcmb-lgarmire/xtyang/CordBlood/03-wenting_results/FlowSorted.CordBloodCombined.450k.compTable.rda")
ref = FlowSorted.CordBloodCombined.450k.compTable
out<-epidish(beta.m=myCombat, ref.m=as.matrix(ref), method="CP")
count <- out$estF
pd = data.frame(cbind(out$estF, pd))
save(pd, file = "pd_all.RData")