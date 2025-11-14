########################################
#01- Ching et al. dataset preprocess
#########################################

#We obtained the dataset from the authors.
# Library:
library(ChAMP)
library(dplyr)
library(lumi)


# 1. Data pre-processing
## 1.1 Load data
myLoad<-champ.load("/nfs/dcmb-lgarmire/liuwent/13-Travers_data/Travers_CB_DNAm_data", arraytype="450k", filterSNPs=FALSE, filterXY=FALSE, filterMultiHit = FALSE)
load("/nfs/dcmb-lgarmire/xtyang/CordBlood/multi_hit_new_450k.RData")#use mask annotation from Zhou et al.
myLoad$beta = myLoad$beta[!rownames(myLoad$beta)%in%multi_hit_new_450k$probeID, ]

newpd<-myLoad$pd[!(myLoad$pd$Sample_Name)=="8",]
newbeta<-myLoad$beta[,!(colnames(myLoad$beta)==8)]

myLoad$pd<-newpd
myLoad$beta<-newbeta

myLoad$pd$Eth2 <- myLoad$pd$EthnicMom1

## Organize ethnicities into bigger categories
for (i in 1:nrow(myLoad$pd)){
  if(myLoad$pd$Eth2[i] %in% c("American Samoan", "Chamorro", "Micronesian","MicronesianChuukese", "MicronesianKosraean", "samoan", "Samoan", "Tongan", "Hawaiian", "Tahitian")){
    myLoad$pd$Eth2[i] = "PacificIslander"}
  
  else if(myLoad$pd$Eth2[i] %in% c("Asian", "Cambodian", "Chinese", "Filipino", "Japanese", "Korean", "Laotian", "Okinawan", "Vietnamese","American Indian")){
    myLoad$pd$Eth2[i] = "Asian"}
  
  else if(myLoad$pd$Eth2[i] %in% c("Caucasian", "English", "German", "Hungarian", "Irish", "Italian", "Portuguese", "Scottish","Spanish")){
    myLoad$pd$Eth2[i] = "Caucasian"}
  
  else{myLoad$pd$Eth2[i] = "Others"}
}

table(myLoad$pd$Eth2)

## Remove two cases which are overlapping with CB_DNAm_PE project:
myLoad$pd <- myLoad$pd[which(myLoad$pd$Same_sample_with_CB != 'y'), ]

myLoad$beta <- myLoad$beta[, which(colnames(myLoad$beta)%in%myLoad$pd$Sample_Name)]


#load("/home/liuwent/13-Travers_data/myLoad.RData")

myLoad$pd$Sample_Name <- as.character(myLoad$pd$Slide)
myLoad$pd$Pool_ID <- as.character(myLoad$pd$Pool_I)
myLoad$pd$Slide <- as.character(myLoad$pd$Slide)
myLoad$pd$MomAge <- as.numeric(myLoad$pd$MomAge)
myLoad$pd$GAWeek <- as.numeric(myLoad$pd$GAWeek)
myLoad$pd$BMI <- as.numeric(myLoad$pd$BMI)


# Normalization(BMIQ)
myNorm<-champ.norm(beta=as.data.frame(myLoad$beta), arraytype="450k", cores = 8)

champ.SVD(beta = as.data.frame(myNorm), PDFplot=TRUE,resultsDir = "/nfs/dcmb-lgarmire/xtyang/CordBlood/03-wenting_results/02-Travis_data")

# save file
save(myLoad, file='/nfs/dcmb-lgarmire/xtyang/CordBlood/03-wenting_results/02-Travis_data/myLoad2.RData')
save(myNorm, file = '/nfs/dcmb-lgarmire/xtyang/CordBlood/03-wenting_results/02-Travis_data/myNorm2.RData')
pd = myLoad$pd
save(pd, file = '/nfs/dcmb-lgarmire/xtyang/CordBlood/03-wenting_results/02-Travis_data/pd2.RData')
