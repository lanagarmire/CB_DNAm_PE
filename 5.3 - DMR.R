library(bumphunter)
library(minfi)
library(dplyr)
library(doParallel)

#load myCombat from 1 - Raw data pre-processing.Rmd 
#load pd_all from 2.0 - cell_type_deconvolution_Gervin.R
load("/nfs/dcmb-lgarmire/liuwent/04-Full_Model/myCombat.RData")
load("/nfs/dcmb-lgarmire/xtyang/CordBlood/03-wenting_results/pd_all.RData")

#design matrix
pd = df
allConfounders = pd%>%dplyr::select("Sample_Group", "GA", "Age", "Parity", "Eth2", 
                                    "CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran","Smoker")
ind <- sapply(allConfounders, is.numeric)
#scale numerical variables
f = function(x){scale(x, center = FALSE)}
allConfounders[ind] <- lapply(allConfounders[ind],f)
formula = paste0(names(allConfounders), collapse = ' + ')
formula = paste0("~", formula)
formula = formula(formula)
design = model.matrix(formula, data = allConfounders)

registerDoParallel(cores = 8)

#get annotation
RSobject <- RatioSet(myCombat, annotation = c(array = "IlluminaHumanMethylationEPIC",annotation = "ilm10b4.hg19"))
probe.features <- getAnnotation(RSobject)
cpg.idx <- intersect(rownames(myCombat),rownames(probe.features))
Anno <- probe.features[cpg.idx,]
Anno <- Anno[order(Anno$chr,Anno$pos),]
cpg.idx <- rownames(Anno)


#identify clusters
minProbes=7
maxGap=300
cl <- clusterMaker(Anno$chr,Anno$pos,maxGap=maxGap)
names(cl) <- cpg.idx
bumphunter.idx <- cpg.idx[which(cl %in% names(which(table(cl)>minProbes)))]
sum(table(cl)>minProbes)

#convert beta to mvalues
Beta <- myCombat[bumphunter.idx,]
Beta <- replace(Beta,which(Beta <= 0.001),0.001)
Beta <- replace(Beta,which(Beta >= 0.999),0.999)
Y <- log((Beta/(1-Beta)),2)

dmr_results <- bumphunter(
  Y,
  design = design,
  chr=Anno[bumphunter.idx,]$chr,
  pos=Anno[bumphunter.idx,]$pos,
  cluster=cl[bumphunter.idx],
  cutoff=0.2,
  B=250,
  nullMethod = "bootstrap"
)
 
save(dmr_results, file = "/nfs/dcmb-lgarmire/xtyang/CordBlood/03-wenting_results/01-new_result_Garvin_ref/dmr_results.RData")

load("01-new_result_Garvin_ref/dmr_results")
adjPvalDmr = 0.05
DMR <- dmr_results$table[which(dmr_results$table$p.valueArea <= adjPvalDmr),]

rownames(DMR) <- paste("DMR",1:nrow(DMR),sep="_")
DMR <- data.frame(DMR[,1:3],width=DMR[,3]-DMR[,2],strand="*",DMR[,4:14])
colnames(DMR)[1:3] <- c("seqnames","start","end") 

write.csv(DMR, file = "/nfs/dcmb-lgarmire/xtyang/CordBlood/03-wenting_results/01-new_result_Garvin_ref/DMR.csv")
