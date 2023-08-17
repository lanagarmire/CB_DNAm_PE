# Library:
library(minfi)
library(limma)
library(lumi)

# Data Shortcut:
load("/home/liuwent/13-Travers_data/pd.RData")
load("/home/liuwent/13-Travers_data/m.RData")
load("/home/liuwent/13-Travers_data/beta.RData")


path <- "/home/liuwent/13-Travers_data/Travers_CB_DNAm_data"
list.files(path)

targets <- read.csv("/home/liuwent/13-Travers_data/Travers_CB_DNAm_data/Trav_pd.csv", as.is = TRUE)
dim(targets) # 21 17

# Remove one case due to the outlier birth weight (mentioned in Travers paper)
targets<-targets[!(targets$Sample_Name)=="8",]
dim(targets) # 20 17

## Organize ethnicities into bigger categories
for (i in 1:nrow(targets)){
  if(targets$EthnicMom1[i] %in% c("American Samoan", "Chamorro", "Micronesian","MicronesianChuukese", "MicronesianKosraean", "samoan", "Samoan", "Tongan", "Hawaiian", "Tahitian")){
    targets$EthnicMom1[i] = "PacificIslander"}
  
  else if(targets$EthnicMom1[i] %in% c("Asian", "Cambodian", "Chinese", "Filipino", "Japanese", "Korean", "Laotian", "Okinawan", "Vietnamese","American Indian")){
    targets$EthnicMom1[i] = "Asian"}
  
  else if(targets$EthnicMom1[i] %in% c("Caucasian", "English", "German", "Hungarian", "Irish", "Italian", "Portuguese", "Scottish","Spanish")){
    targets$EthnicMom1[i] = "Caucasian"}
  
  else{targets$EthnicMom1[i] = "Others"}
}

table(targets$EthnicMom1)

targets$Basename <- gsub(" ", "", paste(targets$Sentrix_ID, "_", targets$Sentrix_Position))

targets$Basename <- file.path(path, targets$Basename)
rgset <- read.metharray(targets$Basename, verbose = TRUE)

pd<-targets
dim(pd) # 20 18
save(pd, file = "/home/liuwent/13-Travers_data/pd.RData")
load("/home/liuwent/13-Travers_data/pd.RData")

dim(getRed(rgset)) # 622,399     20
dim(getGreen(rgset)) # 622,399     20

qcReport(rgset, sampNames=pd$Sample_Name, sampGroups=pd$Sample_Group, pdf="qcReport.pdf")

MSet.raw<-preprocessRaw(rgset)

MSet.norm<-preprocessIllumina(rgset, bg.correct=TRUE, normalize="controls", reference=2)

getMeth(MSet.raw)[1:4,1:3]

getUnmeth(MSet.raw)[1:4,1:3]

beta <- getBeta(MSet.raw,type="Illumina")
dim(beta) # 485512     20
sum(is.na(beta)) # 0
save(beta, file = "/home/liuwent/13-Travers_data/beta.RData")
load("/home/liuwent/13-Travers_data/beta.RData")

m <- getM(MSet.raw)
dim(m) # 485512     20
sum(is.na(m)) # 501
m <- na.omit(m)
dim(m) # 485278     20
save(m, file = "/home/liuwent/13-Travers_data/m.RData")
load("/home/liuwent/13-Travers_data/m.RData")

# 2. regress only on sample group-----------------------------
pd$Sample_Group = relevel(factor(pd$Sample_Group), ref="CONTROL")
design0 = model.matrix(~pd$Sample_Group)
fit0 = lmFit(m, design0)
fit00 = eBayes(fit0)

## Get log_2 fold change and p value
log2fc = fit00$coefficients[,2]
##pval = p.adjust(fit$p.value[,2], method="bonferroni")
pval0 = p.adjust(fit00$p.value[,2], "BH")
length(pval0[pval0<0.05]) # 68,458
cpg <- data.frame(logFC=log2fc, pval=pval0)
sig_hyper <- rownames(cpg[cpg$logFC>0 & cpg$pval<0.05,])
length(sig_hyper) # 16,256
sig_hypo <- rownames(cpg[cpg$logFC<0 & cpg$pval<0.05,])
length(sig_hypo) # 52,202
table_fit00 = topTable(fit00, num=Inf, coef=2, adjust.method = "BH")
table_fit00$Type = ifelse(table_fit00$logFC>0, "Hyper", "Hypo")
sum(table_fit00$adj.P.Val<0.05) # 68,458

## Make a volcano plot (label DMP with different colors)
col_list = sapply(seq(1, dim(fit00$coefficients)[1]), function(i) {
  if (pval0[i]<0.05) {return("red")}
  else {return("black")} })
table(col_list)

pdf(file = "/home/liuwent/13-Travers_data/01-Trav_fit00_minfi.pdf")
par(cex.axis = 1.5, cex.lab = 1.5) 

volcanoplot(fit00, coef=2, col=col_list, highlight=5,
            names=rownames(fit00$coefficients))
# abline(h=-log10(0.05), lty=2, col="red")
abline(v=-1, lty=2, col="grey")
abline(v=1, lty=2, col="grey")
text(-1, 5, "x=-1", cex=1.2)
text(1, 5, "x=1", cex=1.2)
dev.off()

# 3. Cell type deconvolution
## 3.1 Regular version:
library(EpiDISH)
library(ggplot2)
library(reshape2)
library(dplyr)
library(limma)
library(ggpubr)
library(ggsignif)

load("/home/liuwent/04b-cell_type_deconvolution/newnewSetofMarkers_final.RData")

out<-epidish(beta.m=beta, ref.m=as.matrix(newnewSetofMarkers_final), method="CP")
Trav_estF <- out$estF

df = data.frame(cbind(Trav_estF, pd))
Trav_pd_all_minfi <- df
dim(Trav_pd_all_minfi) # 20 25
save(Trav_pd_all_minfi, file = "/home/liuwent/13-Travers_data/Trav_pd_all_minfi.RData")
load("/home/liuwent/13-Travers_data/Trav_pd_all_minfi.RData")

Sample_Group = df$Sample_Group
df1 = data.frame(cbind(df[, 1:7], Sample_Group))

df1_long = melt(df1, id.vars = "Sample_Group")

colnames(df1_long)[1] = "Sample_Group"
colnames(df1_long)[2] = "Cell_Type"
colnames(df1_long)[3] = "Cell_Type_Proportion"

# Define the order of the factor levels for the x-axis
df1_long$Cell_Type <- factor(df1_long$Cell_Type, 
                             levels=c("Gran", "NK", "Bcell", "nRBC", "CD4T", "CD8T", "Mono"))

pdf("/home/liuwent/13-Travers_data/02-Trav_newnewSetofMarkers_CP(151,794)_minfi.pdf")
ggplot(df1_long, aes(Cell_Type, Cell_Type_Proportion, fill=Sample_Group)) + 
  geom_boxplot() +
  # labs(title="Cell Type Proportion by Sample Group") + 
  stat_compare_means(aes(label = after_stat(p.signif))) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))+
  theme(text = element_text(size = 15))
dev.off()

## 3.2 Relationship between phenotypes and cell types
# Library
library(dplyr)
library(limma)
library(ggplot2)
library(tidyverse)

# fit for all clinical variables by using Trav_estF:
pd = pd%>%dplyr::select("Sample_Group", "GAWeek", "MomAge", "BabySex")
pdnames = names(pd)
formstr <- paste0(pdnames, collapse = ' + ')
formstr <- paste0('~', formstr)
formstr <- as.formula(formstr)
design = model.matrix(formstr, data=pd)
fit_estF_cv = lmFit(t(Trav_estF), design)
fit_estF_cv = eBayes(fit_estF_cv)

hm_data <- data.frame(t(-log10(fit_estF_cv$p.value[, -1])))

hm_data <- hm_data %>% 
  rownames_to_column() %>%
  gather(colname, value, -rowname) %>% 
  setNames(c("Clinical_Variable", "Cell_Type", "Value"))
hm_data$Stars <- cut(hm_data$Value, breaks=c(-Inf, -log10(0.05), -log10(0.01), -log10(0.001), Inf), label=c(" ", "*", "**", "***"))

# Define the order of the factor levels for the x-axis
hm_data$Cell_Type <- factor(hm_data$Cell_Type, 
                            levels=c("Gran", "NK", "Bcell", "nRBC", "CD4T", "CD8T", "Mono"))

pdf("/home/liuwent/13-Travers_data/03a-Trav_hm_estF_cv_minfi.pdf")
ggplot(data = hm_data, aes(x=Cell_Type, y=Clinical_Variable, fill=Value)) +
  geom_tile()+
  scale_fill_gradient(low = "white", high = "#FF3300")+
  geom_text(aes(label = Stars), col = "#000000")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))+
  guides(fill=guide_legend(title="-log10(pval)"))+
  theme(text = element_text(size = 15))
dev.off()

## 3.3 Residual version:
# fit for all significant confounders except PE by using Trav_estF:
pd_noPE = pd%>%dplyr::select("GAWeek", "MomAge", "BabySex")
pdnames = names(pd_noPE)
formstr <- paste0(pdnames, collapse = ' + ')
formstr <- paste0('~', formstr)
formstr <- as.formula(formstr)
design = model.matrix(formstr, data=pd_noPE)
fit_noPE = lmFit(t(Trav_estF), design)
fit_noPE = eBayes(fit_noPE)

Residuals_newnew_markers <- residuals(fit_noPE, t(Trav_estF))
# Add mean:
Residuals_newnew_markers <- Residuals_newnew_markers+matrix(apply(t(Trav_estF), 1, mean),
                                                            nrow=ncol(Trav_estF),
                                                            ncol=nrow(Trav_estF))
count <- t(Residuals_newnew_markers)

df = data.frame(cbind(count, pd))

Sample_Group = as.factor(df$Sample_Group)
df1 = data.frame(cbind(df[, 1:7], Sample_Group))

df1_long = melt(df1, id.vars = "Sample_Group")

colnames(df1_long)[1] = "Sample_Group"
colnames(df1_long)[2] = "Cell_Type"
colnames(df1_long)[3] = "Cell_Type_Proportion_Residuals"

# Define the order of the factor levels for the x-axis
df1_long$Cell_Type <- factor(df1_long$Cell_Type, 
                             levels=c("Gran", "NK", "Bcell", "nRBC", "CD4T", "CD8T", "Mono"))

pdf("/home/liuwent/13-Travers_data/03b-Trav_Residuals_newnewSetofMarkers_CP(151,794)_minfi.pdf")
ggplot(df1_long, aes(Cell_Type, Cell_Type_Proportion_Residuals, fill=Sample_Group)) + 
  geom_boxplot() + 
  # labs(title="Cell Type Residuals by Sample Group") + 
  # stat_compare_means(aes(label = after_stat(p.signif))) +
  annotate("text", label = c("ns","ns","ns","ns","ns","ns","ns"), size = 4, col = "black",  x=c(1.1,2.1,3.1,4.1,5.1,6.1,7.1), y = 0.45, vjust=1, hjust=1)+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))+
  theme(text = element_text(size = 15))
dev.off()

# 4. SOV
# Library:
library(shiny)
library(ggplot2)
library(parallel)
library(dplyr)
library(lumi)
library(minfi)
library(RColorBrewer)
library(MASS)
library(limma)
library(ggcorrplot)
library(carData)
library(car)

# Data:
load("/home/liuwent/13-Travers_data/Trav_pd_all_minfi.RData")
load("/home/liuwent/13-Travers_data/m.RData")
load("/home/liuwent/13-Travers_data/beta.RData")

#extract confounders
pd = Trav_pd_all_minfi
pd = pd%>%dplyr::select("Sample_Group", "GAWeek", "MomAge", "BMI", "EthnicMom1", "BabySex", "BabyWgtGram",
                        "CD8T", "CD4T", "NK", "Bcell", "Mono", "Gran", "nRBC")
pd$Sample_Group = as.numeric(pd$Sample_Group)
pd$EthnicMom1 = as.numeric(as.factor(pd$EthnicMom1))

# sovdat = data.frame(t(m))
sovdat = data.frame(t(beta))
# sovdatlist <- list()
# for(i in 1:ncol(sovdat)){
#   sovdatlist[[i]] <- sovdat[,i]
# }

Ftab <- data.frame(Sample_Group = numeric(), stringsAsFactors = FALSE)

for(i in 2:ncol(pd)){
  varname <- names(pd)[i]
  Ftab[varname] <- numeric()
}

calF <- function(probe = probecol){
  newdata <- pd
  pdnames <- names(newdata)
  newdata$beta <- probe
  formstr <- paste0(pdnames, collapse = ' + ')
  formstr <- paste0('beta ~ ', formstr)
  formstr <- as.formula(formstr)
  fit <- lm(formstr, data = newdata)
  aovfit <- Anova(fit, type = 3, singular.ok = TRUE)
  F <- aovfit$`F value`
  F <- F[2:(length(F)-1)]
  names(F) <- pdnames
  F <- as.data.frame(F, stringsAsFactors = FALSE)
  F <- as.data.frame(t(F))
  row.names(F) <- 1
  Ftab <- rbind(Ftab, F)
  return(Ftab)
}

Ftab <- mclapply(X = sovdat, FUN = calF, mc.cores = 40)
Ftab <- do.call(rbind, Ftab)
Fmean <- colMeans(Ftab)
Fmean <- Fmean[order(-Fmean)]
Fmean <- data.frame(Factor = names(Fmean), Fstat = as.vector(Fmean), stringsAsFactors = FALSE)
finalvars <- unique(c('Sample_Group', Fmean$Factor[Fmean$Fstat > 1]))

save(Ftab, Fmean, finalvars, file = "/home/liuwent/13-Travers_data/Trav_res_minfi.RData")

load("/home/liuwent/13-Travers_data/Trav_res_minfi.RData")

sovplot <- function(restab = MSSmean, clustername = 'Case', plottype = 'MSS',
                    textsize = 20){
  resmean <- restab
  samplegroupidx <- match('Sample_Group', resmean$Factor)
  resmean$Factor[samplegroupidx] <- paste0(clustername, '_Control')
  if(plottype == 'MSS'){
    ytitle <- 'Mean Square'
    resmean <- resmean[order(-resmean$MSSstat),]
    resmean$Factor <- factor(resmean$Factor, levels = resmean$Factor, ordered = TRUE)
    p <- ggplot(data = resmean, mapping = aes(x = Factor, y = MSSstat, fill = Factor))
    print(
      p + geom_bar(stat = 'identity') +
        ggtitle('Source of Variance (Type 3 Anova)') +
        ylab(ytitle) +
        xlab('') +
        scale_fill_discrete(guide = FALSE) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 60, hjust = 1, size = textsize)) +
        theme(axis.text.y = element_text(size = textsize))
    )
  }else if(plottype == 'pval'){
    ytitle <- '-log2(p-val)'
    resmean <- resmean[order(-resmean$logpval),]
    resmean$Factor <- factor(resmean$Factor, levels = resmean$Factor, ordered = TRUE)
    
    p <- ggplot(data = resmean, mapping = aes(x = Factor, y = logpval, fill = Factor))
    print(
      p + geom_bar(stat = 'identity') +
        ggtitle('Source of Variance (Type 3 Anova)') +
        ylab(ytitle) +
        xlab('') +
        scale_fill_discrete(guide = FALSE) +
        geom_hline(yintercept = -log2(0.05), color = 'red', size = 1) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 60, hjust = 1, size = textsize)) +
        theme(axis.text.y = element_text(size = textsize))
    )
  }else{
    ytitle <- 'F statistic'
    resmean <- resmean[order(-resmean$Fstat),]
    resmean$Factor <- factor(resmean$Factor, levels = resmean$Factor, ordered = TRUE)
    p <- ggplot(data = resmean, mapping = aes(x = Factor, y = Fstat, fill = Factor))
    print(
      p + geom_bar(stat = 'identity') +
        ggtitle('Source of Variance (Type 3 Anova)') +
        ylab(ytitle) +
        xlab('') +
        scale_fill_discrete(guide = FALSE) +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 60, hjust = 1, size = textsize)) +
        theme(axis.text.y = element_text(size = textsize))+
        geom_hline(yintercept=1, linetype="dashed", color = "red")
    )
  }
}

pdf("/home/liuwent/13-Travers_data/04_Trav_sov_plot_for_all_minfi.pdf")
sovplot(restab = Fmean, plottype = 'F', textsize = 7)
dev.off()

Fmean
#The variables that are over 1 are the confounding factors you will address.

# 5. Adjust for clinical confounders
pd = Trav_pd_all_minfi
allConfounders = pd%>%dplyr::select("Sample_Group", "GAWeek", "MomAge", "BabySex")

ind <- sapply(allConfounders, is.numeric)
#scale numerical variables
f = function(x){scale(x, center = FALSE)}
allConfounders[ind] <- lapply(allConfounders[ind],f)
formula = paste0(names(allConfounders), collapse = ' + ')
formula = paste0("~", formula)
formula = formula(formula)

#design matrix
design2 = model.matrix(formula, data = allConfounders)
#???#default setting control = 1, need to set case = 1???#
# design2[,2] = 1-design2[,2]
colnames(design2)[2] = "Cases"
#fit linear model
fitConfounders = lmFit(m, design2)
fitConfounders = eBayes(fitConfounders)

#make a toptable on second column sample_group
allg.limma <- topTable(fitConfounders, coef="Cases", n=dim(fitConfounders)[1], adjust.method="fdr")
sigg.limma <- subset(allg.limma, adj.P.Val < 0.05)
nonsigg.limma <- subset(allg.limma, adj.P.Val >= 0.05)
dim(sigg.limma) # 0 # 51,486 (reported by Travers paper)

# Make a volcano plot (label DMP with different colors)
pval2 = p.adjust(fitConfounders$p.value[,2], method = "fdr")
col_list2 = sapply(seq(1, dim(fitConfounders$coefficients)[1]), function(i) {
  if (pval2[i]<0.05) {return("red")}
  else {return("black")} })
cpg <- data.frame(logFC=fitConfounders$coefficients[,2],pval=pval2)
sig_hyper <- rownames(cpg[cpg$logFC>0 & cpg$pval<0.05,])
length(sig_hyper) # 0 # 12,563 (reported by Travers paper)
sig_hypo <- rownames(cpg[cpg$logFC<0 & cpg$pval<0.05,])
length(sig_hypo) # 0

pdf("/home/liuwent/13-Travers_data/05-Trav_fitConfounders_cv_minfi.pdf")
par(cex.axis = 1.5, cex.lab = 1.5)

volcanoplot(fitConfounders,coef=2,col=col_list2,highlight=1,
            names=rownames(fitConfounders$coefficients))
# abline(h=-log10(0.05), lty=2, col="red")
abline(v=-1, lty=2, col="grey")
abline(v=1, lty=2, col="grey")
text(-1, 5, "x=-1", cex=1.2)
text(1, 5, "x=1", cex=1.2)
dev.off()

# 6. Adjust for all confounders
pd = Trav_pd_all_minfi
allConfounders = pd%>%dplyr::select("Sample_Group", "GAWeek", "MomAge", "BabySex",
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
#???#default setting control = 1, need to set case = 1???#
# design2[,2] = 1-design2[,2]
colnames(design2)[2] = "Cases"
#fit linear model
fitConfounders = lmFit(m, design2)
fitConfounders = eBayes(fitConfounders)

#make a toptable on second column sample_group
allg.limma <- topTable(fitConfounders, coef="Cases", n=dim(fitConfounders)[1], adjust.method="fdr")
sigg.limma <- subset(allg.limma, adj.P.Val < 0.05)
nonsigg.limma <- subset(allg.limma, adj.P.Val >= 0.05)
dim(sigg.limma) # 1 # 51,486 (reported by Travers paper)

# Make a volcano plot (label DMP with different colors)
pval2 = p.adjust(fitConfounders$p.value[,2], method = "fdr")
col_list2 = sapply(seq(1, dim(fitConfounders$coefficients)[1]), function(i) {
  if (pval2[i]<0.05) {return("red")}
  else {return("black")} })
cpg <- data.frame(logFC=fitConfounders$coefficients[,2],pval=pval2)
sig_hyper <- rownames(cpg[cpg$logFC>0 & cpg$pval<0.05,])
length(sig_hyper) # 1 (cg18827503) # 12,563 (reported by Travers paper)
sig_hypo <- rownames(cpg[cpg$logFC<0 & cpg$pval<0.05,])
length(sig_hypo) # 0

pdf("/home/liuwent/13-Travers_data/06-Trav_fitConfounders_cvct_minfi.pdf")
par(cex.axis = 1.5, cex.lab = 1.5)

volcanoplot(fitConfounders,coef=2,col=col_list2,highlight=1,
            names=rownames(fitConfounders$coefficients))
# abline(h=-log10(0.05), lty=2, col="red")
abline(v=-1, lty=2, col="grey")
abline(v=1, lty=2, col="grey")
text(-1, 5, "x=-1", cex=1.2)
text(1, 5, "x=1", cex=1.2)
dev.off()