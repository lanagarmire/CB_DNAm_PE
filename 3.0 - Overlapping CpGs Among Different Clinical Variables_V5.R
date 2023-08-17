# Library
library(limma)

# Shortcut
load("/home/liuwent/04-Full_Model/pd.RData")
load("/home/liuwent/04-Full_Model/Mvalues.RData")

# Fit on all significant confounders
allConfounders = pd[, c("Sample_Group", "GA", "BMI", "Eth2", "Parity", "Age")]
ind <- sapply(allConfounders, is.numeric)

#scale numerical variables
f = function(x){scale(x, center = FALSE)}
allConfounders[ind] <- lapply(allConfounders[ind],f)

# # Change Eth2 reference
# allConfounders$Eth2 = relevel(factor(allConfounders$Eth2), ref="Others")

#design matrix
formula = paste0(names(allConfounders), collapse = ' + ')
formula = paste0("~", formula)
formula = formula(formula)

design = model.matrix(formula, data = allConfounders)

colnames(design)[2] = "Cases"

# Fit on all significant confounders
fit = lmFit(Mvalues, design)
fit = eBayes(fit)

# PE related significant CpGs:
pval_PE <- p.adjust(fit$p.value[,2], "BH")
pval_PE_sig <- pval_PE[pval_PE<0.05]
length(pval_PE_sig)
# 67
write.csv(pval_PE_sig, file='pval_PE_sig.csv')

# GA related significant CpGs:
pval_GA <- p.adjust(fit$p.value[,3], "BH")
pval_GA_sig <- pval_GA[pval_GA<0.05]
length(pval_GA_sig)
# 63,946
write.csv(pval_GA_sig, file='pval_GA_sig.csv')

# BMI related significant CpGs:
pval_BMI <- p.adjust(fit$p.value[,4], "BH")
pval_BMI_sig <- pval_BMI[pval_BMI<0.05]
length(pval_BMI_sig)
# 1
write.csv(pval_BMI_sig, file='pval_BMI_sig.csv')

# Eth2 related significant CpGs:
## Eth2Caucasian related significant CpGs:
pval_Eth2Caucasian <- p.adjust(fit$p.value[,5], "BH")
pval_Eth2Caucasian_sig <- pval_Eth2Caucasian[pval_Eth2Caucasian<0.05]
length(pval_Eth2Caucasian_sig)
# 0

## Eth2Pacific Islander related significant CpGs:
pval_Eth2PacificIslander <- p.adjust(fit$p.value[,6], "BH")
pval_Eth2PacificIslander_sig <- pval_Eth2PacificIslander[pval_Eth2PacificIslander<0.05]
length(pval_Eth2PacificIslander_sig)
# 0

# Parity related significant CpGs:
pval_Parity <- p.adjust(fit$p.value[,7], "BH")
pval_Parity_sig <- pval_Parity[pval_Parity<0.05]
length(pval_Parity_sig)
# 0

# Age related significant CpGs:
pval_Age <- p.adjust(fit$p.value[,8], "BH")
pval_Age_sig <- pval_Age[pval_Age<0.05]
length(pval_Age_sig)
# 0

