# Library:
library("IlluminaHumanMethylationEPICanno.ilm10b4.hg19")

# Data:
annEPIC = getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# Shrink the data info to Chr6:
annEPIC_chr6 <- annEPIC[which(annEPIC$chr == "chr6"), ]

# split the gene names
for(i in 1:length(annEPIC_chr6$UCSC_RefGene_Name)){annEPIC_chr6$UCSC_RefGene_Name_new[i] = strsplit(annEPIC_chr6$UCSC_RefGene_Name[i], ";")[[1]][1]}
annEPIC_chr6$UCSC_RefGene_Name_new

# split the gene features
for(i in 1:length(annEPIC_chr6$UCSC_RefGene_Group)){annEPIC_chr6$UCSC_RefGene_Group_new[i] = strsplit(annEPIC_chr6$UCSC_RefGene_Group[i], ";")[[1]][1]}
annEPIC_chr6$UCSC_RefGene_Group_new

# Eliminate NAs in UCSC_RefGene_Name_new
annEPIC_chr6_sub <- annEPIC_chr6[which(is.na(annEPIC_chr6$UCSC_RefGene_Name_new) == FALSE), ]

# Shrink the data info correspond to pos in the range of chr6:27486812-27487069 based on the "Islands_Name":
annEPIC_chr6_sub1 <- annEPIC_chr6[annEPIC_chr6$pos >= 27486811 & annEPIC_chr6$pos <= 27487070, ]

# Shrink the data info correspond to pos in the range of 6:27594908-27595149 based on the "HMM_Island":
annEPIC_chr6_sub2 <- annEPIC_chr6[annEPIC_chr6$pos >= 27594907 & annEPIC_chr6$pos <= 27595149, ]

# Shrink the data info correspond to pos in the range of 6:27486079-27487598 based on the "Regulatory_Feature_Name":
annEPIC_chr6_sub3 <- annEPIC_chr6[annEPIC_chr6$pos >= 27486078 & annEPIC_chr6$pos <= 27487599, ]

# Shrink the data info correspond to pos in the range of chr6:27486625-27487515 based on the "DNase_Hypersensitivity_NAME":
annEPIC_chr6_sub4 <- annEPIC_chr6[annEPIC_chr6$pos >= 27486624 & annEPIC_chr6$pos <= 27487516, ]

