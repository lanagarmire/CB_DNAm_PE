######################################
#Analyze in-house DNA methylation data
######################################
# Include cell type deconvolution, differential analysis, cell proportion trend analysis
# Include code to generate figure 2b,2c,3b,3c,4a,4b,5a,5b


# set up -------------------------------------------------------------------------------------
# Library:
library(lumi)
library(limma)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(EpiDISH)
library(reshape2)
library(ggpubr)

# see data preprocessing script
load("/nfs/dcmb-lgarmire/liuwent/04-Full_Model/pd.RData")
load("/nfs/dcmb-lgarmire/xtyang/CordBlood/03-wenting_results/01-new_result_Garvin_ref/myNorm_new.RData")

#cell type deconvolution using Gervin's reference
load("/nfs/dcmb-lgarmire/xtyang/CordBlood/03-wenting_results/FlowSorted.CordBloodCombined.450k.compTable.rda")
ref = FlowSorted.CordBloodCombined.450k.compTable
out<-epidish(beta.m=myNorm0, ref.m=as.matrix(ref), method="CP")
count <- out$estF
pd = data.frame(cbind(out$estF, pd))
save(pd, file = "/nfs/dcmb-lgarmire/xtyang/CordBlood/03-wenting_results/01-new_result_Garvin_ref/pd_all2.RData")

#Figure 2:cell proportion boxplot by sample group-------------------------------------------------------

#figure 2b
cell_type = data.frame(cbind(pd[, 1:7], pd$Sample_Group))
colnames(cell_type)[8] = "Sample_Group"
cell_type_long = melt(cell_type, id.vars = "Sample_Group")

colnames(cell_type_long)[1] = "Sample_Group"
colnames(cell_type_long)[2] = "Cell_Type"
colnames(cell_type_long)[3] = "Cell_Type_Proportion"

# Define the order of the factor levels for the x-axis
cell_type_long$Cell_Type <- factor(cell_type_long$Cell_Type)
beta_tbl_unadjusted = data.frame(matrix(NA, nrow = 2, ncol=7))
for (i in 1:7) {
  fit <- lm(pd[, i] ~ Sample_Group , data = pd)
  beta_tbl_unadjusted[,i] = summary(fit)$coefficients[,1]
}
colnames(beta_tbl_unadjusted)=colnames(pd[,1:7])
rownames(beta_tbl_unadjusted) = rownames(summary(fit)$coefficient)
write.csv(beta_tbl_unadjusted, file = "/nfs/dcmb-lgarmire/xtyang/CordBlood/03-wenting_results/01-new_result_Garvin_ref/cell_prop_regression_beta_tbl_unadjusted.csv")

#figure 2b
tiff("/nfs/dcmb-lgarmire/xtyang/CordBlood/03-wenting_results/01-new_result_Garvin_ref/cell_prop_before_adj2.tiff", res = 200, width = 25, height = 20, units = "cm")
ggplot(cell_type_long, aes(Cell_Type, Cell_Type_Proportion, fill=Sample_Group)) + 
  geom_boxplot() +
  stat_compare_means(aes(label = paste0("p = ", after_stat(p.format)))) + 
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))+
  theme(text = element_text(size = 15))+
  xlab("Cell Type")+
  ylab("Cell Type Proportion")+
  labs(fill='Sample Group') 
dev.off()

#figure 2c
# number of rows = nrow(pd), number of cols = 7
residual_df <- matrix(NA, nrow = nrow(pd), ncol = 7)
pval= c()
pval_tbl = data.frame(matrix(NA, nrow = 9, ncol=7))
beta_tbl = data.frame(matrix(NA, nrow = 9, ncol=7))
for (i in 1:7) {
  fit <- lm(pd[, i] ~ Sample_Group+GA + BMI + Age + Parity + Eth2+Smoker , data = pd)
  pval_tbl[,i]=summary(fit)$coefficients[,4]
  beta_tbl[,i] = summary(fit)$coefficients[,1]
  pval[i] = summary(fit)$coefficient[2,4]
  residual_df[, i] <- fit$residuals+mean(pd[,i], na.rm=T)
  
}
colnames(pval_tbl)=colnames(pd[,1:7])
rownames(pval_tbl) = rownames(summary(fit)$coefficient)
colnames(beta_tbl)=colnames(pd[,1:7])
rownames(beta_tbl) = rownames(summary(fit)$coefficient)
write.csv(pval_tbl, file = "/nfs/dcmb-lgarmire/xtyang/CordBlood/03-wenting_results/01-new_result_Garvin_ref/cell_prop_regression_pval_tbl.csv")
write.csv(beta_tbl, file = "/nfs/dcmb-lgarmire/xtyang/CordBlood/03-wenting_results/01-new_result_Garvin_ref/cell_prop_regression_beta_tbl.csv")

# make it a data.frame and preserve column names
residual_df <- as.data.frame(residual_df)
colnames(residual_df) <- colnames(pd)[1:7]
residual_df$Sample_Group = pd$Sample_Group

residual_long = melt(residual_df, id.vars = "Sample_Group")

colnames(residual_long)[1] = "Sample_Group"
colnames(residual_long)[2] = "Cell_Type"
colnames(residual_long)[3] = "Cell_Type_Proportion"

# Define the order of the factor levels for the x-axis
residual_long$Cell_Type <- factor(residual_long$Cell_Type)

tiff("/nfs/dcmb-lgarmire/xtyang/CordBlood/03-wenting_results/01-new_result_Garvin_ref/cell_prop_after_adj2.tiff", res = 200, width = 25, height = 20, units = "cm")
ggplot(residual_long, aes(Cell_Type, Cell_Type_Proportion, fill=Sample_Group)) + 
  geom_boxplot() +
  annotate("text", label = paste0("p = ",round(pval, 3)), 
           size = 4, col = "black",  
           x=c(1.1,2.1,3.1,4.1,5.1,6.1,7.1), y = 0.8, vjust=1, hjust=1)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))+
  theme(text = element_text(size = 15))+
  xlab("Cell Type")+
  ylab("Cell Type Residual")+
  labs(fill='Sample Group') 
dev.off()

#Figure 3---------------------------------------------------------------
#Null model (figure 3b)
design = model.matrix(formula("~Sample_Group+Slide+Array"), data = pd)
Mvalues = beta2m(myNorm0)
fitNull = lmFit(Mvalues, design)
fitNull = eBayes(fitNull)
toptable_null = limma::topTable(fitNull, n=dim(fitNull)[1], adjust.method="BH", coef = 2)
sum(toptable_null$adj.P.Val<0.05)#9239
pval0 = p.adjust(fitNull$p.value[,2], "BH")

col_list2 = sapply(seq(1, dim(fitNull$coefficients)[1]), function(i) {
  if (pval0[i]<0.05) {return("red")}
  else {return("black")} })

png("/nfs/dcmb-lgarmire/xtyang/CordBlood/03-wenting_results/01-new_result_Garvin_ref/volcano_3b.png")
par(cex.axis = 1.5, cex.lab = 1.5)

volcanoplot(fitNull,coef=2,col=col_list2,highlight=0,
            names=rownames(fitNull$coefficients), xlim=c(-4,4))
# abline(h=-log10(0.05), lty=2, col="red")
abline(v=-1, lty=2, col="grey")
abline(v=1, lty=2, col="grey")
text(-1, 5, "x=-1", cex=1.2)
text(1, 5, "x=1", cex=1.2)
dev.off()


# all confounder model 
#transform compositional cell proportion data with ILR
library(compositions)
library(zCompositions)
cell_type = pd[,1:7]
cell_type[cell_type < 0] <- 0          # guard against negatives
#cell_type_zr <- cmultRepl(cell_type, method = "CZM", output = "prop")  # compositional zero replacement
cell_ilr <- ilr(acomp(cell_type))
colnames(cell_ilr) = c("ct1", "ct2", "ct3", "ct4", "ct5","ct6")

#save(pd, file = "pd_all.RData")
pd$Slide = as.character(pd$Slide)
allConfounders = pd%>%dplyr::select("Sample_Group", "Age", "Parity", "Eth2", "GA", "Smoker","Slide", "Array")%>%cbind(cell_ilr)

formula = paste0(names(allConfounders), collapse = ' + ')
formula = paste0("~", formula)
formula = formula(formula)
design2 = model.matrix(formula, data = allConfounders)
colnames(design2)[2] = "Cases"

#fit linear model
fitConfounders = lmFit(Mvalues, design2)
fitConfounders = eBayes(fitConfounders)
topTable_fitConfounders = limma::topTable(fitConfounders, coef = "Cases",n=dim(fitConfounders)[1], adjust.method="BH")
pval = p.adjust(fitConfounders$p.value[,2], "BH")

#calculate bacon inflation score
t_stats <- fitConfounders$t[,"Cases"] #random shuffled group
df<- fitConfounders$df.total[1]
z_score2 = zscoreT(t_stats, df, approx=FALSE, method = "bailey")
bc2 <- bacon(z_score2)
limma_inflation_score2 = inflation(bc2)#0.922

#make a toptable on second column sample_group
col_list2 = sapply(seq(1, dim(fitConfounders$coefficients)[1]), function(i) {
  if (pval[i]<0.05) {return("red")}
  else {return("black")} })
cpg <- data.frame(logFC=fitConfounders$coefficients[,2],pval=pval)
sig_hyper <- rownames(cpg[cpg$logFC>0 & cpg$pval<0.05,])
length(sig_hyper)#0
sig_hypo <- rownames(cpg[cpg$logFC<0 & cpg$pval<0.05,])
length(sig_hypo)

# figure 3c
png("/nfs/dcmb-lgarmire/xtyang/CordBlood/03-wenting_results/01-new_result_Garvin_ref/volcano_3c.png")
par(cex.axis = 1.5, cex.lab = 1.5)

volcanoplot(fitConfounders,coef=2,col=col_list2,highlight=0,
            names=rownames(fitConfounders$coefficients), xlim=c(-4,4))
# abline(h=-log10(0.05), lty=2, col="red")
abline(v=-1, lty=2, col="grey")
abline(v=1, lty=2, col="grey")
text(-1, 5, "x=-1", cex=1.2)
text(1, 5, "x=1", cex=1.2)
dev.off()

# number of significant cpg associate with each variable---------------------------------------------------

# toptable_ga = limma::topTable(fitConfounders, coef="GA", n=dim(fitConfounders)[1], adjust.method="BH")
# dim(toptable_ga[toptable_ga$adj.P.Val<0.05, ])
# toptable_age = limma::topTable(fitConfounders, coef="Age", n=dim(fitConfounders)[1], adjust.method="BH")
# dim(toptable_age[toptable_age$adj.P.Val<0.05, ])
# toptable_parity = limma::topTable(fitConfounders, coef="Parity", n=dim(fitConfounders)[1], adjust.method="BH")
# dim(toptable_parity[toptable_parity$adj.P.Val<0.05, ])
# toptable_Caucasian = limma::topTable(fitConfounders, coef="Eth2Caucasian", n=dim(fitConfounders)[1], adjust.method="BH")
# dim(toptable_Caucasian[toptable_Caucasian$adj.P.Val<0.05, ])
# toptable_PacificIslander = limma::topTable(fitConfounders, coef="Eth2PacificIslander", n=dim(fitConfounders)[1], adjust.method="BH")
# dim(toptable_PacificIslander[toptable_PacificIslander$adj.P.Val<0.05, ])
# toptable_CD8T = topTable(fitConfounders, coef="ct1", n=dim(fitConfounders)[1], adjust.method="BH")
# dim(toptable_CD8T[toptable_CD8T$adj.P.Val<0.05, ])
# toptable_CD4T = topTable(fitConfounders, coef="ct2", n=dim(fitConfounders)[1], adjust.method="BH")
# dim(toptable_CD8T[toptable_CD4T$adj.P.Val<0.05, ])
# toptable_NK = topTable(fitConfounders, coef="ct3", n=dim(fitConfounders)[1], adjust.method="BH")
# dim(toptable_NK[toptable_NK$adj.P.Val<0.05, ])
# toptable_Bcell = topTable(fitConfounders, coef="ct4", n=dim(fitConfounders)[1], adjust.method="BH")
# dim(toptable_Bcell[toptable_Bcell$adj.P.Val<0.05, ])
# toptable_Mono = topTable(fitConfounders, coef="ct5", n=dim(fitConfounders)[1], adjust.method="BH")
# dim(toptable_Mono[toptable_Mono$adj.P.Val<0.05, ])
# toptable_Gran = topTable(fitConfounders, coef="ct6", n=dim(fitConfounders)[1], adjust.method="BH")
# dim(toptable_Gran[toptable_Gran$adj.P.Val<0.05, ])

# #fit model with only Sample group and GA
# design3 = model.matrix(formula("~Sample_Group+GA+Slide+Array"), data = pd)
# colnames(design3)[2] = "Cases"
# #fit linear model
# Mvalues = beta2m(myNorm0)
# fitConfounders_ga_only = lmFit(Mvalues, design3)
# fitConfounders_ga_only = eBayes(fitConfounders_ga_only)
# toptable_pe = topTable(fitConfounders_ga_only, coef="Cases", n=dim(fitConfounders_ga_only)[1], adjust.method="BH")
# dim(toptable_pe[toptable_pe$adj.P.Val<0.05, ])
# toptable_ga = topTable(fitConfounders_ga_only, coef="GA", n=dim(fitConfounders_ga_only)[1], adjust.method="BH")
# dim(toptable_ga[toptable_ga$adj.P.Val<0.05, ])

# Figure 4---------------------------------------------------------------

# Figure 4a
library(dplyr)
library(ggplot2)
library(ggpubr)   # ggarrange
library(broom)    # tidy()

load("/nfs/dcmb-lgarmire/xtyang/CordBlood/03-wenting_results/01-new_result_Garvin_ref/pd_all2.RData")
data <- pd

# Choose the cell types you want to plot
cell_types <- c("Bcell", "CD4T", "CD8T", "Gran", "Mono", "NK", "nRBC")

# Helper: compute (beta, p) for GA -> cell_type and format a label
ga_label <- function(df, cell) {
  fml <- reformulate("GA", response = cell)     # e.g., Bcell ~ GA
  fit <- lm(fml, data = df)
  row_ga <- tidy(fit) |> filter(term == "GA")
  p <- row_ga$p.value[1]
  beta <- row_ga$estimate[1]
  # label: show beta and p (use NS if not sig)
  if (is.na(p)) {
    lbl <- sprintf("β = %.3f\np = NA", beta)
  } else {
    lbl <- sprintf("β = %.3f\np = %.2e", beta, p)
  }
  lbl
}

# Helper: make a single scatter + lm line + corner label
make_scatter <- function(df, cell) {
  lbl <- ga_label(df, cell)
  ggplot(df, aes(x = GA, y = .data[[cell]])) +
    geom_point(alpha = 0.85) +
    geom_smooth(method = "lm", se = TRUE) +
    scale_y_continuous(breaks = scales::breaks_extended(n = 5)) +
    annotate("text", label = lbl, size = 4, color = "black",
             x = Inf, y = Inf, vjust = 1, hjust = 1) +
    labs(x = "GA", y = cell) +
    theme_bw(base_size = 16) +
    theme(plot.margin = margin(5.5, 14, 5.5, 5.5))  # a bit more right margin for label
}

# Build all plots
plots <- lapply(cell_types, function(ct) make_scatter(data, ct))

# Save in one row (tweak width/height as needed)
png("/nfs/dcmb-lgarmire/xtyang/CordBlood/03-wenting_results/01-new_result_Garvin_ref/cell_prop_by_GA2.png", width = 3000, height = 600, res = 150)
ggarrange(plotlist = plots, nrow = 1, common.legend = FALSE)
dev.off()



#figure 4b------------------------------------------------------------
# Packages
library(dplyr)
library(ggplot2)
library(ggpubr)   # ggarrange
library(broom)    # tidy()
library(rlang)    # .data


load("/nfs/dcmb-lgarmire/xtyang/CordBlood/10-preterm_data/merged_pd2.RData")# this include idiopathic preterm controls
# --- Data (assumes 'merged_pd' already loaded) ---
# Pick your cell-type columns and "other" covariates once here:
cell_types   <- c("Bcell","CD4T","CD8T","Gran","Mono","NK","nRBC")
other_covars <- c("baby_sex", "dataset")  

# Helper: build formula "cell ~ PE + GA + PE:GA + other_covars"
make_formula <- function(cell, others) {
  reformulate(termlabels = c("Sample_Group", "GA", "Sample_Group:GA", others), response = cell)
}

# Helper: extract (beta, p) for the PE:GA interaction and format a label
interaction_label <- function(df, cell, others) {
  fml <- make_formula(cell, others)
  fit <- lm(fml, data = df)
  tt  <- tidy(fit)
  # Be robust to term name order ("PE:GA" vs "GA:PE")
  row <- tt %>% filter(term %in% c("Sample_GroupDisease:GA")) %>% dplyr::slice(1)

  
  if (nrow(row) == 0) {
    return("β(PE×GA) = NA\np = NA")
  } else {
    beta <- row$estimate[1]
    pval <- row$p.value[1]
    sprintf("β(PE×GA) = %.3f\np = %s",
            beta,
            ifelse(is.na(pval), "NA",
                   ifelse(pval < 0.001, format(pval, scientific = TRUE, digits = 2),
                          sprintf("%.3f", pval))))
  }
}

# Helper: one plot per cell type, colored by PE, with lm lines per PE and label
make_plot <- function(df, cell, others) {
  lbl <- interaction_label(df, cell, others)
  ggplot(df, aes(x = GA, y = .data[[cell]], color = Sample_Group)) +
    geom_point(alpha = 0.9, size = 1.2) +
    geom_smooth(aes(group = Sample_Group), method = "lm", se = FALSE) +
    scale_y_continuous(breaks = scales::breaks_extended(5)) +
    annotate("text", label = lbl, x = Inf, y = Inf, hjust = 1, vjust = 1,
             size = 4, color = "black") +
    labs(x = "GA (weeks)", y = cell, color = "Sample Group") +
    theme_bw(base_size = 16) +
    theme(legend.position = "top",
          plot.margin = margin(5.5, 14, 5.5, 5.5))
}

# Build all panels and save in a single-row figure with shared legend
plots <- lapply(cell_types, function(ct) make_plot(merged_pd, ct, other_covars))

png("/nfs/dcmb-lgarmire/xtyang/CordBlood/03-wenting_results/01-new_result_Garvin_ref/cell_prop_by_GA_PE_interaction2.png", width = 3000, height = 600, res = 150)
ggarrange(plotlist = plots, nrow = 1, common.legend = TRUE, legend = "top")
dev.off()


#Figure 5a:---------------------------------------------------------------

# Packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(broom)

# --- Load data ---
load("/nfs/dcmb-lgarmire/xtyang/CordBlood/03-wenting_results/01-new_result_Garvin_ref/pd_all.RData") # has pd with Sample_Group, GA
load("/nfs/dcmb-lgarmire/liuwent/04b-cell_type_deconvolution/pd_PBMC.RData")      # has pd_PBMC with Sample_Group, GA
#load("/nfs/dcmb-lgarmire/liuwent/04b-cell_type_deconvolution/estF_nogran.RData")  # WB (without gran)   -> estF_nogran
load("/nfs/dcmb-lgarmire/liuwent/04b-cell_type_deconvolution/ctp_PBMC_nogran.RData") # PBMC (without gran) -> ctp_PBMC_nogran

#recalculate cell proportion with no granulocyte
cell_mat <- pd[, c("Bcell","CD4T","CD8T","Mono","NK","nRBC")]
cell_row_sum <- rowSums(cell_mat)
estF_nogran <- sweep(cell_mat, 1, cell_row_sum, "/")

# --- Build tidy data ---
# whole blood
data_wb <- pd[, c("Sample_Group","GA","CD8T","CD4T","NK","Bcell","Mono","Gran","nRBC")]|>
  as.data.frame() |>
  mutate(Group = "WB")

# PBMC
data_pbmc <- cbind(pd_PBMC[, c("Sample_Group","GA")], ctp_PBMC_nogran) |>
  as.data.frame() |>
  mutate(Group = "PBMC")

# Combine
data_all <- bind_rows(data_wb, data_pbmc)

# Choose cell types present in the *nogran* sets
cell_types <- intersect(colnames(estF_nogran), colnames(ctp_PBMC_nogran))

# Long (tidy) format
df_long <- data_all %>%
  mutate(
    GA = as.numeric(GA),
    Sample_Group = factor(Sample_Group, levels = c("Controls","Disease")),
    Group = factor(Group, levels = c("WB","PBMC"))
  ) %>%
  pivot_longer(cols = all_of(cell_types), names_to = "Cell", values_to = "Prop") %>%
  mutate(Prop = suppressWarnings(as.numeric(Prop)))

# --- Filter to PE only ---
df_pe <- df_long %>% filter(Sample_Group == "Disease")

# --- Build per-facet labels: β/p for GA within each Group (WB, PBMC) ---
labels_df <- df_pe %>%
  group_by(Cell) %>%
  do({
    fit <- lm(Prop ~ GA*Group, data = .)
    row <- tidy(fit) %>% filter(term == "GA:GroupPBMC")
    beta <- row$estimate[1]
    p    <- row$p.value[1]
    tibble(lbl = sprintf(ifelse(is.na(p), "p=NA", sprintf("interaction p-value =%.2e", p))))
  }) %>%
  summarise(label = paste(lbl, collapse = "\n"), .groups = "drop") %>%
  mutate(x = Inf, y = Inf, hjust = 1, vjust = 1)

# --- One-panel multi-facet plot (PE only) ---
p <- ggplot(df_pe, aes(x = GA, y = Prop, color = Group)) +
  geom_point(alpha = 0.85, size = 1.8) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_color_manual(values = c("WB" = "#d95f02", "PBMC" = "#7570b3")) +
  scale_y_continuous(breaks = scales::breaks_extended(n = 5)) +
  labs(x = "GA (Disease)", y = "Estimated proportion", color = "Source") +
  facet_wrap(~ Cell, ncol = 6, scales = "free_y") +
  geom_text(data = labels_df,
            aes(x = x, y = y, label = label, hjust = hjust, vjust = vjust),
            inherit.aes = FALSE, size = 4) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "white"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(5.5, 14, 5.5, 5.5)
  )

# Save a single figure (PDF/PNG)
ggsave("/nfs/dcmb-lgarmire/xtyang/CordBlood/03-wenting_results/01-new_result_Garvin_ref/PE_WB_vs_PBMC_GA_celltypes_faceted2.png",
       plot = p, width = 3000/150, height = 700/150, units = "in", dpi = 150)

#Figure 5b:---------------------------------------------------------------
cell_types <- c("Bcell","CD4T","CD8T","Mono","NK","nRBC")
df_ctrl <- df_long %>% filter(Sample_Group == "Controls")

labels_df <- df_ctrl %>%
  group_by(Cell) %>%
  do({
    fit <- lm(Prop ~ GA*Group, data = .)
    row <- tidy(fit) %>% filter(term == "GA:GroupPBMC")
    beta <- row$estimate[1]
    p    <- row$p.value[1]
    tibble(lbl = sprintf(ifelse(is.na(p), "p=NA", sprintf("interaction p-value =%.2e", p))))
  }) %>%
  summarise(label = paste(lbl, collapse = "\n"), .groups = "drop") %>%
  mutate(x = Inf, y = Inf, hjust = 1, vjust = 1)

# ---- One multi-facet plot (Controls only) ----
p_ctrl <- ggplot(df_ctrl, aes(x = GA, y = Prop, color = Group)) +
  geom_point(alpha = 0.85, size = 1.8) +
  geom_smooth(method = "lm", se = TRUE) +
  scale_color_manual(values = c("WB" = "#d95f02", "PBMC" = "#7570b3")) +
  scale_y_continuous(breaks = scales::breaks_extended(n = 5)) +
  labs(x = "GA (Control)", y = "Estimated proportion", color = "Source") +
  facet_wrap(~ Cell, ncol = 6, scales = "free_y") +
  geom_text(data = labels_df,
            aes(x = x, y = y, label = label, hjust = hjust, vjust = vjust),
            inherit.aes = FALSE, size = 4) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "white"),
    panel.grid.minor = element_blank(),
    plot.margin = margin(5.5, 14, 5.5, 5.5)
  )

# ---- Save one figure (PDF/PNG) ----
ggsave("/nfs/dcmb-lgarmire/xtyang/CordBlood/03-wenting_results/01-new_result_Garvin_ref/CTRL_WB_vs_PBMC_GA_celltypes_faceted2.png",
       plot = p_ctrl, width = 3000/150, height = 700/150, units = "in", dpi = 150)
