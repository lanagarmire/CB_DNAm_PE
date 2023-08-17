# Library:
library(ggplot2)

# Data:
Gene_df <- data.frame(
  Adjustment = c('No adjustment', 'adjust for CVs', 'adjust for CVs and CTs'),
  Gene = c(4767, 2, 0)
)

Path_df <- data.frame(
  Adjustment = c('No adjustment', 'adjust for CVs', 'adjust for CVs and CTs'),
  Pathway = c(200, 0, 0)
)

# Barplot:
pdf("/home/liuwent/04d-DM_analysis_after_adj/Gene_bar.pdf", width = 10, height = 10)
ggplot(Gene_df, aes(x = factor(Adjustment, levels = c('No adjustment', 'adjust for CVs', 'adjust for CVs and CTs')), y = Gene)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = Gene), vjust = -0.5, size = 8) +
  labs(title = "Gene by Adjustment", x = "Adjustment", y = "Gene") +
  theme(axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), 
        plot.title = element_text(size = 22), axis.title = element_text(size = 22))
dev.off()

pdf("/home/liuwent/04d-DM_analysis_after_adj/Path_bar.pdf", width = 10, height = 10)
ggplot(Path_df, aes(x = factor(Adjustment, levels = c('No adjustment', 'adjust for CVs', 'adjust for CVs and CTs')), y = Pathway)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = Pathway), vjust = -0.5, size = 8) +
  labs(title = "Pathway by Adjustment", x = "Adjustment", y = "Pathway") +
  theme(axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15),
        plot.title = element_text(size = 22), axis.title = element_text(size = 22))
dev.off()