# Library:
library(dplyr)
library(limma)
library(ggplot2)
library(tidyverse)
library(gplots)

# Data from 1 - Raw data pre-processing.Rmd
load("/home/liuwent/04-Full_Model/pd.RData")
load("/home/liuwent/04b-cell_type_deconvolution/estF.RData")

estF <- as.data.frame(estF)
rownames(estF) <- paste(pd$Sample_Name, pd$Sample_Group, sep = "-")

# Heatmap.2:
pdf("/home/liuwent/04b-cell_type_deconvolution/Fig_3A_heatmap.pdf")

# Set the margins to increase space for row names and dendrogram
par(mar=c(10, 10, 2, 2) + 0.1, cex.main=1, cex.lab=0.7, cex.axis=0.7)

# create a vector indicating whether each row corresponds to a case or control sample
status_vec <- ifelse(grepl("Disease", rownames(estF)), "Disease", "Control")

# create a named vector of colors for each status value
status_colors <- c("Disease" = "red", "Control" = "blue")

col_fun <- colorRampPalette(c("blue", "green", "yellow", "red"))

heatmap_plot <- heatmap.2(as.matrix(estF), 
                          Rowv=TRUE, 
                          Colv=TRUE,
                          margins=c(10,5),
                          cexRow=0.5,
                          cexCol=1.5,
                          RowSideColors = status_colors[status_vec],
                          trace = "none", 
                          dendrogram = "both", 
                          col = col_fun(100),
                          scale = "none",
                          # add expression to add a legend
                          add.expr = {
                            # set the fill colors of the legend based on whether the row name contains "Disease" or not
                            legend("topleft",
                                   legend = c("Disease", "Control"),
                                   fill = status_colors[status_vec])
                          },
                          denscol='black',
                          lmat=rbind(c(0,0,4),c(3,1,2),c(0,0,5)),
                          lwid=c(0.1, 0.02, 0.5),
                          lhei=c(0.1, 0.55, 0.15)
                          )

# Define the color key properties
custom_colors <- col_fun(100)
custom_breaks <- c(-Inf, 0, 0.05, 0.1, 0.3, 0.5, Inf)
custom_values <- seq(0, 1, length.out = length(custom_breaks))

heatmap_scale <- scale_color_gradientn(colors = custom_colors, 
                                       values = custom_values, 
                                       breaks = custom_breaks)

heatmap_plot + 
  heatmap_scale

dev.off()

