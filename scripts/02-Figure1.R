
library(ggplot2)
library(ggpubr)
library(ggforce)
library(dplyr)
library(Biobase)
library(variancePartition)
library(rcartocolor)
library(ComplexHeatmap)
library(cowplot)
library(here)

data_dir <- here("data", "processed")

load(file.path(data_dir, "NeoTRIP_singscore_BL_D1C2.RData"))

#*******************************
# Panel C
#*******************************

ct <- scoreSS[, scoreSS$Arm == "CT"]
cta <- scoreSS[, scoreSS$Arm == "CT/A"]

# PCA in CT
pca.ct <- prcomp(t(exprs(ct)))
df.ct <- data.frame(Arm = "CT", 
          pData(ct)[, c("pCR", "Timepoint")], 
          pca.ct$x[, 1:2])
ve.ct <- summary(pca.ct)

# PCA in CT/A
pca.cta <- prcomp(t(exprs(cta)))
df.cta <- data.frame(Arm = "CT/A", 
           pData(cta)[, c("pCR", "Timepoint")], 
           pca.cta$x[, 1:2])
ve.cta <- summary(pca.cta)

# plots
g1 <- ggplot(data = df.ct, aes(x = PC1, y = PC2, shape = pCR)) +
 geom_point(aes(color = Timepoint), size = 3, alpha = 0.6) +
 scale_shape_manual(name = "pCR status", values = c(15, 17)) +
 theme_pubr(legend = "right", border = T) +
 scale_color_manual(name = "Timepoint", values = carto_pal(12, "Bold")[2:1]) +
 xlab(paste0("PC1 (", round(ve.ct$importance[2, 1]*100, 1), "% explained variance)")) +
 ylab(paste0("PC2 (", round(ve.ct$importance[2, 2]*100, 1), "% explained variance)")) +
 labs(title = "CT") +
 theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))

g2 <- ggplot(data = df.cta, aes(x = PC1, y = PC2, shape = pCR)) +
 geom_point(aes(color = Timepoint), size = 3, alpha = 0.6) +
 scale_shape_manual(name = "pCR status", values = c(15, 17)) +
 theme_pubr(legend = "right", border = T) +
 scale_color_manual(name = "Timepoint", values = carto_pal(12, "Bold")[2:1]) +
 xlab(paste0("PC1 (", round(ve.cta$importance[2, 1]*100, 1), "% explained variance)")) +
 ylab(paste0("PC2 (", round(ve.cta$importance[2, 2]*100, 1), "% explained variance)")) +
 labs(title = "CT/A") +
 theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))

g12 <- ggarrange(g1, g2, ncol = 2, nrow = 1, common.legend = T)

#*******************************
# Panel D
#*******************************

meta <- pData(ct)[, c("pCR", "Timepoint", "Patient_ID")]
meta$Timepoint <- factor(meta$Timepoint)
meta$pCR <- factor(meta$pCR)
meta$Patient_ID <- factor(meta$Patient_ID)
form <- ~ (1|pCR) + (1|Timepoint) + (1|Patient_ID)
varPart.ct <- fitExtractVarPartModel(exprs(ct), form, meta)
varPart.ct <- reshape2::melt(varPart.ct)
varPart.ct$Arm <- "CT"


meta <- pData(cta)[, c("pCR", "Timepoint", "Patient_ID")]
meta$Timepoint <- factor(meta$Timepoint)
meta$pCR <- factor(meta$pCR)
meta$Patient_ID <- factor(meta$Patient_ID)
form <- ~ (1|pCR) + (1|Timepoint) + (1|Patient_ID)
varPart.cta <- fitExtractVarPartModel(exprs(cta), form, meta)
varPart.cta <- reshape2::melt(varPart.cta)
varPart.cta$Arm <- "CT/A"


varPart <- rbind(varPart.ct, varPart.cta)

g3 <- ggplot(data = varPart, aes(x = variable, y = value, fill = Arm)) +
 geom_boxplot(coef = NULL) +
 theme_pubr(legend = "top") +
 xlab("") +
 ylab("Proportion variance explained") +
 scale_x_discrete(labels = c("Patient", "pCR status", "Timepoint", "Residuals")) +
 scale_fill_manual(name = "", values = carto_pal(12, "Safe")[1:2]) +
 theme(axis.text.x = element_text(angle = 45, hjust = 1), 
    panel.grid.major.y = element_line(linetype = "dotted", color = "grey80"))

g123 <- cowplot::plot_grid(g12, g3, ncol = 2, nrow = 1, rel_widths = c(1, 0.4), labels = c("C", "D"))

#*******************************
# Panel B
#*******************************

info <- unique(pData(scoreSS)[, c("Patient_ID", "Arm", "pCR")])
info$Baseline <- ifelse(info$Patient_ID %in% scoreSS$Patient_ID[scoreSS$Timepoint == "Baseline"], 1, 0)
info$D1C2 <- ifelse(info$Patient_ID %in% scoreSS$Patient_ID[scoreSS$Timepoint=="D1C2"], 1, 0)
info <- info[order(info$Arm, info$pCR, -info$Baseline, -info$D1C2),]
info$Patient_ID <- factor(info$Patient_ID, levels = info$Patient_ID)
info$Count <- "Paired"
info$Count[info$Baseline == 1 & info$D1C2 == 0] <- "Baseline"
info$Count[info$Baseline == 0 & info$D1C2 == 1] <- "D1C2"
info$Count <- factor(info$Count, levels = c("Paired", "Baseline", "D1C2"))


t2x2 <- data.frame(table(info$Count, paste(info$Arm, info$pCR)))
start <- c(1, cumsum(t2x2$Freq)[-length(cumsum(t2x2$Freq))] + 1)
end <- cumsum(t2x2$Freq)
align_to <- vector("list", length = 4*3)
names(align_to) <- LETTERS[1:length(align_to)]
for(i in 1:length(align_to)){
  print(i)
  align_to[[i]] <- seq(start[i], end[i], 1)
}
panel_fun = function(index, nm) {
  grid.text(length(index), 0.5, 0.5, gp = gpar(fontsize = 10))
}

ba<-HeatmapAnnotation(" " = info$Count,
                      foo = anno_block(
                        align_to = align_to,
                        panel_fun = panel_fun),
                      col = list(" " = structure(carto_pal(12, "Vivid")[1:3] , names = c("Paired", "Baseline", "D1C2"))),
                      show_annotation_name = F,
                      annotation_legend_param = list(" " = list(direction = "horizontal")))
ra<-rowAnnotation("No. of samples "= anno_numeric(colSums(info[, 4:5]),
                                                  rg = c(0, max(colSums(info[, 4:5]))),
                                                  bg_gp = gpar(fill = "#9EB9F3", col = "#9EB9F3")))

hm<-Heatmap(t(info[, 4:5]),
            col=c("grey90", "grey20"),
            show_heatmap_legend = T,
            name = "  ",
            heatmap_legend_param = list(at = c("1", "0"), labels = c("RNA-Seq yes", "RNA-Seq no")),
            cluster_rows = F,
            cluster_columns = F,
            rect_gp = gpar(col = "white", lwd = 0.2),
            show_column_names = F,
            row_names_side = "left",
            bottom_annotation = ba,
            right_annotation = ra,
            column_split = paste(info$Arm,info$pCR),
            height = unit(0.4,units = "in"))

g0 = grid.grabExpr(draw(hm,
                        annotation_legend_side="bottom",
                        heatmap_legend_side="bottom",
                        merge_legends=T))

combined <- plot_grid(g0, g123, ncol = 1, nrow = 2, rel_heights = c(0.5, 1), labels = c("B", ""))

title <- ggdraw() +
  draw_label("Figure 1",
             x = 0, hjust = 0,
             fontface = "bold", size = 14)
final_plot <- plot_grid(title, combined,
                        ncol = 1,
                        rel_heights = c(0.1, 1))
ggsave(file.path(here::here("results", "figures"), "Figure1.eps"),
       plot = final_plot,
       device = cairo_ps,
       width = 11, height = 7)

