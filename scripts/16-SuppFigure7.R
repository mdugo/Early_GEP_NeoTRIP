
library(Biobase)
library(rcartocolor)
library(Biobase)
library(ggplot2)
library(ggpubr)
library(cowplot)

#*******************************
# Dynamics plot
#*******************************

data_dir <- here::here("data", "processed")

load(file.path(data_dir, "NeoTRIP_singscore_BL_D1C2.RData"))
paired <- names(which(table(scoreSS$Patient_ID) == 2))
scoreSSp <- scoreSS[, scoreSS$Patient_ID%in%paired]
scoreSSp <- scoreSSp[, order(scoreSSp$Arm, scoreSSp$Patient_ID)]
scoreSSp <- scoreSSp[, scoreSSp$Arm == "CT/A"]


df <- data.frame(Var1 = exprs(scoreSSp)["Iron_utilization", ], 
        Var2 = exprs(scoreSSp)["cTME_Cytotoxic_cells", ], 
        Timepoint = scoreSSp$Timepoint, 
        Arm = scoreSSp$Arm, 
        pCR = factor(scoreSSp$pCR, levels = c("pCR", "RD")), 
        Patient = scoreSSp$Patient_ID)
df$Group <- paste(df$Arm, df$pCR, sep = " - ")
df$Group <- factor(df$Group, levels = unique(df$Group)[c(1, 2, 4, 3)])
tmp <- df[df$Timepoint == "D1C2", ]
tmp$Var2_cat <- ifelse(tmp$Var2 >= median(tmp$Var2), "cTME_Cytotoxic_cells\nHigh", "cTME_Cytotoxic_cells\nLow")
df <- merge(df, tmp[, c("Patient", "Var2_cat")], by = "Patient")
df$Var2_cat <- factor(df$Var2_cat, levels = c("cTME_Cytotoxic_cells\nLow", "cTME_Cytotoxic_cells\nHigh"))

g1 <- ggplot(data = df, aes(x = Timepoint, y = Var1)) +
 geom_boxplot(fill = "grey80", coef = NULL, width = 0.3) +
 geom_line(aes(group = Patient, color = pCR)) +
 facet_wrap(~Group, nrow = 1) +
 scale_color_manual(name = "", values = carto_pal(12, "Pastel")[c(2, 1)]) +
 xlab("") +
 ylab("Score") +
 labs(title = "Iron_utilization") +
 theme_pubr(border = T, legend = "none") +
 theme(panel.grid.major.y = element_line(colour = "grey80", linewidth = 0.6, linetype = "dotted"), 
    strip.text = element_text(size = 12, face = "bold"), 
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5), 
    axis.text.x = element_text(angle = 90, hjust = 1))

g2 <- ggplot(data = df, aes(x = Timepoint, y = Var2)) +
 geom_boxplot(fill = "grey80", coef = NULL, width = 0.3) +
 geom_line(aes(group = Patient, color = pCR)) +
 facet_wrap(~Group, nrow = 1) +
 scale_color_manual(name = "", values = carto_pal(12, "Pastel")[c(2, 1)]) +
 xlab("") +
 ylab("Score") +
 labs(title = "cTME_Cytotoxic_cells") +
 theme_pubr(border = T, legend = "none") +
 theme(panel.grid.major.y = element_line(colour = "grey80", linewidth = 0.6, linetype = "dotted"), 
    strip.text = element_text(size = 12, face = "bold"), 
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5), 
    axis.text.x = element_text(angle = 90, hjust = 1))

g3 <- ggplot(data = df, aes(x = Timepoint, y = Var1)) +
 geom_boxplot(fill = "grey80", coef = NULL, width = 0.3) +
 geom_line(aes(group = Patient, color = pCR)) +
 facet_wrap(Var2_cat~Group, nrow = 1) +
 scale_color_manual(name = "", values = carto_pal(12, "Pastel")[c(2, 1)]) +
 xlab("") +
 ylab("Iron_utilization score") +
 theme_pubr(border = T, legend = "none") +
 theme(panel.grid.major.y = element_line(colour = "grey80", linewidth = 0.6, linetype = "dotted"), 
    strip.text = element_text(size = 12, face = "bold"), 
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5), 
    axis.text.x = element_text(angle = 90, hjust = 1))

g12 <- plot_grid(g1, g2, ncol = 2, nrow = 1)

#**************************
# Panel assembly
#**************************

combined <- plot_grid(g12, g3, ncol = 1, nrow = 2, labels = c("A", "B"))
title <- ggdraw() +
  draw_label("Supplementary Figure 7",
             x = 0, hjust = 0,
             fontface = "bold", size = 14)
final_plot <- plot_grid(title, combined,
                        ncol = 1,
                        rel_heights = c(0.1, 1))  # spazio sopra
ggsave(file.path(here::here("results", "figures"), "SuppFigure7.eps"),
       plot = final_plot,
       device = cairo_ps,
       width = 12, height = 9)

