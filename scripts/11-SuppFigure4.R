
library(ggplot2)
library(ggpubr)
library(rcartocolor)
library(cowplot)
library(Biobase)
source(file.path(here::here("scripts"), "02-helper_functions.R"))

data_dir <- here::here("data", "processed")

load(file.path(data_dir, "NeoTRIP_singscore_BL_D1C2.RData"))
paired <- names(which(table(scoreSS$Patient_ID) == 2))
scoreSS <- scoreSS[, scoreSS$Patient_ID %in% paired]
scoreSS <- scoreSS[, order(scoreSS$Arm, scoreSS$Patient_ID)]

#*****************************************
# sankey plot of TNBCtypes and timepoints
#*****************************************

df <- pData(scoreSS)[, c("Patient_ID", "Arm", "Timepoint", "TNBCtype", "pCR", "pCR.num")]
df <- reshape(df, idvar = "Patient_ID", timevar = "Timepoint", direction = "wide")
df <- df[, c(1, 2, 4, 5, 3, 7)]
colnames(df)[2:6] <- c("Arm", "pCR", "pCR.num", "Baseline", "D1C2")
df$Combo <- paste(df$Baseline, df$D1C2, sep = "|")
df$Baseline <- factor(df$Baseline, levels = rev(sort(unique(df$Baseline))))
df$D1C2 <- factor(df$D1C2, levels = rev(sort(unique(df$D1C2))))

g3 <- sankey_plot(df = df[df$Arm == "CT", 5:6], plot_title = "CT", fill_colors = c("grey60", carto_pal(12, "Vivid")[c(1:6)]), factor.levels = levels(df$Baseline)) +
 coord_cartesian(xlim = c(0.91, 2.1), clip = 'off') +
 theme(plot.margin = margin(t = 0.5, r = 0, b = 0.5, l = 0, unit = "cm")) 
g4 <- sankey_plot(df = df[df$Arm == "CT/A", 5:6], plot_title = "CT/A", fill_colors = c("grey60", carto_pal(12, "Vivid")[c(1:6)]), factor.levels = levels(df$Baseline)) +
 coord_cartesian(xlim = c(0.91, 2.1), clip = 'off') +
 theme(plot.margin = margin(t = 0.5, r = 0, b = 0.5, l = 0, unit = "cm"))
g34 <- ggarrange(g3, g4, ncol = 2, nrow = 1)

df <- df[-grep("UNS", df$D1C2), ]
df$D1C2 <- droplevels(df$D1C2)
t2x2 <- table(df$D1C2[df$Arm == "CT"], df$pCR[df$Arm == "CT"])
ntot <- rowSums(t2x2)
perc <- data.frame(round(prop.table(t2x2, margin = 1)*100, 1))
t2x2 <- data.frame(t2x2)
t2x2.2 <- table(df$D1C2[df$Arm == "CT/A"], df$pCR[df$Arm == "CT/A"])
ntot2 <- rowSums(t2x2.2)
perc.2 <- data.frame(round(prop.table(t2x2.2, margin = 1)*100, 1))
t2x2.2 <- data.frame(t2x2.2)
df.perc <- data.frame(rbind(perc, perc.2), N = c(t2x2$Freq, t2x2.2$Freq), Tot = c(ntot, ntot, ntot2, ntot2), Arm = rep(c("CT", "CT/A"), each = nrow(perc)))
df.perc$Var2 <- factor(df.perc$Var2, levels = c("RD", "pCR"))
df.perc$Label <- gsub("NaN", "0", paste0(df.perc$Freq, "%\n(N = ", df.perc$N, ")"))
df.perc$Label2 <- gsub("NaN", "0", paste0(df.perc$Freq, "% ", df.perc$N, "/", df.perc$Tot))
df.perc <- df.perc[df.perc$Var2 == "pCR", ]
df.perc$Freq[is.na(df.perc$Freq)] <- 0
g6 <- ggplot(data = df.perc, aes(x = Var1, y = Freq, fill = Arm)) +
 geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
 geom_text(aes(x = Var1, y = Freq+1, label = Label2), color = "black", size = 3.5, position = position_dodge(width = 0.9), angle = 90, hjust = 0) +
 theme_pubr() +
 xlab("") +
 ylab("Frequency (%)") +
 scale_y_continuous(limits = c(0, 80), breaks = seq(0, 80, 20)) +
 scale_fill_manual(name = "", values = carto_pal(12, "Safe")[c(1, 2)]) +
 theme(panel.grid.major.y = element_line(colour = "grey80", linewidth = 0.6, linetype = "dotted"))

#**************************
# Panel assembly
#**************************

combined <- plot_grid(g34, g6, ncol = 2, rel_widths = c(1, 0.5), labels = c("A", "B"), label_size = 16)
title <- ggdraw() +
  draw_label("Supplementary Figure 4",
             x = 0, hjust = 0,
             fontface = "bold", size = 16)
final_plot <- plot_grid(title, combined,
                        ncol = 1,
                        rel_heights = c(0.1, 1)) 
ggsave(file.path(here("results", "figures"), "SuppFigure4.eps"),
       plot = final_plot,
       device = cairo_ps,
       width = 13, height = 5.5)


