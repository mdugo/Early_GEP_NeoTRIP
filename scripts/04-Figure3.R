
library(Biobase)
library(nnet)
library(ggplot2)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)
library(rcartocolor)
library(RColorBrewer)
library(showtext)
library(cowplot)
library(binom)
library(here)
source(file.path(here("scripts"), "00a-helper_functions.R"))

data_dir <- here("data", "processed")

load(file.path(data_dir, "NeoTRIP_singscore_BL_D1C2.RData"))

scoreSS <- scoreSS[, scoreSS$Timepoint == "Baseline"]
sum(is.na(scoreSS$Early.pCR))
scoreSS <- scoreSS[, !is.na(scoreSS$Early.pCR)]
scoreSS$Early.pCR <- gsub("e-pCR", "early-pCR", scoreSS$Early.pCR)
scoreSS$Early.pCR <- factor(scoreSS$Early.pCR, levels = c("Sustained early-pCR", "Transient early-pCR", "pCR", "RD"))

#**********************************************
# Panel A
#**********************************************

clinical <- pData(scoreSS)[, c(2, 4:7)]
clinical$Tumor_cells <- ifelse(clinical$Early.pCR %in% c("Sustained early-pCR", "Transient early-pCR"), "No TC in core biopsy", "TC in core biopsy")
clinical$Tumor_cells <- factor(clinical$Tumor_cells, levels = c("TC in core biopsy", "No TC in core biopsy"))
clinical$Arm <- factor(clinical$Arm, levels = c("CT", "CT/A"))
df <- data.frame(Variable = clinical$Tumor_cells, 
  pCR_cat = clinical$pCR, 
  pCR = clinical$pCR.num, 
  Arm = clinical$Arm)
df$Group <- paste(df$Arm, df$Variable)
df$Group <- factor(df$Group, levels = c("CT TC in core biopsy", "CT No TC in core biopsy", "CT/A TC in core biopsy", "CT/A No TC in core biopsy"))
t2x2 <- table(df$Group, df$pCR)
perc <- round(t2x2/rowSums(t2x2)*100, 1)

pp <- table(df$Variable, df$Arm)
chisq.test(pp)
ppp <- round(pp/colSums(pp)*100, 1)

df.perc <- data.frame(Perc = perc[, 2], 
   N = t2x2[, 2], 
   Variable = gsub("CT/A ", "", gsub("CT ", "", rownames(t2x2))), 
   Arm = gsub(" .*", "", rownames(t2x2)), 
   Group = rownames(t2x2))
df.perc$Arm <- factor(df.perc$Arm, levels = c("CT", "CT/A"))
df.perc$Group <- factor(df.perc$Group, levels = df.perc$Group)
df.perc$Position <- c(1, 2, 4, 5)
df.perc$Variable <- factor(df.perc$Variable, levels = c("TC in core biopsy", "No TC in core biopsy"))

fit1 <- glm(pCR ~ Variable, family = binomial, data = df[df$Arm == "CT", ])
fit2 <- glm(pCR ~ Variable, family = binomial, data = df[df$Arm == "CT/A", ])
fit.int <- glm(pCR ~ Variable*Arm, family = binomial, data = df)

pvals <- data.frame(Arm = c("CT", "CT/A"), 
   pval = c(round(summary(fit1)$coefficients[2, 4], 3), 
    round(summary(fit2)$coefficients[2, 4], 4)), 
   PosLab = c(1.5, 4.5))

g1 <- ggplot(data = df.perc, aes(x = Position, y = Perc, fill = Variable)) +
 geom_bar(stat = "identity", width = 0.99, position = position_dodge(width = 0.9)) +
 xlab("") +
 ylab("pCR rate (%)") +
 scale_x_continuous(breaks = c(1.5, 4.5), labels = c("CT\nN = 116", "CT/A\nN = 114")) +
 scale_y_continuous(breaks = seq(0, 100, 20), limits = c(0, 103)) +
 geom_text(aes(x = Position, y = Perc+5, label = paste0(round(Perc, 1), "")), color = "black", size = 5) +
 scale_fill_manual(name = "", values = carto_pal(7, "Tropic")[c(3, 5)]) +
 theme_pubr(border = F) +
 theme(axis.title = element_text(size = 12), 
 axis.text.x = element_text(size = 12), 
 axis.text.y = element_text(size = 12), 
 title = element_text(size = 18), 
 legend.position = "top") +
 guides(fill = guide_legend(ncol = 2)) +
 annotate('text', 1.5, 92, 
  label = paste("italic(P) == ", pvals$pval[1]), 
  parse = TRUE, 
  hjust = 0.5, size = 4.5) +
 annotate('text', 4.5, 92, 
  label = paste("italic(P) == ", pvals$pval[2]), 
  parse = TRUE, 
  hjust = 0.5, size = 4.5) +
 annotate('text', x = 3, y = 102.5, 
  label = paste("italic(P)[interaction] == ", round(summary(fit.int)$coefficients[4, 4], 3)), 
  parse = TRUE, 
  hjust = 0.5, size = 4.5) +
 annotate("rect", xmin = 1.5, xmax = 4.5, ymin = 98, ymax = 98, alpha = 1, colour = "black") +
 annotate("rect", xmin = 1.5, xmax = 1.5, ymin = 96, ymax = 98, alpha = 1, colour = "black") +
 annotate("rect", xmin = 4.5, xmax = 4.5, ymin = 96, ymax = 98, alpha = 1, colour = "black")

# PPV/NPV calculation for CT

t2x2 <- table(df$Variable[df$Arm == "CT"], df$pCR_cat[df$Arm == "CT"])
TP <- t2x2["No TC in core biopsy", "pCR"]
FP <- t2x2["No TC in core biopsy", "RD"]
FN <- t2x2["TC in core biopsy", "pCR"]
TN <- t2x2["TC in core biopsy", "RD"]
PPV <- TP / (TP + FP)
NPV <- TN / (TN + FN)
ppv_ci <- binom.confint(TP, TP + FP, methods = "exact")
npv_ci <- binom.confint(TN, TN + FN, methods = "exact")
cat("PPV:", round(PPV*100, 1), "%", 
 " (95% CI:", round(ppv_ci$lower*100, 1), "-", round(ppv_ci$upper*100, 1), "%)\n")
cat("NPV:", round(NPV*100, 1), "%", 
 " (95% CI:", round(npv_ci$lower*100, 1), "-", round(npv_ci$upper*100, 1), "%)\n")

# PPV/NPV calculation for CT/A

t2x2 <- table(df$Variable[df$Arm == "CT/A"], df$pCR_cat[df$Arm == "CT/A"])
TP <- t2x2["No TC in core biopsy", "pCR"]
FP <- t2x2["No TC in core biopsy", "RD"]
FN <- t2x2["TC in core biopsy", "pCR"]
TN <- t2x2["TC in core biopsy", "RD"]
PPV <- TP / (TP + FP)
NPV <- TN / (TN + FN)
ppv_ci <- binom.confint(TP, TP + FP, methods = "exact")
npv_ci <- binom.confint(TN, TN + FN, methods = "exact")
cat("PPV:", round(PPV*100, 1), "%", 
 " (95% CI:", round(ppv_ci$lower*100, 1), "-", round(ppv_ci$upper*100, 1), "%)\n")
cat("NPV:", round(NPV*100, 1), "%", 
 " (95% CI:", round(npv_ci$lower*100, 1), "-", round(npv_ci$upper*100, 1), "%)\n")

#**********************************************
# Panel B
#**********************************************

res <- data.frame(Variable = rownames(scoreSS),
         OR_CT = 0,
         P_CT = 0,
         OR_CTA = 0,
         P_CTA = 0,
         ORR_Interaction = 0,
         P_Interaction = 0)

for(i in 1:nrow(scoreSS)){
 df <- data.frame(pCR = scoreSS$Early.pCR, Arm = scoreSS$Arm, Variable = exprs(scoreSS)[i, ])

 # CT
 model <- multinom(pCR ~ Variable, data = df[df$Arm == "CT", ], Hess = TRUE)
 z <- summary(model)$coefficients/summary(model)$standard.errors
 p <- (1 - pnorm(abs(z), 0, 1)) * 2
 res$OR_CT[i] <- exp(summary(model)$coefficients[2, 2])
 res$P_CT[i] <- p[2, 2]
 
 # CT/A
 model <- multinom(pCR ~ Variable, data = df[df$Arm == "CT/A", ], Hess = TRUE)
 z <- summary(model)$coefficients/summary(model)$standard.errors
 p <- (1 - pnorm(abs(z), 0, 1)) * 2
 res$OR_CTA[i] <- exp(summary(model)$coefficients[2, 2])
 res$P_CTA[i] <- p[2, 2]
 
 # Interaction
 model <- multinom(pCR ~ Variable*Arm, data = df, Hess = TRUE)
 z <- summary(model)$coefficients/summary(model)$standard.errors
 p <- (1 - pnorm(abs(z), 0, 1)) * 2
 res$ORR_Interaction[i] <- exp(summary(model)$coefficients[2, 4])
 res$P_Interaction[i] <- p[2, 4]
}
res$FDR_CT <- p.adjust(res$P_CT, method = "BH")
res$FDR_CTA <- p.adjust(res$P_CTA, method = "BH")
res <- res[, c(1, 2, 3, 8, 4, 5, 9, 6, 7)]
write.table(res, file = file.path(here("results", "tables"), "02-Early_pCR_prediciton_baseline.txt"), sep = "\t", row.names=F, quote=F)

# heatmap

sig <- res$Variable[res$P_CTA < 0.05 | res$P_CT < 0.05]
scaledMat <- t(scale(t(exprs(scoreSS))))
medie <- t(apply(scaledMat, 1, function(x)tapply(x, paste(scoreSS$Arm, scoreSS$Early.pCR), mean)))
medie <- medie[sig, ]

intPvalNum <- res$P_Interaction
names(intPvalNum) <- res$Variable
intPvalNum <- intPvalNum[rownames(medie)]
intPvalNum <- round(intPvalNum, 2)

pMat <- res[match(rownames(medie), res$Variable), ]
rownames(pMat) <- pMat$Variable

annotCol <- structure(c(carto_pal(12, "Vivid")[1:9], "grey75"), names = names(table(fData(scoreSS)$Annotation)))
fdata <- fData(scoreSS)
fdata <- fdata[rownames(medie), ]
fdata$Annotation <- droplevels(fdata$Annotation)
la <- rowAnnotation(Annotation = fdata$Annotation, 
   col = list(Annotation = annotCol), 
   annotation_legend_param = list(Annotation = list(at = levels(fdata$Annotation)[c(5, 1, 7, 6, 2, 3, 4)], 
       labels_gp = gpar(fontsize = 11), title_gp = gpar(fontsize = 11, fontface = "bold"))))


ra <- rowAnnotation("P" = anno_text(round(pMat$P_CT, 2), show_name = TRUE, gp = gpar(fontsize = 12, fontface = ifelse(pMat$P_CT<0.05, "bold", "plain"))))

hm1 <- Heatmap(medie[, c(3,4,1,2)], 
  name = "Z-score", 
  heatmap_legend_param = list(title_gp = gpar(fontsize = 12, fontface = "bold"), 
     labels_gp = gpar(fontsize = 12)), 
  cluster_columns = F, 
  row_split = fdata$Annotation, 
  row_title = NULL, 
  cluster_rows = T, 
  clustering_method_rows = "ward.D2", 
  show_row_dend = F, 
  row_names_side = "left", 
  row_names_gp = gpar(fontsize = 12), 
  column_names_gp = gpar(fontsize = 10), 
  column_labels = levels(scoreSS$Early.pCR), 
  column_title = "CT", 
  column_title_gp = gpar(fontsize = 14), 
  rect_gp = gpar(col = "grey40", lwd = 0.6), 
  right_annotation = ra, 
  width = unit(0.7, "in"), 
  height = unit(2.9, "in"), 
  left_annotation = la)


ra2 <- rowAnnotation("P CT/A" = anno_text(paste0(round(pMat$P_CTA, 2), " "), show_name = TRUE, gp = gpar(fontsize = 12, fontface = ifelse(pMat$P_CTA<0.05, "bold", "plain"))), 
   "P interaction" = anno_text(intPvalNum, show_name = TRUE, gp = gpar(fontsize = 12, fontface = ifelse(intPvalNum<0.05, "bold", "plain"))), 
   annotation_label = list("P", expression(P[interaction])))

hm2 <- Heatmap(medie[, c(7,8,5,6)], 
  name = "Z-score2", 
  show_heatmap_legend = FALSE, 
  heatmap_legend_param = list(title_gp = gpar(fontsize = 12, fontface = "bold"), 
      labels_gp = gpar(fontsize = 12)), 
  cluster_columns = F, 
  row_title = NULL, 
  show_row_dend = F, 
  row_names_side = "left", 
  row_names_gp = gpar(fontsize = 12), 
  column_names_gp = gpar(fontsize = 10), 
  clustering_method_rows = "ward.D2", 
  column_labels = levels(scoreSS$Early.pCR), 
  column_title = "CT/A", 
  column_title_gp = gpar(fontsize = 14), 
  rect_gp = gpar(col = "grey40", lwd = 0.6), 
  right_annotation = ra2, 
  width = unit(0.7, "in"), 
  height = unit(2.9, "in"))


g2 = grid.grabExpr(draw(hm1+hm2, merge_legend = TRUE, legend_grouping = "original", annotation_legend_side = "right", heatmap_legend_side = "right"))

#**********************************************
# Panel C
#**********************************************

geneset = "cTME_Neutrophils"
df <- data.frame(Var1 = exprs(scoreSS)[geneset, ], 
  Timepoint = scoreSS$Timepoint, 
  Arm = scoreSS$Arm, 
  pCR = scoreSS$Early.pCR)
nsamp <- data.frame(table(df$pCR, df$Arm))
colnames(nsamp) <- c("pCR", "Arm", "N")
rng <- (max(df$Var1) - min(df$Var1))/20
nsamp$Ypos <- min(df$Var1) - rng
nsamp$Label <- paste0("N = ", nsamp$N)
signif <- data.frame(Arm = c("CT", "CT/A"), Ypos = max(df$Var1) + rng, Label = paste0("P = ", c(round(res$P_CT[res$Variable == geneset], 3), round(res$P_CTA[res$Variable == geneset], 3))))

g3 <- ggplot(data = df, aes(x = Arm, y = Var1, fill = pCR)) +
 geom_boxplot(coef = NULL, width = 0.5, position = position_dodge(width = 0.8)) +
 geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), size = 0.6) +
 scale_fill_manual(name = "", values = carto_pal(12, "Pastel")[c(4, 3, 2, 1)]) +
 xlab("") +
 ylab("Score") +
 labs(title = geneset) +
 theme_pubr(border = T, legend = "right") +
 theme(panel.grid.major.y = element_line(colour = "grey80", linewidth = 0.6, linetype = "dotted"), 
 axis.text = element_text(size = 12), 
 axis.title = element_text(size = 12), 
 plot.title = element_text(size = 12, face = "bold", hjust = 0.5), 
 legend.text = element_text(size = 13)) +
 geom_text(data = nsamp, aes(y = Ypos, label = N), position = position_dodge(width = 0.8), size = 4.5) +
 geom_text(data = signif, aes(x = Arm, y = Ypos, label = Label), size = 4.1, inherit.aes = F)


geneset = "HM_GLYCOLYSIS"
df <- data.frame(Var1 = exprs(scoreSS)[geneset, ], 
  Timepoint = scoreSS$Timepoint, 
  Arm = scoreSS$Arm, 
  pCR = scoreSS$Early.pCR)
nsamp <- data.frame(table(df$pCR, df$Arm))
colnames(nsamp) <- c("pCR", "Arm", "N")
rng <- (max(df$Var1) - min(df$Var1))/20
nsamp$Ypos <- min(df$Var1) - rng
nsamp$Label <- paste0("N = ", nsamp$N)
signif <- data.frame(Arm = c("CT", "CT/A"), Ypos = max(df$Var1) + rng, Label = paste0("P = ", c(round(res$P_CT[res$Variable == geneset], 3), round(res$P_CTA[res$Variable == geneset], 3))))

g4 <- ggplot(data = df, aes(x = Arm, y = Var1, fill = pCR)) +
 geom_boxplot(coef = NULL, width = 0.5, position = position_dodge(width = 0.8)) +
 geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), size = 0.6) +
 scale_fill_manual(name = "", values = carto_pal(12, "Pastel")[c(4, 3, 2, 1)]) +
 xlab("") +
 ylab("Score") +
 labs(title = geneset) +
 theme_pubr(border = T, legend = "right") +
 theme(panel.grid.major.y = element_line(colour = "grey80", linewidth = 0.6, linetype = "dotted"), 
 axis.text = element_text(size = 12), 
 axis.title = element_text(size = 12), 
 plot.title = element_text(size = 12, face = "bold", hjust = 0.5), 
 legend.text = element_text(size = 13)) +
 geom_text(data = nsamp, aes(y = Ypos, label = N), position = position_dodge(width = 0.8), size = 4.5) +
 geom_text(data = signif, aes(x = Arm, y = Ypos, label = Label), size = 4.1, inherit.aes = F)


empty_plot <- ggplot() +
 theme_void() +
 theme(
 plot.margin = margin(0, 0, 0, 0)
 )

#**********************************************
# Panel assembly
#**********************************************

g34 <- ggarrange(g3, g4, nrow = 1, ncol = 2, common.legend = T)
g12 <- plot_grid(g1, g2, nrow = 1, ncol = 2, labels = c("A", "B"), label_size = 16, rel_widths = c(0.42, 1))
g34e <- plot_grid(g34, empty_plot, nrow = 1, ncol = 2, labels = c("C", "D"), label_size = 16)

combined <- plot_grid(g12, g34e, nrow = 2, rel_heights = c(1, 0.8))

title <- ggdraw() +
  draw_label("Figure 3",
             x = 0, hjust = 0,
             fontface = "bold", size = 16)
final_plot <- plot_grid(title, combined,
                        ncol = 1,
                        rel_heights = c(0.1, 1)) 
ggsave(file.path(here("results", "figures"), "Figure3.eps"),
       plot = final_plot,
       device = cairo_ps,
       width = 13.6, height = 9.5)



