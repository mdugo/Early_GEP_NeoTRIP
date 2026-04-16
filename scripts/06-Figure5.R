
library(Biobase)
library(ComplexHeatmap)
library(circlize)
library(rcartocolor)
library(RColorBrewer)
library(Biobase)
library(limma)
library(ggplot2)
library(ggpubr)
library(showtext)
library(cowplot)
library(dplyr)
source(file.path(here("scripts"), "00a-helper_functions.R"))

#**********************************************
#* Analysis of delta with pCR status
#**********************************************

data_dir <- here::here("data", "processed")

load(file.path(data_dir, "NeoTRIP_singscore_BL_D1C2.RData"))
paired <- names(which(table(scoreSS$Patient_ID) == 2))
scoreSSp <- scoreSS[, scoreSS$Patient_ID %in% paired]
scoreSSp <- scoreSSp[, order(scoreSSp$Arm, scoreSSp$Patient_ID, scoreSSp$Timepoint)]

deltaMat <- matrix(0, nrow = nrow(scoreSSp), ncol = length(unique(scoreSSp$Patient_ID)))
rownames(deltaMat) <- rownames(scoreSSp)
colnames(deltaMat) <- unique(scoreSSp$Patient_ID)

for(i in 1:ncol(deltaMat)){
 deltaMat[, i] <- exprs(scoreSSp)[, scoreSSp$Patient_ID == colnames(deltaMat)[i] & scoreSSp$Timepoint == "D1C2"] - exprs(scoreSSp)[, scoreSSp$Patient_ID == colnames(deltaMat)[i] & scoreSSp$Timepoint == "Baseline"]
}
pdata <- unique(pData(scoreSSp)[, c("Patient_ID", "Arm", "pCR", "pCR.num")])
rownames(pdata) <- pdata$Patient_ID
fdata <- fData(scoreSSp)
deltas <- ExpressionSet(assayData = deltaMat, 
           phenoData = new("AnnotatedDataFrame", pdata), 
           featureData = new("AnnotatedDataFrame", fdata))


res.ct.d <- logistic_regression(response = "pCR.num", 
               data = deltas[, deltas$Arm == "CT"], 
               sort = F)
res.cta.d <- logistic_regression(response = "pCR.num", 
                data = deltas[, deltas$Arm == "CT/A"], 
                sort = F)
res.int.d <- interaction_logistic_regression(response = "pCR.num", 
                      covariate_interaction = "Arm", 
                      data = deltas, 
                      sort = F)

fullRes <- merge(res.ct.d, res.cta.d, by = "Variable")
colnames(fullRes) <- gsub("_Variable.y", "_CTA", gsub("_Variable.x", "_CT", colnames(fullRes)))
fullRes <- merge(fullRes, res.int.d[, c("Variable", "P_Interaction", "FDR_Interaction")], by = "Variable")
write.table(fullRes, file = file.path(here::here("results", "tables"), "06-pCR_prediciton_delta.txt"), sep = "\t", row.names = F, quote = F)

#**********************************************
#* Panel A
#**********************************************

res.ct.d$Arm <- "CT"
res.cta.d$Arm <- "CT/A"
sig <- unique(c(res.ct.d$Variable[res.ct.d$FDR_Variable<0.1], res.cta.d$Variable[res.cta.d$FDR_Variable<0.1]))

annotCol <- structure(c(carto_pal(12, "Vivid")[1:9], "grey75"), names = names(table(fData(scoreSS)$Annotation)))

orMat <- cbind(res.ct.d$OR_Variable, 
       res.cta.d$OR_Variable)
rownames(orMat) <- res.ct.d$Variable
colnames(orMat) <- c("BL_CT", "BL_CTA")

fdrMat <- cbind(res.ct.d$FDR_Variable, 
       res.cta.d$FDR_Variable)
rownames(fdrMat) <- res.ct.d$Variable
colnames(fdrMat) <- c("BL_CT", "BL_CTA")

orMat <- orMat[sig, ]
fdrMat <- fdrMat[sig, ]
fdrMat2 <- ifelse(fdrMat<0.05, "\U2217\U2217", ifelse(fdrMat<0.1, "\U2217", ""))
logOrMat <- log10(orMat)


intPval <- ifelse(res.int.d$P_Interaction<0.01, "p < 0.01", ifelse(res.int.d$P_Interaction<0.05, "p < 0.05", "n.s."))
intPval <- factor(intPval, levels = c("p < 0.01", "p < 0.05", "n.s."))
names(intPval) <- res.int.d$Variable
intPval <- intPval[rownames(logOrMat)]
intPvalNum <- res.int.d$P_Interaction
names(intPvalNum) <- res.int.d$Variable
intPvalNum <- intPvalNum[rownames(logOrMat)]
intPvalNum <- round(intPvalNum, 3)

ra <- rowAnnotation("Interaction (Arm)" = intPval, 
         col = list("Interaction (Arm)" = structure(c(carto_pal(7, "SunsetDark")[c(6, 3)], "grey80"), names = levels(intPval))), 
         annotation_legend_param = list(title = "Interaction test"))

textra <- rowAnnotation("P interaction" = anno_text(intPvalNum, show_name = TRUE, gp = gpar(fontsize = 12, fontface = ifelse(intPvalNum<0.05, "bold", "plain"))), 
           annotation_label = expression(P[interaction]))

fdata <- fData(deltas)
fdata <- fdata[rownames(logOrMat), ]
fdata$Annotation <- droplevels(fdata$Annotation)
la <- rowAnnotation(Annotation = fdata$Annotation, 
         col = list(Annotation = annotCol), 
         annotation_legend_param = list(Annotation = list(at = levels(fdata$Annotation)[c(1, 2, 5, 7, 3, 4, 6)])))


orCol <- colorRamp2(c(min(logOrMat), 0, max(logOrMat)), brewer.pal(11, "RdBu")[c(10, 6, 2)])
showtext.auto()
hm.d <- Heatmap(logOrMat, 
       name = "OR", 
       col = orCol, 
       row_split = fdata$Annotation, 
       row_title = NULL, 
       heatmap_legend_param = list(title = expression("log"[10](OR)), 
                     title_gp = gpar(fontsize = 10, fontface = "bold")), 
       cluster_columns = F, 
       cluster_rows = T, 
       clustering_method_rows = "ward.D2", 
       show_row_dend = F, 
       row_names_side = "left", 
       row_names_gp = gpar(fontsize = 13), 
       column_names_gp = gpar(fontsize = 13), 
       column_labels = c("CT", "CT/A"), 
       rect_gp = gpar(col = "grey40", lwd = 0.6), 
       right_annotation = textra, 
       left_annotation = la, 
       width = unit(0.7, "in"), 
       cell_fun = function(j, i, x, y, width, height, fill) {
        grid.text(fdrMat2[i, j], x, y, gp = gpar(fontsize = 11, fontface = "bold"))}
)
lgd_list = list(
 Legend(labels = c("*: FDR < 0.1", "**: FDR < 0.05"), title = "Gene set significance", type = "grid"))
g1 = grid.grabExpr(draw(hm.d, annotation_legend_list = lgd_list, merge_legend = TRUE, legend_grouping = "original"))

#**********************************************
#* Panel B
#**********************************************

mygenesets <- c(
  "HM_MITOTIC_SPINDLE",
 "HM_IL2_STAT5_SIGNALING", 
 "HM_APICAL_SURFACE", 
 "HM_P53_PATHWAY"
 )

plotList <- vector("list", length = length(mygenesets))

for(i in 1:length(mygenesets)){
 
 df <- data.frame(Var1 = exprs(scoreSSp)[mygenesets[i], ], 
         Timepoint = scoreSSp$Timepoint, 
         Arm = scoreSSp$Arm, 
         pCR = factor(scoreSSp$pCR, levels = c("pCR", "RD")), 
         Patient = scoreSSp$Patient_ID)
 df$Group <- paste(df$Arm, df$pCR, sep = " - ")
 df$Group <- factor(df$Group, levels = unique(df$Group)[c(1, 2, 4, 3)])
 
 df.mean <- df %>%
  group_by(Timepoint, Group) %>%
  summarise(
   n = n(), 
   mean_Var1 = mean(Var1, na.rm = TRUE), 
   .groups = "drop"
  ) %>%
  ungroup()
 
 
 ga <- ggplot(data = df, aes(x = Timepoint, y = Var1)) +
  geom_boxplot(fill = "grey80", coef = NULL, width = 0.3) +
  geom_line(aes(group = Patient, color = pCR)) +
  facet_wrap(~Group, nrow = 1) +
  scale_color_manual(name = "", values = carto_pal(12, "Pastel")[c(2, 1)]) +
  xlab("") +
  ylab("Score") +
  theme_pubr(border = T, legend = "none") +
  theme(panel.grid.major.y =  element_line(colour = "grey80", linewidth = 0.6, linetype = "dotted"), 
     strip.text = element_text(size = 12, face = "bold"), 
     plot.title = element_text(size = 14, face = "bold", hjust = 0.5), 
     axis.text.x = element_text(angle = 90, hjust = 1), 
     plot.margin = margin(t = 4, r = 3, b = 6, l = -100, unit = "pt")) +
  geom_point(data = df.mean, aes(x = Timepoint, y = mean_Var1), size = 3, color = "black", fill = "red3", pch = 21, inherit.aes = F) +
  geom_line(data = df.mean, aes(x = Timepoint, y = mean_Var1, group = Group), color = "red3", linewidth = 1.2, inherit.aes = F)
 
 dfp <- data.frame(Delta = exprs(deltas)[mygenesets[i], ], 
         Arm = deltas$Arm, 
         pCR = factor(deltas$pCR, levels = c("pCR", "RD")))
 testRes <- data.frame(Arm = c("CT", "CT/A"), rbind(res.ct.d[res.ct.d$Variable == mygenesets[i], ], 
                        res.cta.d[res.cta.d$Variable == mygenesets[i], ]))
 rng <- (max(dfp$Delta)-min(dfp$Delta))/8
 testRes$Ypos <- max(dfp$Delta) + rng
 testRes$Label <- paste0("OR = ", round(testRes$OR_Variable, 2), "\nFDR = ", round(testRes$FDR_Variable, 3))
 deltaInt <- round(res.int.d$P_Interaction[res.int.d$Variable == mygenesets[i]], 3)
 
 gb <- ggplot(data = dfp, aes(x = Arm, y = Delta, fill = pCR)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_boxplot(coef = NULL, width = 0.3) +
  scale_fill_manual(name = "", values = carto_pal(12, "Pastel")[c(2, 1)]) +
  xlab("") +
  ylab("\u0394 score") +
   ylim(min(dfp$Delta), max(dfp$Delta) + (rng * 1.8)) +
  theme_pubr(border = T, legend = "none") +
  theme(panel.grid.major.y =  element_line(colour = "grey80", linewidth = 0.6, linetype = "dotted"), 
     plot.title = element_text(face = "bold", hjust = 0.5, size = 12), 
     plot.margin = margin(t = 4, r = 3, b = -10, l = 0, unit = "pt")) +
  geom_text(data = testRes, aes(x = Arm, y = Ypos, label = Label), inherit.aes = F, size = 3.4, parse = F)
 
 plotList[[i]] <- plot_grid(ga, gb, ncol = 2, nrow = 1, rel_widths = c(1, 0.7))
 plotList[[i]] <- annotate_figure(plotList[[i]], top = text_grob(mygenesets[i], 
                                 color = "black", face = "bold", size = 14))
}

#**********************************************
#* Panel assembly
#**********************************************

leg <- get_legend(ggplot(data = dfp, aes(x = Arm, y = Delta, fill = pCR)) +
 geom_boxplot(coef = NULL, width = 0.3) +
 scale_fill_manual(name = "", values = carto_pal(12, "Pastel")[c(2, 1)]) +
 theme_pubr(border = T, legend = "top"))
gbox <- ggarrange(plotlist = plotList, nrow = length(mygenesets), ncol = 1, common.legend = T, legend.grob = leg)

combined <- plot_grid(g1, gbox, labels = c("A", "B"), label_size = 20, rel_widths = c(1, 0.7), label_x = c(0, -0.25))

title <- ggdraw() +
  draw_label("Figure 5",
             x = 0, hjust = 0,
             fontface = "bold", size = 20)
final_plot <- plot_grid(title, combined,
                        ncol = 1,
                        rel_heights = c(0.1, 1)) 
ggsave(file.path(here("results", "figures"), "Figure5.eps"),
       plot = final_plot,
       device = cairo_ps,
       width = 16, height = 13.5)
