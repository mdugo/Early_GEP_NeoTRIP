

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

#**********************************************
#* differential expression D1C2 vs Baseline
#**********************************************

data_dir <- here::here("data", "processed")

load(file.path(data_dir, "NeoTRIP_singscore_BL_D1C2.RData"))
paired <- names(which(table(scoreSS$Patient_ID) == 2))
scoreSSp <- scoreSS[, scoreSS$Patient_ID%in%paired]
scoreSSp <- scoreSSp[, order(scoreSSp$Arm, scoreSSp$Patient_ID)]

f <- factor(paste(scoreSSp$Timepoint, gsub("\\/", "", scoreSSp$Arm), scoreSSp$pCR, sep = "_"))
pat <- factor(paste0("p", scoreSSp$Patient_ID))
design <- model.matrix(~ 0 + f + pat)
design <- design[, !colnames(design) %in% c("patpTNBC87556", "patpTNBC88958", "patpTNBC89008")]

fit <- lmFit(scoreSSp, design)
contrast.matrix <- makeContrasts(
  fD1C2_CT_RD-fBaseline_CT_RD, 
         fD1C2_CT_pCR-fBaseline_CT_pCR, 
         fD1C2_CTA_RD-fBaseline_CTA_RD, 
         fD1C2_CTA_pCR-fBaseline_CTA_pCR, 
         (fD1C2_CT_pCR-fBaseline_CT_pCR)-(fD1C2_CT_RD-fBaseline_CT_RD), 
         (fD1C2_CTA_pCR-fBaseline_CTA_pCR)-(fD1C2_CTA_RD-fBaseline_CTA_RD), 
         (((fD1C2_CTA_pCR+fD1C2_CTA_RD)/2)-((fBaseline_CTA_pCR+fBaseline_CTA_RD)/2)) - (((fD1C2_CT_pCR+fD1C2_CT_RD)/2)-((fBaseline_CT_pCR+fBaseline_CT_RD)/2)), 
         ((fD1C2_CT_pCR+fD1C2_CT_RD)/2)-((fBaseline_CT_pCR+fBaseline_CT_RD)/2), 
         ((fD1C2_CTA_pCR+fD1C2_CTA_RD)/2)-((fBaseline_CTA_pCR+fBaseline_CTA_RD)/2), 
         levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
res1 <- topTable(fit2, coef = 1, number = nrow(scoreSS), adjust = "BH", sort.by = "none")
res2 <- topTable(fit2, coef = 2, number = nrow(scoreSS), adjust = "BH", sort.by = "none")
res3 <- topTable(fit2, coef = 3, number = nrow(scoreSS), adjust = "BH", sort.by = "none")
res4 <- topTable(fit2, coef = 4, number = nrow(scoreSS), adjust = "BH", sort.by = "none")
res5 <- topTable(fit2, coef = 5, number = nrow(scoreSS), adjust = "BH", sort.by = "none")
res6 <- topTable(fit2, coef = 6, number = nrow(scoreSS), adjust = "BH", sort.by = "none")
res7 <- topTable(fit2, coef = 7, number = nrow(scoreSS), adjust = "BH", sort.by = "none")
res8 <- topTable(fit2, coef = 8, number = nrow(scoreSS), adjust = "BH", sort.by = "none")
res9 <- topTable(fit2, coef = 9, number = nrow(scoreSS), adjust = "BH", sort.by = "none")

sig.all <- unique(c(res1$GeneSet[res1$adj.P.Val<0.05], 
     res2$GeneSet[res2$adj.P.Val<0.05], 
     res3$GeneSet[res3$adj.P.Val<0.05], 
     res4$GeneSet[res4$adj.P.Val<0.05], 
     res8$GeneSet[res8$adj.P.Val<0.05], 
     res9$GeneSet[res9$adj.P.Val<0.05]))

fcMat <- cbind(res8$logFC, res9$logFC, res1$logFC, res2$logFC, res3$logFC, res4$logFC)
fdrMat <- cbind(res8$adj.P.Val, res9$adj.P.Val, res1$adj.P.Val, res2$adj.P.Val, res3$adj.P.Val, res4$adj.P.Val)
colnames(fcMat) <- c("CT", "CTA", "CT_RD", "CT_pCR", "CTA_RD", "CTA_pCR")
rownames(fcMat) <- res1$GeneSet
colnames(fdrMat) <- c("CT", "CTA", "CT_RD", "CT_pCR", "CTA_RD", "CTA_pCR")
rownames(fdrMat) <- res1$GeneSet
fcMat <- fcMat[sig.all, ]
fdrMat <- fdrMat[sig.all, ]
fdrMat2 <- ifelse(fdrMat < 0.05, "\U2217\U2217", ifelse(fdrMat < 0.1, "\U2217", ""))

intPvalNum <- structure(res7$P.Value, names = res7$GeneSet)
intPvalNum <- intPvalNum[rownames(fcMat)]
textra <- rowAnnotation("P interaction" = anno_text(round(intPvalNum, 2), show_name = TRUE, gp = gpar(fontsize = 8, fontface = ifelse(intPvalNum<0.05, "bold", "plain"))), 
       annotation_label = expression(P[interaction]))

annotCol <- structure(c(carto_pal(12, "Vivid")[1:9], "grey75"), names = names(table(fData(scoreSS)$Annotation)))
fdata <- fData(scoreSS)
fdata <- fdata[rownames(fcMat), ]
fdata$Annotation <- droplevels(fdata$Annotation)
la <- rowAnnotation(Annotation = fdata$Annotation, 
     col = list(Annotation = annotCol), 
     annotation_legend_param = list(Annotation = list(at = levels(fdata$Annotation)[c(1, 5, 6, 8, 3, 9, 10, 4, 2, 7)])))

fcCol <- colorRamp2(c(min(fcMat), 0, max(fcMat)), brewer.pal(11, "RdBu")[c(10, 6, 2)])

hm.all <- Heatmap(fcMat[, 1:2], 
    name = "\U0394 score", 
    col = fcCol, 
    cluster_columns = F, 
    clustering_method_rows = "ward.D2", 
    row_split = fdata$Annotation, 
    row_title = NULL, 
    show_row_dend = F, 
    row_names_side = "left", 
    row_names_gp = gpar(fontsize = 7), 
    column_labels = c("CT", "CT/A"), 
    column_title = "CT+CT/A", 
    rect_gp = gpar(col = "grey40", lwd = 0.6), 
    right_annotation = textra, 
    left_annotation = la, 
    width = unit(0.55, "in"), 
    heatmap_legend_param = list(title_gp = gpar(fontsize = 10, fontface = "bold")), 
    cell_fun = function(j, i, x, y, width, height, fill) {
     grid.text(fdrMat2[, 1:2][i, j], x, y, gp = gpar(fontsize = 9, fontface = "bold"))})


intPvalNum.ct <- structure(res5$P.Value, names = res5$GeneSet)
intPvalNum.ct <- intPvalNum.ct[rownames(fcMat)]
textra.ct <- rowAnnotation("P interaction" = anno_text(round(intPvalNum.ct, 2), show_name = TRUE, gp = gpar(fontsize = 8, fontface = ifelse(intPvalNum.ct<0.05, "bold", "plain"))), 
       annotation_label = expression(P[interaction]))

intPvalNum.cta <- structure(res6$P.Value, names = res6$GeneSet)
intPvalNum.cta <- intPvalNum.cta[rownames(fcMat)]
textra.cta <- rowAnnotation("P interaction" = anno_text(round(intPvalNum.cta, 2), show_name = TRUE, gp = gpar(fontsize = 8, fontface = ifelse(intPvalNum.cta<0.05, "bold", "plain"))), 
       annotation_label = expression(P[interaction]))

annotCol <- structure(c(carto_pal(12, "Vivid")[1:9], "grey75"), names = names(table(fData(scoreSS)$Annotation)))
fdata <- fData(scoreSS)
fdata <- fdata[rownames(fcMat), ]
fdata$Annotation <- droplevels(fdata$Annotation)
la <- rowAnnotation(Annotation = fdata$Annotation, 
     col = list(Annotation = annotCol), 
     annotation_legend_param = list(Annotation = list(at = levels(fdata$Annotation)[c(2, 7, 1, 9, 10, 4, 5, 6, 8, 3)])))

fcCol <- colorRamp2(c(min(fcMat), 0, max(fcMat)), brewer.pal(11, "RdBu")[c(10, 6, 2)])

hm.ct <- Heatmap(fcMat[, 3:4], 
    name = "\U0394 score", 
    show_heatmap_legend = T,
    col = fcCol, 
    cluster_columns = F, 
   clustering_method_rows = "ward.D2", 
    row_split = fdata$Annotation, 
    row_title = NULL, 
    show_row_dend = F, 
    row_names_side = "left", 
    row_names_gp = gpar(fontsize = 7), 
    column_labels = c("RD", "pCR"), 
    column_title = "CT", 
    rect_gp = gpar(col = "grey40", lwd = 0.6), 
    right_annotation = textra.ct, 
    left_annotation = la, 
    width = unit(0.55, "in"), 
    heatmap_legend_param = list(title_gp = gpar(fontsize = 10, fontface = "bold")), 
    cell_fun = function(j, i, x, y, width, height, fill) {
     grid.text(fdrMat2[, 3:4][i, j], x, y, gp = gpar(fontsize = 9, fontface = "bold"))})

hm.cta <- Heatmap(fcMat[, 5:6], 
    name = "\U0394 score2",
    show_heatmap_legend = F,
    col = fcCol, 
    cluster_columns = F, 
    clustering_method_rows = "ward.D2", 
    row_split = fdata$Annotation, 
    row_title = NULL, 
    show_row_dend = F, 
    row_names_side = "left", 
    row_names_gp = gpar(fontsize = 7), 
    column_labels = c("RD", "pCR"), 
    column_title = "CT/A", 
    rect_gp = gpar(col = "grey40", lwd = 0.6), 
    right_annotation = textra.cta, 
    width = unit(0.55, "in"), 
    heatmap_legend_param = list(title_gp = gpar(fontsize = 10, fontface = "bold")), 
    cell_fun = function(j, i, x, y, width, height, fill) {
     grid.text(fdrMat2[, 5:6][i, j], x, y, gp = gpar(fontsize = 9, fontface = "bold"))})

lgd_list = list(
  Legend(labels = c("*: FDR < 0.1", "**: FDR < 0.05"), title = "Gene set significance", type = "grid"))

#**********************************************
#* Panel assembly
#**********************************************

showtext.auto()
g1 <- grid.grabExpr(draw(hm.all, annotation_legend_list = lgd_list, merge_legend = TRUE))
g2 <- grid.grabExpr(draw(hm.ct + hm.cta, gap = unit(7, "mm"), annotation_legend_list = lgd_list, merge_legend = TRUE))
combined <- plot_grid(g1, g2, labels = c("A", "B"), label_size = 16)

title <- ggdraw() +
  draw_label("Supplementary Figure 3",
             x = 0, hjust = 0,
             fontface = "bold", size = 16)
final_plot <- plot_grid(title, combined,
                        ncol = 1,
                        rel_heights = c(0.1, 1)) 
ggsave(file.path(here("results", "figures"), "SuppFigure3.eps"),
       plot = final_plot,
       device = cairo_ps,
       width = 12, height = 12.5)
