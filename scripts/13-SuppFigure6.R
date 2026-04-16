
library(Biobase)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(limma)
library(showtext)
library(rcartocolor)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
source(file.path(here::here("scripts"), "00a-helper_functions.R"))

data_dir <- here::here("data", "processed")

#***********************************************************************
# Panel A-B
#***********************************************************************

# Not adjusted

load(file.path(data_dir, "instaprism_deconvolution.RData"))
load(file.path(data_dir, "NeoTRIP_singscore_BL_D1C2.RData"))
deconv_tpm@Post.ini.ct@theta <- deconv_tpm@Post.ini.ct@theta*100
identical(colnames(scoreSS), colnames(deconv_tpm@Post.ini.ct@theta))
pData(scoreSS) <- cbind(pData(scoreSS), t(deconv_tpm@Post.ini.ct@theta))
colnames(pData(scoreSS))[17] <- "Cancer_Epithelial"
scoreSS$Tumor_fraction <- scoreSS$Cancer_Epithelial
scoreSS$Stromal_fraction <- rowSums(pData(scoreSS)[, c(10:12)])
scoreSS$Immune_fraction <- rowSums(pData(scoreSS)[, c(13:15)])

bl <- scoreSS[, scoreSS$Timepoint == "Baseline"]
res.ct.bl <- logistic_regression(response = "pCR.num", 
                data = bl[, bl$Arm == "CT"], 
                sort = F)
res.cta.bl <- logistic_regression(response = "pCR.num", 
                data = bl[, bl$Arm == "CT/A"], 
                sort = F)
res.int.bl <- interaction_logistic_regression(response = "pCR.num", 
                      covariate_interaction = "Arm", 
                      data = bl, 
                      sort = F)
d1c2 <- scoreSS[, scoreSS$Timepoint == "D1C2"]
res.ct.d1c2 <- logistic_regression(response = "pCR.num", 
                 data = d1c2[, d1c2$Arm == "CT"], 
                 sort = F)
res.cta.d1c2 <- logistic_regression(response = "pCR.num", 
                 data = d1c2[, d1c2$Arm == "CT/A"], 
                 sort = F)
res.int.d1c2 <- interaction_logistic_regression(response = "pCR.num", 
                       covariate_interaction = "Arm", 
                       data = d1c2, 
                       sort = F)


# Adjusted

res.ct.bl2 <- logistic_regression(response = "pCR.num", 
                covariates = colnames(pData(scoreSS))[19:20], 
                data = bl[, bl$Arm == "CT"], 
                sort = F)
res.cta.bl2 <- logistic_regression(response = "pCR.num", 
                 covariates = colnames(pData(scoreSS))[19:20], 
                 data = bl[, bl$Arm == "CT/A"], 
                 sort = F)
res.int.bl2 <- interaction_logistic_regression(response = "pCR.num", 
                       covariates = colnames(pData(scoreSS))[19:20], 
                       covariate_interaction = "Arm", 
                       data = bl, 
                       sort = F)

res.ct.d1c22 <- logistic_regression(response = "pCR.num", 
                 covariates = colnames(pData(scoreSS))[19:20], 
                 data = d1c2[, d1c2$Arm == "CT"], 
                 sort = F)
res.cta.d1c22 <- logistic_regression(response = "pCR.num", 
                  covariates = colnames(pData(scoreSS))[19:20], 
                  data = d1c2[, d1c2$Arm == "CT/A"], 
                  sort = F)
res.int.d1c22 <- interaction_logistic_regression(response = "pCR.num", 
                        covariates = colnames(pData(scoreSS))[19:20], 
                        covariate_interaction = "Arm", 
                        data = d1c2, 
                        sort = F)


# heatmap panel A

sig.bl <- unique(c(res.ct.bl$Variable[res.ct.bl$FDR_Variable<0.05], 
         res.cta.bl$Variable[res.cta.bl$FDR_Variable<0.05], 
         res.ct.bl2$Variable[res.ct.bl2$FDR_Variable<0.05], 
         res.cta.bl2$Variable[res.cta.bl2$FDR_Variable<0.05]))
sig.d1c2 <- unique(c(res.ct.d1c2$Variable[res.ct.d1c2$FDR_Variable<0.05], 
          res.cta.d1c2$Variable[res.cta.d1c2$FDR_Variable<0.05], 
          res.ct.d1c22$Variable[res.ct.d1c22$FDR_Variable<0.05], 
          res.cta.d1c22$Variable[res.cta.d1c22$FDR_Variable<0.05]))

annotCol <- structure(c(carto_pal(12, "Vivid")[1:9], "grey75"), names = names(table(fData(bl)$Annotation)))

orMat <- cbind(res.ct.bl$OR_Variable, 
       res.cta.bl$OR_Variable, 
       res.ct.bl2$OR_Variable, 
       res.cta.bl2$OR_Variable)
rownames(orMat) <- res.ct.bl$Variable
colnames(orMat) <- c("BL_CT", "BL_CTA", "BL_CT_Dec", "BL_CTA_Dec")

fdrMat <- cbind(res.ct.bl$FDR_Variable, 
       res.cta.bl$FDR_Variable, 
       res.ct.bl2$FDR_Variable, 
       res.cta.bl2$FDR_Variable)
rownames(fdrMat) <- res.ct.bl$Variable
colnames(fdrMat) <- c("BL_CT", "BL_CTA", "BL_CT_Dec", "BL_CTA_Dec")
fdrMat <- fdrMat[sig.bl, ]
orMat <- orMat[rownames(fdrMat), ]

fdrMat2 <- ifelse(fdrMat<0.05, "\U2217\U2217", ifelse(fdrMat<0.1, "\U2217", ""))
logOrMat <- log10(orMat+0.01)


intPvalNum <- res.int.bl$P_Interaction
names(intPvalNum) <- res.int.bl$Variable
intPvalNum <- intPvalNum[rownames(logOrMat)]
intPvalNum <- round(intPvalNum, 2)
textra <- rowAnnotation("P interaction" = anno_text(intPvalNum, show_name = TRUE, gp = gpar(fontsize = 12, fontface = ifelse(intPvalNum<0.05, "bold", "plain"))), 
           annotation_label = expression(P[interaction]))

intPvalNum2 <- res.int.bl2$P_Interaction
names(intPvalNum2) <- res.int.bl2$Variable
intPvalNum2 <- intPvalNum2[rownames(logOrMat)]
intPvalNum2 <- round(intPvalNum2, 2)
textra2 <- rowAnnotation("P interaction" = anno_text(intPvalNum2, show_name = TRUE, gp = gpar(fontsize = 12, fontface = ifelse(intPvalNum2<0.05, "bold", "plain"))), 
           annotation_label = expression(P[interaction]))

fdata <- fData(bl)
fdata <- fdata[rownames(logOrMat), ]
fdata$Annotation <- droplevels(fdata$Annotation)
la <- rowAnnotation(Annotation = fdata$Annotation, 
         col = list(Annotation = annotCol), 
         annotation_legend_param = list(Annotation = list(at = levels(fdata$Annotation)[c(6, 8, 7, 4, 3, 5, 2, 1)])))

orCol <- colorRamp2(c(min(logOrMat), 0, max(logOrMat)), brewer.pal(11, "RdBu")[c(10, 6, 2)])
hm.bl <- Heatmap(logOrMat[, 1:2], 
        name = "OR", 
        col = orCol, 
        row_split = fdata$Annotation, 
        row_title = NULL, 
        heatmap_legend_param = list(title = expression("log"[10](OR)), 
                      title_gp = gpar(fontsize = 10, fontface = "bold")), 
        cluster_columns = F, 
        cluster_column_slices = T, 
        cluster_rows = T, 
        clustering_method_rows = "ward.D2", 
        show_row_dend = F, 
        row_names_side = "left", 
        row_names_gp = gpar(fontsize = 11), 
        column_names_gp = gpar(fontsize = 13), 
        column_labels = c("CT", "CT/A"), 
        column_title = "Not adjusted", 
        rect_gp = gpar(col = "grey40", lwd = 0.6), 
        right_annotation = textra, 
        left_annotation = la, 
        width = unit(0.5, "in"), 
        cell_fun = function(j, i, x, y, width, height, fill) {
         grid.text(fdrMat2[, 1:2][i, j], x, y, gp = gpar(fontsize = 9, fontface = "bold"))}
)

hm.bl2 <- Heatmap(logOrMat[, 3:4], 
        name = "OR2", 
        col = orCol, 
        show_heatmap_legend = F, 
        row_split = fdata$Annotation, 
        row_title = NULL, 
        heatmap_legend_param = list(title = expression("log"[10](OR)), 
                      title_gp = gpar(fontsize = 10, fontface = "bold")), 
        cluster_columns = F, 
        cluster_rows = F, 
        show_row_dend = F, 
        show_row_names = F, 
        row_names_side = "left", 
        row_names_gp = gpar(fontsize = 11), 
        column_names_gp = gpar(fontsize = 13), 
        column_labels = c("CT", "CT/A"), 
        column_title = "Adjusted", 
        rect_gp = gpar(col = "grey40", lwd = 0.6), 
        right_annotation = textra2, 
        width = unit(0.5, "in"), 
        cell_fun = function(j, i, x, y, width, height, fill) {
         grid.text(fdrMat2[, 3:4][i, j], x, y, gp = gpar(fontsize = 9, fontface = "bold"))}
)

lgd_list = list(
 Legend(labels = c("*: FDR < 0.1", "**: FDR < 0.05"), title = "Gene set significance", type = "grid"))
showtext.auto()
g1 = grid.grabExpr(draw(hm.bl+hm.bl2, annotation_legend_list = lgd_list, gap = unit(12, "mm"), column_title = "Baseline", merge_legend = T))

# heatmap panel B

orMat <- cbind(res.ct.d1c2$OR_Variable, 
       res.cta.d1c2$OR_Variable, 
       res.ct.d1c22$OR_Variable, 
       res.cta.d1c22$OR_Variable)
rownames(orMat) <- res.ct.d1c2$Variable
colnames(orMat) <- c("CT", "CTA", "CT_Dec", "CTA_Dec")

fdrMat <- cbind(res.ct.d1c2$FDR_Variable, 
       res.ct.d1c2$FDR_Variable, 
       res.ct.d1c22$FDR_Variable, 
       res.ct.d1c22$FDR_Variable)
rownames(fdrMat) <- res.ct.d1c2$Variable
colnames(fdrMat) <- c("CT", "CTA", "CT_Dec", "CTA_Dec")
fdrMat <- fdrMat[sig.d1c2, ]
orMat <- orMat[rownames(fdrMat), ]

fdrMat2 <- ifelse(fdrMat<0.05, "\U2217\U2217", ifelse(fdrMat<0.1, "\U2217", ""))
logOrMat <- log10(orMat+0.01)


intPvalNum <- res.int.d1c2$P_Interaction
names(intPvalNum) <- res.int.d1c2$Variable
intPvalNum <- intPvalNum[rownames(logOrMat)]
intPvalNum <- round(intPvalNum, 2)
textra <- rowAnnotation("P interaction" = anno_text(intPvalNum, show_name = TRUE, gp = gpar(fontsize = 12, fontface = ifelse(intPvalNum<0.05, "bold", "plain"))), 
           annotation_label = expression(P[interaction]))

intPvalNum2 <- res.int.d1c22$P_Interaction
names(intPvalNum2) <- res.int.d1c22$Variable
intPvalNum2 <- intPvalNum2[rownames(logOrMat)]
intPvalNum2 <- round(intPvalNum2, 2)
textra2 <- rowAnnotation("P interaction" = anno_text(intPvalNum2, show_name = TRUE, gp = gpar(fontsize = 12, fontface = ifelse(intPvalNum2<0.05, "bold", "plain"))), 
            annotation_label = expression(P[interaction]))

fdata <- fData(d1c2)
fdata <- fdata[rownames(logOrMat), ]
fdata$Annotation <- droplevels(fdata$Annotation)
la <- rowAnnotation(Annotation = fdata$Annotation, 
         col = list(Annotation = annotCol), 
         annotation_legend_param = list(Annotation = list(at = levels(fdata$Annotation)[c(3, 1, 4, 2, 5)])))

orCol <- colorRamp2(c(min(logOrMat), 0, max(logOrMat)), brewer.pal(11, "RdBu")[c(10, 6, 2)])
hm.d1c2 <- Heatmap(logOrMat[, 1:2], 
        name = "OR", 
        col = orCol, 
        row_split = fdata$Annotation, 
        row_title = NULL, 
        heatmap_legend_param = list(title = expression("log"[10](OR)), 
                      title_gp = gpar(fontsize = 10, fontface = "bold")), 
        cluster_columns = F, 
        cluster_column_slices = T, 
        cluster_rows = T, 
        clustering_method_rows = "ward.D2", 
        show_row_dend = F, 
        row_names_side = "left", 
        row_names_gp = gpar(fontsize = 11), 
        column_names_gp = gpar(fontsize = 13), 
        column_labels = c("CT", "CT/A"), 
        column_title = "Not adjusted", 
        rect_gp = gpar(col = "grey40", lwd = 0.6), 
        right_annotation = textra, 
        left_annotation = la, 
        width = unit(0.6, "in"), 
        #height = unit(6, "in"), 
        cell_fun = function(j, i, x, y, width, height, fill) {
         grid.text(fdrMat2[, 1:2][i, j], x, y, gp = gpar(fontsize = 9, fontface = "bold"))}
)

orCol <- colorRamp2(c(min(logOrMat[, 3:4]), 0, max(logOrMat[, 3:4])), brewer.pal(11, "RdBu")[c(10, 6, 2)])
hm.d1c22 <- Heatmap(logOrMat[, 3:4], 
        name = "OR2", 
        col = orCol, 
        show_heatmap_legend = F, 
        row_split = fdata$Annotation, 
        row_title = NULL, 
        heatmap_legend_param = list(title = expression("log"[10](OR)), 
                      title_gp = gpar(fontsize = 10, fontface = "bold")), 
        cluster_columns = F, 
        cluster_rows = F, 
        show_row_dend = F, 
        show_row_names = F, 
        row_names_side = "left", 
        row_names_gp = gpar(fontsize = 11), 
        column_names_gp = gpar(fontsize = 13), 
        column_labels = c("CT", "CT/A"), 
        column_title = "Adjusted", 
        rect_gp = gpar(col = "grey40", lwd = 0.6), 
        right_annotation = textra2, 
        width = unit(0.6, "in"), 
        cell_fun = function(j, i, x, y, width, height, fill) {
         grid.text(fdrMat2[, 3:4][i, j], x, y, gp = gpar(fontsize = 9, fontface = "bold"))}
)

showtext.auto()
g2 = grid.grabExpr(draw(hm.d1c2+hm.d1c22, annotation_legend_list = lgd_list, gap = unit(12, "mm"), column_title = "D1C2", merge_legend = T))
g12 <- plot_grid(g1, g2, ncol = 1, labels = c("A", "B"), label_size = 20)


#***********************************************************************
# Panel C
#***********************************************************************

paired <- names(which(table(scoreSS$Patient_ID) == 2))
scoreSSp <- scoreSS[, scoreSS$Patient_ID %in% paired]
f <- factor(paste(scoreSSp$Timepoint, gsub("\\/", "", scoreSSp$Arm), scoreSSp$pCR, sep = "_"))
pat <- factor(scoreSSp$Patient_ID)
stromalfrac <- scoreSSp$Stromal_fraction
immunefrac <- scoreSSp$Immune_fraction
design <- model.matrix(~0+f+pat)
design <- design[, !colnames(design) %in% c("patTNBC87556", "patTNBC88958", "patTNBC89008")]

fit <- lmFit(scoreSSp, design)
contrast.matrix <- makeContrasts(((fD1C2_CT_pCR+fD1C2_CT_RD)/2)-((fBaseline_CT_pCR+fBaseline_CT_RD)/2), 
                 ((fD1C2_CTA_pCR+fD1C2_CTA_RD)/2)-((fBaseline_CTA_pCR+fBaseline_CTA_RD)/2), 
                 (((fD1C2_CTA_pCR+fD1C2_CTA_RD)/2)-((fBaseline_CTA_pCR+fBaseline_CTA_RD)/2)) - (((fD1C2_CT_pCR+fD1C2_CT_RD)/2)-((fBaseline_CT_pCR+fBaseline_CT_RD)/2)), 
                 levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
res1 <- topTable(fit2, coef = 1, number = nrow(scoreSSp), adjust = "BH", sort.by = "none")
res2 <- topTable(fit2, coef = 2, number = nrow(scoreSSp), adjust = "BH", sort.by = "none")
res3 <- topTable(fit2, coef = 3, number = nrow(scoreSSp), adjust = "BH", sort.by = "none")

sig.all <- unique(c(res1$GeneSet[res1$adj.P.Val<0.05], 
         res2$GeneSet[res2$adj.P.Val<0.05]))
sig.all <- sig.all[-grep("cTME", sig.all)]

fcMat <- cbind(res1$logFC, res2$logFC)
fdrMat <- cbind(res1$adj.P.Val, res2$adj.P.Val)
colnames(fcMat) <- c("CT", "CTA")
rownames(fcMat) <- res1$GeneSet
colnames(fdrMat) <- c("CT", "CTA")
rownames(fdrMat) <- res1$GeneSet
fcMat <- fcMat[sig.all, ]
fdrMat <- fdrMat[sig.all, ]
fdrMat2 <- ifelse(fdrMat<0.05, "\U2217\U2217", ifelse(fdrMat<0.1, "\U2217", ""))

intPvalNum <- structure(res3$P.Value, names = res3$GeneSet)
intPvalNum <- intPvalNum[sig.all]

textra <- rowAnnotation("P interaction" = anno_text(round(intPvalNum, 3), show_name = TRUE, gp = gpar(fontsize = 10, fontface = ifelse(intPvalNum<0.05, "bold", "plain"))), 
             annotation_label = expression(P[interaction]))

annotCol <- structure(c(carto_pal(12, "Vivid")[1:9], "grey75"), names = names(table(fData(scoreSSp)$Annotation)))
fdata <- fData(scoreSSp)
fdata <- fdata[rownames(fcMat), ]
fdata$Annotation <- droplevels(fdata$Annotation)
la <- rowAnnotation(Annotation = fdata$Annotation, 
         col = list(Annotation = annotCol), 
         annotation_legend_param = list(Annotation = list(at = levels(fdata$Annotation)[c(1, 2, 10, 4, 9, 5, 6, 8, 3, 7)])))

fcCol <- colorRamp2(c(min(fcMat), 0, max(fcMat)), brewer.pal(11, "RdBu")[c(10, 6, 2)])

hm1 <- Heatmap(fcMat, 
        name = "\U0394 score", 
        col = fcCol, 
        cluster_columns = F, 
        clustering_method_rows = "ward.D2", 
        row_split = fdata$Annotation, 
        row_title = NULL, 
        show_row_dend = F, 
        row_names_side = "left", 
        row_names_gp = gpar(fontsize = 10), 
        column_labels = c("CT", "CT/A"), 
        column_title = "Not adjusted", 
        rect_gp = gpar(col = "grey40", lwd = 0.6), 
        right_annotation = textra, 
        left_annotation = la, 
        width = unit(0.55, "in"), 
        heatmap_legend_param = list(title_gp = gpar(fontsize = 9, fontface = "bold")), 
        cell_fun = function(j, i, x, y, width, height, fill) {
         grid.text(fdrMat2[i, j], x, y, gp = gpar(fontsize = 9, fontface = "bold"))})



###

design <- model.matrix(~0+f+pat+stromalfrac + immunefrac)
design <- design[, !colnames(design) %in% c("patTNBC87556", "patTNBC88958", "patTNBC89008")]

fit <- lmFit(scoreSSp, design)
contrast.matrix <- makeContrasts(((fD1C2_CT_pCR+fD1C2_CT_RD)/2)-((fBaseline_CT_pCR+fBaseline_CT_RD)/2), 
                 ((fD1C2_CTA_pCR+fD1C2_CTA_RD)/2)-((fBaseline_CTA_pCR+fBaseline_CTA_RD)/2), 
                 (((fD1C2_CTA_pCR+fD1C2_CTA_RD)/2)-((fBaseline_CTA_pCR+fBaseline_CTA_RD)/2)) - (((fD1C2_CT_pCR+fD1C2_CT_RD)/2)-((fBaseline_CT_pCR+fBaseline_CT_RD)/2)), 
                 levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
res1 <- topTable(fit2, coef = 1, number = nrow(scoreSSp), adjust = "BH", sort.by = "none")
res2 <- topTable(fit2, coef = 2, number = nrow(scoreSSp), adjust = "BH", sort.by = "none")
res3 <- topTable(fit2, coef = 3, number = nrow(scoreSSp), adjust = "BH", sort.by = "none")


fcMat <- cbind(res1$logFC, res2$logFC)
fdrMat <- cbind(res1$adj.P.Val, res2$adj.P.Val)
colnames(fcMat) <- c("CT", "CTA")
rownames(fcMat) <- res1$GeneSet
colnames(fdrMat) <- c("CT", "CTA")
rownames(fdrMat) <- res1$GeneSet
fcMat <- fcMat[sig.all, ]
fdrMat <- fdrMat[sig.all, ]
fdrMat2 <- ifelse(fdrMat<0.05, "\U2217\U2217", ifelse(fdrMat<0.1, "\U2217", ""))

intPvalNum <- structure(res3$P.Value, names = res3$GeneSet)
intPvalNum <- intPvalNum[sig.all]

textra <- rowAnnotation("P interaction" = anno_text(round(intPvalNum, 3), show_name = TRUE, gp = gpar(fontsize = 10, fontface = ifelse(intPvalNum<0.05, "bold", "plain"))), 
           annotation_label = expression(P[interaction]))

fcCol <- colorRamp2(c(min(fcMat), 0, max(fcMat)), brewer.pal(11, "RdBu")[c(10, 6, 2)])

hm2 <- Heatmap(fcMat, 
        name = "\U0394 score\nadjusted", 
        col = fcCol, 
        cluster_columns = F, 
        clustering_method_rows = "ward.D2", 
        row_title = NULL, 
        show_row_dend = F, 
        row_names_side = "left", 
        row_names_gp = gpar(fontsize = 10), 
        column_labels = c("CT", "CT/A"), 
        column_title = "Adjusted", 
        rect_gp = gpar(col = "grey40", lwd = 0.6), 
        right_annotation = textra, 
        width = unit(0.55, "in"), 
        heatmap_legend_param = list(title_gp = gpar(fontsize = 9, fontface = "bold")), 
        cell_fun = function(j, i, x, y, width, height, fill) {
         grid.text(fdrMat2[i, j], x, y, gp = gpar(fontsize = 9, fontface = "bold"))})


showtext.auto()
g4 = grid.grabExpr(draw(hm1+hm2, annotation_legend_list = lgd_list, gap = unit(15, "mm"), merge_legend = TRUE))


#**************************
# Panel assembly
#**************************

combined <- plot_grid(g12, g4, ncol = 2, nrow = 1, labels = c("", "C"), label_size = 20)
title <- ggdraw() +
  draw_label("Supplementary Figure 6",
             x = 0, hjust = 0,
             fontface = "bold", size = 20)
final_plot <- plot_grid(title, combined,
                        ncol = 1,
                        rel_heights = c(0.1, 1))  # spazio sopra
ggsave(file.path(here::here("results", "figures"), "SuppFigure6.eps"),
       plot = final_plot,
       device = cairo_ps,
       width = 19, height = 15)


