
library(Biobase)
library(singscore)
library(ComplexHeatmap)
library(openxlsx)
library(igraph)
library(rcartocolor)
library(data.table)

data_dir <- here::here("data", "processed")

load(file.path(data_dir, "NeoTRIP_singscore_BL_D1C2.RData"))

correl <- cor(t(exprs(scoreSS)), method="pearson")
correl.txt <- ifelse(correl>=0.8, "+", ifelse(correl>=0.7, "-", ""))

corrs <- cor(t(exprs(scoreSS)), method = 'pearson')
corrs[upper.tri(corrs, diag = TRUE)] <- NA
corrs <- as.data.table(corrs)[, Var1 := rownames(corrs)]
corrs <- melt(corrs, id.vars = 'Var1', variable = 'Var2', value = 'corr')[!is.na(corr)]
corrs[, corr := abs(corr)] # high positive or negative correlation
corrs <- corrs[, TooHighCorr := corr >= 0.8]
corrGrph <- graph_from_data_frame(corrs[(TooHighCorr)], directed = FALSE)

clusters <- membership(cluster_edge_betweenness(corrGrph))
clusters <- data.frame(GeneSet=names(clusters), Cluster=paste0("C", as.vector(clusters)))
fdata <- fData(scoreSS)
fdata <- merge(fdata, clusters, by="GeneSet", all.x=T)
rownames(fdata) <- fdata$GeneSet
fdata <- fdata[rownames(scoreSS), ]
identical(rownames(fdata), rownames(scoreSS))
fData(scoreSS) <- fdata


annotCol <- structure(c(carto_pal(12, "Vivid")[1:9], "grey75"), names=names(table(fData(scoreSS)$Annotation)))

ra <- rowAnnotation(Annotation=fData(scoreSS)$Annotation, 
         col=list(Annotation=annotCol), 
         show_legend=F)

ta <- HeatmapAnnotation(Annotation=fData(scoreSS)$Annotation, 
         col=list(Annotation=annotCol), 
         annotation_legend_param = list(Annotation = list(at=levels(fData(scoreSS)$Annotation)[c(7, 2, 5, 10, 4, 9, 8, 6, 3, 1)])))


hm<-Heatmap(correl, 
    name = "Spearman\ncorrelation", 
    row_names_gp = gpar(fontsize = 6), 
    column_names_gp = gpar(fontsize = 6), 
    clustering_method_rows = "ward.D2", 
    clustering_method_columns = "ward.D2", 
    right_annotation = ra, 
    top_annotation = ta, 
    row_split = fData(scoreSS)$Annotation, 
    column_split = fData(scoreSS)$Annotation, 
    column_title = NULL, 
    row_title = NULL, 
    cell_fun = function(j, i, x, y, width, height, fill) {
     grid.text(correl.txt[i, j], x, y, gp = gpar(fontsize = 6))
    })

g1 = grid.grabExpr(draw(hm, annotation_legend_list = lgd_list, merge_legend=TRUE, legend_grouping="original"))

title <- ggdraw() +
  draw_label("Supplementary Figure 1",
             x = 0, hjust = 0,
             fontface = "bold", size = 16)
final_plot <- plot_grid(title, g1,
                        ncol = 1,
                        rel_heights = c(0.1, 1)) 
ggsave(file.path(here("results", "figures"), "SuppFigure1.eps"),
       plot = final_plot,
       device = cairo_ps,
       width = 14, height = 11.5)
