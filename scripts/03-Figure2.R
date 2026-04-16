

library(Biobase)
library(ComplexHeatmap)
library(circlize)
library(rcartocolor)
library(ggplot2)
library(ggpubr)
library(ggplotify)
library(cowplot)
library(RColorBrewer)
library(showtext)
library(colorspace)
library(here)
source(file.path(here("scripts"), "00a-helper_functions.R"))

data_dir <- here("data", "processed")

load(file.path(data_dir, "NeoTRIP_singscore_BL_D1C2.RData"))

#*******************************
# Panel A
#*******************************

scoreSS <- scoreSS[, scoreSS$Timepoint   ==  "Baseline"]

res.ct.bl <- logistic_regression(response = "pCR.num", 
                data = scoreSS[, scoreSS$Arm  == "CT"], 
                sort = F)
res.cta.bl <- logistic_regression(response = "pCR.num", 
                data = scoreSS[, scoreSS$Arm  == "CT/A"], 
                sort = F)
res.int.bl <- interaction_logistic_regression(response = "pCR.num", 
                      covariate_interaction = "Arm", 
                      data = scoreSS, 
                      sort=F)

fullRes <- merge(res.ct.bl, res.cta.bl, by="Variable")
colnames(fullRes) <- gsub("_Variable.y", "_CTA", gsub("_Variable.x", "_CT", colnames(fullRes)))
fullRes <- merge(fullRes, res.int.bl[, c("Variable", "P_Interaction", "FDR_Interaction")], by="Variable")
write.table(fullRes, file = file.path(here("results", "tables"), "01-pCR_prediciton_baseline.txt"))

# heatmap

res.ct.bl$Arm <- "CT"
res.cta.bl$Arm <- "CT/A"
sig <- unique(c(res.ct.bl$Variable[res.ct.bl$FDR_Variable < 0.05],
                res.cta.bl$Variable[res.cta.bl$FDR_Variable < 0.05]))

annotCol <- structure(c(carto_pal(12, "Vivid")[1:9], "grey75"), names = names(table(fData(scoreSS)$Annotation)))

orMat <- cbind(res.ct.bl$OR_Variable, 
       res.cta.bl$OR_Variable)
rownames(orMat) <- res.ct.bl$Variable
colnames(orMat) <- c("BL_CT", "BL_CTA")

fdrMat <- cbind(res.ct.bl$FDR_Variable, 
       res.cta.bl$FDR_Variable)
rownames(fdrMat) <- res.ct.bl$Variable
colnames(fdrMat) <- c("BL_CT", "BL_CTA")

orMat <- orMat[sig, ]
fdrMat <- fdrMat[sig, ]
fdrMat2 <- ifelse(fdrMat < 0.05, "\U2217\U2217", ifelse(fdrMat < 0.1, "\U2217", ""))
logOrMat <- log10(orMat)


intPval <- ifelse(res.int.bl$P_Interaction<0.01, "p < 0.01", ifelse(res.int.bl$P_Interaction<0.05, "p < 0.05", "n.s."))
intPval <- factor(intPval, levels=c("p < 0.01", "p < 0.05", "n.s."))
names(intPval) <- res.int.bl$Variable
intPval <- intPval[rownames(logOrMat)]
intPvalNum <- res.int.bl$P_Interaction
names(intPvalNum) <- res.int.bl$Variable
intPvalNum <- intPvalNum[rownames(logOrMat)]
intPvalNum <- round(intPvalNum, 2)

ra <- rowAnnotation("Interaction (Arm)"=intPval, 
         col=list("Interaction (Arm)"=structure(c(carto_pal(7, "SunsetDark")[c(6, 3)], "grey80"), names=levels(intPval))), 
         annotation_legend_param=list(title="Interaction test"))

textra <- rowAnnotation("P interaction" = anno_text(intPvalNum, show_name = TRUE, gp=gpar(fontsize=12, fontface=ifelse(intPvalNum<0.05, "bold", "plain"))), 
           annotation_label=expression(P[interaction]))

fdata <- fData(scoreSS)
fdata <- fdata[rownames(logOrMat), ]
fdata$Annotation <- droplevels(fdata$Annotation)
la <- rowAnnotation(Annotation=fdata$Annotation, 
         col=list(Annotation=annotCol), 
         annotation_legend_param = list(Annotation = list(at=levels(fdata$Annotation)[c(6, 8, 7, 4, 3, 5, 2, 1)])))

orCol <- colorRamp2(c(min(logOrMat), 0, max(logOrMat)), brewer.pal(11, "RdBu")[c(10, 6, 2)])
showtext.auto()
hm.bl <- Heatmap(logOrMat, 
        name="OR", 
        col=orCol, 
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
        row_names_gp = gpar(fontsize=11), 
        column_names_gp = gpar(fontsize=13), 
        column_labels = c("CT", "CT/A"), 
        column_title = "Baseline", 
        rect_gp = gpar(col = "grey40", lwd = 0.6), 
        right_annotation = textra, 
        left_annotation = la, 
        width = unit(0.55, "in"), 
        #height = unit(6, "in"), 
        cell_fun = function(j, i, x, y, width, height, fill) {
         grid.text(fdrMat2[i, j], x, y, gp = gpar(fontsize = 9, fontface = "bold"))}
)
lgd_list = list(
 Legend(labels = c("*: FDR < 0.1", "**: FDR < 0.05"), title = "Gene set significance", type = "grid"))
g1 = grid.grabExpr(draw(hm.bl, annotation_legend_list = lgd_list, merge_legend = TRUE, legend_grouping = "original"))

#*******************************
# Panel B
#*******************************

res.int.bl2 <- res.int.bl[match(rownames(logOrMat), res.int.bl$Variable), ]
intSig <- res.int.bl2$Variable[round(res.int.bl2$P_Interaction, 2)<0.05]
df <- data.frame(GeneSet=rep(intSig, each=ncol(scoreSS)), 
        SingScore=c(exprs(scoreSS)[intSig[1], ], exprs(scoreSS)[intSig[2], ], exprs(scoreSS)[intSig[3], ], exprs(scoreSS)[intSig[4], ], exprs(scoreSS)[intSig[5], ], exprs(scoreSS)[intSig[6], ], exprs(scoreSS)[intSig[7], ]), 
        Arm=rep(scoreSS$Arm, 7), 
        pCR=rep(scoreSS$pCR, 7))
df$GeneSet <- factor(df$GeneSet, levels = unique(df$GeneSet)[c(5, 2, 4, 1, 6, 7, 3)])
df <- df[order(df$GeneSet, df$Arm, df$SingScore), ]
df$SingScore_cat <- unlist(tapply(df$SingScore, df$GeneSet, function(x)ifelse(x>=median(x), "High", "Low")))
df$Arm <- factor(df$Arm, levels=c("CT", "CT/A"))
df$pCR.num <- ifelse(df$pCR  == "pCR", 1, 0)
df$pCR <- factor(df$pCR, levels=c("RD", "pCR"))
df$Group <- paste(df$SingScore_cat, df$Arm)
df$Group <- factor(df$Group, levels=c("Low CT", "Low CT/A", "High CT", "High CT/A"))


test <- data.frame(GeneSet=intSig, pvalue=0)
test$pvalue[1] <- summary(glm(pCR.num ~ SingScore_cat*Arm, data=df[df$GeneSet  == test$GeneSet[1], ]))$coefficients[4, 4]
test$pvalue[2] <- summary(glm(pCR.num ~ SingScore_cat*Arm, data=df[df$GeneSet  == test$GeneSet[2], ]))$coefficients[4, 4]
test$pvalue[3] <- summary(glm(pCR.num ~ SingScore_cat*Arm, data=df[df$GeneSet  == test$GeneSet[3], ]))$coefficients[4, 4]
test$pvalue[4] <- summary(glm(pCR.num ~ SingScore_cat*Arm, data=df[df$GeneSet  == test$GeneSet[4], ]))$coefficients[4, 4]
test$pvalue[5] <- summary(glm(pCR.num ~ SingScore_cat*Arm, data=df[df$GeneSet  == test$GeneSet[5], ]))$coefficients[4, 4]
test$pvalue[6] <- summary(glm(pCR.num ~ SingScore_cat*Arm, data=df[df$GeneSet  == test$GeneSet[6], ]))$coefficients[4, 4]
test$pvalue <- round(test$pvalue, 3)
test$Label <- paste("italic(P)[interaction]  == ", test$pvalue)

t2x2 <- table(df$Group[df$GeneSet  == intSig[1]], df$pCR[df$GeneSet  == intSig[1]])
perc <- round(t2x2/rowSums(t2x2)*100, 1)
df.perc <- data.frame(Perc=perc[, 2], 
          N=t2x2[, 2], 
          Variable=gsub(" .*", "", rownames(t2x2)), 
          Arm=gsub(".* ", "", rownames(t2x2)), 
          Group=rownames(t2x2))
df.perc$Arm <- factor(df.perc$Arm, levels=c("CT", "CT/A"))
df.perc$xLab <- paste0(df.perc$Variable, "\nN=", rep(c(sum(t2x2[1:2, ]), sum(t2x2[3:4, ])), each=2))
df.perc$Group <- factor(df.perc$Group, levels=df.perc$Group)
df.perc$Position <- c(1, 2, 4, 5)
df.perc$GeneSet <- intSig[1]
for(i in 2:length(intSig)){
 t2x2 <- table(df$Group[df$GeneSet  == intSig[i]], df$pCR[df$GeneSet  == intSig[i]])
 perc <- round(t2x2/rowSums(t2x2)*100, 1)
 tmp <- data.frame(Perc=perc[, 2], 
         N=t2x2[, 2], 
         Variable=gsub(" .*", "", rownames(t2x2)), 
         Arm=gsub(".* ", "", rownames(t2x2)), 
         Group=rownames(t2x2))
 tmp$Arm <- factor(tmp$Arm, levels=c("CT", "CT/A"))
 tmp$xLab <- paste0(tmp$Variable, "\nN=", rep(c(sum(t2x2[1:2, ]), sum(t2x2[3:4, ])), each=2))
 tmp$Group <- factor(tmp$Group, levels=tmp$Group)
 tmp$Position <- c(1, 2, 4, 5)
 tmp$GeneSet <- intSig[i]
 df.perc <- rbind(df.perc, tmp)
}
df.perc <- merge(df.perc, test, by="GeneSet")
df.perc$GeneSet <- factor(df.perc$GeneSet, levels=levels(df$GeneSet))

g2 <- ggplot(data=df.perc, aes(x=Position, y=Perc, fill=Arm)) +
 geom_bar(stat="identity", width=0.99, position = position_dodge(width=0.9)) +
 facet_wrap(~GeneSet, ncol=2) +
 xlab("") +
 ylab("pCR rate (%)") +
 scale_x_continuous(breaks=c(1.5, 4.5), labels=unique(df.perc$xLab)) +
 scale_y_continuous(breaks=seq(0, 100, 20), limits = c(0, 100)) +
 geom_text(aes(x=Position, y=Perc+6, label=paste0(round(Perc, 1), "")), color="black", size=3) +
 scale_fill_manual(name="", values = carto_pal(12, "Safe")[1:2]) +
 theme_pubr(border=T) +
 theme(plot.margin = margin(l=-0.7, r=0.1, t=0, b=0, unit="in"), 
    axis.title = element_text(size=12), 
    axis.text.x=element_text(size=12), 
    axis.text.y=element_text(size=10), 
    legend.text = element_text(size=12), 
    title = element_text(size=16), 
    legend.key.size = unit(0.5, "cm"), 
    legend.position="top") +
 guides(fill=guide_legend(ncol=2)) 

#*******************************
# Panel assembly
#*******************************

combined <- plot_grid(g1, g2, ncol = 2, rel_widths = c(1, 0.4), labels = c("A", "B"), label_size = 16, label_x = c(0, -0.3))

title <- ggdraw() +
  draw_label("Figure 2",
             x = 0, hjust = 0,
             fontface = "bold", size = 16)
final_plot <- plot_grid(title, combined,
                        ncol = 1,
                        rel_heights = c(0.1, 1)) 
ggsave(file.path(here("results", "figures"), "Figure2.eps"),
       plot = final_plot,
       device = cairo_ps,
       width = 12, height = 8.5)

