

library(Biobase)
library(limma)
library(rcartocolor)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(RColorBrewer)
library(showtext)

#***************************************************
# differential expression between arms at baseline
#***************************************************

data_dir <- here::here("data", "processed")

load(file.path(data_dir, "NeoTRIP_singscore_BL_D1C2.RData"))
scoreSS <- scoreSS[, scoreSS$Timepoint == "Baseline"]
f <- factor(gsub("\\/", "", scoreSS$Arm))
design<-model.matrix(~0+f)
fit <- lmFit(scoreSS, design)
contrast.matrix <- makeContrasts(fCT - fCTA,   
                                 levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
res <- topTable(fit2, number = nrow(scoreSS), adjust="BH")
res <- res[order(res$logFC, decreasing = T),]
res$GeneSet <- factor(res$GeneSet, levels = res$GeneSet)

annotCol <- structure(c(carto_pal(12, "Vivid")[1:9], "grey75"), names=names(table(fData(scoreSS)$Annotation)))
res$Annotation <- factor(res$Annotation, levels = unique(res$Annotation))
annotCol <- annotCol[levels(res$Annotation)]
g1 <- ggplot(data = res, aes(x = GeneSet, y = logFC, fill = Annotation)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(name = "Annotation", values = annotCol) +
  theme_pubr(legend = "right") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
  xlab("") +
  ylab("\U0394 score CT vs CT/A")

res<-res[,c(1,2,3,5,6,7)]
write.table(res,file=file.path(here::here("results", "tables"), "05-DE_CTvsCTA_BL.txt"), sep = "\t", row.names = F, quote = F)

#*******************************
#* Panel assembly
#*******************************

showtext.auto()
title <- ggdraw() +
  draw_label("Supplementary Figure 2",
             x = 0, hjust = 0,
             fontface = "bold", size = 16)
final_plot <- plot_grid(title, g1,
                        ncol = 1,
                        rel_heights = c(0.1, 1)) 
ggsave(file.path(here("results", "figures"), "SuppFigure2.eps"),
       plot = final_plot,
       device = cairo_ps,
       width = 15, height = 5.5)




