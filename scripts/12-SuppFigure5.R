
library(Biobase)
library(InstaPrism)
library(ggplot2)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)
library(rcartocolor)
library(RColorBrewer)
library(splines)
library(cowplot)
library(showtext)
library(nnet)
library(here)
library(GSEABase)
library(singscore)

# InstaPrism on tpm

data_dir <- here::here("data", "processed")

# expData <- read.table(file.path(data_dir, "NeoTRIP_baseline_D1C2_TPM_ComBat_all_samples.txt"),
#                       header = T,
#                       sep = "\t",
#                       as.is = T)
# rownames(expData) <- expData$HGNC_approved_symbol
# expData <- expData[, -c(1:3)]
# ref = InstaPrism_reference('BRCA')
# deconv_tpm = InstaPrism(bulk_Expr = expData, refPhi_cs = ref)
# save(deconv_tpm, file=file.path(data_dir, "instaprism_deconvolution.RData"))



load(file.path(data_dir, "NeoTRIP_singscore_BL_D1C2.RData"))
load(file.path(data_dir, "instaprism_deconvolution.RData"))
estimated_frac = t(deconv_tpm@Post.ini.ct@theta)
colnames(estimated_frac)<-gsub(" ","_",colnames(estimated_frac))
estimated_frac<-estimated_frac[colnames(scoreSS),]
identical(colnames(scoreSS), rownames(estimated_frac))

#**********************************************
#* Panel A
#**********************************************

# correlation with cTME signatures

ctme<-exprs(scoreSS)[grep("cTME",rownames(scoreSS)),scoreSS$Timepoint=="Baseline"]
correl<-cor(t(ctme),estimated_frac[colnames(ctme),])
colnames(correl)<-gsub("_E"," e", colnames(correl))

h1<-Heatmap(correl,
        name="Pearson",
        show_row_dend = F,
        clustering_method_columns = "ward.D2",
        clustering_method_rows = "ward.D2",
        column_title = "Baseline",
        row_names_side = "left",
        row_names_gp = gpar(fontsize=9),
        rect_gp = gpar(col = "grey40", lwd = 0.6),
        width = unit(2.6, "in"),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(round(correl,2)[i, j], x, y, gp = gpar(fontsize = 7, fontface = "bold"))})

ctme<-exprs(scoreSS)[grep("cTME",rownames(scoreSS)),scoreSS$Timepoint=="D1C2"]
correl<-cor(t(ctme),estimated_frac[colnames(ctme),])
colnames(correl)<-gsub("_E"," e", colnames(correl))

h2<-Heatmap(correl,
            name="Pearson2",
            show_row_dend = F,
            show_heatmap_legend = F,
            clustering_method_columns = "ward.D2",
            clustering_method_rows = "ward.D2",
            column_title = "D1C2",
            row_names_side = "left",
            row_names_gp = gpar(fontsize=9),
            rect_gp = gpar(col = "grey40", lwd = 0.6),
            width = unit(2.6, "in"),
            cell_fun = function(j, i, x, y, width, height, fill) {
              grid.text(round(correl,2)[i, j], x, y, gp = gpar(fontsize = 7, fontface = "bold"))})

ghm1 = grid.grabExpr(draw(h1+h2,gap = unit(10, "mm"), merge_legend=TRUE))


#**********************************************
#* Panel B
#**********************************************

df<-reshape2::melt(estimated_frac)
colnames(df)<-c("Sample_ID","cellType","Fraction")
pdata<-pData(scoreSS)[,c("Sample_ID","Arm","pCR","Timepoint")]
df<-merge(pdata,df, by="Sample_ID")
df$Sample_ID<-factor(df$Sample_ID)
df_summary <- df %>%
  group_by(Timepoint, Arm, pCR, cellType) %>%
  summarise(
    mean_frac = mean(Fraction),
    sd_frac   = sd(Fraction),
    .groups = "drop"
  )
df_summary$cellType<-gsub("_"," ",df_summary$cellType)

gbar<-ggplot(df_summary, aes(x = pCR, y = mean_frac, fill = cellType)) +
  geom_bar(stat = "identity", width = 0.8,color = "black") +
  facet_wrap(Timepoint~Arm, nrow=1) +
  scale_fill_manual(name="Cell type", values=carto_pal(12,"Prism")[1:8]) +
  xlab("") +
  ylab("Cell fraction") +
  theme_pubr(border=T, legend="right") +
  theme(plot.margin = margin(t=30, r = 6, b = 15, l = 6))


# association with pCR at baseline

pdata<-pData(scoreSS)[scoreSS$Timepoint=="Baseline",c("Arm","pCR","pCR.num","Timepoint")]
bl<-estimated_frac[rownames(pdata),]*10
pdata<-data.frame(bl, pdata)

res.bl<-data.frame(Variable=colnames(estimated_frac),OR_CT=0,CI2.5_CT=0,CI97.5_CT=0,P_CT=0,FDR_CT=0,Pcub_CT=0,
                   OR_CTA=0,CI2.5_CTA=0,CI97.5_CTA=0,P_CTA=0,FDR_CTA=0,Pcub_CTA=0,P_Interaction=0,FDR_Interaction=0)

for(i in 1:nrow(res.bl)){
  tmp<-pdata[pdata$Arm=="CT",]
  fit<-glm(pCR.num~tmp[,i], family="binomial", data=tmp)
  fitCubic<-glm(pCR.num ~ ns(tmp[,i], df = 3), family = "binomial", data=tmp)
  lrtest<-anova(fit, fitCubic, test = "LRT")
  res.bl$Pcub_CT[i]<-lrtest$`Pr(>Chi)`[2]
  res.bl$OR_CT[i]<-exp(summary(fit)$coef[2,1])
  res.bl$CI2.5_CT[i]<-exp(confint(fit)[2,1])
  res.bl$CI97.5_CT[i]<-exp(confint(fit)[2,2])
  res.bl$P_CT[i]<-summary(fit)$coef[2,4]
  tmp<-pdata[pdata$Arm=="CT/A",]
  fit<-glm(pCR.num~tmp[,i], family="binomial", data=tmp)
  fitCubic<-glm(pCR.num ~ ns(tmp[,i], df = 3), family = "binomial", data=tmp)
  lrtest<-anova(fit, fitCubic, test = "LRT")
  res.bl$Pcub_CTA[i]<-lrtest$`Pr(>Chi)`[2]
  res.bl$OR_CTA[i]<-exp(summary(fit)$coef[2,1])
  res.bl$CI2.5_CTA[i]<-exp(confint(fit)[2,1])
  res.bl$CI97.5_CTA[i]<-exp(confint(fit)[2,2])
  res.bl$P_CTA[i]<-summary(fit)$coef[2,4]
  fit<-glm(pCR.num~pdata[,i]*Arm, family="binomial", data=pdata)
  res.bl$P_Interaction[i]<-summary(fit)$coef[4,4]
}
res.bl$FDR_CT<-p.adjust(res.bl$P_CT,method="BH")
res.bl$FDR_CTA<-p.adjust(res.bl$P_CTA,method="BH")
res.bl$FDR_Interaction<-p.adjust(res.bl$P_Interaction,method="BH")


# association with pCR at D1C2

pdata<-pData(scoreSS)[scoreSS$Timepoint=="D1C2",c("Arm","pCR","pCR.num","Timepoint")]
d1c2<-estimated_frac[rownames(pdata),]*10
pdata<-data.frame(d1c2, pdata)

res.d1c2<-data.frame(Variable=colnames(estimated_frac),OR_CT=0,CI2.5_CT=0,CI97.5_CT=0,P_CT=0,FDR_CT=0,Pcub_CT=0,
                   OR_CTA=0,CI2.5_CTA=0,CI97.5_CTA=0,P_CTA=0,FDR_CTA=0,Pcub_CTA=0,P_Interaction=0,FDR_Interaction=0)

for(i in 1:nrow(res.d1c2)){
  tmp<-pdata[pdata$Arm=="CT",]
  fit<-glm(pCR.num~tmp[,i], family="binomial", data=tmp)
  fitCubic<-glm(pCR.num ~ ns(tmp[,i], df = 3), family = "binomial", data=tmp)
  lrtest<-anova(fit, fitCubic, test = "LRT")
  res.d1c2$Pcub_CT[i]<-lrtest$`Pr(>Chi)`[2]
  res.d1c2$OR_CT[i]<-exp(summary(fit)$coef[2,1])
  res.d1c2$CI2.5_CT[i]<-exp(confint(fit)[2,1])
  res.d1c2$CI97.5_CT[i]<-exp(confint(fit)[2,2])
  res.d1c2$P_CT[i]<-summary(fit)$coef[2,4]
  tmp<-pdata[pdata$Arm=="CT/A",]
  fit<-glm(pCR.num~tmp[,i], family="binomial", data=tmp)
  fitCubic<-glm(pCR.num ~ ns(tmp[,i], df = 3), family = "binomial", data=tmp)
  lrtest<-anova(fit, fitCubic, test = "LRT")
  res.d1c2$Pcub_CTA[i]<-lrtest$`Pr(>Chi)`[2]
  res.d1c2$OR_CTA[i]<-exp(summary(fit)$coef[2,1])
  res.d1c2$CI2.5_CTA[i]<-exp(confint(fit)[2,1])
  res.d1c2$CI97.5_CTA[i]<-exp(confint(fit)[2,2])
  res.d1c2$P_CTA[i]<-summary(fit)$coef[2,4]
  fit<-glm(pCR.num~pdata[,i]*Arm, family="binomial", data=pdata)
  res.d1c2$P_Interaction[i]<-summary(fit)$coef[4,4]
}
res.d1c2$FDR_CT<-p.adjust(res.d1c2$P_CT,method="BH")
res.d1c2$FDR_CTA<-p.adjust(res.d1c2$P_CTA,method="BH")
res.d1c2$FDR_Interaction<-p.adjust(res.d1c2$P_Interaction,method="BH")


# association with pCR for delta

pdata<-pData(scoreSS)[c("Arm","pCR","pCR.num","Timepoint","Patient_ID")]
pdata<-pdata[pdata$Patient_ID%in%names(which(table(pdata$Patient_ID)==2)),]
pdata<-pdata[order(pdata$Patient_ID,pdata$Timepoint),]
ef.p<-estimated_frac[rownames(pdata),]*10
pdata2<-data.frame(ef.p, pdata)
deltaFrac<-matrix(0,nrow = ncol(ef.p),ncol=length(unique(pdata2$Patient_ID)))
rownames(deltaFrac)<-colnames(ef.p)
colnames(deltaFrac)<-unique(pdata2$Patient_ID)
for(i in 1:nrow(deltaFrac)){
  deltaFrac[i,]<-pdata2[pdata2$Timepoint=="D1C2",i] - pdata2[pdata2$Timepoint=="Baseline",i]
}

pdata<-data.frame(t(deltaFrac), pdata[pdata$Timepoint=="Baseline",])

res.delta<-data.frame(Variable=colnames(estimated_frac),OR_CT=0,CI2.5_CT=0,CI97.5_CT=0,P_CT=0,FDR_CT=0,Pcub_CT=0,
                     OR_CTA=0,CI2.5_CTA=0,CI97.5_CTA=0,P_CTA=0,FDR_CTA=0,Pcub_CTA=0,P_Interaction=0,FDR_Interaction=0)

for(i in 1:nrow(res.delta)){
  tmp<-pdata[pdata$Arm=="CT",]
  fit<-glm(pCR.num~tmp[,i], family="binomial", data=tmp)
  fitCubic<-glm(pCR.num ~ ns(tmp[,i], df = 3), family = "binomial", data=tmp)
  lrtest<-anova(fit, fitCubic, test = "LRT")
  res.delta$Pcub_CT[i]<-lrtest$`Pr(>Chi)`[2]
  res.delta$OR_CT[i]<-exp(summary(fit)$coef[2,1])
  res.delta$CI2.5_CT[i]<-exp(confint(fit)[2,1])
  res.delta$CI97.5_CT[i]<-exp(confint(fit)[2,2])
  res.delta$P_CT[i]<-summary(fit)$coef[2,4]
  tmp<-pdata[pdata$Arm=="CT/A",]
  fit<-glm(pCR.num~tmp[,i], family="binomial", data=tmp)
  fitCubic<-glm(pCR.num ~ ns(tmp[,i], df = 3), family = "binomial", data=tmp)
  lrtest<-anova(fit, fitCubic, test = "LRT")
  res.delta$Pcub_CTA[i]<-lrtest$`Pr(>Chi)`[2]
  res.delta$OR_CTA[i]<-exp(summary(fit)$coef[2,1])
  res.delta$CI2.5_CTA[i]<-exp(confint(fit)[2,1])
  res.delta$CI97.5_CTA[i]<-exp(confint(fit)[2,2])
  res.delta$P_CTA[i]<-summary(fit)$coef[2,4]
  fit<-glm(pCR.num~pdata[,i]*Arm, family="binomial", data=pdata)
  res.delta$P_Interaction[i]<-summary(fit)$coef[4,4]
}
res.delta$FDR_CT<-p.adjust(res.delta$P_CT,method="BH")
res.delta$FDR_CTA<-p.adjust(res.delta$P_CTA,method="BH")
res.delta$FDR_Interaction<-p.adjust(res.delta$P_Interaction,method="BH")



orMat.bl<-res.bl[,grep("OR_",colnames(res.bl))]
rownames(orMat.bl)<-res.bl$Variable
colnames(orMat.bl)<-c("BL_CT","BL_CTA")
fdrMat.bl<-res.bl[,grep("FDR_",colnames(res.bl))]
rownames(fdrMat.bl)<-res.bl$Variable
colnames(fdrMat.bl)<-c("BL_CT","BL_CTA")
fdrMat2.bl<-ifelse(fdrMat.bl<0.05,"\U2217\U2217",ifelse(fdrMat.bl<0.1,"\U2217",""))
logOrMat.bl<-log10(orMat.bl)
intPvalNum.bl<-res.bl$P_Interaction
names(intPvalNum.bl)<-res.bl$Variable
intPvalNum.bl<-round(intPvalNum.bl,2)
textra<-rowAnnotation("P interaction" = anno_text(intPvalNum.bl, show_name = TRUE,gp=gpar(fontsize=12, fontface=ifelse(intPvalNum.bl<0.05,"bold","plain"))),
                      annotation_label=expression(P[interaction]))
orCol<-colorRamp2(c(min(logOrMat.bl), 0, max(logOrMat.bl)), brewer.pal(11,"RdBu")[c(10,6,2)])
showtext.auto()
hm.bl<-Heatmap(logOrMat.bl,
              name="OR",
              col=orCol,
              heatmap_legend_param = list(title = expression("log"[10](OR)),
                                          title_gp = gpar(fontsize = 10, fontface = "bold")),
              cluster_columns = F,
              cluster_rows = T,
              clustering_method_rows = "ward.D2",
              show_row_dend = F,
              row_names_side = "left",
              row_names_gp = gpar(fontsize=11),
              column_names_gp = gpar(fontsize=13),
              column_labels = c("CT","CT/A"),
              column_title = "Baseline",
              rect_gp = gpar(col = "grey40", lwd = 0.6),
              right_annotation = textra,
              width = unit(0.5, "in"),
              #height = unit(6, "in"),
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(fdrMat2.bl[i, j], x, y, gp = gpar(fontsize = 9, fontface = "bold"))}
)
lgd_list = list(
  Legend(labels = c("*: FDR < 0.1", "**: FDR < 0.05"), title = "Gene set significance", type = "grid"))

orMat.d1c2<-res.d1c2[,grep("OR_",colnames(res.d1c2))]
rownames(orMat.d1c2)<-res.d1c2$Variable
colnames(orMat.d1c2)<-c("D1C2_CT","D1C2_CTA")
fdrMat.d1c2<-res.d1c2[,grep("FDR_",colnames(res.d1c2))]
rownames(fdrMat.d1c2)<-res.d1c2$Variable
colnames(fdrMat.d1c2)<-c("D1C2_CT","D1C2_CTA")
fdrMat2.d1c2<-ifelse(fdrMat.d1c2<0.05,"\U2217\U2217",ifelse(fdrMat.d1c2<0.1,"\U2217",""))
logOrMat.d1c2<-log10(orMat.d1c2)
intPvalNum.d1c2<-res.d1c2$P_Interaction
names(intPvalNum.d1c2)<-res.d1c2$Variable
intPvalNum.d1c2<-round(intPvalNum.d1c2,2)
textra<-rowAnnotation("P interaction" = anno_text(intPvalNum.d1c2, show_name = TRUE,gp=gpar(fontsize=12, fontface=ifelse(intPvalNum.d1c2<0.05,"bold","plain"))),
                      annotation_label=expression(P[interaction]))
orCol<-colorRamp2(c(min(logOrMat.d1c2), 0, max(logOrMat.d1c2)), brewer.pal(11,"RdBu")[c(10,6,2)])
showtext.auto()
hm.d1c2<-Heatmap(logOrMat.d1c2,
               name="OR2",
               col=orCol,
               show_heatmap_legend = F,
               heatmap_legend_param = list(title = expression("log"[10](OR)),
                                           title_gp = gpar(fontsize = 10, fontface = "bold")),
               cluster_columns = F,
               cluster_rows = T,
               clustering_method_rows = "ward.D2",
               show_row_dend = F,
               row_names_side = "left",
               row_names_gp = gpar(fontsize=11),
               column_names_gp = gpar(fontsize=13),
               column_labels = c("CT","CT/A"),
               column_title = "D1C2",
               rect_gp = gpar(col = "grey40", lwd = 0.6),
               right_annotation = textra,
               width = unit(0.5, "in"),
               #height = unit(6, "in"),
               cell_fun = function(j, i, x, y, width, height, fill) {
                 grid.text(fdrMat2.d1c2[i, j], x, y, gp = gpar(fontsize = 9, fontface = "bold"))}
)

orMat.delta<-res.delta[,grep("OR_",colnames(res.delta))]
rownames(orMat.delta)<-res.delta$Variable
colnames(orMat.delta)<-c("BL_CT","BL_CTA")
fdrMat.delta<-res.delta[,grep("FDR_",colnames(res.delta))]
rownames(fdrMat.delta)<-res.delta$Variable
colnames(fdrMat.delta)<-c("BL_CT","BL_CTA")
fdrMat2.delta<-ifelse(fdrMat.delta<0.05,"\U2217\U2217",ifelse(fdrMat.delta<0.1,"\U2217",""))
logOrMat.delta<-log10(orMat.delta)
intPvalNum.delta<-res.delta$P_Interaction
names(intPvalNum.delta)<-res.delta$Variable
intPvalNum.delta<-round(intPvalNum.delta,2)
textra<-rowAnnotation("P interaction" = anno_text(intPvalNum.delta, show_name = TRUE,gp=gpar(fontsize=12, fontface=ifelse(intPvalNum.delta<0.05,"bold","plain"))),
                      annotation_label=expression(P[interaction]))
orCol<-colorRamp2(c(min(logOrMat.delta), 0, max(logOrMat.delta)), brewer.pal(11,"RdBu")[c(10,6,2)])
showtext.auto()
hm.delta<-Heatmap(logOrMat.delta,
               name="OR3",
               col=orCol,
               show_heatmap_legend = F,
               heatmap_legend_param = list(title = expression("log"[10](OR)),
                                           title_gp = gpar(fontsize = 10, fontface = "bold")),
               cluster_columns = F,
               cluster_rows = T,
               clustering_method_rows = "ward.D2",
               show_row_dend = F,
               row_names_side = "left",
               row_names_gp = gpar(fontsize=11),
               column_names_gp = gpar(fontsize=13),
               column_labels = c("CT","CT/A"),
               column_title = "\u0394",
               rect_gp = gpar(col = "grey40", lwd = 0.6),
               right_annotation = textra,
               width = unit(0.5, "in"),
               #height = unit(6, "in"),
               cell_fun = function(j, i, x, y, width, height, fill) {
                 grid.text(fdrMat2.delta[i, j], x, y, gp = gpar(fontsize = 9, fontface = "bold"))}
)
ghm2 = grid.grabExpr(draw(hm.bl+hm.d1c2+hm.delta,gap = unit(10, "mm"),annotation_legend_list = lgd_list, merge_legend=TRUE))

#************************************
# Panel D
#************************************

pdata<-pData(scoreSS)[scoreSS$Timepoint=="Baseline" & !is.na(scoreSS$Early.pCR),]
bl<-estimated_frac[rownames(pdata),]*10
pdata<-data.frame(bl, pdata)
pdata$Early.pCR <- gsub("e-pCR", "early-pCR", pdata$Early.pCR)
pdata$Early.pCR <- factor(pdata$Early.pCR, levels = c("Sustained early-pCR", "Transient early-pCR", "pCR", "RD"))

res.epcr<-data.frame(Variable=colnames(pdata)[1:8],OR_CT=0,P_CT=0,OR_CTA=0,P_CTA=0,ORR_Interaction=0,P_Interaction=0)

for(i in 1:8){
  df<-data.frame(pCR=pdata$Early.pCR,Arm=pdata$Arm, Variable=pdata[,i])

  # CT
  model<-multinom(pCR ~ Variable, data = df[df$Arm=="CT",], Hess = TRUE)
  z <- summary(model)$coefficients/summary(model)$standard.errors
  p <- (1 - pnorm(abs(z), 0, 1)) * 2
  res.epcr$OR_CT[i]<-exp(summary(model)$coefficients[2,2])
  res.epcr$P_CT[i]<-p[2,2]
  
  # CT/A
  model<-multinom(pCR ~ Variable, data = df[df$Arm=="CT/A",], Hess = TRUE)
  z <- summary(model)$coefficients/summary(model)$standard.errors
  p <- (1 - pnorm(abs(z), 0, 1)) * 2
  res.epcr$OR_CTA[i]<-exp(summary(model)$coefficients[2,2])
  res.epcr$P_CTA[i]<-p[2,2]
  
  # Interaction
  model<-multinom(pCR ~ Variable*Arm, data = df, Hess = TRUE)
  z <- summary(model)$coefficients/summary(model)$standard.errors
  p <- (1 - pnorm(abs(z), 0, 1)) * 2
  res.epcr$ORR_Interaction[i]<-exp(summary(model)$coefficients[2,4])
  res.epcr$P_Interaction[i]<-p[2,4]
}
res.epcr$FDR_CT<-p.adjust(res.epcr$P_CT,method="BH")
res.epcr$FDR_CTA<-p.adjust(res.epcr$P_CTA,method="BH")
res.epcr<-res.epcr[,c(1,2,3,8,4,5,9,6,7)]


nsamp<-data.frame(table(pdata$Early.pCR,pdata$Arm))
colnames(nsamp)<-c("Early.pCR","Arm","N")
rng<-(max(pdata$B.cells/10)-min(pdata$B.cells/10))/20
nsamp$Ypos<-min(pdata$B.cells/10)-rng
nsamp$Label<-paste0("N=",nsamp$N)
signif<-data.frame(Arm=c("CT","CT/A"), Ypos=max(pdata$B.cells/10)+rng,Label=paste0("P=",c(round(res.epcr$P_CT[res.epcr$Variable=="B.cells"],3),round(res.epcr$P_CTA[res.epcr$Variable=="B.cells"],3))))
g1<-ggplot(data=pdata, aes(x=Arm, y=B.cells/10, fill=Early.pCR)) +
  geom_boxplot(coef=NULL, width=0.5, position = position_dodge(width=0.8)) +
  geom_jitter(position = position_jitterdodge(jitter.width=0.15,dodge.width = 0.8),size=0.6) +
  scale_fill_manual(name="", values=carto_pal(12,"Pastel")[c(4,3,2,1)]) +
  xlab("") +
  ylab("Cell fraction") +
  labs(title="B-cells") +
  theme_pubr(border=T, legend="right") +
  theme(panel.grid.major.y =  element_line(colour = "grey80", linewidth = 0.6, linetype="dotted"),
        axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        plot.title = element_text(size=12, face="bold", hjust=0.5),
        legend.text = element_text(size=13)) +
  geom_text(data=nsamp,aes(y=Ypos, label=N),position = position_dodge(width = 0.8),size=4.5) +
  geom_text(data=signif,aes(x=Arm,y=Ypos,label=Label),size=4.1, inherit.aes = F)

nsamp<-data.frame(table(pdata$Early.pCR,pdata$Arm))
colnames(nsamp)<-c("Early.pCR","Arm","N")
rng<-(max(pdata$T.cells/10)-min(pdata$T.cells/10))/20
nsamp$Ypos<-min(pdata$T.cells/10)-rng
nsamp$Label<-paste0("N=",nsamp$N)
signif<-data.frame(Arm=c("CT","CT/A"), Ypos=max(pdata$T.cells/10)+rng,Label=paste0("P=",c(round(res.epcr$P_CT[res.epcr$Variable=="T.cells"],3),round(res.epcr$P_CTA[res.epcr$Variable=="T.cells"],3))))
g2<-ggplot(data=pdata, aes(x=Arm, y=T.cells/10, fill=Early.pCR)) +
  geom_boxplot(coef=NULL, width=0.5, position = position_dodge(width=0.8)) +
  geom_jitter(position = position_jitterdodge(jitter.width=0.15,dodge.width = 0.8),size=0.6) +
  scale_fill_manual(name="", values=carto_pal(12,"Pastel")[c(4,3,2,1)]) +
  xlab("") +
  ylab("Cell fraction") +
  labs(title="T-cells") +
  theme_pubr(border=T, legend="right") +
  theme(panel.grid.major.y =  element_line(colour = "grey80", linewidth = 0.6, linetype="dotted"),
        axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        plot.title = element_text(size=12, face="bold", hjust=0.5),
        legend.text = element_text(size=13)) +
  geom_text(data=nsamp,aes(y=Ypos, label=N),position = position_dodge(width = 0.8),size=4.5) +
  geom_text(data=signif,aes(x=Arm,y=Ypos,label=Label),size=4.1, inherit.aes = F)

nsamp<-data.frame(table(pdata$Early.pCR,pdata$Arm))
colnames(nsamp)<-c("Early.pCR","Arm","N")
rng<-(max(pdata$Cancer_Epithelial/10)-min(pdata$Cancer_Epithelial/10))/20
nsamp$Ypos<-min(pdata$Cancer_Epithelial/10)-rng
nsamp$Label<-paste0("N=",nsamp$N)
signif<-data.frame(Arm=c("CT","CT/A"), Ypos=max(pdata$Cancer_Epithelial/10)+rng,Label=paste0("P=",c(round(res.epcr$P_CT[res.epcr$Variable=="Cancer_Epithelial"],3),round(res.epcr$P_CTA[res.epcr$Variable=="Cancer_Epithelial"],3))))
g3<-ggplot(data=pdata, aes(x=Arm, y=Cancer_Epithelial/10, fill=Early.pCR)) +
  geom_boxplot(coef=NULL, width=0.5, position = position_dodge(width=0.8)) +
  geom_jitter(position = position_jitterdodge(jitter.width=0.15,dodge.width = 0.8),size=0.6) +
  scale_fill_manual(name="", values=carto_pal(12,"Pastel")[c(4,3,2,1)]) +
  xlab("") +
  ylab("Cell fraction") +
  labs(title="Cancer epithelial") +
  theme_pubr(border=T, legend="right") +
  theme(panel.grid.major.y =  element_line(colour = "grey80", linewidth = 0.6, linetype="dotted"),
        axis.text = element_text(size=12),
        axis.title = element_text(size=12),
        plot.title = element_text(size=12, face="bold", hjust=0.5),
        legend.text = element_text(size=13)) +
  geom_text(data=nsamp,aes(y=Ypos, label=N),position = position_dodge(width = 0.8),size=4.5) +
  geom_text(data=signif,aes(x=Arm,y=Ypos,label=Label),size=4.1, inherit.aes = F)

gboxp<-ggarrange(g1, g2, g3, nrow=1, ncol=3, common.legend = T)

topPlot<-plot_grid(ghm1, gbar, nrow=1, labels=c("A","B"), label_size = 20)
bottomPlot<-plot_grid(ghm2, gboxp, nrow=1, labels=c("C","D"), label_size = 20)

combined <- plot_grid(topPlot, bottomPlot, nrow=2, rel_heights = c(1, 0.7), label_size = 20)

title <- ggdraw() +
  draw_label("Supplementary Figure 5",
             x = 0, hjust = 0,
             fontface = "bold", size = 16)
final_plot <- plot_grid(title, combined,
                        ncol = 1,
                        rel_heights = c(0.1, 1)) 
ggsave(file.path(here("results", "figures"), "SuppFigure5.eps"),
       plot = final_plot,
       device = cairo_ps,
       width = 16.5, height = 8.5)

