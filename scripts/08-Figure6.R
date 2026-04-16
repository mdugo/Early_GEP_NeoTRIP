

library(xgboost)
library(caret)
library(ggplot2)
library(ggpubr)
library(Biobase)
library(SHAPforxgboost)
library(pdp)
library(pROC)
library(dplyr)
library(doParallel)
library(foreach)
library(openxlsx)
library(cowplot)
library(rcartocolor)
library(colorspace)

data_dir <- here::here("data", "processed")

load(file.path(data_dir, "NeoTRIP_singscore_BL_D1C2.RData"))

bl.ct <- scoreSS[, scoreSS$Timepoint == "Baseline" & scoreSS$Arm == "CT"]
bl.cta <- scoreSS[, scoreSS$Timepoint == "Baseline" & scoreSS$Arm == "CT/A"]
on.ct <- scoreSS[, scoreSS$Timepoint == "D1C2" & scoreSS$Arm == "CT"]
on.cta <- scoreSS[, scoreSS$Timepoint == "D1C2" & scoreSS$Arm == "CT/A"]

paired <- names(which(table(scoreSS$Patient_ID) == 2))
scoreSSp <- scoreSS[, scoreSS$Patient_ID %in% paired]
scoreSSp <- scoreSSp[, order(scoreSSp$Arm, scoreSSp$Patient_ID)]
colnames(scoreSSp) <- scoreSSp$Sample_ID

# expression set for combined timepoints

bltmp <- exprs(scoreSSp)[, scoreSSp$Timepoint == "Baseline"]
rownames(bltmp) <- paste0("BL_", rownames(bltmp))
ontmp <- exprs(scoreSSp)[, scoreSSp$Timepoint == "D1C2"]
rownames(ontmp) <- paste0("D1C2_", rownames(ontmp))
both <- ExpressionSet(assayData  =  rbind(bltmp, ontmp), 
                    phenoData  =  new("AnnotatedDataFrame", pData(scoreSSp)[scoreSSp$Timepoint == "Baseline", ]))

##### WARNING: REPEATED CV TAKES A LONG TIME TO RUN #####
# 
# expList <- list(bl.ct, bl.cta, on.ct, on.cta, both[, both$Arm == "CT"], both[, both$Arm == "CT/A"])
# names(expList) <- c("Baseline CT", "Baseline CT/A", "D1C2 CT", "D1C2 CT/A", "Both CT", "Both CT/A")
# 
# resList <- vector("list", lengt = length(expList))
# names(resList) <- names(expList)
# 
# best_params_list <- vector("list", lengt = length(expList))
# names(best_params_list) <- names(expList)
# 
# # --- GLOBAL PARAMETERS ---
# n_repeats <- 100
# n_folds <- 10
# seed_base <- 123
# nrounds <- 200
# 
# # --- GRID SEARCH ---
# 
# param_grid <- expand.grid(
#   nrounds  =  c(100, 200, 300),
#   max_depth  =  c(3, 6, 9),
#   eta  =  c(0.01, 0.1, 0.3),
#   gamma  =  c(0, 1),
#   colsample_bytree  =  c(0.6, 0.8, 1),
#   min_child_weight  =  c(1, 5),
#   subsample  =  c(0.6, 0.8, 1)
# )
# 
# 
# # --- EXTERNAL LOOP ITERATING OVER DATASETS ---
# for (j in 1:length(expList)) {
#   cat("\n Dataset:", names(expList)[j], "\n")
# 
#   data <- t(exprs(expList[[j]]))
#   label <- expList[[j]]$pCR.num
# 
#   # --- GRID SEARCH FOR HYPERPARAMETERS OPTIMIZATION ---
#   best_auc <- 0
#   best_params <- NULL
# 
#   cat("Hyperparameters optimization...\n")
#   for (i in 1:nrow(param_grid)) {
#     params <- param_grid[i, ]
# 
#     xgb_cv <- tryCatch({
#       xgb.cv(
#         data  =  data,
#         label  =  label,
#         nfold  =  5,
#         nrounds  =  nrounds,
#         metrics  =  "auc",
#         objective  =  "binary:logistic",
#         eta  =  params$eta,
#         max_depth  =  params$max_depth,
#         subsample  =  params$subsample,
#         colsample_bytree  =  params$colsample_bytree,
#         verbose  =  0,
#         stratified  =  TRUE,
#         early_stopping_rounds  =  20
#       )
#     }, error  =  function(e) NULL)
# 
#     if (!is.null(xgb_cv)) {
#       mean_auc <- max(xgb_cv$evaluation_log$test_auc_mean)
#       if (mean_auc > best_auc) {
#         best_auc <- mean_auc
#         best_params <- params
#       }
#     }
#   }
# 
#   if (is.null(best_params)) {
#     warning("No valid model for ", names(expList)[j])
#     next
#   }
# 
#   cat("Optimal hyperparameters found for", names(expList)[j], ":\n")
#   print(best_params)
#   cat("Mean AUC:", round(best_auc, 4), "\n")
# 
#   best_params_list[[names(expList)[j]]] <- best_params
# 
#   # --- PARALLEL REPEATED CV ---
#   n_cores <- parallel::detectCores() - 1
#   cl <- makeCluster(n_cores)
#   registerDoParallel(cl)
#   cat("Using", n_cores, "cores in parallel\n")
# 
#   results_list <- foreach(
#     iteration  =  1:n_repeats,
#     .packages  =  c("xgboost", "dplyr"),
#     .export  =  c("data", "label", "n_folds", "nrounds", "seed_base", "best_params")
#   ) %dopar% {
#     set.seed(seed_base + iteration)
# 
#     tryCatch({
#       xgb_cv <- xgb.cv(
#         data  =  data,
#         label  =  label,
#         nfold  =  n_folds,
#         nrounds  =  nrounds,
#         metrics  =  c("auc", "error"),
#         objective  =  "binary:logistic",
#         eta  =  best_params$eta,
#         max_depth  =  best_params$max_depth,
#         subsample  =  best_params$subsample,
#         colsample_bytree  =  best_params$colsample_bytree,
#         verbose  =  0,
#         stratified  =  TRUE,
#         early_stopping_rounds  =  20
#       )
# 
#       best_iter <- which.max(xgb_cv$evaluation_log$test_auc_mean)
#       best_auc <- xgb_cv$evaluation_log$test_auc_mean[best_iter]
#       best_acc <- 1 - xgb_cv$evaluation_log$test_error_mean[best_iter]
# 
#       # Final model to extract feature importance
#       xgb_model <- xgboost(
#         data  =  data,
#         label  =  label,
#         nrounds  =  best_iter,
#         eta  =  best_params$eta,
#         max_depth  =  best_params$max_depth,
#         subsample  =  best_params$subsample,
#         colsample_bytree  =  best_params$colsample_bytree,
#         objective  =  "binary:logistic",
#         verbose  =  0,
#         nthread  =  1
#       )
# 
#       importance_matrix <- xgb.importance(feature_names  =  colnames(data), model  =  xgb_model)
# 
#       list(
#         iteration  =  iteration,
#         accuracy  =  best_acc,
#         auc  =  best_auc,
#         importance  =  importance_matrix
#       )
#     }, error  =  function(e) {
#       list(iteration  =  iteration, accuracy  =  NA, auc  =  NA, importance  =  NULL)
#     })
#   }
# 
#   stopCluster(cl)
# 
#   # --- RESULTS LIST ---
#   results <- bind_rows(lapply(results_list, function(x)
#     data.frame(iteration  =  x$iteration, accuracy  =  x$accuracy, auc  =  x$auc)))
#   results$dataset <- names(expList)[j]
# 
#   feature_importance_list <- lapply(results_list, function(x) x$importance)
#   feature_importance_list <- feature_importance_list[!sapply(feature_importance_list, is.null)]
# 
#   summary_results <- results %>%
#     filter(!is.na(accuracy)) %>%
#     summarise(
#       mean_accuracy  =  mean(accuracy),
#       sd_accuracy  =  sd(accuracy),
#       mean_auc  =  mean(auc),
#       sd_auc  =  sd(auc)
#     ) %>%
#     mutate(dataset  =  names(expList)[j])
# 
#   if (length(feature_importance_list) > 0) {
#     importance_df <- bind_rows(feature_importance_list, .id  =  "iteration")
#     importance_df$dataset <- names(expList)[j]
# 
#     mean_importance <- importance_df %>%
#       group_by(Feature) %>%
#       summarise(
#         meanGain  =  mean(Gain, na.rm  =  TRUE),
#         stdvGain  =  sd(Gain, na.rm  =  TRUE),
#         meanCover  =  mean(Cover, na.rm  =  TRUE),
#         meanFrequency  =  mean(Frequency, na.rm  =  TRUE)
#       ) %>%
#       arrange(desc(meanGain)) %>%
#       mutate(dataset  =  names(expList)[j])
#   }
# 
#   # --- SAVE IN LIST ---
#   resList[[j]]$results <- results
#   resList[[j]]$summary_results <- summary_results
#   resList[[j]]$importance_df <- importance_df
#   resList[[j]]$mean_importance <- mean_importance
# }
# save(resList, best_params_list, file = file.path(here::here("results"), "06-XGBoost_repeatedCV_results.RData"))
# 

#**********************************************
#* Panel A
#**********************************************

load(file.path(here::here("results"), "06-XGBoost_repeatedCV_results.RData"))

perfs <- lapply(resList, function(z)as.data.frame(z$results))
df <- Reduce(rbind, perfs)
df$dataset <- factor(df$dataset, levels = unique(df$dataset))
df$Arm <- gsub(".* ", "", df$dataset)
df$Timepoint <- gsub(" .*", "", df$dataset)
df$Timepoint <- factor(df$Timepoint, levels = unique(df$Timepoint))

median_acc <- round(tapply(df$accuracy, df$dataset, median, na.rm = T), 2)
median_acc <- data.frame(dataset = names(median_acc), value = median_acc)
median_acc$Timepoint <- gsub(" .*", "", median_acc$dataset)
median_acc$Arm <- gsub(".* ", "", median_acc$dataset)
median_acc$Timepoint <- factor(median_acc$Timepoint, levels = levels(df$Timepoint))

g0 <- ggplot(data = df, aes(x = Arm, y = accuracy, fill = Arm)) +
  geom_boxplot(coef = NULL, width = 0.55) +
  scale_fill_manual(name = "", values = carto_pal(12, "Safe")[1:2]) +
  geom_jitter(width = 0.15, size = 1) +
  geom_text(data = median_acc, aes(x = Arm, label = value), y = 0.81, size  =  3.5) +
  facet_wrap(~Timepoint) +
  ylim(0.48, 0.82) +
  theme_pubr(border = T, legend = "none", base_size  =  9) +
  theme(panel.grid.major.y  =  element_line(color  =  "grey80", linetype = "dotted"), 
        strip.text.x  =  element_text(face  =  "bold", size  =  10)) +
  xlab("") +
  ylab("Mean accuracy")

#**********************************************
#* Panel B
#**********************************************

annot <- fData(scoreSS)
colnames(annot)[1] <- "Feature"
annotCol <- structure(c(carto_pal(12, "Vivid")[1:9], "grey75"), names = names(table(annot$Annotation)))

imp1 <- resList$`Baseline CT/A`$importance_df
imp1 <- imp1 %>%
  group_by(iteration) %>%
  mutate(rank  =  rank(-Gain)) %>%
  ungroup()
imp1 <- merge(imp1, annot, by = "Feature")
mediane1 <- sort(tapply(imp1$rank, imp1$Feature, median), decreasing = F)[1:10]
imp1 <- imp1[imp1$Feature %in% names(mediane1), ]

imp2 <- resList$`D1C2 CT/A`$importance_df
imp2 <- imp2 %>%
  group_by(iteration) %>%
  mutate(rank  =  rank(-Gain)) %>%
  ungroup()
imp2 <- merge(imp2, annot, by = "Feature")
mediane2 <- sort(tapply(imp2$rank, imp2$Feature, median), decreasing = F)[1:10]
imp2 <- imp2[imp2$Feature %in% names(mediane2), ]

imp <- rbind(imp1, imp2)
imp$Feature <- factor(imp$Feature, levels = unique(c(names(mediane1), names(mediane2))))

g1 <- ggplot(data = imp, aes(x = Feature, y = Gain, fill = Annotation)) +
  geom_boxplot(coef = NULL) +
  scale_fill_manual(values = annotCol) +
  coord_flip() +
  xlab("") +
  ylab("Gain") +
  theme_pubr(legend = "top", border = T, base_size  =  9) +
  theme(
    #axis.text.y  =  element_text(size = 11), 
        strip.text.x  =  element_text(face  =  "bold", size  =  10)) +
  facet_wrap(~dataset, ncol = 2, scales  =  "free") +
  guides(fill = guide_legend(ncol = 3))

#**********************************************
#* Panel C
#**********************************************

dtrain <- xgb.DMatrix(data  =  t(exprs(on.cta)), label  =  on.cta$pCR.num)
xgb_model <- xgb.train(
  params  =  as.list(best_params_list$`D1C2 CT/A`[-1]), 
  data  =  dtrain, 
  nrounds = 200, 
  objective  =  "binary:logistic", 
  verbose  =  1
)


myexp <- data.frame(Iron_utilization = exprs(on.cta)["Iron_utilization", ], cTME_Cytotoxic_cells = exprs(on.cta)["cTME_Cytotoxic_cells", ], pCR = on.cta$pCR, pCR.num = on.cta$pCR.num)
myexp$cTME_Cytotoxic_cells_cat <- ifelse(myexp$cTME_Cytotoxic_cells >= median(myexp$cTME_Cytotoxic_cells), "High", "Low")
myexp$Iron_utilization_cat <- ifelse(myexp$Iron_utilization >= median(myexp$Iron_utilization), "High", "Low")
myexp$Combination <- interaction(myexp$Iron_utilization_cat, myexp$cTME_Cytotoxic_cells_cat, sep = "-")
myexp$Combination <- factor(myexp$Combination, levels = c("High-High", "High-Low", "Low-High", "Low-Low"))
myexp$cTME_Cytotoxic_cells_cat <- factor(myexp$cTME_Cytotoxic_cells_cat, levels = c("Low", "High"))
myexp$Iron_utilization_cat <- factor(myexp$Iron_utilization_cat, levels = c("Low", "High"))

p1 <- pdp::partial(xgb_model, 
                   pred.var  =  c("Iron_utilization", "cTME_Cytotoxic_cells"), 
                   train  =  t(exprs(on.cta)), 
                   type  =  "classification")
sigmoid <- function(x) 1 / (1 + exp(-x))
p1$yhat_prob <- sigmoid(p1$yhat)
p1$yhat_scaled <- (p1$yhat - min(p1$yhat)) / (max(p1$yhat) - min(p1$yhat))


g2 <- ggplot(p1, aes(x  =  Iron_utilization, y  =  cTME_Cytotoxic_cells, fill  =  yhat_scaled)) +
  geom_tile() +
  scale_fill_continuous_sequential(name = "P(pCR)", palette = "Greens", alpha = 0.7) +
  ggnewscale::new_scale_fill() +
  geom_point(data = myexp, aes(x = Iron_utilization, y = cTME_Cytotoxic_cells, fill = pCR), pch = 21, inherit.aes  =  F, size = 2.3) +
  scale_fill_manual(name = "pCR status", values  =  carto_pal(12, "Pastel")[2:1]) +
  theme_pubr(legend = "right", base_size  =  9) +
  labs(title = "D1C2 CT/A") +
  theme(plot.title  =  element_text(face = "bold", size = 9)) +
  geom_vline(xintercept  =  median(p1$Iron_utilization), linetype = "dashed") +
  geom_hline(yintercept  =  median(p1$cTME_Cytotoxic_cells), linetype = "dashed")

#**********************************************
#* Panel D
#**********************************************

df1 <- myexp %>%
  group_by(cTME_Cytotoxic_cells_cat) %>%
  summarise(
    n  =  n(), 
    n_pCR  =  sum(pCR  ==  "pCR"), 
    perc_pCR  =  100 * n_pCR / n, 
    label  =  paste0("\nN = ", n)
  ) %>%
  ungroup()
df1$label <- paste0(df1$cTME_Cytotoxic_cells_cat, df1$label)

df2 <- myexp %>%
  group_by(Iron_utilization_cat) %>%
  summarise(
    n  =  n(), 
    n_pCR  =  sum(pCR  ==  "pCR"), 
    perc_pCR  =  100 * n_pCR / n, 
    label  =  paste0("\nN = ", n)
  ) %>%
  ungroup()
df2$label <- paste0(df2$Iron_utilization_cat, df2$label)

df3 <- myexp %>%
  group_by(cTME_Cytotoxic_cells_cat, Iron_utilization_cat) %>%
  summarise(
    n  =  n(), 
    n_pCR  =  sum(pCR  ==  "pCR"), 
    perc_pCR  =  100 * n_pCR / n, 
    label  =  paste0("\nN = ", n)
  ) %>%
  ungroup()


g3a <- ggplot(df1, 
            aes(x  =  label, 
                y  =  perc_pCR)) +
  geom_bar(stat  =  "identity") +
  ylab("pCR rate (%)") +
  xlab("") +
  theme_pubr(base_size  =  9) +
  theme(panel.grid.major.y  =  element_line(linewidth = 0.6, linetype = "dotted", color = "grey20"), 
        plot.title  =  element_text(size = 9, face = "bold")) +
  ylim(0, 100) +
  labs(title = "cTME_Cytotoxic_cells")
g3b <- ggplot(df2, 
            aes(x  =  label, 
                y  =  perc_pCR)) +
  geom_bar(stat  =  "identity") +
  ylab("pCR rate (%)") +
  xlab("") +
  theme_pubr(base_size  =  9) +
  theme(panel.grid.major.y  =  element_line(linewidth = 0.6, linetype = "dotted", color = "grey20"), 
        plot.title  =  element_text(size = 9, face = "bold")) +
  ylim(0, 100) +
  labs(title = "Iron_utilization")
g3c <- ggplot(df3, 
            aes(x  =  cTME_Cytotoxic_cells_cat, 
                y  =  perc_pCR, fill  =  Iron_utilization_cat)) +
  geom_bar(stat  =  "identity", position = position_dodge()) +
  scale_fill_manual(name = "Iron_utilization\ngroup", values = carto_pal(7, "SunsetDark")[c(1, 3)]) +
  ylab("pCR rate (%)") +
  xlab("cTME_Cytotoxic_cells group") +
  theme_pubr(legend = "right", base_size  =  9) +
  theme(panel.grid.major.y  =  element_line(linewidth = 0.6, linetype = "dotted", color = "grey20"), 
        plot.title  =  element_text(size = 9, face = "bold")) +
  ylim(0, 100) +
  labs(title = "Iron_utilization + cTME_Cytotoxic_cells")
g3ab <- plot_grid(g3a, g3b, ncol = 2)
g3 <- plot_grid(g3ab, g3c, ncol = 1)

#**********************************************
#* Panel E
#**********************************************

mod_A <- glm(pCR.num ~ Iron_utilization_cat, data = myexp, family = "binomial")
mod_B <- glm(pCR.num ~ cTME_Cytotoxic_cells_cat, data = myexp, family = "binomial")
mod_AB <- glm(pCR.num ~ Iron_utilization_cat + cTME_Cytotoxic_cells_cat , data = myexp, family = "binomial")
mod_interaction <- glm(pCR.num ~ Iron_utilization_cat * cTME_Cytotoxic_cells_cat, data  =  myexp, family  =  "binomial")

extract_OR <- function(model) {
  est <- coef(summary(model))
  
  if (is.vector(est)) {
    est <- t(as.matrix(est))
    colnames(est) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    rownames(est) <- names(coef(model))
  }
  
  est <- est[rownames(est) !=  "(Intercept)", , drop = FALSE]
  
  OR <- exp(est[, "Estimate"])
  CI_lower <- exp(est[, "Estimate"] - 1.96 * est[, "Std. Error"])
  CI_upper <- exp(est[, "Estimate"] + 1.96 * est[, "Std. Error"])
  p_val <- est[, "Pr(>|z|)"]
  aic_val <- AIC(model)
  
  data.frame(
    Variable  =  rownames(est), 
    OR  =  round(OR, 2), 
    CI_lower  =  round(CI_lower, 2), 
    CI_upper  =  round(CI_upper, 2), 
    p_value  =  signif(p_val, 3), 
    AIC  =  round(aic_val, 1), 
    row.names  =  NULL
  )
}

# Generate a table for each model

tab_A  <- extract_OR(mod_A)  %>% mutate(Model = "Univariable")
tab_B  <- extract_OR(mod_B)  %>% mutate(Model = "Univariable")
tab_AB <- extract_OR(mod_AB) %>% mutate(Model = "Bivariable")
tab_int <-  extract_OR(mod_interaction) %>% mutate(Model = "Bivariable + interaction")

# Combine all tables

final_table <- bind_rows(tab_A, tab_B, tab_AB, tab_int) %>%
  dplyr::select(Model, Variable, OR, CI_lower, CI_upper, p_value, AIC)
final_table <- final_table[!(final_table$Model == "Bivariable + interaction" & final_table$Variable %in% c("cTME_Cytotoxic_cells_catHigh", "Iron_utilization_catHigh")), ]
final_table$Variable[final_table$Model == "Bivariable + interaction"] <- "Interaction"
final_table <- final_table %>%
  mutate(
    Variable  =  factor(Variable, levels = rev(unique(Variable))), 
    Model  =  factor(Model, levels = unique(Model))
  )
final_table$p_value_lab = paste("P  = ", round(final_table$p_value, 4))
final_table$xpos <- 0.1
final_table$ypos <- c(1, 2, 1.5, 1.5, 1)
final_table$Variable <- gsub("_catHigh", "", final_table$Variable)
final_table$Model2 <- paste0(final_table$Model, "\n(AIC  =  ", final_table$AIC, ")")
final_table$Model2 <- factor(final_table$Model2, levels = unique(final_table$Model2))
aic_table <- unique(final_table[, c(1, 7, 9, 10)])
aic_table$AIC_label <- paste("AIC  == ", aic_table$AIC)

# Model comparison by LRT test
lrt_AB_A <- anova(mod_A, mod_AB, test = "LRT")
lrt_AB_B <- anova(mod_B, mod_AB, test = "LRT")
lrt_inter <- anova(mod_AB, mod_interaction, test = "LRT")
extract_LRT <- function(lrt) {
  data.frame(
    p_value   =  format(lrt$`Pr(>Chi)`[2], scientific  =  T, digits  =  3)
  )
}

# LRT table
lrt_table <- rbind(
  cbind(Comparison = "Iron_utilization vs\nBivariable", extract_LRT(lrt_AB_A)), 
  cbind(Comparison = "cTME_Cytotoxic_cells vs\nBivariable", extract_LRT(lrt_AB_B)), 
  cbind(Comparison = "Bivariable vs\nBivariable + interaction", extract_LRT(lrt_inter))
)

# AIC
lrt_table$delta_AIC <- c(
  round(AIC(mod_AB)-AIC(mod_A), 1), 
  round(AIC(mod_AB)-AIC(mod_B), 1), 
  round(AIC(mod_interaction)-AIC(mod_AB), 1)
)
tab_plot <- ggtexttable(lrt_table, rows  =  NULL, theme  =  ttheme("classic", base_size = 8), 
                        cols  =  c("Model comparison", "LRT p-value", "\u0394 AIC"))
g4a <- ggplot(final_table, aes(x = OR, y = Variable, xmin  =  CI_lower, xmax  =  CI_upper)) +
  geom_point(size = 3, position = position_dodge(width = -0.7)) +
  geom_errorbarh(height = 0.2, position = position_dodge(width = -0.7)) +
  geom_vline(xintercept  =  1, linetype = "dashed") +
  scale_x_continuous(trans = 'log10', limits  =  c(0.1, 270)) +
  theme_pubr(border = T, base_size  =  9) +
  xlab("Odds Ratio (log10)") +
  ylab("") +
  theme(panel.grid.major  =  element_line(linewidth = 0.5, linetype = "dotted", color = "grey30"), 
        strip.text.x  =  element_text(face = "bold"), 
        plot.margin  =  margin(t  =  5, r  =  5, l  =  12, b  =  5)) +
  facet_wrap(~Model2, nrow = 4, scales  =  "free_y", strip.position  =  "top")
showtext::showtext_auto()
g4 <- plot_grid(g4a, tab_plot, ncol  =  1, rel_heights  =  c(1, 0.4), align  =  "v")

#**********************************************
#* Panel assembly
#**********************************************

g01 <- plot_grid(g0, g1, ncol = 2, labels  =  c("A", "B"), rel_widths  =  c(0.6, 1))
g34 <- plot_grid(g3, g4a, ncol = 2, labels  =  c("D", "E"), rel_widths   =  c(1, 0.9))
g234 <- plot_grid(g2, g34 , ncol = 2, labels  =  c("C", ""), rel_widths  =  c(0.6, 1))

combined <- plot_grid(g01, g234, ncol = 1, rel_widths  =  c(1, 1))

title <- ggdraw() +
  draw_label("Figure 6",
             x = 0, hjust = 0,
             fontface = "bold", size = 16)
final_plot <- plot_grid(title, combined,
                        ncol = 1,
                        rel_heights = c(0.1, 1)) 
ggsave(file.path(here("results", "figures"), "Figure6.eps"),
       plot = final_plot,
       device = cairo_ps,
       width = 11, height = 7.5)


