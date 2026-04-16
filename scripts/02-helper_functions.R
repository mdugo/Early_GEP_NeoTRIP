

#*******************************
# Loop for univariable logistic regression.
# Response must be a numeric vector where 0 indicates failure and 1 success.
# It an also be specified as a factor (when the first level denotes failure
# and all others success).
#*******************************

logistic_regression <- function(response, covariates = NULL, data, rounding_factor = 3, sort = FALSE){
 res.df <- data.frame(Variable = rownames(data), matrix(0, nrow = nrow(data), ncol = (length(covariates)+1)*5))
 colnames(res.df)[-1] <- paste(rep(c("OR", "CI2.5", "CI97.5", "Z", "P"), length(covariates)+1), rep(c("Variable", covariates), each = 5), sep = "_")
 
 Y <- pData(data)[, response]
 covs <- data.frame(pData(data)[, covariates])
 if(is.null(covariates)){
  resTemp <- apply(exprs(data), 1, function(x){
   z <- glm(Y ~ x, family = binomial)
   ee <- c(round(exp(coef(z))[2], rounding_factor), 
      round(exp(confint.default(z, level = 0.95))[2, 1], rounding_factor), 
      round(exp(confint.default(z, level = 0.95))[2, 2], rounding_factor), 
      summary(z)$coefficients[2, 3], 
      summary(z)$coefficients[2, 4])
  })
  res.df[, -1] <- t(resTemp)
  res.df$FDR_Variable <- p.adjust(res.df$P_Variable, method = "BH")
 }
 if(!is.null(covariates)){
  colnames(covs) <- covariates
  resTemp <- apply(exprs(data), 1, function(x){
   df <- data.frame(Y, 
           Gene = x, 
           covs, 
           stringsAsFactors = F)
   f1 <- as.formula(paste("Y ~ Gene +", paste(colnames(df)[3:ncol(df)], collapse = "+")))
   z <- glm(f1, family = binomial, data = df)
   ee <- matrix(as.vector(rbind(round(exp(coef(z))[-1], rounding_factor), 
                 round(exp(confint.default(z, level = 0.95))[-1, 1], rounding_factor), 
                 round(exp(confint.default(z, level = 0.95))[-1, 2], rounding_factor), 
                 summary(z)$coefficients[-1, 3], 
                 summary(z)$coefficients[-1, 4])), nrow = 1)
  })
  res.df[, -1] <- t(resTemp)
  fdr <- apply(res.df[, grep("^P_", colnames(res.df))], 2, p.adjust, method = "BH")
  colnames(fdr) <- gsub("P_", "FDR_", colnames(fdr))
  res.df <- as.data.frame(cbind(res.df, fdr))
  index <- unname(unlist(sapply(c("Variable", covariates), function(x){
   grep(x, colnames(res.df))})))
  res.df <- res.df[, index]
 }
 if(sort){
  res.df <- res.df[order(res.df$P_Variable), ]
 }
 return(res.df)
}

#*******************************
# Loop for univariable logistic regression with interaction test.
# Response must be a numeric vector where 0 indicates failure and 1 success.
# It an also be specified as a factor (when the first level denotes failure
# and all others success).
#*******************************

interaction_logistic_regression <- function(response, covariate_interaction, covariates_adjustment = NULL, data, rounding_factor = 3, sort = FALSE){
  res.df <- data.frame(Variable = rownames(data), matrix(0, nrow = nrow(data), ncol = (length(covariates_adjustment)+3)*5))
  colnames(res.df)[-1] <- paste(rep(c("OR", "CI2.5", "CI97.5", "Z", "P"), length(covariates_adjustment)+3), rep(c("Variable", covariate_interaction, "Interaction", covariates_adjustment), each = 5), sep = "_")
  
  for(i in 1:nrow(data)){
    df <- data.frame(Y = pData(data)[, response], 
                     Gene = exprs(data)[i, ], 
                     Covariate_interaction = pData(data)[, covariate_interaction], 
                     pData(data)[, covariates_adjustment], 
                     stringsAsFactors = F)
    colnames(df)[grep("pData", colnames(df))] <- covariates_adjustment
    if(is.null(covariates_adjustment)){
      
      fit <- glm(Y ~ Gene*Covariate_interaction, family = binomial, data = df)
      res.df[i, grep("_Variable", colnames(res.df))] <- c(round(exp(coef(fit))[2], rounding_factor), 
                                                          round(exp(confint.default(fit, level = 0.95))[2, 1], rounding_factor), 
                                                          round(exp(confint.default(fit, level = 0.95))[2, 2], rounding_factor), 
                                                          summary(fit)$coefficients[2, 3], 
                                                          summary(fit)$coefficients[2, 4])
      res.df[i, grep(covariate_interaction, colnames(res.df))] <- c(round(exp(coef(fit))[3], rounding_factor), 
                                                                    round(exp(confint.default(fit, level = 0.95))[3, 1], rounding_factor), 
                                                                    round(exp(confint.default(fit, level = 0.95))[3, 2], rounding_factor), 
                                                                    summary(fit)$coefficients[3, 3], 
                                                                    summary(fit)$coefficients[3, 4])
      res.df[i, grep("Interaction", colnames(res.df))] <- c(round(exp(coef(fit))[4], rounding_factor), 
                                                            round(exp(confint.default(fit, level = 0.95))[4, 1], rounding_factor), 
                                                            round(exp(confint.default(fit, level = 0.95))[4, 2], rounding_factor), 
                                                            summary(fit)$coefficients[4, 3], 
                                                            summary(fit)$coefficients[4, 4])
    }
    if(!is.null(covariates_adjustment)){
      f1 <- as.formula(paste("Y ~ Gene * Covariate_interaction +", paste(colnames(df)[4:ncol(df)], collapse = "+")))
      fit <- glm(f1, family = binomial, data = df)
      res.df[i, grep("_Variable", colnames(res.df))] <- c(round(exp(coef(fit))[2], rounding_factor), 
                                                          round(exp(confint.default(fit, level = 0.95))[2, 1], rounding_factor), 
                                                          round(exp(confint.default(fit, level = 0.95))[2, 2], rounding_factor), 
                                                          summary(fit)$coefficients[2, 3], 
                                                          summary(fit)$coefficients[2, 4])
      res.df[i, grep(covariate_interaction, colnames(res.df))] <- c(round(exp(coef(fit))[3], rounding_factor), 
                                                                    round(exp(confint.default(fit, level = 0.95))[3, 1], rounding_factor), 
                                                                    round(exp(confint.default(fit, level = 0.95))[3, 2], rounding_factor), 
                                                                    summary(fit)$coefficients[3, 3], 
                                                                    summary(fit)$coefficients[3, 4])
      res.df[i, grep("Interaction", colnames(res.df))] <- c(round(exp(coef(fit))[4+length(covariates_adjustment)], rounding_factor), 
                                                            round(exp(confint.default(fit, level = 0.95))[4+length(covariates_adjustment), 1], rounding_factor), 
                                                            round(exp(confint.default(fit, level = 0.95))[4+length(covariates_adjustment), 2], rounding_factor), 
                                                            summary(fit)$coefficients[4+length(covariates_adjustment), 3], 
                                                            summary(fit)$coefficients[4+length(covariates_adjustment), 4])
      k = 4
      for(j in 1:length(covariates_adjustment)){
        res.df[i, grep(covariates_adjustment[j], colnames(res.df))] <- c(round(exp(coef(fit))[k], rounding_factor), 
                                                                         round(exp(confint.default(fit, level = 0.95))[k, 1], rounding_factor), 
                                                                         round(exp(confint.default(fit, level = 0.95))[k, 2], rounding_factor), 
                                                                         summary(fit)$coefficients[k, 3], 
                                                                         summary(fit)$coefficients[k, 4])
        k <- k+1
      }
    }
  }
  fdr <- res.df[, grep("^P_", colnames(res.df))]
  colnames(fdr) <- gsub("^P_", "FDR_", colnames(fdr))
  fdr <- apply(fdr, 2, p.adjust, "BH")
  res.df <- as.data.frame(cbind(res.df, fdr))
  colOrder <- grep("Variable", colnames(res.df), value = T)
  colOrder <- c(colOrder, grep(covariate_interaction, colnames(res.df), value = T))
  colOrder <- c(colOrder, grep("Interaction", colnames(res.df), value = T))
  
  if(!is.null(covariates_adjustment)){
    for(z in 1:length(covariates_adjustment)){
      colOrder <- c(colOrder, grep(covariates_adjustment[z], colnames(res.df), value = T))
    }
  }
  res.df <- res.df[, colOrder]
  if(sort){
    res.df <- res.df[order(res.df$P_Interaction), ]
  }
  return(res.df)
}


#*******************************
# Sankey plot
#*******************************

sankey_plot<-function(df, node_width=0.2, node_space=10, node_color="black", node_alpha=0.4,
                   flow_smooth=8, node_label_size=4, node_label_color="black", n_size=4, x_label_size=13, x_label_color="black",
                   legend_name="", fill_colors=c("blue","red"), legend_text_size=12, legend_title_size=14,
                   plot_title="", plot_subtitle_size=13,plot_title_size=16,print_numbers=TRUE, factor.levels=NULL, simPVal=TRUE){
  
  library(dplyr)
  library(ggplot2)
  library(ggsankey)
  
  df <- df %>%
    make_long(colnames(df))
  if(is.null(factor.levels)){
    factor.levels<-union(levels(factor(df$node)), levels(factor(df$next_node)))
  }
  df$node<-factor(df$node,levels=factor.levels)
  df$next_node<-factor(df$next_node,levels=factor.levels)
  
  
  # prepare labels and labels positions
  
  stages<-levels(df$x)
  count.stages<-vector("list", length = length(stages))
  names(count.stages)<-stages
  for(i in 1:length(count.stages)){
    count.stages[[i]]<-table(df$node[df$x%in%stages[i]])
  }
  count.stages<-lapply(count.stages, function(x)x[x>0])
  
  lab.pos.stages<-lapply(count.stages, cumsum)
  for(i in 1:length(lab.pos.stages)){
    lab.pos.stages[[i]][1]<-lab.pos.stages[[i]][1]/2
    if(length(lab.pos.stages[[i]])>1){
      for(j in 2:length(lab.pos.stages[[i]])){
        lab.pos.stages[[i]][j]<-(lab.pos.stages[[i]][j] - (count.stages[[i]][j]/2)) + (node_space*(j-1))
      }
    }
  }
  
  # calculate p-value
  
  tnxn<-matrix(0,ncol=length(stages), nrow=length(factor.levels))
  rownames(tnxn)<-factor.levels
  colnames(tnxn)<-stages
  for(i in 1:length(stages)){
    tnxn[,i]<-as.vector(table(df$node[df$x==stages[i]]))
  }
  if(sum(tnxn<=5)>0){
    if(ncol(tnxn)<3){
      p.test<-fisher.test(tnxn, simulate.p.value=simPVal)$p.value
    } else {
      p.test<-fisher.test(tnxn, simulate.p.value=simPVal)$p.value
    }
  } else {
    p.test<-chisq.test(tnxn)$p.value
  }
  
  if(p.test<0.001){
    p.test<-format(p.test, digits=3, scientific=T)
  } else {
    p.test<-round(p.test,3)
  }
  
  # plot
  
  p<-ggplot(data=df, aes(x = x,
                         next_x = next_x,
                         node = node,
                         next_node = next_node,
                         fill=node)) +
    geom_sankey(smooth=flow_smooth,type="alluvial",
                node.color = node_color,
                alpha=node_alpha,
                space=node_space,
                width=node_width,
                show.legend = F) +
    theme_sankey() +
    theme(axis.text.x = element_text(size=x_label_size, color=x_label_color),
          legend.text = element_text(size=legend_text_size),
          legend.title = element_text(size=legend_title_size),
          plot.title = element_text(face="bold", size=plot_title_size, hjust=0.5),
          plot.subtitle = element_text(size=plot_subtitle_size, hjust=0.5)) +
    scale_fill_manual(name=legend_name,values=fill_colors) +
    xlab("") +
    ylab("") +
    labs(title=plot_title, subtitle = paste0("P = ",p.test))
  
  xs<-c(0.87,2.13)
  hjusts<-c(1,0)
  for(i in 1:length(lab.pos.stages)){
    for(j in 1:length(lab.pos.stages[[i]])){
      p <- p + geom_text(x = xs[i], y = lab.pos.stages[[i]][j], label = paste0(names(lab.pos.stages[[i]])[j]," (N = ",count.stages[[i]][j],")"), inherit.aes = F, size = node_label_size, hjust=hjusts[i])
    }
  }
  
  
  return(p)
  print(p)
}

