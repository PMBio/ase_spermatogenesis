##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
library(randomForest)
require(pROC)

rf_cv_auc <- function(data_fit, folds = 5){     
  n_chunks <- floor(nrow(data_fit) / folds)
  data_fit <- data_fit[sample(nrow(data_fit)), ]
  df_folds <- split(data_fit, (seq(nrow(data_fit))-1) %/% n_chunks)
  lapply(1:(length(df_folds) - 1), function(i){
    train_data <- do.call('rbind', df_folds[-i])
    test_data <- df_folds[[i]]
    rf <- randomForest(
      objective ~ .,
      data = train_data,
      importance = T, 
    )
    prediction <- stats::predict(rf, test_data, type = "vote")
    rf.roc <- roc(test_data$objective, prediction[,2])
    return(auc(rf.roc))
  })
  
}

genes_test <- readRDS("~/Desktop/genes_test.rds")
all_genes <- rownames(genes_test)

genes_test <- readRDS("~/Desktop/genes_test.rds")
genes_save <- rownames(genes_test)
genes_test <- data.frame(apply(genes_test, 2, as.numeric))
rownames(genes_test) <- genes_save
genes_test <- genes_test[!is.na(rownames(genes_test)), ]
genes_test <- genes_test[!is.na(genes_test$chromosome), ]
genes_test <- genes_test[!is.na(genes_test$BF_trans), ]
genes_test <- genes_test[!is.na(genes_test$BF_cis_dynamic), ]
genes_test <- genes_test[!is.na(genes_test$BF_trans_static), ]
genes_test <- genes_test[!is.na(genes_test$BF_trans_dynamic), ]

genes_test$BF_cis_static <- genes_test$NLL_null - genes_test$NLL_alt
genes_test[genes_test$NLL_null < 0 | genes_test$NLL_alt < 0, ]$BF_cis_static <- 0
genes_test[genes_test$BF_cis_static < 0, ]$BF_cis_static <- 0

BF_cutoff <- 10
BF_cutoff_lower <- 2
genes_cis_static <- rownames(genes_test[genes_test$BF_cis_static > BF_cutoff & genes_test$BF_cis_dynamic, ])
genes_cis_dynamic <- rownames(genes_test[genes_test$BF_cis_dynamic > BF_cutoff, ])
genes_cis <- union(genes_cis_static, genes_cis_dynamic)
genes_trans_static <- rownames(genes_test[genes_test$BF_trans_static > BF_cutoff & genes_test$BF_trans_dynamic, ])
genes_trans_dynamic <- rownames(genes_test[genes_test$BF_trans_dynamic > BF_cutoff, ])
genes_trans <- union(genes_trans_static, genes_trans_dynamic)
genetic_effects <- union(genes_cis, genes_trans)

##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### COMPARISON CIS GENERAL GENES
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 

data_fit <- readRDS("~/Desktop/rforest_analysis_files/all_features.rds")

data_fit[is.na(data_fit)] <- 0
data_fit <- data.frame(scale(data_fit))

# comparison: cis general
data_fit <- data_fit[all_genes, ]
data_fit$objective <- rownames(data_fit) %in% genes_cis
data_fit$objective <- factor( data_fit$objective )
data_fit <- data_fit[!is.na(data_fit$objective), ]
table(data_fit$objective)

set.seed(1234)
rf <- randomForest(
  objective ~ .,
  data=data_fit,
  importance = T, 
)

rf.roc <- roc(data_fit$objective, rf$votes[,2])
auc_training <- auc(rf.roc)
auc_test <- rf_cv_auc(data_fit, folds = 5)

importance_df <- data.frame(rf$importance)
importance_df$Feature <- rownames(importance_df)
importance_df <- importance_df[order(importance_df$MeanDecreaseAccuracy, decreasing = T), ]

auc_training_cis <- auc_training
auc_test_cis <- auc_test
importance_df_cis <- importance_df

##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### COMPARISON TRANS GENERAL GENES
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 

data_fit <- readRDS("~/Desktop/rforest_analysis_files/all_features.rds")

data_fit[is.na(data_fit)] <- 0
data_fit <- data.frame(scale(data_fit))

# comparison: all trans
data_fit <- data_fit[all_genes, ]
data_fit$objective <- rownames(data_fit) %in% genes_trans
data_fit$objective <- factor( data_fit$objective )
data_fit <- data_fit[!is.na(data_fit$objective), ]
table(data_fit$objective)

set.seed(1234)
rf <- randomForest(
  objective ~ .,
  data=data_fit,
  importance = T, 
)

rf.roc <- roc(data_fit$objective, rf$votes[,2])
auc_training <- auc(rf.roc)
auc_test <- rf_cv_auc(data_fit, folds = 5)

importance_df <- data.frame(rf$importance)
importance_df$Feature <- rownames(importance_df)
importance_df <- importance_df[order(importance_df$MeanDecreaseAccuracy, decreasing = T), ]

auc_training_trans <- auc_training
auc_test_trans <- auc_test
importance_df_trans <- importance_df

##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### COMPARISON CIS DYNAMIC GENES
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 

data_fit <- readRDS("~/Desktop/rforest_analysis_files/all_features.rds")

data_fit[is.na(data_fit)] <- 0
data_fit <- data.frame(scale(data_fit))

# comparison: dynamic cis effects
data_fit <- data_fit[genes_cis, ]
data_fit$objective <- rownames(data_fit) %in% genes_cis_dynamic
data_fit$objective <- factor( data_fit$objective )
data_fit <- data_fit[!is.na(data_fit$objective), ]
table(data_fit$objective)

set.seed(1234)
rf <- randomForest(
  objective ~ .,
  data=data_fit,
  importance = T, 
)

rf.roc <- roc(data_fit$objective, rf$votes[,2])
auc_training <- auc(rf.roc)
auc_test <- rf_cv_auc(data_fit, folds = 5)

importance_df <- data.frame(rf$importance)
importance_df$Feature <- rownames(importance_df)
importance_df <- importance_df[order(importance_df$MeanDecreaseAccuracy, decreasing = T), ]

auc_training_cis_dynamic <- auc_training
auc_test_cis_dynamic <- auc_test
importance_df_cis_dynamic <- importance_df

##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### COMPARISON TRANS DYNAMIC GENES
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 

# comparison: dynamic trans effects
data_fit <- readRDS("~/Desktop/rforest_analysis_files/all_features.rds")

data_fit[is.na(data_fit)] <- 0
data_fit <- data.frame(scale(data_fit))

data_fit <- data_fit[genes_trans, ]
data_fit$objective <- rownames(data_fit) %in% genes_trans_dynamic
data_fit$objective <- factor( data_fit$objective )
data_fit <- data_fit[!is.na(data_fit$objective), ]
table(data_fit$objective)

set.seed(1234)
rf <- randomForest(
  objective ~ .,
  data=data_fit,
  importance = T, 
)

rf.roc <- roc(data_fit$objective, rf$votes[,2])
auc_training <- auc(rf.roc)
auc_test <- rf_cv_auc(data_fit, folds = 5)

importance_df <- data.frame(rf$importance)
importance_df$Feature <- rownames(importance_df)
importance_df <- importance_df[order(importance_df$MeanDecreaseAccuracy, decreasing = T), ]

auc_training_trans_dynamic <- auc_training
auc_test_trans_dynamic <- auc_test
importance_df_trans_dynamic <- importance_df

##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
##### COMPARISON CIS vs TRANS
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 

# comparison: dynamic trans effects
data_fit <- readRDS("~/Desktop/rforest_analysis_files/all_features.rds")

data_fit[is.na(data_fit)] <- 0
data_fit <- data.frame(scale(data_fit))

data_fit <- data_fit[union(genes_cis, genes_trans), ]
data_fit$objective <- rownames(data_fit) %in% genes_trans
data_fit$objective <- factor( data_fit$objective )
data_fit <- data_fit[!is.na(data_fit$objective), ]
table(data_fit$objective)

set.seed(1234)
rf <- randomForest(
  objective ~ .,
  data=data_fit,
  importance = T, 
)

rf.roc <- roc(data_fit$objective, rf$votes[,2])
auc_training <- auc(rf.roc)
auc_test <- rf_cv_auc(data_fit, folds = 5)

importance_df <- data.frame(rf$importance)
importance_df$Feature <- rownames(importance_df)
importance_df <- importance_df[order(importance_df$MeanDecreaseAccuracy, decreasing = T), ]

auc_training_trans_vs_cis <- auc_training
auc_test_trans_trans_vs_cis <- auc_test
importance_df_trans_vs_cis <- importance_df

####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### 
## Plotting 
####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### 

# First plot all AUCs of all comparisons together

all_ucs_test <- data.frame(
  "AllCis" = unlist(auc_test_cis), 
  "AllTrans" = unlist(auc_test_trans), 
  "DynamicCis" = unlist(auc_test_cis_dynamic), 
  "DynamicTrans" = unlist(auc_test_trans_dynamic), 
  "TransvsCis" = unlist(auc_test_trans_trans_vs_cis)
) %>% pivot_longer(-c())

all_ucs_test_means <- all_ucs_test %>% dplyr::group_by(name) %>%
  summarise(value = mean(value)) %>%
  add_column(col = "AUC_CV")

all_training_aucs <- data.frame(
  "AllCis" = as.numeric(unlist(auc_training_cis)), 
  "AllTrans" = as.numeric(unlist(auc_training_trans)), 
  "DynamicCis" = as.numeric(unlist(auc_training_cis_dynamic)), 
  "DynamicTrans" = as.numeric(unlist(auc_training_trans_dynamic)), 
  "TransvsCis" = as.numeric(unlist(auc_training_trans_vs_cis))
) %>% pivot_longer(-c()) %>%
  add_column(col = "AUC_train")

rbind(all_ucs_test_means, all_training_aucs) %>%
  ggplot(aes(x = name, y = value, fill = col)) + 
  geom_bar(stat = "identity", position = "dodge")  + 
  scale_y_continuous(limits = c(0, 1), expan = c(0, 0)) + 
  scale_fill_manual(values = c("grey", "black")) + 
  theme_classic() + 
  xlab("") + ylab("") + 
  theme(text = element_text(size = 30)) + 
  ggtitle("Model performance") + 
  geom_hline(col = "red", linetype = "dashed", yintercept = 0.5)
ggsave("~/Desktop/PaperWriting/Fig6/AllModelAUCs.pdf")

# Now we plot feature importances for the individual models
# first cis effects

plot_df_1 <- importance_df_cis %>% add_column(Effect = "General") %>% 
  mutate(Feature = gsub("^VAR_", "", Feature)) %>%
  mutate(Feature = gsub("^evolution_|^structural_|^REGULATION_", "", Feature)) %>%
  arrange(MeanDecreaseAccuracy) %>% mutate(Feature = factor(Feature, levels = reorder(Feature, MeanDecreaseAccuracy)))
plot_df_2 <- importance_df_cis_dynamic %>% 
  mutate(Feature = gsub("^VAR_", "", Feature)) %>%
  mutate(Feature = gsub("^evolution_|^structural_|^REGULATION_", "", Feature)) %>%
  add_column(Effect = "Dynamic")

rbind(plot_df_1, plot_df_2) %>%
  ggplot(aes(Feature, MeanDecreaseAccuracy, col = Effect)) + 
    geom_point(size = 3) + coord_flip() + ylim(-0.001, 0.015) + 
    geom_hline(yintercept = 0, linetype = "dashed") + 
    ggtitle("Predictors of Cis effects") + 
    theme_classic() + 
    xlab("") + 
    theme(text = element_text(size=30)) + 
    scale_color_manual(values = c("darkred", "darkgrey"))
ggsave("~/Desktop/PaperWriting/Fig6/Fig6_FeatureImportance_Cis.pdf")

# reduced plots
plot_df_1 <- importance_df_cis %>% add_column(Effect = "General") %>% 
  mutate(Feature = gsub("^VAR_", "", Feature)) %>%
  mutate(Feature = gsub("^evolution_|^structural_|^REGULATION_", "", Feature)) %>%
  arrange(MeanDecreaseAccuracy) %>% mutate(Feature = factor(Feature, levels = reorder(Feature, MeanDecreaseAccuracy)))
plot_df_2 <- importance_df_cis_dynamic %>% 
  mutate(Feature = gsub("^VAR_", "", Feature)) %>%
  mutate(Feature = gsub("^evolution_|^structural_|^REGULATION_", "", Feature)) %>%
  add_column(Effect = "Dynamic")

features_show <- setNames(importance_df_cis$MeanDecreaseAccuracy + importance_df_cis_dynamic[rownames(importance_df_cis), ]$MeanDecreaseAccuracy, 
                          rownames(importance_df_cis))
features_show <- names(sort(features_show, decreasing = T)[1:12])

convert_names <- setNames(c("Exp. Level", "PhastCons (exon)", "PhyloP (exon)", "%SNP (exon)", "3' UTR length", "%SNP 3' UTR", "#Exons", 
                            "Exp. Variability (inter-celltype)", "Exp. Variability (intra-celltype)", "5' UTR length", "%SNP (5' UTR)", "Gene Length"), 
                          features_show)

plot_df_1 <- plot_df_1[features_show, ]
plot_df_2 <- plot_df_2[features_show, ]

plot_df_1$order <- 13 - 1:nrow(plot_df_1)
plot_df_2$order <- 13 - 1:nrow(plot_df_2)

plot_df_1$FeaturePretty <- convert_names[rownames(plot_df_1)]
plot_df_2$FeaturePretty <- convert_names[rownames(plot_df_2)]

features_highlight <- c("VAR_evolution_phastcons_exon", "VAR_expression_level", "VAR_evolution_snpdensity_F1_exon")

features_highlight_index <- plot_df_1[features_highlight, ]$order

rbind(plot_df_1, plot_df_2) %>%
  ggplot(aes(reorder(FeaturePretty, order),  MeanDecreaseAccuracy, col = Effect)) + 
  geom_point(size = 7) + coord_flip() + ylim(-0.001, 0.015) + 
  ggtitle("Predictors of Cis effects") + 
  theme_classic() + 
  xlab("") + 
  theme(text = element_text(size=30)) + 
  scale_color_manual(values = c("red", "darkgrey")) + 
  ylab("Feature Importance") + theme(legend.position = "None") + 
  scale_y_continuous(expand = c(0, 0)) + 
  annotate(x  = 2, y = 0.010, geom = "text", label = paste0("Cis (gen) AUC: ", round(all_ucs_test_means[1, ]$value, digits = 3)), size = 10, col = "grey40") + 
  annotate(x  = 1, y = 0.010, geom = "text", label = paste0("Cis (gen) AUC: ", round(all_ucs_test_means[1, ]$value, digits = 3)), size = 10, col = "red") + 
  annotate(geom = "rect", xmin = features_highlight_index[[1]] + 0.4, xmax = features_highlight_index[[1]] - 0.4, ymin = 0, ymax = 0.015, 
           col = "black", fill = "yellow", alpha = 0.1) + 
  annotate(geom = "rect", xmin = features_highlight_index[[2]] + 0.4, xmax = features_highlight_index[[2]] - 0.4, ymin = 0, ymax = 0.015, 
           col = "black", fill = "yellow", alpha = 0.1) + 
  annotate(geom = "rect", xmin = features_highlight_index[[3]] + 0.4, xmax = features_highlight_index[[3]] - 0.4, ymin = 0, ymax = 0.015, 
           col = "black", fill = "yellow", alpha = 0.1)
ggsave("~/Desktop/PaperWriting/Fig6/Fig6_FeatureImportance_Cis_Reduced.pdf")

# Now we plot feature importances for the individual models
# then trans effects

plot_df_1 <- importance_df_trans %>% add_column(Effect = "General") %>% 
  mutate(Feature = gsub("^VAR_", "", Feature)) %>%
  mutate(Feature = gsub("^evolution_|^structural_|^REGULATION_", "", Feature)) %>%
  arrange(MeanDecreaseAccuracy) %>% mutate(Feature = factor(Feature, levels = reorder(Feature, MeanDecreaseAccuracy)))
plot_df_2 <- importance_df_trans_dynamic %>% 
  mutate(Feature = gsub("^VAR_", "", Feature)) %>%
  mutate(Feature = gsub("^evolution_|^structural_|^REGULATION_", "", Feature)) %>%
  add_column(Effect = "Dynamic")

rbind(plot_df_1, plot_df_2) %>%
  ggplot(aes(Feature, MeanDecreaseAccuracy, col = Effect)) + 
  geom_point(size = 3) + coord_flip() + ylim(-0.001, 0.015) + 
  geom_hline(yintercept = 0, linetype = "dashed") + 
  ggtitle("Predictors of Trans effects") + 
  theme_classic() + 
  xlab("") + 
  theme(text = element_text(size=30)) + 
  scale_color_manual(values = c("darkred", "darkgrey"))
ggsave("~/Desktop/PaperWriting/Fig6/Fig6_FeatureImportance_Trans.pdf")

## Now we plot features with predictive differences

data_fit <- readRDS("~/Desktop/rforest_analysis_files/all_features.rds")

data_fit[is.na(data_fit)] <- 0
data_fit <- data.frame(scale(data_fit))

features_show <- c("VAR_evolution_phastcons_exon", "VAR_expression_level", "VAR_evolution_snpdensity_F1_exon")

all_genes <- setNames(rep("None", nrow(data_fit)), rownames(data_fit))
all_genes[genes_cis] <- "cis"
all_genes[genes_cis_dynamic] <- "dynamic cis"
all_genes <- factor(all_genes, levels = c("None", "cis", "dynamic cis"))


df_plot <- 
  data_fit[,features_show] %>%
  data.frame() %>%
  add_column(effect = all_genes) %>%
  pivot_longer(-effect) %>%
  mutate(name = convert_names[as.character(name)]) %>%
  mutate(name = factor(name, levels = c("Exp. Level", "PhastCons (exon)", "%SNP (exon)"))) %>%
  mutate(effect = factor(effect, levels = c("dynamic cis", "cis", "None")))

df_plot %>%
 ggplot(aes(x = effect, y = value)) + 
  geom_violin(fill = "grey") +
  stat_summary(size = 1) + 
  # geom_boxplot(width = 0.1) +
  facet_wrap(~name, scales = "free", nrow = 3) +
  coord_flip() + 
   ggpubr::stat_compare_means(comparisons = list(c("None", "cis"), c("None", "dynamic cis"), c("cis", "dynamic cis")), 
                              label = "p.signif", coord.flip = TRUE, size = 10) + 
  theme_classic() + 
  theme(text = element_text(size = 50)) + 
  xlab("") + ylab("")
ggsave("~/Desktop/PaperWriting/Fig6/Fig6_CisEffectsViolin.pdf")

## Now we plot features with predictive difference

