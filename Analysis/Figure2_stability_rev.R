library(SingleCellExperiment)

theme_paper <- function(textsize = 30){
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1), 
        text = element_text(size = textsize))
}

path = "~/Desktop/Projects/ASE_Spermatogenesis_Paper_rerun_Revisions/"
setwd(path)

# Read scDALI results 
data_results <- readRDS("~/Desktop/Projects/ASE_Spermatogenesis_Paper_rerun_Revisions/Data/processed/genes_test.rds")

### read processed data and generate allelic ratios
list_data <- readRDS("~/Desktop/Projects/ASE_Spermatogenesis_Paper_rerun_Revisions//Data/processed/binned_data.rds")

data_f0_b6 <- list_data[[1]]
data_f0_cast <- list_data[[2]]
data_f1_b6 <- list_data[[3]]
data_f1_cast <- list_data[[4]]

sfs = c(sum(data_f0_b6), sum(data_f0_cast)) / sum(data_f0_b6)
data_f0_b6_norm <- data_f0_b6 / sfs[[1]]
data_f0_cast_norm <- data_f0_cast / sfs[[2]]

data_f1_b6_scaled <- t(t(data_f1_b6) / colSums(data_f1_b6))
data_f1_cast_scaled <- t(t(data_f1_cast) / colSums(data_f1_cast))

genes_use <- rownames(data_f0_b6_norm)
rs_f1 <- (data_f1_b6[genes_use, ]) / (data_f1_b6[genes_use, ] + data_f1_cast[genes_use, ])
#rs_f0 <- (data_f0_b6_norm[genes_use, ] + 1) / (data_f0_b6_norm[genes_use, ] + data_f0_cast_norm[genes_use, ] + 2)
#rs_f1 <- (data_f1_b6[genes_use, ] + 1) / (data_f1_b6[genes_use, ] + data_f1_cast[genes_use, ] + 2)

exp_cutoff <- .5

# will be useful to look at results in different clusters of genes
clustering_df <- readRDS("./Data/processed/clustering_results.rds")
genes_with_cis <- rownames(clustering_df)

# 

## Partition genes into expressed, then up- and down-regulated
total_expression <- data_f1_b6 + data_f1_cast

normalization_factors <- colMeans(log(total_expression + 1))
normalized_centered_matrix <- t(scale(t(log(total_expression[genes_with_cis, ] + 1)) / normalization_factors))
#normalized_centered_matrix <- normalized_centered_matrix[genes_with_cis, ]

dim(normalized_centered_matrix) # we look at 1116 genes with dynamic cis-effects

## Analysis overview: 
# We want to address how ASE relates to up- and downregulation of genes

# First classify genes by their expression patterns
# Get genes with dynamic cis effects
#genes_with_cis <- data_results[p.adjust(data_results$dali_pval_polynomial) < 0.01, ]

# ComplexHeatmap::Heatmap(log(total_expression[rownames(genes_with_cis), ] + 1), cluster_columns = F)
# ComplexHeatmap::Heatmap(t(scale(t(log(total_expression[rownames(genes_with_cis), ] + 1)))), cluster_columns = F)

loess_smooth_matrix <- function(mat){
  smoothed_matrix <- t(apply(mat, 1, function(x){
    df_here = data.frame(Index = 1:ncol(mat), Values = x)
    loess_fit = loess(data = df_here, Values ~ Index)$fitted
  }))
  smoothed_matrix
}

# Run loess fits for total expression
normalized_centered_matrix_loess <- loess_smooth_matrix(normalized_centered_matrix)
left_annotation = ComplexHeatmap::rowAnnotation(df = clustering_df[order(clustering_df$Cluster_Expression), ])
ComplexHeatmap::Heatmap(normalized_centered_matrix_loess[rownames(clustering_df[order(clustering_df$Cluster_Expression), ]), ], cluster_columns = F, left_annotation = left_annotation, cluster_rows = F)

# find where maximum expression is along trajectory
peaks <- apply(normalized_centered_matrix_loess, 1, function(x){which(x == max(x))})

# compute average ASE during up- and down-regulation --> plot
average_ase_updown <- lapply(rownames(normalized_centered_matrix_loess), function(gene){
  peak_here <- peaks[[gene]]
  expressed_intervals <- total_expression[gene, ] > 0.1
  ases = rs_f1[gene, ]
  prev_ase = mean(ases[1:peak_here][expressed_intervals[1:peak_here]], na.rm = T)
  post_ase = mean(ases[peak_here:length(ases)][expressed_intervals[peak_here:length(ases)]], na.rm = T)
  return(c("Pre" = prev_ase, "Post" = post_ase))
})

average_ase_updown <- lapply(rownames(normalized_centered_matrix_loess), function(gene){
  peak_here <- peaks[[gene]]
  expressed_intervals <- total_expression[gene, ] > 0.1
  ases = rs_f1[gene, ]
  prev_ase = mean(ases[1:peak_here][expressed_intervals[1:peak_here]], na.rm = T)
  post_ase = mean(ases[peak_here:length(ases)][expressed_intervals[peak_here:length(ases)]], na.rm = T)
  return(c("Pre" = prev_ase, "Post" = post_ase))
})

df_test <- do.call("rbind", average_ase_updown) %>% 
  data.frame() %>% 
  add_column(Gene = rownames(normalized_centered_matrix_loess)) %>%
  add_column(Peak = peaks) %>%
  dplyr::filter(Peak > 1 & Peak < 100)

## get ase in quantiles of expression
genes_expressed <- rownames(normalized_centered_matrix_loess)[rowSums(total_expression[rownames(normalized_centered_matrix_loess), ] > 1) > 0]
genes_expressed <- df_test$Gene
average_ase_updown <- lapply(genes_expressed, function(gene){
  expressed_intervals <- total_expression[gene, ] > .1
  ases = rs_f1[gene, ][expressed_intervals]
  return(unlist(lapply(split(ases, cut(seq_along(ases), 5, labels = F)), mean)))
})
average_ase_updown <- average_ase_updown[unlist(lapply(average_ase_updown, length)) == 5]
average_ase_updown_dfs <- do.call("rbind", average_ase_updown)

# In which percentile is ASE the strongest?
test <- abs(do.call("rbind", average_ase_updown) - 0.5)
scaled_test <- data.frame(t(apply(test, 1, function(x){x - min(x)})))
colnames(test) <- paste0("Quantile_", 1:5)
colnames(scaled_test) <- paste0("Quantile_", 1:5)

color_vector <- setNames(c("green", "grey", "grey", "grey", "cyan"), nm = paste0("Quantile_", 1:5))
test %>%
  data.frame() %>% pivot_longer(-c()) %>%
  ggplot(aes(x = name, y = value, fill = name)) + geom_violin() + 
    geom_boxplot(width = 0.1, fill = "white", outlier.color = NA) +  theme_paper(textsize = 40) + 
    ylab("Absolute effect size (Allelic ratio - 0.5)") + xlab("") + 
    scale_fill_manual(values = color_vector) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
    theme(legend.position = "None") + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, 0.51))
ggsave("~/Desktop/Projects/ASE_Spermatogenesis_Paper_rerun_Revisions/Plots/FigureS3_rev//FigS3_effect_size_per_quantile.pdf", width = 14, height = 12)

scaled_test$max_quantile <- as.numeric(unlist(lapply(1:nrow(scaled_test), function(i){x = scaled_test[i, ]; which(x == max(x))[[1]]})))
scaled_test %>%
  data.frame() %>%
  mutate(max_quantile = paste0("Quantile_", max_quantile)) %>%
  ggplot(aes(x = max_quantile, fill = max_quantile)) + geom_bar(col = "black") + theme_paper(textsize = 40) + 
    xlab("") + ylab("Number of genes") +  scale_fill_manual(values = color_vector) + theme(legend.position = "None") + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, 390))
ggsave("~/Desktop/Projects/ASE_Spermatogenesis_Paper_rerun_Revisions/Plots/FigureS3_rev/FigS3_effect_size_per_quantile_barplot.pdf", width = 14, height = 12)

average_ase_updown_df <- do.call("rbind", average_ase_updown)
df_here <- data.frame(
  Pre = average_ase_updown_df[,1], 
  Post = average_ase_updown_df[,5], 
  Gene = genes_expressed, 
  Peak = peaks[genes_expressed]
)

eps = .05

df_test %>%
  dplyr::filter(Peak > 1 & Peak < 100) %>%
  mutate(Pre_cent = abs(Pre - 0.5)) %>%
  mutate(Post_cent = abs(Post - 0.5)) %>%
  mutate(Pre_cent = replace_na(Pre_cent, 0.5)) %>%
  mutate(Post_cent = replace_na(Post_cent, 0.5)) %>%
  mutate(dASE = (Post_cent * 1) - (Pre_cent * 1)) %>%
  dplyr::filter(!is.na(dASE)) %>%
  mutate(Category = ifelse(dASE > 0, "Late stronger", "Early stronger")) %>%
  mutate(Category = ifelse(sign(Pre_cent) != sign(Post_cent), "Switch", Category)) %>%
  mutate(Category = ifelse(abs(dASE) < 0.05, "No change", Category)) %>%
  add_column(expression_cluster = clustering_df[.$Gene, ]$Cluster_Expression) -> for_plot

for_plot %>%
  ggplot(aes(Pre, Post, col = Category, size = abs(dASE))) + geom_point()

for_plot %>%
  ggplot(aes(Pre, Post, col = Category)) + 
    geom_hline(yintercept = 0.5, linetype = 'dashed') + 
    geom_vline(xintercept = 0.5, linetype = 'dashed') + 
    geom_hline(yintercept = 0.5 + eps, linetype = 'dashed') + 
    geom_hline(yintercept = 0.5 - eps, linetype = 'dashed') + 
    geom_vline(xintercept = 0.5 + eps, linetype = 'dashed') + 
    geom_vline(xintercept = 0.5 - eps, linetype = 'dashed') + 
    geom_point(size = 3) + geom_abline(linetype = "dashed") + 
    xlab("AI during up-regulation") + ylab("AI during down-regulation") + theme_paper(textsize = 40) + 
    scale_color_manual(values = c("green", "cyan", "grey"))
ggsave("~/Desktop/Projects/ASE_Spermatogenesis_Paper_rerun_Revisions/Plots/FigureS3_rev/FigS3_early_late_categories_scatter.pdf")

df_here %>%
  dplyr::filter(Peak > 5 & Peak < 100) %>%
  mutate(Pre = abs(Pre - 0.5)) %>%
  mutate(Post = abs(Post - 0.5)) %>%
  mutate(dASE = Post - Pre) %>%
  dplyr::filter(!is.na(dASE)) %>%
  mutate(Category = ifelse(dASE > 0, "Late stronger", "Early stronger")) %>%
  mutate(Category = ifelse(sign(Pre) != sign(Post), "Switch", Category)) %>%
  mutate(Category = ifelse(abs(dASE) < eps, "No change", Category)) %>%
  add_column(expression_cluster = clustering_df[.$Gene, ]$Cluster_Expression) -> testy

round(table(testy$Category) / sum(table(testy$Category)), digits = 2)







