# This script generates the supplement for figure 2 that isnt generated in the script for the main figure 2

# Functions
library(ggplot2)
library(viridis)

setwd("~/Desktop/Projects/ASE_Spermatogenesis_Paper_rerun_Revisions//")

library_colors <- list(
  "B6_Sample1" = "black",
  "B6_Sample2" = "grey",
  "CAST_Sample3" = "chocolate4",
  "CAST_Sample4" = "chocolate",
  "F1_B6_CAST_Sample5" = "purple4",
  "F1_B6_CAST_Sample6" = "purple"
)

sample_colors <- list(
  "B6" = "black", 
  "CAST" = "chocolate", 
  "F1_B6_CAST" = "purple"
)

theme_paper <- function(textsize = 20){
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1), 
        text = element_text(size = textsize))
}

calc_SS <- function(df) sum(as.matrix(dist(df)^2)) / (2 * nrow(df))
calc_SS_partition <- function(df, clusters){
  data_here = lapply(unique(clusters), function(cl){
    dd = df[clusters == cl, ]
    calc_SS(dd)
  })
  sum(unlist(data_here))
}

make_pseudotime_smooth <- function(data, pt, n_partitions = 100){
  pseudotime_here <- pt
  pseudotime_here <- pseudotime_here - min(pseudotime_here)
  pseudotime_here <- pseudotime_here / max(pseudotime_here) # between 0 and 1
  intervals <- c(0, 1:n_partitions / n_partitions)
  
  expression_smoothed <- do.call("cbind", 
                                 lapply(1:n_partitions, function(i){
                                   p_lower <- intervals[i]
                                   p_upper <-  intervals[i + 1]
                                   ix <- which(p_lower <= pseudotime_here & pseudotime_here <= p_upper)
                                   rowMeans(data[,ix])
                                 }))
  expression_smoothed
}

# read in data from persistent allelic imbalance tests
summary_data_meanTest <- readRDS("~/Desktop/PhD/Projects/ASE_Spermatogenesis_Paper_rerun_Revisions//Data/processed//Dali_LRT_calls.csv")


# exclude genes with mapping effects (see get mapping effects) and genes with no reads
summary_data_meanTest <- summary_data_meanTest[!summary_data_meanTest$Mapping_effect & 
                                                 summary_data_meanTest$coverage_ref > 1 & summary_data_meanTest$coverage_alt > 1, ]

# read results from dynamic tests and add to summary df
results_linear_kernel <- readRDS("~/Desktop/PhD/Projects/ASE_Spermatogenesis_Paper_rerun_Revisions/Data/processed//Dali_scoreTest_linear.csv")
results_polynomial_kernel <- readRDS("~/Desktop/PhD/Projects/ASE_Spermatogenesis_Paper_rerun_Revisions/Data/processed//Dali_scoreTest_poly.csv")
results_cluster_kernel <- readRDS("~/Desktop/PhD/Projects/ASE_Spermatogenesis_Paper_rerun_Revisions/Data/processed//Dali_scoreTest_discrete.csv")

summary_data_meanTest$pval_dali_linear <- results_linear_kernel
summary_data_meanTest$pval_dali_polynomial <- results_polynomial_kernel
summary_data_meanTest$pval_dali_discrete <- results_cluster_kernel

# exclude genes with failed convergence for persistent tests
table(-log10(summary_data_meanTest$padj) < 5000000)
summary_data_meanTest_here <- summary_data_meanTest[-log10(summary_data_meanTest$padj) < 5000000, ]

# get genes with corrected pvalue < 0.01 for persistent and dynamic tests
pval_cutoff <- 0.01
corrected_pvals_x <- data.frame(pval = unlist(summary_data_meanTest$pval), padj = summary_data_meanTest$padj)
corrected_pvals_x <- corrected_pvals_x[order(corrected_pvals_x$padj), ]
corrected_pvals_x_cutoff <- corrected_pvals_x[corrected_pvals_x$padj > 0.01, 1][[1]]
corrected_pvals_y <- data.frame(pval = (results_linear_kernel), padj = p.adjust(results_linear_kernel))
corrected_pvals_y <- corrected_pvals_y[order(corrected_pvals_y$padj), ]
corrected_pvals_y_cutoff <- corrected_pvals_y[corrected_pvals_y$padj > 0.01, 1][[1]]

#### Read and smooth data
# add heatmap of total expression
data_total_expression <- readRDS("~/Desktop/PhD/Projects/ASE_Spermatogenesis_Paper_rerun_Revisions/Data/processed/final_sce_f1_dataset_for_supplement_pseudotime.rds")
data_total_expression <- data_total_expression[,data_total_expression$Library != "Sample7"]
data_total_expression_save <- data_total_expression
data_total_expression <- data_total_expression[,data_total_expression$AnnotatedClusters %in% c("ES", "RS", "SC", "SG")]
data_total_expression <- data_total_expression[,data_total_expression$Library %in% c("Sample5", "Sample6")]

# show pseudotime for S2.1
pt <- data_total_expression_save$Pseudotime
# pt <- - (pt - max(pt))

# recompute umap dimensionality reduction only on germ cells
set.seed(12345)
corrected <- metadata(data_total_expression_save)$corrected[colnames(data_total_expression_save), ]
umap <- umap::umap(corrected)
reducedDims(data_total_expression_save)$UMAP <- umap$layout

data.frame(
  Umap1 = umap$layout[,1], 
  Umap2 = umap$layout[,2], 
  Species = unlist(data_total_expression_save$Sample), 
  Pseudotime = pt
) %>%
  ggplot(aes(Umap1, Umap2, col = Pseudotime)) + geom_point() + 
  scale_color_viridis() + theme_paper(textsize = 35) + facet_wrap(~Species) 
ggsave("./Plots/FigureS2/FigS2_1_pseudotime_umap.pdf", width = 20, height = 8)

# show sample distribution across pseudotime
# data.frame(
#   Library = paste0(data_total_expression_save$Sample, "_", data_total_expression_save$Library), 
#   Pseudotime = data_total_expression_save$Pseudotime
# ) %>%
#   ggplot(aes(x = Pseudotime, fill = Library)) + geom_density(alpha = 0.5) + 
#   scale_color_viridis() + theme_paper(textsize = 50) + scale_fill_manual(values = library_colors)
# ggsave("./Plots/FigureS2/FigS2_2_pseudotime_distribution.pdf")

data.frame(
  Library = paste0(data_total_expression_save$Sample, "_", data_total_expression_save$Library), 
  Species = unlist(data_total_expression_save$Sample), 
  Pseudotime = data_total_expression_save$Pseudotime
) %>%
  ggplot(aes(x = Pseudotime, fill = Library)) + geom_density(alpha = 0.5) + 
  scale_color_viridis() + theme_paper(textsize = 35) + scale_fill_manual(values = library_colors) + 
  facet_wrap(~Species) + ylab("Density of cells") + 
  theme(legend.position = "bottom")
ggsave("./Plots/FigureS2/FigS2_3_pseudotime_distribution_sep.pdf")

# gp_results_total <- readRDS("./Data/processed/GP_fits_all.rds")
# names(gp_results_total) <- rownames(summary_data_meanTest)
# gp_results_total <- gp_results_total[!unlist(lapply(gp_results_total, function(x){any(is.na(x))}))]
# 
# n_partitions = 100
# pseudotime_here <- data_total_expression$Pseudotime
# pseudotime_here <- pseudotime_here - min(pseudotime_here)
# max_pseudotime <- max(pseudotime_here)
# interval_width = max(max_pseudotime) / 50
# pseudotime_partition <- 1:n_partitions * max(max_pseudotime) / n_partitions
# 
# # asdasdasdasd
# expression_smoothed <- do.call("cbind", 
#                                lapply(pseudotime_partition, function(p){
#                                  p_lower <- max(0, p - interval_width)
#                                  p_upper <- min(max_pseudotime, p + interval_width)
#                                  ix <- which(p_lower <= pseudotime_here & pseudotime_here <= p_upper)
#                                  a = rowMeans(logcounts(data_total_expression[,names(ix)]))
#                                  return(a)
#                                }))
# 
# gp_latent_matrix <- t(do.call("cbind", lapply(gp_results_total, function(x){x$posterior_mean})))
# gp_latent_matrix <- gp_latent_matrix[,order(data_total_expression$Pseudotime)]
# rownames(gp_latent_matrix) <- names(gp_results_total)
# 
# dynamic_cis_effects <- rownames(summary_data_meanTest[p.adjust(summary_data_meanTest$pval_dali_polynomial) < 0.01, ])
# dynamic_cis_effects <- intersect(dynamic_cis_effects, rownames(gp_latent_matrix))
# cell_use <- 1:ncol(gp_latent_matrix) %% 10 == 0
# gp_latent_matrix <- gp_latent_matrix[dynamic_cis_effects, ]
# 
# gp_latent_matrix_scaled <- t(scale(t(gp_latent_matrix)))
# 
# colnames(gp_latent_matrix) <- names(pseudotime_here[order(pseudotime_here)])
# 
# ase_smoothed <- do.call("cbind", 
#                         lapply(pseudotime_partition, function(p){
#                           p_lower <- max(0, p - interval_width)
#                           p_upper <- min(max_pseudotime, p + interval_width)
#                           ix <- which(p_lower <= pseudotime_here & pseudotime_here <= p_upper)
#                           a = rowMeans(gp_latent_matrix[,names(ix)])
#                           return(a)
#                         }))
# 
# gp_latent_matrix <- gp_latent_matrix[,cell_use]
# 
# ## compute clusters on scaled gp fits for ASE
# clusterings <- lapply(1:30, function(i){cutree(hclust(dist(gp_latent_matrix_scaled)), k = i)})
# 
# SSDs <- lapply(clusterings, function(x){calc_SS_partition(gp_latent_matrix, x)}) # 7 looks ok
# plot(unlist(SSDs))
# ggplot(data.frame(Clusters = 1:30, SSD = unlist(SSDs)), aes(Clusters, SSD)) + geom_point() + geom_line() + theme_classic() + 
#   theme(text = element_text(size = 20)) + geom_hline(yintercept = SSDs[[7]], col = "red", linetype = "dashed") + 
#   ggtitle("ASE cluster choice")
# 
# clustering_choose <-  clusterings[[7]]
# 
# annotation_row_expression <- data.frame(
#   Cluster = clustering_choose,
#   Dummy = 1
# )
# 
# rownames(annotation_row_expression) <- rownames(gp_latent_matrix)
# annotation_row_expression <- annotation_row_expression[order(annotation_row_expression$Cluster),]
# 
# row_gaps <- cumsum(table(annotation_row_expression))
# points_plot <-  1:ncol(gp_latent_matrix) %% 10 == 0
# 
# ## compute clusters on scaled expression patterns
# expression_smoothed_cut <- expression_smoothed[rownames(gp_latent_matrix), ]
# 
# expression_smoothed_scaled <- t(scale(t(expression_smoothed_cut)))
# 
# clusterings_expression <- lapply(1:30, function(i){cutree(hclust(dist(expression_smoothed_scaled)), k = i)})
# SSDs_expression <- lapply(clusterings_expression, function(x){calc_SS_partition(expression_smoothed_scaled, x)}) # 7 looks ok
# 
# ggplot(data.frame(Clusters = 1:30, SSD = unlist(SSDs_expression)), aes(Clusters, SSD)) + geom_point() + geom_line() + theme_classic() + 
#   theme(text = element_text(size = 20)) + geom_hline(yintercept = SSDs_expression[[7]], col = "red", linetype = "dashed") + 
#   ggtitle("Expression cluster choice")
# 
# clusterings_reordered <- rank(168 - table(c(clusterings[[7]], 7)), ties.method = c("first"))[clustering_choose]
# 
# clusters_raw_ase <- lapply(1:30, function(i){cutree(hclust(dist(gp_latent_matrix)), k = i)})
# SSDs_raw_ase <- lapply(clusters_raw_ase, function(x){calc_SS_partition(gp_latent_matrix, x)}) # 7 looks ok
# plot(unlist(SSDs_raw_ase))
# 
# clustering_df <- data.frame(
#   Cluster_Expression = clusterings_expression[[7]],
#   Cluster_ASE = clusterings_reordered
# )
# 
# clustering_df <- clustering_df[order(clustering_df$Cluster_ASE, clustering_df$Cluster_Expression), ]

# get cluster meta-profiles
# ase_smoothed_scaled <- t(scale(t(ase_smoothed)))
# ase_smoothed_scaled <- ase_smoothed
# 
# meta_profiles_ase <- do.call("rbind", lapply(1:7, function(i){
#   genes = rownames(clustering_df[clustering_df$Cluster_ASE == i, ])
#   colMeans(ase_smoothed_scaled[genes, ])
# }))
# 
# ggplot(reshape2::melt(meta_profiles_ase), aes(x = Var2, y = value)) + geom_line() + facet_wrap(~Var1) + 
#   geom_hline(yintercept = 0.5, linetype = "dashed") + ggtitle("Cluster profiles ASE")
# 
# meta_profiles_expression <- do.call("rbind", lapply(1:7, function(i){
#   genes = rownames(clustering_df[clustering_df$Cluster_Expression == i, ])
#   colMeans(expression_smoothed_cut[genes, ])
# }))
# 
# ggplot(reshape2::melt(meta_profiles_expression), aes(x = Var2, y = value)) + geom_line() + facet_wrap(~Var1) + ggtitle("Cluster profiles Expression")

# confusion matrix
# 
# confusion_matrix <- table(clustering_df)
# row_sums <- rowSums(confusion_matrix)
# col_sums <- colSums(confusion_matrix)
# confusion_matrix_expected <- row_sums %o% col_sums / sum(confusion_matrix)
# 
# ggplot(data.frame(confusion_matrix), aes(Cluster_Expression, Cluster_ASE, fill = Freq)) + geom_tile() + 
#   scale_fill_gradientn(colours = c("white", "red")) + theme_classic() + theme(text = element_text(size = 20)) + 
#   geom_text(aes(label = Freq))
# 
# ggplot(data.frame(confusion_matrix / confusion_matrix_expected), aes(Cluster_Expression, Cluster_ASE, fill = Freq)) + geom_tile() + 
#   scale_fill_gradientn(colours = c("white", "red")) + theme_classic() + theme(text = element_text(size = 20)) + 
#   geom_text(aes(label = round(Freq, 2)))
# 
# clustering_df <- clustering_df[order(clustering_df$Cluster_Expression),]
# pheatmap::pheatmap(expression_smoothed_cut[rownames(clustering_df), ],
#                    cluster_cols = F, cluster_rows = F, show_rownames = F, show_colnames = F, 
#                    annotation_row = clustering_df, scale = "row",
#                    gaps_row = cumsum(table(clustering_df$Cluster_Expression)), 
#                    color = rev(RColorBrewer::brewer.pal(name = "RdBu", n = 10)))

# 
# gp_latent_matrix_pt <- make_pseudotime_smooth(gp_latent_matrix, sort(pseudotime_here[cell_use]), 100)

# colfunc <- colorRampPalette(c("black", "white", "chocolate"))
# clustering_df <- clustering_df[order(clustering_df$Cluster_ASE), ]
# pheatmap::pheatmap(gp_latent_matrix_pt[rownames(clustering_df), ],
#                    cluster_cols = F, cluster_rows = F, show_rownames = F, show_colnames = F, 
#                    annotation_row = clustering_df, scale = "row",
#                    gaps_row = cumsum(table(clustering_df$Cluster_Expression)),
#                    color = colfunc(10))

# clustering_df <- clustering_df[order(clustering_df$Cluster_raw_ASE), ]
# pheatmap::pheatmap(gp_latent_matrix_pt[rownames(clustering_df), ],
#                    cluster_cols = F, cluster_rows = F, show_rownames = F, show_colnames = F, 
#                    annotation_row = clustering_df, 
#                    gaps_row = cumsum(table(clustering_df$Cluster_Expression)),
#                    color = colfunc(10))

# pheatmap::pheatmap(gp_latent_matrix_pt[rownames(clustering_df), ],
#                    cluster_cols = F, cluster_rows = F, show_rownames = F, show_colnames = F, 
#                    annotation_row = clustering_df, scale = "row",
#                    gaps_row = cumsum(table(clustering_df$Cluster_Expression)),
#                    color = colfunc(10))

#cis_effects_binned <- make_pseudotime_smooth(gp_latent_matrix, sort(data_use[,data_use$Species == "CASTxB6"]$Pseudotime), 
#                                             n_partitions = 100)

#pheatmap::pheatmap(abs(cis_effects_binned - 0.5),
#                   cluster_cols = F, cluster_rows = T, show_rownames = F, show_colnames = F, 
#                   scale = "row",
#                   color = colfunc(10))

####

# source("~/Desktop/PhD/RCode/Multipurpose/R_gsea_functions.R")
# 
# gs_msigdb_h <- get_gs_database_msigdb()
# gs_msigdb_c5 <- get_gs_database_msigdb(which_set = "C5")
# #gs_kegg <- get_gs_database_kegg()
# 
# gsea_hallmark_res <- lapply(unique(clustering_df$Cluster_ASE), 
# function(i){
#   print(i)
#   genes = rownames(clustering_df[clustering_df$Cluster_ASE == i, ])
#   go_overrepresentation(genes, gs_list = gs_msigdb_c5, univ = rownames(data_f1), is_human = F)
# })
# 
# head(gsea_hallmark_res[[7]], n = 10)
# 
# gsea_total <- go_overrepresentation(rownames(clustering_df), gs_list = gs_msigdb_c5, 
#                                     univ = rownames(data_f1), is_human = F)
# 
# head(gsea_total, n = 10)
# 
# simpleCap <- function(x) {
#   s <- strsplit(x, " ")[[1]]
#   paste(toupper(substring(s, 1, 1)), substring(s, 2),
#         sep = "", collapse = " ")
# }
# 
# head(gsea_hallmark_res[[6]])
# 
# gs_here <- unlist(lapply(tolower(gs_msigdb_c5[["GO_RIBONUCLEOPROTEIN_COMPLEX"]]), simpleCap))
# 
# intersect(gs_here, rownames(clustering_df[clustering_df$Cluster_ASE == 5, ]))
# intersect(gs_here, rownames(clustering_df))

# gene = "Nusap1"
# plot_gene_GP(data_f1, gene, 
#              latent_ase = gp_results_total[[gene]], remove_zero = F) + 
#   ggtitle(gene)

# ggplot(data.frame(
#   Pseudotime = data_total_expression$Pseudotime, 
#   Expression = as.numeric(log(counts(data_total_expression["Tnp2", ])) + 1)
# ), aes(Pseudotime, Expression)) + geom_point()
# 
# colnames(gp_latent_matrix) <- paste0("Cell_", 1:ncol(gp_latent_matrix))
# column_annotation = data.frame(
#   Pseudotime = data_total_expression$Pseudotime, 
#   Total_Expression = log10(colSums(counts(data_total_expression)) + 1), 
#   Protamine_Expression = as.numeric(log(counts(data_total_expression["Tnp2", ]) + 1))
# )
# column_annotation <- column_annotation[order(column_annotation$Pseudotime), ]
# rownames(column_annotation) <- colnames(gp_latent_matrix)
# 
# column_annotation_collapsed <- data.frame(
#   t(make_pseudotime_smooth(t(column_annotation), column_annotation$Pseudotime, n_partitions = 100))
# )
# colnames(gp_latent_matrix_pt) <- paste0("Cell_", 1:100)
# rownames(column_annotation_collapsed) <- paste0("Cell_", 1:100)
# 
# pdf("~/Desktop/PaperWriting/Fig3/ASE_matrix_wrong.pdf")
# pheatmap::pheatmap(gp_latent_matrix[rownames(clustering_df[order(clustering_df$Cluster_ASE), ]), points_plot],
#                    cluster_cols = F, cluster_rows = F, show_rownames = F, show_colnames = F, 
#                    annotation_row = clustering_df, scale = 'row', 
#                    annotation_col = column_annotation, 
#                    gaps_row = cumsum(table(clustering_df$Cluster_Expression)),
#                    color = colfunc(10))
# dev.off()
# 
# pdf("~/Desktop/PaperWriting/Fig3/ASE_matrix.pdf")
# pheatmap::pheatmap(gp_latent_matrix_pt[rownames(clustering_df[order(clustering_df$Cluster_ASE), ]), ],
#                    cluster_cols = F, cluster_rows = F, show_rownames = F, show_colnames = F, 
#                    annotation_row = clustering_df, scale = 'row', 
#                    annotation_col = data.frame(column_annotation_collapsed), 
#                    gaps_row = cumsum(table(clustering_df$Cluster_ASE)),
#                    color = colfunc(10))
# dev.off()
# 
# column_annotation_collapsed_2 <- column_annotation_collapsed
# 
# scale_min_max <- function(x){
#   return((x - min(x)) / (max(x) - min(x)))
# }
# 
# column_annotation_collapsed_2$Total_Expression <- scale_min_max(column_annotation_collapsed_2$Total_Expression)
# column_annotation_collapsed_2$Protamine_Expression <-  scale_min_max(column_annotation_collapsed_2$Protamine_Expression)
# 
# ggplot(data.frame(column_annotation_collapsed_2)) + 
#   geom_smooth(aes(Pseudotime, Total_Expression), method = "loess", span = 0.2, col = "grey", size = 2) + 
#   geom_smooth(aes(Pseudotime, Protamine_Expression), method = "loess", span = 0.2, col = "#3cb371", size = 2) + 
#   theme_classic() + 
#   theme(text = element_text(size = 30), 
#         axis.ticks.x = element_blank(), 
#         axis.text.x = element_blank()) + 
#   scale_x_continuous(expand = c(0, 0)) +
#   scale_y_continuous(expand = c(0, 0)) +
#   scale_y_continuous(breaks = c(0, 1), limits = c(-0.05, 1.1)) + 
#   xlab("") + ylab("") + 
#   geom_text(x = 40, y = 0.9, label = "Total expression", size = 10, col = "grey") + 
#   geom_text(x = 80, y = 0.1, label = "Protamine expression", size = 10, col = "#3cb371") + 
#   theme(plot.margin=grid::unit(c(0.1,0.1,0.1,0.1), "mm"))
# 
# ggsave("~/Desktop/PaperWriting/Fig3/Top_of_ase_matrix.pdf", width = 10, height = 1.5)
# 
# pheatmap::pheatmap(gp_latent_matrix[rownames(clustering_df[order(clustering_df$Cluster_ASE), ]), ],
#                    cluster_cols = F, cluster_rows = T, show_rownames = F, show_colnames = F, 
#                    annotation_row = clustering_df, 
#                    annotation_col = column_annotation, 
#                    gaps_row = cumsum(table(clustering_df$Cluster_Expression)),
#                    color = colfunc(10))
# 
# test <- gp_latent_matrix[rownames(clustering_df[order(clustering_df$Cluster_ASE), ]), points_plot]
# test <- t(apply(test,1 , function(x){(x - 0)}))
# #test[test > 0.2] <- 0.2
# #test[test < -0.2] <- -0.2
# 
# clustering_df$mean_ase <- rowMeans(gp_latent_matrix)
# 
# pheatmap::pheatmap((test),
#                    cluster_cols = F, cluster_rows = F, show_rownames = F, show_colnames = F, 
#                    annotation_row = clustering_df, 
#                    annotation_col = column_annotation, 
#                    gaps_row = cumsum(table(clustering_df$Cluster_Expression)))
# 
# pheatmap::pheatmap((test),
#                    cluster_cols = F, cluster_rows = F, show_rownames = F, show_colnames = F, 
#                    annotation_row = clustering_df, scale = "row", 
#                    annotation_col = column_annotation, 
#                    gaps_row = cumsum(table(clustering_df$Cluster_Expression)))

