# This script generates the supplement for figure 2 that isnt generated in the script for the main figure 2

# Functions
library(ggplot2)
library(viridis)

setwd("~/Desktop/Projects/ASE_Spermatogenesis_Paper_rerun//")

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
summary_data_meanTest <- readRDS("./Data/processed/Dali_LRT_calls.csv")

# exclude genes with mapping effects (see get mapping effects) and genes with no reads
summary_data_meanTest <- summary_data_meanTest[!summary_data_meanTest$Mapping_effect & 
                                                 summary_data_meanTest$coverage_ref > 1 & summary_data_meanTest$coverage_alt > 1, ]

# read results from dynamic tests and add to summary df
results_linear_kernel <- readRDS("./Data/processed/Dali_scoreTest_linear.csv")
results_polynomial_kernel <- readRDS("./Data/processed/Dali_scoreTest_poly.csv")
results_cluster_kernel <- readRDS("./Data/processed/Dali_scoreTest_discrete.csv")

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
data_total_expression <- readRDS("./Data/processed/sce_merged_new.rds")
data_total_expression_save <- data_total_expression
data_total_expression <- data_total_expression[,data_total_expression$AnnotatedClusters %in% c("ES", "RS", "SC", "SG")]
data_total_expression <- data_total_expression[,data_total_expression$Library %in% c("Sample5", "Sample6")]

# show pseudotime for S2.1
pt <- data_total_expression_save$Pseudotime
#pt <- - (pt - max(pt))

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
  scale_color_viridis() + theme_paper(textsize = 40) + facet_wrap(~Species)
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
  scale_color_viridis() + theme_paper(textsize = 40) + scale_fill_manual(values = library_colors) + 
  facet_wrap(~Species) + ylab("Density of cells")
ggsave("./Plots/FigureS2/FigS2_3_pseudotime_distribution_sep.pdf")

