### Load libraries
setwd("~/Desktop/PhD/Projects/ASE_Spermatogenesis_Paper_rerun_Revisions/")

library(scran)
library(scater)
library(viridis)
library(pheatmap)
library(zoo)
library(cowplot)
source("./Scripts/General/auxiliary.R")
source("./Scripts/General/reuse_functions.R")
source("./Scripts/General/ase_functions.R")

## ---------------------- dali functions ---------------------- 

# initialize python path
library(reticulate)
use_python("/Users/jasper/Library/r-miniconda/bin/python", required = T)
dali_path = "~/Desktop/PhD/Projects/ASE_Spermatogenesis/Scripts/dali/"

# define R wrappers around scDALI functions
test_regions_R <- function(A, D, 
                           cell_state, 
                           n_cores = 1L)
{
  source_python(paste0(dali_path, "/dali/my_test_regions.py"))
  test_regions(np_array(A), 
               np_array(D), 
               np_array(cell_state), 
               n_cores = n_cores)
}

test_mean_R <- function(a, d, mean_mu = 0.5)
{
  source_python(paste0(dali_path, "/dali/my_test_mean.py"))
  res = test_mean(np_array(a), 
                  np_array(d), 
                  mean_mu = mean_mu)
  names(res) = c("Mean", "Theta", "NLL_null", "NLL_alt", "pval")
  res
}

run_gp_R <- function(A, D, 
                     cell_state, 
                     kernel = "Linear",
                     n_cores = 1L)
{
  source_python(paste0(dali_path, "/dali/my_run_gp.py"))
  res = run_gp(np_array(A), 
               np_array(D), 
               np_array(cell_state), 
               kernel = kernel,
               n_cores = n_cores)
  res
}


## ---------------------- R functions ---------------------- 

# function to plot ase trajectory for a gene
plot_gene <- function(sce, gene, remove_zero = F){
  data_test <- data.frame(
    pt = sce$Pseudotime,
    exp_ref = as.numeric(counts_reference(sce[gene,])),
    exp_alt = as.numeric(counts_alternative(sce[gene,])),
    ase = as.numeric(allelic_ratios(sce[gene,]))
  )
  if (remove_zero){
    data_test <- data_test[data_test$exp_ref + data_test$exp_alt > 0,]
  }
  data_test$exp_total = data_test$exp_ref + data_test$exp_alt
  data_test <- data_test[order(data_test$pt),]
  
  p1 <- ggplot(data_test) + 
    geom_smooth(aes(pt, exp_ref, color = "B6")) + 
    geom_smooth(aes(pt, exp_alt, color = "Cast")) + 
    theme_classic() + xlim(c(0, 4)) + xlab("") + ylab("Expression") +
    theme(legend.position="top") + 
    scale_color_manual(values = c("B6" = "black", "Cast" = "brown"))
  
  p2 <- ggplot(data_test, aes(pt, ase)) +
    ylim(c(0, 1)) + geom_jitter(width = 0.1, height = 0.02, size = 0.1) + 
    scale_color_viridis() + 
    theme_classic() + geom_smooth(color = "purple") + xlim(c(0, 4)) + 
    xlab("Pseudotime") + ylab("Allelic Bias")
  
  cowplot::plot_grid(p1, p2, 
                     nrow = 2, align = "v", 
                     rel_heights = c(0.4, 0.6))
}

# function to plot ase trajectory for a gene
plot_gene_GP <- function(sce, gene, latent_ase, remove_zero = F, scale_var = 1.96){
  data_test <- data.frame(
    pt_here = sce$Pseudotime,
    exp_ref = as.numeric(counts_reference(sce[gene,])),
    exp_alt = as.numeric(counts_alternative(sce[gene,])),
    latent_ase = latent_ase$posterior_mean,
    latent_var_lower = latent_ase$posterior_mean - scale_var * sqrt(latent_ase$posterior_var),
    latent_var_upper = latent_ase$posterior_mean + scale_var * sqrt(latent_ase$posterior_var)
  )
  data_test$ase <- data_test$exp_ref / (data_test$exp_alt + data_test$exp_ref)
  
  if (remove_zero){
    data_test <- data_test[data_test$exp_ref + data_test$exp_alt > 0,]
  }
  
  data_test$exp_total = data_test$exp_ref + data_test$exp_alt
  data_test <- data_test[order(data_test$pt),]
  
  print(summary(data_test$pt_here))
  
  p1 <- ggplot(data_test) + 
    geom_smooth(aes(pt_here, exp_ref, color = "B6")) + 
    geom_smooth(aes(pt_here, exp_alt, color = "Cast")) + 
    theme_classic() + xlim(c(0, max(data_test$pt_here + 1))) + xlab("") + ylab("Expression") +
    theme(legend.position="top") + 
    scale_color_manual(values = c("B6" = "black", "Cast" = "brown"))
  
  p2 <- ggplot(data_test, aes(pt_here, ase)) +
    ylim(c(0, 1)) + geom_point(width = 0.1, height = 0.02, size = 0.1) + 
    scale_color_viridis() + 
    geom_line(aes(pt_here, latent_ase), color = "green") + 
    geom_ribbon(aes(x = pt_here, ymin = latent_var_lower, ymax = latent_var_upper), 
                color = "green", alpha = 0.2) + 
    theme_classic() + xlim(c(0, max(data_test$pt_here + 1))) + 
    xlab("Pseudotime") + ylab("Allelic Bias")
  
  cowplot::plot_grid(p1, p2, 
                     nrow = 2, align = "v", 
                     rel_heights = c(0.4, 0.6))
  
}

counts_reference <- function(sce){
  assays(sce)[["counts_reference"]]
}

counts_alternative <- function(sce){
  assays(sce)[["counts_alternative"]]
}

allelic_ratios <- function(sce){
  assays(sce)[["allelic_ratios"]]
}

compute_quantile_diff <- function(d, q = 0.05){
  as.numeric(abs(quantile(d, q) - quantile(d, 1 - q)))
}

theme_tsne <- function(){
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_line(color = "grey"),
        plot.background=element_blank())
}

theme_paper <- function(){
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_line(color = "grey"),
        plot.background=element_blank())
}

## ------------------------------------------------------- 
## ---------------------- read data ---------------------- 
## ------------------------------------------------------- 

data <- readRDS("./Data/processed/sce_merged_new.rds")

# subset on germ cells and F1
data <- data[,data$CellType %in% c("ES", "RS", "SC", "SG")]
data_f1 <- data[,data$Library %in% c("Sample5", "Sample6")]
data_f1$Species <- list(
  "Sample1" = "B6",
  "Sample2" = "B6",
  "Sample3" = "CAST",
  "Sample4" = "CAST",
  "Sample5" = "CASTxB6",
  "Sample6" = "CASTxB6"
)[data_f1$Library]

# get gene statistics to see which ones make sense to even analyze
data_f1 <- data_f1[,data_f1$CellType %in% c("SG", "SC", "RS", "ES")]

summary_data <- data.frame(
  coverage_total = rowSums(counts(data_f1)),
  coverage_ref = rowSums(counts_reference(data_f1)),
  coverage_alt = rowSums(counts_alternative(data_f1))
)

# only use genes with at least 1000 allelic reads
summary_data$genes_use <- summary_data$coverage_ref + summary_data$coverage_alt > 1000 & 
  summary_data$coverage_total > 10

genes.test <- summary_data[summary_data$genes_use,]
genes.show <- rownames(genes.test)

#### Analysis scDALI
## Part1: LLR test if gene has allelic bias

genes_test <- rownames(genes.test)
data_test_ref <- counts_reference(data_f1[genes_test,])
data_test_alt <- counts_alternative(data_f1[genes_test,])

## here we actually run the tests - comment out if needed
#run beta-binomial test to detect genes with persistent ASE
# llr_pvals <- data.frame(do.call("rbind", lapply(genes_test, function(gene){
#  print(gene)
#  a = data_test_ref[gene,]
#  d = data_test_alt[gene,] + data_test_ref[gene,]
#  pval = test_mean_R(a, d)
#  return(pval)
# })))

# rownames(llr_pvals) <- genes_test
# summary_data_meanTest <- cbind(summary_data[genes_test,], llr_pvals[genes_test, ])
# summary_data_meanTest[is.na(summary_data_meanTest$pval),]$pval <- 1
# summary_data_meanTest$padj <- p.adjust(summary_data_meanTest$pval)
# summary_data_meanTest$ase <- summary_data_meanTest$coverage_ref /
#  (summary_data_meanTest$coverage_alt + summary_data_meanTest$coverage_ref)

## find mapping effects
# data_parental <- readRDS("./Data/processed/sce_merged_new.rds")
# data_f1 <- annotate_chromosome_sce(data_f1)
# av.counts.b6 <- rowSums(counts(data_parental[,data_parental$Library %in% c("Sample1", "Sample2")]))
# av.counts.cast <- rowSums(counts(data_parental[,data_parental$Library %in% c("Sample3", "Sample4")]))
# sfs <- c(sum(av.counts.b6), sum(av.counts.cast)) / sum(av.counts.b6)
# 
# summary_parental <- data.frame(
#   CountsB6 = av.counts.b6 / sfs[[1]],
#   CountsCast = av.counts.cast / sfs[[2]],
#   LogFC = log(av.counts.b6 + 1) - log(av.counts.cast + 1)
# )
# 
# summary_parental <- summary_parental[rownames(summary_data_meanTest),]
# summary_parental$LogFC_filial <- log(summary_data_meanTest$coverage_ref + 1) - log(summary_data_meanTest$coverage_alt + 1)
# summary_parental$Chromosome <- rowData(data_f1)[rownames(summary_parental),]$chromosome_name
# summary_parental$Chromosome_X_MT <- summary_parental$Chromosome %in% c("X", "MT")
# 
# FC_cutoff <- 8
# df_genes_show <- summary_parental[(abs(summary_parental$LogFC) > FC_cutoff |
#                                      abs(summary_parental$LogFC_filial) > FC_cutoff),]
# 
# genes_exclude <- abs(summary_parental$LogFC) < 1.5 & summary_parental$LogFC_filial < -3
# 
# summary_data_meanTest$Mapping_effect <- genes_exclude

#saveRDS(summary_data_meanTest, "./Data/processed/Dali_LRT_calls.csv")

summary_data_meanTest <- readRDS("./Data/processed/Dali_LRT_calls.csv")
summary_data_meanTest$Mean <- unlist(summary_data_meanTest$Mean)
summary_data_meanTest$Theta <- unlist(summary_data_meanTest$Theta)
summary_data_meanTest$NLL_null <- unlist(summary_data_meanTest$NLL_null)
summary_data_meanTest$NLL_alt <- unlist(summary_data_meanTest$NLL_alt)
summary_data_meanTest$pval <- unlist(summary_data_meanTest$pval)

#### Analysis dynamic scASE
## Part1:  Run score test with linear kernel

genes_test <- summary_data_meanTest[!summary_data_meanTest$Mapping_effect,]
genes.test <- genes_test[genes_test$coverage_ref > 1 & genes_test$coverage_alt > 1,]

A <- t(counts_reference(data_f1[rownames(genes.test),]))
D <- t(counts_reference(data_f1[rownames(genes.test),]) + 
       counts_alternative(data_f1[rownames(genes.test),]))
pseudotime <- data_f1$Pseudotime

# results_linear_kernel <- test_regions_R(A, D, pseudotime, n_cores = 4L)
# saveRDS(results_linear_kernel, "./Data/processed/Dali_scoreTest_linear.csv")
results_linear_kernel <- readRDS("./Data/processed/Dali_scoreTest_linear.csv")

genes.test$dali_pval_linear <- results_linear_kernel
genes.test$chromosome <- rowData(data_f1)[rownames(genes.test),]$chromosome_name

table(genes.test$dali_pval_linear < .01)

# we also test for dependence on celltype specifically using a celltype-kernel

make_onehot_cluster <- function(clusters){
  mm = do.call('cbind', lapply(unique(clusters), function(x){as.numeric(x == clusters)}))
  colnames(mm) = unique(clusters)
  rownames(mm) = names(clusters)
  mm
}

cluster_kernel <- make_onehot_cluster(data_f1$CellType)

# results_cluster_kernel <- test_regions_R(A, D, cluster_kernel, n_cores = 4L)
# names(results_cluster_kernel) <- rownames(genes.test)
# saveRDS(results_cluster_kernel, "./Data/processed/Dali_scoreTest_discrete.csv")

make_polynomial_kernel <- function(cell_state, degree = 1){
  cell_state = cell_state / max(cell_state)
  t(do.call("rbind", lapply(1:degree, function(i){
    cell_state ** i
  })))
}

# 
# pseudotime_polynomial <- make_polynomial_kernel(pseudotime, degree = 3)

# results_polynomial_kernel <- test_regions_R(A, D, pseudotime_polynomial)
# saveRDS(results_polynomial_kernel, "./Data/processed/Dali_scoreTest_poly.csv")
# names(results_polynomial_kernel) <- colnames(A)

#results_polynomial_kernel <- readRDS("./Data/processed/Dali_scoreTest_poly.csv")
#results_polynomial_kernel[is.na(results_polynomial_kernel)] <- 1
#genes.test$dali_pval_polynomial <- results_polynomial_kernel

#saveRDS(genes.test, "./Data/processed/Dali_full_results.rds")

genes.test <- readRDS("./Data/processed/Dali_full_results.rds")

### ------ downstream analysis --------------------------
### calculate smoothed profiles for each significant gene

genes_sig <- rownames(genes.test[p.adjust(genes.test$dali_pval_polynomial) < 0.01,])
genes_sig <- rownames(genes.test)

A_here <- A[, genes_sig]
D_here <- D[, genes_sig]

pseudotime_poly <- make_polynomial_kernel(data_f1$Pseudotime, degree = 3)

# run GP over all, doing it gene by gene makes it work even if some fail

sample_indices <- sample(1:nrow(A_here), nrow(A_here))
sample_indices <- sample_indices[order(sample_indices)]

gp_results_total <- lapply((1:ncol(A_here)), function(i){
  result = tryCatch({
    print(paste0("Current gene: ", colnames(D_here)[[i]]))
    if (sum(D_here[sample_indices, i]) < 1000){
      return(list(NA, NA))
    }
    gp_results_total <- run_gp_R(A_here[sample_indices, i],
                                 D_here[sample_indices, i],
                                 cell_state = pseudotime[sample_indices],
                                 kernel = "RBF")

    gp_results_total

  }, error = function(error_condition) {
    return(list(NA, NA))
  })
})

saveRDS(gp_results_total, "./Data/processed/GP_fits_all.rds")
saveRDS(gp_results_total, "./Data/processed/GP_fits.rds")

gp_results_total <- readRDS("./Data/processed/GP_fits.rds")




