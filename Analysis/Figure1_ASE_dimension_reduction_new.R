library(scran)
library(scater)
library(viridis)
library(pheatmap)
library(zoo)
source("./Scripts/General/auxiliary.R")
source("./Scripts/General/reuse_functions.R")

## Read data and subset on F1 samples
setwd("~/Desktop/PhD/Projects/ASE_Spermatogenesis_Paper_rerun//")

data <- readRDS("./Data/processed/sce_merged_new.rds")
data <- annotate_chromosome_sce(data)

# amount of reads with allele-specific signal?

data.frame(
  counts_total = colSums(counts(data)), 
  counts_ae = colSums(counts_reference(data)) + colSums(counts_alternative(data)), 
  sample = data$Library
) %>% group_by(sample) %>%
  summarize(counts_total = sum(counts_total), counts_ae = sum(counts_ae)) %>%
  add_column(ratio = .$counts_ae / .$counts_total)

#data <- data[,data$CellType %in% c("ES", "RS", "SC", "SG")]
data_f1 <- data[,data$Library %in% c("Sample5", "Sample6")]
data_f1 <- data_f1[!is.na(rowData(data_f1)$chromosome_name), ]

# get genes with high counts in both reference and alternative
ref_data <- counts_reference(data_f1)
total_data <- counts_reference(data_f1) + counts_alternative(data_f1)
total_reference <- rowSums(counts_reference(data_f1))
total_alternative <- rowSums(counts_alternative(data_f1))
genes_check <- total_reference > 500 & total_alternative > 500 & !rowData(data_f1)$chromosome_name %in% c("MT", "X")
genes_use <- names(which(genes_check))
#genes_check <- total_reference + total_alternative > 0

# fit BB model to look for genes with high over-dispersion

fit_bb_model <- function(gene){
  tryCatch({
    k = ref_data[gene, ]
    N = total_data[gene, ]
    ll <- function(alpha, beta) {
      -sum(dbetabinom.ab(k, N, alpha, beta, log = TRUE))
    }
    m <- mle(ll, start = list(alpha = 1, beta = 10), method = "L-BFGS-B")
    coefs = coef(m)
    alpha = coefs[[1]]
    beta = coefs[[2]]
    c(gene, "mean" = alpha / (alpha + beta), "rho" = 1 / (alpha + beta + 1))
  }, error=function(cond) {
    return(c(gene, NA, NA))
  })
}

bb_fits <- lapply(genes_use, fit_bb_model)
saveRDS(bb_fits, "~/Desktop/BB_fits_dimred.rds")
#bb_fits <- readRDS("~/Desktop/BB_fits_dimred.rds")
bb_fits_df <- data.frame(do.call("rbind", bb_fits))
bb_fits_df$Gene <- genes_use
bb_fits_df$mean <- as.numeric(bb_fits_df$mean)
bb_fits_df$rho <- as.numeric(bb_fits_df$rho)
bb_fits_df <- bb_fits_df[!is.na(bb_fits_df$mean), ]
bb_fits_df <- bb_fits_df[order(bb_fits_df$rho, decreasing = T), ]

head(bb_fits_df, n = 50)

genes_check <- bb_fits_df[bb_fits_df$rho > 0.05, ]$Gene

counts_ref <- counts_reference(data_f1[genes_check, ])
counts_alt <- counts_alternative(data_f1[genes_check, ])
counts_total <- counts_ref + counts_alt

gene = "Tex21"
data.frame(
  celltype = data_f1$CellType, 
  value = counts_ref[gene, ] / (counts_ref[gene, ] + counts_alt[gene, ])
) %>% ggplot(aes(x = celltype, y = value)) +
  geom_violin() + 
  stat_summary()

ratios_stabilized <- logit((counts_ref + 1) / (counts_ref + counts_alt + 2))
pca_here <- prcomp(t(ratios_stabilized))
pca_here_genes <- data.frame(pca_here$rotation)

ggplot(data.frame(
  PC1 = pca_here$x[,1], 
  PC2 = pca_here$x[,2], 
  celltype = data_f1$CellType
), aes(PC1, PC2, col = celltype)) + geom_point()

umap_here <- umap::umap(pca_here$x[,1:50])

ggplot(data.frame(
  umap1 = umap_here$layout[,1], 
  umap2 = umap_here$layout[,2], 
  celltype = data_f1$CellType
), aes(umap1, umap2, col = celltype)) + geom_point(size = 3)

input_ref = counts_ref[rowSums(counts_ref) > 0, ]

my_binom_deviance <- function(X, N, p_null = 0.5){
  # implement deviance score for binomial data
  # we see that for np = const, n --> inf the score grows and for k/n ~ p_null it is smallest
  output = X * log(X / (N * p_null)) + (N - X) * log((N - X)/(N*(1 - p_null)))
  # output[is.na(output)] = 0 # x * log(x) --> 0
  output
}

#N = 100
#plot(my_binom_deviance(0:N, N))

rates = counts_ref / (counts_ref + counts_alt)
rates_stab = (counts_ref + 1) / (counts_ref + counts_alt + 2)
stand_errors = rates_stab * (1 - rates_stab) / (counts_ref + counts_alt + 1)
stand_errors_inv = 1 / stand_errors

xpand <- function(d) do.call("expand.grid", rep(list(1:ncol(d)), 2))
euc_norm <- function(x) sqrt(sum(x^2))

euc_dist <- function(mat, weights) {
  iter <- xpand(mat)
  vec <- mapply(function(i,j) {
    weights = sum(sqrt(weights[, i] * weights[, j]))
    weights = weights / sum(weights)
    sum(weights * (mat[,i] - mat[,j]) ^ 2)
  }, iter[,1], iter[,2])
  matrix(vec,ncol(mat), ncol(mat))
}

n_cells = 100
testy <- euc_dist(rates_stab[, 1:n_cells], stand_errors_inv[, 1:n_cells])
pheatmap(testy)

testy <- cor(as.matrix(rates_stab[,1:500]))

pca_here <- prcomp(testy)$x[,1:50]
umap_here <- umap::umap(pca_here)

data.frame(
  CellType = data$CellType[1:500], 
  umap1 = umap_here$layout[,1], 
  umap2 = umap_here$layout[,2]
) %>% 
  ggplot(aes(umap1, umap2, col = CellType)) + geom_point()

distances_check <- do.call("rbind", lapply(1:ncol(rates)[1:10], function(i){
  unlist(lapply(1:ncol(rates)[1:10], function(j){
    distance_measure(rates_stab[,i], rates_stab[,j], stand_errors_inv[,i], stand_errors_inv[,j])
  }))
}))

distance_measure(rates_stab[,1], rates_stab[,2], stand_errors_inv[,1], stand_errors_inv[,2])



average_ps <- rowMeans(counts_ref / (counts_ref + counts_alt), na.rm = T)

#deviance_data <- my_binom_deviance(counts_ref, counts_total)
deviance_data_relative <- do.call("rbind", lapply(1:nrow(counts_ref), function(i){
  X = counts_ref[i, ]
  N = counts_ref[i, ] + counts_alt[i, ]
  my_binom_deviance(X, N, p_null = average_ps[i])
}))

deviance_data_relative <- deviance_data
rownames(deviance_data_relative) <- rownames(counts_ref)

deviance_data_relative <- deviance_data_relative[,!colnames(deviance_data_relative) %in% c("Sample5_CTTGATTAGTCTAGAA-1")]

# deviance_data_nona <- deviance_data
# deviance_data_nona[is.na(deviance_data_nona)] <- 0
# deviance_data_nona <- deviance_data_nona[rowSums(deviance_data_nona) > 0, ]

# deviance_data[deviance_data > 5] <- 5

# only include features with at least some non-NA values
num_na <- apply(deviance_data_relative, 1, function(x){sum(!is.na(x))})
which_genes <- names( which(num_na > ncol(deviance_data_relative) / 10 ))
deviance_data_relative <- deviance_data_relative[which_genes, ]
rownames(deviance_data_relative) <- which_genes

deviance_data_input <- sign(deviance_data_relative) * sqrt(deviance_data_relative)
deviance_data_input <- t(scale(t(deviance_data_input)))
rownames(deviance_data_input) <- which_genes
# testytest <- prcomp(t(deviance_data), na.action = na.omit)

testytest <- rowQuantiles(deviance_data_relative, na.rm = T)[,4] - rowQuantiles(deviance_data_relative, na.rm = T)[,2]
names(testytest) <- rownames(deviance_data_relative)
sort(testytest, decreasing = T)[1:100]

data.frame(
  celltype = data_f1$CellType, 
  value = deviance_data_input["Ddt", ]
) %>% ggplot(aes(x = celltype, y = value)) +
  geom_jitter() + 
  geom_boxplot()

meta_ases <- data.frame(do.call("cbind", lapply(unique(data_f1$CellType), function(ct){
  xx <- deviance_data_input_scaled[,data_f1$CellType == ct]
  return(rowMeans(xx, na.rm = T))
})))
colnames(meta_ases) <- unique(data_f1$CellType)


## look at independence between gene expression and ase using hoeffers D

# hoeff_d_per_gene <- function(gene, filtering_cutoff = - 1){
#   counts_ref_gene <- counts_ref[gene,]
#   counts_total_gene <- counts_total[gene,]
#   deviance_gene <- deviance_data[gene, ]
#   
#   test_data_1 <- data.frame(
#     TotalExp = counts_total_gene, 
#     RefExp = counts_ref_gene, 
#     Deviance = deviance_gene
#   )
#   
#   test_data_1 <- test_data_1[!is.na(test_data_1$Deviance), ]
#   test_data_1 <- test_data_1[test_data_1$Deviance > 0, ]
#   
#   Hmisc::hoeffd(test_data_1[test_data_1$TotalExp > filtering_cutoff, ]$TotalExp, 
#                 test_data_1[test_data_1$TotalExp > filtering_cutoff, ]$Deviance)$P[1, 2]
# }
# 
# gene = rownames(deviance_data)[[2]]
# gene = "Tnp2"
# counts_ref_gene <- counts_ref[gene,]
# counts_alt_gene <- counts_alt[gene,]
# counts_total_gene <- counts_total[gene,]
# deviance_gene <- deviance_data[gene, ]
# 
# test_data_1 <- data.frame(
#   TotalExp = counts_total_gene, 
#   RefExp = counts_ref_gene, 
#   AltExp = counts_alt_gene, 
#   Deviance = deviance_gene
# )
# 
# test_data_1 <- test_data_1[!is.na(test_data_1$Deviance), ]
# test_data_1 <- test_data_1[!(test_data_1$RefExp == 0 | test_data_1$AltExp == 0), ]
# 
# p1 <- ggplot(test_data_1, aes(log10(TotalExp + 1), log10(RefExp + 1))) + geom_jitter() + ylim(0, 8) + ggtitle(gene) + 
#   geom_smooth()
# p2 <- ggplot(test_data_1, aes(log10(TotalExp + 1), log10(Deviance + 1))) + geom_jitter() + ylim(0, 5) + 
#   geom_smooth()
# gridExtra::grid.arrange(p1, p2, nrow = 2)
# hoeff_d_per_gene(gene)
# 
# testytest <- lm(Deviance ~ TotalExp, test_data_1)

# compute genome-wide

# check with multiple cutoffs
# hoeff_results <- unlist(lapply(rownames(deviance_data)[1:100], function(g){hoeff_d_per_gene(g)}))
# hist(hoeff_results)
# 
# hoeff_results <- unlist(lapply(rownames(deviance_data)[1:100], function(g){hoeff_d_per_gene(g, filtering_cutoff = 5)}))
# hist(hoeff_results)

## use MOFA to compute dimensionality reduction

library(MOFA2)

# choose interesting features
deviance_data_here <- deviance_data_relative[names(sort(testytest, decreasing = T)), ]

# deviance_data_here <- deviance_data_here[rownames(deviance_data_here) %in% c("Hbb-bs", "Cyp17a1")]
# deviance_data_here <- 
n_na_cells <- apply(deviance_data_here, 1, function(x){sum(!is.na(x))})
feature_mean <- apply(deviance_data_here, 1, function(x){mean(x, na.rm = T)})
feature_variance <- apply(deviance_data_here, 1, function(x){sd(x, na.rm = T)})

ggplot(data.frame(ncells = n_na_cells, var = feature_variance), aes(ncells, var)) + geom_point() + scale_x_log10()
ggplot(data.frame(ncells = n_na_cells, mean = feature_mean, var = feature_variance), aes(mean, var)) + geom_point() + scale_x_log10()
ggplot(data.frame(ncells = n_na_cells, mean = feature_mean, var = feature_variance / feature_mean), aes(mean, var)) + geom_point() + scale_x_log10()

# choose variable features

deviance_data_here <- ratios_stabilized

# cells_subset <- data_f1$CellType %in% c("ES", "RS", "SC")
#MOFAobject <- create_mofa(data = list("View1" = as.matrix(deviance_data_here[, cells_subset])))
MOFAobject <- create_mofa(data = list("View1" = as.matrix(deviance_data_here)))
MOFAobject <- prepare_mofa(MOFAobject)
MOFAobject <- run_mofa(MOFAobject)
MOFAobject <- run_umap(MOFAobject)

MOFA2::plot_factor_cor(MOFAobject)

# plot_dimred(MOFAobject, method = "UMAP", color_by = )
plot_weights(MOFAobject, factors = 1:3)

data.frame(
  umap1 = MOFAobject@dim_red$UMAP[,2], 
  umap2 = MOFAobject@dim_red$UMAP[,3],
  celltype = data_f1$CellType
) %>%
  ggplot(aes(umap1, umap2, col = celltype)) + geom_point() + theme_classic()

factor_values <- get_factors(MOFAobject)$group1

data.frame(
  factor1 = factor_values[,1], 
  factor2 = factor_values[,2]
) %>%
  ggplot(aes(factor1, factor2)) + geom_point()

data.frame(
  umap1 = MOFAobject@dim_red$UMAP[,2], 
  umap2 = MOFAobject@dim_red$UMAP[,3],
  factor = factor_values[,3]
) %>%
  ggplot(aes(umap1, umap2, col = factor)) + geom_point() + theme_classic() + scale_color_viridis(limits = c(-5, 10))

features_choose

gene = "Piwil1"
gene = features_choose[15]
print(gene)
data.frame(
  umap1 = MOFAobject@dim_red$UMAP[,2], 
  umap2 = MOFAobject@dim_red$UMAP[,3],
  dev_value =  as.numeric(deviance_data[gene, ])[cells_subset], 
  is_na = is.na(deviance_data[gene, ])
) %>%
  ggplot(aes(umap1, umap2, col = log(dev_value + 1), shape = is_na)) + geom_point() + theme_classic() + scale_color_viridis()

pca_output <- prcomp(t(deviance_data_nona), center = T, scale = T)
pca_output_x <- pca_output$x[,1:100]
pca_output_x <- pca_output_x[!duplicated(pca_output_x), ]
tsne_rates <- Rtsne::Rtsne(pca_output_x)
umap_rates <- umap::umap(pca_output_x)

genes_df <- pca_output$rotation
genes_df <- genes_df[order(rowSums(genes_df[,1:10] ** 2), decreasing = T), ]
head(genes_df[,1:5], n = 30)

colour_vector <- c("SG" = "#046735",
                   "SC" = "#D31281",
                   "RS" = "#282B69",
                   "ES" = "grey60",
                   "Sertoli" = "#643F18",
                   "Leydig" = "#985D25",
                   "Immune" = "#FFDFB2")

ggplot(data.frame(
  Umap1 = umap_rates$layout[,1], 
  Umap2 = umap_rates$layout[,2], 
  CellType = data_f1[genes_check,rownames(pca_output_x)]$CellType
), aes(Umap1, Umap2, col = CellType)) + geom_point() + 
  theme_tsne() + scale_color_manual(values = colour_vector) + 
  ggtitle("Embedding based on allelic rates") + theme(text = element_text(size = 20))

ggsave("~/Desktop/Projects/ASE_Spermatogenesis/Manuscript/Plots/Figure2/umap_ase.pdf")


### try glmpca
library(glmpca)

#create a simple dataset with two clusters
mu<-rep(c(.5,3),each=10)
mu<-matrix(rnorm(100*20),nrow=100)
mu[1:50,1:10] <- mu[1:50,1:10] + 2
mu <- rev_logit(mu)
clust<-rep(c("red","black"),each=10)
Y<-matrix(rbinom(n = prod(dim(mu)), prob = mu , size = 10), nrow=nrow(mu))

#visualize the latent structure
res<-glmpca(Y, L = 2, fam = "binom", )
factors<-res$factors
plot(factors[,1],factors[,2],col=clust,pch=19)

input_data <- counts_ref[,1:500]
input_data <- input_data[rowSums(input_data) > 0, ]
res <- glmpca(input_data, L = 2, fam = "binom")

factors<-res$factors
plot(factors[,1],factors[,2],col=clust,pch=19)

