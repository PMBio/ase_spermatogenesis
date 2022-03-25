# Idea from https://github.com/willtownes/scrna2019
# read data

# seurat also uses similar ideas
# https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1874-1

library(scran)
library(scater)
library(viridis)
library(pheatmap)
library(zoo)
source("~/Desktop/PhD/RCode/Multipurpose/auxiliary.R")
source("~/Desktop/PhD/Projects/ASE_Spermatogenesis/Scripts/ase_functions.R")

setwd("~/Desktop/PhD/Projects/ASE_Spermatogenesis/")

data <- readRDS("./ProcessedData/sce_merged_new_sparse.rds")

#data <- data[,data$CellType %in% c("ES", "RS", "SC", "SG")]
data_f1 <- data[,data$Library %in% c("Sample5", "Sample6")]
data_f1 <- data_f1[!is.na(rowData(data_f1)$chromosome_name), ]

# get genes with high counts in both reference and alternative
total_reference <- rowSums(counts_reference(data_f1))
total_alternative <- rowSums(counts_alternative(data_f1))
genes_check <- total_reference > 1000 & total_alternative > 1000 & !rowData(data_f1)$chromosome_name %in% c("MT", "X")
#genes_check <- total_reference + total_alternative > 0

counts_ref <- counts_reference(data_f1[genes_check,1:1000])
counts_alt <- counts_alternative(data_f1[genes_check,1:1000])
counts_total <- counts_ref + counts_alt

input_ref = counts_ref[rowSums(counts_ref) > 0, ]

my_binom_deviance <- function(X, N, p_null = 0.5){
  # implement deviance score for binomial data
  # we see that for np = const, n --> inf the score grows and for k/n ~ p_null it is smallest
  output = X * log(X / (N * p_null)) + (N - X) * log((N - X)/(N*(1-p_null)))
  output[is.na(output)] = 0 # x * log(x) --> 0
  output
}

#X = outer(1:10, 1:10)
#N = matrix(rep(100, 100), nrow = 10)

#pheatmap(my_binom_deviance(X, N), cluster_rows = F, cluster_cols = F)

deviance_data <- my_binom_deviance(counts_ref, counts_total)
deviance_data <- deviance_data[rowSums(deviance_data) > 0, ]

## look at deviance scores across genes

mean_deviance <- sort(rowMeans(deviance_data), decreasing = T)
deviance_detected <- sort(apply(deviance_data, 1, function(x){sum(x > 0)}), decreasing = T)
top_deviant_features <- names(deviance_detected[1:1000])

ggplot2::qplot(mean_deviance, deviance_detected) + scale_x_log10()

hist(log10(mean_deviance + 1), breaks = 100)

pca_output <- prcomp(t(deviance_data[top_deviant_features, ]))
pca_output_x <- pca_output$x[,1:20]
pca_output_x <- pca_output_x[!duplicated(pca_output_x), ]
tsne_rates <- Rtsne::Rtsne(pca_output_x)

ggplot(data.frame(
  Tsne1 = tsne_rates$Y[,1], 
  Tsne2 = tsne_rates$Y[,2], 
  CellType = data_f1[genes_check,rownames(pca_output_x)]$CellType
), aes(Tsne1, Tsne2, col = CellType)) + geom_point() + 
  ggtitle("ASE, top 1000 features")

ggplot(data.frame(
  Tsne1 = tsne_rates$Y[,1], 
  Tsne2 = tsne_rates$Y[,2], 
  Total = log(colSums(counts_total)[rownames(pca_output_x)] + 1)
), aes(Tsne1, Tsne2, col = Total)) + geom_point()

ggplot(data.frame(
  PC1 = pca_output_x[,1], 
  PC2 = pca_output_x[,2], 
  CellType = data_f1[genes_check,rownames(pca_output_x)]$CellType
), aes(PC1, PC2, col = CellType)) + geom_point()

ggplot(data.frame(
  Tsne1 = reducedDims(data_f1[,rownames(pca_output_x)])$TSNE[,1], 
  Tsne2 = reducedDims(data_f1[,rownames(pca_output_x)])$TSNE[,2], 
  CellType = data_f1[genes_check,rownames(pca_output_x)]$CellType
), aes(Tsne1, Tsne2, col = CellType)) + geom_point() + 
  ggtitle("Gene expression, top 1000 features")

### worried that this will be biased by celltype-specific expression -- look at constitutively expressed genes

ncells_expressed <- rowSums(counts_total > 0)
variance <- rowVars(as.matrix(counts_total))

ggplot(data.frame(
  ncells_expressed, 
  variance
), aes(ncells_expressed, variance)) + geom_point() + scale_x_log10() + scale_y_log10()

genes_broad_expressed <- names(ncells_expressed[ncells_expressed > 500])
genes_broad_expressed <- intersect(genes_broad_expressed, rownames(deviance_data))

pca_output <- prcomp(t(deviance_data[genes_broad_expressed, ]))
pca_output_x <- pca_output$x[,1:20]
pca_output_x <- pca_output_x[!duplicated(pca_output_x), ]
genes_df <- pca_output$rotation
genes_df <- genes_df[order(rowSums(genes_df[,1:10] ** 2), decreasing = T), ]
head(genes_df[,1:5], n = 20)
tsne_rates <- Rtsne::Rtsne(pca_output_x)

ggplot(data.frame(
  Tsne1 = tsne_rates$Y[,1], 
  Tsne2 = tsne_rates$Y[,2], 
  CellType = data_f1[genes_check,rownames(pca_output_x)]$CellType
), aes(Tsne1, Tsne2, col = CellType)) + geom_point() + 
  ggtitle("ASE, broadly_expressed genes (>50% of cells)")

ggplot(data.frame(
  Tsne1 = tsne_rates$Y[,1], 
  Tsne2 = tsne_rates$Y[,2], 
  Residual = deviance_data["Ubc", rownames(pca_output_x)]
), aes(Tsne1, Tsne2, col = Residual)) + geom_point() + 
  ggtitle("ASE, broadly_expressed genes (>50% of cells)") + scale_color_viridis()

ggplot(data.frame(
  Tsne1 = tsne_rates$Y[,1], 
  Tsne2 = tsne_rates$Y[,2], 
  Residual = colSums(deviance_data[, rownames(pca_output_x)])
), aes(Tsne1, Tsne2, col = Residual)) + geom_point() + 
  ggtitle("ASE, broadly_expressed genes (>50% of cells)") + scale_color_viridis()


#### finally on all cells

counts_ref <- counts_reference(data_f1[genes_check,])
counts_alt <- counts_alternative(data_f1[genes_check,])
counts_total <- counts_ref + counts_alt

input_ref = counts_ref[rowSums(counts_ref) > 0, ]

my_binom_deviance <- function(X, N, p_null = 0.5){
  # implement deviance score for binomial data
  # we see that for np = const, n --> inf the score grows and for k/n ~ p_null it is smallest
  output = X * log(X / (N * p_null)) + (N - X) * log((N - X)/(N*(1-p_null)))
  output[is.na(output)] = 0 # x * log(x) --> 0
  output
}

#X = outer(1:10, 1:10)
#N = matrix(rep(100, 100), nrow = 10)

#pheatmap(my_binom_deviance(X, N), cluster_rows = F, cluster_cols = F)

deviance_data <- my_binom_deviance(counts_ref, counts_total)
deviance_data <- deviance_data[rowSums(deviance_data) > 0, ]

## look at deviance scores across genes

mean_deviance <- sort(rowMeans(deviance_data), decreasing = T)
var_deviance <- sort(rowVars(as.matrix(deviance_data)), decreasing = T)
top_deviant_features <- names(mean_deviance[1:length(mean_deviance)])

hist(log10(mean_deviance + 1), breaks = 100)

pca_output <- prcomp(t(deviance_data[top_deviant_features, ]), center = T, scale = T)
pca_output_x <- pca_output$x[,1:100]
pca_output_x <- pca_output_x[!duplicated(pca_output_x), ]
tsne_rates <- Rtsne::Rtsne(pca_output_x)
umap_rates <- umap::umap(pca_output_x)

genes_df <- pca_output$rotation
genes_df <- genes_df[order(rowSums(genes_df[,1:10] ** 2), decreasing = T), ]
head(genes_df[,1:5], n = 30)

ggplot(data.frame(
  Umap1 = umap_rates$layout[,1], 
  Umap2 = umap_rates$layout[,2], 
  Expression = log(counts_total["Cox8a", rownames(pca_output_x)] + 1)
), aes(Umap1, Umap2, col = Expression)) + geom_point(aes(alpha = Expression)) + 
  ggtitle("ASE, broadly_expressed genes (>50% of cells)") + scale_color_gradient(low = "grey", high = "red")

ggplot(data.frame(
  Umap1 = umap_rates$layout[,1], 
  Umap2 = umap_rates$layout[,2], 
  Residual = log(deviance_data["Cox8a", rownames(pca_output_x)] + 1)
), aes(Umap1, Umap2, col = Residual)) + geom_point(aes(alpha = Residual)) + 
  ggtitle("ASE, broadly_expressed genes (>50% of cells)") + scale_color_gradient(low = "grey", high = "red")

ggplot(data.frame(
  Dim1 = umap_rates$layout[,1], 
  Dim2 = umap_rates$layout[,2], 
  CellType = data_f1[genes_check,rownames(pca_output_x)]$CellType
), aes(Dim1, Dim2, col = CellType)) + geom_point() + 
  ggtitle("ASE, 2000 features, all cells")

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

ggplot(data.frame(
  Tsne1 = tsne_rates$Y[,1], 
  Tsne2 = tsne_rates$Y[,2], 
  CellType = data_f1[genes_check,rownames(pca_output_x)]$CellType
), aes(Tsne1, Tsne2, col = CellType)) + geom_point() + 
  ggtitle("ASE, 2000 features, all cells")

total_deviance <- data.frame(
  TotalDeviance = colMeans(deviance_data),
  TotalExpression = colMeans(counts_total), 
  CellType= data_f1$CellType
)

ggplot(total_deviance, aes(TotalDeviance, TotalExpression, col = CellType)) + geom_point() + scale_x_log10() + scale_y_log10()
ggplot(total_deviance, aes(CellType, TotalExpression, col = CellType)) + geom_boxplot() + scale_y_log10() + 
  ggtitle("Library Size")
ggplot(total_deviance, aes(CellType, TotalDeviance, col = CellType)) + geom_boxplot() + scale_y_log10() + 
  ggtitle("Total Deviance")

#### 

gene = 'Pou5f2'
counts_ref_gene <- counts_ref[gene,]
counts_total_gene <- counts_total[gene,]
deviance_gene <- deviance_data[gene, ]

test_data_1 <- data.frame(
  TotalExp = counts_total_gene, 
  RefExp = counts_ref_gene, 
  Deviance = deviance_gene
)

p1 <- ggplot(test_data_1, aes(log(TotalExp + 1), log(RefExp + 1))) + geom_jitter() + ylim(0, 8) + ggtitle(gene) + 
  geom_smooth()
p2 <- ggplot(test_data_1, aes(log(TotalExp + 1), log(Deviance + 1))) + geom_jitter() + ylim(0, 5) + 
  geom_smooth()
gridExtra::grid.arrange(p1, p2, nrow = 2)

# check 



Hmisc::hoeffd(test_data_1[test_data_1$TotalExp > 10, ]$TotalExp, test_data_1[test_data_1$TotalExp > 10, ]$Deviance)
Hmisc::hoeffd(test_data_1$TotalExp, test_data_1$Deviance)



#### look at how top contributing features are wrt seq depth

df_check <- data.frame(
  Gene = rownames(pca_output$rotation), 
  PCA_contribution = rowSums(pca_output$rotation[,1:20] ** 2), 
  SeqDepth = rowSums(counts_total[rownames(pca_output$rotation), ]), 
  InfCells = apply(deviance_data[rownames(pca_output$rotation), ], 1, function(x){sum(x > 0)})
)

df_check <- df_check[order(df_check$PCA_contribution, decreasing = T), ]

ggplot(df_check, aes(SeqDepth, PCA_contribution)) + geom_point() + scale_x_log10() + geom_smooth()
ggplot(df_check, aes(InfCells, PCA_contribution)) + geom_point() + scale_x_log10() + geom_smooth()

### Residuals per celltype

data_gene_plot <- data.frame(
  CellType = data_f1$CellType, 
  Deviance = deviance_data["Rnaseh2a", ]
)

ggplot(data_gene_plot, aes(CellType, log(Deviance + 1), col = CellType)) + geom_boxplot() + geom_jitter(alpha = 0.1)

##

ratios = counts_ref / counts_total
ratios <- abs(ratios - 0.5)

gene = "Aldh1a1"
expression_gene <- counts_total[gene, ]
residuals_gene <- deviance_data[gene, ]
ase_gene <- ratios[gene, ]

ggplot(data.frame("Exp" = expression_gene, "ASE" = ase_gene), 
       aes(Exp, ASE)) + geom_jitter() + geom_smooth()+ scale_x_log10()

ggplot(data.frame("Exp" = expression_gene, "Res" = residuals_gene), 
       aes(Exp, Res)) + geom_jitter() + geom_smooth() + scale_y_log10() + scale_x_log10()

ggplot(data.frame(
  Exp = rowMeans(counts_total), 
  ASE = rowMeans(ratios, na.rm = T), 
  Res = rowMeans(deviance_data)
), aes(Exp, Res)) + geom_point() + scale_x_log10() + scale_y_log10() + geom_smooth()

# fit curve through this
loess.fit <- loess(data = data.frame(Exp = log10(rowMeans(counts_total)), Res = log10(rowMeans(deviance_data))), Res ~ Exp)
loess.fit.resids <- loess.fit$residuals
loess.fit.resids[order(loess.fit.resids, decreasing = T)][1:100]

ratios_fixed <- ratios
ratios_fixed[is.na(ratios_fixed)] <- 0.5

ggplot(data.frame(
  Exp = rowMeans(counts_total), 
  ASE = rowMeans(ratios, na.rm = T), 
  ASE_prior = rowMeans(ratios_fixed),
  Res = rowMeans(deviance_data)
), aes(Exp, ASE)) + geom_point() + scale_x_log10() + geom_smooth()

ggplot(data.frame(
  Exp = rowMeans(counts_total), 
  ASE = rowMeans(ratios, na.rm = T), 
  ASE_prior = rowMeans(ratios_fixed),
  Res = rowMeans(deviance_data)
), aes(Exp, ASE_prior)) + geom_point() + scale_x_log10() + geom_smooth()

##

testy <- data.frame(
  Exp = rowMeans(counts_total), 
  ASE = rowMeans(ratios, na.rm = T), 
  Res = rowMeans(deviance_data)
)

head(testy[order(testy$Res, decreasing = T), ])





####  MOFA stuff #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
data_f1_norm <- logNormCounts(data_f1)
hvgs <- HVG(data_f1_norm)

expression_mofa <- logcounts(data_f1_norm[hvgs, ])
expression_mofa <- t(scale(t(expression_mofa)))
rownames(expression_mofa) <- paste0(rownames(expression_mofa), "_expression")
deviance_mofa <- deviance_data[top_deviant_features[1:1000], ]
deviance_mofa <- apply(deviance_mofa, 2, function(x){x / sum(x)})
deviance_mofa <- t(scale(t(deviance_mofa)))
rownames(deviance_mofa) <- paste0(rownames(deviance_mofa), "_ase")

library(data.table)
library(purrr)
library(MOFA2)
library(reticulate)
use_python("/usr/local/anaconda3/bin/python3.8", required = TRUE)
mofa <- import("mofapy2")
mofa_entrypoint <- mofa$run.entry_point$entry_point()

source("~/Desktop/PhD/RCode/Multipurpose/auxiliary.R")

which_cells <- data_f1$CellType
which_cells <- sample(which(which_cells %in% c("ES", "RS", "SG", "SC")), 5000)

MOFAobject <- create_mofa(list("View1" = as.matrix(expression_mofa[,which_cells]), 
                               "View2" = log(as.matrix(deviance_mofa[,which_cells] + 1))))

plot_data_overview(MOFAobject)

MOFAobject <- prepare_mofa(MOFAobject)
MOFAobject.trained <- run_mofa(MOFAobject, use_basilisk = T)

plot_variance_explained(MOFAobject.trained, max_r2 = 20)
plot_variance_explained(MOFAobject.trained, y="factor", plot_total = T)[[2]]
plot_factor_cor(MOFAobject.trained)

plot_factor(MOFAobject.trained, 
            factor = 1:6)

plot_weights(MOFAobject.trained,
             view = "View2",
             factor = 3,
             nfeatures = 10,     # Number of features to highlight
             scale = T,          # Scale weights from -1 to 1
             abs = F             # Take the absolute value?
)
plot_top_weights(MOFAobject.trained, view = "View1", factor = 5)
plot_top_weights(MOFAobject.trained, view = "View2", factor = 5)

#MOFAobject.trained <- run_tsne(MOFAobject.trained, seed = 222)
MOFAobject.trained <- run_umap(MOFAobject.trained)

plot_variance_explained(MOFAobject.trained, max_r2 = 10)
MOFAobject.trained@samples_metadata <- cbind(MOFAobject.trained@samples_metadata, "CellType" = data_f1$CellType[which_cells])
plot_dimred(MOFAobject.trained,
            method = "UMAP", 
            color_by = "CellType"
)
plot_dimred(MOFAobject.trained,
            method = "UMAP", 
            color_by = "Factor5"
)

# analysis:
# 1) to which extent does ASE influence gene expression factors? (might point to non-independence)
# 2) are there ASE-factors and do they tell us more about structure?
# 3) does including ASE improve clustering?


