setwd("~/Desktop/Projects/ASE_Spermatogenesis_Paper/")

library(scran)
library(scater)
library(viridis)
library(pheatmap)
library(zoo)
library(tidyverse)
source("./Scripts/General/auxiliary.R")
source("./Scripts/General//reuse_functions.R")

data_evo <- readRDS("./data/processed/final_sce_evo_dataset.rds")
colnames(data_evo) <- paste0(data_evo$Library, "_", data_evo$Barcode)
corrected_data <- metadata(data_evo)$corrected
rownames(corrected_data) <- colnames(data_evo)
data_evo <- data_evo[,!data_evo$CellType %in% c("Immune", "Leydig", "Sertoli")]
corrected_data <- corrected_data[colnames(data_evo), ]
dim(corrected_data)

# load other dataset for reference 
data_ref <- readRDS("~/Desktop/Projects/ASE_Spermatogenesis_Paper/Data/processed/data_with_pseudotime.rds")
data_ref <- logNormCounts(data_ref)
data_ref <- data_ref[,!data_ref$CellType %in% c("Immune", "Leydig", "Sertoli")]

# combine sces
genes_both <- intersect(rownames(data_ref), rownames(data_evo))

coldata_here <- colData(data_evo)[,c("Sample", "Library")]
colnames(coldata_here) <- c("Sample", "Library")
to_combine_evo <- SingleCellExperiment(
  assays = list("counts" = counts(data_evo)), 
  colData = coldata_here
)

coldata_here <- colData(data_ref)[,c("Species", "Library")]
colnames(coldata_here) <- c("Sample", "Library")
to_combine_ref <- SingleCellExperiment(
  assays = list("counts" = counts(data_ref)), 
  colData = coldata_here
)

split.sce <- function(sce, groups, colData.name = "SubCluster"){
  # List to collect individual single cell experiments
  list.out <- list()
  for(i in groups){
    print(i)
    cur_sce <- sce[,as.character(colData(sce)[[colData.name]]) == as.character(i)]
    cur_sce <- computeSumFactors(cur_sce)
    cur_sce <- logNormCounts(cur_sce)
    list.out[[i]] <- cur_sce
  }
  names(list.out) <- groups
  list.out
}

data_combine <- cbind(to_combine_evo[genes_both, ], to_combine_ref[genes_both, ])
data_combine <- logNormCounts(data_combine)

set.seed(123)
data_combine_split <- split.sce(data_combine, unique(data_combine$Library), "Library")
batch_corrected <- batch.correction(data_combine_split)
batch_corrected_knn_graph <- buildSNNGraph(t(batch_corrected))
princomp_here <- prcomp(batch_corrected)

experiments <- c(rep("exp1", ncol(to_combine_evo)), rep("exp2", ncol(to_combine_ref)))

# find k NNs and get average pseudotime for each cell
ref_cells <- t(princomp_here$x[experiments == "exp2", 1:50])
query_cells <- t(princomp_here$x[experiments == "exp1", 1:50])
pseudotime_ref <- data_ref$Pseudotime

k = 50 # 50 was good
cors <- cor(ref_cells, query_cells)
nns <- apply(cors, 2, function(x){order(x, decreasing = T)[1:k]})
pt_assigned <- apply(nns, 2, function(x){mean(pseudotime_ref[x])}) # median was good

ggplot() + 
  geom_histogram(data = data.frame(x = pseudotime_ref), aes(x), fill = "red", bins = 100, alpha = 0.3, col = "grey") + 
  geom_histogram(data = data.frame(x = c(1, pt_assigned), 
                                   species = factor(c("REF", data_evo$Sample), levels = c("REF", "B6", "CAST", "CAROLI"))), 
                 aes(- (x - max(x)), fill = species), bins = 100, alpha = 1) + 
  theme_classic() + xlab("Pseudotime") + ylab("Number of Cells") + 
  theme(text = element_text(size = 30)) + 
  scale_fill_manual(values = c("REF" = "grey", "B6" = "black", "CAST" = "chocolate", "CAROLI" = "brown")) + 
  facet_wrap(~species, nrow = 1)
ggsave("./Plots/FigureS5/FigS5_7.pdf")

# check pseudotime on dimensionality reduction
data.frame(
  umap1 = reducedDims(data_evo)$UMAP[,1], 
  umap2 = reducedDims(data_evo)$UMAP[,2],
  species = factor(data_evo$Sample, levels = c("B6", "CAST", "CAROLI")), 
  pseudotime = pt_assigned[colnames(data_evo)]
) %>% ggplot(aes(umap1, umap2, col = - (pseudotime - max(pseudotime)))) + 
  geom_point() + facet_wrap(~species) + 
  theme_tsne() + scale_color_viridis() + 
  theme(text = element_text(size = 30)) + labs(color="Pseudotime")
ggsave("./Plots/FigureS5/FigS5_8.pdf", width = 15, height = 8)

pt <- pt_assigned

n_partitions <- 100
pseudotime_here <- pt
pseudotime_here <- max(pseudotime_here) - pseudotime_here
pseudotime_here <- pseudotime_here / max(pseudotime_here) # between 0 and 1
intervals <- c(0, 1:n_partitions / n_partitions)

make_pseudotime_smooth <- function(data, pt, n_partitions = 100){
  pseudotime_here <- pt
  #pseudotime_here <- pseudotime_here - min(pseudotime_here)
  pseudotime_here <- pseudotime_here / max(pseudotime_here) # between 0 and 1
  intervals <- c(0, 1:n_partitions / n_partitions)
  
  expression_smoothed <- do.call("cbind", 
                                 lapply(1:n_partitions, function(i){
                                   p_lower <- intervals[i]
                                   p_upper <-  intervals[i + 1]
                                   ix <- which(p_lower <= pseudotime_here & pseudotime_here <= p_upper)
                                   if (length(ix) == 0){
                                     return(rep(0, nrow(data)))
                                   }
                                   if (length(ix) == 1){
                                     print("test")
                                     return(data[,ix])
                                   }
                                   rowMeans(data[,ix])
                                 }))
  expression_smoothed
}

data_f0_b6 <- make_pseudotime_smooth(counts(data_evo[,data_evo$Sample == "B6"]),
                                     pseudotime_here[data_evo$Sample == "B6"], n_partitions = 100)
data_f0_cast <- make_pseudotime_smooth(counts(data_evo[,data_evo$Sample == "CAST"]),
                                     pseudotime_here[data_evo$Sample == "CAST"], n_partitions = 100)
data_f0_caroli <- make_pseudotime_smooth(counts(data_evo[,data_evo$Sample == "CAROLI"]),
                                         pseudotime_here[data_evo$Sample == "CAROLI"], n_partitions = 100)

total_counts_per_gene <- rowSums(counts(data_evo))
genes_test <- names(total_counts_per_gene[total_counts_per_gene > 1000])

data_f0_b6 <- data_f0_b6[genes_test, ]
data_f0_cast <- data_f0_cast[genes_test, ]
data_f0_caroli <- data_f0_caroli[genes_test, ]

# save for fitting 
saveRDS(list(data_f0_b6, data_f0_cast, data_f0_caroli), file = "~/Desktop/Projects/ASE_Spermatogenesis_Paper/Data/processed/data_evo_for_fitting.rds")



