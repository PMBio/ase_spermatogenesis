library(DropletUtils)
library(scran)
library(reshape2)
library(scater)
library(ggplot2)
library(Rtsne)
library(irlba)
library(RColorBrewer)
library(viridis)
library(pheatmap)

setwd("~/Desktop/Projects/ASE_Spermatogenesis_Paper_rerun/")
source("./Scripts/General/auxiliary.R")

#### Read individual libraries

sce.sample1 <- read10xCounts("./Data/raw/f1_dataset/Sample1/")
colData(sce.sample1)$Sample <- rep("B6", ncol(sce.sample1))
colData(sce.sample1)$Library <- rep("Sample1", ncol(sce.sample1))

sce.sample2 <- read10xCounts("./Data/raw/f1_dataset/Sample2/")
colData(sce.sample2)$Sample <- rep("B6", ncol(sce.sample2))
colData(sce.sample2)$Library <- rep("Sample2", ncol(sce.sample2))

sce.sample3 <- read10xCounts("./Data/raw/f1_dataset/Sample3/")
colData(sce.sample3)$Sample <- rep("CAST", ncol(sce.sample3))
colData(sce.sample3)$Library <- rep("Sample3", ncol(sce.sample3))

sce.sample4 <- read10xCounts("./Data/raw/f1_dataset/Sample4/")
colData(sce.sample4)$Sample <- rep("CAST", ncol(sce.sample4))
colData(sce.sample4)$Library <- rep("Sample4", ncol(sce.sample4))

sce.sample5 <- read10xCounts("./Data/raw/f1_dataset/Sample5/")
colData(sce.sample5)$Sample <- rep("F1_B6_CAST", ncol(sce.sample5))
colData(sce.sample5)$Library <- rep("Sample5", ncol(sce.sample5))

sce.sample6 <- read10xCounts("./Data/raw/f1_dataset/Sample6/")
colData(sce.sample6)$Sample <- rep("F1_B6_CAST", ncol(sce.sample6))
colData(sce.sample6)$Library <- rep("Sample6", ncol(sce.sample6))

#### merge dataset 

sce.all <- cbind(sce.sample1, sce.sample2,
                 sce.sample3, sce.sample4,
                 sce.sample5, sce.sample6)
rm(list = c("sce.sample1", "sce.sample2", "sce.sample3", "sce.sample4", "sce.sample5", "sce.sample6"))

colnames(sce.all) <- paste0(sce.all$Library, "_", sce.all$Barcode)

# calculate number of detected cells per library
cur_stats <- melt(table(colData(sce.all)$Sample, colData(sce.all)$Library))
cur_stats <- cur_stats[cur_stats$value > 0,]
cur_stats <- cur_stats[order(cur_stats$Var1),]
stats.df <- data.frame(row.names = cur_stats$Var2,
                       Sample = cur_stats$Var1,
                       Library = cur_stats$Var2,
                       No_cells = cur_stats$value)

# ---------------------------------------------------------------------------------------------
ggplot(stats.df, aes(x = Library, y = No_cells, fill = Sample)) + 
  geom_bar(stat = "identity") + theme_classic() + 
  coord_flip() + ylim(0, 6000) + ggtitle("Number of cells before filtering")
# ---------------------------------------------------------------------------------------------

QC.metrics <- perCellQCMetrics(sce.all)
colData(sce.all) <- cbind(colData(sce.all), QC.metrics)
sce.all <- sce.all[,colData(sce.all)$detected > 500]
sce.all <- sce.all[,colData(sce.all)$total > 500]
# Remove genes that are not expressed
sce.all <- sce.all[Matrix::rowSums(counts(sce.all)) > 0,]
# Add to stats data frame
cur_stats <- melt(table(colData(sce.all)$Sample, colData(sce.all)$Library))
cur_stats <- cur_stats[cur_stats$value > 0,]
cur_stats <- cur_stats[order(cur_stats$Var1),]
stats.df$AfterFiltering <- cur_stats$value

write.csv(stats.df, "./Data/processed/filtering_stats.csv")

# ---------------------------------------------------------------------------------------------
ggplot(stats.df, aes(x = Library, y = AfterFiltering, fill = Sample)) + 
  geom_bar(stat = "identity") + theme_classic() + 
  coord_flip() + ylim(0, 6000) + ggtitle("Number of cells after filtering")
# ---------------------------------------------------------------------------------------------

rownames(sce.all) <- rowData(sce.all)$Symbol

clusters <- quickCluster(sce.all, method = "igraph",use.ranks=FALSE, min.size = 100)
sce.all <- computeSumFactors(sce.all, clusters=clusters)
sce.all <- logNormCounts(sce.all)

set.seed(123)
HVgenes <- HVG(sce = sce.all, n = 1000)
HVgenes[1:50]

pca <- prcomp_irlba(t(logcounts(sce.all[HVgenes,])), n = 50)
tsne <- Rtsne(pca$x, pca = FALSE)
reducedDims(sce.all)$TSNE <- tsne$Y

batches <- unique(paste(colData(sce.all)$Sample, colData(sce.all)$Library))
col_vector_batch <- vector(length = length(batches))
names(col_vector_batch) <- batches
col_vector_batch["B6 Sample1"] <- colorRampPalette(c("white", "black"))(10)[10]
col_vector_batch["B6 Sample2"] <- colorRampPalette(c("white", "black"))(10)[4]
col_vector_batch["CAST Sample3"] <- colorRampPalette(c("white", "chocolate"))(10)[10]
col_vector_batch["CAST Sample4"] <- colorRampPalette(c("white", "chocolate"))(10)[4]
col_vector_batch["F1_B6_CAST Sample5"] <- colorRampPalette(c("white", "violet"))(10)[10]
col_vector_batch["F1_B6_CAST Sample6"] <- colorRampPalette(c("white", "violet"))(10)[4]

col_vector_species <- vector(length = 3)
names(col_vector_species) <- c("B6", "CAST", "F1")
col_vector_species["B6"] <- colorRampPalette(c("white", "black"))(10)[8]
col_vector_species["CAST"] <- colorRampPalette(c("white", "chocolate"))(10)[8]
col_vector_species["F1"] <- colorRampPalette(c("white", "violet"))(10)[8]

# ---------------------------------------------------------------------------------------------
ggplot(data.frame(tsne1 = reducedDims(sce.all)$TSNE[,1],
                  tsne2 = reducedDims(sce.all)$TSNE[,2],
                  batch = paste(colData(sce.all)$Sample, 
                                colData(sce.all)$Library))) +
  geom_point(aes(tsne1, tsne2, colour = batch), alpha = 0.6) +
  scale_color_manual(values = col_vector_batch)
# ---------------------------------------------------------------------------------------------
sce.single <- split.sce(sce = sce.all, groups = unique(colData(sce.all)$Library), 
                        colData.name = "Library")

# Batch correction
set.seed(123)
corrected <- batch.correction(sce.single)
rm(sce.single)
# Save batch corrected matrix in sce object
metadata(sce.all)$corrected <- corrected

# dimension reductions
set.seed(12345)
tsne <- Rtsne(corrected, pca = FALSE, perplexity = 350)
reducedDims(sce.all)$TSNE <- tsne$Y
umap <- umap::umap(corrected)
reducedDims(sce.all)$UMAP <- umap$layout

# ---------------------------------------------------------------------------------------------

# Batches
ggplot(data.frame(tsne1 = reducedDims(sce.all)$TSNE[,1],
                  tsne2 = reducedDims(sce.all)$TSNE[,2],
                  batch = paste(colData(sce.all)$Sample,
                                colData(sce.all)$Library))) +
  geom_point(aes(tsne1, tsne2, colour = batch)) +
  scale_color_manual(values = col_vector_batch) + 
  ggtitle("Batchcorrected")

# Number genes expressed
ggplot(data.frame(tsne1 = reducedDims(sce.all)$TSNE[,1],
                  tsne2 = reducedDims(sce.all)$TSNE[,2],
                  gene = colData(sce.all)$total)) +
  geom_point(aes(tsne1, tsne2, colour = gene)) + scale_colour_viridis() + 
  ggtitle("Number of genes expressed")

# Plot gene expression
gene = "Stra8"
ggplot(data.frame(tsne1 = reducedDims(sce.all)$TSNE[,1],
                  tsne2 = reducedDims(sce.all)$TSNE[,2],
                  gene = logcounts(sce.all)[gene,])) +
  geom_point(aes(tsne1, tsne2, colour = gene)) + scale_colour_viridis() + 
  ggtitle(gene)

gene = "Hormad1"
ggplot(data.frame(tsne1 = reducedDims(sce.all)$TSNE[,1],
                  tsne2 = reducedDims(sce.all)$TSNE[,2],
                  gene = logcounts(sce.all)[gene,])) +
  geom_point(aes(tsne1, tsne2, colour = gene)) + scale_colour_viridis() + 
  ggtitle(gene)

gene = "Piwil1"
ggplot(data.frame(tsne1 = reducedDims(sce.all)$TSNE[,1],
                  tsne2 = reducedDims(sce.all)$TSNE[,2],
                  gene = logcounts(sce.all)[gene,])) +
  geom_point(aes(tsne1, tsne2, colour = gene)) + scale_colour_viridis() + 
  ggtitle(gene)

gene = "Pou5f2"
ggplot(data.frame(tsne1 = reducedDims(sce.all)$TSNE[,1],
                  tsne2 = reducedDims(sce.all)$TSNE[,2],
                  gene = logcounts(sce.all)[gene,])) +
  geom_point(aes(tsne1, tsne2, colour = gene)) + scale_colour_viridis() + 
  ggtitle(gene)

gene = "Tex21"
ggplot(data.frame(tsne1 = reducedDims(sce.all)$TSNE[,1],
                  tsne2 = reducedDims(sce.all)$TSNE[,2],
                  gene = logcounts(sce.all)[gene,])) +
  geom_point(aes(tsne1, tsne2, colour = gene)) + scale_colour_viridis() + 
  ggtitle(gene)

gene = "Prm1"
ggplot(data.frame(tsne1 = reducedDims(sce.all)$TSNE[,1],
                  tsne2 = reducedDims(sce.all)$TSNE[,2],
                  gene = logcounts(sce.all)[gene,])) +
  geom_point(aes(tsne1, tsne2, colour = gene)) + scale_colour_viridis() + 
  ggtitle(gene)

gene = "Cst12"
ggplot(data.frame(tsne1 = reducedDims(sce.all)$TSNE[,1],
                  tsne2 = reducedDims(sce.all)$TSNE[,2],
                  gene = logcounts(sce.all)[gene,])) +
  geom_point(aes(tsne1, tsne2, colour = gene)) + scale_colour_viridis() + 
  ggtitle(gene)

gene = "Insl3"
ggplot(data.frame(tsne1 = reducedDims(sce.all)$TSNE[,1],
                  tsne2 = reducedDims(sce.all)$TSNE[,2],
                  gene = logcounts(sce.all)[gene,])) +
  geom_point(aes(tsne1, tsne2, colour = gene)) + scale_colour_viridis() + 
  ggtitle(gene)

gene = "Acta2"
ggplot(data.frame(tsne1 = reducedDims(sce.all)$TSNE[,1],
                  tsne2 = reducedDims(sce.all)$TSNE[,2],
                  gene = logcounts(sce.all)[gene,])) +
  geom_point(aes(tsne1, tsne2, colour = gene)) + scale_colour_viridis() + 
  ggtitle(gene)

gene = "Acta2"
ggplot(data.frame(tsne1 = reducedDims(sce.all)$UMAP[,1],
                  tsne2 = reducedDims(sce.all)$UMAP[,2],
                  gene = logcounts(sce.all)[gene,])) +
  geom_point(aes(tsne1, tsne2, colour = gene)) + scale_colour_viridis() + 
  ggtitle(gene)

# ---------------------------------------------------------------------------------------------

# clusters batch corrected data
set.seed(1234)
g <- buildSNNGraph(t(metadata(sce.all)$corrected), k = 4)
set.seed(1234)
clusters <- igraph::cluster_louvain(g)$membership
table(clusters)

##### annotate clusters
marker.genes <- findMarkers(sce.all, clusters)
marker.genes <- lapply(marker.genes, function(x){x[x$summary.logFC > 0,]})

mito.proportion <- colCounts(as.matrix(counts(sce.all[grepl("^mt-", rownames(sce.all)),]))) / 
  colCounts(as.matrix(counts(sce.all)))

ggplot(data.frame(tsne1 = reducedDims(sce.all)$TSNE[,1],
                  tsne2 = reducedDims(sce.all)$TSNE[,2],
                  mito_prop = mito.proportion)) +
  geom_point(aes(tsne1, tsne2, colour = mito_prop)) + 
  scale_color_viridis()

# ---------------------------------------------------------------------------------------------
ggplot(data.frame(tsne1 = reducedDims(sce.all)$TSNE[,1],
                  tsne2 = reducedDims(sce.all)$TSNE[,2],
                  cluster = as.factor(clusters))) +
  geom_point(aes(tsne1, tsne2, colour = cluster))
# ---------------------------------------------------------------------------------------------

p <- ggplot(data.frame(tsne1 = reducedDims(sce.all)$TSNE[,1],
                       tsne2 = reducedDims(sce.all)$TSNE[,2],
                       cluster = as.factor(clusters))) +
  geom_point(aes(tsne1, tsne2, colour = cluster))

plotly::ggplotly(p)

cluster_annotation <- list(
  "1" = "RS",
  "2" = "ES",
  "3" = "SG",
  "4" = "Sertoli",
  "5" = "RS",
  "6" = "RS",
  "7" = "Outliers",
  "8" = "Sertoli",
  "9" = "ES",
  "10" = "ES",
  "11" = "SC",
  "12" = "ES",
  "13" = "Outliers",
  "14" = "ES",
  "15" = "RS",
  "16" = "Outliers",
  "17" = "RS",
  "18" = "SC",
  "19" = "Immune/Endothelial",
  "20" = "Outliers",
  "21" = "Outliers",
  "22" = "SC",
  "23" = "ES",
  "24" = "ES",
  "25" = "Leydig"
)

cluster_annotation_fine <- list(
  "1" = "Round Spermatids (2)",
  "2" = "Elongating Spermatids (5)",
  "3" = "Spermatogonia",
  "4" = "Sertoli",
  "5" = "Round Spermatids (3)",
  "6" = "Round Spermatids (5)",
  "7" = "Outliers",
  "8" = "Sertoli",
  "9" = "Elongating Spermatids (2)",
  "10" = "Elongating Spermatids (3)",
  "11" = "Diplotene",
  "12" = "Elongating Spermatids (4)",
  "13" = "Outliers", 
  "14" = "Elongating Spermatids (1)",
  "15" = "Round Spermatids (4)",
  "16" = "Outliers", 
  "17" = "Round Spermatids (1)",
  "18" = "Spermatogonia / Leptotene",
  "19" = "Immune/Endothelial",
  "20" = "Outliers",
  "21" = "Outliers",
  "22" = "Meiosis",
  "23" = "Elongating Spermatids (6)",
  "24" = "Elongating Spermatids (7)",
  "25" = "Leydig"
)

sce.all$AnnotatedClusters <- unlist(cluster_annotation[as.character(clusters)])
sce.all$AnnotatedClustersFine <- unlist(cluster_annotation_fine[as.character(clusters)])

ggplot(data.frame(tsne1 = reducedDims(sce.all)$UMAP[,1],
                  tsne2 = reducedDims(sce.all)$UMAP[,2],
                  cluster = clusters)) +
  geom_point(aes(tsne1, tsne2, colour = cluster))

ggplot(data.frame(tsne1 = reducedDims(sce.all)$UMAP[,1],
                  tsne2 = reducedDims(sce.all)$UMAP[,2],
                  cluster = sce.all$AnnotatedClusters)) +
  geom_point(aes(tsne1, tsne2, colour = cluster))

pp <- ggplot(data.frame(tsne1 = reducedDims(sce.all)$UMAP[,1],
                  tsne2 = reducedDims(sce.all)$UMAP[,2],
                  cluster = unlist(cluster_annotation_fine[clusters]))) +
  geom_point(aes(tsne1, tsne2, colour = cluster))

plotly::ggplotly(pp)

##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
sce.all.no.outliers <- sce.all[,sce.all$AnnotatedClusters != "Outliers" & !sce.all$AnnotatedClusters %in% c("Leydig", "Immune/Endothelial", "Sertoli")]

simple_sfs <- librarySizeFactors(sce.all.no.outliers)

sce.all.no.outliers <- computeSumFactors(sce.all.no.outliers)
sce.all.no.outliers <- logNormCounts(sce.all.no.outliers, size_factors = simple_sfs)

set.seed(123)
HVGs_here <- HVG(sce.all.no.outliers)
pca_here <- prcomp(t(logcounts(sce.all.no.outliers[HVGs_here, ])))
prcurve <- principal_curve(pca_here$x[,1:2])

pseudotime = prcurve$lambda
sce.all.no.outliers$Pseudotime <- pseudotime

# add pseudotime to main object
sce.all$Pseudotime <-  NA
sce.all[,colnames(sce.all.no.outliers)]$Pseudotime <- sce.all.no.outliers$Pseudotime

ggplot(data.frame(
  umap1 = reducedDims(sce.all.no.outliers)$UMAP[,1], 
  umap2 = reducedDims(sce.all.no.outliers)$UMAP[,2], 
  Library = sce.all.no.outliers$Library,
  PT = sce.all.no.outliers$Pseudotime,
  CellType = sce.all.no.outliers$AnnotatedClusters
), aes(umap1, umap2, col = PT, shape = CellType)) + geom_point()

ggplot(data.frame(
  umap1 = reducedDims(sce.all)$UMAP[,1], 
  umap2 = reducedDims(sce.all)$UMAP[,2], 
  Library = sce.all$Library,
  PT = sce.all$Pseudotime,
  CellType = sce.all$AnnotatedClusters
), aes(umap1, umap2, col = PT)) + geom_point()

# temporary
saveRDS(sce.all, "./Data/processed/final_sce_f1_dataset_for_supplement.rds")
saveRDS(sce.all.no.outliers, "./Data/processed/final_sce_f1_dataset_for_supplement_pseudotime.rds")

# final
saveRDS(sce.all, "./Data/processed/final_sce_f1_dataset.rds")
saveRDS(sce.all.no.outliers, "./Data/processed/final_sce_f1_dataset_pseudotime.rds")

