library(ggplot2)
library(scran)
library(scater)

setwd("~/Desktop/PhD/Projects/ASE_Spermatogenesis_Paper_rerun/")

source("./Scripts/General/auxiliary.R")
source("./Scripts/General/ase_functions.R")
source("./Scripts/General/reuse_functions.R")

# we get the "true" cells from the standard pipeline
sce.all <- readRDS("./Data/processed/final_sce_f1_dataset.rds")
corrected.expression <- metadata(sce.all)$corrected
colnames(sce.all) <- paste0(sce.all$Library, "_", sce.all$Barcode)
rownames(corrected.expression) <- colnames(sce.all)
sce.filtered <- sce.all[,!grepl("Outliers", sce.all$AnnotatedClusters)]
rm(sce.all)

### read ase-data
data.sample1 <- read_10x_ase("./data/raw/f1_dataset/ase_feature_matrix_sample1/")
data.sample2 <- read_10x_ase("./data/raw/f1_dataset/ase_feature_matrix_sample2/")
data.sample3 <- read_10x_ase("./data/raw/f1_dataset/ase_feature_matrix_sample3/")
data.sample4 <- read_10x_ase("./data/raw/f1_dataset/ase_feature_matrix_sample4/")
data.sample5 <- read_10x_ase("./data/raw/f1_dataset/ase_feature_matrix_sample5/")
data.sample6 <- read_10x_ase("./data/raw/f1_dataset/ase_feature_matrix_sample6/")

add_empty_rows <- function(data, vector){
  data = data[rownames(data) %in% vector,]
  new_rows = setdiff(vector, rownames(data))
  add_zeros = matrix(rep(0, length(new_rows) * ncol(data)), 
                     nrow = length(new_rows),
                     ncol = ncol(data))
  rownames(add_zeros) = new_rows
  dd = rbind(data, add_zeros)
  #rownames(dd) = c(rownames(data), new_rows)
  dd[vector, ]
  # dd_out <- do.call("rbind", lapply(vector, function(x){
  #   print(x)
  #   dd[x,]
  # }))
  #dd
}

make_ase_sce <- function(sce_full, data_reference, data_alternative){
  colnames(sce_full) <- sce_full$Barcode
  measured.genes = union(rownames(data_reference), 
                         rownames(data_alternative))
  measured.genes = union(rownames(sce_full), measured.genes)
  data_full = add_empty_rows(as.matrix(counts(sce_full)), measured.genes)
  data_reference = add_empty_rows(data_reference, measured.genes)
  data_alternative = add_empty_rows(data_alternative, measured.genes)
  
  #return(list(data_full, data_reference, data_alternative))
  
  barcodes_use <- sce_full$Barcode
  barcodes_use <- intersect(barcodes_use, 
                            intersect(colnames(data_reference), 
                                      colnames(data_alternative)))
  
  sce <- SingleCellExperiment(
    assays = list("counts" = as(data_full[,barcodes_use], "dgTMatrix"),
                  "counts_reference" = as(data_reference[,barcodes_use], "dgTMatrix"), 
                  "counts_alternative" = as(data_alternative[,barcodes_use], "dgTMatrix"))
  )
  sce
}

# Sample1
sce.sample1 <- make_ase_sce(sce.filtered[ ,sce.filtered$Library == "Sample1"], 
                            data.sample1$reference,
                            data.sample1$alternative)
sample1.celltypes <- sce.filtered[,sce.filtered$Library == "Sample1"]$AnnotatedClusters
names(sample1.celltypes) <- sce.filtered[,sce.filtered$Library == "Sample1"]$Barcode
sce.sample1$CellType <- sample1.celltypes[colnames(sce.sample1)]
sce.sample1$Library <- "Sample1"
sce.sample1$Species <- "B6"
rm(data.sample1)

# Sample2
sce.sample2 <- make_ase_sce(sce.filtered[ ,sce.filtered$Library == "Sample2"], 
                            data.sample2$reference,
                            data.sample2$alternative)
sample2.celltypes <- sce.filtered[,sce.filtered$Library == "Sample2"]$AnnotatedClusters
names(sample2.celltypes) <- sce.filtered[,sce.filtered$Library == "Sample2"]$Barcode
sce.sample2$CellType <- sample2.celltypes[colnames(sce.sample2)]
sce.sample2$Library <- "Sample2"
sce.sample2$Species <- "B6"
rm(data.sample2)

# Sample3
sce.sample3 <- make_ase_sce(sce.filtered[ ,sce.filtered$Library == "Sample3"], 
                            data.sample3$reference,
                            data.sample3$alternative)
sample3.celltypes <- sce.filtered[,sce.filtered$Library == "Sample3"]$AnnotatedClusters
names(sample3.celltypes) <- sce.filtered[,sce.filtered$Library == "Sample3"]$Barcode
sce.sample3$CellType <- sample3.celltypes[colnames(sce.sample3)]
sce.sample3$Library <- "Sample3"
sce.sample3$Species <- "CAST"
rm(data.sample3)

# Sample4
sce.sample4 <- make_ase_sce(sce.filtered[ ,sce.filtered$Library == "Sample4"], 
                            data.sample4$reference,
                            data.sample4$alternative)
sample4.celltypes <- sce.filtered[,sce.filtered$Library == "Sample4"]$AnnotatedClusters
names(sample4.celltypes) <- sce.filtered[,sce.filtered$Library == "Sample4"]$Barcode
sce.sample4$CellType <- sample4.celltypes[colnames(sce.sample4)]
sce.sample4$Library <- "Sample4"
sce.sample4$Species <- "CAST"
rm(data.sample4)

# Sample5
sce.sample5 <- make_ase_sce(sce.filtered[ ,sce.filtered$Library == "Sample5"], 
                            data.sample5$reference,
                            data.sample5$alternative)
sample5.celltypes <- sce.filtered[,sce.filtered$Library == "Sample5"]$AnnotatedClusters
names(sample5.celltypes) <- sce.filtered[,sce.filtered$Library == "Sample5"]$Barcode
sce.sample5$CellType <- sample5.celltypes[colnames(sce.sample5)]
sce.sample5$Library <- "Sample5"
sce.sample5$Species <- "B6xCAST"
rm(data.sample5)

# Sample6
sce.sample6 <- make_ase_sce(sce.filtered[ ,sce.filtered$Library == "Sample6"], 
                            data.sample6$reference,
                            data.sample6$alternative)
sample6.celltypes <- sce.filtered[,sce.filtered$Library == "Sample6"]$AnnotatedClusters
names(sample6.celltypes) <- sce.filtered[,sce.filtered$Library == "Sample6"]$Barcode
sce.sample6$CellType <- sample6.celltypes[colnames(sce.sample6)]
sce.sample6$Library <- "Sample6"
sce.sample6$Species <- "B6xCAST"
rm(data.sample6)

joint.genes <- Reduce(intersect, 
                      list(rownames(sce.sample1),
                           rownames(sce.sample2),
                           rownames(sce.sample3),
                           rownames(sce.sample4),
                           rownames(sce.sample5),
                           rownames(sce.sample6)))

sce.merged <- Reduce(cbind, 
                     list(sce.sample1[joint.genes,],
                          sce.sample2[joint.genes,],
                          sce.sample3[joint.genes,],
                          sce.sample4[joint.genes,],
                          sce.sample5[joint.genes,],
                          sce.sample6[joint.genes,]))

rowData(sce.merged) <- rowData(sce.filtered)[rownames(sce.merged), ]
sce.merged_2 <- annotate_chromosome_sce(sce.merged)
colnames(sce.merged) <- paste0(sce.merged$Library, "_", colnames(sce.merged))

rm(list = c("sce.sample1", "sce.sample2", 
            "sce.sample3", "sce.sample4",
            "sce.sample5", "sce.sample6"))

# # Do some QC
# total_counts_reference <- lapply(
#   list(assays(sce.sample1)[['counts_reference']],
#        assays(sce.sample2)[['counts_reference']],
#        assays(sce.sample3)[['counts_reference']],
#        assays(sce.sample4)[['counts_reference']],
#        assays(sce.sample5)[['counts_reference']],
#        assays(sce.sample6)[['counts_reference']]),
#   colSums
# )
# names(total_counts_reference) <- paste0("Sample_", 1:6)
# 
# total_counts_alternative <- lapply(
#   list(assays(sce.sample1)[['counts_alternative']],
#        assays(sce.sample2)[['counts_alternative']],
#        assays(sce.sample3)[['counts_alternative']],
#        assays(sce.sample4)[['counts_alternative']],
#        assays(sce.sample5)[['counts_alternative']],
#        assays(sce.sample6)[['counts_alternative']]),
#   colSums
# )
# 
# 
# 
# df_total_counts_reference <- do.call("rbind", lapply(1:6, function(i){
#   x = total_counts_reference[[i]]
#   data.frame(
#     Sample = paste0("Sample_", i),
#     Allele = "Reference",
#     Counts = x
#   )
# }))
# 
# df_total_counts_alternative <- do.call("rbind", lapply(1:6, function(i){
#   x = total_counts_alternative[[i]]
#   data.frame(
#     Sample = paste0("Sample_", i),
#     Allele = "Alternative",
#     Counts = x
#   )
# }))
# 
# df_total_counts <- rbind(df_total_counts_reference, 
#                          df_total_counts_alternative)
# 
# df_total_counts$Sample <- unlist(list("Sample_1" = "B6parental_Sample1", 
#                                       "Sample_2" = "B6parental_Sample2",
#                                       "Sample_3" = "Castparental_Sample3",
#                                       "Sample_4" = "Castparental_Sample4",
#                                       "Sample_5" = "B6xCast_Sample5",
#                                       "Sample_6" = "B6xCast_Sample6"))[df_total_counts$Sample]
# 
# df_total_counts$Sample <- factor(df_total_counts$Sample, levels = c("B6parental_Sample1", 
#                                                                     "B6parental_Sample2",
#                                                                     "Castparental_Sample3",
#                                                                     "Castparental_Sample4",
#                                                                     "B6xCast_Sample5",
#                                                                     "B6xCast_Sample6"))
# 
# # Plot total counts per cell and highlight fold change over it:
# average.reference <- aggregate(df_total_counts[df_total_counts$Allele == "Reference", 3], 
#                                list(df_total_counts[df_total_counts$Allele == "Reference",]$Sample), mean)$x
# average.alternative <- aggregate(df_total_counts[df_total_counts$Allele == "Alternative", 3], 
#                                  list(df_total_counts[df_total_counts$Allele == "Alternative",]$Sample), mean)$x
# ratios <- round(average.reference / (average.reference + average.alternative), digits = 3) * 100
# 
# p_top <- ggplot(df_total_counts, aes(x = Sample, y = Counts + 1, fill = Allele)) + 
#   geom_boxplot() + 
#   ylab("log(ReadCounts + 1) per cell") + 
#   xlab("") + 
#   scale_y_log10(limits = c(1, 1.5e5)) + 
#   annotate("text", label = paste0(as.character(ratios), "%"), 
#            x = 1:6, y = 60000) + 
#   annotate("text", label = "Percent reference reads: ", 
#            x = 1.6, y = 140000) + 
#   theme_classic() + 
#   theme(axis.text.x = element_text(angle = 45, hjust=1))
# 
# # Plot aggregated data across cells to compare on gene level
# 
# sce_sample1 <- sce.merged[,sce.merged$Library == "Sample1"]
# df_plot_1 <- data.frame(
#   Counts_Reference = rowMeans(counts_reference(sce_sample1)),
#   Counts_Alternative = 0
# )
# df_plot_1[names(rowMeans(counts_reference(sce_sample1))),]$Counts_Alternative <- 
#   rowMeans(counts_alternative(sce_sample1))
# 
# sce_sample3 <- sce.merged[,sce.merged$Library == "Sample3"]
# df_plot_3 <- data.frame(
#   Counts_Reference = rowMeans(counts_reference(sce_sample3)),
#   Counts_Alternative = 0
# )
# df_plot_3[names(rowMeans(counts_reference(sce_sample3))),]$Counts_Alternative <- 
#   rowMeans(counts_alternative(sce_sample3))
# 
# sce_sample5 <- sce.merged[,sce.merged$Library == "Sample5"]
# sce_sample5 <- annotate_chromosome(sce_sample5)
# df_plot_5 <- data.frame(
#   Counts_Reference = rowMeans(counts_reference(sce_sample5)),
#   Counts_Alternative = 0,
#   Chromosome = rowData(sce_sample5)$chromosome_name %in% c("MT", "X")
# )
# df_plot_5[names(rowMeans(counts_reference(sce_sample5))),]$Counts_Alternative <- 
#   rowMeans(counts_alternative(sce_sample5))
# 
# p1 <- ggplot(df_plot_1, aes(log(Counts_Reference + 1), log(Counts_Alternative + 1))) + 
#   geom_point(size = 0.5) + 
#   xlim(0, 4) + ylim(0, 4) + 
#   coord_fixed() + 
#   theme_classic() + 
#   ggtitle("Reads mapping to reference/alternative \n alleles in B6")
# 
# p2 <- ggplot(df_plot_3, aes(log(Counts_Reference + 1), log(Counts_Alternative + 1))) + 
#   geom_point(size = 0.5) + 
#   xlim(0, 4) + ylim(0, 4) + 
#   coord_fixed() + 
#   theme_classic() + 
#   ggtitle("Reads mapping to reference/alternative \n alleles in CAST")
# 
# p3 <- ggplot(df_plot_5, aes(log(Counts_Reference + 1), log(Counts_Alternative + 1), 
#                             col = Chromosome)) + 
#   geom_point(size = 0.5) + 
#   xlim(0, 4) + ylim(0, 4) + 
#   coord_fixed() + 
#   theme_classic() + 
#   ggtitle("Reads mapping to reference/alternative \n alleles in F1") + 
#   scale_color_manual(values = list("FALSE" = "black", "TRUE" = "red"),
#                      name = "X/MT-chromosome")
# 
# rm(sce_sample1)
# rm(sce_sample3)
# rm(sce_sample5)

colnames(sce.merged) <- colnames(sce.merged)
cells.both <- intersect(colnames(sce.filtered), colnames(sce.merged))
sce.merged <- sce.merged[,cells.both]
colData(sce.merged) <- cbind(colData(sce.merged), colData(sce.filtered[,cells.both]))
reducedDims(sce.merged)$TSNE <- reducedDims(sce.filtered[,cells.both])$TSNE
reducedDims(sce.merged)$UMAP <- reducedDims(sce.filtered[,cells.both])$UMAP
metadata(sce.merged)$corrected <- metadata(sce.filtered)$corrected
sce.merged <- logNormCounts(sce.merged)

saveRDS(sce.merged, "./Data/processed/sce_merged_new.rds")

# Do some QC
# total_counts_reference <- lapply(
#   unique(sce.merged$Library),
#   function(x){
#     sce.filtered <- assays(sce.merged[,sce.merged$Library == x])[["counts_reference"]]
#     colSums(sce.filtered)
#   }
# )
# names(total_counts_reference) <- unique(sce.merged$Library)
# 
# total_counts_alternative <- lapply(
#   unique(sce.merged$Library),
#   function(x){
#     sce.filtered <- assays(sce.merged[,sce.merged$Library == x])[["counts_alternative"]]
#     colSums(sce.filtered)
#   }
# )
# names(total_counts_alternative) <- unique(sce.merged$Library)
# 
# df_total_counts_reference <- do.call("rbind", lapply(1:6, function(i){
#   x = total_counts_reference[[i]]
#   data.frame(
#     Sample = paste0("Sample_", i),
#     Allele = "Reference",
#     Counts = x
#   )
# }))
# 
# df_total_counts_alternative <- do.call("rbind", lapply(1:6, function(i){
#   x = total_counts_alternative[[i]]
#   data.frame(
#     Sample = paste0("Sample_", i),
#     Allele = "Alternative",
#     Counts = x
#   )
# }))
# 
# df_total_counts <- rbind(df_total_counts_reference,
#                          df_total_counts_alternative)
# 
# df_total_counts$Sample <- unlist(list("Sample_1" = "B6parental_Sample1",
#                                       "Sample_2" = "B6parental_Sample2",
#                                       "Sample_3" = "Castparental_Sample3",
#                                       "Sample_4" = "Castparental_Sample4",
#                                       "Sample_5" = "B6xCast_Sample5",
#                                       "Sample_6" = "B6xCast_Sample6"))[df_total_counts$Sample]
# 
# df_total_counts$Sample <- factor(df_total_counts$Sample, levels = c("B6parental_Sample1",
#                                                                     "B6parental_Sample2",
#                                                                     "Castparental_Sample3",
#                                                                     "Castparental_Sample4",
#                                                                     "B6xCast_Sample5",
#                                                                     "B6xCast_Sample6"))
# 
# # Plot total counts per cell and highlight fold change over it:
# average.reference <- aggregate(df_total_counts[df_total_counts$Allele == "Reference", 3],
#                                list(df_total_counts[df_total_counts$Allele == "Reference",]$Sample), mean)$x
# average.alternative <- aggregate(df_total_counts[df_total_counts$Allele == "Alternative", 3],
#                                  list(df_total_counts[df_total_counts$Allele == "Alternative",]$Sample), mean)$x
# ratios <- round(average.reference / (average.reference + average.alternative), digits = 3) * 100
# 
# p_top <- ggplot(df_total_counts, aes(x = Sample, y = Counts + 1, fill = Allele)) +
#   geom_boxplot() +
#   ylab("log10(ReadCounts + 1) per cell") +
#   xlab("") +
#   scale_y_log10(limits = c(1, 1.5e5)) +
#   annotate("text", label = paste0(as.character(ratios), "%"),
#            x = 1:6, y = 60000, size = 5) +
#   annotate("text", label = "Percent reference reads: ",
#            x = 1.6, y = 140000, size = 5) +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 45, hjust=1),
#         text = element_text(size = 20))
# 
# p_top
# 
# # Plot aggregated data across cells to compare on gene level
# 
# sce_sample1 <- sce.merged[,sce.merged$Library == "Sample1"]
# df_plot_1 <- data.frame(
#   Counts_Reference = rowMeans(counts_reference(sce_sample1)),
#   Counts_Alternative = 0
# )
# df_plot_1[names(rowMeans(counts_reference(sce_sample1))),]$Counts_Alternative <-
#   rowMeans(counts_alternative(sce_sample1))
# 
# sce_sample3 <- sce.merged[,sce.merged$Library == "Sample3"]
# df_plot_3 <- data.frame(
#   Counts_Reference = rowMeans(counts_reference(sce_sample3)),
#   Counts_Alternative = 0
# )
# df_plot_3[names(rowMeans(counts_reference(sce_sample3))),]$Counts_Alternative <-
#   rowMeans(counts_alternative(sce_sample3))
# 
# sce_sample5 <- sce.merged[,sce.merged$Library == "Sample5"]
# sce_sample5 <- annotate_chromosome(sce_sample5)
# df_plot_5 <- data.frame(
#   Counts_Reference = rowMeans(counts_reference(sce_sample5)),
#   Counts_Alternative = 0,
#   Chromosome = rowData(sce_sample5)$chromosome_name %in% c("MT", "X")
# )
# df_plot_5[names(rowMeans(counts_reference(sce_sample5))),]$Counts_Alternative <-
#   rowMeans(counts_alternative(sce_sample5))
# 
# 
# p1 <- ggplot(df_plot_1, aes(log(Counts_Reference + 1), log(Counts_Alternative + 1))) +
#   geom_point(size = 0.5) +
#   xlim(0, 4) + ylim(0, 4) +
#   coord_fixed() +
#   theme_classic() +
#   geom_abline() +
#   ggtitle("Reads mapping to reference/alternative \n alleles in B6") +
#   theme(text = element_text(size = 20))
# 
# 
# p2 <- ggplot(df_plot_3, aes(log(Counts_Reference + 1), log(Counts_Alternative + 1))) +
#   geom_point(size = 0.5) +
#   xlim(0, 4) + ylim(0, 4) +
#   coord_fixed() +
#   theme_classic() +
#   geom_abline() +
#   ggtitle("Reads mapping to reference/alternative \n alleles in CAST") +
#   theme(text = element_text(size = 20))
# 
# 
# p3 <- ggplot(df_plot_5, aes(log(Counts_Reference + 1), log(Counts_Alternative + 1),
#                             col = Chromosome)) +
#   geom_point(size = 0.5) +
#   xlim(0, 4) + ylim(0, 4) +
#   coord_fixed() +
#   theme_classic() +
#   geom_abline() +
#   ggtitle("Reads mapping to reference/alternative \n alleles in F1") +
#   scale_color_manual(values = list("FALSE" = "black", "TRUE" = "red"),
#                      name = "X/MT-chromosome") +
#   theme(text = element_text(size = 20))
# 
# p1
# p2
# p3
# 
# # Look at ratios per cell type
# 
# df_total_counts$CellType <- rep(sce.merged$CellType, 2)
# df_total_counts_f1 <- df_total_counts[df_total_counts$Sample %in% c("B6xCast_Sample5", "B6xCast_Sample6"), ]
# 
# ggplot(df_total_counts_f1, aes(x = CellType, y = Counts + 1, fill = Allele)) +
#   geom_boxplot() +
#   ylab("log10(ReadCounts + 1) per cell") +
#   xlab("") +
#   scale_y_log10() +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 45, hjust=1),
#         text = element_text(size = 20))
# 
# # How many allele-specific reads per gene in F1?
# 
# total_counts_per_gene <- rowSums(counts_reference(sce.merged[,sce.merged$Library %in% c("Sample5", "Sample6")])) +
#   rowSums(counts_alternative(sce.merged[,sce.merged$Library %in% c("Sample5", "Sample6")]))
# 
# total_counts_per_gene_df <- data.frame(total_counts_per_gene)
# total_counts_per_gene_df$include <- total_counts_per_gene_df$total_counts_per_gene > 1000
# total_counts_per_gene_df[total_counts_per_gene_df$total_counts_per_gene == 0, ]$include <- "no_quantification"
# 
# n_included_genes <- dim(total_counts_per_gene_df[total_counts_per_gene_df$include == "TRUE", ])[[1]]
# percent_included_genes <- n_included_genes / dim(total_counts_per_gene_df)[[1]]
# 
# ggplot(total_counts_per_gene_df, aes(x = total_counts_per_gene + 1, fill = include)) + geom_histogram() + scale_x_log10() +
#   theme_classic() +
#   annotate("text", label = paste0("Included genes: ", n_included_genes) , x = 1000, y = 15000, size = 5) +
#   annotate("text", label = paste0("Percent Included genes: ", round(percent_included_genes, digits = 4)) , x = 1000, y = 14000, size = 5)

