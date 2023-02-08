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

colnames(sce.merged) <- colnames(sce.merged)
cells.both <- intersect(colnames(sce.filtered), colnames(sce.merged))
sce.merged <- sce.merged[,cells.both]
colData(sce.merged) <- cbind(colData(sce.merged), colData(sce.filtered[,cells.both]))
reducedDims(sce.merged)$TSNE <- reducedDims(sce.filtered[,cells.both])$TSNE
reducedDims(sce.merged)$UMAP <- reducedDims(sce.filtered[,cells.both])$UMAP
metadata(sce.merged)$corrected <- metadata(sce.filtered)$corrected
sce.merged <- logNormCounts(sce.merged)

saveRDS(sce.merged, "./Data/processed/sce_merged_new.rds")

