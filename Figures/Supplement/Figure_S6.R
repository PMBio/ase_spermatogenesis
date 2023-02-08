library(scran)
library(scater)
library(tidyverse)

setwd("~/Desktop/Projects/ASE_Spermatogenesis_Paper_rerun/Revisions/")

source("../Scripts/General/auxiliary.R")
source("../Scripts/General/reuse_functions.R")

colours_cis_trans <- c(
  "conserved" =  "grey",
  "cis" =  "orange",
  "trans" =  "cyan",
  "cis_trans" =  "Purple"
)

data <- readRDS("./Data/processed/sce_merged_new.rds")
data$IndividualSamples <- paste0(data$Dataset, "_", data$Library)
data <- data[,data$Library != "Sample7"]

aggregate_test <- aggregateAcrossCells(data, data$IndividualSamples, use.assay.type = c("counts_reference", "counts_alternative"))

normalize_sf <- function(counts){
  col.save <- colnames(counts)
  row.save <- rownames(counts)
  sfs <- colSums(counts) / colSums(counts)[[1]]
  counts <- do.call("cbind", lapply(1:ncol(counts), function(i){
    return(counts[,i] / sfs[[i]])
  }))
  colnames(counts) <- col.save
  rownames(counts) <- row.save
  return(round(counts))
}
make_bulk_celltype <- function(data){
  
  aggregate_sce <- aggregateAcrossCells(data, data$IndividualSamples, use.assay.type = c("counts_reference", "counts_alternative"))
  
  bulk_par_b6 <- counts_reference(aggregate_sce[,aggregate_sce$Species == "B6"])
  colnames(bulk_par_b6) <- paste0(colnames(bulk_par_b6), "_Reference")
  bulk_par_cast <- counts_alternative(aggregate_sce[,aggregate_sce$Species == "CAST"])
  colnames(bulk_par_cast) <- paste0(colnames(bulk_par_cast), "_Alternative")
  
  sce_par_cur_perLibrary <- cbind(bulk_par_b6, bulk_par_cast)
  
  bulk_fil_b6 <- counts_reference(aggregate_sce[,aggregate_sce$Species == "B6xCAST"])
  colnames(bulk_fil_b6) <- paste0(colnames(bulk_fil_b6), "_Reference")
  bulk_fil_cast <- counts_alternative(aggregate_sce[,aggregate_sce$Species == "B6xCAST"])
  colnames(bulk_fil_cast) <- paste0(colnames(bulk_fil_cast), "_Alternative")
  
  sce_fil_cur_perLibrary <- cbind(bulk_fil_b6, bulk_fil_cast)
  
  # sce_par_cur_perLibrary = data.frame(
  #   "Sample1_Reference" = rowSums(counts(data[,data$Library == "Sample1"])),
  #   "Sample2_Reference" = rowSums(counts(data[,data$Library == "Sample2"])),
  #   "Sample3_Alternative" = rowSums(counts(data[,data$Library == "Sample3"])),
  #   "Sample4_Alternative" = rowSums(counts(data[,data$Library == "Sample4"]))
  # )
  # 
  # sce_fil_cur_perLibrary = data.frame(
  #   "Sample5_Reference" = rowSums(counts_reference(data[,data$Library == "Sample5"])),
  #   "Sample6_Reference" = rowSums(counts_reference(data[,data$Library == "Sample6"])),
  #   "Sample5_Alternative" = rowSums(counts_alternative(data[,data$Library == "Sample5"])),
  #   "Sample6_Alternative" = rowSums(counts_alternative(data[,data$Library == "Sample6"]))
  # )
  # 
  genes.both <- intersect(rownames(sce_par_cur_perLibrary),
                          rownames(sce_fil_cur_perLibrary))
  
  sce_par_cur_perLibrary <- sce_par_cur_perLibrary[genes.both,]
  sce_fil_cur_perLibrary <- sce_fil_cur_perLibrary[genes.both,]
  
  # disps.par.b6 <- estimateDisp(sce_par_cur_perLibrary[,1:2])$tagwise.dispersion
  # names(disps.par.b6) <- rownames(sce_par_cur_perLibrary)
  # disps.par.cast <- estimateDisp(sce_par_cur_perLibrary[,3:4])$tagwise.dispersion
  # names(disps.par.cast) <- rownames(sce_par_cur_perLibrary)
  
  sce_par_cur_perLibrary_norm <- normalize_sf(sce_par_cur_perLibrary)
  sce_fil_cur_perLibrary_norm <- normalize_sf(sce_fil_cur_perLibrary)
  
  return(list(sce_par_cur_perLibrary_norm, sce_fil_cur_perLibrary_norm))
}

bulk_sce <- aggregateAcrossCells(data, data$IndividualSamples, use.assay.type = c("counts", "counts_reference", "counts_alternative"))
bulk_sce$sizeFactor <- calculateSumFactors(bulk_sce)
bulk_sce <- logNormCounts(bulk_sce)
reducedDims(bulk_sce)[["PCA"]] <- calculatePCA(bulk_sce)
plotPCA(bulk_sce, colour_by = "Dataset", text_by = "Species")

bulk_all <- make_bulk_celltype(data)
bulk_sc <- make_bulk_celltype(data[,data$CellType == "SC"])
bulk_rs <- make_bulk_celltype(data[,data$CellType == "RS"])
bulk_es <- make_bulk_celltype(data[,data$CellType == "ES"])





model.selection.results.bulk <- readRDS("./Data/processed/model_selection_bulk.rds")
model.selection.results.sc <- readRDS("./Data/processed/model_selection_sc.rds")
model.selection.results.rs <- readRDS("./Data/processed/model_selection_rs.rds")
model.selection.results.es <- readRDS("./Data/processed/model_selection_es.rds")

make_fc_strat_plot <- function(data){
  data %>%
    mutate(fold_change_interval = cut(sqrt(FC_filial ** 2 + FC_parental ** 2), c(0, 0.05, 0.1, 0.2, 0.5, 1, 1.5, 2, Inf))) %>%
    group_by(fold_change_interval, CategoryNew) %>%
    summarize(n_genes = n()) %>%
    mutate(n_total = sum(n_genes)) %>%
    dplyr::filter(!is.na(fold_change_interval)) %>%
    ggplot(aes(x = fold_change_interval, y = n_genes, fill = CategoryNew)) + geom_bar(position = "fill", stat = "identity") + 
    xlab("log aFC interval") + ylab("Proportion of effects") + scale_fill_manual(values = colours_cis_trans) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + theme_paper(textsize = 50) + 
    geom_text(aes(label = n_total, y = 1), vjust = 1.5, size = 10)
}
make_explevel_strat_plot <- function(data){
  data %>%
    mutate(fold_change_interval = cut(Expr_Par_B6, breaks = c(0, 10, 100, 1000, 10000, 100000, Inf))) %>%
    group_by(fold_change_interval, CategoryNew) %>%
    summarize(n_genes = n()) %>%
    mutate(n_total = sum(n_genes)) %>%
    dplyr::filter(!is.na(fold_change_interval)) %>%
    ggplot(aes(x = fold_change_interval, y = n_genes, fill = CategoryNew)) + geom_bar(position = "fill", stat = "identity") + 
    xlab("Gene expression interval") + ylab("Proportion of effects") + scale_fill_manual(values = colours_cis_trans) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + theme_paper(textsize = 50) + 
    geom_text(aes(label = n_total, y = 1), vjust = 1.5, size = 10)
}
make_effectsize_cis_plot <- function(data){
  data %>%
    ggplot(aes(x = CategoryNew, y = abs((FC_filial)))) + geom_jitter() + geom_boxplot(aes(fill = CategoryNew), outlier.colour = NA) + 
    theme_paper(textsize = 50) + coord_flip() + 
    ylab("magnitude cis-effect \n (|aFC (F1)|)") + xlab("") + scale_fill_manual(values = colours_cis_trans) + 
    geom_hline(yintercept = 0.5, linetype = "dashed") + 
    geom_hline(yintercept = 1, linetype = "dashed") + 
    theme(legend.position = "None")
}
make_effectsize_trans_plot <- function(data){
  data %>%
    ggplot(aes(x = CategoryNew, y = abs((FC_filial) - (FC_parental)))) + geom_jitter() + geom_boxplot(aes(fill = CategoryNew), outlier.colour = NA) +
    theme_paper(textsize = 50) + coord_flip() + 
    ylab("magnitude trans-effect \n (|aFC (F1) - aFC (F0)|)") + xlab("") + scale_fill_manual(values = colours_cis_trans) + 
    geom_hline(yintercept = 0.5, linetype = "dashed") + 
    geom_hline(yintercept = 1, linetype = "dashed") + 
    theme(legend.position = "None")
}

## get number of genes with high FCs in cis and trans:
table(abs(model.selection.results.bulk[model.selection.results.bulk$CategoryNew %in% c("cis", "cis_trans"), ]$FC_filial) > 1) / 
  sum(table(abs(model.selection.results.bulk[model.selection.results.bulk$CategoryNew %in% c("cis", "cis_trans"), ]$FC_filial) > 1))
table(abs(model.selection.results.sc[model.selection.results.sc$CategoryNew %in% c("cis", "cis_trans"), ]$FC_filial) > 0.5) / 
  sum(table(abs(model.selection.results.sc[model.selection.results.sc$CategoryNew %in% c("cis", "cis_trans"), ]$FC_filial) > 0.5))
table(abs(model.selection.results.sc[model.selection.results.sc$CategoryNew %in% c("cis", "cis_trans"), ]$FC_filial) > 1) / 
  sum(table(abs(model.selection.results.sc[model.selection.results.sc$CategoryNew %in% c("cis", "cis_trans"), ]$FC_filial) > 1))
table(abs(model.selection.results.rs[model.selection.results.rs$CategoryNew %in% c("cis", "cis_trans"), ]$FC_filial) > 1) / 
  sum(table(abs(model.selection.results.rs[model.selection.results.rs$CategoryNew %in% c("cis", "cis_trans"), ]$FC_filial) > 1))
table(abs(model.selection.results.es[model.selection.results.es$CategoryNew %in% c("cis", "cis_trans"), ]$FC_filial) > 1) / 
  sum(table(abs(model.selection.results.es[model.selection.results.es$CategoryNew %in% c("cis", "cis_trans"), ]$FC_filial) > 1))

table(abs(model.selection.results.bulk[model.selection.results.bulk$CategoryNew %in% c("trans", "cis_trans"), ]$FC_parental - 
            model.selection.results.bulk[model.selection.results.bulk$CategoryNew %in% c("trans", "cis_trans"), ]$FC_filial) > 1) / 
  sum(table(abs(model.selection.results.bulk[model.selection.results.bulk$CategoryNew %in% c("trans", "cis_trans"), ]$FC_parental - 
                  model.selection.results.bulk[model.selection.results.bulk$CategoryNew %in% c("trans", "cis_trans"), ]$FC_filial) > 1))
table(abs(model.selection.results.sc[model.selection.results.sc$CategoryNew %in% c("trans", "cis_trans"), ]$FC_parental - 
            model.selection.results.sc[model.selection.results.sc$CategoryNew %in% c("trans", "cis_trans"), ]$FC_filial) > 1) / 
  sum(table(abs(model.selection.results.sc[model.selection.results.sc$CategoryNew %in% c("trans", "cis_trans"), ]$FC_parental - 
                  model.selection.results.sc[model.selection.results.sc$CategoryNew %in% c("trans", "cis_trans"), ]$FC_filial) > 1))
table(abs(model.selection.results.rs[model.selection.results.rs$CategoryNew %in% c("trans", "cis_trans"), ]$FC_parental - 
            model.selection.results.rs[model.selection.results.rs$CategoryNew %in% c("trans", "cis_trans"), ]$FC_filial) > 1) / 
  sum(table(abs(model.selection.results.rs[model.selection.results.rs$CategoryNew %in% c("trans", "cis_trans"), ]$FC_parental - 
                  model.selection.results.rs[model.selection.results.rs$CategoryNew %in% c("trans", "cis_trans"), ]$FC_filial) > 1))
table(abs(model.selection.results.es[model.selection.results.es$CategoryNew %in% c("trans", "cis_trans"), ]$FC_parental - 
            model.selection.results.es[model.selection.results.es$CategoryNew %in% c("trans", "cis_trans"), ]$FC_filial) > 1) / 
  sum(table(abs(model.selection.results.es[model.selection.results.es$CategoryNew %in% c("trans", "cis_trans"), ]$FC_parental - 
                  model.selection.results.es[model.selection.results.es$CategoryNew %in% c("trans", "cis_trans"), ]$FC_filial) > 1))

table(abs(model.selection.results.sc[model.selection.results.sc$CategoryNew %in% c("cis", "cis_trans"), ]$FC_filial) > 1) / 
  sum(table(abs(model.selection.results.sc[model.selection.results.sc$CategoryNew %in% c("cis", "cis_trans"), ]$FC_filial) > 1))
table(abs(model.selection.results.rs[model.selection.results.rs$CategoryNew %in% c("cis", "cis_trans"), ]$FC_filial) > 1) / 
  sum(table(abs(model.selection.results.rs[model.selection.results.rs$CategoryNew %in% c("cis", "cis_trans"), ]$FC_filial) > 1))
table(abs(model.selection.results.es[model.selection.results.es$CategoryNew %in% c("cis", "cis_trans"), ]$FC_filial) > 1) / 
  sum(table(abs(model.selection.results.es[model.selection.results.es$CategoryNew %in% c("cis", "cis_trans"), ]$FC_filial) > 1))


make_fc_strat_plot(model.selection.results.bulk) + ggtitle("Bulk")
ggsave("~/Desktop/Projects/ASE_Spermatogenesis_Paper_rerun_Revisions/Plots/FigureS1_rev/fold_change_stratification_bulk.pdf")
make_fc_strat_plot(model.selection.results.sc) + ggtitle("SC")
ggsave("~/Desktop/Projects/ASE_Spermatogenesis_Paper_rerun_Revisions/Plots/FigureS1_rev/fold_change_stratification_sc.pdf")
make_fc_strat_plot(model.selection.results.rs) + ggtitle("RS")
ggsave("~/Desktop/Projects/ASE_Spermatogenesis_Paper_rerun_Revisions/Plots/FigureS1_rev/fold_change_stratification_rs.pdf")
make_fc_strat_plot(model.selection.results.es) + ggtitle("ES")
ggsave("~/Desktop/Projects/ASE_Spermatogenesis_Paper_rerun_Revisions/Plots/FigureS1_rev/fold_change_stratification_es.pdf")

make_explevel_strat_plot(model.selection.results.bulk) + ggtitle("Bulk")
ggsave("~/Desktop/Projects/ASE_Spermatogenesis_Paper_rerun_Revisions/Plots/FigureS1_rev/explevel_stratification_bulk.pdf")
make_explevel_strat_plot(model.selection.results.sc) + ggtitle("SC")
ggsave("~/Desktop/Projects/ASE_Spermatogenesis_Paper_rerun_Revisions/Plots/FigureS1_rev/explevel_stratification_sc.pdf")
make_explevel_strat_plot(model.selection.results.rs) + ggtitle("RS")
ggsave("~/Desktop/Projects/ASE_Spermatogenesis_Paper_rerun_Revisions/Plots/FigureS1_rev/explevel_stratification_rs.pdf")
make_explevel_strat_plot(model.selection.results.es) + ggtitle("ES")
ggsave("~/Desktop/Projects/ASE_Spermatogenesis_Paper_rerun_Revisions/Plots/FigureS1_rev/explevel_stratification_es.pdf")

make_effectsize_cis_plot(model.selection.results.bulk) + ggtitle("Bulk")
ggsave("~/Desktop/Projects/ASE_Spermatogenesis_Paper_rerun_Revisions/Plots/FigureS1_rev/effect_size_cis_bulk.pdf")
make_effectsize_cis_plot(model.selection.results.sc) + ggtitle("SC")
ggsave("~/Desktop/Projects/ASE_Spermatogenesis_Paper_rerun_Revisions/Plots/FigureS1_rev/effect_size_cis_sc.pdf")
make_effectsize_cis_plot(model.selection.results.rs) + ggtitle("RS")
ggsave("~/Desktop/Projects/ASE_Spermatogenesis_Paper_rerun_Revisions/Plots/FigureS1_rev/effect_size_cis_rs.pdf")
make_effectsize_cis_plot(model.selection.results.es) + ggtitle("ES")
ggsave("~/Desktop/Projects/ASE_Spermatogenesis_Paper_rerun_Revisions/Plots/FigureS1_rev/effect_size_cis_es.pdf")

# check how many genes show allelic fold changes above certain thresholds
table(abs(model.selection.results.sc[model.selection.results.sc$Category == "cis", ]$FC_filial) > 0.5)[[2]] / 
  sum(table(abs(model.selection.results.sc[model.selection.results.sc$Category == "cis", ]$FC_filial) > 0.5))

table(abs(model.selection.results.sc[model.selection.results.sc$Category == "cis", ]$FC_filial) > 1)[[2]] / 
  sum(table(abs(model.selection.results.sc[model.selection.results.sc$Category == "cis", ]$FC_filial) > 0.5))

make_effectsize_trans_plot(model.selection.results.bulk) + ggtitle("Bulk")
ggsave("~/Desktop/Projects/ASE_Spermatogenesis_Paper_rerun_Revisions/Plots/FigureS1_rev/effect_size_trans_bulk.pdf")
make_effectsize_trans_plot(model.selection.results.sc) + ggtitle("SC")
ggsave("~/Desktop/Projects/ASE_Spermatogenesis_Paper_rerun_Revisions/Plots/FigureS1_rev/effect_size_trans_sc.pdf")
make_effectsize_trans_plot(model.selection.results.rs) + ggtitle("RS")
ggsave("~/Desktop/Projects/ASE_Spermatogenesis_Paper_rerun_Revisions/Plots/FigureS1_rev/effect_size_trans_rs.pdf")
make_effectsize_trans_plot(model.selection.results.es) + ggtitle("ES")
ggsave("~/Desktop/Projects/ASE_Spermatogenesis_Paper_rerun_Revisions/Plots/FigureS1_rev/effect_size_trans_es.pdf")

##### 
# Idea: Replace categorization with BIC > X categorization
# Look at number and effect sizes per cell type using (beta-)-binomial model

library(DESeq2)

# check DESeq2 on parents results against binomial LM
deseq.test <- DESeqDataSetFromMatrix(bulk_all[[1]],
                                     colData = DataFrame(row.names = colnames(bulk_all[[1]]), Strain = c(rep("B6", 6), rep("CAST", 6))), 
                                     design = ~Strain)
deseq.test <- DESeq(deseq.test)
degs.test <- results(deseq.test)

degs.test %>% data.frame() %>%
  dplyr::filter(baseMean > 50) -> test
table(test$padj < 0.1)
perc <- paste0(round(table(test$padj < 0.1)[[2]] / sum(table(test$padj < 0.1)), digits = 4) * 100, "%")
perc_2 <- paste0(round(table(test$padj < 0.1 & abs(test$log2FoldChange > 0.5))[[2]] / sum(table(test$padj < 0.1)), digits = 4) * 100, "%")

make_volcano_plot <- function(data){
  data %>% data.frame() %>%
    dplyr::filter(baseMean > 50) %>%
    ggplot(aes(x = log2FoldChange, y = -log10(pvalue), col = padj < 0.1)) + geom_point() + theme_paper(textsize = 50) + 
    annotate(geom = "text", x = 0, y = -log10(min(degs.test$pvalue, na.rm = T)) + 45, size = 12, 
             label = paste0("DEGs (p_adj < 0.1): ", table(test$padj < 0.1)[[2]], " / ", sum(table(test$padj < 0.1)), 
                            " (", perc, ")")) + 
    annotate(geom = "text", x = 0, y = -log10(min(degs.test$pvalue, na.rm = T)) + 20, size = 12, 
             label = paste0("|LogFC| > 0.5: ", table(test$padj < 0.1 & abs(test$log2FoldChange > 0.5))[[2]], " / ", sum(table(test$padj < 0.1)), 
                            " (", perc_2, ")")) + 
    scale_color_manual(values = c("black", "red"))
}

make_volcano_plot(degs.test) + ggtitle("Bulk")
ggsave("~/Desktop/Projects/ASE_Spermatogenesis_Paper_rerun_Revisions/Plots/FigureS1_rev/volcano_plot_all.pdf")

### SC
deseq.test <- DESeqDataSetFromMatrix(bulk_sc[[1]],
                                     colData = DataFrame(row.names = colnames(bulk_sc[[1]]), Strain = c(rep("B6", 6), rep("CAST", 6))), 
                                     design = ~Strain)
deseq.test <- DESeq(deseq.test)
degs.test <- results(deseq.test)

degs.test %>% data.frame() %>%
  dplyr::filter(baseMean > 50) -> test
table(test$padj < 0.1)
perc <- paste0(round(table(test$padj < 0.1)[[2]] / sum(table(test$padj < 0.1)), digits = 4) * 100, "%")
perc_2 <- paste0(round(table(test$padj < 0.1 & abs(test$log2FoldChange > 0.5))[[2]] / sum(table(test$padj < 0.1)), digits = 4) * 100, "%")

make_volcano_plot(degs.test) + ggtitle("SC")
ggsave("~/Desktop/Projects/ASE_Spermatogenesis_Paper_rerun_Revisions/Plots/FigureS1_rev/volcano_plot_sc.pdf")

### RS
deseq.test <- DESeqDataSetFromMatrix(bulk_rs[[1]],
                                     colData = DataFrame(row.names = colnames(bulk_rs[[1]]), Strain = c(rep("B6", 6), rep("CAST", 6))), 
                                     design = ~Strain)
deseq.test <- DESeq(deseq.test)
degs.test <- results(deseq.test)

degs.test %>% data.frame() %>%
  dplyr::filter(baseMean > 50) -> test
table(test$padj < 0.1)
perc <- paste0(round(table(test$padj < 0.1)[[2]] / sum(table(test$padj < 0.1)), digits = 4) * 100, "%")
perc_2 <- paste0(round(table(test$padj < 0.1 & abs(test$log2FoldChange > 0.5))[[2]] / sum(table(test$padj < 0.1)), digits = 4) * 100, "%")

make_volcano_plot(degs.test) + ggtitle("RS")
ggsave("~/Desktop/Projects/ASE_Spermatogenesis_Paper_rerun_Revisions/Plots/FigureS1_rev/volcano_plot_rs.pdf")

### ES
deseq.test <- DESeqDataSetFromMatrix(bulk_es[[1]],
                                     colData = DataFrame(row.names = colnames(bulk_es[[1]]), Strain = c(rep("B6", 6), rep("CAST", 6))), 
                                     design = ~Strain)
deseq.test <- DESeq(deseq.test)
degs.test <- results(deseq.test)

degs.test %>% data.frame() %>%
  dplyr::filter(baseMean > 50) -> test
table(test$padj < 0.1)
perc <- paste0(round(table(test$padj < 0.1)[[2]] / sum(table(test$padj < 0.1)), digits = 4) * 100, "%")
perc_2 <- paste0(round(table(test$padj < 0.1 & abs(test$log2FoldChange > 0.5))[[2]] / sum(table(test$padj < 0.1)), digits = 4) * 100, "%")

make_volcano_plot(degs.test) + ggtitle("ES")
ggsave("~/Desktop/Projects/ASE_Spermatogenesis_Paper_rerun_Revisions/Plots/FigureS1_rev/volcano_plot_es.pdf")


### 
angela_categorization <- readxl::read_excel("~/Downloads/SuppTables.xls", sheet = 6) %>%
  column_to_rownames("geneID")
angela_categorization <- angela_categorization[!duplicated(angela_categorization$geneName), ]
rownames(angela_categorization) <- angela_categorization$geneName

model.selection.results.bulk <- readRDS("./Data/processed/model_selection_bulk.rds")
rownames(model.selection.results.bulk) <- model.selection.results.bulk$Gene

genes_in_both <- intersect(angela_categorization$geneName, rownames(model.selection.results.bulk))
genes_from_liver <- angela_categorization$geneName[angela_categorization$geneName %in% genes_in_both]

# how many genes in testis dataset are categorized?
model.selection.results.bulk$in_liver <- model.selection.results.bulk$Gene %in% genes_in_both

model.selection.results.bulk %>%
  ggplot(aes(x = "", fill = in_liver)) + geom_bar(col = "black", position = "fill") + theme_paper() + scale_fill_manual(values = c("white", "grey")) + 
  ylab("% of genes detected in liver (Goncalves et al)") + scale_y_continuous(expand = c(0, 0)) + xlab("")

model.selection.results.bulk_inliver <- model.selection.results.bulk[model.selection.results.bulk$in_liver, ]

angela_categorization$in_testis <- ifelse(angela_categorization$geneName %in% genes_in_both, "Detected", "NotDetected")
angela_categorization[genes_from_liver, ]$in_testis <- as.character(model.selection.results.bulk_inliver[genes_from_liver, ]$CategoryNew)

number_of_effects <- data.frame(table(angela_categorization$category))

angela_categorization %>%
  mutate(in_testis = factor(in_testis, levels = c("NotDetected", "conserved", "cis", "trans", "cis_trans"))) %>%
  mutate(category = factor(category, levels = c("CONS", "CIS", "TRANS", "CIS&TRANS"))) %>%
  ggplot() + geom_bar(aes(x = category, fill = in_testis), col = "black", position = "fill") + theme_paper(textsize = 30) + 
  #scale_fill_manual(values = c("white", "grey")) + 
  ylab("% of genes detected in liver (Goncalves et al)") + scale_y_continuous(expand = c(0, 0)) + xlab("") + 
  geom_text(aes(x = Var1, y = 0.95, label = Freq), data = number_of_effects, size = 20) + 
  scale_fill_manual(values = c(colours_cis_trans, "NotDetected" = "white"))
ggsave("~/Desktop/Projects/ASE_Spermatogenesis_Paper_rerun_Revisions/Plots/FigureS1_rev/comparison_to_goncalves.pdf")

number_of_effects <- angela_categorization %>% 
  dplyr::filter(in_testis != "NotDetected") %>% 
  .$category %>% table() %>% data.frame()

angela_categorization %>%
  dplyr::filter(!in_testis == "NotDetected") %>%
  mutate(in_testis = factor(in_testis, levels = c("conserved", "cis", "trans", "cis_trans"))) %>%
  mutate(category = factor(category, levels = c("CONS", "CIS", "TRANS", "CIS&TRANS"))) %>%
  ggplot() + geom_bar(aes(x = category, fill = in_testis), col = "black", position = "fill") + theme_paper(textsize = 30) + 
  #scale_fill_manual(values = c("white", "grey")) + 
  geom_text(aes(x = ., y = 0.95, label = Freq), data = number_of_effects, size = 20) + 
  ylab("% of genes detected in liver (Goncalves et al)") + scale_y_continuous(expand = c(0, 0)) + xlab("") + 
  scale_fill_manual(values = c(colours_cis_trans))
ggsave("~/Desktop/Projects/ASE_Spermatogenesis_Paper_rerun_Revisions/Plots/FigureS1_rev/comparison_to_goncalves_only_detected.pdf")

