### ATAC data analysis for the paper

####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  
####  Read data, integrate peaks, check scatters
####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  ####  

library(tidyverse)
library(GenomicRanges)
library(scran)

theme_paper <- function(textsize = 20){
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1), 
        text = element_text(size = textsize))
}


# Read in count files
get_annotated_data <- function(path1, path2){
  data_ref <- read_table(path1, col_names = F)
  data_alt <- read_table(path2, col_names = F)
  
  data_full <- data.frame(
    peak = data_ref$X1, 
    data_ref = data_ref$X2, 
    data_alt = data_alt$X2
  )
  rownames(data_full) <- data_full$peak
  
  ranges <- data.frame(do.call('rbind', lapply(data_full[,1], function(x){unlist(str_split(x, pattern = "_"))})))
  ranges_df <- data.frame(
    Chr = ranges$X1, 
    Start = as.numeric(ranges$X2), 
    End = as.numeric(ranges$X3)
  )
  rownames(ranges_df) <- data_full$peak
  ranges_granges <- makeGRangesFromDataFrame(ranges_df[!is.na(ranges_df$Start), ])
  
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(org.Mm.eg.db)
  txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
  ranges_granges_annotated <- ChIPseeker::annotatePeak(ranges_granges, TxDb = txdb)
  ranges_granges_annotated <- ranges_granges_annotated@anno
  ranges_granges_annotated$CountsRef = data_full[names(ranges_granges_annotated), ][,2]
  ranges_granges_annotated$CountsAlt = data_full[names(ranges_granges_annotated), ][,3]
  #ranges_granges_annotated$CountsRef = data_full[,2]
  #ranges_granges_annotated$CountsAlt = data_full[,3]
  ranges_granges_annotated
}

# data_846663 <- get_annotated_data("~/Desktop/ATAC_test/AS-846663-LR-64674_counts_allele1.txt.gz", "~/Desktop/ATAC_test/AS-846663-LR-64674_counts_allele2.txt.gz") # 1N
data_846665 <- get_annotated_data("~/Desktop/ATAC_test/AS-846665-LR-64674_counts_allele1.txt.gz", "~/Desktop/ATAC_test/AS-846665-LR-64674_counts_allele2.txt.gz") # 4N
# data_846667 <- get_annotated_data("~/Desktop/ATAC_test/AS-846667-LR-64674_counts_allele1.txt.gz", "~/Desktop/ATAC_test/AS-846667-LR-64674_counts_allele2.txt.gz") # 1N
data_846669 <- get_annotated_data("~/Desktop/ATAC_test/AS-846669-LR-64674_counts_allele1.txt.gz", "~/Desktop/ATAC_test/AS-846669-LR-64674_counts_allele2.txt.gz") # 4N

# compare SCs and STs - we don't need this necessarily
overlapping_peaks <- findOverlaps(data_846665, data_846667, minoverlap = 0.8)

compare_st_sc <- data.frame(
  SC_ref = data_846665[from(overlapping_peaks), ]$CountsRef, 
  SC_alt = data_846665[from(overlapping_peaks), ]$CountsAlt, 
  ST_ref = data_846667[to(overlapping_peaks), ]$CountsRef, 
  ST_alt = data_846667[to(overlapping_peaks), ]$CountsAlt
)

size_factors <- c(sum(compare_st_sc$ST_ref + compare_st_sc$ST_alt), sum(compare_st_sc$SC_ref + compare_st_sc$SC_alt))
size_factors <- size_factors / mean(size_factors)

compare_st_sc %>%
  ggplot(aes(x = (SC_ref + SC_alt), y = (ST_ref + ST_alt))) + geom_point(size = 0.1) + scale_x_log10() + scale_y_log10() + theme_paper() + geom_abline() + 
  xlab("Spermatocytes") + ylab("Spermatids")

compare_st_sc %>%
  ggplot(aes(x = (SC_ref + SC_alt) / size_factors[[1]], y = (ST_ref + ST_alt) / size_factors[[2]])) + 
    geom_point(size = 0.1) + scale_x_log10() + scale_y_log10() + theme_paper() + geom_abline() + 
    xlab("Spermatocytes") + ylab("Spermatids")

# compare SCs and SCs
overlapping_peaks <- findOverlaps(data_846665, data_846669, minoverlap = 0.9)

compare_sc_sc <- data.frame(
  SC_1_ref = data_846665[from(overlapping_peaks), ]$CountsRef, 
  SC_1_alt = data_846665[from(overlapping_peaks), ]$CountsAlt, 
  SC_2_ref = data_846669[to(overlapping_peaks), ]$CountsRef, 
  SC_2_alt = data_846669[to(overlapping_peaks), ]$CountsAlt, 
  chromosome = seqnames(data_846665[from(overlapping_peaks), ])
) %>% mutate(chromosome %in% c("chrX", "chrMT"))

size_factors <- c(sum(compare_sc_sc$SC_1_ref + compare_sc_sc$SC_1_alt), sum(compare_sc_sc$SC_2_ref + compare_sc_sc$SC_2_alt))
size_factors <- size_factors / mean(size_factors)

compare_sc_sc %>%
  ggplot(aes(x = (SC_1_ref + SC_1_alt) / size_factors[[1]], y = (SC_2_ref + SC_2_alt) / size_factors[[2]])) + 
  geom_point(size = 0.1) + scale_x_log10() + scale_y_log10() + theme_paper(textsize = 40) + geom_abline() + 
  xlab("Spermatocytes (Rep 1) (log10 reads)") + ylab("Spermatocytes (Rep 2) (log10 reads)") + coord_fixed()
ggsave("~/Desktop/Projects/ASE_Spermatogenesis_Paper_rerun_Revisions/Plots/FigureS3_rev/ATAC_scatter_replicates.pdf")

compare_sc_sc %>%
  dplyr::filter((SC_1_ref + SC_1_alt + SC_2_ref + SC_2_alt) / 4 > 50) %>%
  ggplot(aes(x = SC_1_ref / (SC_1_ref + SC_1_alt), y = SC_2_ref / (SC_2_ref + SC_2_alt))) + 
  geom_point(size = 0.1) + theme_paper(textsize = 40) + geom_abline() + 
  xlab("Spermatocytes (Rep 1) (allelic ratio)") + ylab("Spermatocytes (Rep 2) (allelic ratio)") + coord_fixed()
ggsave("~/Desktop/Projects/ASE_Spermatogenesis_Paper_rerun_Revisions/Plots/FigureS3_rev/ATAC_scatter_replicates_ase.pdf")

# Define joint sites between two samples
compare_sc_sc <- mergeByOverlaps(data_846665, data_846669, minoverlap = 0.9)

data_full <- data.frame(
  chromosome = paste0("chr", compare_sc_sc$geneChr), 
  start = compare_sc_sc$geneStart, 
  end = compare_sc_sc$geneEnd,
  Gene = compare_sc_sc$geneId, 
  Annotation = compare_sc_sc$annotation,
  DistanceToTSS = compare_sc_sc$distanceToTSS, 
  Ref_1 = compare_sc_sc$data_846665$CountsRef, 
  Ref_2 = compare_sc_sc$data_846669$CountsRef, 
  Alt_1 = compare_sc_sc$data_846665$CountsAlt, 
  Alt_2 = compare_sc_sc$data_846669$CountsAlt
) %>%
  mutate(data_ref = Ref_1 + Ref_2) %>%
  mutate(data_alt = Alt_1 + Alt_2) %>%
  mutate(total = data_ref + data_alt) %>%
  mutate(ASE = data_ref / (data_ref + data_alt))
data_full <- data_full[!duplicated(data_full), ]

rbind(
  data.frame(
    Annotation = gsub(" \\(.*\\)", "", data_846665$annotation), 
    Count = data_846665$CountsRef + data_846665$CountsRef,
    Sample = "Rep1" 
  ), 
  data.frame(
    Annotation = gsub(" \\(.*\\)", "", data_846669$annotation), 
    Count = data_846669$CountsRef + data_846669$CountsRef,
    Sample = "Rep2"
  )
) %>%
  dplyr::filter(Count > 2 * 50) %>%
  ggplot(aes(x = Sample, fill = Annotation)) + geom_bar(position = "fill") + theme_paper(textsize = 40) + 
  ylab("Number of Peaks")
ggsave("~/Desktop/Projects/ASE_Spermatogenesis_Paper_rerun_Revisions/Plots/FigureS3_rev/peak_distribution_per_sample.pdf", width = 10)

data.frame(
  Annotation = gsub(" \\(.*\\)", "", data_full$Annotation), 
  Count = data_full$total
) %>%
  dplyr::filter(Count > 2 * 50) %>%
  ggplot(aes(x = 1, fill = Annotation)) + geom_bar(position = "fill") + theme_paper(textsize = 30) + 
  xlab("Binned expression (log10(Count + 1))") + ylab("Number of Peaks")

data_full$annotation <- gsub("Intron.*", "Intron", data_full$Annotation)
data_full$annotation <- gsub("Exon.*", "Exon", data_full$Annotation)
data_full$symbol <- mapIds(org.Mm.eg.db, keys = gsub("\\.[0-9]*", "", data_full$Gene), column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")


# For the following analysis, exclude peaks more than 20kb away from any gene
data_full <- data_full %>%
  dplyr::filter(abs(DistanceToTSS) < 20000)

#data_full <- data_full[!is.na(data_full$ASE), ]
#data_full <- data_full[data_full$data_ref + data_full$data_alt > 100, ]

data_results <- readRDS("~/Desktop/Projects/ASE_Spermatogenesis_Paper_rerun_Revisions//Data/processed/genes_test.rds")

data_full <- data_full %>%
  # dplyr::filter(!grepl("Distal", Annotation)) %>%
  dplyr::filter(total > 100)

test <- data_full %>%
  data.frame() %>%
  dplyr::select(c(symbol, ASE)) %>%
  mutate(ASE = ASE) %>%
  mutate(ASE_effect = abs(ASE - 0.5)) %>%
  group_by(symbol) %>%
  summarize(ASE = ASE[which(max(ASE_effect) == ASE_effect)[[1]]]) %>%
  mutate(symbol = unlist(symbol)) %>%
  dplyr::filter(!is.na(symbol)) %>%
  column_to_rownames("symbol") %>%
  mutate(ATAC_effect = abs(ASE - 0.5) > 0.1)

data_results$ATAC_effect <- test[rownames(data_results), ]$ATAC_effect
data_results$ATAC_ASE <- test[rownames(data_results), ]$ASE

# how many sites are close to genes we are interested in?
dim(data_full[data_full$symbol %in% rownames(data_results), ])


data_results_here <- data_results %>%
  dplyr::filter(!is.na(data_results$ATAC_ASE))

dim(data_results_here)

set.seed(12345)
data_results_here$ATAC_shuffled <- data_results_here$ATAC_ASE[sample(1:nrow(data_results_here), nrow(data_results_here))]

data_results_here <- data_results_here[data_results_here$padj < 0.1, ]

lapply(c(0.1, 0.2, 0.3, 0.4), function(cutoff){
  
  tt_shuff <- table(abs(data_results_here[abs(data_results_here$ase - 0.5) > cutoff, ]$ATAC_shuffled - 0.5) > 0.1) %>% data.frame()
  tt_shuff$Cov <- "Shuffled"
  tt_shuff$Group <- c("No ATAC", "ATAC")
  tt_shuff$Cutoff <- cutoff
  tt_shuff$Freq <- tt_shuff$Freq / sum(tt_shuff$Freq)
    
  tt <- table(abs(data_results_here[abs(data_results_here$ase - 0.5) > cutoff, ]$ATAC_ASE - 0.5) > 0.1) %>% data.frame()
  tt$Cov <- "Real"
  tt$Group <-  c("No ATAC", "ATAC")
  tt$Cutoff <- cutoff
  tt$Freq <- tt$Freq / sum(tt$Freq)
  
  rbind(tt_shuff, tt)
  
}) %>% do.call("rbind", .) -> results_df

results_df %>%
  dplyr::filter(!grepl("No ATAC", Group)) %>%
  mutate(Cov = factor(Cov, levels = c("Shuffled", "Real"))) %>%
  mutate(Cutoff = paste0(">", Cutoff)) %>%
  ggplot(aes(x = Cutoff, fill = Cov, y = Freq)) + 
  geom_bar(stat = "identity", position = "dodge", col = "black") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + theme_paper(textsize = 35) + xlab("Cutoff AI RNA") + ylab("Fraction of genes with AI in CA > 0.1") + 
    scale_fill_manual(values = c("white", "black")) + scale_y_continuous(expand = c(0, 0), limits = c(0, 1.1)) + labs("") + labs(fill = "")
ggsave("~/Desktop/Projects/ASE_Spermatogenesis_Paper_rerun_Revisions/Plots/FigureS3_rev/ATAC_enrichment_rna_bulk.pdf")

NULL


## get ratios per bin 
setwd("~/Desktop/Projects/ASE_Spermatogenesis_Paper_rerun/")
data <- readRDS("./Data/processed/sce_merged_new.rds")
data$IndividualSamples <- paste0(data$Dataset, "_", data$Library)
data <- data[,data$Library != "Sample7"]
aggregate_test <- aggregateAcrossCells(data, DataFrame(Individual = data$IndividualSamples, CellType = data$CellType), 
                                       use.assay.type = c("counts_reference", "counts_alternative"))

# ask whether enrichment is stronger in scs?

lapply(c("SC", "RS", "ES"), function(ct){
  aggregate_test_celltype <- aggregate_test[,aggregate_test$CellType == ct]
  ases_here <- rowSums(assays(aggregate_test_celltype)[["counts_reference"]]) / 
    (rowSums(assays(aggregate_test_celltype)[["counts_reference"]] + assays(aggregate_test_celltype)[["counts_alternative"]]))
}) %>% do.call("cbind", .) -> ratios_per_celltype
colnames(ratios_per_celltype) <- c("SC", "RS", "ES")

data_results_here <- data_results %>%
  dplyr::filter(!is.na(data_results$ATAC_ASE))

set.seed(1234)
data_results_here$ATAC_shuffled <- data_results_here$ATAC_ASE[sample(1:nrow(data_results_here), nrow(data_results_here))]
data_results_here$sign <- sign(data_results_here$ATAC_shuffled - 0.5) == sign(data_results_here$ase - 0.5)

data_results_here <- cbind(data_results_here, ratios_per_celltype[rownames(data_results_here), ])

get_per_celltype <- function(celltype){lapply(c(0.1, 0.2, 0.3, 0.4), function(cutoff){
  
  data_results_here$ase_here <- data_results_here[,celltype]
  
  tt_shuff <- table(abs(data_results_here[abs(data_results_here$ase_here - 0.5) > cutoff, ]$ATAC_shuffled - 0.5) > 0.1) %>% data.frame()
  tt_shuff$Cov <- "Shuffled"
  tt_shuff$Group <- c("No ATAC", "ATAC")
  tt_shuff$Cutoff <- cutoff
  tt_shuff$Freq <- tt_shuff$Freq / sum(tt_shuff$Freq)
  tt_shuff$CellType <- celltype
  
  tt <- table(abs(data_results_here[abs(data_results_here$ase_here - 0.5) > cutoff, ]$ATAC_ASE - 0.5) > 0.1) %>% data.frame()
  tt$Cov <- "Real"
  tt$Group <-  c("No ATAC", "ATAC")
  tt$Cutoff <- cutoff
  tt$Freq <- tt$Freq / sum(tt$Freq)
  tt$CellType <- celltype
  
  rbind(tt_shuff, tt)
  
}) %>% do.call("rbind", .)}

results_df_sc <- get_per_celltype("SC")
results_df_rs <- get_per_celltype("RS")
results_df_es <- get_per_celltype("ES")

rbind(cbind(results_df, CellType = "Bulk"), results_df_sc, results_df_rs, results_df_es)  %>%
  dplyr::filter(!grepl("No ATAC", Group)) %>%
  mutate(Cov = factor(Cov, levels = c("Shuffled", "Real"))) %>%
  mutate(Cutoff = paste0(">", Cutoff)) %>%
  ggplot(aes(x = Cutoff, fill = interaction(Cov, CellType), y = Freq)) + 
  geom_bar(stat = "identity", position = "dodge", col = "black", width = 0.8) + 
  #facet_wrap(~Cutoff, nrow = 1) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + theme_paper(textsize=  40) + ylab("Fraction of genes with CA > 0.1") + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) + 
  labs(fill = "") + 
  scale_fill_manual(values = c("white", "black", "grey92", "grey60", "#bcbee8", "#282B69", "#f0bdda", "#D31281"))

rbind(results_df_sc, results_df_rs, results_df_es)  %>%
  dplyr::filter(!grepl("No ATAC", Group)) %>%
  mutate(Cov = factor(Cov, levels = c("Shuffled", "Real"))) %>%
  mutate(Cutoff = paste0(">", Cutoff)) %>%
  ggplot(aes(x = Cutoff, fill = interaction(Cov, CellType), y = Freq)) + 
    geom_bar(stat = "identity", position = "dodge", col = "black", width = 0.8) + 
    #facet_wrap(~Cutoff, nrow = 1) + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + theme_paper(textsize=  40) + ylab("Fraction of genes with CA > 0.1") + 
    scale_y_continuous(expand = c(0, 0), limits = c(0, 1)) + 
    labs(fill = "") + 
    scale_fill_manual(values = c("grey92", "grey60", "#bcbee8", "#282B69", "#f0bdda", "#D31281"))
ggsave("~/Desktop/Projects/ASE_Spermatogenesis_Paper_rerun_Revisions/Plots/FigureS3_rev/ATAC_enrichment_rna_per_celltype.pdf")

## We now look at allelic imbalance of dynamic genes with AI_ATAC in spermatocytes
### read processed data and generate allelic ratios
list_data <- readRDS("~/Desktop/Projects/ASE_Spermatogenesis_Paper_rerun_Revisions//Data/processed/binned_data.rds")

data_f0_b6 <- list_data[[1]]
data_f0_cast <- list_data[[2]]
data_f1_b6 <- list_data[[3]]
data_f1_cast <- list_data[[4]]

sfs = c(sum(data_f0_b6), sum(data_f0_cast)) / sum(data_f0_b6)
data_f0_b6_norm <- data_f0_b6 / sfs[[1]]
data_f0_cast_norm <- data_f0_cast / sfs[[2]]

data_f1_b6_scaled <- t(t(data_f1_b6) / colSums(data_f1_b6))
data_f1_cast_scaled <- t(t(data_f1_cast) / colSums(data_f1_cast))

genes_use <- rownames(data_f0_b6_norm)
rs_f1 <- (data_f1_b6[genes_use, ]) / (data_f1_b6[genes_use, ] + data_f1_cast[genes_use, ])
rs_f1_dynamic <- rs_f1[rownames(data_results[p.adjust(data_results$dali_pval_linear) < 0.1, ]), ]

data.frame(
  ase_sc = (rowMeans(abs(rs_f1_dynamic[, 5:20] - 0.5), na.rm = T)), 
  ase_non_sc = (rowMeans(abs(rs_f1_dynamic[, 21:100] - 0.5), na.rm = T))
) %>% 
  rownames_to_column("Gene") %>% 
  add_column(ATAC_ASE = abs(data_results[.$Gene, ]$ATAC_ASE - 0.5)) %>%
  mutate(ase_difference = ase_sc - ase_non_sc) %>%
  mutate(ase_difference_interval = cut(ase_difference, 5)) -> test

colors_here <- c("grey", rev(colorRampPalette(c("#D31281", "grey"))(4)))
counts_df <- data.frame(table(test$ase_difference_interval))

my_comparisons <- lapply(unique(test$ase_difference_interval)[unique(test$ase_difference_interval) != "(0.213,0.347]"], function(x){c(as.character(x), "(0.213,0.347]")})
test %>%
  ggplot(aes(ase_difference_interval, y = ATAC_ASE, fill = ase_difference_interval)) + geom_violin() + theme_paper(40) + stat_summary() + 
    xlab("AI increase in Spermatocytes") + ylab("ATAC AI (Spermatocytes)") + 
    # ggpubr::stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + 
    scale_fill_manual(values = c(colors_here)) + 
    theme(legend.position = "None") + 
    annotate(label = counts_df$Freq, x = counts_df$Var1, y = 0.55, geom = "text", size = 15) + 
    coord_flip() + ylim(0, 0.6)
ggsave("~/Desktop/Projects/ASE_Spermatogenesis_Paper_rerun_Revisions/Plots/FigureS3_rev/ATAC_dynamic_analysis.pdf")

pvals <- lapply(unique(test$ase_difference_interval), function(x){
  data_up = test[test$ase_difference_interval == "(0.213,0.347]", ]
  data_other = test[test$ase_difference_interval == x, ]
  wilcox.test(data_up$ATAC_ASE, data_other$ATAC_ASE)$p.value
})
pvals <- pvals[pvals < 1]

unique(test$ase_difference_interval)
p.adjust(unlist(pvals))

