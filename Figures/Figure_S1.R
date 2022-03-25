## Supplementary plots Figure 1

library(tidyverse)
library(reshape2)
library(scran)
library(scater)

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

theme_paper <- function(textsize = 20){
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1), 
        text = element_text(size = textsize))
}

counts_reference <- function(sce){
  assays(sce)[["counts_reference"]]
}

counts_alternative <- function(sce){
  assays(sce)[["counts_alternative"]]
}

sce.all <- readRDS("~/Desktop/Projects/ASE_Spermatogenesis_Paper/Data/processed/final_sce_f1_dataset_for_supplement.rds")
sce.all <- sce.all[,!grepl("Outliers", sce.all$AnnotatedClusters)]

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

# 

# recompute a clean umap without outlier cells
set.seed(12345)
corrected <- metadata(sce.all)$corrected[colnames(sce.all), ]
umap <- umap::umap(corrected)
reducedDims(sce.all)$UMAP <- umap$layout

# Supp 1_1:
stats.df <- read.csv("~/Desktop/Projects/ASE_Spermatogenesis_Paper/Data/processed/filtering_stats.csv", row.names = 1)

stats.df %>% 
  dplyr::rename(BeforeFiltering = No_cells) %>%
  pivot_longer(-c(Sample, Library)) %>%
  ggplot(aes(x = Library, y = value, fill = Sample, group = name)) + 
    geom_bar(stat = "identity", position = "dodge") + 
    geom_point(aes(x = Library, y = value + 300, col = name, group = name), position = position_dodge(width = 1), size = 10) + 
    theme_classic() + 
    coord_flip() + 
    scale_color_manual(values = c("black", "grey")) + ylab("Number of cells") + 
    scale_fill_manual(values = sample_colors) + 
    theme_paper(textsize = 50) + xlab("")

ggsave("~/Desktop/Projects/ASE_Spermatogenesis_Paper/Plots/FigureS1/FigS1_1.pdf")

# Supp 1_2: 

data_coverage <- data.frame(
  Library = colData(sce.all)$Library, 
  Sample = colData(sce.all)$Sample, 
  Library_Sample = paste0(colData(sce.all)$Sample, "_", colData(sce.all)$Library), 
  detected = colData(sce.all)$detected
)

ggplot(data_coverage, aes(y = detected, x = Library, fill = Library_Sample)) + 
  geom_violin() + scale_fill_manual(values = library_colors) + theme_classic() + 
  theme(text = element_text(size = 20)) + ylab("Number of detected genes") + xlab("") + 
  coord_flip() + theme_paper(textsize = 50) + 
  theme(legend.position = "None") + theme(aspect.ratio = 1.5)

ggsave("~/Desktop/Projects/ASE_Spermatogenesis_Paper/Plots/FigureS1/FigS1_2.pdf")


# Supp 1_3:

data_tsne <- data.frame(
  Library = colData(sce.all)$Library, 
  Sample = colData(sce.all)$Sample, 
  Library_Sample = paste0(colData(sce.all)$Sample, "_", colData(sce.all)$Library), 
  Umap1 = reducedDims(sce.all)$UMAP[,1], 
  Umap2 = reducedDims(sce.all)$UMAP[,2]
)

ggplot(data_tsne, aes(Umap1, Umap2, col = Library_Sample)) + 
  geom_point(alpha = 0.2) + theme_tsne() + facet_wrap(~Sample, nrow = 2) +
  scale_color_manual(values = library_colors) + coord_fixed() + 
  theme(text = element_text(size= 50), legend.position = "None")

ggsave("~/Desktop/Projects/ASE_Spermatogenesis_Paper/Plots/FigureS1/FigS1_3.pdf")

# Supp 1_4: 

colour_vector <- c("Spermatogonia" = "#046735",
                   "Spermatogonia / Leptotene" = "#FBDCEA", 
                   "Pachytene" = "#D31281",
                   "Diplotene" = "#7E2678",
                   "Meiosis" = "#502269",
                   "Round Spermatids (1)" = "#B3A9D3",
                   "Round Spermatids (2)" = "#D8DAEC",
                   "Round Spermatids (3)" = "#9595C9",
                   "Round Spermatids (4)" = "#6167AF",
                   "Round Spermatids (5)" = "#2A348B",
                   "Elongating Spermatids (1)" = "#171447",
                   "Elongating Spermatids (2)" = "#22232B",
                   "Elongating Spermatids (3)" = "#444A5E",
                   "Elongating Spermatids (4)" = "#606B89",
                   "Elongating Spermatids (5)" = "#737884",
                   "Elongating Spermatids (6)" = "#737884",
                   "Elongating Spermatids (7)" = "#737884",
                   "Sertoli" = "#643F18",
                   "Leydig" = "#985D25",
                   "Immune" = "#C48C58")

data.frame(
  Umap1 = reducedDims(sce.all)$UMAP[,1],
  Umap2 = reducedDims(sce.all)$UMAP[,2],
  Cluster = sce.all$AnnotatedClustersFine
) %>%
  ggplot(aes(Umap1, Umap2, col = Cluster)) + 
  scale_color_manual(values = colour_vector) + 
  geom_point(size = 1) + 
  theme_paper(textsize = 40) + 
  theme(legend.title = element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=10)))

ggsave("~/Desktop/Projects/ASE_Spermatogenesis_Paper/Plots/FigureS1/FigS1_4.pdf")

# Supp 1_5: 

marker_genes <- c("Dazl", "Hormad1", "Piwil1", "Pou5f2", "Tex21", "Prm1", "Cst12", "Insl3", "Acta2")

data_genes <- cbind(2 ** data.frame(as.matrix(t(logcounts(sce.all[marker_genes, ])))), 
                    Sample = colData(sce.all)$Sample, 
                    CellType = colData(sce.all)$AnnotatedClusters)

data_genes_transformed <-
  data_genes %>% pivot_longer(-c(Sample, CellType), names_to = "Gene") %>%
  group_by(Sample, CellType, Gene) %>%
  summarize(value = mean(value)) %>% ungroup() %>%
  group_by(Gene) %>%
  mutate(value = scale(value))

data_genes_transformed$CellType <- factor(data_genes_transformed$CellType, levels = c("SG", "SC", "RS", "ES", "Sertoli", "Leydig", "Immune"))
data_genes_transformed$Gene <- factor(data_genes_transformed$Gene, levels = marker_genes)

ggplot(data_genes_transformed, aes(x = CellType, y = Gene, fill = Sample, size = value, group = Sample)) + 
  geom_point(position=position_dodge(width=0.8), shape = 22, col = "black") + scale_fill_manual(values = sample_colors) + 
  theme_paper(textsize = 50) + scale_size_area() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave("~/Desktop/Projects/ASE_Spermatogenesis_Paper/Plots/FigureS1/FigS1_5.pdf")

# Supp 1_6:

data.frame(
  Library = sce.all$Library, 
  Sample = sce.all$Sample, 
  Library_Sample = paste0(sce.all$Sample, "_", sce.all$Library), 
  CellType = sce.all$AnnotatedClusters
) %>% group_by(Library_Sample, CellType) %>%
  dplyr::count() %>%
  group_by(Library_Sample) %>%
  mutate(per =  100 * n / sum(n)) %>%
ggplot(aes(x = CellType, fill = Library_Sample, y = per)) + 
  geom_bar(position = "dodge", stat = "identity") + 
  ylim(0, 50) +
  scale_fill_manual(values = library_colors) + 
  theme_paper(textsize = 40) + ylab("Cell Type Proportions [%]") + 
  coord_flip() + xlab("")

ggsave("~/Desktop/Projects/ASE_Spermatogenesis_Paper/Plots/FigureS1/FigS1_6.pdf")


# Allelic data

# Supp 1F:

sce.merged <- readRDS("~/Desktop/Projects/ASE_Spermatogenesis_Paper/Data/processed/sce_merged_new_sparse.rds")

df_total_counts <- data.frame(
  Species = sce.merged$Species, 
  Library = sce.merged$Library, 
  Counts_Ref = colSums(counts_reference(sce.merged)), 
  Counts_Alt = colSums(counts_alternative(sce.merged))
) %>% pivot_longer(-c(Species, Library)) %>%
  dplyr::rename(Allele = name) %>% dplyr::rename(Counts = value)

average.reference <- df_total_counts %>%
  dplyr::filter(Allele == "Counts_Ref") %>%
  group_by(Library) %>%
  summarize(sum = sum(Counts))
average.alternative <- df_total_counts %>%
  dplyr::filter(Allele == "Counts_Alt") %>%
  group_by(Library) %>%
  summarize(sum = sum(Counts))
ratios <- round(average.reference$sum / (average.reference$sum + average.alternative$sum), digits = 3) * 100
  
ggplot(df_total_counts, aes(x = Library, y = Counts + 1, fill = Allele)) + 
  geom_boxplot() + 
  ylab("log(ReadCounts + 1) per cell") + 
  xlab("") + 
  scale_y_log10(limits = c(1, 1.5e5)) + 
  annotate("text", label = paste0(as.character(ratios), "%"), 
           x = 1:6, y = 60000, size = 10) + 
  annotate("text", label = "Percent reference reads: ", 
           x = 2.6, y = 140000, size = 12) + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
  theme_paper(textsize = 40)

ggsave("~/Desktop/Projects/ASE_Spermatogenesis_Paper/Plots/FigureS1/FigS1_7.pdf")

# Plot aggregated data across cells to compare on gene level
sce_sample1 <- sce.merged[,sce.merged$Library == "Sample1"]
df_plot_1 <- data.frame(
  Counts_Reference = rowMeans(counts_reference(sce_sample1)),
  Counts_Alternative = 0
)
df_plot_1[names(rowMeans(counts_reference(sce_sample1))),]$Counts_Alternative <- 
  rowMeans(counts_alternative(sce_sample1))

sce_sample3 <- sce.merged[,sce.merged$Library == "Sample3"]
df_plot_3 <- data.frame(
  Counts_Reference = rowMeans(counts_reference(sce_sample3)),
  Counts_Alternative = 0
)
df_plot_3[names(rowMeans(counts_reference(sce_sample3))),]$Counts_Alternative <- 
  rowMeans(counts_alternative(sce_sample3))

sce_sample5 <- sce.merged[,sce.merged$Library == "Sample5"]
sce_sample5 <- annotate_chromosome(sce_sample5)
df_plot_5 <- data.frame(
  Counts_Reference = rowMeans(counts_reference(sce_sample5)),
  Counts_Alternative = 0,
  Chromosome = rowData(sce_sample5)$chromosome_name %in% c("MT", "X")
)
df_plot_5[names(rowMeans(counts_reference(sce_sample5))),]$Counts_Alternative <- 
  rowMeans(counts_alternative(sce_sample5))

p1 <- ggplot(df_plot_1, aes(log(Counts_Reference + 1), log(Counts_Alternative + 1))) + 
  geom_point(size = 0.5) + 
  xlim(0, 4) + ylim(0, 4) + 
  coord_fixed() + 
  theme_classic() + 
  annotate("text", label = "B6 (F0)", 
           x = 0.7, y = 3.9, size = 15) + 
  xlab("log( counts B6 + 1)") + 
  ylab("log( counts CAST + 1)") + 
  theme_paper(textsize = 50)

p2 <- ggplot(df_plot_3, aes(log(Counts_Reference + 1), log(Counts_Alternative + 1))) + 
  geom_point(size = 0.5) + 
  xlim(0, 4) + ylim(0, 4) + 
  coord_fixed() + 
  theme_classic() + 
  annotate("text", label = "CAST (F0)", 
           x = 1, y = 3.9, size = 15) + 
  xlab("log( counts B6 + 1)") + 
  ylab("log( counts CAST + 1)") + 
  theme_paper(textsize = 50)

p3 <- ggplot(df_plot_5, aes(log(Counts_Reference + 1), log(Counts_Alternative + 1))) + 
  geom_point(size = 0.5) + 
  xlim(0, 4) + ylim(0, 4) + 
  coord_fixed() + 
  theme_classic() + 
  annotate("text", label = "B6 x CAST (F1)", 
           x = 1.4, y = 3.9, size = 15) + 
  xlab("log( counts B6 + 1)") + 
  ylab("log( counts CAST + 1)") + 
  theme_paper(textsize = 50)

ggsave("~/Desktop/Projects/ASE_Spermatogenesis_Paper/Plots/FigureS1/FigS1_8a.pdf", p1)
ggsave("~/Desktop/Projects/ASE_Spermatogenesis_Paper/Plots/FigureS1/FigS1_8b.pdf", p2)
ggsave("~/Desktop/Projects/ASE_Spermatogenesis_Paper/Plots/FigureS1/FigS1_8c.pdf", p3)

# Look at ratios per cell type

sce.merged.cut <- sce.merged[!grepl("^mt-", rownames(sce.merged)), ]
sce.merged.cut <- sce.merged[!rowData(sce.merged)$chromosome_name %in% c("MT"), ]

df_total_counts_f1 <- data.frame(
  Library = sce.merged.cut$Library, 
  CellType = factor(sce.merged.cut$CellType, levels = c("SG", "SC", "RS", "ES", "Sertoli", "Leydig", "Immune")), 
  Counts_Ref = colSums(counts_reference(sce.merged.cut)), 
  Counts_Alt = colSums(counts_alternative(sce.merged.cut))
) %>% 
  dplyr::filter(Library %in% c("Sample5", "Sample6")) %>%
  add_column(ratio = .$Counts_Ref / (.$Counts_Ref + .$Counts_Alt))

ggplot(df_total_counts_f1, aes(x = CellType, y = ratio)) + 
  geom_boxplot() + 
  ylab("B6 / (B6 + CAST)") + 
  annotate("text", label = "Allelic ratios per cell", 
           x = 2, y = 0.95, size = 10) + 
  xlab("") + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, hjust=1)) + 
  geom_hline(yintercept = 0.5, linetype = "dashed") + 
  ylim(0, 1) + 
  theme_paper(textsize = 50)

ggsave("~/Desktop/Projects/ASE_Spermatogenesis_Paper/Plots/FigureS1/FigS1_9.pdf")

## goncalves fitting
# analysis of goncalves fitting procedure

library(ggplot2)

fits_bulk <- readRDS("~/Desktop/Projects/ASE_Spermatogenesis_Paper/Data/processed/model_selection_bulk.rds")
fits_sc <- readRDS("~/Desktop/Projects/ASE_Spermatogenesis_Paper/Data/processed/model_selection_sc.rds")
fits_rs <- readRDS("~/Desktop/Projects/ASE_Spermatogenesis_Paper/Data/processed/model_selection_rs.rds")
fits_es <- readRDS("~/Desktop/Projects/ASE_Spermatogenesis_Paper/Data/processed/model_selection_es.rds")

# first look at bayes factor distributions and see which effects are strongest

ggplot(reshape2::melt(data.frame(
  bulk = fits_bulk$LL.conserved, 
  fits_sc = fits_bulk$LL.conserved, 
  fits_rs = fits_bulk$LL.conserved, 
  fits_es = fits_bulk$LL.conserved
)), aes(fill = variable, x = value)) + geom_histogram() + facet_wrap(~variable) + scale_x_log10()

ggplot(reshape2::melt(data.frame(
  bulk = fits_bulk$LL.cis, 
  fits_sc = fits_bulk$LL.cis, 
  fits_rs = fits_bulk$LL.cis, 
  fits_es = fits_bulk$LL.cis
)), aes(fill = variable, x = value)) + geom_histogram() + facet_wrap(~variable) + scale_x_log10()

ggplot(reshape2::melt(data.frame(
  bulk = fits_bulk$LL.trans, 
  fits_sc = fits_bulk$LL.trans, 
  fits_rs = fits_bulk$LL.trans, 
  fits_es = fits_bulk$LL.trans
)), aes(fill = variable, x = value)) + geom_histogram() + facet_wrap(~variable) + scale_x_log10()

ggplot(reshape2::melt(data.frame(
  bulk = fits_bulk$LL.cis_trans, 
  fits_sc = fits_bulk$LL.cis_trans, 
  fits_rs = fits_bulk$LL.cis_trans, 
  fits_es = fits_bulk$LL.cis_trans
)), aes(fill = variable, x = value)) + geom_histogram() + facet_wrap(~variable) + scale_x_log10()

# first look at bayes factor distributions and see which effects are strongest

ggplot(data.frame(
  conserved = fits_bulk$LL.conserved, 
  effect = fits_bulk$LL.cis, 
  col = fits_bulk$Category
), aes(conserved, effect, col = col)) + geom_point() + scale_x_log10() + scale_y_log10()

ggplot(data.frame(
  conserved = fits_bulk$LL.conserved, 
  effect = fits_bulk$LL.cis, 
  col = fits_bulk$Category
), aes(conserved - effect, fill = col)) + geom_histogram() + scale_y_log10() + facet_wrap(~col)

ggplot(data.frame(
  conserved = fits_bulk$LL.conserved, 
  effect = fits_bulk$LL.trans, 
  col = fits_bulk$Category, 
  fc = fits_bulk$FC_parental
), aes(conserved - effect, -abs(fc),  fill = col)) + geom_point() + facet_wrap(~col)

# quantify number of effects

effects_combined <- data.frame(
  Gene  = Reduce("union", list(fits_bulk$Gene, fits_sc$Gene, fits_rs$Gene, fits_es$Gene))
)
rownames(effects_combined) <- effects_combined$Gene
effects_combined$cat_bulk <- NA
effects_combined$cat_sc <- NA
effects_combined$cat_rs <- NA
effects_combined$cat_es <- NA

effects_combined[fits_bulk$Gene, ]$cat_bulk <- as.character(fits_bulk$Category)
effects_combined[fits_sc$Gene, ]$cat_sc <- as.character(fits_sc$Category)
effects_combined[fits_rs$Gene, ]$cat_rs <- as.character(fits_rs$Category)
effects_combined[fits_es$Gene, ]$cat_es <- as.character(fits_es$Category)

# compute number of genes with effects
effect_counts_bulk <- table(effects_combined$cat_bulk)
effect_counts_sc <- table(effects_combined$cat_sc)
effect_counts_rs <- table(effects_combined$cat_rs)
effect_counts_es <- table(effects_combined$cat_es)

sum(effect_counts_bulk) - effect_counts_bulk[["conserved"]]
sum(effect_counts_sc) - effect_counts_sc[["conserved"]]
sum(effect_counts_rs) - effect_counts_rs[["conserved"]]
sum(effect_counts_es) - effect_counts_es[["conserved"]]

effects_combined$has_effect_bulk <- effects_combined$cat_bulk != "conserved"
effects_combined$has_effect_celltype <- apply(effects_combined[,c("cat_sc", "cat_rs", "cat_es")], 1, function(x){any(x != "conserved")})

table(effects_combined$has_effect_bulk)
table(effects_combined$has_effect_celltype)

rbind(
  data.frame(
    Analysis = "bulk", 
    table(effects_combined$has_effect_bulk)
  ), 
  data.frame(
    Analysis = "celltype", 
    table(effects_combined$has_effect_celltype)
  )
) %>%
  ggplot(aes(x = Analysis, y = Freq, fill = Var1)) + geom_bar(stat = "identity") + 
  theme_paper()

# by this analysis, we see more effects -- but this has a multiple testing problem as oli mentioned

# do the same thing only on strong effects
# get bayes factors (-LL_alt - (-LL_hyp))

fits_bulk$bf_cis <- fits_bulk$LL.conserved - fits_bulk$LL.cis
fits_bulk$bf_trans <- fits_bulk$LL.conserved - fits_bulk$LL.trans
fits_bulk$bf_cis_trans <- fits_bulk$LL.conserved - fits_bulk$LL.cis_trans
fits_bulk$combined_bf <- apply(fits_bulk[,c("bf_cis", "bf_trans", "bf_cis_trans")], 1, function(x){max(x)})

fits_sc$bf_cis <- fits_sc$LL.conserved - fits_sc$LL.cis
fits_sc$bf_trans <- fits_sc$LL.conserved - fits_sc$LL.trans
fits_sc$bf_cis_trans <- fits_sc$LL.conserved - fits_sc$LL.cis_trans
fits_sc$combined_bf <- apply(fits_sc[,c("bf_cis", "bf_trans", "bf_cis_trans")], 1, function(x){max(x)})

fits_rs$bf_cis <- fits_rs$LL.conserved - fits_rs$LL.cis
fits_rs$bf_trans <- fits_rs$LL.conserved - fits_rs$LL.trans
fits_rs$bf_cis_trans <- fits_rs$LL.conserved - fits_rs$LL.cis_trans
fits_rs$combined_bf <- apply(fits_rs[,c("bf_cis", "bf_trans", "bf_cis_trans")], 1, function(x){max(x)})

fits_es$bf_cis <- fits_es$LL.conserved - fits_es$LL.cis
fits_es$bf_trans <- fits_es$LL.conserved - fits_es$LL.trans
fits_es$bf_cis_trans <- fits_es$LL.conserved - fits_es$LL.cis_trans
fits_es$combined_bf <- apply(fits_es[,c("bf_cis", "bf_trans", "bf_cis_trans")], 1, function(x){max(x)})

bf_cutoff <- 10
effects_combined_strict <- effects_combined[,1:5]

effects_combined_strict[fits_bulk$combined_bf < bf_cutoff, ]$cat_bulk <- "conserved"
effects_combined_strict[fits_sc$combined_bf < bf_cutoff, ]$cat_sc <- "conserved"
effects_combined_strict[fits_rs$combined_bf < bf_cutoff, ]$cat_rs <- "conserved"
effects_combined_strict[fits_es$combined_bf < bf_cutoff, ]$cat_es <- "conserved"

effects_combined_strict$has_bulk_effect <- effects_combined_strict$cat_bulk != "conserved"
effects_combined_strict$has_celltype_effect <- apply(effects_combined_strict[,c("cat_sc", "cat_rs", "cat_es")], 1, function(x){any(x != "conserved")})

table(effects_combined_strict$has_bulk_effect)
table(effects_combined_strict$has_celltype_effect)

rbind(
  data.frame(
    Analysis = "bulk", 
    table(effects_combined_strict$has_bulk_effect)
  ), 
  data.frame(
    Analysis = "celltype", 
    table(effects_combined_strict$has_celltype_effect)
  )
) %>%
  ggplot(aes(x = Analysis, y = Freq, fill = Var1)) + geom_bar(stat = "identity") + 
  theme_paper(textsize = 50)

ggsave("~/Desktop/")

# by this analysis, we see more effects -- but this has a multiple testing problem as oli mentioned

# do the same thing only on strong effects
# get bayes factors (-LL_alt - (-LL_hyp))

fits_bulk$bf_cis <- fits_bulk$LL.conserved - fits_bulk$LL.cis
fits_bulk$bf_trans <- fits_bulk$LL.conserved - fits_bulk$LL.trans
fits_bulk$bf_cis_trans <- fits_bulk$LL.conserved - fits_bulk$LL.cis_trans
fits_bulk$combined_bf <- apply(fits_bulk[,c("bf_cis", "bf_trans", "bf_cis_trans")], 1, function(x){max(x)})

fits_sc$bf_cis <- fits_sc$LL.conserved - fits_sc$LL.cis
fits_sc$bf_trans <- fits_sc$LL.conserved - fits_sc$LL.trans
fits_sc$bf_cis_trans <- fits_sc$LL.conserved - fits_sc$LL.cis_trans
fits_sc$combined_bf <- apply(fits_sc[,c("bf_cis", "bf_trans", "bf_cis_trans")], 1, function(x){max(x)})

fits_rs$bf_cis <- fits_rs$LL.conserved - fits_rs$LL.cis
fits_rs$bf_trans <- fits_rs$LL.conserved - fits_rs$LL.trans
fits_rs$bf_cis_trans <- fits_rs$LL.conserved - fits_rs$LL.cis_trans
fits_rs$combined_bf <- apply(fits_rs[,c("bf_cis", "bf_trans", "bf_cis_trans")], 1, function(x){max(x)})

fits_es$bf_cis <- fits_es$LL.conserved - fits_es$LL.cis
fits_es$bf_trans <- fits_es$LL.conserved - fits_es$LL.trans
fits_es$bf_cis_trans <- fits_es$LL.conserved - fits_es$LL.cis_trans
fits_es$combined_bf <- apply(fits_es[,c("bf_cis", "bf_trans", "bf_cis_trans")], 1, function(x){max(x)})

bf_cutoff <- 10
effects_combined_strict <- effects_combined[,1:5]

effects_combined_strict[fits_bulk$combined_bf < bf_cutoff, ]$cat_bulk <- "conserved"
effects_combined_strict[fits_sc$combined_bf < bf_cutoff, ]$cat_sc <- "conserved"
effects_combined_strict[fits_rs$combined_bf < bf_cutoff, ]$cat_rs <- "conserved"
effects_combined_strict[fits_es$combined_bf < bf_cutoff, ]$cat_es <- "conserved"

effects_combined_strict$has_bulk_effect <- effects_combined_strict$cat_bulk != "conserved"
effects_combined_strict$has_celltype_effect <- apply(effects_combined_strict[,c("cat_sc", "cat_rs", "cat_es")], 1, function(x){any(x != "conserved")})

table(effects_combined_strict$has_bulk_effect)
table(effects_combined_strict$has_celltype_effect)

rbind(
  data.frame(
    Analysis = "bulk", 
    table(effects_combined_strict$has_bulk_effect)
  ), 
  data.frame(
    Analysis = "celltype", 
    table(effects_combined_strict$has_celltype_effect)
  )
) %>%
  ggplot(aes(x = Analysis, y = Freq, fill = Var1)) + geom_bar(stat = "identity") + 
  theme_paper(textsize = 40) + theme(legend.position = "None") + 
  annotate(geom = "text", x = c(1, 2), y = 1000, label = "with genetic effect", size = 10) +
  annotate(geom = "text", x = c(1, 2), y = 4000, label = "without genetic effect", size = 10) + 
  xlab("") + ylab("Number of genes")

ggsave("~/Desktop/Projects/ASE_Spermatogenesis_Paper/Plots/FigureS1/FigS1_11.pdf")

effect_counts_bulk <- table(effects_combined_strict$cat_bulk)
effect_counts_sc <- table(effects_combined_strict$cat_sc)
effect_counts_rs <- table(effects_combined_strict$cat_rs)
effect_counts_es <- table(effects_combined_strict$cat_es)

sum(effect_counts_bulk) - effect_counts_bulk[["conserved"]]
sum(effect_counts_sc) - effect_counts_sc[["conserved"]]
sum(effect_counts_rs) - effect_counts_rs[["conserved"]]
sum(effect_counts_es) - effect_counts_es[["conserved"]]

ggplot(rbind(
  data.frame(value = fits_bulk$combined_bf, comp = "bulk"), 
  data.frame(value = fits_sc$combined_bf, comp = "sc"), 
  data.frame(value = fits_rs$combined_bf, comp = "rs"), 
  data.frame(value = fits_es$combined_bf, comp = "es")
), aes(x = value, fill = comp)) + geom_histogram(bins = 100) + facet_wrap(~comp) + geom_vline(xintercept = bf_cutoff) + 
  theme_classic()

plot_df <- rbind(
  data.frame(value = fits_bulk$combined_bf, comp = "bulk"), 
  data.frame(value = fits_sc$combined_bf, comp = "sc"), 
  data.frame(value = fits_rs$combined_bf, comp = "rs"), 
  data.frame(value = fits_es$combined_bf, comp = "es")
)

plot_df[plot_df$value > 100, ]$value <- 100

ggplot(plot_df, aes(x = value, fill = comp)) + geom_histogram(bins = 100) + facet_wrap(~comp) + geom_vline(xintercept = bf_cutoff) + 
  theme_classic() + theme(text = element_text(size = 20))

ggplot(fits_bulk, aes(bf_cis, abs(FC_filial), col = Category)) + geom_point() + geom_vline(xintercept = 10) + theme_classic() + 
  theme(text = element_text(size = 20))

# finally look at changes between the bayes factors for individual models

full_df <- rbind(
  cbind(fits_sc, CellType = "SC"), 
  cbind(fits_rs, CellType = "RS"), 
  cbind(fits_es, CellType = "ES")
)
full_df <- full_df[,c("Gene", "Category", "bf_cis", "bf_trans", "bf_cis_trans", "CellType")]

library(tidyverse)
library(dplyr)

full_df %>% dplyr::select(Gene, Category, bf_cis, CellType) %>%
  pivot_wider(names_from = CellType, values_from = bf_cis) %>%
  ggplot(aes(SC, RS)) + geom_point() + geom_abline() + 
  theme_classic() + geom_vline(xintercept = 10, linetype = "dashed") + geom_hline(yintercept = 10, linetype = "dashed") + 
  xlim(0, 100) + ylim(0, 100) + ggtitle("Comparison of cis effects") + 
  theme_paper(textsize = 40) + 
  xlab("Likelihood ratio (Spermatocytes)") + 
  ylab("Likelihood ratio (Round Spermatids)")
ggsave("~/Desktop/Projects/ASE_Spermatogenesis_Paper/Plots/FigureS1/FigS1_11a.pdf")

full_df %>% dplyr::select(Gene, Category, bf_trans, CellType) %>%
  pivot_wider(names_from = CellType, values_from = bf_trans) %>%
  ggplot(aes(SC, RS)) + geom_point() + geom_abline() + 
  theme_classic() + geom_vline(xintercept = 10, linetype = "dashed") + geom_hline(yintercept = 10, linetype = "dashed") + 
  xlim(0, 100) + ylim(0, 100) + ggtitle("Comparison of trans effects") + 
  theme_paper(textsize = 40) + 
  xlab("Likelihood ratio (Spermatocytes)") + 
  ylab("Likelihood ratio (Round Spermatids)")
ggsave("~/Desktop/Projects/ASE_Spermatogenesis_Paper/Plots/FigureS1/FigS1_11b.pdf")


