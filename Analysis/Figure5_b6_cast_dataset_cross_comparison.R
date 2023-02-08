### 

library(scran)
library(scater)
library(tidyverse)

setwd("~/Desktop/Projects/ASE_Spermatogenesis_Paper_rerun/Revisions/")

source("../Scripts/General/auxiliary.R")
source("../Scripts/General/reuse_functions.R")

theme_paper <- function(textsize = 30){
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        text = element_text(size = textsize))
}

colour_vector <- c("SG" = "#046735",
                   "SC" = "#D31281",
                   "RS" = "#282B69",
                   "ES" = "grey60",
                   "Sertoli" = "#643F18",
                   "Leydig" = "#985D25",
                   "Immune" = "#FFDFB2")


data <- readRDS("./Data/processed/sce_merged_new.rds")
data$IndividualSamples <- paste0(data$Dataset, "_", data$Library)
data <- data[,data$Library != "Sample7"]

bulk_all <- aggregateAcrossCells(data, ids = data$Sample)
bulk_sc <- aggregateAcrossCells(data[,data$CellType == "SC"], ids = data[,data$CellType == "SC"]$Sample)
bulk_rs <- aggregateAcrossCells(data[,data$CellType == "RS"], ids = data[,data$CellType == "RS"]$Sample)
bulk_es <- aggregateAcrossCells(data[,data$CellType == "ES"], ids = data[,data$CellType == "ES"]$Sample)

data <- readRDS("~/Desktop/Projects/ASE_Spermatogenesis_Paper_rerun/Data/processed/final_sce_evo_dataset.rds")
data$IndividualSamples <- paste0(data$Dataset, "_", data$Library)
#data <- data[,data$Library != "Sample7"]

bulk_all_evo <- aggregateAcrossCells(data, ids = data$Sample)
bulk_sc_evo <- aggregateAcrossCells(data[,data$CellType %in% c("Pachytene", "Diplotene", "Meiosis")], ids = data[,data$CellType %in% c("Pachytene", "Diplotene", "Meiosis")]$Sample)
bulk_rs_evo <- aggregateAcrossCells(data[,grepl("Round Spermatids", data$CellType)], ids = data[,grepl("Round Spermatids", data$CellType)]$Sample)
bulk_es_evo <- aggregateAcrossCells(data[,grepl("Elongating Spermatids", data$CellType)], ids = data[,grepl("Elongating Spermatids", data$CellType)]$Sample)

# Get fold changes CAST / B6 dataset 1

joint_genes <- intersect(rownames(bulk_all), rownames(bulk_all_evo))

# data.frame(
#   b6_sc_ds1 = rowSums(counts(bulk_sc[joint_genes,"B6"])), 
#   b6_sc_ds2 = rowSums(counts(bulk_sc_evo[joint_genes,"B6"]))
# ) %>% ggplot(aes(b6_sc_ds1, b6_sc_ds2)) + geom_point() + scale_x_log10() + scale_y_log10() + geom_abline()

joint_data_sc <- data.frame(
  b6_sc_ds1 = rowSums(counts(bulk_sc[joint_genes,"B6"])), 
  b6_sc_ds2 = rowSums(counts(bulk_sc_evo[joint_genes,"B6"])), 
  cast_sc_ds1 = rowSums(counts(bulk_sc[joint_genes,"CAST"])), 
  cast_sc_ds2 = rowSums(counts(bulk_sc_evo[joint_genes,"CAST"]))
)

joint_data_sc <- joint_data_sc[rowMeans(joint_data_sc) > 50, ]
 
joint_data_sc <- data.frame(t(t(joint_data_sc) / colSums(joint_data_sc)))

joint_data_sc <- joint_data_sc %>%
  mutate(logfc_ds1 = log2(b6_sc_ds1 / cast_sc_ds1)) %>%
  mutate(logfc_ds2 = log2(b6_sc_ds2 / cast_sc_ds2)) %>%
  add_column(chromosome = rowData(bulk_sc)[rownames(joint_data_sc), ]$chromosome_name)
  
joint_data_sc %>%
  dplyr::filter(!(is.na(logfc_ds1) | is.na(logfc_ds2))) %>% 
  dplyr::filter(!(is.infinite(logfc_ds1) | is.infinite(logfc_ds2))) %>% 
  dplyr::select(c("logfc_ds1", "logfc_ds2")) %>% cor() -> cors_here 

joint_data_sc %>% ggplot(aes(logfc_ds1, logfc_ds2)) + geom_point(size = 0.5, col = colour_vector[["SC"]]) + geom_abline() + theme_paper(textsize = 30) +
  xlab("log2(FC) B6 / Cast (Spermatids, Dataset1)") + ylab("log2(FC) B6 / Cast (Spermatids, Dataset2)") + 
  annotate(x = 8, y = -8, geom = "text", label = paste0("Pearson correlation: \n", round(cors_here[1, 2], digits = 4)), size = 10)
ggsave("~/Desktop/revision_plots/CrossComparisonDatasets_SC_FCs.pdf")



joint_data_rs <- data.frame(
  b6_sc_ds1 = rowSums(counts(bulk_rs[joint_genes,"B6"])), 
  b6_sc_ds2 = rowSums(counts(bulk_rs_evo[joint_genes,"B6"])), 
  cast_sc_ds1 = rowSums(counts(bulk_rs[joint_genes,"CAST"])), 
  cast_sc_ds2 = rowSums(counts(bulk_rs_evo[joint_genes,"CAST"]))
)

joint_data_rs <- joint_data_rs[rowMeans(joint_data_rs) > 50, ]

joint_data_rs <- data.frame(t(t(joint_data_rs) / colSums(joint_data_rs)))

joint_data_rs <- joint_data_rs %>%
  mutate(logfc_ds1 = log2(b6_sc_ds1 / cast_sc_ds1)) %>%
  mutate(logfc_ds2 = log2(b6_sc_ds2 / cast_sc_ds2)) %>%
  add_column(chromosome = rowData(bulk_sc)[rownames(joint_data_rs), ]$chromosome_name)

joint_data_rs %>%
  dplyr::filter(!(is.na(logfc_ds1) | is.na(logfc_ds2))) %>% 
  dplyr::filter(!(is.infinite(logfc_ds1) | is.infinite(logfc_ds2))) %>% 
  dplyr::select(c("logfc_ds1", "logfc_ds2")) %>% cor() -> cors_here

joint_data_rs %>% ggplot(aes(logfc_ds1, logfc_ds2)) + geom_point(size = 0.5, col = colour_vector[["RS"]]) + geom_abline() + theme_paper(textsize = 30) +
  xlab("log2(FC) B6 / Cast (Spermatids, Dataset1)") + ylab("log2(FC) B6 / Cast (Spermatids, Dataset2)") + 
  annotate(x = 8, y = -8, geom = "text", label = paste0("Pearson correlation: \n", round(cors_here[1, 2], digits = 4)), size = 10)
ggsave("~/Desktop/revision_plots/CrossComparisonDatasets_RS_FCs.pdf")

# 
genes_joint_here <- intersect(rownames(joint_data_sc), rownames(joint_data_rs))
colnames(joint_data_sc) <- paste0(colnames(joint_data_sc), "_SC")
colnames(joint_data_rs) <- paste0(colnames(joint_data_rs), "_RS")
joint_data_joint <- cbind(joint_data_sc[genes_joint_here, ], joint_data_rs[genes_joint_here, ])

#cors_here <- cor(testytest$logfc_ds1_SC - testytest$logfc_ds1_RS, testytest$logfc_ds2_SC - testytest$logfc_ds2_RS)

joint_data_joint %>%
  ggplot(aes(logfc_ds1_SC - logfc_ds1_RS, logfc_ds2_SC - logfc_ds2_RS)) + geom_point(size = 0.1) + 
  geom_point(size = 0.5, col = "black") + geom_abline() + theme_paper(textsize = 30) +
  xlab("Fold change difference SC - RS (dataset_1)") + ylab("Fold change difference SC - RS (dataset2)") + 
  annotate(x = 5, y = -5, geom = "text", label = paste0("Pearson correlation: \n", round(cors_here, digits = 4)), size = 10)
ggsave("~/Desktop/revision_plots/CrossComparisonDatasets_SC_vs_RS_FCs.pdf")



