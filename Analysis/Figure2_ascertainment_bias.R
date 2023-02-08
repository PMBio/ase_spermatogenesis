setwd("~/Desktop/Projects/ASE_Spermatogenesis_Paper_rerun_Revisions/")

data_exp1_sample5 <- read_delim("~/Desktop/Projects/ASE_Spermatogenesis_Paper_rerun_Revisions/Revisions/ascertainment_bias_data/Exp1_Sample5_per_snp_coverage.tsv", delim = "\t", 
                              col_names = c("Gene", "Position", "REF", "ALT"))
data_exp1_sample6 <- read_delim("~/Desktop/Projects/ASE_Spermatogenesis_Paper_rerun_Revisions/Revisions/ascertainment_bias_data/Exp1_Sample6_per_snp_coverage.tsv", delim = "\t", 
                                col_names = c("Gene", "Position", "REF", "ALT"))

data_exp2_sample5 <- read_delim("~/Desktop/Projects/ASE_Spermatogenesis_Paper_rerun_Revisions/Revisions/ascertainment_bias_data/Exp2_Sample5_per_snp_coverage.tsv", delim = "\t", 
                                col_names = c("Gene", "Position", "REF", "ALT"))
data_exp2_sample6 <- read_delim("~/Desktop/Projects/ASE_Spermatogenesis_Paper_rerun_Revisions/Revisions/ascertainment_bias_data/Exp2_Sample6_per_snp_coverage.tsv", delim = "\t", 
                                col_names = c("Gene", "Position", "REF", "ALT"))

data_exp3_sample5 <- read_delim("~/Desktop/Projects/ASE_Spermatogenesis_Paper_rerun_Revisions/Revisions/ascertainment_bias_data/Exp3_Sample5_per_snp_coverage.tsv", delim = "\t", 
                                col_names = c("Gene", "Position", "REF", "ALT"))
data_exp3_sample6 <- read_delim("~/Desktop/Projects/ASE_Spermatogenesis_Paper_rerun_Revisions/Revisions/ascertainment_bias_data/Exp3_Sample6_per_snp_coverage.tsv", delim = "\t", 
                                col_names = c("Gene", "Position", "REF", "ALT"))

merge_files <- function(x, y){
  merged_df = merge(x, y, by.x = c("Gene", "Position"), by.y = c("Gene", "Position"), all = T)
  merged_df[is.na(merged_df)] <- 0
  merged_df = data.frame("Gene" = merged_df$Gene, "Position" = merged_df$Position, "REF" = merged_df$REF.x + merged_df$REF.y, "ALT" = merged_df$ALT.x + merged_df$ALT.y)
  merged_df
}

data_exp1 <- merge_files(data_exp1_sample5, data_exp1_sample6)
data_exp2 <- merge_files(data_exp2_sample5, data_exp2_sample6)
data_exp3 <- merge_files(data_exp3_sample5, data_exp3_sample6)

data_full_1 <- merge_files(data_exp1, data_exp2)
data <-  merge_files(data_full_1, data_exp3)

data.sce <- rowData(readRDS("~/Desktop/Projects/ASE_Spermatogenesis_Paper_rerun_Revisions/Data/processed/final_sce_f1_dataset.rds"))

genetic_effects_all <- readRDS("~/Desktop/Projects/ASE_Spermatogenesis_Paper_rerun_Revisions/Data/processed/genes_test.rds")

genetic_effects <- genetic_effects_all %>%
  mutate(static_effect = padj < 0.01) %>%
  mutate(dynamic_effect = p.adjust(dali_pval_polynomial) < 0.01) %>%
  dplyr::select(c(static_effect, dynamic_effect))

genetic_effects$ensembl <- data.sce[rownames(genetic_effects), ]$ID
genetic_effects$gene <- rownames(genetic_effects)

genetic_effects <- merge(genetic_effects, data, by.x = "ensembl", by.y = "Gene")
  
genetic_effects %>% mutate(Total = REF + ALT) %>%
  dplyr::filter(Total > 200) %>%
  group_by(ensembl) %>%
  summarize(n_snps = n(), static_effect = unique(static_effect), dynamic_effect = unique(dynamic_effect), gene = unique(gene)) -> test

ltm::biserial.cor(log2(test$n_snps), test$static_effect, level = 2)
ltm::biserial.cor(log2(test$n_snps), test$dynamic_effect, level = 2)

n_buckets = 10
log_base = 2
test %>%
  dplyr::filter(n_snps > 0) %>%
  mutate(bucket = cut(n_snps, breaks = c(0, 2^(1:n_buckets)), labels = c(paste0('>= ',scales::comma(c(2^(1:n_buckets- 1))))))) %>%
  ggplot(aes(x = bucket, fill = static_effect)) + geom_bar(position = "fill") + 
    theme_paper(textsize = 40) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + xlab("") + scale_fill_manual(values = c("grey", "black")) + 
    labs(fill = "Gene has persistent \n cis-effect") + ylab("Proportion") + scale_y_continuous(expand = c(0, 0))
ggsave("./Plots/FigureS9_rev/static_effect_proportions.pdf")

table(test$n_snps >= 2) / sum(table(test$n_snps >= 2))
table(test$n_snps >= 5) / sum(table(test$n_snps >= 5))

test %>%
  dplyr::filter(n_snps > 0) %>%
  mutate(bucket = cut(n_snps, breaks = c(0, 2^(1:n_buckets)), labels = c(paste0('>= ',scales::comma(c(2^(1:n_buckets- 1))))))) %>%
  ggplot(aes(x = bucket, fill = dynamic_effect)) + geom_bar(position = "fill") + 
    theme_paper(textsize = 40) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + xlab("") + scale_fill_manual(values = c("grey", "black")) + 
    labs(fill = "Gene has dynamic \n cis-effect") + ylab("Proportion") + scale_y_continuous(expand = c(0, 0))
ggsave("./Plots/FigureS9_rev/dynamic_effect_proportions.pdf")

## check if genes with different SNP densities have different scores between SNP-loadings

library(ggplot2)

features_prepared <- data.frame(readRDS("./Data/processed/all_features.rds"))
exon_df_results <- setNames(features_prepared$VAR_evolution_phastcons_exon, rownames(features_prepared))

test$Score_exon <- exon_df_results[test$gene]

test %>%
  mutate(bucket = cut(n_snps, breaks = c(0, 2^(1:n_buckets)), labels = c(paste0('>= ',scales::comma(c(2^(1:n_buckets- 1))))))) %>%
  mutate(conserved = !static_effect & !dynamic_effect) %>%
  mutate(static = static_effect & !dynamic_effect) %>%
  mutate(dynamic = dynamic_effect) %>%
  dplyr::select(c(n_snps, Score_exon, bucket, conserved, static, dynamic)) %>%
  pivot_longer(-c(n_snps, Score_exon, bucket)) %>%
  dplyr::filter(value) %>%
  mutate(name = factor(name, levels = c("conserved", "static", "dynamic"))) %>%
  ggplot(aes(x = name, y = Score_exon, fill = bucket)) + theme_paper(textsize = 25) + 
    facet_wrap(~bucket) + geom_boxplot(fill = "grey") + ylab("Exonic conservation (PhastCons)") + xlab("")
ggsave("./Plots/FigureS9_rev/phastcons_across_levels.pdf")
                                                                                    