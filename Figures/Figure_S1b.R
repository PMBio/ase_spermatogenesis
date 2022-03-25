# analysis of goncalves fitting procedure

theme_paper <- function(){
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1), 
        text = element_text(size = 20))
}


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
  theme_paper()

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
  #replace_na(list(SC = 1, RS = 1, ES = 1)) %>%
  ggplot(aes(SC, RS)) + geom_point() + geom_abline() + 
  theme_classic() + geom_vline(xintercept = 10, linetype = "dashed") + geom_hline(yintercept = 10, linetype = "dashed") + 
  xlim(0, 100) + ylim(0, 100) + ggtitle("Comparison of cis effects")

full_df %>% dplyr::select(Gene, Category, bf_trans, CellType) %>%
  pivot_wider(names_from = CellType, values_from = bf_trans) %>%
  #replace_na(list(SC = 1, RS = 1, ES = 1)) %>%
  ggplot(aes(SC, RS)) + geom_point() + geom_abline() + 
  theme_classic() + geom_vline(xintercept = 10, linetype = "dashed") + geom_hline(yintercept = 10, linetype = "dashed") + 
  xlim(0, 100) + ylim(0, 100) + ggtitle("Comparison of trans effects")

full_df %>% dplyr::select(Gene, Category, bf_cis_trans, CellType) %>%
  pivot_wider(names_from = CellType, values_from = bf_cis_trans) %>%
  #replace_na(list(SC = 1, RS = 1, ES = 1)) %>%
  ggplot(aes(SC, RS)) + geom_point() + geom_abline() + 
  theme_classic() + geom_vline(xintercept = 10, linetype = "dashed") + geom_hline(yintercept = 10, linetype = "dashed") + 
  xlim(0, 100) + ylim(0, 100)


