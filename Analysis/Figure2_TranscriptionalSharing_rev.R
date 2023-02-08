### Re comment: 
# Need to rule out that transcriptional sharing induces bias

# Transc sharing will not change the allelic ratio in spermatids, but might impact the variance
# We obtain list of genes which have been shown to be subject to transc sharing

# 1) show that these genes have higher over-dispersion in RS in our data
# 2) show increased correlation between neighbouring genes for unshared transcripts
# 3) show that there is no co-loc between genes with sharing and dynamic / non-dynamic genetic effects

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

table_path <- "~/Desktop/Projects/ASE_Spermatogenesis_Paper_rerun_Revisions/Revisions/abb1723-Bhutani-SM-Tables-S1-to-S11-Zipped/table_s1_mouse_genoinformativity.csv"
genoinformative <- read.csv(table_path, skip = 1, sep = ",")

ggplot(genoinformative, aes(GIM_status, Posterior_mean_genoinformativity)) + geom_boxplot()

genes_test <- readRDS("~/Desktop/Projects/ASE_Spermatogenesis_Paper_rerun_Revisions/Data/processed/genes_test.rds")

library(EnsDb.Mmusculus.v79)
geneIDs1 <- ensembldb::select(EnsDb.Mmusculus.v79, keys= genoinformative$Ensembl_gene_ID, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
genoinformative$GeneName <- geneIDs1$SYMBOL
genoinformative <- genoinformative[!duplicated(genoinformative$GeneName), ]
rownames(genoinformative) <- genoinformative$GeneName

genes_test$genoinformativity <- genoinformative[rownames(genes_test), ]$Posterior_mean_genoinformativity
genes_test[!is.na(genes_test$genoinformativity), ]
genes_test$GIM_status <- genoinformative[rownames(genes_test), ]$GIM_status
genes_test <- genes_test[!is.na(genes_test$GIM_status), ]

ggplot(genes_test, aes(GIM_status, abs(Mean - 0.5))) + geom_boxplot() + theme_paper()
ggplot(genes_test, aes(GIM_status, qdiff_cis)) + geom_boxplot() + theme_paper() + scale_y_sqrt()
ggplot(genes_test, aes(GIM_status, -log10(pval))) + geom_boxplot() + theme_paper()
ggplot(genes_test, aes(GIM_status, -log10(dali_pval_linear))) + geom_boxplot() + theme_paper() + scale_y_sqrt()

rename_columns <- list(
  "Mean" = "(1) Absolute AI - 0.5", 
  "qdiff_cis" = "(2) Interquartile range of AI", 
  "pval" = "(3) scDALI persistent p-value", 
  "dali_pval_linear" = "(4) scDALI dynamic p-value"
)

genes_test %>%
  dplyr::select(c(GIM_status, Mean, qdiff_cis, pval, dali_pval_linear)) %>%
  dplyr::rename() %>%
  mutate(Mean = abs(Mean - 0.5)) %>%
  mutate(pval = -log10(pval)) %>%
  mutate(dali_pval_linear = -log10(dali_pval_linear)) %>%
  pivot_longer(-c(GIM_status)) %>%
  mutate(nice_names = unlist(rename_columns[name])) %>%
  ggplot(aes(x = GIM_status, y = value)) + facet_wrap(~nice_names, scales = "free_x") + geom_boxplot(fill = "grey") + scale_y_sqrt() + theme_paper(35) + coord_flip() + 
    xlab("") + ylab("")
ggsave("~/Desktop/Projects/ASE_Spermatogenesis_Paper_rerun_Revisions/Plots/FigureS4_rev/genoinformative_genes_comparison.pdf")

fisher.test(table(abs(genes_test$Mean - 0.5) > 0.1, genes_test$GIM_status))
fisher.test(table(genes_test$qdiff_cis > 0.1, genes_test$GIM_status))

fisher.test(table(p.adjust(genes_test$pval) < 0.1, genes_test$GIM_status))
fisher.test(table(p.adjust(genes_test$dali_pval_polynomial) < 0.1, genes_test$GIM_status))



