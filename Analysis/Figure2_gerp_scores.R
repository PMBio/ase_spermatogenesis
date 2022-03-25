### gerp scores

library(GenomicScores)

# get gene coordinates

library("ensembldb")
library("EnsDb.Mmusculus.v79")
genome <- EnsDb.Mmusculus.v79

# get mapping from genes to transcripts
all_genes <- genes(EnsDb.Mmusculus.v79, columns = c("tx_name", "symbol"))
transcript_to_symbol <- data.frame(tx = all_genes$tx_name, symbol = all_genes$symbol) %>%
  column_to_rownames("tx")

# get exon coordinates

exons_full_collapsed <- ensembldb::exons(EnsDb.Mmusculus.v79, columns = c("gene_name", "symbol"))

introns_full <- intronsByTranscript(EnsDb.Mmusculus.v79, use.names = T)
introns_full_collapsed <- unlist(introns_full)
introns_full_collapsed$tx_name <- names(introns_full_collapsed)
introns_full_collapsed$symbol <- transcript_to_symbol[introns_full_collapsed$tx_name, ]

five_utrs_full <- fiveUTRsByTranscript(EnsDb.Mmusculus.v79)
five_utrs_full_collapsed <- unlist(five_utrs_full)
five_utrs_full_collapsed$tx_name <- names(five_utrs_full_collapsed)
five_utrs_full_collapsed$symbol <- transcript_to_symbol[five_utrs_full_collapsed$tx_name, ]

three_utrs_full <- threeUTRsByTranscript(EnsDb.Mmusculus.v79)
three_utrs_full_collapsed <- unlist(three_utrs_full)
three_utrs_full_collapsed$tx_name <- names(three_utrs_full_collapsed)
three_utrs_full_collapsed$symbol <- transcript_to_symbol[three_utrs_full_collapsed$tx_name, ]

promoters_full <- promoters(EnsDb.Mmusculus.v79, use.names = T, upstream = 2000, downstream = 200, columns = "symbol")
promoters_full_collapsed <- promoters_full

all_coordinates <- list(
  exons_full_collapsed, 
  introns_full_collapsed, 
  five_utrs_full_collapsed, 
  three_utrs_full_collapsed, 
  promoters_full_collapsed
)
names(all_coordinates) <- c("Exon", "Intron", "FiveUTR", "ThreeUTR", "Promoter")
saveRDS(all_coordinates, "~/Desktop/all_gene_coordinates_annotated.rds")

# cons_scores <- getGScores("phyloP60way.UCSC.mm10")
# scores_exons <- gscores(cons_scores, exons_full_collapsed)
# scores_introns <- gscores(cons_scores, introns_full_collapsed)
# scores_five_utr <- gscores(cons_scores, five_utrs_full_collapsed)
# scores_three_utr <- gscores(cons_scores, three_utrs_full_collapsed)
# scores_promoter <- gscores(cons_scores, promoters_full_collapsed[1:(length(promoters_full_collapsed) - 1000)])
# # this crashes for one or more entries at the end, so i just cut off the last 1000 entries -- 
# # this only excludes Y-linked, chrM and unplaced scaffolds that we don't look at anyways (we can exclude in the beginning)
# 
# scores_exons_new <- scores_exons
# scores_introns_new <- scores_introns
# scores_five_utr_new <- scores_five_utr
# scores_three_utr_new <- scores_three_utr
# #scores_exons_new <- scores_exons
# 
# ### assign genes to transcript and compute average, min and max scores
# # exons
# 
# library(EnsDb.Mmusculus.v79)
# 
# scores_exons_new$gene_symbol <- mapIds(org.Mm.eg.db, scores_exons_new$name, 'SYMBOL', 'ENTREZID')
# scores_exons_new <- scores_exons_new[scores_exons_new$gene_symbol %in% rownames(genes_test), ]
# exon_df_results <- data.frame(
#   gene_symbols = scores_exons_new$gene_symbol, 
#   score = scores_exons_new$default
# )
# exon_df_results <- exon_df_results[!is.na(exon_df_results$score), ]
# exon_df_results <- aggregate(exon_df_results, list(exon_df_results$gene_symbols), function(x){c(mean(x), min(x), max(x))})
# 
# # introns
# txmap <- ensembldb::select(EnsDb.Mmusculus.v79, keys=scores_introns_new$name, keytype = "TXID", columns = c("SYMBOL","TXID"))
# rownames(txmap) <- txmap$TXID
# txmap <- txmap[txmap$SYMBOL %in% rownames(genes_test), ]
# scores_introns_new <- scores_introns_new[scores_introns_new$name %in% txmap$TXID, ]
# scores_introns_new$gene_symbol <- txmap[scores_introns_new$name, ]$SYMBOL
# intron_df_results <- data.frame(
#   gene_symbols = scores_introns_new$gene_symbol, 
#   score = scores_introns_new$default
# )
# intron_df_results <- intron_df_results[!is.na(intron_df_results$score), ]
# intron_df_results <- aggregate(intron_df_results, list(intron_df_results$gene_symbols), 
#                                function(x){c(mean(x), min(x), max(x))})
# 
# # 5 utrs
# txmap <- ensembldb::select(EnsDb.Mmusculus.v79, keys=scores_five_utr$name, keytype = "TXID", columns = c("SYMBOL","TXID"))
# rownames(txmap) <- txmap$TXID
# txmap <- txmap[txmap$SYMBOL %in% rownames(genes_test), ]
# scores_five_utr <- scores_five_utr[scores_five_utr$name %in% txmap$TXID, ]
# scores_five_utr$gene_symbol <- txmap[scores_five_utr$name, ]$SYMBOL
# five_utr_df_results <- data.frame(
#   gene_symbols = scores_five_utr$gene_symbol, 
#   score = scores_five_utr$default
# )
# five_utr_df_results <- five_utr_df_results[!is.na(five_utr_df_results$score), ]
# five_utr_df_results <- aggregate(five_utr_df_results, list(five_utr_df_results$gene_symbols), 
#                                function(x){c(mean(x), min(x), max(x))})
# 
# # promoters
# scores_promoter$name <- gsub(".[0-9]*$", "", scores_promoter$tx_name)
# txmap <- ensembldb::select(EnsDb.Mmusculus.v79, keys=scores_promoter$name, keytype = "TXID", columns = c("SYMBOL","TXID"))
# rownames(txmap) <- txmap$TXID
# txmap <- txmap[txmap$SYMBOL %in% rownames(genes_test), ]
# scores_promoter <- scores_promoter[scores_promoter$name %in% txmap$TXID, ]
# scores_promoter$gene_symbol <- txmap[scores_promoter$name, ]$SYMBOL
# promoter_df_results <- data.frame(
#   gene_symbols = scores_promoter$gene_symbol, 
#   score = scores_promoter$default
# )
# promoter_df_results <- promoter_df_results[!is.na(promoter_df_results$score), ]
# promoter_df_results <- aggregate(promoter_df_results, list(promoter_df_results$gene_symbols), 
#                                  function(x){c(mean(x), min(x), max(x))})
# 
# # three utrs
# txmap <- ensembldb::select(EnsDb.Mmusculus.v79, keys=scores_three_utr$name, keytype = "TXID", columns = c("SYMBOL","TXID"))
# rownames(txmap) <- txmap$TXID
# txmap <- txmap[txmap$SYMBOL %in% rownames(genes_test), ]
# scores_three_utr <- scores_three_utr[scores_three_utr$name %in% txmap$TXID, ]
# scores_three_utr$gene_symbol <- txmap[scores_three_utr$name, ]$SYMBOL
# three_utr_df_results <- data.frame(
#   gene_symbols = scores_three_utr$gene_symbol, 
#   score = scores_three_utr$default
# )
# three_utr_df_results <- three_utr_df_results[!is.na(three_utr_df_results$score), ]
# three_utr_df_results <- aggregate(three_utr_df_results, list(three_utr_df_results$gene_symbols), 
#                                   function(x){c(mean(x), min(x), max(x))})
# 
# ### Summary 
# 
# rownames(exon_df_results) <- exon_df_results$Group.1
# rownames(intron_df_results) <- intron_df_results$Group.1
# rownames(five_utr_df_results) <- five_utr_df_results$Group.1
# rownames(three_utr_df_results) <- three_utr_df_results$Group.1
# rownames(promoter_df_results) <- promoter_df_results$Group.1

#saveRDS(list(exon_df_results, intron_df_results, five_utr_df_results, three_utr_df_results, promoter_df_results), 
#        "~/Desktop/phastcons_scores.rds")
# phast_cons_scores <- readRDS("~/Desktop/Projects/ASE_Spermatogenesis_Paper/Data/processed/phastcons_scores.rds")
