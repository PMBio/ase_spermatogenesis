library("GenomicFeatures")
library("TxDb.Mmusculus.UCSC.mm10.knownGene")
genome <- TxDb.Mmusculus.UCSC.mm10.knownGene

genes_test <- readRDS("~/Desktop/genes_test.rds")

 # get exon coordinates

exons_full <- exonsBy(TxDb.Mmusculus.UCSC.mm10.knownGene, by = c("gene"))
names(exons_full) <- gsub("\\.[0-9]*$", "", names(exons_full))
exons_full <- exons_full[lapply(exons_full, length) > 0]
exons_full <- lapply(1:length(exons_full), function(i){
  x = exons_full[[i]]
  name = names(exons_full)[[i]]
  x$name = name
  x
})
exons_full_collapsed <- do.call("c", exons_full)

introns_full <- intronsByTranscript(TxDb.Mmusculus.UCSC.mm10.knownGene, use.names = T)
names(introns_full) <- gsub("\\.[0-9]*$", "", names(introns_full))
introns_full <- introns_full[lapply(introns_full, length) > 0]
introns_full <- lapply(1:length(introns_full), function(i){
  x = introns_full[[i]]
  name = names(introns_full)[[i]]
  x$name = name
  x
})
introns_full_collapsed <- do.call("c", introns_full)

five_utrs_full <- fiveUTRsByTranscript(TxDb.Mmusculus.UCSC.mm10.knownGene, use.names = T)
names(five_utrs_full) <- gsub("\\.[0-9]*$", "", names(five_utrs_full))
five_utrs_full <- five_utrs_full[lapply(five_utrs_full, length) > 0]
five_utrs_full <- lapply(1:length(five_utrs_full), function(i){
  x = five_utrs_full[[i]]
  name = names(five_utrs_full)[[i]]
  x$name = name
  x
})
five_utrs_full_collapsed <- do.call("c", five_utrs_full)

three_utrs_full <- threeUTRsByTranscript(TxDb.Mmusculus.UCSC.mm10.knownGene, use.names = T)
names(three_utrs_full) <- gsub("\\.[0-9]*$", "", names(three_utrs_full))
three_utrs_full <- three_utrs_full[lapply(three_utrs_full, length) > 0]
three_utrs_full <- lapply(1:length(three_utrs_full), function(i){
  x = three_utrs_full[[i]]
  name = names(three_utrs_full)[[i]]
  x$name = name
  x
})
three_utrs_full_collapsed <- do.call("c", three_utrs_full)

promoters_full <- promoters(TxDb.Mmusculus.UCSC.mm10.knownGene, use.names = T)
names(promoters_full) <- gsub("\\.[0-9]*$", "", names(promoters_full))
promoters_full_collapsed <- promoters_full

# get mapping from genes to transcripts
library(org.Mm.eg.db)
library(annotate)
names(exons_full) <- unlist(lapply(exons_full, function(x){x$name[[1]]}))
gene_symbols <- getSYMBOL(names(exons_full), data = "org.Mm.eg.db")

exons_full <- exons_full[gene_symbols %in% rownames(genes_test)]

# transcripts to genes
library(biomaRt)
mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
res <- getBM(attributes = c('ensembl_transcript_id_version',
                            'ensembl_gene_id',
                            'external_transcript_name',
                            'external_gene_name'), mart = mart)
res <- res[res$external_gene_name %in% rownames(genes_test), ]
res$ensembl_transcript_id_version <- gsub("\\.[0-9]*$", "", res$ensembl_transcript_id_version)
rownames(res) <- res$ensembl_transcript_id_version
introns_full <- introns_full[names(introns_full) %in% res$ensembl_transcript_id_version]
five_utrs_full <- five_utrs_full[names(five_utrs_full) %in% res$ensembl_transcript_id_version]
three_utrs_full <- three_utrs_full[names(three_utrs_full) %in% res$ensembl_transcript_id_version]
promoters_full <- promoters_full[names(promoters_full) %in% res$ensembl_transcript_id_version]

# read snp files

read_snp_file <- function(x){
  snps.chr1 <- read.csv(x, skip = 1, sep = "\t", header = F)
  colnames(snps.chr1) <- c("ID", "Chr", "Start", "Strand", "Alleles", "Info")
  snps.chr1$Strand <- "+"
  snps.chr1$Chr <- paste0("chr", snps.chr1$Chr)
  snps.chr1 <- makeGRangesFromDataFrame(snps.chr1, seqnames.field = "Chr", start.field = "Start", end.field = "Start")
  snps.chr1
}

all_snps <- lapply(list.files("~/Desktop/Cast_SNPs/", full.names = T), read_snp_file)
all_snps_combined <- do.call("c", all_snps)

promoters_full$n_snps <- countOverlaps(promoters_full, all_snps_combined)
promoters_full$snp_rate <- promoters_full$n_snps / width(promoters_full)
promoters_full$gene_name <- res[gsub("\\..*$", "", promoters_full$tx_name), ]$external_gene_name
promoters_full_agg <- data.frame("rate" = promoters_full$snp_rate, "name" = promoters_full$gene_name)
promoters_full_agg <- data.frame(aggregate(promoters_full_agg, list(promoters_full_agg$name), mean))
rownames(promoters_full_agg) <- promoters_full_agg$Group.1

gene_names_save <- rownames(genes_test)
genes_test <- data.frame(apply(genes_test, 2, as.numeric))
rownames(genes_test) <- gene_names_save
genes_test$promoter_snps <- promoters_full_agg[rownames(genes_test), ]$rate

ggplot(genes_test, aes(BF_trans_dynamic > 10, promoter_snps)) + geom_boxplot() + ggpubr::stat_compare_means()
