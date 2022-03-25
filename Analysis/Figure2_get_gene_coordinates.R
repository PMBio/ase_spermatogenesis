library("GenomicFeatures")
library("TxDb.Mmusculus.UCSC.mm10.knownGene")
genome <- TxDb.Mmusculus.UCSC.mm10.knownGene

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

all_coordinates <- list(
  exons_full, 
  introns_full, 
  five_utrs_full, 
  three_utrs_full, 
  promoters_full
)

all_coordinates_collapsed <- list(
  exons_full_collapsed, 
  introns_full_collapsed, 
  five_utrs_full_collapsed, 
  three_utrs_full_collapsed, 
  promoters_full_collapsed
)

saveRDS(all_coordinates, "~/Desktop/all_gene_coordinates.rds")
saveRDS(all_coordinates_collapsed, "~/Desktop/all_gene_coordinates_collapsed.rds")

# get mapping from genes to transcripts
library(org.Mm.eg.db)
# gene_symbols <- getSYMBOL(names(exons_full), data = "org.Mm.eg.db")

# exons_full <- exons_full[gene_symbols %in% rownames(genes_test)]

names(exons_full) <- unique(exons_full_collapsed$name)
names(introns_full) <- unique(introns_full_collapsed$name)
names(five_utrs_full) <- unique(five_utrs_full_collapsed$name)
names(three_utrs_full) <- unique(three_utrs_full_collapsed$name)

# transcripts to genes
library(biomaRt)
mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
res <- getBM(attributes = c('ensembl_transcript_id_version',
                            'ensembl_gene_id',
                            'external_transcript_name',
                            'external_gene_name'), mart = mart)
res <- res[res$external_gene_name %in% rownames(genes_test), ]
res$ensembl_transcript_id_version <- gsub("\\.[0-9]*$", "", res$ensembl_transcript_id_version)
introns_full <- introns_full[names(introns_full) %in% res$ensembl_transcript_id_version]
five_utrs_full <- five_utrs_full[names(five_utrs_full) %in% res$ensembl_transcript_id_version]
three_utrs_full <- three_utrs_full[names(three_utrs_full) %in% res$ensembl_transcript_id_version]
promoters_full <- promoters_full[names(promoters_full) %in% res$ensembl_transcript_id_version]

library(EnsDb.Mmusculus.v79)
exons_full_collapsed$gene_name <- mapIds(org.Mm.eg.db, exons_full_collapsed$name, 'SYMBOL', 'ENTREZID')
introns_full_collapsed$gene_name <- mapIds(org.Mm.eg.db, introns_full_collapsed$name, 'SYMBOL', 'ENSEMBLTRANS')
five_utrs_full_collapsed$gene_name <- mapIds(org.Mm.eg.db, five_utrs_full_collapsed$name, 'SYMBOL', 'ENSEMBLTRANS')
three_utrs_full_collapsed$gene_name <- mapIds(org.Mm.eg.db, three_utrs_full_collapsed$name, 'SYMBOL', 'ENSEMBLTRANS')
promoters_full_collapsed$gene_name <- mapIds(org.Mm.eg.db, names(promoters_full_collapsed), 'SYMBOL', 'ENSEMBLTRANS')

all_coordinates_annotated <- list(
  exons_full_collapsed, 
  introns_full_collapsed, 
  five_utrs_full_collapsed, 
  three_utrs_full_collapsed, 
  promoters_full_collapsed
)

saveRDS(all_coordinates_annotated, "./Data/processed/all_gene_coordinates_annotated.rds")
