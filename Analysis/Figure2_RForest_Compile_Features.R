# This script collects the features for the random forest analysis in Figure 2

library(GenomeInfoDb)
setwd("~/Desktop/Projects/ASE_Spermatogenesis_Paper_rerun_Revisions/")

genes_test <- readRDS("./Data/processed/genes_test.rds")
genes_test_vector <- rownames(genes_test)

# get coordinates of exons, introns, promoters, utrs 
gene_coordinates <- readRDS("./Data/processed/all_gene_coordinates_annotated.rds")

convert_granges <- function(x){
  seqlevelsStyle(x) <- "UCSC"
  x <- x[seqnames(x) %in% c(paste0("chr", 1:19)),  ]
  seqlevels(x) <- c(paste0("chr", 1:19))
  x
}
gene_coordinates <- lapply(gene_coordinates, convert_granges)

# ---- ## ---- ## ---- ## ---- #  
# ---- # Expression based features
# ---- ## ---- ## ---- ## ---- #
# Read data

library(scran)
library(scater)
library(stringr)
library(tidyverse)

data <- readRDS("./Data/processed/sce_merged_new.rds")
data <- logNormCounts(data)

# Get total expression level
VAR_expression_level <- rowMeans(counts(data[genes_test_vector, ]))

# Get pseudobulk expression, compute coefficient of var to estimate celltype-dependence (var / mean ^ 2)
collapsed.expression <- lapply(unique(data$CellType), function(cl){
  sce.here <- data[genes_test_vector, data$CellType == cl]
  cl_lib_means <- rowMeans(logcounts(sce.here))
  cl_lib_means
})
collapsed.expression <- do.call("cbind", collapsed.expression)
VAR_expression_variability <- sqrt(rowVars(collapsed.expression)) / (VAR_expression_level)

# Get single-cell level variability as squared coefficient of variation (var / mean ^ 2)
# This is just an estimate here... 
collapsed.variability <- lapply(unique(data$CellType), function(cl){
  sce.here <- data[genes_test_vector, data$CellType == cl]
  cl_lib_means <- rowVars(as.matrix(logcounts(sce.here)))
  cl_lib_means
})
collapsed.variability <- do.call("cbind", collapsed.variability)
VAR_expression_variability_post_celltype <- sqrt(rowMeans(collapsed.variability)) / (VAR_expression_level)

get_loess_residuals <- function(x, y){
  # first set up a lowess fit:
  lfit <- lowess(x,y)
  # create a functional version of the lowess fit
  lfun <- approxfun(lfit)
  fitted <- lfun(x)
  resid <- y-fitted
  resid
}

VAR_expression_variability <- get_loess_residuals(VAR_expression_level, VAR_expression_variability)
VAR_expression_variability_post_celltype <- get_loess_residuals(VAR_expression_level, VAR_expression_variability_post_celltype)

# house keeping genes
convertHumanGeneList <- function(x){
  library("biomaRt")
  human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
  mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human,
                   attributesL = c("mgi_symbol"), martL = mouse, uniqueRows = T)
  humanx <- unique(genesV2[, 2])
  return(humanx)
}

hks <- readxl::read_xlsx("./Data/processed/classifier/gkaa609_supplemental_files/Supplementary_Table1.xlsx", skip = 1)
hks <- convertHumanGeneList(as.character(hks$Gene))
VAR_expression_isHousekeeper <- genes_test_vector %in% hks

# ---- ## ---- ## ---- ## ---- #  
# ---- # Promoter architecture features
# ---- ## ---- ## ---- ## ---- #
# Which genes have TATA boxes or Initiator elements?
# Need to kick out genes for which there is no annotation
universe <- read.table("./Data/processed/classifier/Promoter_features/mouse_epdnew_QsG7N.ALL.bed")
all_genes_promoter <- gsub("_.*$", "", universe$V4)
tata.data <- read.table("./Data/processed/classifier/Promoter_features/mouse_epdnew_Ow9Mk.TATA.bed", sep = "\t")
tata.data$gene <- gsub("_.*$", "", tata.data$V4)

VAR_promoter_hasTATA <- rownames(genes_test) %in% tata.data$gene
VAR_promoter_hasTATA[!rownames(genes_test) %in% all_genes_promoter] <- NA

inr.data <- read.table("./Data/processed/classifier/Promoter_features/mouse_epdnew_7Qe23.INR.bed", sep = "\t")
inr.data$gene <- gsub("_.*$", "", inr.data$V4)

VAR_promoter_hasINR <- rownames(genes_test) %in% inr.data$gene
VAR_promoter_hasINR[!rownames(genes_test) %in% all_genes_promoter] <- NA

# GC content

#Get the sequences and compute the GC content
library(BSgenome.Mmusculus.UCSC.mm10)
GetGC <- function(bsgenome, gr){
  seqs <- BSgenome::getSeq(bsgenome, gr)
  return(as.numeric(Biostrings::letterFrequency(x = seqs, letters = "GC", as.prob = TRUE)))
}

promoter_targetgenes <- gene_coordinates$Promoter
promoter_targetgenes <- promoter_targetgenes[promoter_targetgenes$symbol %in% rownames(genes_test), ]
promoter_targetgenes_here <- keepStandardChromosomes(promoter_targetgenes, pruning.mode = "coarse")
promoter_targetgenes_here <- promoter_targetgenes_here[seqnames(promoter_targetgenes_here) != "MT", ]
promoter_targetgenes_here <- dropSeqlevels(promoter_targetgenes_here, c("MT"))
seqlevelsStyle(promoter_targetgenes_here) <- "UCSC"
promoter_targetgenes_here$gc_content_promoters <- GetGC(BSgenome.Mmusculus.UCSC.mm10, promoter_targetgenes_here)

df_in <- data.frame(gene = promoter_targetgenes_here$symbol, gc = promoter_targetgenes_here$gc_content_promoters)
df_in <- aggregate(df_in, list(df_in$gene), mean) %>% 
  column_to_rownames("Group.1")
VAR_promoter_gcContent <- df_in[rownames(genes_test), ]$gc

# CpG islands overlaying promoters
library(rtracklayer)
mySession = browserSession("UCSC")
genome(mySession) <- "mm10"
promoter_targetgenes <- gene_coordinates$Promoter
promoter_targetgenes <- promoter_targetgenes[promoter_targetgenes$symbol %in% rownames(genes_test), ]
tbl.cpg_island <- getTable(ucscTableQuery(mySession, track="cpgIslandExt", table="cpgIslandExt"))
tbl.cpg_island <- makeGRangesFromDataFrame(tbl.cpg_island)
seqlevelsStyle(promoter_targetgenes) <- "UCSC"
promoter_targetgenes <- keepStandardChromosomes(promoter_targetgenes, pruning.mode = "coarse")
genes_with_CpG <- unique(subsetByOverlaps(promoter_targetgenes, tbl.cpg_island)$symbol)
VAR_promoter_has_CpG <- rownames(genes_test) %in% genes_with_CpG

# ---- ## ---- ## ---- ## ---- #  
# ---- # Other structural features
# ---- ## ---- ## ---- ## ---- #

# get number of exons
exons_here <- gene_coordinates$Exon
introns_here <- gene_coordinates$Intron
exon_number <- table(exons_here$symbol)
VAR_structural_nExons <- exon_number[genes_test_vector]

# get gene lengths 
exons_lengths <- aggregate(width(exons_here), by = list(exons_here$symbol), sum)
intron_lengths <- aggregate(width(introns_here), by = list(introns_here$symbol), sum)
both_lengths <- merge(exons_lengths, intron_lengths, by = "Group.1", all = T)
both_lengths[is.na(both_lengths)] <- 0
both_lengths$total_length <- both_lengths$x.x + both_lengths$x.y
both_lengths <- both_lengths[both_lengths$Group.1 %in% genes_test_vector, ]
rownames(both_lengths) <- both_lengths$Group.1
both_lengths <- both_lengths[genes_test_vector, ]

VAR_structural_gene_length <- both_lengths$total_length
VAR_structural_exon_proportion <- both_lengths$x.x / VAR_structural_gene_length
VAR_structural_intron_proportion <- both_lengths$x.y / VAR_structural_gene_length

# GC content exons / introns
exons_here_cut <- exons_here[exons_here$symbol %in% rownames(genes_test), ]
introns_here_cut <- introns_here[introns_here$symbol %in% rownames(genes_test), ]

seqlevelsStyle(exons_here_cut) <- "UCSC"
seqlevelsStyle(introns_here_cut) <- "UCSC"
exons_here_cut <- keepStandardChromosomes(dropSeqlevels(exons_here_cut, value = "chrM"), pruning.mode = "coarse")
introns_here_cut <- keepStandardChromosomes(dropSeqlevels(introns_here_cut, value = "chrM"), pruning.mode = "coarse")

exons_here_cut$gc <- GetGC(BSgenome.Mmusculus.UCSC.mm10, exons_here_cut)
introns_here_cut$gc <- GetGC(BSgenome.Mmusculus.UCSC.mm10, introns_here_cut)

gc_exon_aggregate <- exons_here_cut %>% 
  data.frame() %>%
  as_tibble() %>%
  dplyr::select(symbol, gc) %>%
  group_by(symbol) %>%
  summarize(gc = mean(gc))

gc_intron_aggregate <- introns_here_cut %>% 
  data.frame() %>%
  as_tibble() %>%
  dplyr::select(symbol, gc) %>%
  group_by(symbol) %>%
  summarize(gc = mean(gc))

VAR_structural_ExonGC <- rep(0.5, length(genes_test_vector))
names(VAR_structural_ExonGC) <- genes_test_vector
VAR_structural_ExonGC[gc_exon_aggregate$symbol] <- gc_exon_aggregate$gc

VAR_structural_IntronGC <- rep(0.5, length(genes_test_vector))
names(VAR_structural_IntronGC) <- genes_test_vector
VAR_structural_IntronGC[gc_intron_aggregate$symbol] <- gc_intron_aggregate$gc

# number of transcripts
library(ensembldb)
library(EnsDb.Mmusculus.v79)

all_transcripts <- transcripts(EnsDb.Mmusculus.v79)
gene_symbols <- ensembldb::mapIds(EnsDb.Mmusculus.v79, all_transcripts$tx_name, "SYMBOL", "TXNAME")
all_transcripts$gene_symbol <- gene_symbols
df_out <- data.frame(tx = all_transcripts$tx_name, symbols = all_transcripts$gene_symbol)
df_out <- unique(df_out)
VAR_structural_nTranscripts <- table(df_out$symbols)[rownames(genes_test)]

# UTRs
all_five_utrs <- gene_coordinates$FiveUTR
five_utr_lengths <- data.frame(symbol = all_five_utrs$symbol, utr_length = width(all_five_utrs)) %>%
  group_by(symbol) %>% summarize(utr_length = mean(utr_length)) %>% column_to_rownames("symbol")
VAR_structural_fiveUTRlength <- five_utr_lengths[genes_test_vector, ]

all_three_utrs <- gene_coordinates$ThreeUTR
three_utr_lengths <- data.frame(symbol = all_three_utrs$symbol, utr_length = width(all_three_utrs)) %>%
  group_by(symbol) %>% summarize(utr_length = mean(utr_length)) %>% column_to_rownames("symbol")
VAR_structural_threeUTRlength <- three_utr_lengths[genes_test_vector, ]

############## 

library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(AnnotationHub)

# hgdownload.cse.ucsc.edu/goldenPath/mm10/database/cpgIslandExt.txt.gz
#cpgs <- read_delim("~/Downloads/cpgIslandExt(1).txt.gz", delim  = "\t", col_names = c("Number", "Chromosome", "Start", "End", "Annotation", "NA", "NA", "NA", "NA", "NA", "NA", "NA"))
#cpgs <- makeGRangesFromDataFrame(cpgs)

# promoters_genes <- promoters(TxDb.Mmusculus.UCSC.mm10.knownGene, upstream = 500, downstream = 100)
# cpg_promoter_overlaps <- subsetByOverlaps(promoters_genes, cpgs)
# 
# ## CTCF binding sites
# hub <- AnnotationHub()
# ctcf_sites <- query(hub, c("CTCF_mm10"))
# ctcf_sites <- hub["AH95568"]

# ---- ## ---- ## ---- ## ---- #  
# ---- # Evolutionary features
# ---- ## ---- ## ---- ## ---- #

library(GenomicScores)

# PhyloP60 scores
cons_scores <- getGScores("phyloP60way.UCSC.mm10")
phastcons_scores <-  getGScores('phastCons60way.UCSC.mm10')

cons_input <- lapply(gene_coordinates, function(x){x[x$symbol %in% genes_test_vector, ]})

cons_input <- lapply(cons_input, function(x){
  seqlevelsStyle(x) <- "UCSC"
  x <- keepStandardChromosomes(x, pruning.mode = "coarse"); x})
scores_exons <- gscores(phastcons_scores, cons_input$Exon)
scores_introns <- gscores(phastcons_scores, cons_input$Intron)
scores_five_utr <- gscores(phastcons_scores, cons_input$FiveUTR)
scores_three_utr <- gscores(phastcons_scores, cons_input$ThreeUTR)
scores_promoter <- gscores(phastcons_scores, cons_input$Promoter)

phastcons_scores_exon <- scores_exons %>% data.frame() %>%
  dplyr::select(symbol, default) %>%
  group_by(symbol) %>% summarize(value = mean(default, na.rm = T)) %>%
  column_to_rownames("symbol")
phastcons_scores_intron <- scores_introns %>% data.frame() %>%
  dplyr::select(symbol, default) %>%
  group_by(symbol) %>% summarize(value = mean(default, na.rm = T)) %>%
  column_to_rownames("symbol")
phastcons_scores_five_utr <- scores_five_utr %>% data.frame() %>%
  dplyr::select(symbol, default) %>%
  group_by(symbol) %>% summarize(value = mean(default, na.rm = T)) %>%
  column_to_rownames("symbol")
phastcons_scores_three_utr <- scores_three_utr %>% data.frame() %>%
  dplyr::select(symbol, default) %>%
  group_by(symbol) %>% summarize(value = mean(default, na.rm = T)) %>%
  column_to_rownames("symbol")
phastcons_scores_promoter <- scores_promoter %>% data.frame() %>%
  dplyr::select(symbol, default) %>%
  group_by(symbol) %>% summarize(value = mean(default, na.rm = T)) %>%
  column_to_rownames("symbol")

VAR_evolution_phastcons_exon <- phastcons_scores_exon[genes_test_vector, ]
VAR_evolution_phastcons_intron <- phastcons_scores_intron[genes_test_vector, ]
VAR_evolution_phastcons_five_utr <- phastcons_scores_five_utr[genes_test_vector, ]
VAR_evolution_phastcons_three_utr <- phastcons_scores_three_utr[genes_test_vector, ]
VAR_evolution_phastcons_promoter <- phastcons_scores_promoter[genes_test_vector, ]

phylop_scores_exons <- gscores(cons_scores, cons_input$Exon)
phylop_scores_introns <- gscores(cons_scores, cons_input$Intron)
phylop_scores_five_utr <- gscores(cons_scores, cons_input$FiveUTR)
phylop_scores_three_utr <- gscores(cons_scores, cons_input$ThreeUTR)
phylop_scores_promoter <- gscores(cons_scores, cons_input$Promoter)

phylop_scores_exon <- phylop_scores_exons %>% data.frame() %>%
  dplyr::select(symbol, default) %>%
  group_by(symbol) %>% summarize(value = mean(default, na.rm = T)) %>%
  column_to_rownames("symbol")
phylop_scores_intron <- phylop_scores_introns %>% data.frame() %>%
  dplyr::select(symbol, default) %>%
  group_by(symbol) %>% summarize(value = mean(default, na.rm = T)) %>%
  column_to_rownames("symbol")
phylop_scores_five_utr <- phylop_scores_five_utr %>% data.frame() %>%
  dplyr::select(symbol, default) %>%
  group_by(symbol) %>% summarize(value = mean(default, na.rm = T)) %>%
  column_to_rownames("symbol")
phylop_scores_three_utr <- phylop_scores_three_utr %>% data.frame() %>%
  dplyr::select(symbol, default) %>%
  group_by(symbol) %>% summarize(value = mean(default, na.rm = T)) %>%
  column_to_rownames("symbol")
phylop_scores_promoter <- phylop_scores_promoter %>% data.frame() %>%
  dplyr::select(symbol, default) %>%
  group_by(symbol) %>% summarize(value = mean(default, na.rm = T)) %>%
  column_to_rownames("symbol")

VAR_evolution_phylop_exon <- phylop_scores_exon[genes_test_vector, ]
VAR_evolution_phylop_intron <- phylop_scores_intron[genes_test_vector, ]
VAR_evolution_phylop_five_utr <- phylop_scores_five_utr[genes_test_vector, ]
VAR_evolution_phylop_three_utr <- phylop_scores_three_utr[genes_test_vector, ]
VAR_evolution_phylop_promoter <- phylop_scores_promoter[genes_test_vector, ]

# SNP density in CAST / across mouse strains
read_snp_file <- function(x){
  snps.chr1 <- read.csv(x, skip = 1, sep = "\t", header = F)
  colnames(snps.chr1) <- c("ID", "Chr", "Start", "Strand", "Alleles", "Info")
  snps.chr1$Strand <- "+"
  snps.chr1$Chr <- paste0("chr", snps.chr1$Chr)
  snps.chr1 <- makeGRangesFromDataFrame(snps.chr1, seqnames.field = "Chr", start.field = "Start", end.field = "Start")
  snps.chr1
}

all_snps <- lapply(list.files("./Stuff/rforest_analysis_files/Cast_SNPs/", full.names = T), read_snp_file)
all_snps_combined <- do.call("c", all_snps)

compute_snp_density <- function(df_here, genes){
  df_here$n_snps <- as.numeric(countOverlaps(df_here, all_snps_combined))
  df_here$snp_rate <- df_here$n_snps / width(df_here)
  df_here <- df_here %>% data.frame() %>% dplyr::select(symbol, snp_rate) %>% group_by(symbol) %>%
    summarize(snp_rate = mean(snp_rate, na.rm = T)) %>% column_to_rownames("symbol")
  df_here[genes, ]
}

VAR_evolution_snpdensity_F1_exon <- compute_snp_density(cons_input$Exon, genes_test_vector)
VAR_evolution_snpdensity_F1_intron <- compute_snp_density(cons_input$Intron, genes_test_vector)
VAR_evolution_snpdensity_F1_fiveUTR <- compute_snp_density(cons_input$FiveUTR, genes_test_vector)
VAR_evolution_snpdensity_F1_threeUTR <- compute_snp_density(cons_input$ThreeUTR, genes_test_vector)
VAR_evolution_snpdensity_F1_promoter <- compute_snp_density(cons_input$Promoter, genes_test_vector)

############## 
# Epigenomic features (testis)
############## 

testis_specific_promoter <- readRDS("./Data/processed/classifier/testis_specific_enhancer_active.rds") %>%
  mutate(chr = paste0("chr", chr)) %>% makeGRangesFromDataFrame()
testis_specific_enhancer_primed <- readRDS("./Data/processed/classifier/testis_specific_enhancer_active.rds") %>%
  mutate(chr = paste0("chr", chr)) %>% makeGRangesFromDataFrame()
testis_specific_enhancer_active <- readRDS("./Data/processed/classifier/testis_specific_enhancer_active.rds") %>%
  mutate(chr = paste0("chr", chr)) %>% makeGRangesFromDataFrame()
tissue_shared_promoter <- readRDS("./Data/processed/classifier/tissue_shared_promoter.rds") %>%
  mutate(chr = paste0("chr", chr)) %>% makeGRangesFromDataFrame()
tissue_shared_enhancer_primed <- readRDS("./Data/processed/classifier/tissue_shared_promoter.rds") %>%
  mutate(chr = paste0("chr", chr)) %>% makeGRangesFromDataFrame()
tissue_shared_enhancer_active <- readRDS("./Data/processed/classifier/tissue_shared_promoter.rds") %>%
  mutate(chr = paste0("chr", chr)) %>% makeGRangesFromDataFrame()

library(tidyverse)

# find overlaps with promoters
# testis_specific
input_here <- cons_input$Promoter
input_here$counts <- countOverlaps(input_here, testis_specific_promoter)
VAR_REGULATION_has_testis_specific_promoter <- data.frame(symbol = input_here$symbol, counts = input_here$counts) %>%
  group_by(symbol) %>% summarize(has_feature = sum(counts) > 0) %>% column_to_rownames("symbol")

# tissue_shared
input_here <- cons_input$Promoter
input_here$counts <- countOverlaps(input_here, tissue_shared_promoter)
VAR_REGULATION_has_tissue_shared_promoter <- data.frame(symbol = input_here$symbol, counts = input_here$counts) %>%
  group_by(symbol) %>% summarize(has_feature = sum(counts) > 0) %>% column_to_rownames("symbol")

# find overlaps with testis specific enhancers
# active
enhancer_range <- 1e5

input_here <- cons_input$Promoter
start(input_here) <- start(input_here) - enhancer_range
end(input_here) <- end(input_here) + enhancer_range
input_here$counts <- countOverlaps(input_here, testis_specific_enhancer_active)
VAR_REGULATION_n_testis_specific_active_enhancers <- data.frame(symbol = input_here$symbol, counts = input_here$counts) %>%
  group_by(symbol) %>% summarize(n_enhancer = sum(counts)) %>% column_to_rownames("symbol")
# primed
input_here <- cons_input$Promoter
start(input_here) <- start(input_here) - enhancer_range
end(input_here) <- end(input_here) + enhancer_range
input_here$counts <- countOverlaps(input_here, testis_specific_enhancer_primed)
VAR_REGULATION_n_testis_specific_primed_enhancers <- data.frame(symbol = input_here$symbol, counts = input_here$counts) %>%
  group_by(symbol) %>% summarize(n_enhancer = sum(counts)) %>% column_to_rownames("symbol")

# find overlaps with tissue shared enhancers
# active
input_here <- cons_input$Promoter
start(input_here) <- start(input_here) - enhancer_range
end(input_here) <- end(input_here) + enhancer_range
input_here$counts <- countOverlaps(input_here, tissue_shared_enhancer_active)
VAR_REGULATION_n_tissue_shared_active_enhancers <- data.frame(symbol = input_here$symbol, counts = input_here$counts) %>%
  group_by(symbol) %>% summarize(n_enhancer = sum(counts)) %>% column_to_rownames("symbol")
# primed
input_here <- cons_input$Promoter
start(input_here) <- start(input_here) - enhancer_range
end(input_here) <- end(input_here) + enhancer_range
input_here$counts <- countOverlaps(input_here, tissue_shared_enhancer_primed)
VAR_REGULATION_n_tissue_shared_primed_enhancers <- data.frame(symbol = input_here$symbol, counts = input_here$counts) %>%
  group_by(symbol) %>% summarize(n_enhancer = sum(counts)) %>% column_to_rownames("symbol")

# Build Model
# which genes in which groups?

genes_test <- readRDS("./Data/processed/genes_test.rds")
genes_save <- rownames(genes_test)
genes_test <- data.frame(apply(genes_test, 2, as.numeric))
rownames(genes_test) <- genes_save
genes_test <- genes_test[!is.na(rownames(genes_test)), ]
genes_test <- genes_test[!is.na(genes_test$chromosome), ]
genes_test <- genes_test[!is.na(genes_test$BF_trans), ]
genes_test <- genes_test[!is.na(genes_test$BF_cis_dynamic), ]
genes_test <- genes_test[!is.na(genes_test$BF_trans_static), ]
genes_test <- genes_test[!is.na(genes_test$BF_trans_dynamic), ]

genes_test$BF_cis_static <- genes_test$NLL_null - genes_test$NLL_alt
genes_test[genes_test$NLL_null < 0 | genes_test$NLL_alt < 0, ]$BF_cis_static <- 0
genes_test[genes_test$BF_cis_static < 0, ]$BF_cis_static <- 0

BF_cutoff <- 10
BF_cutoff_lower <- 2
genes_cis_static <- rownames(genes_test[genes_test$BF_cis_static > BF_cutoff & genes_test$BF_cis_dynamic, ])
genes_cis_dynamic <- rownames(genes_test[genes_test$BF_cis_dynamic > BF_cutoff, ])
genes_cis <- union(genes_cis_static, genes_cis_dynamic)
genes_trans_static <- rownames(genes_test[genes_test$BF_trans_static > BF_cutoff & genes_test$BF_trans_dynamic, ])
genes_trans_dynamic <- rownames(genes_test[genes_test$BF_trans_dynamic > BF_cutoff, ])
genes_trans <- union(genes_trans_static, genes_trans_dynamic)
genetic_effects <- union(genes_cis, genes_trans)

# Build dataset
all_genes <- names(VAR_expression_level)

zero_log <- function(x){
  x <- x - min(x, na.rm = T)
  log( x + min(x[x > 0], na.rm = T) )
}

data_fit_list <- list(
  # promoter features
  VAR_promoter_hasTATA = as.numeric(VAR_promoter_hasTATA),
  VAR_promoter_hasINR = as.numeric(VAR_promoter_hasINR), 
  VAR_promoter_gcContent = VAR_promoter_gcContent, 
  VAR_promoter_has_CpG = VAR_promoter_has_CpG, 
  # gene context features
  # distal regulator features
  # gene body features
  VAR_structural_gene_length = log(VAR_structural_gene_length + 1), 
  VAR_structural_exon_proportion = VAR_structural_exon_proportion, 
  VAR_structural_exon_number = log( as.numeric(VAR_structural_nExons) + 1), 
  VAR_structural_fiveUTRlength = log( VAR_structural_fiveUTRlength + 1), 
  VAR_structural_threeUTRlength = log( VAR_structural_threeUTRlength + 1) , 
  VAR_structural_ExonGC = VAR_structural_ExonGC, 
  VAR_structural_IntronGC = VAR_structural_IntronGC, 
  # UTR features
  # evolution / genetics features
  VAR_evolution_phylop_exon = VAR_evolution_phylop_exon, 
  VAR_evolution_phylop_intron = VAR_evolution_phylop_intron, 
  VAR_evolution_phylop_five_utr = VAR_evolution_phylop_five_utr, 
  VAR_evolution_phylop_three_utr = VAR_evolution_phylop_three_utr, 
  VAR_evolution_phylop_promoter = VAR_evolution_phylop_promoter, 
  VAR_evolution_phastcons_exon = VAR_evolution_phastcons_exon, 
  VAR_evolution_phastcons_intron = VAR_evolution_phastcons_intron, 
  VAR_evolution_phastcons_five_utr = VAR_evolution_phastcons_five_utr, 
  VAR_evolution_phastcons_three_utr = VAR_evolution_phastcons_three_utr, 
  VAR_evolution_phastcons_promoter = VAR_evolution_phastcons_promoter, 
  VAR_evolution_snpdensity_F1_exon = zero_log(VAR_evolution_snpdensity_F1_exon), 
  VAR_evolution_snpdensity_F1_intron = zero_log(VAR_evolution_snpdensity_F1_intron), 
  VAR_evolution_snpdensity_F1_fiveUTR = zero_log(VAR_evolution_snpdensity_F1_fiveUTR), 
  VAR_evolution_snpdensity_F1_threeUTR = zero_log(VAR_evolution_snpdensity_F1_threeUTR), 
  VAR_evolution_snpdensity_F1_promoter = zero_log(VAR_evolution_snpdensity_F1_promoter), 
  # expression features
  VAR_expression_level = log( VAR_expression_level),
  VAR_expression_variability = zero_log(VAR_expression_variability),
  VAR_expression_variability_post_celltype = zero_log(VAR_expression_variability_post_celltype), 
  VAR_expression_isHousekeeper = VAR_expression_isHousekeeper, 
  # gene regulatory features
  VAR_REGULATION_has_testis_specific_promoter = VAR_REGULATION_has_testis_specific_promoter[all_genes, 1], 
  VAR_REGULATION_has_tissue_shared_promoter = VAR_REGULATION_has_tissue_shared_promoter[all_genes, 1], 
  VAR_REGULATION_n_testis_specific_active_enhancers = log(VAR_REGULATION_n_testis_specific_active_enhancers[all_genes, 1] + 1), 
  VAR_REGULATION_n_testis_specific_primed_enhancers = log(VAR_REGULATION_n_testis_specific_primed_enhancers[all_genes, 1] + 1), 
  VAR_REGULATION_n_tissue_shared_active_enhancers = log(VAR_REGULATION_n_tissue_shared_active_enhancers[all_genes, 1] + 1), 
  VAR_REGULATION_n_tissue_shared_primed_enhancers = log(VAR_REGULATION_n_tissue_shared_primed_enhancers[all_genes, 1] + 1)
)

data_fit <- do.call("cbind", data_fit_list)

saveRDS(data_fit, "./Data/processed/all_features.rds")

