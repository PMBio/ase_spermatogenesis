### ASE_paper functions

plot_colours <- function(colours){
  # mod from https://stackoverflow.com/questions/51867716/plot-colors-with-hex-values-in-r
  # take colors, name the column, then cbind
  
  color_df = data.frame(
    ColorName = names(colours),
    ColorHex = colours
  )
  
  # plot it, using scale_color_identity to interpret the hex as a color
  ggplot(color_df, aes(x = ColorName, y = 1, fill = ColorHex)) +
    geom_tile() +
    geom_text(aes(label = ColorName)) +
    scale_fill_identity() + 
    coord_flip() + theme_void()
}

annotate_chromosome <- function(sce){
  library(biomaRt)
  ensembl <- useMart("ensembl")
  ensembl <- useDataset(dataset = "mmusculus_gene_ensembl", 
                        mart = ensembl)
  genenames <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name"),
                     mart = ensembl)
  rownames(genenames) <- make.unique(genenames$external_gene_name)
  genenames <- genenames[!duplicated(genenames),]
  rowData(sce) <- genenames[rownames(sce),]
  return(sce)
}

annotate_chromosome_list <- function(genes){
  library(biomaRt)
  ensembl <- useMart("ensembl")
  ensembl <- useDataset(dataset = "mmusculus_gene_ensembl", 
                        mart = ensembl)
  genenames <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name"),
                     mart = ensembl)
  rownames(genenames) <- make.unique(genenames$external_gene_name)
  genenames <- genenames[!duplicated(genenames),]
  return(genenames[genes,])
}

counts_reference <- function(sce){
  assays(sce)[["counts_reference"]]
}

counts_alternative <- function(sce){
  assays(sce)[["counts_alternative"]]
}

theme_paper <- function(textsize = 20){
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1), 
        text = element_text(size = textsize))
}

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

annotate_chromosome_new <- function(to_annotate, by_identifier = "gene_symbol"){
  library(EnsDb.Mmusculus.v79)
  all_genes <- genes(EnsDb.Mmusculus.v79)
  if ( by_identifier == "gene_symbol" ){
    df_here <- data.frame(
      chromosome_name = seqnames(all_genes), 
      genes = all_genes$symbol
    )
    df_here <- df_here[!duplicated(df_here$genes), ]
    rownames(df_here) <- df_here$genes
  }
  df_here[to_annotate, ]
}

annotate_chromosome_sce <- function(to_annotate, by_identifier = "gene_symbol"){
  library(EnsDb.Mmusculus.v79)
  all_genes <- genes(EnsDb.Mmusculus.v79)
  if ( by_identifier == "gene_symbol" ){
    df_here <- data.frame(
      chromosome_name = seqnames(all_genes), 
      genes = all_genes$symbol
    )
    df_here <- df_here[!duplicated(df_here$genes), ]
    rownames(df_here) <- df_here$genes
  }
  new_coldata <- cbind(rowData(to_annotate), "chromosome_name" = df_here[rownames(to_annotate), 1])
  rowData(to_annotate) <- new_coldata
  to_annotate
}

plot_gene_GP <- function(sce, gene, latent_ase, remove_zero = F, scale_var = 1.96, textsize = 20){
  data_test <- data.frame(
    pt_here = sce$Pseudotime,
    exp_ref = as.numeric(counts_reference(sce[gene,])),
    exp_alt = as.numeric(counts_alternative(sce[gene,])),
    #ase = as.numeric(allelic_ratios(sce[gene,])),
    latent_ase = latent_ase$posterior_mean,
    latent_var_lower = latent_ase$posterior_mean - scale_var * sqrt(latent_ase$posterior_var),
    latent_var_upper = latent_ase$posterior_mean + scale_var * sqrt(latent_ase$posterior_var)
  )
  data_test$ase <- data_test$exp_ref / (data_test$exp_alt + data_test$exp_ref)
  
  if (remove_zero){
    data_test <- data_test[data_test$exp_ref + data_test$exp_alt > 0,]
  }
  
  data_test$exp_total = data_test$exp_ref + data_test$exp_alt
  data_test <- data_test[order(data_test$pt),]
  
  print(summary(data_test$pt_here))
  
  p1 <- ggplot(data_test) + 
    geom_jitter(aes(pt_here, exp_ref, color = "B6"), alpha = 0.2, size = 0.2, height = 0.01) + 
    geom_jitter(aes(pt_here, exp_alt, color = "Cast"), alpha = 0.2, size = 0.2, height = 0.01) + 
    geom_smooth(aes(pt_here, exp_ref, color = "B6"), size = 1) + 
    geom_smooth(aes(pt_here, exp_alt, color = "Cast"), size = 1) + 
    theme_classic() + xlim(c(0, max(data_test$pt_here + 1))) + xlab("") + ylab("Expression") +
    theme(legend.position="top") + 
    scale_y_log10(expand = c(0, 0)) + 
    scale_color_manual(values = c("B6" = "black", "Cast" = "chocolate")) + 
    theme(text = element_text(size = textsize)) +
    theme(axis.line.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank()) + 
    ggtitle(gene)
  p1
  
  p2 <- ggplot(data_test, aes(pt_here, ase)) +
    ylim(c(0, 1)) + geom_point(width = 0.1, height = 0.02, size = 0.1) + 
    scale_color_viridis() + 
    geom_line(aes(pt_here, latent_ase), color = "green") + 
    geom_ribbon(aes(x = pt_here, ymin = latent_var_lower, ymax = latent_var_upper), 
                color = "green", alpha = 0.2) + 
    theme_classic() + xlim(c(0, max(data_test$pt_here + 1))) + 
    xlab("Pseudotime") + ylab("Allelic Bias") + 
    theme(text = element_text(size = textsize)) + 
    scale_y_continuous(breaks = c(0, 0.5, 1), expand = c(0, 0)) + 
    theme(axis.line.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank())
  
  
  cowplot::plot_grid(p1, p2, 
                     nrow = 2, align = "v", 
                     rel_heights = c(0.4, 0.6)) 
}

logit <- function(x){
  log(x / (1 - x))
}

rev_logit <- function(x){
  1 / (1 + exp(-x))
}

convert_ensembl_symbol <- function(data){
  library(EnsDb.Mmusculus.v79)
  txdb <- EnsDb.Mmusculus.v79
  k <- keys(txdb, keytype = "GENEID")
  tx2gene <- AnnotationDbi::select(txdb, k, "SYMBOL", "GENEID")
  joint.genes <- intersect(rownames(data), tx2gene$GENEID)
  tx2gene.filt <- tx2gene[tx2gene$GENEID %in% joint.genes,]
  data.filt <- data[rownames(data) %in% joint.genes,]
  rownames(tx2gene.filt) <- tx2gene.filt$GENEID
  rownames(data.filt) <- make.unique(tx2gene.filt[joint.genes,]$SYMBOL)
  data.filt
}

# function to read allele specific counts and returns as list
read_10x_ase <- function(sample.dir){
  barcodes <- read.csv(paste0(sample.dir, "/barcodes.reference.tsv"), header = F)
  genes <- read.csv(paste0(sample.dir, "/genes.reference.tsv"), header = F)
  matrix <- read.csv(paste0(sample.dir, "/matrix.reference.tsv"), header = F, sep = "\t")
  matrix.wide <- reshape2::dcast(matrix, V1 ~ V2, fill = 0)
  rownames(matrix.wide) <- matrix.wide[,1]
  matrix.wide.reference <- t(matrix.wide[,-1])
  matrix.wide.reference <- convert_ensembl_symbol(matrix.wide.reference)
  
  barcodes <- read.csv(paste0(sample.dir, "/barcodes.alternative.tsv"), header = F)
  genes <- read.csv(paste0(sample.dir, "/genes.alternative.tsv"), header = F)
  matrix <- read.csv(paste0(sample.dir, "/matrix.alternative.tsv"), header = F, sep = "\t")
  matrix.wide <- reshape2::dcast(matrix, V1 ~ V2, fill = 0)
  rownames(matrix.wide) <- matrix.wide[,1]
  matrix.wide.alternative <- t(matrix.wide[,-1])
  matrix.wide.alternative <- convert_ensembl_symbol(matrix.wide.alternative)
  
  joint.bcs <- union(colnames(matrix.wide.reference), colnames(matrix.wide.alternative))
  
  to.add.reference <- matrix(rep(0, (length(joint.bcs) - ncol(matrix.wide.reference)) * 
                                   nrow(matrix.wide.reference)),
                             nrow = nrow(matrix.wide.reference),
                             ncol = length(joint.bcs) - ncol(matrix.wide.reference))
  colnames(to.add.reference) <- joint.bcs[!joint.bcs %in% colnames(matrix.wide.reference)]
  
  to.add.alternative <- matrix(rep(0, (length(joint.bcs) - ncol(matrix.wide.alternative)) * 
                                     nrow(matrix.wide.alternative)),
                               nrow = nrow(matrix.wide.alternative),
                               ncol = length(joint.bcs) - ncol(matrix.wide.alternative))
  colnames(to.add.alternative) <- joint.bcs[!joint.bcs %in% colnames(matrix.wide.alternative)]
  
  matrix.wide.reference <- cbind(matrix.wide.reference, to.add.reference)
  matrix.wide.alternative <- cbind(matrix.wide.alternative, to.add.alternative)
  
  matrix.wide.reference = matrix.wide.reference[rownames(matrix.wide.reference) != "",]
  matrix.wide.alternative = matrix.wide.alternative[rownames(matrix.wide.alternative) != "",]
  
  return(list("reference" = matrix.wide.reference, "alternative" = matrix.wide.alternative))
}

annotate_chromosome <- function(sce){
  library(biomaRt)
  ensembl <- useMart("ensembl")
  ensembl <- useDataset(dataset = "mmusculus_gene_ensembl", 
                        mart = ensembl)
  genenames <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name"),
                     mart = ensembl)
  rownames(genenames) <- make.unique(genenames$external_gene_name)
  genenames <- genenames[!duplicated(genenames),]
  rowData(sce) <- genenames[rownames(sce),]
  return(sce)
}

annotate_chromosome_list <- function(genes){
  library(biomaRt)
  ensembl <- useMart("ensembl")
  ensembl <- useDataset(dataset = "mmusculus_gene_ensembl", 
                        mart = ensembl)
  genenames <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name"),
                     mart = ensembl)
  rownames(genenames) <- make.unique(genenames$external_gene_name)
  genenames <- genenames[!duplicated(genenames),]
  return(genenames[genes,])
}

counts_reference <- function(sce){
  assays(sce)[["counts_reference"]]
}

counts_alternative <- function(sce){
  assays(sce)[["counts_alternative"]]
}


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

annotate_chromosome <- function(sce){
  library(biomaRt)
  ensembl <- useMart("ensembl")
  ensembl <- useDataset(dataset = "mmusculus_gene_ensembl", 
                        mart = ensembl)
  genenames <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name"),
                     mart = ensembl)
  rownames(genenames) <- make.unique(genenames$external_gene_name)
  genenames <- genenames[!duplicated(genenames),]
  rowData(sce) <- genenames[rownames(sce),]
  return(sce)
}
annotate_chromosome_list <- function(genes){
  library(biomaRt)
  ensembl <- useMart("ensembl")
  ensembl <- useDataset(dataset = "mmusculus_gene_ensembl", 
                        mart = ensembl)
  genenames <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name"),
                     mart = ensembl)
  rownames(genenames) <- make.unique(genenames$external_gene_name)
  genenames <- genenames[!duplicated(genenames),]
  return(genenames[genes,])
}
counts_reference <- function(sce){
  assays(sce)[["counts_reference"]]
}
counts_alternative <- function(sce){
  assays(sce)[["counts_alternative"]]
}