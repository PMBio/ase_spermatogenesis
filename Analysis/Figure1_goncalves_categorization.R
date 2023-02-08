library(scran)
library(scater)
library(tidyverse)

setwd("~/Desktop/Projects/ASE_Spermatogenesis_Paper_rerun/Revisions/")

source("../Scripts/General/auxiliary.R")
source("../Scripts/General/reuse_functions.R")

colours_cis_trans <- c(
  "conserved" =  "grey",
  "cis" =  "orange",
  "trans" =  "cyan",
  "cis_trans" =  "Purple"
)

data <- readRDS("./Data/processed/sce_merged_new.rds")
data$IndividualSamples <- paste0(data$Dataset, "_", data$Library)
data <- data[,data$Library != "Sample7"]

aggregate_test <- aggregateAcrossCells(data, data$IndividualSamples, use.assay.type = c("counts_reference", "counts_alternative"))

normalize_sf <- function(counts){
  col.save <- colnames(counts)
  row.save <- rownames(counts)
  sfs <- colSums(counts) / colSums(counts)[[1]]
  counts <- do.call("cbind", lapply(1:ncol(counts), function(i){
    return(counts[,i] / sfs[[i]])
  }))
  colnames(counts) <- col.save
  rownames(counts) <- row.save
  return(round(counts))
}
make_bulk_celltype <- function(data){
  
  aggregate_sce <- aggregateAcrossCells(data, data$IndividualSamples, use.assay.type = c("counts_reference", "counts_alternative"))
  
  bulk_par_b6 <- counts_reference(aggregate_sce[,aggregate_sce$Species == "B6"])
  colnames(bulk_par_b6) <- paste0(colnames(bulk_par_b6), "_Reference")
  bulk_par_cast <- counts_alternative(aggregate_sce[,aggregate_sce$Species == "CAST"])
  colnames(bulk_par_cast) <- paste0(colnames(bulk_par_cast), "_Alternative")
  
  sce_par_cur_perLibrary <- cbind(bulk_par_b6, bulk_par_cast)
  
  bulk_fil_b6 <- counts_reference(aggregate_sce[,aggregate_sce$Species == "B6xCAST"])
  colnames(bulk_fil_b6) <- paste0(colnames(bulk_fil_b6), "_Reference")
  bulk_fil_cast <- counts_alternative(aggregate_sce[,aggregate_sce$Species == "B6xCAST"])
  colnames(bulk_fil_cast) <- paste0(colnames(bulk_fil_cast), "_Alternative")
  
  sce_fil_cur_perLibrary <- cbind(bulk_fil_b6, bulk_fil_cast)
  
  # sce_par_cur_perLibrary = data.frame(
  #   "Sample1_Reference" = rowSums(counts(data[,data$Library == "Sample1"])),
  #   "Sample2_Reference" = rowSums(counts(data[,data$Library == "Sample2"])),
  #   "Sample3_Alternative" = rowSums(counts(data[,data$Library == "Sample3"])),
  #   "Sample4_Alternative" = rowSums(counts(data[,data$Library == "Sample4"]))
  # )
  # 
  # sce_fil_cur_perLibrary = data.frame(
  #   "Sample5_Reference" = rowSums(counts_reference(data[,data$Library == "Sample5"])),
  #   "Sample6_Reference" = rowSums(counts_reference(data[,data$Library == "Sample6"])),
  #   "Sample5_Alternative" = rowSums(counts_alternative(data[,data$Library == "Sample5"])),
  #   "Sample6_Alternative" = rowSums(counts_alternative(data[,data$Library == "Sample6"]))
  # )
  # 
  genes.both <- intersect(rownames(sce_par_cur_perLibrary),
                          rownames(sce_fil_cur_perLibrary))
  
  sce_par_cur_perLibrary <- sce_par_cur_perLibrary[genes.both,]
  sce_fil_cur_perLibrary <- sce_fil_cur_perLibrary[genes.both,]
  
  # disps.par.b6 <- estimateDisp(sce_par_cur_perLibrary[,1:2])$tagwise.dispersion
  # names(disps.par.b6) <- rownames(sce_par_cur_perLibrary)
  # disps.par.cast <- estimateDisp(sce_par_cur_perLibrary[,3:4])$tagwise.dispersion
  # names(disps.par.cast) <- rownames(sce_par_cur_perLibrary)
  
  sce_par_cur_perLibrary_norm <- normalize_sf(sce_par_cur_perLibrary)
  sce_fil_cur_perLibrary_norm <- normalize_sf(sce_fil_cur_perLibrary)
  
  return(list(sce_par_cur_perLibrary_norm, sce_fil_cur_perLibrary_norm))
}

bulk_sce <- aggregateAcrossCells(data, data$IndividualSamples, use.assay.type = c("counts", "counts_reference", "counts_alternative"))
bulk_sce$sizeFactor <- calculateSumFactors(bulk_sce)
bulk_sce <- logNormCounts(bulk_sce)
reducedDims(bulk_sce)[["PCA"]] <- calculatePCA(bulk_sce)
plotPCA(bulk_sce, colour_by = "Dataset", text_by = "Species")

bulk_all <- make_bulk_celltype(data)
bulk_sc <- make_bulk_celltype(data[,data$CellType == "SC"])
bulk_rs <- make_bulk_celltype(data[,data$CellType == "RS"])
bulk_es <- make_bulk_celltype(data[,data$CellType == "ES"])

# downsample parentals 
# downsample_parentals <- function(x){
#   sfs <- estimateSizeFactorsForMatrix(cbind(x$, par_cast))
#   sfs <- min(sfs) / sfs
#   
#   par_b6_norm <- downsampleMatrix(par_b6, prop = sfs[1:ncol(par_b6)], bycol = T)
#   par_cast_norm <- downsampleMatrix(par_cast, prop = sfs[(ncol(par_cast)+1):(ncol(par_cast) + ncol(par_cast))], bycol = T)
#   
# }

bulk_all <- list(
  bulk_all[[1]][rowMeans(bulk_all[[1]]) > 50 & rowMeans(bulk_all[[2]]) > 50, ],
  bulk_all[[2]][rowMeans(bulk_all[[1]]) > 50 & rowMeans(bulk_all[[2]]) > 50, ]
)

# Check samples with PCA at this stage
pca_parental <- prcomp(t(log(bulk_all[[1]] + 1)))
plot(pca_parental$sdev / sum(pca_parental$sdev) * 100)

data.frame(
  PC1 = pca_parental$x[,1], 
  PC2 = pca_parental$x[,2], 
  Exp = gsub("_Sample[0-9]_(Reference|Alternative)", "", colnames(bulk_all[[1]]))
) %>%
  ggplot(aes(x = PC1, y = PC2)) + geom_point(aes(col = Exp))

pca_parental <- prcomp(t(log(bulk_all[[2]] + 1)))
plot(pca_parental$sdev / sum(pca_parental$sdev) * 100)

data.frame(
  PC1 = pca_parental$x[,1], 
  PC2 = pca_parental$x[,2], 
  Exp = gsub("_Sample[0-9]_(Reference|Alternative)", "", colnames(bulk_all[[1]]))
) %>%
  ggplot(aes(x = PC1, y = PC2)) + geom_point(aes(col = Exp))

plot_fc_df <- data.frame(
  FC_parental = log(rowMeans(bulk_all[[1]][,1:6]) + 1) - log(rowMeans(bulk_all[[1]][,7:12]) + 1),
  FC_filial = log(rowMeans(bulk_all[[2]][,1:6]) + 1) - log(rowMeans(bulk_all[[2]][,7:12]) + 1)
)

chromosome_annotation <- annotate_chromosome_list(rownames(bulk_all[[1]]))
is_maternal <- chromosome_annotation$chromosome_name %in% c("MT", "X")

plot_fc_df$is_mito <- is_maternal

bulk_correlation <- cor(plot_fc_df$FC_parental[!plot_fc_df$is_mito ], plot_fc_df$FC_filial[!plot_fc_df$is_mito ])
bulk_correlation <- round(bulk_correlation, digits = 2)

plot_fc_df %>%
  dplyr::filter(!is_maternal) %>%
  ggplot(aes(x = FC_filial, y = FC_parental)) + geom_point(size = 0.5, col = "grey")  + theme_paper() + 
  geom_abline() + geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed")

convert_fc_ase <- function(fc){
  exp(-fc) / (exp(-fc) + 1)
}

plot_fc_df %>%
  dplyr::filter(!is_maternal) %>%
  ggplot(aes(x = convert_fc_ase(FC_filial), y = convert_fc_ase(FC_parental))) + geom_point(size = 0.5, col = "grey")  + theme_paper() + 
  geom_abline() + geom_hline(yintercept = 0.5, linetype = "dashed") + geom_vline(xintercept = 0.5, linetype = "dashed")

########################################################################################################################################################
# define functions for model selection
########################################################################################################################################################
library(emdbook)

LL_conserved <- function(mu1, disp, a){
  
  # conserved: mu1 = mu2, a = b
  
  R_p_b6 = dnbinom(data_gene$p_b6, mu = mu1, size = disp)
  R_p_cast = dnbinom(data_gene$p_cast, mu = mu1, size = disp)
  
  R_f1 = unlist(
    lapply(1:length(data_gene$f_b6), function(i){
      n_a = data_gene$f_b6[i]
      n_b = data_gene$f_cast[i]
      n_i = n_a + n_b
      
      emdbook::dbetabinom(n_a, size = n_i, shape1 = a, shape2 = a)
    })
  )
  
  -( sum(log(R_p_b6)) + sum(log(R_p_cast)) + sum(log(R_f1)))
}
LL_cis <- function(mu1, mu2, disp, a){
  
  # conserved: mu1 = mu2, a = b
  
  R_p_b6 = dnbinom(data_gene$p_b6, mu = mu1, size = disp)
  R_p_cast = dnbinom(data_gene$p_cast, mu = mu2, size = disp)
  
  R_f1 = unlist(
    lapply(1:length(data_gene$f_b6), function(i){
      n_a = data_gene$f_b6[i]
      n_b = data_gene$f_cast[i]
      n_i = n_a + n_b
      
      emdbook::dbetabinom(n_a, size = n_i, shape1 = a, shape2 = a * mu2 / mu1)
    })
  )
  
  -( sum(log(R_p_b6)) + sum(log(R_p_cast)) + sum(log(R_f1)))
}
LL_trans <- function(mu1, mu2, disp, a){
  
  # conserved: mu1 = mu2, a = b
  
  R_p_b6 = dnbinom(data_gene$p_b6, mu = mu1, size = disp)
  R_p_cast = dnbinom(data_gene$p_cast, mu = mu2, size = disp)
  
  R_f1 = unlist(
    lapply(1:length(data_gene$f_b6), function(i){
      n_a = data_gene$f_b6[i]
      n_b = data_gene$f_cast[i]
      n_i = n_a + n_b
      
      emdbook::dbetabinom(n_a, size = n_i, shape1 = a, shape2 = a)
    })
  )
  
  -( sum(log(R_p_b6)) + sum(log(R_p_cast)) + sum(log(R_f1)))
}
LL_cis_trans <- function(mu1,mu2, disp, a, b){
  
  # conserved: mu1 = mu2, a = b
  
  R_p_b6 = dnbinom(data_gene$p_b6, mu = mu1, size = disp)
  R_p_cast = dnbinom(data_gene$p_cast, mu = mu2, size = disp)
  
  R_f1 = unlist(
    lapply(1:length(data_gene$f_b6), function(i){
      n_a = data_gene$f_b6[i]
      n_b = data_gene$f_cast[i]
      n_i = n_a + n_b
      
      emdbook::dbetabinom(n_a, size = n_i, shape1 = a, shape2 = b)
    })
  )
  
  -( sum(log(R_p_b6)) + sum(log(R_p_cast)) + sum(log(R_f1)))
}

mom_bb <- function(x1, x2){
  xx = x1
  ns = x1 + x2
  k = 1
  m_1 = (sum(ns*xx**k))/(sum(ns)) 
  
  k = 2
  m_2 = (sum(ns*xx**k))/(sum(ns))
  
  n = mean(ns)
  
  alpha = (n*m_1-m_2)/(n*(m_2/m_1-m_1-1)+m_1)
  beta = (n-m_1)*(n-m_2/m_1)/(n*(m_2/m_1-m_1-1)+m_1)
  return(c(alpha, beta))
}
model_selection_per_gene <- function(gene, curr_disp){
  
  print(gene)
  
  # get the relevant data
  # data_gene <<- list(
  #   p_b6 = as.numeric((t(t(par_b6[gene, ]) / sfs_total[1:6]))),
  #   p_cast = as.numeric((t(t(par_cast[gene, ]) / sfs_total[7:12]))),
  #   f_b6 = as.numeric(f1_b6_init[gene, ]),
  #   f_cast = as.numeric(f1_cast_init[gene, ])
  # )
  
  data_gene <<- list(
    p_b6 = as.numeric(sce_par_cur_perLibrary_norm[gene ,1:6]),
    p_cast = as.numeric(sce_par_cur_perLibrary_norm[gene ,7:12]),
    f_b6 = as.numeric(sce_fil_cur_perLibrary_norm[gene, 1:6]),
    f_cast = as.numeric(sce_fil_cur_perLibrary_norm[gene, 7:12])
  )
  
  # data_gene <<- list(
  #   p_b6 = as.numeric(par_b6[gene, ]),
  #   p_cast = as.numeric(par_cast[gene, ]),
  #   f_b6 = as.numeric(f1_b6_init[gene, ]),
  #   f_cast = as.numeric(f1_cast_init[gene, ])
  # )
  
  # data_gene <<- list(
  #   p_b6 = round(rnorm(6, sd = 10) + 200),
  #   p_cast =  round(rnorm(6, sd = 10) + 200),
  #   f_b6 =  round(rnorm(12, sd = 10) + 200),
  #   f_cast =  round(rnorm(12, sd = 10) + 200)
  # )
  
  # 
  #print(c(log2(mean(data_gene$p_b6) / mean(data_gene$p_cast)), log2(mean(data_gene$f_b6) / mean(data_gene$f_cast))))
  
  #plot(unlist(data_gene))
  
  # set starting parameters as empirical means
  
  mu1_start = mean(data_gene$p_b6)
  mu2_start = mean(data_gene$p_cast)
  
  # what to set a_start, b_start as ???
  # 
  
  mom_estimate = mom_bb(data_gene$f_b6, data_gene$f_cast)
  
  a_start = max(mom_estimate[[1]], 0.1)
  b_start = max(mom_estimate[[2]], 0.1)
  
  #a_start = 1
  #b_start = 1
  
  # get disp parameters
  
  # curr_disp1 = 1 / disps.par.b6[gene]
  # curr_disp2 = 1 / disps.par.cast[gene]
  
  curr_disp = 1 / disps.par.joint[gene]
  
  #curr_disp1 = 1000
  #curr_disp2 = 1000
  
  # print("Empirical means:\n")
  # print(mu1_start)
  # print(mu2_start)
  # 
  # print("EdgeR dispersions means:\n")
  # print(curr_disp1)
  # print(curr_disp2)
  
  tryCatch({
    suppressWarnings({
      #print("con")
      
      mle.conserved <- mle(LL_conserved, start = list(mu1 = mu1_start, a = a_start), fixed = list(disp = curr_disp))
      mle.cis <- mle(LL_cis, start = list(mu1 = mu1_start, mu2 = mu2_start, a = a_start), fixed = list(disp = curr_disp))
      mle.trans <- mle(LL_trans, start = list(mu1 = mu1_start, mu2 = mu2_start, a = a_start), fixed = list(disp = curr_disp))
      mle.trans.cis <- mle(LL_cis_trans, start = list(mu1 = mu1_start, mu2 = mu2_start, a = a_start, b = b_start), fixed = list(disp = curr_disp))
      
      mle.conserved@nobs <- length(unlist(data_gene[1:3]))
      mle.cis@nobs <- length(unlist(data_gene[1:3]))
      mle.trans@nobs <- length(unlist(data_gene[1:3]))
      mle.trans.cis@nobs <- length(unlist(data_gene[1:3]))
      
      # mle.conserved <- mle(LL_conserved, start = list(mu1 = mu1_start, a = a_start, disp = 1))
      # mle.conserved@nobs <- length(unlist(data_gene[1:3]))
      # #print("cis")
      # 
      # mle.cis <- mle(LL_cis, start = list(mu1 = mu1_start, mu2 = mu2_start, a = a_start, disp = 1))
      # mle.cis@nobs <- length(unlist(data_gene[1:3]))
      # #print("trans")
      # 
      # mle.trans <- mle(LL_trans, start = list(mu1 = mu1_start, mu2 = mu2_start, a = a_start, disp = 1))
      # mle.trans@nobs <- length(unlist(data_gene[1:3]))
      # #print("cis-trans")
      # 
      # mle.trans.cis <- mle(LL_cis_trans, start = list(mu1 = mu1_start, mu2 = mu2_start, a = a_start, b = b_start, disp = 1))
      # mle.trans.cis@nobs <- length(unlist(data_gene[1:3]))
      
      lls <- c("conserved" = logLik(mle.conserved),
               "cis" = logLik(mle.cis),
               "trans" = logLik(mle.trans),
               "cis_trans" = logLik(mle.trans.cis))
      
      aics <- c("conserved" = BIC(mle.conserved),
                "cis" = BIC(mle.cis),
                "trans" = BIC(mle.trans),
                "cis_trans" = BIC(mle.trans.cis))
      
      # aics <- c("conserved" = AIC(mle.conserved),
      #           "cis" = AIC(mle.cis),
      #           "trans" = AIC(mle.trans),
      #           "cis_trans" = AIC(mle.trans.cis))
      aics <- round(aics, digits = 4)
      print(aics - min(aics))
      category_ll <- names(lls[lls == max(lls)])
      category <- names(aics[aics == min(aics)])
    })
  },
  error=function(cond) {
    print("Error!!!")
    lls <<- rep(NA, 4)
    aics <<- rep(NA, 4)
    category <<- "NotConverged"
    category_ll <<- "NotConverged"
    return(aics)
  })
  
  expr_p_b6 <- mean(data_gene$p_b6)
  expr_p_cast <- mean(data_gene$p_cast)
  expr_f_b6 <- mean(data_gene$f_b6)
  expr_f_cast <- mean(data_gene$f_cast)
  fc_parental <- round(log(expr_p_b6) - log(expr_p_cast), digits = 4)
  fc_filial <- round(log(expr_f_b6) - log(expr_f_cast), digits = 4)
  
  return(c("Gene" = gene, 
           "Expr_Par_B6" = expr_p_b6,
           "Expr_Par_Cast" = expr_p_cast,
           "Expr_Fil_B6" = expr_f_b6,
           "Expr_Fil_Cast" = expr_f_cast,
           "FC_parental" = fc_parental,
           "FC_filial" = fc_filial,
           "LL" = lls["conserved"],
           "LL" = lls["cis"],
           "LL" = lls["trans"],
           "LL" = lls["cis_trans"],
           "AIC" = aics["conserved"], 
           "AIC" = aics["cis"], 
           "AIC" = aics["trans"], 
           "AIC" = aics["cis_trans"],
           "Category" = category, 
           "Category_LL" = category_ll))
}

########################################################################################################################################################
# Run MLE across cells
########################################################################################################################################################
sce_par_cur_perLibrary_norm <- bulk_all[[1]]
sce_fil_cur_perLibrary_norm <- bulk_all[[2]]

sce_par_cur_perLibrary_norm <- sce_par_cur_perLibrary_norm[!is_maternal, ]
sce_fil_cur_perLibrary_norm <- sce_fil_cur_perLibrary_norm[!is_maternal, ]

genes.test <- rownames(sce_par_cur_perLibrary_norm)

library(DESeq2)
deseqset.here <- DESeqDataSetFromMatrix(cbind(sce_par_cur_perLibrary_norm, sce_fil_cur_perLibrary_norm), 
                                        colData = DataFrame(Samples = colnames(cbind(sce_par_cur_perLibrary_norm, sce_fil_cur_perLibrary_norm))), design = ~1)
deseqset.here <- estimateSizeFactors(deseqset.here)
deseqset.here <- estimateDispersions(deseqset.here)
disps.par.joint <- rowData(deseqset.here)$dispersion
names(disps.par.joint) <- rownames(sce_par_cur_perLibrary_norm)

model.selection.results.bulk <- data.frame(
  do.call("rbind", lapply(genes.test, function(x){
    model_selection_per_gene(x, disps.par.joint)
  }
)))

model.selection.results.bulk.save <- model.selection.results.bulk
model.selection.results.bulk <- model.selection.results.bulk.save

model.selection.results.bulk[,2:15] <- apply(model.selection.results.bulk[,2:15], 2, as.numeric)
model.selection.results.bulk$chromosome <- annotate_chromosome_list(model.selection.results.bulk$Gene)$chromosome_name
model.selection.results.bulk <- model.selection.results.bulk[!is.na(model.selection.results.bulk$chromosome), ]
model.selection.results.bulk <- model.selection.results.bulk[!model.selection.results.bulk$chromosome %in% c("X", "MT"), ]
model.selection.results.bulk <- model.selection.results.bulk[model.selection.results.bulk$Category != "NotConverged", ]
model.selection.results.bulk$Category <- factor(model.selection.results.bulk$Category, 
                                                levels = c("conserved", "cis", "trans", "cis_trans"))

cf = 3
model.selection.results.bulk[abs(model.selection.results.bulk$FC_filial) > cf, ]$FC_filial <- 
  cf * sign(model.selection.results.bulk[abs(model.selection.results.bulk$FC_filial) > cf, ]$FC_filial)
model.selection.results.bulk[abs(model.selection.results.bulk$FC_parental) > cf, ]$FC_parental <- 
  cf * sign(model.selection.results.bulk[abs(model.selection.results.bulk$FC_parental) > cf, ]$FC_parental)
sigsn <- !(abs(model.selection.results.bulk$FC_filial) == cf | abs(model.selection.results.bulk$FC_parental) == cf)
model.selection.results.bulk$plot_sign <- ifelse(sigsn, "1", "2")

names(colours_cis_trans) <- tolower(names(colours_cis_trans))

apply(model.selection.results.bulk, 1, function(x){
  aic_diffs = as.numeric(x["AIC.conserved"]) - as.numeric(x[c("AIC.conserved", "AIC.cis", "AIC.trans", "AIC.cis_trans")])
  if (all(aic_diffs < 4)){
    return("conserved")
  } else {
    return(c("AIC.conserved", "AIC.cis", "AIC.trans", "AIC.cis_trans")[aic_diffs == max(aic_diffs)])
  }
}) -> CategoryNew

model.selection.results.bulk$CategoryNew <- factor(gsub("AIC\\.", "", CategoryNew), c("conserved", "cis", "trans", "cis_trans"))

scatter_p <- 
  model.selection.results.bulk %>%
  mutate(CategoryNew = factor(CategoryNew, levels = c("conserved", "cis", "cis_trans", "trans"))) %>%
  dplyr::arrange(CategoryNew) %>%
  ggplot(aes(FC_filial, FC_parental, col = CategoryNew)) + 
  geom_abline(linetype = "dashed") + 
  geom_hline(linetype = "dashed", yintercept = 0) + 
  geom_vline(linetype = "dashed", xintercept = 0) + 
  geom_point(size = 2, aes(shape = plot_sign)) + 
  xlab("log2(F1_B6/F1_CAST)") + ylab("log2(F0_B6/F0_CAST)") + 
  coord_fixed() + labs(colour="Regulatory Category") + 
  scale_color_manual(values = colours_cis_trans) + 
  scale_x_continuous(limits = c(-cf, cf), expand = c(0.01, 0.01)) + 
  scale_y_continuous(limits = c(-cf, cf), expand = c(0.01, 0.01)) + 
  theme_classic() + 
  theme(text = element_text(size=30)) +
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 1)) + 
  theme(legend.position = "none") + 
  ggtitle("All celltypes") + theme(plot.title = element_text(hjust = 0.5))
ggsave("~/Desktop/Projects/ASE_Spermatogenesis_Paper_rerun_Revisions/Plots/Figure1/Fig1_model_comparison_bulk_new_scatter.pdf", width = 8, height = 8)

distribution_p <- ggplot(model.selection.results.bulk, 
                         aes(y = 1, fill = forcats::fct_rev(CategoryNew))) + 
  geom_bar(position = "fill", col = "black") + 
  theme_paper(textsize = 30) + 
  scale_fill_manual(values = colours_cis_trans) + 
  theme(axis.line.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y = element_blank(), 
        axis.title.y=element_blank()) + 
  theme(plot.margin=unit(c(0.5, 5, 0.5, 6.5),"cm")) + 
  theme(legend.position = "none") + 
  xlab("") + 
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
ggsave("~/Desktop/Projects/ASE_Spermatogenesis_Paper_rerun_Revisions/Plots/Figure1/Fig1_model_comparison_bulk_new_barplot.pdf", width = 12, height = 2.5)

pp <- gridExtra::grid.arrange(scatter_p, distribution_p, nrow = 2, 
                              heights=c(4,1))

saveRDS(model.selection.results.bulk, "./Data/processed/model_selection_bulk.rds")
model.selection.results.bulk <- readRDS("./Data/processed/model_selection_bulk.rds")

########################################################################################################################################################
# Run MLE for spermatocytes
########################################################################################################################################################
sce_par_cur_perLibrary_norm <- bulk_sc[[1]][rownames(bulk_all[[1]]), ]
sce_fil_cur_perLibrary_norm <- bulk_sc[[2]][rownames(bulk_all[[1]]), ]

sce_par_cur_perLibrary_norm <- sce_par_cur_perLibrary_norm[!is_maternal, ]
sce_fil_cur_perLibrary_norm <- sce_fil_cur_perLibrary_norm[!is_maternal, ]

genes.test <- rownames(sce_par_cur_perLibrary_norm)

library(DESeq2)
deseqset.here <- DESeqDataSetFromMatrix(cbind(sce_par_cur_perLibrary_norm, sce_fil_cur_perLibrary_norm), 
                                        colData = DataFrame(Samples = colnames(cbind(sce_par_cur_perLibrary_norm, sce_fil_cur_perLibrary_norm))), design = ~1)
deseqset.here <- estimateSizeFactors(deseqset.here)
deseqset.here <- estimateDispersions(deseqset.here)
disps.par.joint <- rowData(deseqset.here)$dispersion
names(disps.par.joint) <- rownames(sce_par_cur_perLibrary_norm)

model.selection.results.sc <- data.frame(
  do.call("rbind", lapply(genes.test, function(x){
    model_selection_per_gene(x, disps.par.joint)
  }
  )))

model.selection.results.sc[,2:15] <- apply(model.selection.results.sc[,2:15], 2, as.numeric)
model.selection.results.sc$chromosome <- annotate_chromosome_list(model.selection.results.sc$Gene)$chromosome_name
model.selection.results.sc <- model.selection.results.sc[!is.na(model.selection.results.sc$chromosome), ]
model.selection.results.sc <- model.selection.results.sc[!model.selection.results.sc$chromosome %in% c("X", "MT"), ]
model.selection.results.sc <- model.selection.results.sc[model.selection.results.sc$Category != "NotConverged", ]
model.selection.results.sc$Category <- factor(model.selection.results.sc$Category, 
                                              levels = c("conserved", "cis", "trans", "cis_trans"))

model.selection.results.sc[abs(model.selection.results.sc$FC_filial) > cf, ]$FC_filial <- 
  cf * sign(model.selection.results.sc[abs(model.selection.results.sc$FC_filial) > cf, ]$FC_filial)
model.selection.results.sc[abs(model.selection.results.sc$FC_parental) > cf, ]$FC_parental <- 
  cf * sign(model.selection.results.sc[abs(model.selection.results.sc$FC_parental) > cf, ]$FC_parental)
sigsn <- !(abs(model.selection.results.sc$FC_filial) == cf | abs(model.selection.results.sc$FC_parental) == cf)
model.selection.results.sc$plot_sign <- ifelse(sigsn, "1", "2")

apply(model.selection.results.sc, 1, function(x){
  aic_diffs = as.numeric(x["AIC.conserved"]) - as.numeric(x[c("AIC.conserved", "AIC.cis", "AIC.trans", "AIC.cis_trans")])
  if (all(aic_diffs < 4)){
    return("conserved")
  } else {
    return(c("AIC.conserved", "AIC.cis", "AIC.trans", "AIC.cis_trans")[aic_diffs == max(aic_diffs)])
  }
}) -> CategoryNew

model.selection.results.sc$CategoryNew <- factor(gsub("AIC\\.", "", CategoryNew), c("conserved", "cis", "trans", "cis_trans"))

scatter_p <- 
  model.selection.results.sc %>%
  mutate(CategoryNew = factor(CategoryNew, levels = c("conserved", "cis", "cis_trans", "trans"))) %>%
  dplyr::arrange(CategoryNew) %>%
  ggplot(aes(FC_filial, FC_parental, col = CategoryNew)) + 
  geom_abline(linetype = "dashed") + 
  geom_hline(linetype = "dashed", yintercept = 0) + 
  geom_vline(linetype = "dashed", xintercept = 0) + 
  geom_point(size = 2, aes(shape = plot_sign)) + 
  xlab("log2(F1_B6/F1_CAST)") + ylab("log2(F0_B6/F0_CAST)") + 
  coord_fixed() + labs(colour="Regulatory Category") + 
  scale_color_manual(values = colours_cis_trans) + 
  scale_x_continuous(limits = c(-cf, cf), expand = c(0.01, 0.01)) + 
  scale_y_continuous(limits = c(-cf, cf), expand = c(0.01, 0.01)) + 
  theme_classic() + 
  theme(text = element_text(size=30)) +
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 1)) + 
  theme(legend.position = "none") + 
  ggtitle("Spermatocytes") + theme(plot.title = element_text(hjust = 0.5))
ggsave("~/Desktop/Projects/ASE_Spermatogenesis_Paper_rerun_Revisions/Plots/Figure1/Fig1_model_comparison_sc_new_scatter.pdf", width = 8, height = 8)

distribution_p <- ggplot(model.selection.results.sc, 
                         aes(y = 1, fill = forcats::fct_rev(CategoryNew))) + 
  geom_bar(position = "fill", col = "black") + 
  theme_paper(textsize = 30) + 
  scale_fill_manual(values = colours_cis_trans) + 
  theme(axis.line.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y = element_blank(), 
        axis.title.y=element_blank()) + 
  theme(plot.margin=unit(c(0.5, 5, 0.5, 6.5),"cm")) + 
  theme(legend.position = "none") + 
  xlab("") + 
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
ggsave("~/Desktop/Projects/ASE_Spermatogenesis_Paper_rerun_Revisions/Plots/Figure1/Fig1_model_comparison_sc_new_barplot.pdf", width = 12, height = 2.5)

pp <- gridExtra::grid.arrange(scatter_p, distribution_p, nrow = 2, 
                              heights=c(4,1))

saveRDS(model.selection.results.sc, "./Data/processed/model_selection_sc.rds")
model.selection.results.sc <- readRDS("./Data/processed/model_selection_sc.rds")

########################################################################################################################################################
# Run MLE for round spermatids
########################################################################################################################################################
sce_par_cur_perLibrary_norm <- bulk_rs[[1]][rownames(bulk_all[[1]]), ]
sce_fil_cur_perLibrary_norm <- bulk_rs[[2]][rownames(bulk_all[[1]]), ]

genes.test <- rownames(sce_par_cur_perLibrary_norm)

library(DESeq2)
deseqset.here <- DESeqDataSetFromMatrix(cbind(sce_par_cur_perLibrary_norm, sce_fil_cur_perLibrary_norm), 
                                        colData = DataFrame(Samples = colnames(cbind(sce_par_cur_perLibrary_norm, sce_fil_cur_perLibrary_norm))), design = ~1)
deseqset.here <- estimateSizeFactors(deseqset.here)
deseqset.here <- estimateDispersions(deseqset.here)
disps.par.joint <- rowData(deseqset.here)$dispersion
names(disps.par.joint) <- rownames(sce_par_cur_perLibrary_norm)

model.selection.results.rs <- data.frame(
  do.call("rbind", lapply(genes.test, function(x){
    model_selection_per_gene(x, disps.par.joint)
  }
  )))

model.selection.results.rs[,2:15] <- apply(model.selection.results.rs[,2:15], 2, as.numeric)
model.selection.results.rs$chromosome <- annotate_chromosome_list(model.selection.results.rs$Gene)$chromosome_name
model.selection.results.rs <- model.selection.results.rs[!is.na(model.selection.results.rs$chromosome), ]
model.selection.results.rs <- model.selection.results.rs[!model.selection.results.rs$chromosome %in% c("X", "MT"), ]
model.selection.results.rs <- model.selection.results.rs[model.selection.results.rs$Category != "NotConverged", ]
model.selection.results.rs$Category <- factor(model.selection.results.rs$Category, 
                                              levels = c("conserved", "cis", "trans", "cis_trans"))

model.selection.results.rs[abs(model.selection.results.rs$FC_filial) > cf, ]$FC_filial <- 
  cf * sign(model.selection.results.rs[abs(model.selection.results.rs$FC_filial) > cf, ]$FC_filial)
model.selection.results.rs[abs(model.selection.results.rs$FC_parental) > cf, ]$FC_parental <- 
  cf * sign(model.selection.results.rs[abs(model.selection.results.rs$FC_parental) > cf, ]$FC_parental)
sigsn <- !(abs(model.selection.results.rs$FC_filial) == cf | abs(model.selection.results.rs$FC_parental) == cf)
model.selection.results.rs$plot_sign <- ifelse(sigsn, "1", "2")

apply(model.selection.results.rs, 1, function(x){
  aic_diffs = as.numeric(x["AIC.conserved"]) - as.numeric(x[c("AIC.conserved", "AIC.cis", "AIC.trans", "AIC.cis_trans")])
  if (all(aic_diffs < 4)){
    return("conserved")
  } else {
    return(c("AIC.conserved", "AIC.cis", "AIC.trans", "AIC.cis_trans")[aic_diffs == max(aic_diffs)])
  }
}) -> CategoryNew

model.selection.results.rs$CategoryNew <- factor(gsub("AIC\\.", "", CategoryNew), c("conserved", "cis", "trans", "cis_trans"))

scatter_p <- 
  model.selection.results.rs %>%
  mutate(CategoryNew = factor(CategoryNew, levels = c("conserved", "cis", "cis_trans", "trans"))) %>%
  dplyr::arrange(CategoryNew) %>%
  ggplot(aes(FC_filial, FC_parental, col = CategoryNew)) + 
  geom_abline(linetype = "dashed") + 
  geom_hline(linetype = "dashed", yintercept = 0) + 
  geom_vline(linetype = "dashed", xintercept = 0) + 
  geom_point(size = 2, aes(shape = plot_sign)) + 
  xlab("log2(F1_B6/F1_CAST)") + ylab("log2(F0_B6/F0_CAST)") + 
  coord_fixed() + labs(colour="Regulatory Category") + 
  scale_color_manual(values = colours_cis_trans) + 
  scale_x_continuous(limits = c(-cf, cf), expand = c(0.01, 0.01)) + 
  scale_y_continuous(limits = c(-cf, cf), expand = c(0.01, 0.01)) + 
  theme_classic() + 
  theme(text = element_text(size=30)) +
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 1)) + 
  theme(legend.position = "none") + 
  ggtitle("Round \n Spermatids") + theme(plot.title = element_text(hjust = 0.5))
ggsave("~/Desktop/Projects/ASE_Spermatogenesis_Paper_rerun_Revisions/Plots/Figure1/Fig1_model_comparison_rs_new_scatter.pdf", width = 8, height = 8)

distribution_p <- ggplot(model.selection.results.rs, 
                         aes(y = 1, fill = forcats::fct_rev(CategoryNew))) + 
  geom_bar(position = "fill", col = "black") + 
  theme_paper(textsize = 30) + 
  scale_fill_manual(values = colours_cis_trans) + 
  theme(axis.line.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y = element_blank(), 
        axis.title.y=element_blank()) + 
  theme(plot.margin=unit(c(0.5, 5, 0.5, 6.5),"cm")) + 
  theme(legend.position = "none") + 
  xlab("") + 
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
ggsave("~/Desktop/Projects/ASE_Spermatogenesis_Paper_rerun_Revisions/Plots/Figure1/Fig1_model_comparison_rs_new_barplot.pdf", width = 12, height = 2.5)

pp <- gridExtra::grid.arrange(scatter_p, distribution_p, nrow = 2, 
                              heights=c(4,1))
ggsave("~/Desktop/Projects/ASE_Spermatogenesis_Paper_rerun_Revisions/Plots/Figure1/Fig1_model_comparison_rs_new.pdf")

saveRDS(model.selection.results.rs, "./Data/processed/model_selection_rs.rds")
model.selection.results.rs <- readRDS("./Data/processed/model_selection_rs.rds")

########################################################################################################################################################
# Run MLE for elongatic spermatids
########################################################################################################################################################
sce_par_cur_perLibrary_norm <- bulk_es[[1]][rownames(bulk_all[[1]]), ]
sce_fil_cur_perLibrary_norm <- bulk_es[[2]][rownames(bulk_all[[1]]), ]

genes.test <- rownames(sce_par_cur_perLibrary_norm)

library(DESeq2)
deseqset.here <- DESeqDataSetFromMatrix(cbind(sce_par_cur_perLibrary_norm, sce_fil_cur_perLibrary_norm), 
                                        colData = DataFrame(Samples = colnames(cbind(sce_par_cur_perLibrary_norm, sce_fil_cur_perLibrary_norm))), design = ~1)
deseqset.here <- estimateSizeFactors(deseqset.here)
deseqset.here <- estimateDispersions(deseqset.here)
disps.par.joint <- rowData(deseqset.here)$dispersion
names(disps.par.joint) <- rownames(sce_par_cur_perLibrary_norm)
model.selection.results.es <- data.frame(
  do.call("rbind", lapply(genes.test, function(x){
    model_selection_per_gene(x, disps.par.joint)
  }
)))

model.selection.results.es[,2:15] <- apply(model.selection.results.es[,2:15], 2, as.numeric)
model.selection.results.es$chromosome <- annotate_chromosome_list(model.selection.results.es$Gene)$chromosome_name
model.selection.results.es <- model.selection.results.es[!is.na(model.selection.results.es$chromosome), ]
model.selection.results.es <- model.selection.results.es[!model.selection.results.es$chromosome %in% c("X", "MT"), ]
model.selection.results.es <- model.selection.results.es[model.selection.results.es$Category != "NotConverged", ]
model.selection.results.es$Category <- factor(model.selection.results.es$Category, 
                                              levels = c("conserved", "cis", "trans", "cis_trans"))

model.selection.results.es[abs(model.selection.results.es$FC_filial) > cf, ]$FC_filial <- 
  cf * sign(model.selection.results.es[abs(model.selection.results.es$FC_filial) > cf, ]$FC_filial)
model.selection.results.es[abs(model.selection.results.es$FC_parental) > cf, ]$FC_parental <- 
  cf * sign(model.selection.results.es[abs(model.selection.results.es$FC_parental) > cf, ]$FC_parental)
sigsn <- !(abs(model.selection.results.es$FC_filial) == cf | abs(model.selection.results.es$FC_parental) == cf)
model.selection.results.es$plot_sign <- ifelse(sigsn, "1", "2")

apply(model.selection.results.es, 1, function(x){
  aic_diffs = as.numeric(x["AIC.conserved"]) - as.numeric(x[c("AIC.conserved", "AIC.cis", "AIC.trans", "AIC.cis_trans")])
  if (all(aic_diffs < 4)){
    return("conserved")
  } else {
    return(c("AIC.conserved", "AIC.cis", "AIC.trans", "AIC.cis_trans")[aic_diffs == max(aic_diffs)])
  }
}) -> CategoryNew

model.selection.results.es$CategoryNew <- factor(gsub("AIC\\.", "", CategoryNew), c("conserved", "cis", "trans", "cis_trans"))

scatter_p <- 
  model.selection.results.es %>%
  mutate(CategoryNew = factor(CategoryNew, levels = c("conserved", "cis", "cis_trans", "trans"))) %>%
  dplyr::arrange(CategoryNew) %>%
  ggplot(aes(FC_filial, FC_parental, col = CategoryNew)) + 
  geom_abline(linetype = "dashed") + 
  geom_hline(linetype = "dashed", yintercept = 0) + 
  geom_vline(linetype = "dashed", xintercept = 0) + 
  geom_point(size = 2, aes(shape = plot_sign)) + 
  xlab("log2(F1_B6/F1_CAST)") + ylab("log2(F0_B6/F0_CAST)") + 
  coord_fixed() + labs(colour="Regulatory Category") + 
  scale_color_manual(values = colours_cis_trans) + 
  scale_x_continuous(limits = c(-cf, cf), expand = c(0.01, 0.01)) + 
  scale_y_continuous(limits = c(-cf, cf), expand = c(0.01, 0.01)) + 
  theme_classic() + 
  theme(text = element_text(size=30)) +
  theme(panel.background = element_rect(fill = "white", colour = "black", size = 1)) + 
  theme(legend.position = "none") + 
  ggtitle("Elongating \n Spermatids") + theme(plot.title = element_text(hjust = 0.5))
ggsave("~/Desktop/Projects/ASE_Spermatogenesis_Paper_rerun_Revisions/Plots/Figure1/Fig1_model_comparison_es_new_scatter.pdf", width = 8, height = 8)

distribution_p <- ggplot(model.selection.results.es, 
                         aes(y = 1, fill = forcats::fct_rev(CategoryNew))) + 
  geom_bar(position = "fill", col = "black") + 
  theme_paper(textsize = 30) + 
  scale_fill_manual(values = colours_cis_trans) + 
  theme(axis.line.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y = element_blank(), 
        axis.title.y=element_blank()) + 
  theme(plot.margin=unit(c(0.5, 5, 0.5, 6.5),"cm")) + 
  theme(legend.position = "none") + 
  xlab("") + 
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
ggsave("~/Desktop/Projects/ASE_Spermatogenesis_Paper_rerun_Revisions/Plots/Figure1/Fig1_model_comparison_es_new_barplot.pdf", width = 12, height = 2.5)

pp <- gridExtra::grid.arrange(scatter_p, distribution_p, nrow = 2, 
                              heights=c(4,1))
ggsave("~/Desktop/Projects/ASE_Spermatogenesis_Paper_rerun_Revisions/Plots/Figure1/Fig1_model_comparison_es_new.pdf")

saveRDS(model.selection.results.es, "./Data/processed/model_selection_es.rds")
model.selection.results.es <- readRDS("./Data/processed/model_selection_es.rds")



