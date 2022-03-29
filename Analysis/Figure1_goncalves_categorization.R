library(scran)
library(scater)

source("./Scripts/General/auxiliary.R")
source("./Scripts/General/reuse_functions.R")

setwd("~/Desktop/Projects/ASE_Spermatogenesis_Paper_rerun/")

data <- readRDS("./Data/processed/sce_merged_new.rds")
colnames(data) <- paste0(data$Library, "_", colnames(data))

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
  
  sce_par_cur_perLibrary = data.frame(
    "Sample1_Reference" = rowSums(counts(data[,data$Library == "Sample1"])),
    "Sample2_Reference" = rowSums(counts(data[,data$Library == "Sample2"])),
    "Sample3_Reference" = rowSums(counts(data[,data$Library == "Sample3"])),
    "Sample4_Reference" = rowSums(counts(data[,data$Library == "Sample4"]))
  )
  
  sce_fil_cur_perLibrary = data.frame(
    "Sample5_Reference" = rowSums(counts_reference(data[,data$Library == "Sample5"])),
    "Sample6_Reference" = rowSums(counts_reference(data[,data$Library == "Sample6"])),
    "Sample5_Alternative" = rowSums(counts_alternative(data[,data$Library == "Sample5"])),
    "Sample6_Alternative" = rowSums(counts_alternative(data[,data$Library == "Sample6"]))
  )
  
  genes.both <- intersect(rownames(sce_par_cur_perLibrary),
                          rownames(sce_fil_cur_perLibrary))
  
  sce_par_cur_perLibrary <- sce_par_cur_perLibrary[genes.both,]
  sce_fil_cur_perLibrary <- sce_fil_cur_perLibrary[genes.both,]
  
  #return(list(sce_par_cur_perLibrary, sce_fil_cur_perLibrary))
  
  disps.par.b6 <- estimateDisp(sce_par_cur_perLibrary[,1:2])$tagwise.dispersion
  names(disps.par.b6) <- rownames(sce_par_cur_perLibrary)
  disps.par.cast <- estimateDisp(sce_par_cur_perLibrary[,3:4])$tagwise.dispersion
  names(disps.par.cast) <- rownames(sce_par_cur_perLibrary)
  
  sce_par_cur_perLibrary_norm <- normalize_sf(sce_par_cur_perLibrary)
  sce_fil_cur_perLibrary_norm <- normalize_sf(sce_fil_cur_perLibrary)
  
  return(list(sce_par_cur_perLibrary_norm, sce_fil_cur_perLibrary_norm))
}

bulk_all <- make_bulk_celltype(data)
bulk_sc <- make_bulk_celltype(data[,data$CellType == "SC"])
bulk_rs <- make_bulk_celltype(data[,data$CellType == "RS"])
bulk_es <- make_bulk_celltype(data[,data$CellType == "ES"])

bulk_all <- list(
  bulk_all[[1]][rowMeans(bulk_all[[1]]) > 100 & rowMeans(bulk_all[[2]]) > 100, ],
  bulk_all[[2]][rowMeans(bulk_all[[1]]) > 100 & rowMeans(bulk_all[[2]]) > 100, ]
)

plot_fc_df <- data.frame(
  FC_parental = log(rowMeans(bulk_all[[1]][,1:2]) + 1) - log(rowMeans(bulk_all[[1]][,3:4]) + 1),
  FC_filial = log(rowMeans(bulk_all[[2]][,1:2]) + 1) - log(rowMeans(bulk_all[[2]][,3:4]) + 1)
)

chromosome_annotation <- annotate_chromosome_list(rownames(bulk_all[[1]]))
is_maternal <- chromosome_annotation$chromosome_name %in% c("MT", "X")

plot_fc_df$is_mito <- is_maternal

bulk_correlation <- cor(plot_fc_df$FC_parental[!plot_fc_df$is_mito ], plot_fc_df$FC_filial[!plot_fc_df$is_mito ])
bulk_correlation <- round(bulk_correlation, digits = 2)

########################################################################################################################################################
# define functions for model selection
########################################################################################################################################################
library(emdbook)

LL_conserved <- function(mu1, disp1, disp2, a){
  
  # conserved: mu1 = mu2, a = b
  
  R_p_b6 = dnbinom(data_gene$p_b6, mu = mu1, size = disp1)
  R_p_cast = dnbinom(data_gene$p_cast, mu = mu1, size = disp2)
  
  R_f1 = sum(
    unlist(
      lapply(1:length(data_gene$f_b6), function(i){
        n_a = data_gene$f_b6[i]
        n_b = data_gene$f_cast[i]
        n_i = n_a + n_b
        
        dbetabinom(n_a, size = n_i, shape1 = a, shape2 = a)
      })
    )
  )
  
  -( sum(log(R_p_b6)) + sum(log(R_p_cast)) + log(R_f1))
}
LL_cis <- function(mu1, mu2, disp1, disp2, a){
  
  # conserved: mu1 = mu2, a = b
  
  R_p_b6 = dnbinom(data_gene$p_b6, mu = mu1, size = disp1)
  R_p_cast = dnbinom(data_gene$p_cast, mu = mu2, size = disp2)
  
  R_f1 = sum(
    unlist(
      lapply(1:length(data_gene$f_b6), function(i){
        n_a = data_gene$f_b6[i]
        n_b = data_gene$f_cast[i]
        n_i = n_a + n_b
        
        dbetabinom(n_a, size = n_i, shape1 = a, shape2 = a * mu2/mu1)
      })
    )
  )
  
  -( sum(log(R_p_b6)) + sum(log(R_p_cast)) + log(R_f1))
}
LL_trans <- function(mu1,mu2, disp1, disp2, a){
  
  # conserved: mu1 = mu2, a = b
  
  R_p_b6 = dnbinom(data_gene$p_b6, mu = mu1, size = disp1)
  R_p_cast = dnbinom(data_gene$p_cast, mu = mu2, size = disp2)
  
  R_f1 = sum(
    unlist(
      lapply(1:length(data_gene$f_b6), function(i){
        n_a = data_gene$f_b6[i]
        n_b = data_gene$f_cast[i]
        n_i = n_a + n_b
        
        dbetabinom(n_a, size = n_i, shape1 = a, shape2 = a)
      })
    )
  )
  
  -( sum(log(R_p_b6)) + sum(log(R_p_cast)) + log(R_f1))
}
LL_cis_trans <- function(mu1,mu2, disp1, disp2, a, b){
  
  # conserved: mu1 = mu2, a = b
  
  R_p_b6 = dnbinom(data_gene$p_b6, mu = mu1, size = disp1)
  R_p_cast = dnbinom(data_gene$p_cast, mu = mu2, size = disp2)
  
  R_f1 = sum(
    unlist(
      lapply(1:length(data_gene$f_b6), function(i){
        n_a = data_gene$f_b6[i]
        n_b = data_gene$f_cast[i]
        n_i = n_a + n_b
        
        dbetabinom(n_a, size = n_i, shape1 = a, shape2 = b)
      })
    )
  )
  
  -( sum(log(R_p_b6)) + sum(log(R_p_cast)) + log(R_f1))
}

model_selection_per_gene <- function(gene, curr_disp1, curr_disp2){
  
  print(gene)
  
  # get the relevant data
  data_gene <<- list(
    p_b6 = as.numeric(sce_par_cur_perLibrary_norm[gene ,1:2]),
    p_cast = as.numeric(sce_par_cur_perLibrary_norm[gene ,3:4]),
    f_b6 = as.numeric(sce_fil_cur_perLibrary_norm[gene, 1:2]),
    f_cast = as.numeric(sce_fil_cur_perLibrary_norm[gene, 3:4])
  )
  
  # set starting parameters as empirical means
  
  mu1_start = mean(data_gene$p_b6)
  mu2_start = mean(data_gene$p_cast)
  
  # what to set a_start, b_start as ???
  
  a_start = (mu1_start + mu2_start) / 2
  b_start = a_start
  
  # get disp parameters
  
  curr_disp1 = 1 / curr_disp1[gene]
  curr_disp2 = 1 / curr_disp2[gene]
  
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
      mle.conserved <- mle(LL_conserved, start = list(mu1 = mu1_start, a = a_start),
                           fixed = list(disp1 = curr_disp1, disp2 = curr_disp2))
      #print("cis")
      
      mle.cis <- mle(LL_cis, start = list(mu1 = mu1_start, mu2 = mu2_start, a = a_start),
                     fixed = list(disp1 = curr_disp1, disp2 = curr_disp2))
      #print("trans")
      
      mle.trans <- mle(LL_trans, start = list(mu1 = mu1_start, mu2 = mu2_start, a = a_start),
                       fixed = list(disp1 = curr_disp1, disp2 = curr_disp2))
      #print("cis-trans")
      
      mle.trans.cis <- mle(LL_cis_trans, start = list(mu1 = mu1_start, mu2 = mu2_start, a = a_start, b = b_start),
                           fixed = list(disp1 = curr_disp1, disp2 = curr_disp2))
      
      aics <- c("conserved" = AIC(mle.conserved),
                "cis" = AIC(mle.cis),
                "trans" = AIC(mle.trans),
                "cis_trans" = AIC(mle.trans.cis))
      lls <- c("conserved" = mle.conserved@min, 
               "cis" = mle.cis@min, 
               "trans" = mle.trans@min, 
               "cis_trans" = mle.trans.cis@min)
      lls <- round(lls, digits = 4)
      print(lls)
      aics <- round(aics, digits = 4)
      category <- names(aics[aics == min(aics)])
    })
  },
  error=function(cond) {
    print("Error!!!")
    aics <<- rep(NA, 4)
    lls <<- rep(NA, 4)
    category <<- "NotConverged"
    return(c(aics, lls))
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
           "AIC" = aics["conserved"], 
           "AIC" = aics["cis"], 
           "AIC" = aics["trans"], 
           "AIC" = aics["cis_trans"],
           "LL" = lls["conserved"], 
           "LL" = lls["cis"], 
           "LL" = lls["trans"], 
           "LL" = lls["cis_trans"],
           "Category" = category))
}


########################################################################################################################################################
# Run MLE across cells
########################################################################################################################################################sce_par_cur_perLibrary_norm <- bulk_all[[1]]
sce_par_cur_perLibrary_norm <- bulk_all[[1]]
sce_fil_cur_perLibrary_norm <- bulk_all[[2]]

sce_par_cur_perLibrary_norm <- sce_par_cur_perLibrary_norm[!is_maternal, ]
sce_fil_cur_perLibrary_norm <- sce_fil_cur_perLibrary_norm[!is_maternal, ]

genes.test <- rownames(sce_par_cur_perLibrary_norm)

disps.par.b6 <- estimateDisp(sce_par_cur_perLibrary_norm[,1:2])$tagwise.dispersion
names(disps.par.b6) <- rownames(sce_par_cur_perLibrary_norm)
disps.par.cast <- estimateDisp(sce_par_cur_perLibrary_norm[,3:4])$tagwise.dispersion
names(disps.par.cast) <- rownames(sce_par_cur_perLibrary_norm)

model.selection.results.bulk <- data.frame(
  do.call("rbind", lapply(genes.test, function(x){
    model_selection_per_gene(x, disps.par.b6, disps.par.cast)
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

saveRDS(model.selection.results.bulk, "./Data/processed/model_selection_bulk.rds")

########################################################################################################################################################
# Run MLE for spermatocytes
########################################################################################################################################################
sce_par_cur_perLibrary_norm <- bulk_sc[[1]][rownames(bulk_all[[1]]), ]
sce_fil_cur_perLibrary_norm <- bulk_sc[[2]][rownames(bulk_all[[1]]), ]

genes.test <- rownames(sce_par_cur_perLibrary_norm)

disps.par.b6 <- estimateDisp(sce_par_cur_perLibrary_norm[,1:2])$tagwise.dispersion
names(disps.par.b6) <- rownames(sce_par_cur_perLibrary_norm)
disps.par.cast <- estimateDisp(sce_par_cur_perLibrary_norm[,3:4])$tagwise.dispersion
names(disps.par.cast) <- rownames(sce_par_cur_perLibrary_norm)

model.selection.results.sc <- data.frame(
  do.call("rbind", lapply(genes.test, function(x){
    model_selection_per_gene(x, disps.par.b6, disps.par.cast)
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

saveRDS(model.selection.results.sc, "./Data/processed/model_selection_sc.rds")

########################################################################################################################################################
# Run MLE for round spermatids
########################################################################################################################################################
sce_par_cur_perLibrary_norm <- bulk_rs[[1]][rownames(bulk_all[[1]]), ]
sce_fil_cur_perLibrary_norm <- bulk_rs[[2]][rownames(bulk_all[[1]]), ]

genes.test <- rownames(sce_par_cur_perLibrary_norm)

disps.par.b6 <- estimateDisp(sce_par_cur_perLibrary_norm[,1:2])$tagwise.dispersion
names(disps.par.b6) <- rownames(sce_par_cur_perLibrary_norm)
disps.par.cast <- estimateDisp(sce_par_cur_perLibrary_norm[,3:4])$tagwise.dispersion
names(disps.par.cast) <- rownames(sce_par_cur_perLibrary_norm)

model.selection.results.rs <- data.frame(
  do.call("rbind", lapply(genes.test, function(x){
    model_selection_per_gene(x, disps.par.b6, disps.par.cast)
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

saveRDS(model.selection.results.rs, "./Data/processed/model_selection_rs.rds")



########################################################################################################################################################
# Run MLE for elongatic spermatids
########################################################################################################################################################

sce_par_cur_perLibrary_norm <- bulk_es[[1]][rownames(bulk_all[[1]]), ]
sce_fil_cur_perLibrary_norm <- bulk_es[[2]][rownames(bulk_all[[1]]), ]

genes.test <- rownames(sce_par_cur_perLibrary_norm)

disps.par.b6 <- estimateDisp(sce_par_cur_perLibrary_norm[,1:2])$tagwise.dispersion
names(disps.par.b6) <- rownames(sce_par_cur_perLibrary_norm)
disps.par.cast <- estimateDisp(sce_par_cur_perLibrary_norm[,3:4])$tagwise.dispersion
names(disps.par.cast) <- rownames(sce_par_cur_perLibrary_norm)

model.selection.results.es <- data.frame(
  do.call("rbind", lapply(genes.test, function(x){
    model_selection_per_gene(x, disps.par.b6, disps.par.cast)
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

saveRDS(model.selection.results.es, "./Data/processed/model_selection_es.rds")
