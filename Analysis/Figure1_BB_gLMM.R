library(VGAM)

source("./Scripts/General/auxiliary.R")
source("./Scripts/General/ase_functions.R")
source("./Scripts/General/reuse_functions.R")

# for gene: get counts
# fit VGAM:glm with and without celltype 
# LRT test

setwd("~/Desktop/PhD/Projects/ASE_Spermatogenesis_Paper_rerun//")

colour_vector <- c("SG" = "#046735",
                   "SC" = "#D31281",
                   "RS" = "#282B69",
                   "ES" = "grey60",
                   "Sertoli" = "#643F18",
                   "Leydig" = "#985D25",
                   "Immune" = "#FFDFB2")

sce <- readRDS("./Data/processed/sce_merged_new.rds")
data_f1 <- sce[,sce$Library %in% c("Sample5", "Sample6")]
sce <- data_f1

model.selection.results.sc <- readRDS("./Data/processed/model_selection_sc.rds")
model.selection.results.rs <- readRDS("./Data/processed/model_selection_rs.rds")
model.selection.results.es <- readRDS("./Data/processed/model_selection_es.rds")

celltype1 = "SC"
celltype2 = "RS"

genes_check <- intersect(model.selection.results.sc$Gene, model.selection.results.rs$Gene)

ds_ref_ct1 <- counts_reference(sce[genes_check, data_f1$CellType == celltype1])
ds_alt_ct1 <- counts_alternative(sce[genes_check, data_f1$CellType == celltype1])
ds_ref_ct2 <- counts_reference(sce[genes_check, data_f1$CellType == celltype2])
ds_alt_ct2 <- counts_alternative(sce[genes_check, data_f1$CellType == celltype2])
library_ct1 <- data_f1[,data_f1$CellType == celltype1]$Library
library_ct2 <- data_f1[,data_f1$CellType == celltype2]$Library

fit_vglm <- function(gene){
  print(gene)
  fit_df <- rbind(
    data.frame(counts_ref = ds_ref_ct1[gene, ], counts_alt = ds_alt_ct1[gene, ], 
               celltype = celltype1, library = library_ct1), 
    data.frame(counts_ref = ds_ref_ct2[gene, ], counts_alt = ds_alt_ct2[gene, ], 
               celltype = celltype2, library = library_ct2)
  )
  
  tryCatch({
    model_celltype <- VGAM::vglm(cbind(counts_ref, counts_alt) ~ celltype + library, 
                                 betabinomial, trace = FALSE, data = fit_df)
    model_constant <- VGAM::vglm(cbind(counts_ref, counts_alt) ~ 1 + library, 
                                 betabinomial, trace = FALSE, data = fit_df)
    
    return(list(gene, model_constant, model_celltype))
    
  }, error=function(cond) {
    message(paste("Not converged"))
    # Choose a return value in case of error
    model_celltype <<- NA
    model_constant <<- NA
    return(c(gene, model_constant, model_celltype))
  })
}

fit_results <- lapply(genes_check, fit_vglm)
saveRDS(fit_results, "./Data/processed/vglm_results_sc_rs.rds")

pvals <- lapply(fit_results, function(x){
  gene = x[[1]]
  model_1 = x[[2]]
  model_2 = x[[3]]
  if (is.na(model_1)){
    return(c(gene, NA))
  }
  else {
    return(c(gene, lrtest(model_1, model_2)@Body[2, 5]))
  }
})

pvals_df = data.frame(do.call("rbind", pvals))
pvals_df$X2 <- as.numeric(pvals_df$X2)
rownames(pvals_df) <- pvals_df$X1

saveRDS(pvals_df, "./Data/processed/vglm_pvals_sc_rs")

# --------------------------------------------------------------------------------------------

celltype1 = "SC"
celltype2 = "ES"

genes_check <- intersect(model.selection.results.sc$Gene, model.selection.results.es$Gene)

ds_ref_ct1 <- counts_reference(sce[genes_check, data_f1$CellType == celltype1])
ds_alt_ct1 <- counts_alternative(sce[genes_check, data_f1$CellType == celltype1])
ds_ref_ct2 <- counts_reference(sce[genes_check, data_f1$CellType == celltype2])
ds_alt_ct2 <- counts_alternative(sce[genes_check, data_f1$CellType == celltype2])
library_ct1 <- data_f1[,data_f1$CellType == celltype1]$Library
library_ct2 <- data_f1[,data_f1$CellType == celltype2]$Library

fit_vglm <- function(gene){
  print(gene)
  fit_df <- rbind(
    data.frame(counts_ref = ds_ref_ct1[gene, ], counts_alt = ds_alt_ct1[gene, ], 
               celltype = celltype1, library = library_ct1), 
    data.frame(counts_ref = ds_ref_ct2[gene, ], counts_alt = ds_alt_ct2[gene, ], 
               celltype = celltype2, library = library_ct2)
  )
  
  tryCatch({
    model_celltype <- VGAM::vglm(cbind(counts_ref, counts_alt) ~ celltype + library, 
                                 betabinomial, trace = FALSE, data = fit_df)
    model_constant <- VGAM::vglm(cbind(counts_ref, counts_alt) ~ 1 + library, 
                                 betabinomial, trace = FALSE, data = fit_df)
    
    return(list(gene, model_constant, model_celltype))
    
  }, error=function(cond) {
    message(paste("Not converged"))
    # Choose a return value in case of error
    model_celltype <<- NA
    model_constant <<- NA
    return(c(gene, model_constant, model_celltype))
  })
}

fit_results <- lapply(genes_check, fit_vglm)
saveRDS(fit_results, "./Data/processed/vglm_results_sc_es.rds")

pvals <- lapply(fit_results, function(x){
  gene = x[[1]]
  model_1 = x[[2]]
  model_2 = x[[3]]
  if (is.na(model_1)){
    return(c(gene, NA))
  }
  else {
    return(c(gene, lrtest(model_1, model_2)@Body[2, 5]))
  }
})

pvals_df = data.frame(do.call("rbind", pvals))
pvals_df$X2 <- as.numeric(pvals_df$X2)
rownames(pvals_df) <- pvals_df$X1

saveRDS(pvals_df, "./Data/processed/vglm_pvals_sc_es")
