###############################################
#### Script to collect auxiliary functions ####
###############################################

# Load libraries needed in this script
library(scater)
library(scran)
library(princurve)
library(edgeR)

#### Split Single cell experiment
split.sce <- function(sce, groups, colData.name = "SubCluster"){
  # List to collect individual single cell experiments
  list.out <- list()
  for(i in groups){
    cur_sce <- sce[,as.character(colData(sce)[[colData.name]]) == as.character(i)]
    cur_sce <- computeSumFactors(cur_sce)
    cur_sce <- logNormCounts(cur_sce)
    list.out[[i]] <- cur_sce
  }
  names(list.out) <- groups
  list.out
}


#### HVG
# Calculate highly variable genes
# sce: one or multiple single cell experiment objects in form of a list
HVG <- function(sce, numberGenes = 1000){
  # One single cell expereiment object
  if(typeof(sce) == "S4"){
    getTopHVGs(modelGeneVar(sce), n = numberGenes)
  }
  # Multiple single cell experiment objects
  else if(typeof(sce) == "list") {
    lapply(sce, function(n){
      getTopHVGs(modelGeneVar(sce), n = numberGenes)
    })
  }
  else if(class(sce) == "matrix" | class(sce) == "data.frame"){
    getTopHVGs(modelGeneVar(sce), n = numberGenes)
  }
  else{
    print("First argument must either be a single cell experiment object \n
          or a list")
  }
}

#### Find specifc marker genes
marker.detection <- function(sce, clusters, blocking = NULL){
  # User scran function findMarkers to perform differential expression
  cur_markers <- findMarkers(sce, clusters, block = blocking)

  # Collect group specific markers
  markers.spec <- lapply(cur_markers, function(n){
    if(!is.na(n$Top[1])){
    cur_n <- n[n$FDR < 0.1 & apply(as.matrix(n[,4:ncol(n)]), 1, function(x){sum(x > 0)}) == ncol(n) - 3,]
      if(nrow(cur_n) > 0){
        cur_n$GeneName <- rowData(sce)$Symbol[match(rownames(cur_n), rowData(sce)$ID)]
      }
    }
    else{
      cur_n <- NULL
    }
    cur_n
  })
  
}

#### Compute pseudorank
PT <- function(rd, clusters, col_vector, 
               exclude = NULL, start = NULL, end = NULL){
    if(!is.null(exclude)){
      cur_rd <- rd[!exclude,]
      
      cur_lin <- principal_curve(cur_rd, start = start)
      
      plot(cur_rd, col = col_vector[clusters[!exclude]], 
           pch = 16, type = "p")
      lines(cur_lin, lwd = 3)
      
      mat.out <- matrix(data = NA, ncol = ncol(cur_rd) + 1, nrow = length(clusters))
      rownames(mat.out) <- names(clusters)
      colnames(mat.out) <- c(colnames(cur_rd), "rank")
      
      mat.out[!exclude,1:ncol(cur_rd)] <- cur_lin$s
      mat.out[!exclude,"rank"] <- order(cur_lin$tag)
      
      mat.out
    }
    else{
      cur_rd <- rd
      
      cur_lin <- principal_curve(cur_rd, start = start)
      
      plot(cur_rd, col = col_vector[clusters], 
           pch = 16, type = "p")
      lines(cur_lin, lwd = 3)
      
      mat.out <- matrix(data = NA, ncol = ncol(cur_rd) + 2, nrow = length(clusters))
      rownames(mat.out) <- names(clusters)
      colnames(mat.out) <- c(colnames(cur_rd), "rank", "lambda")
      
      mat.out[,1:ncol(cur_rd)] <- cur_lin$s
      mat.out[,"rank"] <- order(cur_lin$ord)
      mat.out[,"lambda"] <- cur_lin$lambda
      
      mat.out
    }
}

#### Batch correction
#### Batch correction
batch.correction <- function(sce, number.HVG = 1000){
  library(batchelor)
  # Calculate highly variable genes and merge
  HVG.genes <- lapply(sce, function(sce.part){
    HVG <- modelGeneVar(sce.part)
  })
  
  HVG.df <- do.call("combineVar", HVG.genes)
  rownames(HVG.df) <- rownames(HVG.genes[[1]]) ### added by jasper - make sure this is indeed the correct order
  HVG.df <- HVG.df[order(HVG.df$bio, decreasing = TRUE),]
  genes <- rownames(HVG.df)[1:number.HVG]
  
  # Batch correction
  func <- paste0("fastMNN(", 
                 paste0("as.matrix(logcounts(sce[[", 1:length(sce), "]])[genes,])", collapse=", "), 
                 ")")
  corrected <- eval( parse(text=func) )
  return(as.matrix(t(assays(corrected)[["reconstructed"]])))
#  t(corrected$corrected)
}

batch.correction_old <- function(sce, number.HVG = 1000){
  # Calculate highly variable genes and merge
  HVG.genes <- lapply(sce, function(n){
    HVG <- trendVar(n, use.spikes = FALSE)
    decomposeVar(n, HVG)
  })
  
  HVG.df <- do.call("combineVar", HVG.genes)
  HVG.df <- HVG.df[order(HVG.df$bio, decreasing = TRUE),]
  genes <- rownames(HVG.df)[1:number.HVG]
  
  # Batch correction
  #func <- paste0("mnnCorrect(", 
  #                   paste0("as.matrix(logcounts(sce[[", 1:length(sce), "]])[genes,])", collapse=", "), 
  #                   ", cos.norm.in=TRUE, cos.norm.out=TRUE, sigma=0.1)")
  func <- paste0("mnnCorrect(", 
                     paste0("as.matrix(logcounts(sce[", 1:length(sce), "])[genes,])", collapse=", "), 
                     ", cos.norm.in=TRUE, cos.norm.out=TRUE, sigma=0.1)")
  corrected <- eval( parse(text=func) )
  do.call("cbind", corrected$corrected)
}

# Differnetial expression testing using edgeR
DE.edgeR <- function(sce, conditions, covariate, lfc, FDR){
  # Collect summed counts in matrix
  mat <- matrix(data = NA, 
                ncol = length(unique(paste(covariate, conditions, sep = " "))), 
                nrow = nrow(sce))
  rownames(mat) <- rownames(counts(sce))
  colnames(mat) <- unique(paste(covariate, conditions, sep = " "))
  
  for(j in colnames(mat)){
    cur_covariate <- unlist(strsplit(j, " "))[1]
    cur_condition <- unlist(strsplit(j, " "))[2]
    mat[,j] <- Matrix::rowSums(counts(sce)[,covariate == cur_covariate &
                                             conditions == cur_condition]) 
  }
  
  # Perform differential testing
  print(sapply(colnames(mat), 
               function(n){unlist(strsplit(n, " "))[1]}))
  y <- DGEList(counts=mat,
               group=sapply(colnames(mat), 
                            function(n){unlist(strsplit(n, " "))[1]}))
  y <- calcNormFactors(y)
  
  # Generate design matrix
  design <- model.matrix(~0+sapply(colnames(mat), function(n){unlist(strsplit(n, " "))[2]}))
  colnames(design) <- substring(colnames(design), regexpr("})", colnames(design)) + 2)
  y <- estimateDisp(y,design)
  
  # Fit the model
  fit <- glmQLFit(y,design, robust = TRUE)
  qlf <- glmTreat(fit,coef=2, lfc = lfc, 
                  contrast = eval(parse(text = paste("makeContrasts(", colnames(design)[1],  " - ", 
                                                     colnames(design)[2], ", levels = design)", sep = ""))))
  cur_markers <- topTags(qlf, n = nrow(qlf$table))$table
  
  # Save markers
  cur_out <- list()
  cur_out[[colnames(design)[2]]] <- cur_markers[cur_markers$logFC < 0 & cur_markers$FDR <= FDR,]
  cur_out[[colnames(design)[2]]]$Genename <- rowData(sce)$Symbol[match(rownames(cur_out[[colnames(design)[2]]]),
                                                                       rowData(sce)$ID)]
  cur_out[[colnames(design)[1]]] <- cur_markers[cur_markers$logFC > 0 & cur_markers$FDR <= FDR,]
  cur_out[[colnames(design)[1]]]$Genename <- rowData(sce)$Symbol[match(rownames(cur_out[[colnames(design)[1]]]),
                                                                       rowData(sce)$ID)]
  cur_out
}

# Multi-group DE
multi.DE <- function(sce, conditions, covariate, select.marker = TRUE, lfc, FDR){
  # Select unique conditions
  conditions <- sub("-", "_", conditions)
  cond <- unique(conditions)
  
  # List to store tests
  cur_out <- list()
  
  # Loop through different comparisons
  for(i in seq(1,length(cond)-1)){
    for(j in seq(i+1, length(cond), 1)){
      cur_sce <- sce[,conditions %in% cond[c(i,j)]]
      cur_conditions <- conditions[conditions %in% cond[c(i,j)]]
      cur_covariate <- covariate[conditions %in% cond[c(i,j)]]
      
      # Collect summed counts in matrix
      mat <- matrix(data = NA, 
                    ncol = length(unique(paste(cur_covariate, cur_conditions, sep = " "))), 
                    nrow = nrow(cur_sce))
      rownames(mat) <- rownames(counts(cur_sce))
      colnames(mat) <- unique(paste(cur_covariate, cur_conditions, sep = " "))
      
      for(k in colnames(mat)){
        cur_cov <- unlist(strsplit(k, " "))[1]
        cur_cond <- unlist(strsplit(k, " "))[2]
        if(sum(cur_covariate == cur_cov &
               cur_conditions == cur_cond) < 2){
          next
        }
        mat[,k] <- Matrix::rowSums(counts(cur_sce)[,cur_covariate == cur_cov &
                                                     cur_conditions == cur_cond]) 
      }
      mat <- mat[,!is.na(mat[1,])]
      
      # Perform differential testing
      y <- DGEList(counts=mat,
                   group=sapply(colnames(mat), 
                                function(n){unlist(strsplit(n, " "))[1]}))
      y <- calcNormFactors(y)
      
      # Generate design matrix
      design <- model.matrix(~0+sapply(colnames(mat), function(n){unlist(strsplit(n, " "))[2]}))
      colnames(design) <- substring(colnames(design), regexpr("})", colnames(design)) + 2)
      print(design)
      y <- estimateDisp(y,design)
      
      # Fit the model
      fit <- glmQLFit(y,design, robust = TRUE)
      qlf <- glmTreat(fit,coef=2, lfc = lfc, 
                      contrast = eval(parse(text = paste("makeContrasts(", colnames(design)[1],  " - ", 
                                                         colnames(design)[2], ", levels = design)", sep = ""))))
      # Save results
      cur_out[[paste(colnames(design)[1], colnames(design)[2])]] <- qlf$table
    }
  }
  
  # Combine results into dataframe
  final_out <- list()
  for(i in cond){
    cur_df <- data.frame(row.names = rownames(sce))
    cur_p <- data.frame(row.names = rownames(sce))
    cur_comp <- cur_out[unlist(lapply(strsplit(names(cur_out), " "),
                                      function(n){n[1] == i | n[2] == i}))] 
    for(j in 1:length(cur_comp)){
      if(which(unlist(strsplit(names(cur_comp)[j], " ")) == i) == 1){
        cur_df[[unlist(strsplit(names(cur_comp)[j], " "))[2]]] <- cur_comp[[j]]$logFC
        cur_p[[unlist(strsplit(names(cur_comp)[j], " "))[2]]] <- cur_comp[[j]]$PValue
      }
      else{
        cur_df[[unlist(strsplit(names(cur_comp)[j], " "))[1]]] <- -cur_comp[[j]]$logFC
        cur_p[[unlist(strsplit(names(cur_comp)[j], " "))[1]]] <- cur_comp[[j]]$PValue
      }
    }
    
    # Add P values to data.fram
    cur_df <- cbind(cur_df, cur_p)
    colnames(cur_df) <- paste(colnames(cur_df), c(rep("logFC", ncol(cur_df)/2),
                                                  rep("PValue", ncol(cur_df)/2)), sep = "_")  
    
    # Combine P values
    cur_df$Combined_PValue <- eval(parse(text = paste0("combinePValues(", 
                                                       paste0("cur_p[,", 1:ncol(cur_p), "],", collapse = " "), 
                                                       " method = 'z')")))
    cur_df$FDR <- p.adjust(cur_df$Combined_PValue, method = "fdr")
    
    # Add the gene names
    cur_df$ID <- rownames(cur_df)
    cur_df$Symbol <- rowData(sce)$Symbol[match(rownames(cur_df), rowData(sce)$ID)]
    
    # Order genes pased on P value
    cur_df <- cur_df[order(cur_df$FDR, decreasing = FALSE),]
    
    # Exclude genes with FDR smaller than the specified value
    cur_df <- cur_df[cur_df$FDR <= FDR,]
    
    # Select only genes with throught positive logFC if chosen by user
    if(select.marker){
      cur_df <- cur_df[apply(cur_df[,1:(ncol(cur_df) - 4)], 1, function(n){
        sum(n > 0) == ncol(cur_df) - 4
      }),]
    }
    
    # Save in list
    final_out[[i]] <- cur_df
  }
  # Return object
  final_out
}

# Pairwise distances
pairwiseCor <- function(x = NULL, y = NULL, method = "spearman"){
  if(is.null(x)){
    stop("This function takes one or two SingleCellExperiment objects as input")
  }
  
  if(is.null(y)){
    # Cell-to-cell distance
    cur_cor <- as.dist(cor(as.matrix(logcounts(x)), method = method))
    cur_cor <- as.vector(cur_cor)
  }
  else{
    cur_cor <- cor(as.matrix(logcounts(x)),
                                          as.matrix(logcounts(y)), 
                                          method = method)
    cur_cor <- as.vector(cur_cor)
  }
  cur_cor
}


#### Compute pseudotime with destiny
diffusionPT <- function(sce, HVG, clusters, col_vector,
                        exclude = NULL){
  if(!is.null(exclude)){
    dm <- DiffusionMap(t(as.matrix(logcounts(sce)[HVG,!exclude])), k = 20)
    
    plot(dm, col = col_vector[clusters[!exclude]], 
         pch = 16, type = "p")
    
    dpt <- DPT(dm = dm)
    
    dpt$DPT1
  }
  else{
    dm <- DiffusionMap(t(as.matrix(logcounts(sce)[HVG,])), k = 20)
    
    plot(dm, col = col_vector[clusters], 
         pch = 16, type = "p")
    
    dpt <- DPT(dm = dm)
    
    dpt$DPT1
  }
}
