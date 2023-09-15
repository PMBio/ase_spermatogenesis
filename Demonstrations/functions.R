counts_reference <- function(sce){
  assays(sce)[["counts_reference"]]
}

counts_alternative <- function(sce){
  assays(sce)[["counts_alternative"]]
}

allelic_ratios <- function(sce){
  assays(sce)[["allelic_ratios"]]
}

# initialize python path -- this should point 
library(reticulate)
use_python("/Users/jasper/Library/r-miniconda/bin/python", required = T)
dali_path = "~/Desktop/PhD/Projects/ASE_Spermatogenesis/Scripts/dali/"

# define R wrappers around scDALI functions
test_regions_R <- function(A, D, 
                           cell_state, 
                           n_cores = 1L)
{
  source_python(paste0(dali_path, "/dali/my_test_regions.py"))
  test_regions(np_array(A), 
               np_array(D), 
               np_array(cell_state), 
               n_cores = n_cores)
}

test_mean_R <- function(a, d, mean_mu = 0.5)
{
  source_python(paste0(dali_path, "/dali/my_test_mean.py"))
  res = test_mean(np_array(a), 
                  np_array(d), 
                  mean_mu = mean_mu)
  names(res) = c("Mean", "Theta", "NLL_null", "NLL_alt", "pval")
  res
}

run_gp_R <- function(A, D, 
                     cell_state, 
                     kernel = "Linear",
                     n_cores = 1L)
{
  source_python(paste0(dali_path, "/dali/my_run_gp.py"))
  res = run_gp(np_array(A), 
               np_array(D), 
               np_array(cell_state), 
               kernel = kernel,
               n_cores = n_cores)
  res
}

# This function takes a cell x genes matrix and a pseudotemporal ordering vector and 
# partitions the matrix into n_partitions intervals with even width in pseudotime space, 
# computing average values per gene
make_pseudotime_smooth <- function(data, pt, n_partitions = 100){
  pseudotime_here <- pt
  pseudotime_here <- pseudotime_here - min(pseudotime_here)
  pseudotime_here <- pseudotime_here / max(pseudotime_here) # between 0 and 1
  intervals <- c(0, 1:n_partitions / n_partitions)
  
  expression_smoothed <- do.call("cbind", 
                                 lapply(1:n_partitions, function(i){
                                   p_lower <- intervals[i]
                                   p_upper <-  intervals[i + 1]
                                   ix <- which(p_lower <= pseudotime_here & pseudotime_here <= p_upper)
                                   rowMeans(data[,ix])
                                 }))
  expression_smoothed
}

### plotting functions
plot_gene_GP <- 
  function(sce, gene, latent_ase, remove_zero = F, scale_var = 1.96, textsize = 20){
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

plot_gene <- 
  function(sce, gene, remove_zero = F, scale_var = 1.96, textsize = 20){
    data_test <- data.frame(
      pt_here = sce$Pseudotime,
      exp_ref = as.numeric(counts_reference(sce[gene,])),
      exp_alt = as.numeric(counts_alternative(sce[gene,]))
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



#### 



plot_cis_trans_gene_mean <- function(gene){
  
  rev_logit <- function(x){
    1 / (1 + exp(-x))
  }
  
  asd <- results_per_gene_expressed_intervals[[gene]]
  
  for_plotting <- lapply(1:length(asd[[5]][[2]]), function(i){
    xx = asd[[5]][[2]][[i]]
    xx = data.frame(do.call("cbind", xx))
    xx$Sample = i
    xx$X3 <- sqrt(xx$X3)
    xx})
  
  for_plotting_merged <- do.call("rbind", for_plotting)
  for_plotting_merged$Sample <- as.factor(for_plotting_merged$Sample)
  for_plotting_merged$mean <- for_plotting_merged$X2
  
  data_plot <- data.frame(
    "F0_P" = data_f0_b6[gene, ],
    "F0_M" = data_f0_cast[gene, ],
    "F1_P" = data_f1_b6[gene, ],
    "F1_M" = data_f1_cast[gene, ]
  )
  
  norm_factors <- c(sum(colSums(data_f0_b6)), 
                    sum(colSums(data_f0_cast)))
  norm_factors <- norm_factors / norm_factors[[1]]
  
  data_plot[,1:2] <- data.frame(t(t(data_plot[,1:2]) / norm_factors))
  data_plot$Interval <- 1:nrow(data_plot)
  for_plotting_merged$Sample <- ifelse(for_plotting_merged$Sample == 1, "F0", "F1")
  
  ggplot(for_plotting_merged, aes(X1, y = rev_logit(X2), col = Sample)) + 
    geom_line(linetype = "dashed", size = 1.5) + 
    scale_color_manual(values = c("F0" = "black", "F1" = "darkgrey")) + 
    # geom_ribbon(aes(ymin = rev_logit(error_lower), ymax = rev_logit(error_upper), fill = Sample), 
    #             alpha = 0.2) + 
    theme(legend.position="top") + 
    geom_point(data = data_plot, aes(Interval, (F0_P) / (F0_M + F0_P)), col = "black", size = 2) + 
    geom_point(data = data_plot, aes(Interval, (F1_P) / (F1_M + F1_P)), col = "darkgrey", size = 2) + 
    xlim(0, 100) + scale_fill_manual(values = c("red", "darkgrey")) + ylim(c(0, 1)) + xlab("Pseudotime") + 
    ylab("Allelic Imbalance") + theme_classic() + ggtitle(gene) + 
    geom_hline(yintercept = 0.5, linetype = "dashed") + 
    theme(text = element_text(size = 30), legend.position = "None") + 
    scale_x_continuous(breaks = c(0, 50, 100)) + 
    scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1))
}