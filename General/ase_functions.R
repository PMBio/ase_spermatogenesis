## ---------------------- dali functions ---------------------- 

library(reticulate)
use_python("/Users/jasper/Library/r-miniconda/bin/python", required = T)
dali_path = "~/Desktop/PhD/Projects/ASE_Spermatogenesis/Scripts/dali/"

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


## ---------------------- R functions ---------------------- 

plot_gene <- function(sce, gene, remove_zero = F){
  data_test <- data.frame(
    pt = sce$Pseudotime,
    exp_ref = as.numeric(counts_reference(sce[gene,])),
    exp_alt = as.numeric(counts_alternative(sce[gene,])),
    ase = as.numeric(allelic_ratios(sce[gene,]))
  )
  if (remove_zero){
    data_test <- data_test[data_test$exp_ref + data_test$exp_alt > 0,]
  }
  data_test$exp_total = data_test$exp_ref + data_test$exp_alt
  data_test <- data_test[order(data_test$pt),]
  
  p1 <- ggplot(data_test) + 
    geom_smooth(aes(pt, exp_ref, color = "B6")) + 
    geom_smooth(aes(pt, exp_alt, color = "Cast")) + 
    theme_classic() + xlim(c(0, 4)) + xlab("") + ylab("Expression") +
    theme(legend.position="top") + 
    scale_color_manual(values = c("B6" = "black", "Cast" = "brown"))
  
  
  p2 <- ggplot(data_test, aes(pt, ase)) +
    ylim(c(0, 1)) + geom_jitter(width = 0.1, height = 0.02, size = 0.1) + 
    scale_color_viridis() + 
    theme_classic() + geom_smooth(color = "purple") + xlim(c(0, 4)) + 
    xlab("Pseudotime") + ylab("Allelic Bias")
  
  
  
  cowplot::plot_grid(p1, p2, 
                     nrow = 2, align = "v", 
                     rel_heights = c(0.4, 0.6))
  
  
}

plot_gene_GP <- function(sce, gene, latent_ase, remove_zero = F){
  data_test <- data.frame(
    pt = sce$Pseudotime,
    exp_ref = as.numeric(counts_reference(sce[gene,])),
    exp_alt = as.numeric(counts_alternative(sce[gene,])),
    ase = as.numeric(allelic_ratios(sce[gene,])),
    latent_ase = latent_ase$posterior_mean,
    latent_var_lower = latent_ase$posterior_mean - latent_ase$posterior_var,
    latent_var_upper = latent_ase$posterior_mean + latent_ase$posterior_var
  )
  
  if (remove_zero){
    data_test <- data_test[data_test$exp_ref + data_test$exp_alt > 0,]
  }
  
  data_test$exp_total = data_test$exp_ref + data_test$exp_alt
  data_test <- data_test[order(data_test$pt),]
  
  # return(data_test)
  
  
  p1 <- ggplot(data_test) + 
    geom_smooth(aes(pt, exp_ref, color = "B6")) + 
    geom_smooth(aes(pt, exp_alt, color = "Cast")) + 
    theme_classic() + xlim(c(0, 4)) + xlab("") + ylab("Expression") +
    theme(legend.position="top") + 
    scale_color_manual(values = c("B6" = "black", "Cast" = "brown"))
  
  
  p2 <- ggplot(data_test, aes(pt, ase)) +
    ylim(c(0, 1)) + geom_jitter(width = 0.1, height = 0.02, size = 0.1) + 
    scale_color_viridis() + 
    geom_line(aes(pt, latent_ase), color = "green") + 
    geom_ribbon(aes(x = pt, ymin = latent_var_lower, ymax = latent_var_upper), 
                color = "green", alpha = 0.2) + 
    theme_classic() + xlim(c(0, 4)) + 
    xlab("Pseudotime") + ylab("Allelic Bias")
  
  cowplot::plot_grid(p1, p2, 
                     nrow = 2, align = "v", 
                     rel_heights = c(0.4, 0.6))
  
  
}

counts_reference <- function(sce){
  assays(sce)[["counts_reference"]]
}

counts_alternative <- function(sce){
  assays(sce)[["counts_alternative"]]
}

allelic_ratios <- function(sce){
  assays(sce)[["allelic_ratios"]]
}

compute_quantile_diff <- function(d, q = 0.05){
  as.numeric(abs(quantile(d, q) - quantile(d, 1 - q)))
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

### 


