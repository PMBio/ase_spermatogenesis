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

## from https://github.com/ChrKoenig/R_marginal_plot/blob/master/marginal_plot.R
marginal_plot = function(x, y, group = NULL, data = NULL, lm_show = FALSE, lm_formula = y ~ x, bw = "nrd0", adjust = 1, alpha = 1, plot_legend = T, ...)
{
  require(scales)
  ###############
  # Plots a scatterplot with marginal probability density functions for x and y. 
  # Data may be grouped or ungrouped. 
  # For each group, a linear fit can be plotted. It is hidden by default, but can be shown by providing lm_show = TRUE.
  # The model can be modified using the 'lm_formula' argument.
  # The 'bw' and 'adjust' argument specify the granularity used for estimating probability density functions. See ?density for more information.
  # For large datasets, opacity may be decreased by setting alpha to a value between 0 and 1.
  # Additional graphical parameters are passed to the main plot, so you can customize axis labels, titles etc.
  ###############
  moreargs = eval(substitute(list(...)))
  
  # prepare consistent df
  if(missing(group)){
    if(missing(data)){
      if(length(x) != length(y)){stop("Length of arguments not equal")}
      data = data.frame(x = as.numeric(x), y = as.numeric(y))
    } else {
      data = data.frame(x = as.numeric(data[,deparse(substitute(x))]), 
                        y = as.numeric(data[,deparse(substitute(y))]))
    }
    if(sum(!complete.cases(data)) > 0){
      warning(sprintf("Removed %i rows with missing data", sum(!complete.cases(data))))
      data = data[complete.cases(data),]
    }
    group_colors = "black"
  } else {
    if(missing(data)){
      if(length(x) != length(y) | length(x) != length(group)){stop("Length of arguments not equal")}
      data = data.frame(x = as.numeric(x), y = as.numeric(y), group = as.factor(group))
    } else {
      data = data.frame(x = as.numeric(data[,deparse(substitute(x))]), 
                        y = as.numeric(data[,deparse(substitute(y))]),
                        group = as.factor(data[,deparse(substitute(group))]))
    }
    if(sum(!complete.cases(data)) > 0){
      warning(sprintf("Removed %i rows with missing data", sum(!complete.cases(data))))
      data = data[complete.cases(data),]
    }
    data = subset(data, group %in% names(which(table(data$group) > 5)))
    data$group = droplevels(data$group)
    group_colors = rainbow(length(unique(data$group)))
  } 
  
  # log-transform data (this is need for correct plotting of density functions)
  if(!is.null(moreargs$log)){
    if(!moreargs$log %in% c("y", "x", "yx", "xy")){
      warning("Ignoring invalid 'log' argument. Use 'y', 'x', 'yx' or 'xy.")
    } else {
      data = data[apply(data[unlist(strsplit(moreargs$log, ""))], 1, function(x) !any(x <= 0)), ]
      data[,unlist(strsplit(moreargs$log, ""))] = log10(data[,unlist(strsplit(moreargs$log, ""))])
    }
    moreargs$log = NULL # remove to prevent double logarithm when plotting
  }
  
  # Catch unwanted user inputs
  if(!is.null(moreargs$col)){moreargs$col = NULL}
  if(!is.null(moreargs$type)){moreargs$type = "p"}
  
  # get some default plotting arguments
  if(is.null(moreargs$xlim)){moreargs$xlim = range(data$x)} 
  if(is.null(moreargs$ylim)){moreargs$ylim = range(data$y)}
  if(is.null(moreargs$xlab)){moreargs$xlab = deparse(substitute(x))}
  if(is.null(moreargs$ylab)){moreargs$ylab = deparse(substitute(y))}
  if(is.null(moreargs$las)){moreargs$las = 1} 
  
  # plotting
  tryCatch(expr = {
    ifelse(!is.null(data$group), data_split <- split(data, data$group), data_split <- list(data))
    orig_par = par(no.readonly = T)
    par(mar = c(0.25,5,1,0))
    layout(matrix(1:4, nrow = 2, byrow = T), widths = c(10,3), heights = c(3,10))
    
    # upper density plot
    plot(NULL, type = "n", xlim = moreargs$xlim, ylab = "density",
         ylim = c(0, max(sapply(data_split, function(group_set) max(density(group_set$x, bw = bw)$y)))), main = NA, axes = F)
    axis(2, las = 1)
    mapply(function(group_set, group_color){lines(density(group_set$x, bw = bw, adjust = adjust), col = group_color, lwd = 2)}, data_split, group_colors)
    
    # legend
    par(mar = c(0.25,0.25,0,0))
    plot.new()
    if(!missing(group) & plot_legend){
      legend("center", levels(data$group), fill = group_colors, border = group_colors, bty = "n", title = deparse(substitute(group)), title.adj = 0.1)
    }
    
    # main plot
    par(mar = c(4,5,0,0))
    if(missing(group)){
      do.call(plot, c(list(x = quote(data$x), y = quote(data$y), col = quote(scales::alpha("black", alpha))), moreargs))
    } else {
      do.call(plot, c(list(x = quote(data$x), y = quote(data$y), col = quote(scales::alpha(group_colors[data$group], alpha))), moreargs))
    }
    axis(3, labels = F, tck = 0.01)
    axis(4, labels = F, tck = 0.01)
    box()
    
    if(lm_show == TRUE & !is.null(lm_formula)){
      mapply(function(group_set, group_color){
        lm_tmp = lm(lm_formula, data = group_set)
        x_coords = seq(min(group_set$x), max(group_set$x), length.out = 100)
        y_coords = predict(lm_tmp, newdata = data.frame(x = x_coords))
        lines(x = x_coords, y = y_coords, col = group_color, lwd = 2.5)
      }, data_split, rgb(t(ceiling(col2rgb(group_colors)*0.8)), maxColorValue = 255))
    }
    
    # right density plot
    par(mar = c(4,0.25,0,1))
    plot(NULL, type = "n", ylim = moreargs$ylim, xlim = c(0, max(sapply(data_split, function(group_set) max(density(group_set$y, bw = bw)$y)))), main = NA, axes = F, xlab = "density")
    mapply(function(group_set, group_color){lines(x = density(group_set$y, bw = bw, adjust = adjust)$y, y = density(group_set$y, bw = bw)$x, col = group_color, lwd = 2)}, data_split, group_colors)
    axis(1)
  }, finally = {
    par(orig_par)
  })}

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


