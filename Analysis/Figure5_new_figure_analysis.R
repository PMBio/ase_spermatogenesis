theme_paper <- function(textsize = 20){
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1), 
        text = element_text(size = textsize))
}
setwd("~/Desktop/Projects/ASE_Spermatogenesis_Paper_rerun/")

results_per_gene_expressed_intervals_cast <- readRDS("./Data/processed/gp_fits_DE_cast_new.rds")
results_per_gene_expressed_intervals <- results_per_gene_expressed_intervals_cast
results_per_gene_expressed_intervals <- results_per_gene_expressed_intervals[!unlist(lapply(results_per_gene_expressed_intervals, function(x){any(is.na(x))}))]

likelihood_df <- data.frame(
  likelihood_cons = unlist(lapply(results_per_gene_expressed_intervals, function(x){x[[1]][[3]]})),
  likelihood_stat = unlist(lapply(results_per_gene_expressed_intervals, function(x){x[[2]][[3]]})),
  likelihood_dyn = unlist(lapply(results_per_gene_expressed_intervals, function(x){x[[3]][[3]]}))
)
bf_dynamic_deg_cast <- likelihood_df$likelihood_dyn - likelihood_df$likelihood_stat
names(bf_dynamic_deg_cast) <- names(results_per_gene_expressed_intervals)

results_per_gene_expressed_intervals_caroli <- readRDS("./Data/processed/gp_fits_DE_caroli_new.rds")
results_per_gene_expressed_intervals <- results_per_gene_expressed_intervals_caroli
results_per_gene_expressed_intervals <- results_per_gene_expressed_intervals[!unlist(lapply(results_per_gene_expressed_intervals, function(x){any(is.na(x))}))]

likelihood_df <- data.frame(
  likelihood_cons = unlist(lapply(results_per_gene_expressed_intervals, function(x){x[[1]][[3]]})),
  likelihood_stat = unlist(lapply(results_per_gene_expressed_intervals, function(x){x[[2]][[3]]})),
  likelihood_dyn = unlist(lapply(results_per_gene_expressed_intervals, function(x){x[[3]][[3]]}))
)
bf_dynamic_deg_caroli <- likelihood_df$likelihood_dyn - likelihood_df$likelihood_stat
names(bf_dynamic_deg_caroli) <- names(results_per_gene_expressed_intervals)

ggplot(data.frame( bf_cast = bf_dynamic_deg_cast, bf_caroli = bf_dynamic_deg_caroli[names(bf_dynamic_deg_cast)]), aes(bf_cast, bf_caroli)) + geom_point() + 
  coord_fixed() + xlim(-0.1, 210) + ylim(-0.1, 210) + theme_paper(textsize = 30) + geom_smooth() + 
  xlab("CAST (log BF)") + ylab("CAROLI (log BF)")
ggsave("./Plots/FigureS5/FigS5_9.pdf")

table(
  CAST_vs_B6 = bf_dynamic_deg_cast[intersect(names(bf_dynamic_deg_cast), names(bf_dynamic_deg_caroli))] > 10, 
  CAROLI_vs_B6 = bf_dynamic_deg_caroli[intersect(names(bf_dynamic_deg_cast), names(bf_dynamic_deg_caroli))] > 10
) %>% data.frame() %>%
  ggplot(aes(CAST_vs_B6, CAROLI_vs_B6, fill = Freq)) + geom_tile() + 
  geom_text(aes(label = Freq), size = 20) + xlab("CAST vs B6") + ylab("CAROLI vs B6") + 
  scale_fill_gradient(low = "white", high = "darkblue", limits = c(0, 2500)) + 
  theme_paper() + theme(text = element_text(size = 40))
ggsave("./Plots/FigureS5/FigS5_10.pdf")

cor(bf_dynamic_deg_cast, rank(bf_dynamic_deg_caroli[names(bf_dynamic_deg_cast)]))
cor(bf_dynamic_deg_cast, rank(bf_dynamic_deg_caroli[names(bf_dynamic_deg_cast)]), method = "spearman")

# joint clustering approach here as well

data <- readRDS("./Data/processed/data_evo_for_fitting.rds")

genes_test <- rownames(data[[1]])
genes_test <- genes_test[!genes_test %in% c("Krt82", "Uroc1", "Rax")]

data_f0_b6 <- data[[1]]
data_f0_cast <- data[[2]]
data_f0_caroli <- data[[3]]

sfs_vector_b6_cast <- colSums(data_f0_b6) / colSums(data_f0_cast)
sfs_vector_b6_caroli <- colSums(data_f0_b6) / colSums(data_f0_caroli)
sfs_vector_b6_caroli[is.na(sfs_vector_b6_caroli)] <- 1
data_f0_b6_norm <- data_f0_b6 / 1
data_f0_cast_norm <- t(t(data_f0_cast) * sfs_vector_b6_cast)
data_f0_caroli_norm <- t(t(data_f0_caroli) * sfs_vector_b6_caroli)

rs_f0_noStab_cast <- (data_f0_b6_norm[genes_test, ]) / (data_f0_b6_norm[genes_test, ] + data_f0_cast_norm[genes_test, ])
rs_f0_cast <- (data_f0_b6_norm[genes_test, ] + 1) / (data_f0_b6_norm[genes_test, ] + data_f0_cast_norm[genes_test, ] + 2)
rs_f0_noStab_caroli <- (data_f0_b6_norm[genes_test, ]) / (data_f0_b6_norm[genes_test, ] + data_f0_caroli_norm[genes_test, ])
rs_f0_caroli <- (data_f0_b6_norm[genes_test, ] + 1) / (data_f0_b6_norm[genes_test, ] + data_f0_caroli_norm[genes_test, ] + 2)

total_exp <- data_f0_b6 + data_f0_cast + data_f0_caroli
total_exp <- total_exp[genes_test,]

set.seed(123)
scale_to_zero_one <- function(x){(x - min(x)) / (max(x) - min(x))}
make_meta_profiles <- function(dataset, k = 5){
  dis_mat <- dist(dataset)
  expression_clusters <- cutree(hclust(dis_mat, method = "complete"), k = k)
  cluster_counts <- table(expression_clusters)
  print(cluster_counts)
  res <- lapply(unique(expression_clusters), function(i){
    dd <- dataset[expression_clusters == i, ]
    if (  cluster_counts[[i]] == 1 ){
      return(dd)
    }
    return(colMeans(dd, na.rm = T))
  })
  return(do.call("rbind", res))
}

rs_f0_plot_cast <- rs_f0_cast[names(bf_dynamic_deg_cast[bf_dynamic_deg_cast > 50]), ]
rs_f0_plot_cast <-  data.frame(t(rollapply(t(rs_f0_plot_cast), width = 10, by = 1, FUN = mean, align = "left")))
rs_f0_plot_cast <- abs(rs_f0_plot_cast - 0.5)
rs_f0_plot_unscaled_cast <- rs_f0_plot_cast
rs_f0_plot_cast <- t(apply(rs_f0_plot_cast, 1, scale_to_zero_one))
cast_input <- rs_f0_plot_cast

rs_f0_plot_caroli <- rs_f0_caroli[names(bf_dynamic_deg_caroli[bf_dynamic_deg_caroli > 50]), ]
rs_f0_plot_caroli <-  data.frame(t(rollapply(t(rs_f0_plot_caroli), width = 10, by = 1, FUN = mean, align = "left")))
rs_f0_plot_caroli <- abs(rs_f0_plot_caroli - 0.5)
rs_f0_plot_unscaled_caroli <- rs_f0_plot_caroli
rs_f0_plot_caroli <- t(apply(rs_f0_plot_caroli, 1, scale_to_zero_one))
caroli_input <- rs_f0_plot_caroli

# combine datasets
cast_names <- paste0("CAST_", rownames(cast_input))
caroli_names <- paste0("CAROLI_", rownames(caroli_input))

combined_input <- rbind(cast_input, caroli_input)
rownames(combined_input) <- c(cast_names, caroli_names)

set.seed(123)
n_clusters <- 7
clusters <- cutree(hclust(dist(combined_input), method = "complete"), k = n_clusters)

trace_results <- lapply(1:7, function(i){
  dds = combined_input[clusters == i, ]
  cast_trace = colMeans(dds[grepl("CAST", rownames(dds)), ])
  caroli_trace = colMeans(dds[grepl("CAROLI", rownames(dds)), ])
  dds <- data.frame(do.call("rbind", list(cast_trace, caroli_trace)))
  dds$Cluster = paste0("Cluster_", i)
  dds$View = c(c("cast", "caroli"))
  dds
})
trace_results <- data.frame(do.call("rbind", trace_results))

trace_results %>%
  pivot_longer(-c("Cluster", "View")) %>%
  mutate(name = as.numeric(gsub("X", "", name))) %>%
  ggplot(aes(x = name, y = value, col = View)) + geom_point() + 
  facet_wrap(~Cluster) + 
  theme_classic()

# Trace plots for individual clusters

late_clusters <- paste0("Cluster_", c(1, 3, 4, 5))
cluster_colors <- rep("grey", 7)
names(cluster_colors) <- paste0("Cluster_", 1:7)
cluster_colors[late_clusters] <- RColorBrewer::brewer.pal(name = "Dark2", 4)

# for cast
trace_results %>% 
  dplyr::filter(View == "cast") %>%
  pivot_longer(-c("Cluster", "View")) %>%
  mutate(name = as.numeric(gsub("X", "", name))) %>%
  dplyr::filter(Cluster %in% late_clusters) %>%
ggplot() + 
  geom_smooth(span = 0.5, size = 2, aes(x = name, y = value, col = Cluster)) + 
  theme_classic() +
  theme(text = element_text(size = 30)) +
  xlab("Pseudotime") +
  ylab("Relative effect size") +
  theme(legend.position = "None") +
  scale_color_manual(values = cluster_colors) +
  ggtitle("Differential Expression (F0, CAST)") +
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
ggsave("./Plots/Figure5/CAST_trace_plot.pdf", width = 16, height = 5)

data.frame(
  Cluster = factor(paste0("Cluster_", unique(clusters)), levels = paste0("Cluster_", c(1, 3, 4, 5, 2, 6, 7))),
  nGenes = as.numeric(table(clusters))
) %>% ggplot(aes(x="", y = nGenes, fill = as.factor(Cluster))) +
  geom_bar(stat = "identity", width = 1, color = "black") +
  coord_polar("y", start=0) + 
  scale_fill_manual(values = cluster_colors) + 
  theme_void()
ggsave("./Plots/Figure5/CAST_venn.pdf",  width = 20, height = 20)

# for caroli
trace_results %>% 
  dplyr::filter(View == "caroli") %>%
  pivot_longer(-c("Cluster", "View")) %>%
  mutate(name = as.numeric(gsub("X", "", name))) %>%
  dplyr::filter(Cluster %in% late_clusters) %>%
ggplot() + 
  geom_smooth(span = 0.5, size = 2, aes(x = name, y = value, col = Cluster)) + 
  theme_classic() +
  theme(text = element_text(size = 30)) +
  xlab("Pseudotime") +
  ylab("Relative effect size") +
  theme(legend.position = "None") +
  scale_color_manual(values = cluster_colors) +
  ggtitle("Differential Expression (F0, CAROLI)") +
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
ggsave("./Plots/Figure5/CAROLI_trace_plot.pdf", width = 16, height = 5)

data.frame(
  Cluster = factor(paste0("Cluster_", unique(clusters)), levels = paste0("Cluster_", c(1, 3, 4, 5, 2, 6, 7))),
  nGenes = as.numeric(table(clusters))
) %>% ggplot(aes(x="", y = nGenes, fill = as.factor(Cluster))) +
  geom_bar(stat = "identity", width = 1, color = "black") +
  coord_polar("y", start=0) + 
  scale_fill_manual(values = cluster_colors) + 
  theme_void()
ggsave("./Plots/Figure5/CAROLI_venn.pdf", width = 20, height = 20)

total_effect_sizes_unscaled <- data.frame(
  "Bins" = 1:91,
  "total" = (colSums(data_f0_b6) + colSums(data_f0_cast))[1:91],
  "deg" = colMeans(rs_f0_plot_unscaled_cast)[1:91]
) %>% pivot_longer(-c("Bins"))

ggplot(total_effect_sizes_unscaled[total_effect_sizes_unscaled$name == "deg", ],
       aes(x = Bins, y = value, col = name)) +
  geom_point(alpha = 0.3) +
  geom_smooth(span = 0.2) +
  geom_smooth(span = 1, col = "black", linetype = 'dashed') +
  theme_classic() +
  theme(text = element_text(size = 30),
        legend.position = "None") +
  ylab("") + xlab("") +
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
ggsave("./Plots/Figure5/F5_clustering_CAST_unscaled.pdf",  width = 16, height = 2.5)

total_effect_sizes_unscaled <- data.frame(
  "Bins" = 1:91,
  "total" = (colSums(data_f0_b6) + colSums(data_f0_caroli))[1:91],
  "deg" = colMeans(rs_f0_plot_unscaled_caroli)[1:91]
) %>% pivot_longer(-c("Bins"))

ggplot(total_effect_sizes_unscaled[total_effect_sizes_unscaled$name == "deg", ],
       aes(x = Bins, y = value, col = name)) +
  geom_point(alpha = 0.3) +
  geom_smooth(span = 0.2) +
  geom_smooth(span = 1, col = "black", linetype = 'dashed') +
  theme_classic() +
  theme(text = element_text(size = 30),
        legend.position = "None") +
  ylab("") + xlab("") +
  theme(axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
ggsave("./Plots/Figure5/F5_clustering_CAROLI_unscaled.pdf",  width = 16, height = 2.5)

