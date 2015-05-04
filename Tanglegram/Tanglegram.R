library(dendextend)
library(gplots)
library(reshape)
library(RColorBrewer)
library(ape)
library(colorspace)

#read the similarity matrices into R 
cgmlst <- read.table("../../../../Dropbox/0 - Publications_bh/Thesis/Chapter 03 - EpiQuant/final_analysis/cgMLST_Can274.txt", header = T, sep = '\t', check.names = F)
row.names(cgmlst) <- cgmlst[,1]
cgmlst <- cgmlst[,-c(1,2)]
epi <- read.table("../../../../Dropbox/0 - Publications_bh/Thesis/Chapter 03 - EpiQuant/final_analysis/Epi_Sim_matrix(40_40_20).txt", header = T, sep = '\t', check.names = F)


#### Compute trees for calculating wallace coefficients
hc_cg <- hclust(dist.gene(cgmlst), method = 'single')
hc_epi <- hclust(dist(epi), method = 'single')

cg_clusters <- cutree(hc_cg, h = 0.05*max(hc_cg$height), order_clusters_as_data = F) 
nclus_cg <- length(hc_cg$height) - length(unique(cg_clusters))
hc_cg$height[1:(nclus_cg+1)] <- 0 
hc_cg$height <- (hc_cg$height / max(hc_cg$height))
print(length(unique(hc_cg$height)))

epi_clusters <- cutree(hc_epi, h = 0.44*max(hc_epi$height), order_clusters_as_data = F) 
nclus_epi <- length(hc_epi$height) - length(unique(epi_clusters))
hc_epi$height[1:(nclus_epi+1)] <- 0 
hc_epi$height <- (hc_epi$height / max(hc_epi$height))
print(length(unique(hc_epi$height)))

d_cg <- as.dendrogram(hc_cg)
d_epi <- as.dendrogram(hc_epi)

# plot(hc_cg)
# plot(d_epi)


labels(d_cg) <- as.character(labels(d_cg))
labels(d_epi) <- as.character(labels(d_epi))

dendo_random <- untangle_random_search(d_cg, d_epi, R = 1000, L = 2)
entanglement(dendo_random[[1]], dendo_random[[2]])


dend3 <- dendo_random
dendo2side <- untangle_step_rotate_2side(dend3[[1]], dend3[[2]], max_n_iterations = 5, print_times= T, direction = 'backward', k_seq = 35:2)
entanglement(dendo2side[[1]], dendo2side[[2]])


dend_heights1 <- heights_per_k.dendrogram(dendo_random[[1]])
dend1 <- untangle_step_rotate_1side(dend1 = dendo_random[[1]], dend2_fixed = dendo_random[[2]], L = 1, dend_heights_per_k = dend_heights1, k_seq = 2:40)
tanglegram(dend1$best_dend, dend1$dend2_fixed)

dend_heights2 <- heights_per_k.dendrogram(dend1$dend2_fixed)
dend2 <- untangle_step_rotate_1side(dend1 = dend1$dend2_fixed, dend2_fixed = dend1$best_dend, L = 1, dend_heights_per_k = dend_heights2, k_seq = 2:16)

dend4 <- untangle_step_rotate_1side(dend2$dend2_fixed, dend2_fixed = dend2$best_dend, L = 1, dend_heights_per_k = dend_heights1, k_seq = 2:40)
entanglement(dend4[[1]], dend4[[2]])
dend5 <- untangle_step_rotate_1side(dend2$best_dend, dend2_fixed = dend2$dend2_fixed, L = 1, dend_heights_per_k = dend_heights2, k_seq = 2:20)
entanglement(dend5[[1]], dend5[[2]])






num_k <- length(unique(hc_epi$height))
dend4[[2]] <-  color_branches(dend4[[2]],
                              num_k,
                              col = rainbow_hcl(num_k))
dend4[[1]] <- color_branches(dend4[[1]], k = 1, col = 'black')


col_lines_left2 <- rainbow_hcl(num_k)[cutree(dend4[[2]],
                      num_k, order_clusters_as_data = F,
                      sort_cluster_numbers = T)]


tanglegram(dend4[[2]], dend4[[1]], margin_inner = 7,
           color_lines = col_lines_left2, 
           columns_width = c(4,2,4),
           lwd = 2,
           lab.cex = 0.45, 
           cex.main = 0.95,
           main_left = "Epidemiologic Clustered Tree", 
           main_right = "cgMLST Clustered Tree")






