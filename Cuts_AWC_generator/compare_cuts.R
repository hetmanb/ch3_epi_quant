#### CGF Analysis ####
library(ape)
library(dendextend)
library(gplots)
library(reshape)
library(plyr)

epi <- read.table("../../../../Dropbox/0 - Publications_bh/Thesis/Chapter 03 - EpiQuant/input_datas_274/epi_matrix_176(nohuman).txt", header = T, sep = '\t', check.names = F)
epi_1 <- abs(epi-1)
hc_epi <- hclust(as.dist(epi_1), method = 'single')

cgmlst <- read.table("../../../../Dropbox//0 - Publications_bh/Thesis/Chapter 03 - EpiQuant/input_datas_274/cgMLST_Can274.txt", header = T, sep = '\t', check.names = F)
row.names(cgmlst) <- cgmlst[,1]
cgmlst <- cgmlst[,-c(1:2)]
d_cgmlst <- dist.gene(cgmlst, method = 'percentage')
# d_cgmlst <- d_cgmlst/ncol(cgmlst)
hc_cgmlst <- hclust(d_cgmlst, method = 'single')


mlst <- read.table("../../../../Dropbox//0 - Publications_bh/Thesis/Chapter 03 - EpiQuant/input_datas_274/mlst_274.txt", header = T, sep = '\t', check.names = F)
row.names(mlst) <- mlst[,1]
mlst <- mlst[,-c(1:2)]
d_mlst <- dist.gene(mlst, method = 'percentage')
# d_mlst <- d_mlst/ncol(mlst)
hc_mlst <- hclust(d_mlst, method = 'single')


cgf <- read.table("../../../../Dropbox//0 - Publications_bh/Thesis/Chapter 03 - EpiQuant/input_datas_274/cgf_274.txt", header = T, sep = '\t')
row.names(cgf) <- cgf[,1]
cgf <- cgf[,-c(1:2)]
d_cgf <- dist.gene(cgf, method = 'percentage' )
hc_cgf <- hclust(d_cgf, method = 'single')


 


######### Testing the order of cuts #######
# cgf2 <- cgf[1:20,]
# d2_cgf <- dist.gene(cgf2, method = 'pairwise')
# d2_cgf <- d2_cgf/ncol(cgf2)
# hc2_cgf <- hclust(d2_cgf, method = 'single')
# seq <- seq(0,0.1,0.01)
# cgf2_cuts <- data.frame(matrix(ncol = length(seq), nrow = 0))
# colnames(cgf2_cuts) <- seq
# 
# for(i in 1:length(seq)){
#   x <- cutree(hc2_cgf, #k = seq[i])
#               h= seq[i]*max(hc2_cgf$height),
#               order_clusters_as_data = T)
#   cgf2_cuts[1:20, i] <- x
# #   y$cgf2[i] <- max(x)
# }
# row.names(cgf2_cuts) <- hc2_cgf$labels[hc2_cgf$order]
# 
# hc2_dend <- as.dendrogram(hc2_cgf)
# plot(hc2_dend)
#########  --- END -- Testing the order of cuts #######







seq <- seq(0,0.9,0.01)
# seq <- seq(0.75, 0.25, -0.10)
y <- data.frame(seq)

cgmlst_cuts <- data.frame(matrix(ncol = length(seq), nrow = 0))
colnames(cgmlst_cuts) <- seq
cgf_cuts <- data.frame(matrix(ncol = length(seq), nrow = 0))
colnames(cgf_cuts) <- seq
epicuts <- data.frame(matrix(ncol = length(seq), nrow = 0))
colnames(epicuts) <- seq



for(i in 1:length(seq)){
  x <- cutree(hc_cgmlst, #k = seq[i])
              h= seq[i]*max(hc_cgmlst$height),
              order_clusters_as_data = T)
  cgmlst_cuts[1:274, i] <- x
    y$cgmlst[i] <- max(x)
}
row.names(cgmlst_cuts) <- hc_cgmlst$labels[hc_cgmlst$order]

for(i in 1:length(seq)){
  x <- cutree(hc_cgf, #k = seq[i])
                h= seq[i]*max(hc_cgf$height),
                order_clusters_as_data = T)
  cgf_cuts[1:274, i] <- x
  y$cgf[i] <- max(x)
}
row.names(cgf_cuts) <- hc_cgf$labels[hc_cgf$order]

for(i in 1:length(seq)){
  x <- cutree(hc_epi, #k = seq[i])
                h= seq[i]*max(hc_epi$height),
                order_clusters_as_data = T)
  epicuts[1:176, i] <- x
  y$epi[i] <- max(x)
}
row.names(epicuts) <- hc_epi$labels[hc_epi$order]



# all_cuts_data <- join_all(list(epicuts, cgf_cuts, cgmlst_cuts), by = 0) 
cgf_cuts_data <- merge(epicuts, cgf_cuts, by = 0) 
colnames(cgf_cuts_data) <- gsub("x", "epi", colnames(cgf_cuts_data))
colnames(cgf_cuts_data) <- gsub("y", "cgf", colnames(cgf_cuts_data))

cgMLST_cuts_data <- merge(cgf_cuts, cgmlst_cuts, by = 0) 
colnames(cgMLST_cuts_data) <- gsub("x", "cgf", colnames(cgMLST_cuts_data))
colnames(cgMLST_cuts_data) <- gsub("y", "cgMLST", colnames(cgMLST_cuts_data))




write.table(cgf_cuts_data, "thesis_anal/epi_vs_cgf_AWC.txt", sep = '\t', row.names = T, col.names = T)
write.table(cgMLST_cuts_data, "thesis_anal/cgf_vs_cgMLST_singlelinkage.txt", sep = '\t', row.names = T, col.names = T)
