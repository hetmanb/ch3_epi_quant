library(dplyr)
library(ggplot2)
library(RColorBrewer)
df <- read.table('../../../0 - Publications_bh/Thesis/Chapter 03 - EpiQuant/input_datas_274/cgmlst_RMLST_goe_AWC.txt', sep = '\t', header  = T)
# df 


awc <- ggplot(df, aes((Nclus_A), (AWC_AB)))
awc <- awc + geom_point(aes(colour = (Nclus_B))) 
awc <- awc + scale_color_gradient2(low = "blue", mid = "red", high = "green", midpoint = 30)
awc <- awc + geom_vline(xintercept = c(112,83,60,47,37,25,14), colour = 'purple', linetype = 'dotted')
awc <- awc + annotate("text", x = c(112,83,60,47,37,25,14), y = 0.05,
                      angle= 90, size = 3, color = 'black', vjust = -0.5,
                      label=c("cgMLST_95%", "cgMLST_90%", "cgMLST_80%", "cgMLST_70%", "cgMLST_60%", "cgMLST_50%", "cgMLST_40%"))
awc <- awc + fte_theme()
awc <- awc + scale_x_reverse()
awc <- awc + labs(title = "Adjusted Wallace Concordance of cgMLST versus MLST Clustering",
                  x = "# cgMLST Clusters",
                  #                   x = "cgMLST Threshold (%)",
                  y = "Adjusted Wallace cgMLST to MLST",
                  colour = "# MLST Clusters")
awc
