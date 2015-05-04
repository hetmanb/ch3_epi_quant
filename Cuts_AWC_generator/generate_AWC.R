source("wallace_helper.R")

awc_table <- read.table("../../../0 - Publications_bh/Thesis/Chapter 03 - EpiQuant/input_datas_274/cg", sep = '\t', header = T, check.names=F)
row.names(awc_table) <- awc_table[,1]
awc_table <- awc_table[,-1]
# awc_table <- awc_table[,-(6:11)]

# awc_table <- all_cuts_data
awc_table <- awc_table[,-135]

d <- expand.grid(1:ncol(awc_table), 1:ncol(awc_table))
results_table <- data.frame()
for(i in 1:length(d[,1])){
  print(i)
  a <- d[i,1]
  b <- d[i,2]
  
  aw <- list(adj_wallace(awc_table[,a], awc_table[,b]))
  aw <- aw[[1]]
  awa <- aw$Adjusted_Wallace_A_vs_B
  awb <- aw$Adjusted_Wallace_B_vs_A
  results_table[i,1] <- awa
  results_table[i,2] <- awb
}

col <- colnames(awc_table)
for(i in 1:length(d[,1])){  
results_table[i,3] <- col[d[i,1]]
results_table[i,4] <- col[d[i,2]]
}
results_table <- results_table[,c(3,4,1,2)]
colnames(results_table) <- c("Type_A","Type_B","AWC_AB","AWC_BA")

write.table(results_table, "../../../0 - Publications_bh/Thesis/Chapter 03 - EpiQuant/input_datas_274/cgmlst_RMLST_goe_AWC.txt", sep = '\t', row.names = F)
