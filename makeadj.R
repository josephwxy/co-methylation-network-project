library(readr)
library(tibble)
library(topGO)
library(GSEABase)
library(dplyr)
library(igraph)
library(WGCNA)
library(multtest)

new_testm <- read_rds("new_testm.rds")
adjmake <- function(x, y, quant, pcut){
  cor <- cor(x, y, method = 'spearman')
  rawp <- corPvalueFisher(cor, 38)
  mt <- mt.rawp2adjp(rawp, proc = 'Bonferroni')
  adj <- mt$adjp[,2]
  adjp <- matrix(adj[order(mt$index)], ncol = ncol(cor))
  cor <- abs(cor)
  scc_cut <- quantile(cor, quant)
  cor[adjp > pcut | cor < scc_cut] <- 0
  cor[cor > 0] <- 1
  diag(cor) <- 0
  cor
}
adj_meth <- adjmake(x = new_testm, y = new_testm, quant = 0.8, pcut = 0.05)
write_rds(adj_meth,"adj_meth_85.rds")

#==================================================
adj_meth <- read_rds("adj_meth_85.rds")

gme <- graph.adjacency(adj_meth, mode = 'undirected')

png(filename = 'degree-distribution-85.png', width = 700, height = 550, units = 'px')
dd <- degree_distribution(gme)
plot(log10(seq(length(dd))),log10(dd), 
     xlab = 'log10(degree)', ylab = 'log10(frequency)', 
     main = 'degree distribution of co-methylation network', pch = 20, col = 'red')
avg.degree <- mean(degree(gme, V(gme)))
abline(v = avg.degree, lwd = 1, lty = 2)
dev.off()
