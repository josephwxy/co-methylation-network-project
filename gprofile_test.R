library(readr)
library(tibble)
library(GSEABase)
library(dplyr)
library(igraph)
library(WGCNA)
library(multtest)
library(gProfileR)
m_dict_raw <- read_rds("m_dict_raw.rds")
m_dict_BP <- read_rds("m_dict_BP.rds")
merged_sites <- read_rds("merged_sites.rds")
new_testm <- read_rds("new_testm.rds")
ENSG <- read_rds("ENSG.rds")
load("ENID2GO.Rdata")
m_dict_BP <- arrange(m_dict_BP,ENSEMBL)


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

adj_meth <- adjmake(x = new_testm, y = new_testm, quant = 0.90, pcut = 0.05)
write_rds(adj_meth,"adj_meth.rds")

gme <- graph.adjacency(adj_meth, mode = 'undirected')
write_rds(gme,"gme.rds")

sites <- sort(colnames(new_testm))


find_nb <- function(namelist = sites, g = gme){
  l = list()
  for (i in 1:length(sites)){
    site <- sites[i]
    l[[i]] <- as_ids(neighbors(g, site))
  }
  l
}

m_dict_BP$nb <- find_nb(m_dict_BP$ENSEMBL)

m_dict_BP$neighbor <- lapply(m_dict_BP$nb, function(x) ENSG[match(x, ENSG$cid),]$gene)

write_rds(m_dict_BP, "m_dict_BP_neighbor.rds")

GO <- lapply(m_dict_BP$neighbor, function(x) gprofiler(x, organism = "hsapiens"))

write_rds(GO, "GO.rds")