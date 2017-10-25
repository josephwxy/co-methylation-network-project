library(readr)
library(tibble)
library(topGO)
library(GSEABase)
library(dplyr)
library(igraph)
library(WGCNA)
library(multtest)
library(gProfileR)
library(AnnotationDbi)
library(future)
library(clusterProfiler)
library(preprocessCore)
library(AnnotationDbi)

m_dict_raw <- read_rds("m_dict_raw.rds")
merged_site <- read_rds("merged_sites.rds")
new_testm <- read_rds("new_testm.rds")

new_testm <- new_testm[match(sort(rownames(new_testm)),rownames(new_testm)),]


a <- normalize.quantiles(t(new_testm))
a <- t(a)
colnames(a) <- colnames(new_testm)
rownames(a) <- rownames(new_testm)

new_testm <- a
write_rds(new_testm, "new_testm.rds")
# ENID2GO
ENID2GO_BP <- inverseList(annFUN.org("BP", mapping = "org.Hs.eg.db", ID = "ensembl"))
ENID2GO_CC <- inverseList(annFUN.org("CC", mapping = "org.Hs.eg.db", ID = "ensembl"))
ENID2GO_MF <- inverseList(annFUN.org("MF", mapping = "org.Hs.eg.db", ID = "ensembl"))
save(ENID2GO_BP,ENID2GO_CC,ENID2GO_MF,file = "ENID2GO.Rdata")
load("ENID2GO.Rdata")

merged_sites_unlist = merged_sites[rep(1:nrow(merged_sites),elementNROWS(merged_sites$members)),]
merged_sites_unlist$members_unlist = unlist(merged_sites$members)
merged_sites <- merged_sites_unlist[-3]
merged_sites <- arrange(merged_sites, as.numeric(gsub("m6A_site_", "", merged_sites$members_unlist)))
merged_sites <- merged_sites[-2]
t <- split(merged_sites,merged_sites$cid)
modName <- sapply(t, function(a) a[,3])
names(modName) <- names(t)
cid <- names(t)

ensg <- lapply(split(merged_sites,merged_sites$cid),function(x) unique(x$gene))
a <- NULL
for (i in 1:length(ensg)){
  a <- append(a, ensg[[i]])
}
ENSG <- as.tibble(cbind(names(ensg),a))
colnames(ENSG) <- c("cid","gene")
write_rds(ENSG, "ENSG.rds")

ENID <- ENSG[match(cid,ENSG$cid),]$gene
m_dict <- cbind(cid,ENID)
m_dict <- as.tibble(m_dict)
m_dict$modName <- modName

# delete the rows whose neighbors contain "LRG_"
m_dict <- m_dict[-(20175:20497),]

# mdict_BP
m_dict_BP <- m_dict
m_dict_BP$GO_exact_BP <- ENID2GO_BP[match(m_dict_BP$ENID, names(ENID2GO_BP))]
write_rds(m_dict_BP, 'm_dict_BP.rds')

# mdict_CC
m_dict_CC <- m_dict
m_dict_CC$GO_exact_CC <- ENID2GO_CC[match(m_dict_CC$ENID, names(ENID2GO_CC))]
write_rds(m_dict_CC, 'm_dict_CC.rds')

# mdict_MF
m_dict_MF <- m_dict
m_dict_MF$GO_exact_MF <- ENID2GO_MF[match(m_dict_MF$ENID, names(ENID2GO_MF))]
write_rds(m_dict_MF, 'm_dict_MF.rds')

# build co-methylation network
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
adj_meth <- adjmake(x = new_testm, y = new_testm, quant = 0.80, pcut = 0.05)
write_rds(adj_meth,"adj_meth.rds")
gme <- graph.adjacency(adj_meth, mode = 'undirected')
write_rds(gme,"gme.rds")

# degree distribution of co-methylation network
png(filename = 'degree-distribution.png', width = 700, height = 550, units = 'px')
dd <- degree_distribution(gme)
plot(log10(seq(length(dd))),log10(dd), 
     xlab = 'log10(degree)', ylab = 'log10(frequency)', 
     main = 'degree distribution of co-methylation network', pch = 20, col = 'red')
avg.degree <- mean(degree(gme, V(gme)))
abline(v = avg.degree, lwd = 1, lty = 2)
dev.off()

# find neighbors
have <- names(ENID2GO_BP)
sites <- colnames(new_testm)
find_nb <- function(namelist = sites, g = gme){
  l = list()
  for (i in 1:length(sites)){
    site <- sites[i]
    l[[i]] <- as_ids(neighbors(g, site))
  }
  l
}
m_dict_BP$nb <- find_nb(m_dict_BP$cid)

m_dict_BP$nb_ENID <- lapply(m_dict_BP$nb, function(x) unique(ENSG[match(x, ENSG$cid),]$gene))
write_rds(m_dict_BP,"m_dict_BP.rds")

# only select the gene which have equal and more than  neighbors
m_dict_BP_trim <- filter(m_dict_BP, lengths(nb_ENID) > 4)

# m_dict_BP_trim$entrez <- lapply(m_dict_BP_trim$nb_ENID, function(x) 
#   bitr(x, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db"))
m_dict_BP_trim$entrez <- lapply(m_dict_BP_trim$nb_ENID, function(a) 
  unique(na.omit(AnnotationDbi::select(org.Hs.eg.db, keys = a, columns ="ENTREZID", keytype="ENSEMBL")$ENTREZID)))
write_rds(m_dict_BP_trim,"m_dict_BP_trim.rds")

background <- unique(unlist(m_dict_BP_trim$entrez))
# BP annotation on neighbors
i <- 0
GO_BP_predict <- function(a){
  i <<- i + 1; print(i)
  GO_pre <- enrichGO(gene = a, OrgDb = "org.Hs.eg.db", ont = "BP", pvalueCutoff = 0.01, universe = background, readable = TRUE)
  if (nrow(GO_pre) == 0){
    GO_pre <- NULL
  }
  else {
    GO_pre <- GO_pre[,c(1:6,9)]
  }
}
BP <- lapply(m_dict_BP_trim$entrez, GO_BP_predict)
write_rds(BP, "BP_term.rds")
a <- enrichGO(gene = m_dict_BP_trim$entrez[[1]], OrgDb = "org.Hs.eg.db", ont = "BP", pvalueCutoff = 0.01, readable = TRUE)

#a <- lapply(m_dict_BP_trim$entrez, function(a) enrichGO(gene = a, OrgDb = "org.Hs.eg.db", ont = "BP", readable = TRUE)[,c(1:6,9)])
#unlist(lapply(m_dict_BP_trim$entrez, function(a) length(a)))

#=====Proxy Error=====  
#gpro <- function(a) {
#  if (length(a) == 0){
#    goresult <- NULL
#  }
#  else{
#    goresult <- gprofiler(a, organism = "hsapiens")
#  }
#}
#GO <- lapply(m_dict_BP$nb_ENID, gpro)
#=====================

# module-based cluster
# =======MCL========
adj <- exportNetworkToCytoscape(adj_meth)$edgeData
write_tsv(adj,"adj.tsv")
adj <- read_tsv("adj.tsv")
adj <- adj[,1:3]
adj <- as.matrix(adj)
for (i in 1:nrow(adj)){
  adj [i,] <- paste0(adj[i,1]," ",adj[i,2]," ",adj[i,3])
}
adj <- as.data.frame(adj[,1])
write_tsv(adj,"adjfile.tsv")

# ============ mcl clustering (in terminal) ==================
# mcl adjfile.tsv -I 1.8 --abc -o out.adjfile1.8

adj1.8 <- as.matrix(read_csv("out.adjfile1.8", col_names = FALSE, na = c("NA")))
mclmodule <- NULL
for (i in 1:length(adj1.8)){
  mclmodule[[i]] <- strsplit(adj1.8[i],"\t")
}

mcl_module <- NULL
for (i in 1:length(mclmodule)){
  if (elementNROWS(mclmodule[[i]]) > 4){
    mcl_module <- append(mcl_module, mclmodule[[i]])
  }
}

cluster <- NULL
meth <- NULL
for (j in 1:length(mcl_module)){
  cluster <- append(cluster, rep(j,length(mcl_module[[j]])))
  meth <- append(meth, mcl_module[[j]])
}
  
module_num <- c(1:length(unique(cluster)))
module_size <- unlist(lapply(mcl_module, function(a) length(a)))
module_dict <- as.tibble(cbind(module_num, module_size))
module_dict$cid <- mcl_module
module_dict$ENID <- lapply(module_dict$cid, function(x) unique(ENSG[match(x, ENSG$cid),]$gene))
module_dict$entrez <- lapply(module_dict$ENID, function(a) 
  unique(na.omit(AnnotationDbi::select(org.Hs.eg.db, keys = a, columns ="ENTREZID", keytype="ENSEMBL")$ENTREZID)))
module_dict_trim <- filter(module_dict, lengths(entrez) > 3)
module_background <- unique(unlist(module_dict_trim$entrez))

enrichGO(gene = module_dict_trim$entrez[[6]], OrgDb = "org.Hs.eg.db", ont = "BP", pvalueCutoff = 0.01, universe = module_background, readable = TRUE)

i <- 0
GO_predict <- function(a){
  i <<- i + 1; print(i)
  GO_pre <- enrichGO(gene = a, OrgDb = "org.Hs.eg.db", ont = "BP", pvalueCutoff = 0.01, universe = module_background, readable = TRUE)
  if (nrow(GO_pre) == 0){
    GO_pre <- NULL
  }
  else {
    GO_pre <- GO_pre[,c(1:6,9)]
  }
}
a <- lapply(module_dict_trim$entrez, GO_predict)

# =====================topGO=========================
ifhave <- mapply(function(cv){any(cv %in% have)}, cv = m_dict_BP$nb_ENID)
eval <- m_dict_BP[ifhave,]
eval

test.stat <- new('classicCount', testStatistic = GOFisherTest, name = 'fisher test')
i <- 0
annGO <- function(ids, geneID2GO. = ENID2GO_BP){
  i <<- i+1; print(i)
  geneList <- factor(as.integer(have %in% ids))
  names(geneList) <- have
  GOdata <- new("topGOdata", ontology = 'BP', 
                allGenes = geneList, annot = annFUN.gene2GO,
                gene2GO = geneID2GO.)
  resultFisher <- getSigGroups(GOdata, test.stat)
  t <- sort(score(resultFisher))
  t <- t[t < 0.1] 
  t
}

eval$GO_predict <- mapply(annGO, ids = eval$nb_ENID)