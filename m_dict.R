library(readr)
library(tibble)
library(topGO)
library(org.Hs.eg.db)
library(GSEABase)

m_dict_raw <- read_rds("m_dict_raw.rds")
merged_sites <- read_rds("merged_sites.rds")
new_testm <- read_rds("new_testm.rds")
load("ENID2GO.Rdata")

# m_dict_BP
m_dict_BP <- m_dict_raw
m_dict_BP$GO_exact_BP <- ENID2GO_BP[match(m_dict_BP$ENSEMBL, names(ENID2GO_BP))]
m_dict_BP <- as_tibble(m_dict_BP)
fl <- "http://geneontology.org/ontology/subsets/goslim_generic.obo"
slim <- getOBOCollection(fl)
makeslim <- function(golist , slim. = slim, c = 1){
  ans <- list()
  for (i in seq_along(golist)){
    if (length(golist[[i]]) > 0){
      mygo <- GOCollection(golist[[i]])
      temp <- as.data.frame(goSlim(mygo,slim.,'BP'))
      ans[[i]] <- rownames(temp[temp$Count > c,])
    }
    else {ans[[i]] <- NULL}
  }
  ans
}
m_dict_BP$GO_slim <- makeslim(m_dict_BP$GO_exact_BP)
write_rds(m_dict_BP, 'm_dict_BP.rds')

# m_dict_CC
m_dict_CC <- m_dict_raw
m_dict_CC$GO_exact_CC <- ENID2GO_CC[match(m_dict_CC$ENSEMBL, names(ENID2GO_CC))]
m_dict_CC <- as_tibble(m_dict_CC)
makeslim <- function(golist , slim. = slim, c = 1){
  ans <- list()
  for (i in seq_along(golist)){
    if (length(golist[[i]]) > 0){
      mygo <- GOCollection(golist[[i]])
      temp <- as.data.frame(goSlim(mygo,slim.,'CC'))
      ans[[i]] <- rownames(temp[temp$Count > c,])
    }
    else {ans[[i]] <- NULL}
  }
  ans
}
m_dict_CC$GO_slim <- makeslim(m_dict_CC$GO_exact_CC)
write_rds(m_dict_CC, 'm_dict_CC.rds')

# m_dict_MF
m_dict_MF <- m_dict_raw
m_dict_MF$GO_exact_MF <- ENID2GO_MF[match(m_dict_MF$ENSEMBL, names(ENID2GO_MF))]
m_dict_MF <- as_tibble(m_dict_MF)
makeslim <- function(golist , slim. = slim, c = 1){
  ans <- list()
  for (i in seq_along(golist)){
    if (length(golist[[i]]) > 0){
      mygo <- GOCollection(golist[[i]])
      temp <- as.data.frame(goSlim(mygo,slim.,'MF'))
      ans[[i]] <- rownames(temp[temp$Count > c,])
    }
    else {ans[[i]] <- NULL}
  }
  ans
}
m_dict_MF$GO_slim <- makeslim(m_dict_MF$GO_exact_MF)
write_rds(m_dict_MF, 'm_dict_MF.rds')