#!/home/jia.meng/Zhen/R/bin/Rscript
library(GenomicAlignments)
library(readr)
setwd('/scratch/lhcumt/work/metdb2/align_visualize_expr/Homo_sapiens/human_hg19/')
bin <- read_rds('/home/jia.meng/Xiangyu/bin.rds')
aln <- readGAlignments('p001_HEK293T_S1_SYSY_input/p001_HEK293T_S1_SYSY_input_accept_hits.bam')
rpkm <- countOverlaps(bin, aln)/((200/10^3)*(length(aln)/10^6))
write_tsv(data.frame(p001_HEK293T_S1_SYSY_input = rpkm),'/home/jia.meng/Xiangyu/rpkm/Rpkm_1.txt')
