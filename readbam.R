rm rpkm.sh

cat >>rpkm.sh<<EOF__
#!/home/jia.meng/Zhen/R/bin/Rscript
library(GenomicAlignments)
library(readr)

folders <- readLines("folders.txt")
foldernames <- gsub("^.*\\\/", "", folders)

bamfiles <- paste0(folders, "/", foldernames, "_accept_hits.bam")
bin <- read_rds("Xiangyu/bin.rds")
rpkmgr <- bin

for (i in seq_along(bamfiles)) 
  aln <- readGAlignments(bamfiles[i])
  rpkm <- countOverlaps(bin, aln)/((200/10^3)*(length(aln)/10^6))
  df <- as.data.frame(mcols(rpkmgr))
  df[[foldernames[i]]] <- rpkm 
  mcols(rpkmgr) <- df
  print(paste(folders[i], "finished"))
}

write_rds(rpkmgr, "rpkmgr.rds")
EOF__

chmod u+x rpkm.sh

echo "./rpkm.sh" | qsub -cwd -q bigmem.q -l bigmem -m abe -M 496691875@qq.com -o rpkm.out -N rpkm


