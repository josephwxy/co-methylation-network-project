h = function(df){
  chr <- df$chromosome
  modstart = df$modStart
  modId = df$modId
  strand = df$strand
  supnum = df$supportNum
  gr  = GRanges(seqnames = chr, IRanges(start = modstart, width = 2), strand = strand)
  gr$supportNum = supnum
  gr
}
