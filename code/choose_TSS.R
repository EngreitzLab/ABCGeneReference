suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(optparse)})

option.list <- list(
  make_option("--infile", default="/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/GENCODEv29/gencode.v29.transcripts.level12.basic.protein_coding.bed", type="character", help="input file"),
  make_option("--outfile", default="/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/GENCODEv29/gencode.v29.transcripts.level12.basic.protein_coding.unique_TSS.bed", type="character", help="output file"),
  make_option("--geneOutFile", default="/oak/stanford/groups/engreitz/Users/sheth/hg38_resources/GENCODEv29/gencode.v29.transcripts.level12.basic.protein_coding.uniqe_TSS.gene_bodies.bed", type="character"))
opt <- parse_args(OptionParser(option_list=option.list))

in.file = opt$infile; out.file = opt$outfile; gene.out.file=opt$geneOutFile

TSS_stats <- function() {
  message("Reading in data...")
  ## Read in data
  transcript.list.full = read.table(gzfile(in.file), header=FALSE, sep="\t") %>%
    setNames(c('chr', 'start', 'end', 'gene', 'level', 'strand', 'gene_type'))
  transcript.list = dplyr::select(transcript.list.full, -end) %>% distinct()
  
  genes = dplyr::select(transcript.list, gene) %>% distinct()
  print(nrow(genes))

  for (i in 1:nrow(genes)){
    goi = genes$gene[i]
    sub = dplyr::filter(transcript.list, gene==goi)
    
    genes$nTSS[i]=nrow(sub)
    genes$nL1[i]=nrow(filter(sub, level==1))
    genes$nL2[i]=nrow(filter(sub, level==2))
    }
  
  message("Writing output...")
  readr::write_tsv(genes, file='gene.TSS.numbers.txt', col_names=TRUE, quote_escape=FALSE)
}

more_TSS_stats <- function() {
  ## count how many genes have !(one TSS or one level 1 TSS)-- can we manually curate?
  df = read.table(file='gene.TSS.numbers.txt', header=TRUE, sep="\t")
  n = 0
  for (i in 1:nrow(df)){
    if (df$nTSS[i]==1) {
      n=n+1
      } else
      {
        if (df$nL1[i]==1) {
          n=n+1
        } else if (df$nL1[i]==0 & df$nL2[i]==1) {
          n=n+1
        }
      }
  }
  print(19803-n) # number without easy solution
}

choose_TSS <- function() {
  transcript.list.full = read.table(gzfile(in.file), header=FALSE, sep="\t") %>%
    setNames(c('chr', 'start', 'end', 'gene', 'level', 'strand', 'gene_type'))
  transcript.list.full = transcript.list.full[order(transcript.list.full$gene),]
  transcript.list.full$TSS = as.integer(if_else(transcript.list.full$strand=="+", transcript.list.full$start, transcript.list.full$end))
  transcript.list.full = transcript.list.full[order(transcript.list.full$gene, transcript.list.full$TSS),]
  transcript.list = dplyr::select(transcript.list.full, -c(start, end)) %>% distinct()

  genes = dplyr::select(transcript.list, gene) %>% distinct()
  #print(nrow(genes))
  TSS = data.frame(chr=character(), gene=character(), level=integer(),
                    strand=character(), gene_type=character(), TSS=integer())
  no.TSS = data.frame(chr=character(),start=integer(), end=integer(), gene=character(), 
                      level=integer(), strand=character(), gene_type=character(), TSS=integer(),
                      choice=character())
  geneFile = data.frame(chr=character(), start=integer(), end=integer(), gene=character(), 
                        level=integer(), strand=character(), gene_type=character(),
                        choice=character(), TSS=integer())
  threshold = 500
  
  for (i in 1:nrow(genes)){

    goi = genes$gene[i]
    sub = filter(transcript.list, gene==goi)
    sub.full = filter(transcript.list.full, gene==goi)
    L1 = filter(sub, level==1)
    L2 = filter(sub, level==2)
    L1.full = filter(sub.full, level==1)
    L2.full = filter(sub.full, level==2)

    if (nrow(sub)==1 && !is.element(goi, TSS$gene)){ # if there is only one TSS, use that
      sub$choice="One TSS"
      TSS = rbind(TSS, sub)
      sub.full$choice = "One TSS"
      geneFile = rbind(geneFile, sub.full)
      message("Check 1")
    } 
    
    if (nrow(L1)==1 && nrow(L2>0) && !is.element(goi, TSS$gene)) { # if there are multiple TSSs, but only one L1, use that
      L1$choice = "One L1 TSS"
      TSS = rbind(TSS, L1)
      L1.full$choice = "One L1 TSS"
      geneFile = rbind(geneFile, L1.full)
      message("Check 2")
      
    }
    
    if (nrow(L1.full)>1 && !is.element(goi, TSS$gene)) { # if there are multiple L1 TSSs, but any are within "threshold" bp, take average and use that
      for (j in 1:(nrow(L1.full)-1)){
        diff = as.numeric(abs(L1.full$TSS[j]-L1.full$TSS[j+1]))
        if (diff<=threshold && !(is.element(goi, TSS$gene))){
          TSS = rbind(TSS, c(chr=L1.full$chr[j],  gene=goi, level=1, 
                             strand=L1.full$strand[j], gene_type=L1.full$gene_type[j],
                             TSS=ceiling((L1.full$TSS[j]+L1.full$TSS[j+1])/2),
                             choice="Multiple close L1 TSSs"))
          end.this = if_else(L1.full$strand[j]=="+", ceiling((L1.full$end[j]+L1.full$end[j+1])/2) , ceiling((L1.full$TSS[j]+L1.full$TSS[j+1])/2))
          start.this = if_else(L1.full$strand[j]=="+", ceiling((L1.full$TSS[j]+L1.full$TSS[j+1])/2),ceiling((L1.full$end[j]+L1.full$end[j+1])/2))
          geneFile = rbind(geneFile, c(chr=L1.full$chr[j],  start=ceiling(L1.full$TSS[j]+L1.full$TSS[j+1])/2,
                           TSS=ceiling((L1.full$TSS[j]+L1.full$TSS[j+1])/2),
                           end=end.this, start=start.this, gene=goi, level=1,
                           strand=L1.full$strand[j], gene_type=L1.full$gene_type[j],
                           TSS=ceiling((L1.full$TSS[j]+L1.full$TSS[j+1])/2),
                           choice="Multiple close L1 TSSs"))
        }
      }
      message("Check 3")
      
    }
    
    if (nrow(sub.full)>1 && !is.element(goi, TSS$gene)) { # if there are multiple TSSs, but any are within 250 bp, take average and use that
      for (j in 1:(nrow(sub.full)-1)){
        diff = as.numeric(abs(sub.full$TSS[j]-sub.full$TSS[j+1]))
        if (diff<=threshold && !is.element(goi, TSS$gene)){
          TSS = rbind(TSS, c(chr=sub.full$chr[j], gene=goi, level=min(sub.full$level[j], sub.full$level[j+1]),
                             strand=sub.full$strand[j], gene_type=sub.full$gene_type[j],
                             TSS=ceiling((sub.full$TSS[j]+sub.full$TSS[j+1])/2),
                             choice="Multiple close TSSs"))
          end.this = if_else(sub.full$strand[j]=="+", ceiling((sub.full$end[j]+sub.full$end[j+1])/2),ceiling((sub.full$start[j]+sub.full$start[j+1])/2))
          start.this = if_else(sub.full$strand[j]=="+", ceiling((sub.full$start[j]+sub.full$start[j+1])/2),ceiling((sub.full$end[j]+sub.full$end[j+1])/2))
          geneFile = rbind(geneFile, c(chr=sub.full$chr[j], start = start.this, end=end.this,
                                       gene=goi, level=min(sub.full$level[j], sub.full$level[j+1]),
                                       strand=sub.full$strand[j], gene_type=sub.full$gene_type[j],
                                       TSS=ceiling((sub.full$TSS[j]+sub.full$TSS[j+1])/2),
                                       choice="Multiple close TSSs"))
        }
      }
      message("Check 4")
      
    }
    
    if (!is.element(goi, TSS$gene)) {
      # choose longest transcript's TSS
      sub.full$transcript.length = abs(sub.full$start-sub.full$end)
      longest = filter(sub.full, transcript.length==max(sub.full$transcript.length))
      TSS = rbind(TSS, c(chr=longest$chr[1], gene=goi, level=longest$level[1],
                         strand=longest$strand[1], gene_type=longest$gene_type[1],
                         TSS=longest$TSS[1],choice="Longest transcript"))
      geneFile = rbind(geneFile, c(chr=longest$chr[1], start=longest$start[1], end = longest$end[1],
                                   gene=goi, level=longest$level[1],
                                   strand=longest$strand[1], gene_type=longest$gene_type[1],
                                   TSS=longest$TSS[1],
                                   choice="Longest transcript"))  
      message("Check 5")
      
    }
    
    
    if (!is.element(goi, TSS$gene)) {
      no.TSS = rbind(no.TSS, sub.full)
    }
    message("Check 6")
  }
    
  #print(nrow((TSS)))
  unique.genes = dplyr::select(TSS, gene) %>% distinct()
  #print(nrow(unique.genes))
  #print(head(unique.genes))
  
  TSS$TSS.start = as.numeric(TSS$TSS)-250
  TSS$TSS.end = as.numeric(TSS$TSS)+250
  TSS = TSS[, c('chr', 'TSS', 'TSS.start', 'TSS.end', 'gene', 'level', 'strand', 'gene_type', 'choice')]
    
  # convert to 500 bp regions
  readr::write_tsv(TSS, file=out.file, col_names=TRUE, quote_escape=FALSE)
  readr::write_tsv(geneFile, gene.out.file, col_names=TRUE, quote_escape=FALSE)
  #print(head(no.TSS))
  #readr::write_tsv(no.TSS, file='no.TSS.tsv', col_names=TRUE, quote_escape=FALSE)
}


#TSS_stats()
#more_TSS_stats()
choose_TSS()

