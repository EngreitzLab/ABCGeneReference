suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(optparse)})

main <- function() {
  ## Read in data
  message("Reading in data...")
  
  option.list <- list(
    make_option("--infile", type="character", help="input file"),
    make_option("--outfile", type="character", help="output file"))
  opt <- parse_args(OptionParser(option_list=option.list))
  
  in.file = opt$infile; out.file = opt$outfile
  
  transcripts = read.table(in.file, header=FALSE, sep="\t") %>%
    setNames(c('chr', 'start', 'end', 'gene', 'level', 'strand', 'gene_type'))

  # goal = get these columns (in order) chr, TSS start, TSS end, TSS, score (level), strand, gene, gene start, gene end
  transcripts$TSS = ifelse(transcripts$strand=="+", transcripts$start, transcripts$end)
  transcripts$start_TSS = transcripts$TSS-250
  transcripts$end_TSS = transcripts$TSS+250
  
  fullList = data.frame(chr=transcripts$chr, TSS_start=transcripts$start_TSS, TSS_end=transcripts$end_TSS, 
                        TSS=transcripts$TSS, score=transcripts$level, strand=transcripts$strand, gene=transcripts$gene, 
                        gene_start=transcripts$start, gene_end=transcripts$end)
  
  message("Writing output...")
  readr::write_tsv(fullList, file=out.file, col_names=FALSE, quote_escape=FALSE)
}

main()
