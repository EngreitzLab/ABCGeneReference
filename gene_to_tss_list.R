suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(optparse)})

main <- function() {
  message("Reading in data...")
  ## Read in data
  option.list <- list(
    make_option("--infile", type="character", help="input file"),
    make_option("--outfile", type="character", help="output file"))
  opt <- parse_args(OptionParser(option_list=option.list))
  
  in.file = opt$infile; out.file = opt$outfile
  
  geneList = read.table(gzfile(in.file), header=FALSE, sep="\t") %>%
    setNames(c('chr', 'start', 'end', 'gene', 'level', 'strand', 'gene_type'))
  
  TSS = ifelse(geneList$strand=="+", geneList$start, geneList$end)
  TSSList = geneList
  TSSList$start = TSS-250
  TSSList$end = TSS+250
  
  message("Writing output...")
  readr::write_tsv(TSSList, file=out.file, col_names=FALSE, quote_escape=FALSE)
}

main()
