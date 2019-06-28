

# Importation of packages.
#print(.libPaths())

library("BiocParallel")
register(MulticoreParam(7))
library(DESeq2)

# Parse arguments from terminal.
arguments <- commandArgs(trailingOnly=TRUE)
dock <- as.character(arguments[1])

###############################################################################

# Define functions.

ReadFileTextList <- function(path=NULL) {
  values <- scan(
    file=path,
    what="",
    sep="\n",
  )
  return(values)
}

ReadFileTextTable <- function(path=NULL, index) {
  #rownames(data) <- data$column
  data <- read.table(
    path,
    header=TRUE,
    sep="\t",
    quote="\"'",
    dec=".",
    row.names=index,
    check.names=FALSE,
  )
  return(data)
}

###############################################################################

# Define script's behavior.

ExecuteProcedure <- function(dock=NULL) {
  setwd(dock)
  print(getwd())
  path.tissues <- "./tissue/tissues.txt"
  tissues <- ReadFileTextList(path=path.tissues)
  print(tissues)
  print(tissues[3])

  # Attempt a trial comparison by differential gene expression.

  path.data.sample <- "./tissue/sets/adipose_samples.tsv"
  path.data.gene <- "./tissue/sets/adipose_genes.tsv"

  data.sample <- ReadFileTextTable(path=path.data.sample, index="sample")
  data.gene <- ReadFileTextTable(path=path.data.gene, index="gene")
  print(data.sample[1:10, ])
  print(data.gene[1:10, 1:10])

  # Initialize data in DESeq2.

  data.deseq <- DESeqDataSetFromMatrix(
    countData=data.gene,
    colData=data.sample,
    design= ~ tissue_minor,
  )
  data.deseq <- DESeq(
    data.deseq,
    parallel=TRUE,
  )
  report <- results(data.deseq)
  report.sort <- report[order(report$pvalue), ]
  print(report.sort)
  summary(report.sort)

}

ExecuteProcedure(dock=dock)
