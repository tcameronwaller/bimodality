

# Importation of packages.

#library(DESeq2)

# Parse arguments from terminal.
arguments <- commandArgs(trailingOnly=TRUE)
dock <- as.character(arguments[1])




# Define script's behavior.

ExecuteProcedure <- function(dock=NULL){
    print("Hello Cameron")
    print(dock)

    print(.libPaths())

    print(browseVignettes("DESeq2"))
}

ExecuteProcedure(dock=dock)
