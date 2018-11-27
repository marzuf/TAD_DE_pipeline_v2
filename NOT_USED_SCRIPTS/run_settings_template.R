
# in this file, settings that are specific for a run on a dataset

# gives path to output folder
pipOutFold <- "."

# full path (starting with /mnt/...)
# following format expected for the input
# colnames = samplesID
# rownames = geneID
# !!! geneID are expected not difficulted

# *************************************************************************************************************************
# ************************************ SETTINGS FOR 0_prepGeneData
# *************************************************************************************************************************

rnaseqDT_file <- "/mnt/ed4/marie/other_datasets/GSE90749/GSE90749_genes.fpkm_tracking.combined.batch.all.txt"
my_sep <- "\t"
# input is Rdata or txt file ?
# TRUE if the input is Rdata
inRdata <- FALSE

# can be ensemblID, entrezID, geneSymbol
geneID_format <- "ensemblID"
stopifnot(geneID_format %in% c("ensemblID", "entrezID", "geneSymbol"))

# are geneID rownames ? -> "rn" or numeric giving the column
geneID_loc <- 1
stopifnot(geneID_loc == "rn" | is.numeric(geneID_loc))

# from the start, take only the genes that are not ambiguous
removeDupGeneID <- TRUE

# *************************************************************************************************************************
# ************************************ SETTINGS FOR 1_runGeneDE
# *************************************************************************************************************************

# logFC < 0 if cond2 downregulated with respect to cond1

# labels for conditions
cond1 <- ""
cond2 <- ""

# path to sampleID for each condition - should be Rdata
sample1_file <-  ""
sample2_file <- ""





