
# in this file, settings that are specific for a run on a dataset

# gives path to output folder
pipOutFold <- "OUTPUT_FOLDER/TCGA_brca_lum_bas"

# full path (starting with /mnt/...)
# following format expected for the input
# colnames = samplesID
# rownames = geneID
# !!! geneID are expected not difficulted

# *************************************************************************************************************************
# ************************************ SETTINGS FOR 0_prepGeneData
# *************************************************************************************************************************

rnaseqDT_file <- "/mnt/ed4/marie/other_datasets/TCGA_basal_luminal/rnaseqDT_v2.Rdata"
my_sep <- "\t"
# input is Rdata or txt file ?
# TRUE if the input is Rdata
inRdata <- TRUE

# can be ensemblID, entrezID, geneSymbol
geneID_format <- "entrezID"
stopifnot(geneID_format %in% c("ensemblID", "entrezID", "geneSymbol"))

# are geneID rownames ? -> "rn" or numeric giving the column
geneID_loc <- "rn"
stopifnot(geneID_loc == "rn" | is.numeric(geneID_loc))

removeDupGeneID <- TRUE

# *************************************************************************************************************************
# ************************************ SETTINGS FOR 1_runGeneDE
# *************************************************************************************************************************

# labels for conditions
cond1 <- "luminal"
cond2 <- "basal"

# path to sampleID for each condition - should be Rdata ( ! sample1 for cond1, sample2 for cond2 ! )
sample1_file <- "/mnt/ed4/marie/other_datasets/TCGA_basal_luminal/lumID_rnaseq_v2.Rdata"
sample2_file <- "/mnt/ed4/marie/other_datasets/TCGA_basal_luminal/basalID_rnaseq_v2.Rdata"


minCpmRatio <- 20/888 


