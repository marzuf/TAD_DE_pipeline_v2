startTime <- Sys.time()


SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")


################  OUTPUT
# - rna_geneList.Rdata
# - rna_rnaseqDT.Rdata
# - rna_madnorm_rnaseqDT.Rdata
# - pipeline_geneList.Rdata
# - pipeline_regionList.Rdata
# (+ some pictures)
# - update 04.03.2018: for compatibility with updated 5_, output rna_fpkmDT.Rdata = rna_rnaseqDT.Rdata
################################################################################

# => 04.01 version: do not favor TAD #if(any(grepl("_TAD", tmp_g2t_DT$region))) 
# IN THIS VERSION -> KEEP THE GENE THAT HAS HIGHEST OVERLAP WITH A REGION
# if mapping position for none of the gene is available, take a random one
# otherwise take the one with max overlap (if more than 1 with equal overlap, a random of those ones)

pipScriptDir <- paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2")

# because I sample random genes if there are ties in their # bp overlap
set.seed(20171223)

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1)
settingF <- args[1]
stopifnot(file.exists(settingF))

script_name <- "0EZH2v2clean_prepGeneData"
stopifnot(file.exists(paste0(pipScriptDir, "/", script_name, ".R")))
cat(paste0("> START ", script_name,  "\n"))
# -> rename to 0_ so that saved data can be used in following scripts
script_name <- "0_prepGeneData"

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")
# cat(paste0("setDir = ", setDir, "\n"))
source("main_settings.R") # setDir is the main_settings not in run_settings
#### IMPORTANT !! run_settings.R SHOULD BE AFTER main_settings.R TO OVERWRITE DEFAULT FILES !!!
source(settingF)
source(paste0(pipScriptDir, "/", "TAD_DE_utils.R"))
suppressPackageStartupMessages(library(edgeR, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(ggplot2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(ggpubr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

# if microarray was not set in the settings file -> by default set to  FALSE
if(!exists("microarray")) microarray <- FALSE

stopifnot(microarray)

cat(paste0("... RNA DATA COME FROM MICROARRAY: ", as.character(microarray), "\n"))

plotType <- "svg"

##### IN THE VERSION A_ INSTEAD OF 0_
# take as input a specific gene2tad file
# 

######### 0_prepGene
# ->  according to the user settings, create a geneList object 
# ... with the same genes as rownames in input data
# ... geneList["geneName"] = entrezID
# FOR DEBUG:
# setDir <- "/media/electron/"
# TADpos_file <- paste0(setDir, "/mnt/ed4/marie/scripts/EZH2_final_MAPQ/06_12_50kb_MAPQFILTER/consensus/gene2tad/KARPAS_DMSO_LY19WT_DMSO_c1_r100000_v0_w-1/TopDom/TopDom_KARPAS_DMSO_LY19WT_DMSO_c1_r100000_v0_w-1_TADposDT.txt")
# gene2tadDT_file <- paste0(setDir, "/mnt/ed4/marie/scripts/EZH2_final_MAPQ/06_12_50kb_MAPQFILTER/consensus/gene2tad/KARPAS_DMSO_LY19WT_DMSO_c1_r100000_v0_w-1/TopDom/TopDom_KARPAS_DMSO_LY19WT_DMSO_c1_r100000_v0_w-1_gene2tadDT.txt")
# entrezDT_file <-  paste0(setDir, "/mnt/ed4/marie/scripts/EZH2_final_MAPQ/06_12_50kb_MAPQFILTER/consensus/gene2tad/KARPAS_DMSO_LY19WT_DMSO_c1_r100000_v0_w-1/TopDom/TopDom_KARPAS_DMSO_LY19WT_DMSO_c1_r100000_v0_w-1_entrezDT.txt")
# rnaseqDT_file <- "/mnt/ed4/marie/scripts/EZH2_final_MAPQ/GENE_DATA_MZ/CL/cell.lines.byGene_noDup_MZ.RData"
# my_sep <- "\t"
# inRdata <- TRUE
# geneID_format <- "entrezID"
# stopifnot(geneID_format %in% c("ensemblID", "entrezID", "geneSymbol"))
# geneID_loc <- "rn"
# removeDupGeneID <- TRUE
# cond1 <- "wt"
# cond2 <- "mut"
# sample1_file <- "/mnt/ed4/marie/scripts/EZH2_final_MAPQ/gene_expression/CL/wt_ID.Rdata"
# sample2_file <- "/mnt/ed4/marie/scripts/EZH2_final_MAPQ/gene_expression/CL/mut_ID.Rdata"
#####################

# create the directories
curr_outFold <- paste0(pipOutFold, "/", script_name)
system(paste0("mkdir -p ", curr_outFold))

pipLogFile <- paste0(pipOutFold, "/", script_name, "_logFile.txt")
system(paste0("rm -f ", pipLogFile))

# ADDED 27.11.2018 to check using other files
txt <- paste0("inputDataType\t=\t", inputDataType, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0("gene2tadDT_file\t=\t", gene2tadDT_file, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0("TADpos_file\t=\t", TADpos_file, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0("settingF\t=\t", settingF, "\n")
printAndLog(txt, pipLogFile)

# load the expression data
if(inRdata) {
  x <- load(paste0(setDir, rnaseqDT_file))
  rnaseqDT <- eval(parse(text=x))
} else {
  if(geneID_loc == "rn"){
    rnaseqDT <- read.delim(paste0(setDir, rnaseqDT_file), sep=my_sep, header=T, row.names = 1)
  } else {
    rnaseqDT <- read.delim(paste0(setDir, rnaseqDT_file), sep=my_sep, header=T)
    rownames(rnaseqDT) <- rnaseqDT[, geneID_loc]
    rnaseqDT <- rnaseqDT[,-geneID_loc]
  }
}
# rnaseqDT[1:5,1:5]

gene2tadDT <- read.delim(gene2tadDT_file, header=F, col.names = c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = F)
gene2tadDT$entrezID <- as.character(gene2tadDT$entrezID)


########################################################################################################################## 
########################################################################################################################## ADDITIONAL STEPS HERE FOR EZH2 <<<<--- START
##########################################################################################################################

cat("... PREPROCESSING: transform rows with multiple gene IDs\n")

# the expression DT still contains unique expression values rows (with some rownames = many genes)
# because I need them to be unique for the gene expression DE analysis
# but the madnorm DT contains 1 row per gene per TAD (to not have duplicate)

# before the pipeline, the expression data have already be cleaned
# so that the gene ID they contain are not duplicated
# however some genes from a given row can still map to the same TAD
# which means that I will use a duplicated row for a given TAD => I need to change this for the MAD norm 
# that I will use for the intra-TAD correlation

# => v2: because during the reshuffling I cannot guarantee that duplicate rows cannot end 
# in the same TAD, so keep the rows unique
# so from the many genes, keep only the one with available mapping positions assigned to the largest number of TADs
# I know they all have mapping positions because this was done upstream
# output: expression DT with everywhere 1 gene/row + available mapping positions
rnaseqDT_tmp <- data.frame(rnaseqDT)
rnaseqDT_tmp$rowNbr <- 1:nrow(rnaseqDT_tmp)
singleGene_rnaseqDT <- rnaseqDT_tmp[ ! grepl("///", rownames(rnaseqDT_tmp)),]
manyGenes_rnaseqDT <- rnaseqDT_tmp[  grepl("///", rownames(rnaseqDT_tmp)),]

txt <- paste0(toupper(script_name), "> Number of rows with unique gene ID: ", nrow(singleGene_rnaseqDT), "\n")
printAndLog(txt, pipLogFile)
txt <- paste0(toupper(script_name), "> Number of rows with multiple gene IDs: ", nrow(manyGenes_rnaseqDT), "\n")
printAndLog(txt, pipLogFile)

#  FOR EZH2: processed upstream: no duplicates (but still need to remove then that would induce duplicate rows for a given TAD)
all_genes_init <- unlist(sapply(rownames(rnaseqDT_tmp), function(x) unlist(strsplit(x,split="///"))))
txt <- paste0(toupper(script_name), "> Total number of gene IDs in rownames: ", length(all_genes_init), "\n")
printAndLog(txt, pipLogFile)
stopifnot(!any(duplicated(all_genes_init)))
stopifnot(length(unique(singleGene_rnaseqDT$rowNbr)) + length(unique(manyGenes_rnaseqDT$rowNbr)) == length(unique(rnaseqDT_tmp$rowNbr)))

### in v2 version: where there are multiple genes, retain the 1st one that maps to a TAD
### in v2b version: where there are multiple genes, retain the one that maps to the largest number of TADs
# in the last version -> assignment is unambiguous (1 gene maps only to 1 TAD)
# => CHANGE 04.01 => DO NOT FAVOR TADs !

TADposDT <- read.delim(TADpos_file, header=F, col.names = c("chromo", "region", "TAD_start", "TAD_end"), stringsAsFactors = F)
unlist_genes <- unlist(sapply(rownames(manyGenes_rnaseqDT), function(x) unlist(strsplit(x,split="///"))))
sub_g2t_DT <- left_join(gene2tadDT[gene2tadDT$entrezID %in% unlist_genes,], TADposDT, by=c("chromo", "region"))
sub_g2t_DT$overlap <- pmin(sub_g2t_DT$end, sub_g2t_DT$TAD_end) - pmax(sub_g2t_DT$start, sub_g2t_DT$TAD_start)
sub_g2t_DT$entrezID <- as.character(sub_g2t_DT$entrezID)

one_of_many_rnaseqDT <- foreach(i = 1:dim(manyGenes_rnaseqDT)[1], .combine='rbind') %dopar% {
  gene_list <- rownames(manyGenes_rnaseqDT)[i]
  all_genes <-  unlist(strsplit(gene_list,split="///"))
  tmp_g2t_DT <- sub_g2t_DT[sub_g2t_DT$entrezID %in% all_genes,]
  # => 04.01 version: do not favor TAD #if(any(grepl("_TAD", tmp_g2t_DT$region))) 
  # IN THIS VERSION -> KEEP THE GENE THAT HAS HIGHEST OVERLAP WITH A REGION
  # if mapping position for none of the gene is available, take a random one
  # otherwise take the one with max overlap (if more than 1 with equal overlap, a random of those ones)
  if(nrow(tmp_g2t_DT) == 0) {
    select_gene <- sample(x = all_genes, size=1)
  } else{
    stopifnot(!any(is.na(tmp_g2t_DT$overlap)))
    select_genes <- tmp_g2t_DT$entrezID[tmp_g2t_DT$overlap == max(tmp_g2t_DT$overlap)]  
    select_gene <- sample(x = select_genes, size=1)
  }
  stopifnot(select_gene %in% all_genes)
  newDT <- manyGenes_rnaseqDT[i,,drop=FALSE]
  newDT$genes <- select_gene
  newDT
}

stopifnot(all(manyGenes_rnaseqDT$genes %in% unlist_genes))
stopifnot(nrow(one_of_many_rnaseqDT) == nrow(manyGenes_rnaseqDT))
stopifnot(!any(duplicated(one_of_many_rnaseqDT$genes)))  # this should never be the case, because filtered upstream
rownames(one_of_many_rnaseqDT) <- one_of_many_rnaseqDT$genes
stopifnot(!any(grepl("///", rownames(one_of_many_rnaseqDT))))
one_of_many_rnaseqDT$genes <- NULL
stopifnot(nrow(one_of_many_rnaseqDT) == nrow(manyGenes_rnaseqDT))
stopifnot(length(unique(singleGene_rnaseqDT$rowNbr)) + length(unique(one_of_many_rnaseqDT$rowNbr)) == length(unique(rnaseqDT_tmp$rowNbr)))
# I should not have new expression row numbers
stopifnot(all(c(singleGene_rnaseqDT$rowNbr, one_of_many_rnaseqDT$rowNbr) %in% c(singleGene_rnaseqDT$rowNbr, manyGenes_rnaseqDT$rowNbr)))
stopifnot(!any(one_of_many_rnaseqDT$rowNbr %in% singleGene_rnaseqDT$rowNbr))  # this should never be the case because duplicated filtered upstream
stopifnot(!any(rownames(one_of_many_rnaseqDT) %in% rownames(singleGene_rnaseqDT)))  # this should never be the case because duplicated filtered upstream

# merge the 2 data frames
all_DT <- rbind(one_of_many_rnaseqDT, singleGene_rnaseqDT)
txt <- paste0(toupper(script_name), "> Total number of rows from expression table with one gene/row everywhere: ", nrow(all_DT), "\n")
printAndLog(txt, pipLogFile)
all_genes_after <- unlist(sapply(rownames(all_DT), function(x) unlist(strsplit(x,split="///"))))
txt <- paste0(toupper(script_name), ">  Total number of gene IDs in rownames after one-from-many gene selection:", length(all_genes_after), "\n")
printAndLog(txt, pipLogFile)
stopifnot(nrow(all_DT) == length(all_genes_after))
stopifnot(!any(duplicated(all_DT$rowNbr)))
stopifnot(!any(grepl("///", rownames(all_DT))))
all_DT <- all_DT[order(all_DT$rowNbr),]
all_DT$rowNbr <- NULL
stopifnot(!any(duplicated(all_DT)))

rm_genes_write <- all_genes_init[! all_genes_init %in% all_genes_after]
cat(paste0(rm_genes_write, collapse="\n"), file = paste0(curr_outFold, "/", "log_rm_genes_rmFromMany.txt"), append=FALSE)

########################################################################################
# FILTER 0: DISCARD GENES WITHOUT MAPPING POSITIONS
########################################################################################

# --> THIS WILL NOT WORK FOR EZH2 AS ROWNAMNES AND GENE LIST DO NOT CORRESPOND  -->>> this should be ok for ezh2v2
# in the new version, I systematically discard the ambiguous ones
cat("... FILTER0: discard genes without mapping positions\n")

# take only the genes that I can map 
cat("... take only the ones with mapping positions:\n")
txt <- paste0(toupper(script_name), ">  Total number of rows before mapping position filter: ", nrow(all_DT), "\n")
printAndLog(txt, pipLogFile)
rm_genes_write <- rownames(all_DT)[!rownames(all_DT) %in% gene2tadDT$entrezID]
cat(paste0(rm_genes_write, collapse="\n"), file = paste0(curr_outFold, "/", "log_rm_genes_notMapping.txt"), append=FALSE)
all_DT <- all_DT[rownames(all_DT) %in% gene2tadDT$entrezID,]
txt <- paste0(toupper(script_name), ">  Total number of rows after mapping position filter:", nrow(all_DT), "\n")
printAndLog(txt, pipLogFile)
stopifnot(all(rownames(all_DT) %in% gene2tadDT$entrezID))

### => now at the end I have an expression table that has
### - all rownames with a single gene name
### - all the gene names are not duplicated
### - they all have a mapping position available

rnaseqDT <- all_DT

refGenes <- rownames(rnaseqDT)
stopifnot(!any(grepl("///", refGenes)))

# there should be no duplicated from EZH2 data, and all should have mapping positions
stopifnot(!any(duplicated(refGenes)))
stopifnot(all(refGenes %in% gene2tadDT$entrezID))
cat(paste0("... Number of genes after EZH2 filtering: ", length(refGenes), "\n"))

### ---> NOT DONE FOR EZH2 DATA: filter genes that map to multiple TADs <<<<--- FOR EZH2 DATA, OK IF A GENE MAPS TO > 1 TAD
rna_geneList <- setNames(refGenes, refGenes)
cat("... for EZH2 data, genes allowed to map to multiple TADs")
# check the returned list
stopifnot(all(rna_geneList %in% gene2tadDT$entrezID))
stopifnot(all(names(rna_geneList) %in% refGenes))
stopifnot(!any(duplicated(rna_geneList)))
stopifnot(!any(duplicated(names(rna_geneList))))
stopifnot(names(rna_geneList) == rownames(rnaseqDT))


# rna_geneList will hold the list of genes for which I have 1) genomic positions 2) expression data


########################################################################################
# FILTER 1: DISCARD GENES WITH DUPLICATED GENE ID
# (only if removeDupGeneID == TRUE)
########################################################################################

# --> THIS WILL NOT WORK FOR EZH2 AS ROWNAMNES AND GENE LIST DO NOT CORRESPOND  -->>> this should be ok for ezh2v2
# in the new version, I systematically discard the ambiguous ones
cat("... FILTER1: genes with not duplicated gene ID\n")

if(removeDupGeneID) {
  cat("... remove ambiguous genes with duplicated entrezID\n")
  initLen <- length(rna_geneList)
  stopifnot(initLen == nrow(rnaseqDT))
  dupEntrez <- unique(as.character(rna_geneList[which(duplicated(rna_geneList))]))
  
  rm_genes_write <- rna_geneList[rna_geneList %in% dupEntrez]
  cat(paste0(rm_genes_write, collapse="\n"), file = paste0(curr_outFold, "/", "log_rm_genes_duplicated_id.txt"), append=FALSE)
  
  rna_geneList <- rna_geneList[!rna_geneList %in% dupEntrez]
  rnaseqDT <- rnaseqDT[rownames(rnaseqDT) %in% names(rna_geneList),]
  stopifnot(all(rownames(rnaseqDT) == names(rna_geneList)))
  stopifnot(length(rna_geneList) == nrow(rnaseqDT))
  txt <- paste0(toupper(script_name), "> Remove geneID with ambiguous entrezID, keep : ", length(rna_geneList), "/", initLen, "\n")
  printAndLog(txt, pipLogFile)
  stopifnot(!any(duplicated(rna_geneList)))
}


##########################################################################################################################
##########################################################################################################################
##########################################################################################################################

##### TAKE ONLY THE SAMPLES OF INTEREST
samp1 <- eval(parse(text=load(paste0(setDir, "/", sample1_file))))
samp2 <- eval(parse(text=load(paste0(setDir, "/", sample2_file))))
rnaseqDT <- rnaseqDT[, c(samp1, samp2)]
stopifnot( ncol(rnaseqDT) == (length(samp1) + length(samp2)) )

save(rna_geneList, file = paste0(curr_outFold, "/", "rna_geneList.Rdata"))
cat(paste0("... written: ", paste0(curr_outFold, "/", "rna_geneList.Rdata"), "\n"))

## ---> EZH2 => this is the DT that will be used for the DE analysis, should not contain duplicated rows
rna_rnaseqDT <- rnaseqDT
stopifnot(!any(duplicated(rna_rnaseqDT)))
save(rna_rnaseqDT, file = paste0(curr_outFold, "/", "rna_rnaseqDT.Rdata"))
cat(paste0("... written: ", paste0(curr_outFold, "/", "rna_rnaseqDT.Rdata"), "\n"))

rna_fpkmDT <- rna_rnaseqDT
stopifnot(!any(duplicated(rna_fpkmDT)))
save(rna_fpkmDT, file = paste0(curr_outFold, "/", "rna_fpkmDT.Rdata"))
cat(paste0("... written: ", paste0(curr_outFold, "/", "rna_fpkmDT.Rdata"), "\n"))

#### PREPARE THE QQNORM/MADNORM DATA FOR THE GENES I USED FOR DE ANALYSIS
cat("... normalize the data for other analyses (qqnorm for RNA-seq, MADnorm for microarray)\n")

stopifnot(!any(grepl("///", rownames(rna_rnaseqDT))))
full_rnaseqDT <- rna_rnaseqDT
stopifnot(!any(duplicated(full_rnaseqDT)))

if(microarray) {
  # --> CHANGED HERE FOR EZH2 <<------------- madnorm on the DT that can have multiple rows
  rna_rnaseqDT_tmp <- full_rnaseqDT
  # median center
  all_meds <- apply(rna_rnaseqDT_tmp, 1, median)
  rna_rnaseqDT_tmp <- sweep(rna_rnaseqDT_tmp, 1, all_meds, "-")  # substract from each row their corresponding median
  # mad norm
  # In order to use the MAD as a consistent estimator for the estimation of the standard deviation σ, one takes
  # σ ^ = k ⋅ MAD where k is a constant scale factor, which depends on the distribution.[1]
  # For normally distributed data k is taken to be: ~1.4826
  all_mads <-1.4826 * apply(abs(rna_rnaseqDT_tmp), 1, median)
  rna_madnorm_rnaseqDT <- sweep(rna_rnaseqDT_tmp, 1, all_mads, "/") 
  stopifnot( all ( dim(rna_rnaseqDT) == dim(rna_madnorm_rnaseqDT)))
  rownames(rna_madnorm_rnaseqDT) <- rownames(rna_rnaseqDT)
  colnames(rna_madnorm_rnaseqDT) <- colnames(rna_rnaseqDT)
  save(rna_madnorm_rnaseqDT, file = paste0(curr_outFold, "/", "rna_madnorm_rnaseqDT.Rdata"))
  cat(paste0("... written: ", paste0(curr_outFold, "/", "rna_madnorm_rnaseqDT.Rdata"), "\n"))
} else {
  stop("SHOULD NOT HAPPEN\n")
  rna_qqnorm_rnaseqDT <- t(apply(rna_rnaseqDT, 1, quantNorm))
  stopifnot(all(dim(rna_qqnorm_rnaseqDT) == dim(rna_rnaseqDT)))
  rownames(rna_qqnorm_rnaseqDT) <- rownames(rna_rnaseqDT)
  colnames(rna_qqnorm_rnaseqDT) <- colnames(rna_rnaseqDT)
  save(rna_qqnorm_rnaseqDT, file = paste0(curr_outFold, "/", "rna_qqnorm_rnaseqDT.Rdata"))
  cat(paste0("... written: ", paste0(curr_outFold, "/", "rna_qqnorm_rnaseqDT.Rdata"), "\n"))
}

# => rna_geneList: all the genes for which I have position information
# => countFilter_geneList: all the genes for which I have position information and that are >= CPM threshold

pipeline_geneList <- rna_geneList

########################################################################################
# FILTER 2: FILTER GENES THAT PASS MIN CPM THRESHOLD
# (only if useFilterCountData == TRUE)
########################################################################################
cat("... FILTER2: genes over min CPM threshold (or not)\n")

if(useFilterCountData) {
  stopifnot(is.numeric(rna_rnaseqDT[1,1]))
  # take only the genes for which I have position (filter 1)
  countFilter_rnaseqDT <- rna_rnaseqDT[which(rownames(rna_rnaseqDT) %in% names(pipeline_geneList)),]
  stopifnot(all(rownames(countFilter_rnaseqDT) == names(pipeline_geneList)))
  # ensure the samples are present in the column names
  stopifnot(all(samp1 %in% colnames(countFilter_rnaseqDT)))
  stopifnot(all(samp2 %in% colnames(countFilter_rnaseqDT)))
  countFilter_rnaseqDT <- countFilter_rnaseqDT[,c(samp1, samp2)]
  if(!microarray) {
    stop("SHOULD NOT HAPPEN\n")
    # FILTER THE EXPRESSION DATA TO MIN CPM FILTER
    cpm_exprDT <- cpm(countFilter_rnaseqDT)
    rowsToKeep <- rowSums(cpm_exprDT) >= (minCpmRatio * ncol(countFilter_rnaseqDT))
    txt <- paste0(toupper(script_name), "> useFilterCountData is TRUE -> CPM-filtered geneList; to keep: ", sum(rowsToKeep), "/", length(pipeline_geneList), "\n")
    printAndLog(txt, pipLogFile)
    pipeline_geneList <- pipeline_geneList[rowsToKeep]
    stopifnot(length(pipeline_geneList) == sum(rowsToKeep))
  } else {
    # CANNOT APPLY CPM FILTER  !
    txt <- paste0(toupper(script_name), "> !!! CPM filter not applied !!!", "\n")
    printAndLog(txt, pipLogFile)
    cpm_exprDT <- countFilter_rnaseqDT
    rowsToKeep <- rep(TRUE, nrow(cpm_exprDT))
    pipeline_geneList <- pipeline_geneList[rowsToKeep]
    stopifnot(length(pipeline_geneList) == sum(rowsToKeep))
    stopifnot(all(pipeline_geneList == pipeline_geneList))
  }
}

########################################################################################
# FILTER 3: FILTER TO USE ONLY TAD REGIONS (AND ONLY GENES BELONGING TO) NOT INTER-TAD
# (only if useTADonly == TRUE)
########################################################################################
cat("... FILTER3: use TAD regions only (or not) \n")

gene2tadDT_filtered <- read.delim(gene2tadDT_file, header=F, col.names = c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = F)
gene2tadDT_filtered$entrezID <- as.character(gene2tadDT_filtered$entrezID)
gene2tadDT_intact <- gene2tadDT_filtered

gene2tadDT_filtered <- gene2tadDT_filtered[gene2tadDT_filtered$entrezID %in% pipeline_geneList,]
### FOR EZH2 V2 THIS MIGHT NOT BE TRUE AS GENES CAN MAP TO MULTI TADS
# stopifnot(length(pipeline_geneList) == nrow(gene2tadDT_filtered))
initRegionsLen <- length(unique(gene2tadDT_filtered$region))

#setequal
stopifnot(all(gene2tadDT_filtered$entrezID %in% pipeline_geneList) & all(pipeline_geneList %in% gene2tadDT_filtered$entrezID))

if(useTADonly) {
  initNrow <- nrow(gene2tadDT_filtered)
  gene2tadDT_filtered <- gene2tadDT_filtered[grep("_TAD", gene2tadDT_filtered$region),]
  txt <- paste0(toupper(script_name), "> useTADonly is TRUE -> only TAD regions: ", nrow(gene2tadDT_filtered), "/", initNrow, " regions\n")
  printAndLog(txt, pipLogFile)
  initLen <- length(pipeline_geneList)
  
  cat(paste0(as.character(gene2tadDT_filtered$region), collapse="\n"), file = paste0(curr_outFold, "/", "log_TADs_retained_beforeSizeFilter.txt"), append=FALSE)
  
  rm_genes_write <- pipeline_geneList[!pipeline_geneList %in% gene2tadDT_filtered$entrezID]
  cat(paste0(rm_genes_write, collapse="\n"), file = paste0(curr_outFold, "/", "log_rm_genes_not_in_TAD.txt"), append=FALSE)
  
  pipeline_geneList <- pipeline_geneList[pipeline_geneList %in% gene2tadDT_filtered$entrezID]
  txt <- paste0(toupper(script_name), "> useTADonly is TRUE -> TAD-filtered geneList; to keep: ", length(pipeline_geneList), "/", initLen, "\n")
  printAndLog(txt, pipLogFile)
}

########################################################################################
# BEFORE FILTER 4: LOOK AT TAD GENE NUMBER DISTRIBUTION
# (only if useFilterSizeData == TRUE)
########################################################################################

cat("... FILTER4: filter TADs based on number of genes/TAD\n")

tadGeneNbr_DT <- aggregate(entrezID ~ region, data=gene2tadDT_filtered, FUN=length)
colnames(tadGeneNbr_DT) <- c("regionID", "nbrGenes") 

txt <-  paste0(toupper(script_name), "> Summary number of genes/TAD:", 
               paste0(paste0("... ", names(summary(tadGeneNbr_DT$nbrGenes)), ": ", round(as.numeric(summary(tadGeneNbr_DT$nbrGenes)),2)), collapse="\n"))

tadGeneNbr_DT$region <- ifelse(regexpr("_TAD", tadGeneNbr_DT$regionID) > 0, "TAD", "BOUND")
tadGeneNbr_DT$nbrLog10 <- log10(tadGeneNbr_DT$nbrGenes)

if(useTADonly) {
  stopifnot(all(tadGeneNbr_DT$region == "TAD"))
}

p1 <- ggdensity(tadGeneNbr_DT, x = "nbrLog10",
                title ="Distribution log10(# genes) before filtering",
                # add = "mean",
                rug = TRUE,
                xlab="log10 # of genes/region",
                color = "region", fill = "region",
                palette = c("#00AFBB", "#E7B800"))
p1 <- p1 + theme(plot.title = element_text(hjust = 0.5))

if(! useFilterSizeData) {
  maxQuantGeneTAD <- 1
  minNbrGeneTAD <- 0
}

upperLimit <- as.numeric(quantile(tadGeneNbr_DT$nbrGenes, probs=maxQuantGeneTAD))
txt <- paste0(toupper(script_name), "> Current threshold for TAD size: >= ", minNbrGeneTAD, " and <= ", upperLimit, " (", round(maxQuantGeneTAD*100,2), "%)\n")
printAndLog(txt, pipLogFile)
tadGeneNbr_DT <- tadGeneNbr_DT[tadGeneNbr_DT$nbrGenes >= minNbrGeneTAD & tadGeneNbr_DT$nbrGenes <= upperLimit,]

rangeTADgenes <- c(minNbrGeneTAD, upperLimit)
save(rangeTADgenes, file = paste0(curr_outFold, "/", "rangeTADgenes.Rdata"))
cat(paste0("... written: ", paste0(curr_outFold, "/", "rangeTADgenes.Rdata"), "\n"))

p2 <- ggdensity(tadGeneNbr_DT, x = "nbrLog10",
                title ="Distribution log10(# genes) after filtering",
                # add = "mean", 
                rug = TRUE,
                xlab="log10 # of genes/region",
                color = "region", fill = "region",
                palette = c("#00AFBB", "#E7B800"))
p2 <- p2 + theme(plot.title = element_text(hjust = 0.5))

p1_p2 <- ggarrange(p1, p2, 
                   labels = c("A", "B"),
                   ncol = 2, nrow = 1)

ggsave(filename = paste0(curr_outFold, "/", "TADsize_density.", plotType), plot = p1_p2, height=7, width=10)


########################################################################################
# FILTER 4: FILTER TAD SIZE TO USE ONLY TAD (AND ONLY GENES BELONGING TO) >= LOWER AND <= UPPER BOUND
# (only if useFilterSizeData == TRUE)
########################################################################################

#setequal
stopifnot(all(gene2tadDT_filtered$entrezID %in% pipeline_geneList) & all(pipeline_geneList %in% gene2tadDT_filtered$entrezID))

if(useFilterSizeData){
  
  rm_genes_write <- unique(as.character(gene2tadDT_filtered$region[gene2tadDT_filtered$region %in% as.character(tadGeneNbr_DT$regionID)]))
  cat(paste0(rm_genes_write, collapse="\n"), file = paste0(curr_outFold, "/", "log_rm_regions_tadsize.txt"), append=FALSE)
  
  gene2tadDT_filtered <- gene2tadDT_filtered[gene2tadDT_filtered$region %in% as.character(tadGeneNbr_DT$regionID),]
  initLen <- length(pipeline_geneList)
  
  cat(paste0(as.character(gene2tadDT_filtered$region), collapse="\n"), file = paste0(curr_outFold, "/", "log_TADs_retained_afterSizeFilter.txt"), append=FALSE)
  
  rm_genes_write <- pipeline_geneList[!pipeline_geneList %in% gene2tadDT_filtered$entrezID]
  cat(paste0(rm_genes_write, collapse="\n"), file = paste0(curr_outFold, "/", "log_rm_genes_tadsize.txt"), append=FALSE)
  
  pipeline_geneList <- pipeline_geneList[pipeline_geneList %in% gene2tadDT_filtered$entrezID]
  txt <- paste0(toupper(script_name), "> useFilterSizeData is TRUE -> size-filtered geneList; to keep: ", length(pipeline_geneList), "/", initLen, "\n")
  printAndLog(txt, pipLogFile)
}

txt <- paste0(toupper(script_name), "> pipeline_geneList compared to available mapped genes: ", length(pipeline_geneList), "/", length(rna_geneList), "\n")
printAndLog(txt, pipLogFile)
save(pipeline_geneList, file = paste0(curr_outFold, "/", "pipeline_geneList.Rdata"))
cat(paste0("... written: ", paste0(curr_outFold, "/", "pipeline_geneList.Rdata"), "\n"))

pipeline_regionList <- unique(as.character(gene2tadDT_filtered$region))
txt <- paste0(toupper(script_name), "> pipeline_regionList compared to available regions: ", length(pipeline_regionList), "/", initRegionsLen, "\n")
printAndLog(txt, pipLogFile)
save(pipeline_regionList, file = paste0(curr_outFold, "/", "pipeline_regionList.Rdata"))
cat(paste0("... written: ", paste0(curr_outFold, "/", "pipeline_regionList.Rdata"), "\n"))

########################################################################################################################## 
########################################################### ADDITIONAL STEPS HERE FOR EZH2 - LOOK AT SESN1 AND SESN2
##########################################################################################################################
outFile <-  paste0(curr_outFold, "/", "SESN_TAD_genes_positionList.txt")
system(paste0("rm -f ", outFile))
sink(outFile)

mygenes <- c("FOXO3", "ARMC2", "SESN1", "SESN2")

gene2tadDT <- read.delim(gene2tadDT_file, header=F, col.names = c("entrezID", "chromo", "gene_start", "gene_end", "region"), stringsAsFactors = F)
gene2tadDT <- gene2tadDT[gene2tadDT$entrezID %in% mygenes,]

TADposDT <- read.delim(TADpos_file, header=F, col.names = c("chromo", "region", "TAD_start", "TAD_end"), stringsAsFactors = F)
mytadDT <- left_join(gene2tadDT, TADposDT, by=c("chromo", "region"))
mytadDT <- mytadDT[order(mytadDT$chromo, mytadDT$gene_start),]
print(mytadDT)

for(gene in mygenes) {
  cat(paste0(gene, " retained in pipeline_geneList: ", as.character(gene %in% pipeline_geneList), "\n"))
}
sink()

cat(paste0("... written: ", outFile, "\n"))

########################################################################################################################## 
txt <- paste0(startTime, "\n", Sys.time(), "\n")
printAndLog(txt, pipLogFile)
cat(paste0("*** DONE: ", script_name, "\n"))

