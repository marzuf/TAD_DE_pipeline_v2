startTime <- Sys.time()

################  OUTPUT
# - rna_geneList.Rdata
# - rna_rnaseqDT.Rdata
# - rna_madnorm_rnaseqDT.Rdata
# - pipeline_geneList.Rdata
# - pipeline_regionList.Rdata
# (+ some pictures)
################################################################################


args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1)
settingF <- args[1]
stopifnot(file.exists(settingF))

script_name <- "A_prepGeneData"
stopifnot(file.exists(paste0(script_name, ".R")))
cat(paste0("> START ", script_name,  "\n"))
# -> rename to 0_ so that saved data can be used in following scripts
script_name <- "0_prepGeneData"

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")
# cat(paste0("setDir = ", setDir, "\n"))
source("main_settings.R") # setDir is the main_settings not in run_settings
#### IMPORTANT !! run_settings.R SHOULD BE AFTER main_settings.R TO OVERWRITE DEFAULT FILES !!!
#source("run_settings.R")
source(settingF)
source("TAD_DE_utils.R")
suppressPackageStartupMessages(library(edgeR, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(ggplot2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(ggpubr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

# if microarray was not set in the settings file -> by default set to  FALSE
if(!exists("microarray")) microarray <- FALSE

cat(paste0("... RNA DATA COME FROM MICROARRAY: ", as.character(microarray), "\n"))

plotType <- "svg"

##### IN THE VERSION A_ INSTEAD OF 0_
# take as input a specific gene2tad file
# 

######### 0_prepGene
# ->  according to the user settings, create a geneList object 
# ... with the same genes as rownames in input data
# ... geneList["geneName"] = entrezID
#####################

# create the directories
curr_outFold <- paste0(pipOutFold, "/", script_name)
system(paste0("mkdir -p ", curr_outFold))

pipLogFile <- paste0(pipOutFold, "/", format(Sys.time(), "%Y%d%m%H%M%S"),"_", script_name, "_logFile.txt")
system(paste0("rm -f ", pipLogFile))

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

refGenes <- as.character(rownames(rnaseqDT))
cat(paste0("... initial number of rows (genes): ", length(refGenes), "\n"))


gene2tadDT <- read.delim(gene2tadDT_file, header=F, col.names = c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = F)
gene2tadDT$entrezID <- as.character(gene2tadDT$entrezID)

# remove the genes that cannot be mapped unambigously
dupG2Tentrez <- gene2tadDT$entrezID[duplicated(gene2tadDT$entrezID)]
gene2tadDT <- gene2tadDT[!gene2tadDT$entrezID %in% dupG2Tentrez,]
stopifnot(!any(duplicated(gene2tadDT$entrezID)))

# the list of genes to retain as input
initNrow <- nrow(rnaseqDT)
cat("... remove ambiguous genes mapped to more than one TAD\n")
rnaseqDT <- rnaseqDT[rownames(rnaseqDT) %in% gene2tadDT$entrezID,]
txt <- paste0(toupper(script_name), "> Remove geneID mapping to more than one TAD: ", nrow(rnaseqDT), "/", initNrow, "\n")
printAndLog(txt, pipLogFile)
rna_geneList <- setNames(rownames(rnaseqDT), rownames(rnaseqDT))

# check the returned list
stopifnot(all(rna_geneList %in% gene2tadDT$entrezID))
stopifnot(all(names(rna_geneList) %in% refGenes))
stopifnot(!any(duplicated(rna_geneList)))
stopifnot(!any(duplicated(names(rna_geneList))))

txt <- paste0(toupper(script_name), "> Number of genes retained: ", length(rna_geneList), "/", length(refGenes), "\n")
printAndLog(txt, pipLogFile)
# discard from the rnaseqDT the rows for which I have no information
rnaseqDT <- rnaseqDT[rownames(rnaseqDT) %in% names(rna_geneList),]
stopifnot(all(rownames(rnaseqDT) == names(rna_geneList)))

# in the new version, I systematically discard the ambiguous ones
if(removeDupGeneID) {
  cat("... remove ambiguous genes with duplicated entrezID\n")
  initLen <- length(rna_geneList)
  stopifnot(initLen == nrow(rnaseqDT))
  dupEntrez <- unique(as.character(rna_geneList[which(duplicated(rna_geneList))]))
  rna_geneList <- rna_geneList[!rna_geneList %in% dupEntrez]
  rnaseqDT <- rnaseqDT[rownames(rnaseqDT) %in% names(rna_geneList),]
  stopifnot(all(rownames(rnaseqDT) == names(rna_geneList)))
  stopifnot(length(rna_geneList) == nrow(rnaseqDT))
  txt <- paste0(toupper(script_name), "> Remove geneID with ambiguous entrezID, final # of rows  : ", length(rna_geneList), "/", initLen, "\n")
  printAndLog(txt, pipLogFile)
  stopifnot(!any(duplicated(rna_geneList)))
}


##### TAKE ONLY THE SAMPLES OF INTEREST
samp1 <- eval(parse(text=load(paste0(setDir, "/", sample1_file))))
samp2 <- eval(parse(text=load(paste0(setDir, "/", sample2_file))))
rnaseqDT <- rnaseqDT[, c(samp1, samp2)]
stopifnot( ncol(rnaseqDT) == (length(samp1) + length(samp2)) )

save(rna_geneList, file = paste0(curr_outFold, "/", "rna_geneList.Rdata"))
cat(paste0("... written: ", paste0(curr_outFold, "/", "rna_geneList.Rdata"), "\n"))

rna_rnaseqDT <- rnaseqDT
save(rna_rnaseqDT, file = paste0(curr_outFold, "/", "rna_rnaseqDT.Rdata"))
cat(paste0("... written: ", paste0(curr_outFold, "/", "rna_rnaseqDT.Rdata"), "\n"))

#### PREPARE THE QQNORM/MADNORM DATA FOR THE GENES I USED FOR DE ANALYSIS
cat("... normalize the data for other analyses (qqnorm for RNA-seq, MADnorm for microarray)\n")

if(microarray) {
  rna_rnaseqDT_tmp <- rna_rnaseqDT
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


# 
# => rna_geneList: all the genes for which I have position information
# => countFilter_geneList: all the genes for which I have position information and that are >= CPM threshold

pipeline_geneList <- rna_geneList

########################################################################################
# FILTER 2: FILTER GENES THAT PASS MIN CPM THRESHOLD
# (only if useFilterCountData == TRUE)
########################################################################################
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
  }
}

########################################################################################
# FILTER 3: FILTER TO USE ONLY TAD REGIONS (AND ONLY GENES BELONGING TO) NOT INTER-TAD
# (only if useTADonly == TRUE)
########################################################################################
gene2tadDT <- read.delim(gene2tadDT_file, header=F, col.names = c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = F)
gene2tadDT$entrezID <- as.character(gene2tadDT$entrezID)
gene2tadDT_intact <- gene2tadDT

gene2tadDT <- gene2tadDT[gene2tadDT$entrezID %in% pipeline_geneList,]
stopifnot(length(pipeline_geneList) == nrow(gene2tadDT))
initRegionsLen <- length(unique(gene2tadDT$region))

#setequal
stopifnot(all(gene2tadDT$entrezID %in% pipeline_geneList) & all(pipeline_geneList %in% gene2tadDT$entrezID))

if(useTADonly) {
  initNrow <- nrow(gene2tadDT)
  gene2tadDT <- gene2tadDT[grep("_TAD", gene2tadDT$region),]
  txt <- paste0(toupper(script_name), "> useTADonly is TRUE -> only TAD regions: ", nrow(gene2tadDT), "/", initNrow, " regions\n")
  printAndLog(txt, pipLogFile)
  initLen <- length(pipeline_geneList)
  pipeline_geneList <- pipeline_geneList[pipeline_geneList %in% gene2tadDT$entrezID]
  txt <- paste0(toupper(script_name), "> useTADonly is TRUE -> TAD-filtered geneList; to keep: ", length(pipeline_geneList), "/", initLen, "\n")
  printAndLog(txt, pipLogFile)
}

########################################################################################
# BEFORE FILTER 4: LOOK AT TAD GENE NUMBER DISTRIBUTION
# (only if useFilterSizeData == TRUE)
########################################################################################

tadGeneNbr_DT <- aggregate(entrezID ~ region, data=gene2tadDT, FUN=length)
colnames(tadGeneNbr_DT) <- c("regionID", "nbrGenes") 
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


upperLimit <- as.numeric(quantile(tadGeneNbr_DT$nbrGenes, probs=maxQuantGeneTAD))
txt <- paste0(toupper(script_name), "> Current threhsold for TAD size: >= ", minNbrGeneTAD, " and <= ", upperLimit, " (", maxQuantGeneTAD, "%)\n")
printAndLog(txt, pipLogFile)
tadGeneNbr_DT <- tadGeneNbr_DT[tadGeneNbr_DT$nbrGenes >= minNbrGeneTAD & tadGeneNbr_DT$nbrGenes <= upperLimit,]

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
stopifnot(all(gene2tadDT$entrezID %in% pipeline_geneList) & all(pipeline_geneList %in% gene2tadDT$entrezID))

if(useFilterSizeData){
  gene2tadDT <- gene2tadDT[gene2tadDT$region %in% as.character(tadGeneNbr_DT$regionID),]
  initLen <- length(pipeline_geneList)
  pipeline_geneList <- pipeline_geneList[pipeline_geneList %in% gene2tadDT$entrezID]
  txt <- paste0(toupper(script_name), "> useFilterSizeData is TRUE -> size-filtered geneList; to keep: ", length(pipeline_geneList), "/", initLen, "\n")
  printAndLog(txt, pipLogFile)
}


txt <- paste0(toupper(script_name), "> pipeline_geneList compared to available genes: ", length(pipeline_geneList), "/", length(rna_geneList), "\n")
printAndLog(txt, pipLogFile)
save(pipeline_geneList, file = paste0(curr_outFold, "/", "pipeline_geneList.Rdata"))
cat(paste0("... written: ", paste0(curr_outFold, "/", "pipeline_geneList.Rdata"), "\n"))

pipeline_regionList <- unique(as.character(gene2tadDT$region))
txt <- paste0(toupper(script_name), "> pipeline_regionList compared to available regions: ", length(pipeline_regionList), "/", initRegionsLen, "\n")
printAndLog(txt, pipLogFile)
save(pipeline_regionList, file = paste0(curr_outFold, "/", "pipeline_regionList.Rdata"))
cat(paste0("... written: ", paste0(curr_outFold, "/", "pipeline_regionList.Rdata"), "\n"))

txt <- paste0(startTime, "\n", Sys.time(), "\n")
printAndLog(txt, pipLogFile)

cat(paste0("*** DONE: ", script_name, "\n"))

