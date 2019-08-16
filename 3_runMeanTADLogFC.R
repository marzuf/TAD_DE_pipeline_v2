#!/usr/bin/Rscript

startTime <- Sys.time()

################  USE THE FOLLOWING FILES FROM PREVIOUS STEPS
# - script0: pipeline_regionList.Rdata
# - script0: rna_geneList.Rdata
# - script0: pipeline_geneList.Rdata
# - script0: rna_madnorm_rnaseqDT.Rdata
# - script1: DE_topTable.Rdata
# - script1: DE_geneList.Rdata
################################################################################

################  OUTPUT
# - /all_meanLogFC_TAD.Rdata
################################################################################

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1)
settingF <- args[1]
stopifnot(file.exists(settingF))

pipScriptDir <- paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2")

script0_name <- "0_prepGeneData"
script1_name <- "1_runGeneDE"
script_name <- "3_runMeanTADLogFC"
stopifnot(file.exists(paste0(pipScriptDir, "/", script_name, ".R")))
cat(paste0("> START ", script_name,  "\n"))

source("main_settings.R")
source(settingF)
source(paste0(pipScriptDir, "/", "TAD_DE_utils.R"))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

# create the directories
curr_outFold <- paste0(pipOutFold, "/", script_name)
system(paste0("mkdir -p ", curr_outFold))

pipLogFile <- paste0(pipOutFold, "/", format(Sys.time(), "%Y%d%m%H%M%S"),"_", script_name, "_logFile.txt")
system(paste0("rm -f ", pipLogFile))

registerDoMC(ifelse(SSHFS, 2, nCpu)) # from main_settings.R

# ADDED 16.11.2018 to check using other files
txt <- paste0("inputDataType\t=\t", inputDataType, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0("gene2tadDT_file\t=\t", gene2tadDT_file, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0("TADpos_file\t=\t", TADpos_file, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0("settingF\t=\t", settingF, "\n")
printAndLog(txt, pipLogFile)

################################***********************************************************************************
############ LOAD INPUT DATA
################################***********************************************************************************
gene2tadDT <- read.delim(gene2tadDT_file, header=F, col.names = c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = F)
gene2tadDT$entrezID <- as.character(gene2tadDT$entrezID)

DE_topTable <- eval(parse(text = load(paste0(pipOutFold, "/", script1_name, "/DE_topTable.Rdata"))))
DE_geneList <- eval(parse(text = load(paste0(pipOutFold, "/", script1_name, "/DE_geneList.Rdata"))))

pipeline_geneList <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/pipeline_geneList.Rdata"))))
pipeline_regionList <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/pipeline_regionList.Rdata"))))

if(useTADonly) {
  if(any(grepl("_BOUND", pipeline_regionList))) {
    stop("! data were not prepared for \"useTADonly\" !")
  }
}

stopifnot(all(DE_topTable$genes %in% names(DE_geneList)))
stopifnot(!any(duplicated(names(DE_geneList))))

entrezList <- unlist(sapply(DE_topTable$genes, function(x) DE_geneList[x]))
names(entrezList) <- DE_topTable$genes
stopifnot(length(entrezList) == length(DE_topTable$genes))

# replace the gene symbol rownames by ensemblID rownames
logFC_DT <- data.frame(entrezID =  entrezList,
                       logFC = DE_topTable$logFC, stringsAsFactors = F)

rownames(logFC_DT) <- NULL
initNrow <- nrow(logFC_DT)
logFC_DT <- logFC_DT[logFC_DT$entrezID %in% pipeline_geneList,]
txt <- paste0(toupper(script_name), "> Take only filtered genes: ", nrow(logFC_DT), "/", initNrow, "\n")
printAndLog(txt, pipLogFile)

### take only the filtered data according to initial settings
gene2tadDT <- gene2tadDT[gene2tadDT$entrezID %in% as.character(pipeline_geneList),]
initLen <- length(unique(gene2tadDT$region))
gene2tadDT <- gene2tadDT[gene2tadDT$region %in% pipeline_regionList,]
txt <- paste0(toupper(script_name), "> Take only filtered regions: ", length(unique(gene2tadDT$region)), "/", initLen, "\n")
printAndLog(txt, pipLogFile)

################################***********************************************************************************
################################********************************************* get observed logFC for all regions
################################***********************************************************************************

cat(paste0("... start computing mean logFC by TAD \n"))

head(logFC_DT)

mergedDT <- left_join(logFC_DT, gene2tadDT[,c("entrezID", "region")], by="entrezID")


save(mergedDT, file="mergedDT.Rdata")
save(logFC_DT, file="logFC_DT.Rdata")
save(gene2tadDT, file="gene2tadDT.Rdata")

stopifnot(nrow(mergedDT) == nrow(na.omit(mergedDT)))

mean_DT <- aggregate(logFC ~ region, data=mergedDT, FUN=mean)
all_meanLogFC_TAD <- setNames(mean_DT$logFC, mean_DT$region)
stopifnot(length(all_meanLogFC_TAD) == length(unique(gene2tadDT$region)))
txt <- paste0(toupper(script_name), "> Number of regions for which mean logFC computed: ", length(all_meanLogFC_TAD), "\n")
printAndLog(txt, pipLogFile)

if(useTADonly) {
    initLen <- length(all_meanLogFC_TAD)
    all_meanLogFC_TAD <- all_meanLogFC_TAD[grep("_TAD", names(all_meanLogFC_TAD))]
    txt <- paste0(toupper(script_name), "> Take only the TAD regions: ", length(all_meanLogFC_TAD),"/", initLen, "\n")
    printAndLog(txt, pipLogFile)
}

save(all_meanLogFC_TAD, file= paste0(curr_outFold, "/all_meanLogFC_TAD.Rdata"))
cat(paste0("... written: ", curr_outFold, "/all_meanLogFC_TAD.Rdata", "\n"))

txt <- paste0(startTime, "\n", Sys.time(), "\n")
printAndLog(txt, pipLogFile)

cat(paste0("*** DONE: ", script_name, "\n"))

