#!/usr/bin/Rscript

options(scipen=100)

startTime <- Sys.time()

### !!! HARD CODED corMet "pearson"

################  USE THE FOLLOWING FILES FROM PREVIOUS STEPS
# - script0: pipeline_regionList.Rdata
# - script0: pipeline_geneList.Rdata
# - script0: rna_madnorm_rnaseqDT.Rdata or rna_qqnorm_rnaseqDT.Rdata
# - script5sameNbr: sample_around_TADs_sameNbr.Rdata
################################################################################

################  OUTPUT
# - meanCorr_permDT.Rdata
################################################################################

### !!! TAKE ALL DATA FROM THE DS IN THE FOLDER !!!

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1)
settingF <- args[1]
stopifnot(file.exists(settingF))

pipScriptDir <- paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2")

script0_name <- "0_prepGeneData"
script1_name <- "1_runGeneDE"
script5sameNbr_name <- "5sameNbr_runPermutationsCorr"
script_name <- "6sameNbr_runPermutationsMeanLogFC"
stopifnot(file.exists(paste0(pipScriptDir, "/", script_name, ".R")))
cat(paste0("> START ", script_name,  "\n"))

source("main_settings.R")
source(settingF)
source(file.path(pipScriptDir, "TAD_DE_utils.R"))



suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
registerDoMC(ifelse(SSHFS, 2, nCpu)) # from main_settings.R
suppressPackageStartupMessages(library(reshape2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

# if microarray was not set in the settings file -> by default set to  FALSE
if(!exists("microarray")) microarray <- FALSE

# create the directories
curr_outFold <- paste0(pipOutFold, "/", script_name)
system(paste0("mkdir -p ", curr_outFold))

pipLogFile <- paste0(pipOutFold, "/", format(Sys.time(), "%Y%d%m%H%M%S"),"_", script_name, "_logFile.txt")
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


geneList_file <- file.path(pipOutFold, script0_name, "pipeline_geneList.Rdata")
stopifnot(file.exists(geneList_file))
geneList <- eval(parse(text = load(geneList_file)))

regList_file <- file.path(pipOutFold, script0_name, "pipeline_regionList.Rdata")
stopifnot(file.exists(regList_file))
regionList <- eval(parse(text = load(regList_file)))


DE_topTable <- eval(parse(text = load(paste0(pipOutFold, "/", script1_name, "/DE_topTable.Rdata"))))
DE_topTable$genes <- as.character(DE_topTable$genes)

ds_sample_data_file <- file.path(pipOutFold, script5sameNbr_name, "sample_around_TADs_sameNbr.Rdata")
stopifnot(file.exists(ds_sample_data_file))
ds_sample_data <- eval(parse(text = load(ds_sample_data_file)))

all_regs <- names(ds_sample_data)
stopifnot(setequal(all_regs, regionList))

meanFC_sample_around_TADs_sameNbr <- foreach(reg = all_regs) %dopar% {

  tad_data <- ds_sample_data[[paste0(reg)]]
  
#save(tad_data, file="tad_data.Rdata", version=2)

  if(tad_data$nGenes > 0) {
    ########## => TAKING THE TADs AND SAMPLING ON BOTH SIDES OF THE TADs
    sample_genes <- tad_data$genes

	stopifnot(sample_genes %in% c(tad_data$genes_left, tad_data$genes_right))

#cat(paste0("sample_genes", "\n"))
#cat(paste0(sample_genes ,"\n"))
#cat(paste0("right_genes", "\n"))
#cat(paste0(tad_data$genes_left ,"\n"))
#cat(paste0("left_genes", "\n"))
#cat(paste0(tad_data$genes_right,"\n"))

    tad_genes <- tad_data$tad_genes
    stopifnot(! sample_genes %in% tad_genes)
    stopifnot(sample_genes %in% geneList)
    
    sampleTAD_genes <- names(geneList)[geneList %in% sample_genes]  # the names used in the nrom_rnaseqDT
    
	stopifnot(sampleTAD_genes %in% DE_topTable$genes)

    meanFC <- mean(DE_topTable$logFC[DE_topTable$genes %in% sampleTAD_genes])

	stopifnot(!is.na(meanFC))


  } else  {

	meanFC <- NA

  }

  if(tad_data$nGenes_right > 0) {
    ########## => TAKING THE TADs AND SAMPLING ON BOTH SIDES OF THE TADs
    sample_genes_right <- tad_data$genes_right


    tad_genes <- tad_data$tad_genes
    stopifnot(! sample_genes_right %in% tad_genes)
    stopifnot(sample_genes_right %in% geneList)
    
    sampleTAD_genes_right <- names(geneList)[geneList %in% sample_genes_right]  # the names used in the nrom_rnaseqDT
    
	stopifnot(sampleTAD_genes_right %in% DE_topTable$genes)

    meanFC_right <- mean(DE_topTable$logFC[DE_topTable$genes %in% sampleTAD_genes_right])

	stopifnot(!is.na(meanFC_right))


  } else  {

	meanFC_right<- NA

  }


  if(tad_data$nGenes_left > 0) {
    ########## => TAKING THE TADs AND SAMPLING ON BOTH SIDES OF THE TADs
    sample_genes_left <- tad_data$genes_left


    tad_genes <- tad_data$tad_genes
    stopifnot(! sample_genes_left %in% tad_genes)
    stopifnot(sample_genes_left %in% geneList)
    
    sampleTAD_genes_left <- names(geneList)[geneList %in% sample_genes_left]  # the names used in the nrom_rnaseqDT
    
	stopifnot(sampleTAD_genes_left %in% DE_topTable$genes)

    meanFC_left <- mean(DE_topTable$logFC[DE_topTable$genes %in% sampleTAD_genes_left])

	stopifnot(!is.na(meanFC_left))


  } else  {

	meanFC_left<- NA

  }



  list(meanLogFC = meanFC,meanLogFC_left=meanFC_left, meanLogFC_right=meanFC_right )
} # end iterating over all regions
names(meanFC_sample_around_TADs_sameNbr) <- all_regs 

outFile <- file.path(curr_outFold, "meanFC_sample_around_TADs_sameNbr.Rdata")
save(meanFC_sample_around_TADs_sameNbr, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

########################################################
########################################################
########################################################
txt <- paste0(startTime, "\n", Sys.time(), "\n")
printAndLog(txt, pipLogFile)
cat(paste0("*** DONE: ", script_name, "\n"))






