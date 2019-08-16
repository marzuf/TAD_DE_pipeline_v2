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
script_name <- "7sameNbr_runPermutationsMeanTADCorr"
stopifnot(file.exists(paste0(pipScriptDir, "/", script_name, ".R")))
cat(paste0("> START ", script_name,  "\n"))

source("main_settings.R")
source(settingF)
source(file.path(pipScriptDir, "TAD_DE_utils.R"))
source(file.path(pipScriptDir,  "TAD_DE_utils_meanCorr.R"))

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

corrMeth <- "pearson"
txt <- paste0("corrMeth\t=\t", corrMeth, "\n")
printAndLog(txt, pipLogFile)





geneList_file <- file.path(pipOutFold, script0_name, "pipeline_geneList.Rdata")
stopifnot(file.exists(geneList_file))
geneList <- eval(parse(text = load(geneList_file)))



regionList_file <- file.path(pipOutFold, script0_name, "pipeline_regionList.Rdata")
stopifnot(file.exists(regionList_file))
regionList <- eval(parse(text = load(regionList_file)))




sample1_file <- file.path(setDir, sample1_file)
sample2_file <- file.path(setDir, sample2_file)
stopifnot(file.exists(sample1_file))
stopifnot(file.exists(sample2_file))

cond1_ID <- eval(parse(text = load(sample1_file)))
cond2_ID <- eval(parse(text = load(sample2_file)))

qqnormDTfile <- file.path(pipOutFold,script0_name, "rna_qqnorm_rnaseqDT.Rdata")
stopifnot(file.exists(qqnormDTfile))
qqnormDT <- eval(parse(text = load(qqnormDTfile)))

stopifnot(names(geneList) %in% rownames(qqnormDT))
stopifnot(setequal(colnames(qqnormDT), c(cond1_ID, cond2_ID)))

norm_rnaseqDT <- qqnormDT[names(geneList),]    # !!! ENSURE THAT THE QQNORM IN THE SAME ORDER AS THE GENELIST !!!

stopifnot(rownames(norm_rnaseqDT) == names(geneList))
stopifnot(!duplicated(names(geneList)))

ds_sample_data_file <- file.path(pipOutFold, script5sameNbr_name, "sample_around_TADs_sameNbr.Rdata")
stopifnot(file.exists(ds_sample_data_file))
ds_sample_data <- eval(parse(text = load(ds_sample_data_file)))

all_regs <- names(ds_sample_data)
stopifnot(setequal(all_regs, regionList))

meanCorr_sample_around_TADs_sameNbr <- foreach(reg = all_regs) %dopar% {

  tad_data <- ds_sample_data[[paste0(reg)]]
  
  if(tad_data$nGenes > 0) {
    ########## => TAKING THE TADs AND SAMPLING ON BOTH SIDES OF THE TADs
    sample_genes <- tad_data$genes
    tad_genes <- tad_data$tad_genes
    stopifnot(! sample_genes %in% tad_genes)
    stopifnot(sample_genes %in% geneList)
    
    inTAD_genes <- names(geneList)[geneList %in% tad_genes]      # the names used in the nrom_rnaseqDT
    outTAD_genes <- names(geneList)[geneList %in% sample_genes]  # the names used in the nrom_rnaseqDT
    
    nTotGenes <- length(inTAD_genes) + length(outTAD_genes)
    
    stopifnot(! inTAD_genes %in% outTAD_genes)
    stopifnot(! outTAD_genes %in% inTAD_genes)
    stopifnot(inTAD_genes %in% rownames(norm_rnaseqDT))
    stopifnot(outTAD_genes %in% rownames(norm_rnaseqDT))
    
    sub_normDT <- norm_rnaseqDT[c(inTAD_genes, outTAD_genes),]
    
    stopifnot(nrow(sub_normDT) == nTotGenes)
    stopifnot(rownames(sub_normDT) == c(inTAD_genes, outTAD_genes))
    stopifnot(cond1_ID %in% colnames(sub_normDT))
    stopifnot(cond2_ID %in% colnames(sub_normDT))
    
    sub_normDT_cond1 <- sub_normDT[,cond1_ID]
    sub_normDT_cond2 <- sub_normDT[,cond2_ID]
    
    stopifnot(nrow(sub_normDT) == nrow(sub_normDT_cond1))
    stopifnot(nrow(sub_normDT) == nrow(sub_normDT_cond2))
    stopifnot( ncol(sub_normDT_cond1) + ncol(sub_normDT_cond2) == ncol(sub_normDT))
    stopifnot( ncol(sub_normDT_cond1) == length(cond1_ID))
    stopifnot(ncol(sub_normDT_cond2) == length(cond2_ID))
    
    meanCorr_all <- get_meanCorr_value(
                   exprMatrix = sub_normDT, 
                   inside_genes = inTAD_genes,      # inside_genes and outside_genes should be in rownames of exprMatrix
                   outside_genes = outTAD_genes, 
                   cormet = corrMeth
                   )
    meanCorr_cond1 <- get_meanCorr_value(
      exprMatrix = sub_normDT_cond1, 
      inside_genes = inTAD_genes, 
      outside_genes = outTAD_genes, 
      cormet = corrMeth
    )
    meanCorr_cond2 <- get_meanCorr_value(
      exprMatrix = sub_normDT_cond2, 
      inside_genes = inTAD_genes, 
      outside_genes = outTAD_genes, 
      cormet = corrMeth
    )
    ########## => TAKING THE TADs AND SAMPLING ON THE RIGHT ONLY
    if(tad_data$nGenes_right > 0) {
      sample_genes_right <- tad_data$genes_right
      outTAD_genes_right <- names(geneList)[geneList %in% sample_genes_right]  # the names used in the nrom_rnaseqDT
      stopifnot(! sample_genes_right %in% tad_genes)
      nTotGenes_right <- length(inTAD_genes) + length(outTAD_genes_right)
      
      stopifnot(! inTAD_genes %in% outTAD_genes_right)
      stopifnot(! outTAD_genes_right %in% inTAD_genes)
      stopifnot(outTAD_genes_right %in% rownames(norm_rnaseqDT))
      
      sub_normDT_right <- norm_rnaseqDT[c(inTAD_genes, outTAD_genes_right),]
      
      stopifnot(nrow(sub_normDT_right) == nTotGenes_right)
      stopifnot(rownames(sub_normDT_right) == c(inTAD_genes, outTAD_genes_right))
      stopifnot(cond1_ID %in% colnames(sub_normDT_right))
      stopifnot(cond2_ID %in% colnames(sub_normDT_right))
      
      sub_normDT_cond1_right <- sub_normDT_right[,cond1_ID]
      sub_normDT_cond2_right <- sub_normDT_right[,cond2_ID]
      
      stopifnot(nrow(sub_normDT_right) == nrow(sub_normDT_cond1_right))
      stopifnot(nrow(sub_normDT_right) == nrow(sub_normDT_cond2_right))
      stopifnot( ncol(sub_normDT_cond1_right) + ncol(sub_normDT_cond2_right) == ncol(sub_normDT_right))
      stopifnot( ncol(sub_normDT_cond1_right) == length(cond1_ID))
      stopifnot(ncol(sub_normDT_cond2_right) == length(cond2_ID))
      
      meanCorrRight_all <- get_meanCorr_value(
        exprMatrix = sub_normDT_right, 
        inside_genes = inTAD_genes,      # inside_genes and outside_genes should be in rownames of exprMatrix
        outside_genes = outTAD_genes_right, 
        cormet = corrMeth
      )
      meanCorrRight_cond1 <- get_meanCorr_value(
        exprMatrix = sub_normDT_cond1_right, 
        inside_genes = inTAD_genes,      # inside_genes and outside_genes should be in rownames of exprMatrix
        outside_genes = outTAD_genes_right, 
        cormet = corrMeth
      )
      meanCorrRight_cond2 <- get_meanCorr_value(
        exprMatrix = sub_normDT_cond2_right, 
        inside_genes = inTAD_genes,      # inside_genes and outside_genes should be in rownames of exprMatrix
        outside_genes = outTAD_genes_right, 
        cormet = corrMeth
      )
    } else {
      meanCorrRight_all <- NA
      meanCorrRight_cond1 <- NA
      meanCorrRight_cond2 <- NA
    }
    if(tad_data$nGenes_left > 0) {
      sample_genes_left <- tad_data$genes_left
      stopifnot(! sample_genes_left %in% tad_genes)
      outTAD_genes_left <- names(geneList)[geneList %in% sample_genes_left]  # the names used in the nrom_rnaseqDT
      stopifnot(! sample_genes_left %in% tad_genes)
      
      nTotGenes_left <- length(inTAD_genes) + length(outTAD_genes_left)
      stopifnot(! inTAD_genes %in% outTAD_genes_left)
      stopifnot(! outTAD_genes_left %in% inTAD_genes)
      stopifnot(outTAD_genes_left %in% rownames(norm_rnaseqDT))
      
      sub_normDT_left <- norm_rnaseqDT[c(inTAD_genes, outTAD_genes_left),]
      stopifnot(nrow(sub_normDT_left) == nTotGenes_left)
      stopifnot(rownames(sub_normDT_left) == c(inTAD_genes, outTAD_genes_left))
      stopifnot(cond1_ID %in% colnames(sub_normDT_left))
      stopifnot(cond2_ID %in% colnames(sub_normDT_left))
      
      sub_normDT_cond1_left <- sub_normDT_left[,cond1_ID]
      sub_normDT_cond2_left <- sub_normDT_left[,cond2_ID]
      
      stopifnot(nrow(sub_normDT_left) == nrow(sub_normDT_cond1_left))
      stopifnot(nrow(sub_normDT_left) == nrow(sub_normDT_cond2_left))
      stopifnot( ncol(sub_normDT_cond1_left) + ncol(sub_normDT_cond2_left) == ncol(sub_normDT_left))
      stopifnot( ncol(sub_normDT_cond1_left) == length(cond1_ID))
      stopifnot(ncol(sub_normDT_cond2_left) == length(cond2_ID))
      
      meanCorrLeft_all <- get_meanCorr_value(
        exprMatrix = sub_normDT_left, 
        inside_genes = inTAD_genes,      # inside_genes and outside_genes should be in rownames of exprMatrix
        outside_genes = outTAD_genes_left, 
        cormet = corrMeth
      )
      meanCorrLeft_cond1 <- get_meanCorr_value(
        exprMatrix = sub_normDT_cond1_left, 
        inside_genes = inTAD_genes,      # inside_genes and outside_genes should be in rownames of exprMatrix
        outside_genes = outTAD_genes_left, 
        cormet = corrMeth
      )
      meanCorrLeft_cond2 <- get_meanCorr_value(
        exprMatrix = sub_normDT_cond2_left, 
        inside_genes = inTAD_genes,      # inside_genes and outside_genes should be in rownames of exprMatrix
        outside_genes = outTAD_genes_left, 
        cormet = corrMeth
      )
    } else {
      meanCorrLeft_all <- NA
      meanCorrLeft_cond1 <- NA
      meanCorrLeft_cond2 <- NA
    }
  } else {
    meanCorr_all <- NA
    meanCorr_cond1 <- NA
    meanCorr_cond2 <- NA
    
    meanCorrRight_all <- NA
    meanCorrRight_cond1 <- NA
    meanCorrRight_cond2 <- NA
    
    meanCorrLeft_all <- NA
    meanCorrLeft_cond1 <- NA
    meanCorrLeft_cond2 <- NA
  }
  list(
    nGenes = tad_data$nGenes,
    meanCorr = meanCorr_all,
    meanCorr_cond1 = meanCorr_cond1,
    meanCorr_cond2 = meanCorr_cond2,
    nGenes_right = tad_data$nGenes_right,
    meanCorr_right = meanCorrRight_all,
    meanCorr_cond1_right = meanCorrRight_cond1,
    meanCorr_cond2_right = meanCorrRight_cond2,
    nGenes_left = tad_data$nGenes_left,
    meanCorr_left = meanCorrLeft_all,
    meanCorr_cond1_left = meanCorrLeft_cond1,
    meanCorr_cond2_left = meanCorrLeft_cond2
  )
} # end iterating over all regions
names(meanCorr_sample_around_TADs_sameNbr) <- all_regs 

outFile <- file.path(curr_outFold, "meanCorr_sample_around_TADs_sameNbr.Rdata")
save(meanCorr_sample_around_TADs_sameNbr, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

########################################################
########################################################
########################################################
txt <- paste0(startTime, "\n", Sys.time(), "\n")
printAndLog(txt, pipLogFile)
cat(paste0("*** DONE: ", script_name, "\n"))






