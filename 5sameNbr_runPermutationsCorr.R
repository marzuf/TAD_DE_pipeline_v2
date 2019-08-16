#!/usr/bin/Rscript

options(scipen=100)

startTime <- Sys.time()

#### UPDATE: do not take raw counts but fpkm data !!!

################  USE THE FOLLOWING FILES FROM PREVIOUS STEPS
# - script0: pipeline_regionList.Rdata
# - script0: pipeline_geneList.Rdata
################################################################################

################  OUTPUT
# - sample_around_TADs_sameNbr.Rdata
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
script_name <- "5sameNbr_runPermutationsCorr"
stopifnot(file.exists(paste0(pipScriptDir, "/", script_name, ".R")))
cat(paste0("> START ", script_name,  "\n"))

source("main_settings.R")
source(settingF)
source(paste0(pipScriptDir, "/", "TAD_DE_utils.R"))

# create the directories
curr_outFold <- paste0(pipOutFold, "/", script_name)
system(paste0("mkdir -p ", curr_outFold))

pipLogFile <- paste0(pipOutFold, "/", format(Sys.time(), "%Y%d%m%H%M%S"),"_", script_name, "_logFile.txt")
system(paste0("rm -f ", pipLogFile))

# ADDED 16.11.2018 to check using other files
txt <- paste0("inputDataType\t=\t", inputDataType, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0("gene2tadDT_file\t=\t", gene2tadDT_file, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0("TADpos_file\t=\t", TADpos_file, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0("settingF\t=\t", settingF, "\n")
printAndLog(txt, pipLogFile)

suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
registerDoMC(nCpu)

### RETRIEVE GENE-TO-TAD ASSIGNMENT
stopifnot(file.exists(gene2tadDT_file))
g2t_DT <- read.delim(gene2tadDT_file, header=F, col.names = c("entrezID",  "chromo", "start", "end", "region"), stringsAsFactors = FALSE)
g2t_DT$entrezID <- as.character(g2t_DT$entrezID)

### RETRIEVE THE TAD POSITIONS
stopifnot(file.exists(TADpos_file))
tadpos_DT <- read.delim(TADpos_file, header=F, col.names=c("chromo", "region", "start", "end"), stringsAsFactors = FALSE)
stopifnot(is.numeric(tadpos_DT$start))
stopifnot(is.numeric(tadpos_DT$end))
tadpos_DT <- tadpos_DT[grepl("_TAD", tadpos_DT$region),,drop=FALSE] 

### KEEP ONLY THE TADs USED IN THE PIPELINE
tadListFile <- file.path(pipOutFold, script0_name, "pipeline_regionList.Rdata")
stopifnot(file.exists(tadListFile))
pipeline_tadList <- eval(parse(text = load(tadListFile))) # not adjusted
stopifnot(pipeline_tadList %in% tadpos_DT$region)
tadpos_DT <- tadpos_DT[tadpos_DT$region %in% pipeline_tadList,]
stopifnot(!duplicated(pipeline_tadList))

### RETRIEVE THE GENES USED IN THE PIPELINE - script0
geneListFile <- file.path(pipOutFold, script0_name, "pipeline_geneList.Rdata")
stopifnot(file.exists(geneListFile))
pipeline_geneList <- eval(parse(text = load(geneListFile))) # not adjusted
stopifnot(pipeline_geneList %in% g2t_DT$entrezID)
# stopifnot(names(pipeline_geneList) %in% g2t_DT$entrezID) -> FALSE
g2t_DT <- g2t_DT[as.character(g2t_DT$entrezID) %in% as.character(pipeline_geneList),,drop=FALSE]
stopifnot(length(pipeline_geneList) == nrow(g2t_DT))
stopifnot(g2t_DT$entrezID %in% pipeline_geneList)
g2t_DT$chromo <- as.character(g2t_DT$chromo)

stopifnot(g2t_DT$entrezID %in% pipeline_geneList)
stopifnot(grepl("TAD", g2t_DT$region))

tadpos_DT$mid_pos <- (tadpos_DT$start+tadpos_DT$end)/2
g2t_DT$mid_pos <- (g2t_DT$start+g2t_DT$end)/2

sample_around_TADs <- foreach(reg = pipeline_tadList) %dopar% {
  cat("...... start TAD : \t", reg, "\n")
  
  curr_chromo <- as.character(tadpos_DT$chromo[tadpos_DT$region == reg])
  
  curr_start <- (tadpos_DT$start[tadpos_DT$region == reg])
  stopifnot(is.numeric(curr_start))
  
  curr_end <- (tadpos_DT$end[tadpos_DT$region == reg])
  stopifnot(is.numeric(curr_end))

  stopifnot(length(curr_chromo) == 1)
  stopifnot(length(curr_start) == 1)
  stopifnot(length(curr_end) == 1)
  
  curr_midPos <- (curr_start+curr_end)/2
  stopifnot(curr_midPos == tadpos_DT$mid_pos[tadpos_DT$region == reg])
  
  reg_genes <- g2t_DT$entrezID[g2t_DT$region == reg]
  stopifnot(length(reg_genes) > 0)
  
  curr_nGenes <- length(reg_genes)
  
  # !!! EXTRACT GENES BASED ON START POSITION RELATIVE TO BD
  # !!! SMALLER THAN / GREATER *OR EQUAL* THAN BD POSITION (smaller not equal otherwise genes could come twice)
  
  curr_g2t <- g2t_DT[g2t_DT$chromo == curr_chromo,,drop=FALSE]
  stopifnot(nrow(curr_g2t) > 0)
  
  stopifnot(is.numeric(curr_g2t$start), is.numeric(curr_g2t$end))
  curr_g2t <- curr_g2t[order(curr_g2t$start, curr_g2t$end),,drop=FALSE]
  
  # distance to TAD center
  curr_genesOutsideDT <- curr_g2t[ !curr_g2t$entrezID %in% reg_genes,,drop=FALSE]
  stopifnot(nrow(curr_genesOutsideDT) > 0)
    
  #>>> take the same number of genes, either on left or right
  curr_genesOutsideDT$distToTAD <- abs(curr_genesOutsideDT$mid_pos - curr_midPos)
  curr_genesOutsideDT <- curr_genesOutsideDT[order(curr_genesOutsideDT$distToTAD, decreasing=F),,drop=FALSE]
  
  sample_around_genes <- curr_genesOutsideDT$entrezID[1:curr_nGenes]
  
  all_dist <- curr_genesOutsideDT$distToTAD[1:curr_nGenes]
  stopifnot(!is.na(all_dist))
  
  stopifnot(!is.na(curr_genesOutsideDT))
  
  stopifnot(length(all_dist) == length(sample_around_genes))
  
  #>>> take the same number of genes on the left
  curr_genesOutsideDT_left <- curr_genesOutsideDT[curr_genesOutsideDT$mid_pos < curr_midPos,]  # SMALLER
  curr_nGenes_left <- min(curr_nGenes, nrow(curr_genesOutsideDT_left)) 
  
  if(curr_nGenes_left > 0) {
    sample_around_genes_left <- curr_genesOutsideDT_left$entrezID[1:curr_nGenes_left]
    all_dist_left <- curr_genesOutsideDT_left$distToTAD[1:curr_nGenes_left]  
    stopifnot(!is.na(all_dist_left))
    ### ONLY CONSIDER IN COEXPR GENES USED IN PIPELINE ???
    stopifnot(sample_around_genes_left %in% pipeline_geneList)
    stopifnot(!sample_around_genes_left %in% reg_genes)
    
  } else {
    sample_around_genes_left <- character(0)
    all_dist_left <- c()
  }
  
  curr_genesOutsideDT_right <- curr_genesOutsideDT[curr_genesOutsideDT$mid_pos >= curr_midPos,]  # BIGGER OR EQUAL
  curr_nGenes_right <- min(curr_nGenes, nrow(curr_genesOutsideDT_right)) 
  
  if(curr_nGenes_right > 0) {
    sample_around_genes_right <- curr_genesOutsideDT_right$entrezID[1:curr_nGenes_right]
    all_dist_right <- curr_genesOutsideDT_right$distToTAD[1:curr_nGenes_right]
    stopifnot(!is.na(all_dist_right))
    
    stopifnot(sample_around_genes_right %in% pipeline_geneList)
    stopifnot(!sample_around_genes_right %in% reg_genes)
    
  } else {
    sample_around_genes_right <- character(0)
    all_dist_right <- c()
    
  }
  
  stopifnot(length(sample_around_genes) == curr_nGenes)
  stopifnot(length(sample_around_genes_left) <= curr_nGenes)
  stopifnot(length(sample_around_genes_right) <= curr_nGenes)
  
  stopifnot(length(all_dist) == length(sample_around_genes))
  stopifnot(length(all_dist_right) == length(sample_around_genes_right))
  stopifnot(length(all_dist_left) == length(sample_around_genes_left))
  
  ### ONLY CONSIDER IN COEXPR GENES USED IN PIPELINE ???
  stopifnot(sample_around_genes %in% pipeline_geneList)
  
  list(
      tad_genes = reg_genes,
      
       genes = sample_around_genes,
       nGenes = length(sample_around_genes),
       minDist = min(all_dist),
       maxDist = max(all_dist),
      
      genes_left = sample_around_genes_left,
      nGenes_left = curr_nGenes_left,
      minDist_left = min(all_dist_left),
      maxDist_left = max(all_dist_left),
      
      genes_right = sample_around_genes_right,
      nGenes_right =  curr_nGenes_right,
      minDist_right = min(all_dist_right),
      maxDist_right = max(all_dist_right)
      )
} # end foreach-iterating over TADs
names(sample_around_TADs) <- pipeline_tadList
sample_around_TADs_sameNbr <- sample_around_TADs

outFile <- file.path(curr_outFold, "sample_around_TADs_sameNbr.Rdata")
save(sample_around_TADs_sameNbr, file=outFile)
cat(paste0("... written: ", outFile, "\n"))

########################################################
########################################################
########################################################
txt <- paste0(startTime, "\n", Sys.time(), "\n")
printAndLog(txt, pipLogFile)
cat(paste0("*** DONE: ", script_name, "\n"))















