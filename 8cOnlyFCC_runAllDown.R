#!/usr/bin/Rscript

options(scipen=100)

startTime <- Sys.time()

################  USE THE FOLLOWING FILES FROM PREVIOUS STEPS
# - script0: pipeline_regionList.Rdata
# - script0: pipeline_geneList.Rdata
# - script0: rangeTADgenes.Rdata
# - script1: DE_topTable.Rdata
# - script1: DE_geneList.Rdata
# - script5: permutationsDT.Rdata
# - script5b: permutationsDT_<sizeRatio>.Rdata
# - script5c: permutationsDT_<sizeRatio>.Rdata
# - script5d: permutationsDT_shuffle.Rdata
################################################################################

################  OUTPUT
# - all_obs_<ratio_type>.Rdata 
# - <ratio_type>.Rdata if step8_for_permutGenes
# - <ratio_type>_randomGaussianList_<sizeratio>.Rdata if step8_for_randomTADsGaussian
# - <ratio_type>_randomFixSizeList_<sizeratio>.Rdata if step8_for_randomTADsFix
# - <ratio_type>_randomFixSizeList_<sizeratio>.Rdata if step8_for_randomTADsShuffle
################################################################################

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))  
suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))  

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1)
settingF <- args[1]
stopifnot(file.exists(settingF))

pipScriptDir <- paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2")

script0_name <- "0_prepGeneData"
script1_name <- "1_runGeneDE"
script5_name <- "5_runPermutationsMedian"
script5b_name <- "5b_runPermutationsRandomTADsFixSize"
script5c_name <- "5c_runPermutationsRandomTADsGaussian"
script5d_name <- "5d_runPermutationsRandomTADsShuffle"
script_name <- "8cOnlyFCC_runAllDown"
stopifnot(file.exists(paste0(pipScriptDir, "/", script_name, ".R")))
cat(paste0("> START ", script_name,  "\n"))

source("main_settings.R")
source(settingF)
source(paste0(pipScriptDir, "/", "TAD_DE_utils.R"))
source(paste0(pipScriptDir, "/", "all_down_ratio_fct.R"))

registerDoMC(ifelse(SSHFS,2, nCpu)) # loaded from main_settings.R

# create the directories
curr_outFold <- paste0(pipOutFold, "/", script_name)
system(paste0("mkdir -p ", curr_outFold))

pipLogFile <- paste0(pipOutFold, "/", format(Sys.time(), "%Y%d%m%H%M%S"),"_", script_name, "_logFile.txt")
system(paste0("rm -f ", pipLogFile))

# ADDED 27.11.2018 to check using other files
txt <- paste0("gene2tadDT_file\t=\t", gene2tadDT_file, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0("TADpos_file\t=\t", TADpos_file, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0("settingF\t=\t", settingF, "\n")
printAndLog(txt, pipLogFile)

allDown_fct <- paste0("get_", allDown, "ByRegion_v2")
names(allDown_fct) <- allDown

# if not set: default is to run for gene permutation and not for randomTADs

# ADDED 16.11.2018 to check using other files
txt <- paste0("inputDataType\t=\t", inputDataType, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0("gene2tadDT_file\t=\t", gene2tadDT_file, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0("TADpos_file\t=\t", TADpos_file, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0("settingF\t=\t", settingF, "\n")
printAndLog(txt, pipLogFile)

####****************************************************************** HARD-CODED DEFAULT SETTINGS IF NOT SET IN main_settings.R
if(!exists("step8_for_permutGenes")) step8_for_permutGenes <- TRUE
if(!exists("step8_for_randomTADsFix")) step8_for_randomTADsFix <- FALSE
if(!exists("step8_for_randomTADsGaussian")) step8_for_randomTADsGaussian <- FALSE
if(!exists("step8_for_randomTADsShuffle")) step8_for_randomTADsShuffle <- FALSE
#****************************************************************************************************************************
################################****************************************************************************************
####################################################### PREPARE INPUT DATA
################################****************************************************************************************

rangeTADsize <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/rangeTADgenes.Rdata"))))
stopifnot(rangeTADsize[1] >= minNbrGeneTAD)


# INPUT DATA
gene2tadDT <- read.delim(gene2tadDT_file, header=F, col.names = c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = F)
gene2tadDT$entrezID <- as.character(gene2tadDT$entrezID)

DE_topTable <- eval(parse(text = load(paste0(pipOutFold, "/", script1_name, "/DE_topTable.Rdata"))))
DE_geneList <- eval(parse(text = load(paste0(pipOutFold, "/", script1_name, "/DE_geneList.Rdata"))))
stopifnot(nrow(DE_topTable) == length(DE_geneList))
# stopifnot(all(DE_topTable$genes %in% names(DE_geneList)))
stopifnot(!any(duplicated(names(DE_geneList))))
stopifnot(!any(duplicated(names(DE_topTable$genes))))
### TO USE THE FUNCTIONS FOR RATIO DOWN DE_TOPTABLE$GENES => SHOULD GIVE ENTREZ ID
entrezList <- unlist(sapply(DE_topTable$genes, function(x) DE_geneList[x]))
DE_topTable$genes <- entrezList

pipeline_regionList <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/pipeline_regionList.Rdata"))))
pipeline_geneList <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/pipeline_geneList.Rdata"))))

stopifnot(all(pipeline_geneList %in% DE_topTable$genes))
## !! this is what I should not do here
# initNrow <- nrow(DE_topTable)
# DE_topTable <- DE_topTable[DE_topTable$genes %in% pipeline_geneList,]
# txt <- paste0(toupper(script_name), "> Take only filtered genes: ", nrow(DE_topTable), "/", initNrow, "\n")
# printAndLog(txt, pipLogFile)
# I filter the DE table later based on the input data

gene2tadDT <- gene2tadDT[gene2tadDT$entrezID %in% entrezList,]
initLen <- length(unique(gene2tadDT$region))
gene2tadDT <- gene2tadDT[gene2tadDT$region %in% pipeline_regionList,]
txt <- paste0(toupper(script_name), "> Take only filtered regions: ", length(unique(gene2tadDT$region)), "/", initLen, "\n")
printAndLog(txt, pipLogFile)


### !! FOR THIS 8c version -> only prodSignedRatio
# done because not run intraTADcorr for the 100'000 permutations !!!

allDown <- "prodSignedRatio"

################################****************************************************************************************
####################################################### START THE CALCULATIONS 
################################****************************************************************************************


if(step8_for_permutGenes){
  cat("... load permutation data ...\n")
  permutationsDT <- eval(parse(text = load(paste0(pipOutFold, "/", script5_name, "/permutationsDT.Rdata"))))
  ### TAKE ONLY THE GENES FOR WHICH A LOG FC VALUE IS AVAILABLE (I.E. THE ONES USED FOR DE ANALYSIS)
  initNrow <- nrow(permutationsDT)
  permutationsDT <- permutationsDT[which(rownames(permutationsDT) %in% DE_topTable$genes),, drop=FALSE]
  txt <- paste0(toupper(script_name), "> Take only the genes for which logFC value is available (the ones used in DE analysis): ", nrow(permutationsDT), "/", initNrow, "\n")
  printAndLog(txt, pipLogFile)
  stopifnot(nrow(permutationsDT) == initNrow)
  
  all_regions <- sort(unique(as.character(permutationsDT[,2])))
  
  stopifnot(all(rownames(permutationsDT) %in% pipeline_geneList))
  stopifnot(nrow(permutationsDT) == length(pipeline_geneList))
  
  if(useTADonly) {
    if(! all (regexpr("_TAD",  gene2tadDT$region[gene2tadDT$entrezID %in% rownames(permutationsDT)]) > 0 )) {
      stop("make not sense to filter TAD genes after permutations if permutations were run with genes belonging to non-TAD regions\n")
    }
    initLen <- length(all_regions)
    all_regions <- all_regions[grep("_TAD", all_regions)]
    # if want to use only TADs, make more sens if permutations were run without the TADs
    txt <- paste0(toupper(script_name), "> Take only TAD regions: ", length(all_regions), "/", initLen, "\n")
    printAndLog(txt, pipLogFile)    
    if(length(all_regions) < initLen){
      stop("make not sense to filter TAD regions after permutations if permutations were run with genes belonging to non-TAD regions\n")
    }
  }
}

### START FOR observed and permutGenes
# for permutGene -> I can use the same DE_topTable because it is the same set of genes each time
# for the other random data, filtering of DE_topTable is done within get_statFromShuffle_list_para_list

for(curr_ratio_type in allDown) {
  cat(paste0("... start calculate ", curr_ratio_type, "/TAD\n"))
  DE_topTable_obs <- DE_topTable[DE_topTable$genes %in% pipeline_geneList,]
  all_obs_ratio <- do.call(allDown_fct[curr_ratio_type], list(g2TADdt=gene2tadDT, DEdt=DE_topTable_obs))
  txt <- paste0(toupper(script_name), "> Number of observed ", curr_ratio_type, " computed: ", length(all_obs_ratio), "\n")
  printAndLog(txt, pipLogFile)
  # assign(curr_ratio_type,get("b"))
  save(all_obs_ratio, file = paste0(curr_outFold, "/all_obs_", curr_ratio_type, ".Rdata"))
  cat(paste0("... written: ",  paste0(curr_outFold, "/all_obs_", curr_ratio_type, ".Rdata"), "\n"))
  if(step8_for_permutGenes){
    cat("... start calculate permutations ", curr_ratio_type, "/TAD\n")
    ratio_permDT <- get_statFromShuffle_para(DEdt=DE_topTable_obs, shuffData=permutationsDT, stat_fct=allDown_fct[curr_ratio_type], ncpu=nCpu, TADonly=F)
    txt <- paste0(toupper(script_name), "> Number of permutations for which ", curr_ratio_type," computed: ", ncol(ratio_permDT), "\n")
    printAndLog(txt, pipLogFile)
    txt <- paste0(toupper(script_name), "> Number of regions for which ", curr_ratio_type, " computed: ", nrow(ratio_permDT), "\n")
    printAndLog(txt, pipLogFile)
    save(ratio_permDT, file = paste0(curr_outFold, "/", curr_ratio_type, "_permDT.Rdata"))
    cat(paste0("... written: ",  paste0(curr_outFold, "/", curr_ratio_type, "_permDT.Rdata"), "\n"))
  }
}

### START FOR partition fixSize

if(step8_for_randomTADsFix){
  for(sizeRatio in percentTADsizeFix){
    cat("... load randomTADsFixSize data ...\n")
    randomFixSizeDT <- eval(parse(text = load(paste0(pipOutFold, "/", script5b_name, "/permutationsList_fixSize_", gsub("\\.","", sizeRatio), ".Rdata"))))
    
    ### TAKE ONLY THE GENES FOR WHICH A LOG FC VALUE IS AVAILABLE (I.E. THE ONES USED FOR DE ANALYSIS)
    # initNrow <- nrow(randomFixSizeDT)
    # randomFixSizeDT <- randomFixSizeDT[which(rownames(randomFixSizeDT) %in% DE_topTable$genes),,drop=FALSE]
    # txt <- paste0(toupper(script_name), "> Take only the genes for which logFC value is available (the ones used in DE analysis): ", nrow(randomFixSizeDT), "/", initNrow, "\n")
    # printAndLog(txt, pipLogFile)
    # stopifnot(nrow(randomFixSizeDT) > 0)
    
    # does not make sense to filter the TADs for this kind of permutation !!!
    for(curr_ratio_type in allDown) {
      cat("... start calculate randomTADsFixSize ", curr_ratio_type, "/TAD\n")
      ratio_fixSizeList <- get_statFromShuffle_list_para_list(DEdt=DE_topTable, shuffData=randomFixSizeDT, stat_fct=allDown_fct[curr_ratio_type], ncpu=nCpu, TADonly=F)
      txt <- paste0(toupper(script_name), "> Number of randomTADsFixSize for which ", curr_ratio_type," computed: ", length(ratio_fixSizeList), "\n")
      printAndLog(txt, pipLogFile)
      
      
      # # !!! but I need to remove the regions with too much or too many genes !!!
      # # for(i in 1:length(ratio_fixSizeList)) {
      # filter_ratio_fixSizeList <- lapply(1:length(ratio_fixSizeList), function(i) {
      #   nGenes <- table(permutationsDT[,i])
      #   upperLimit <- as.numeric(quantile(as.numeric(nGenes), probs=maxQuantGeneTAD))
      #   toKeepTADs <- names(nGenes)[as.numeric(nGenes) >= minNbrGeneTAD & as.numeric(nGenes) <= upperLimit ]
      #   beforeFilter <- ratio_fixSizeList[[i]]
      #   stopifnot(toKeepTADs %in% names(beforeFilter))
      #   beforeFilter[toKeepTADs]
      # })
      # ratio_fixSizeList <- filter_ratio_fixSizeList
      
      stopifnot(unlist(lapply(ratio_fixSizeList, function(x) !any(is.na(x)))))
      save(ratio_fixSizeList, file = paste0(curr_outFold, "/", curr_ratio_type, "_randomFixSizeList_", gsub("\\.","", sizeRatio), ".Rdata"))
      cat(paste0("... written: ",  paste0(curr_outFold, "/", curr_ratio_type, "_randomFixSizeList_", gsub("\\.","", sizeRatio), ".Rdata"), "\n"))
    }
  }
}


### START FOR partition gaussian

if(step8_for_randomTADsGaussian){
  for(sizeRatio in percentTADsizeGaussian) {
    cat("... load randomTADsGaussian data ...\n")
    randomGaussianDT <- eval(parse(text = load(paste0(pipOutFold, "/", script5c_name, "/permutationsList_gaussian_", gsub("\\.","", sizeRatio), ".Rdata"))))
    
    ### TAKE ONLY THE GENES FOR WHICH A LOG FC VALUE IS AVAILABLE (I.E. THE ONES USED FOR DE ANALYSIS)
    # initNrow <- nrow(randomGaussianDT)
    # randomGaussianDT <- randomGaussianDT[which(rownames(randomGaussianDT) %in% DE_topTable$genes),,drop=FALSE]
    # txt <- paste0(toupper(script_name), "> Take only the genes for which logFC value is available (the ones used in DE analysis): ", nrow(randomGaussianDT), "/", initNrow, "\n")
    # printAndLog(txt, pipLogFile)
    # stopifnot(nrow(randomGaussianDT) > 0)
    
    # does not make sense to filter the TADs for this kind of permutation !!!
    for(curr_ratio_type in allDown) {
      cat("... start calculate randomTADsGaussian ", curr_ratio_type, "/TAD\n")
      ratio_gaussianList <- get_statFromShuffle_list_para_list(DEdt=DE_topTable, shuffData=randomGaussianDT, stat_fct=allDown_fct[curr_ratio_type], ncpu=nCpu, TADonly=F)
      txt <- paste0(toupper(script_name), "> Number of randomTADsGaussian for which ", curr_ratio_type," computed: ", length(ratio_gaussianList), "\n")
      printAndLog(txt, pipLogFile)
      
      
      # !!! but I need to remove the regions with too much or too many genes !!!
      # for(i in 1:length(ratio_gaussianList)) {
      # filter_ratio_gaussianList <- lapply(1:length(ratio_gaussianList), function(i) {
      #   nGenes <- table(permutationsDT[,i])
      #   upperLimit <- as.numeric(quantile(as.numeric(nGenes), probs=maxQuantGeneTAD))
      #   toKeepTADs <- names(nGenes)[as.numeric(nGenes) >= minNbrGeneTAD & as.numeric(nGenes) <= upperLimit ]
      #   beforeFilter <- ratio_gaussianList[[i]]
      #   stopifnot(toKeepTADs %in% names(beforeFilter))
      #   beforeFilter[toKeepTADs]
      # })
      # ratio_gaussianList <- filter_ratio_gaussianList
      
      
      stopifnot(unlist(lapply(ratio_gaussianList, function(x) !any(is.na(x)))))
      save(ratio_gaussianList, file = paste0(curr_outFold, "/", curr_ratio_type, "_randomGaussianList_", gsub("\\.","", sizeRatio), ".Rdata"))
      cat(paste0("... written: ",  paste0(curr_outFold, "/", curr_ratio_type, "_randomGaussianList_", gsub("\\.","", sizeRatio), ".Rdata"), "\n"))
    }
  }
}



### START FOR partition shuffle

if(step8_for_randomTADsShuffle){
  cat("... load randomTADsShuffle data ...\n")
  randomShuffleDT <- eval(parse(text = load(paste0(pipOutFold, "/", script5d_name, "/permutationsList_shuffle.Rdata"))))
  
  if(useFilterSizeData)
    stopifnot(sapply(randomShuffleDT, function(x) range(table(x))[1] >= rangeTADsize[1] & range(table(x))[2] <= rangeTADsize[2]))
  
  ### TAKE ONLY THE GENES FOR WHICH A LOG FC VALUE IS AVAILABLE (I.E. THE ONES USED FOR DE ANALYSIS)
  # initNrow <- nrow(randomShuffleDT)
  # randomShuffleDT <- randomShuffleDT[which(rownames(randomShuffleDT) %in% DE_topTable$genes),,drop=FALSE]
  # txt <- paste0(toupper(script_name), "> Take only the genes for which logFC value is available (the ones used in DE analysis): ", nrow(randomShuffleDT), "/", initNrow, "\n")
  # printAndLog(txt, pipLogFile)
  # stopifnot(nrow(randomShuffleDT) > 0)
  # # does not make sense to filter the TADs for this kind of permutation !!!
  
  for(curr_ratio_type in allDown) {
    cat("... start calculate randomTADsShuffle ", curr_ratio_type, "/TAD\n")
    # ensure that the regions are filtered correctly
    
    ratio_shuffleList <- get_statFromShuffle_list_para_list(DEdt=DE_topTable, shuffData=randomShuffleDT, stat_fct=allDown_fct[curr_ratio_type], ncpu=nCpu, TADonly=F)
    txt <- paste0(toupper(script_name), "> Number of randomTADsShuffle for which ", curr_ratio_type," computed: ", length(ratio_shuffleList), "\n")
    printAndLog(txt, pipLogFile)
    # !!! but I need to remove the regions with too much or too many genes !!!
    # I don't need to filter anymore because the filtering based on 
    # for(i in 1:length(ratio_shuffleList)) {
    # filter_ratio_shuffleList <- lapply(1:length(ratio_shuffleList), function(i) {
    #   nGenes <- table(permutationsDT[,i])
    #   upperLimit <- as.numeric(quantile(as.numeric(nGenes), probs=maxQuantGeneTAD))
    #   toKeepTADs <- names(nGenes)[as.numeric(nGenes) >= minNbrGeneTAD & as.numeric(nGenes) <= upperLimit ]
    #   beforeFilter <- ratio_shuffleList[[i]]
    #   stopifnot(toKeepTADs %in% names(beforeFilter))
    #   beforeFilter[toKeepTADs]
    # })
    # ratio_shuffleList <- filter_ratio_shuffleList
    
    stopifnot(unlist(lapply(ratio_shuffleList, function(x) !any(is.na(x)))))
    save(ratio_shuffleList, file = paste0(curr_outFold, "/", curr_ratio_type, "_randomShuffleList.Rdata"))
    cat(paste0("... written: ",  paste0(curr_outFold, "/", curr_ratio_type, "_randomShuffleList.Rdata"), "\n"))
  }
}


################################****************************************************************************************
####################################################### WRITE OUTPUT
################################****************************************************************************************


txt <- paste0(startTime, "\n", Sys.time(), "\n")
printAndLog(txt, pipLogFile)

cat(paste0("*** DONE: ", script_name, "\n"))







