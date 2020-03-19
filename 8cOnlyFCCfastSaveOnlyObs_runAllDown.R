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
script_name <- "8cOnlyFCConlyObs_runAllDown"
#stopifnot(file.exists(paste0(pipScriptDir, "/", script_name, ".R")))
cat(paste0("> START ", script_name,  "\n"))

source("main_settings.R")
source(settingF)
source(paste0(pipScriptDir, "/", "TAD_DE_utils.R"))
source(paste0(pipScriptDir, "/", "all_down_ratio_fct.R"))

source(paste0(pipScriptDir, "/", "my_save_pigz.R")) # to use customed fastSave save.pigz()

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


}







################################****************************************************************************************
####################################################### WRITE OUTPUT
################################****************************************************************************************


txt <- paste0(startTime, "\n", Sys.time(), "\n")
printAndLog(txt, pipLogFile)

cat(paste0("*** DONE: ", script_name, "\n"))







