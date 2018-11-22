options(scipen=100)

startTime <- Sys.time()

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1)
settingF <- args[1]
stopifnot(file.exists(settingF))

pipScriptDir <- paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2")

script0_name <- "0_prepGeneData"
script1_name <- "1_runGeneDE"
script8_name <- "8c_runAllDown"
script9_name <- "9_runEmpPvalMeanTADLogFC"
script10_name <- "10_runEmpPvalMeanTADCorr"
script11_name <- "11_runEmpPvalCombined"
script17pre_name <- "170_score_auc_pval_withShuffle"
script_name <- "17ipreRatios_build_score_table_limited_pval_withShuffle"
stopifnot(file.exists(paste0(pipScriptDir, "/", script_name, ".R")))
cat(paste0("> START ", script_name,  "\n"))

# settingF <- "BUILDDT_SETTING_FILES/run_settings_buildDT.R"

source("main_settings.R")
#source("run_settings.R")
source(settingF)
source(paste0(pipScriptDir, "/", "TAD_DE_utils.R"))
# source(paste0(pipScriptDir, "/", "calc_smile_score.R"))
# source(paste0(pipScriptDir, "/", "calc_fit_ks.R"))
suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
suppressPackageStartupMessages(library(ggplot2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
suppressPackageStartupMessages(library(ggpubr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 

# registerDoMC(ifelse(SSHFS,2,20))
registerDoMC(ifelse(SSHFS,2, nCpu)) # loaded from main_settings.R

# create the directories
curr_outFold <- paste0(pipOutFold, "/", script_name)
system(paste0("mkdir -p ", curr_outFold))

pipLogFile <- paste0(pipOutFold, "/", format(Sys.time(), "%Y%d%m%H%M%S"),"_", script_name, "_logFile.txt")
system(paste0("rm -f ", pipLogFile))

processDT <- read.delim(process_file,header=T, sep="\t", stringsAsFactors=F, col.names=c("dataset", "process"))
processDT_withSamp <- read.delim(process_file_with_samp,header=T, sep="\t", stringsAsFactors=F)

stopifnot(all(processDT$dataset == processDT_withSamp$dataset))

plotType <- "svg"

ref_to_rank <- "ratioDown"

# retrieved the thresholds for which we will retrieve the % of TADs (defined in main_settings.R)
stopifnot(exists("ratioDown_thresh"))
stopifnot(exists("prodSignedRatio_thresh"))

# if the threshold for ratioDown is passed as 0.8 -> convert to 0.2 because later I used >= (1-thresh)
if(ratioDown_thresh > 0.5) ratioDown_thresh <- 1-ratioDown_thresh

stopifnot(ratioDown_thresh >= 0 & ratioDown_thresh <= 0.5)
stopifnot(prodSignedRatio_thresh >= -1 & prodSignedRatio_thresh <= 1)

# "permThresh" the quantile of permutations to take is loaded from main_settings.R

# allDown retrieved from main_settings.R
# toPlotRanked retrieved from run_settings_<>.R

#*********************************************************
# build table with following columns:

# dataset_name | smile score | slope_obs | slope_meanPerm | mean_slopeAllPerm | slope_obs/slope_meanPerm | KS pval sigFit | KS D stat sigFit | KS pval cumsum | KS D stat cumsum

#*********************************************************

txt <- paste0(toupper(script_name), "> Build score table for following datasets: ", paste0(comp_folders, collapse = ", "), "\n")
printAndLog(txt, pipLogFile)

########################################### INTERSECT FOR THE MEAN TAD CORR

#*************************** FOR DEBUG ******************************************************
# allDown = allDown[1:2]
# comp_folders <- c( "OUTPUT_FOLDER/TCGA_brca_lum_bas", "OUTPUT_FOLDER/TCGA_crc_msi_mss",
#                    "OUTPUT_FOLDER/TCGA_stad_msi_gs", "OUTPUT_FOLDER/TCGA_ucec_msi_cnl")
#********************************************************************************************

all_ratio_DT <-  foreach(curr_ds = comp_folders, .combine = 'rbind') %do% {
  
  # retrieve the ratio of ratioDown > thresh
  all_obs_ratioDown <- eval(parse(text = load(paste0(curr_ds, "/", script8_name, "/", "all_obs_ratioDown.Rdata"))))
  
  # prodSignedRatio
  all_obs_prodSignedRatio <- eval(parse(text = load(paste0(curr_ds, "/", script8_name, "/", "all_obs_prodSignedRatio.Rdata"))))
  
  ratio_ratioDown <- length(all_obs_ratioDown[all_obs_ratioDown <= ratioDown_thresh | all_obs_ratioDown >= (1-ratioDown_thresh)]) / length(all_obs_ratioDown)
  
  ratio_prodSignedRatio <- length(all_obs_prodSignedRatio[all_obs_prodSignedRatio >= prodSignedRatio_thresh]) / length(all_obs_prodSignedRatio)
  
  data.frame(dataset = basename(curr_ds),
             ratio_ratioDown = ratio_ratioDown,
             ratio_prodSignedRatio = ratio_prodSignedRatio,
           stringsAsFactors = F)
}
colnames(all_ratio_DT) <- c("dataset", paste0("ratio_ratioDown_", ratioDown_thresh), paste0("ratio_prodSignedRatio_", prodSignedRatio_thresh))

stopifnot(!any(duplicated(all_ratio_DT$dataset)))
rownames(all_ratio_DT) <- all_ratio_DT$dataset

# all_scores_DT <- foreach(curr_ds = comp_folders, .combine = 'rbind') %dopar% {
# run with do because parallelize inside plot_cumsumDiff05_AUC2()
all_scores_DT <- foreach(curr_ds = comp_folders, .combine = 'rbind') %do% {
  
  cat(paste0("**** START computing scores for: ", basename(curr_ds), "\n"))
  
  curr_ds_allratio <- eval(parse(text = load(paste0(curr_ds, "/", script17pre_name, "/", "allratio_auc_pval.Rdata"))))
  
  
  # 2) then iterate over the ratio type to get the AUC ratio
  allratio_auc <- foreach(curr_ratio_type = allDown, .combine = 'c') %do% {
    
    cat(paste0("*** START - retrieve ", curr_ratio_type, "\n"))
    
    ############################### GET THE AUC AND THE PVAL - FOR GENE PERMUTATION DATA
    toretrieve <- paste0(curr_ratio_type, "_auc_permGenes")
    stopifnot(toretrieve %in% names(curr_ds_allratio))
    auc_ratio_stat_dep05 <-   curr_ds_allratio[toretrieve]
    
    toretrieve <- paste0(curr_ratio_type, "_pval_permGenes")
    stopifnot(toretrieve %in% names(curr_ds_allratio))
    auc_ratio_stat_dep05_pval <-  curr_ds_allratio[toretrieve]
    
    toretrieve <- paste0(curr_ratio_type, "_auc_shuffTADs")
    stopifnot(toretrieve %in% names(curr_ds_allratio))
    auc_ratio_stat_dep05_shuffTADs <-  curr_ds_allratio[toretrieve]
    
    toretrieve <- paste0(curr_ratio_type, "_pval_shuffTADs")
    stopifnot(toretrieve %in% names(curr_ds_allratio))
    auc_ratio_stat_dep05_pval_shuffTADs <-  curr_ds_allratio[toretrieve]
    
    c(auc_ratio_stat_dep05, auc_ratio_stat_dep05_pval, auc_ratio_stat_dep05_shuffTADs, auc_ratio_stat_dep05_pval_shuffTADs)
    
  }
  ds_vect <- c(basename(curr_ds),  allratio_auc)
  names(ds_vect) <- c("dataset", paste0(rep(allDown,each=4), c("_auc_permGenes", "_pval_permGenes", "_auc_shuffTADs", "_pval_shuffTADs")))
  ds_vect
}

# cat(paste0("... DT built, assemble and rank. Ranking will be down according to: ", ref_to_rank, "\n"))
# dtFile <- paste0(curr_outFold, "/all_scores_DT.Rdata")
# save(all_scores_DT, file = dtFile)
# cat(paste0("... written: ", dtFile, "\n"))
# 
# dtFile <- paste0(curr_outFold, "/all_scores_DT.Rdata")
# load(dtFile)

all_scores_DT_s <- all_scores_DT
all_scores_DT <- data.frame(all_scores_DT, stringsAsFactors = F)
rownames(all_scores_DT) <- all_scores_DT$dataset
all_scores_DT$dataset <- NULL
# ensure this is not a characer
stopifnot(is.character(all_scores_DT[1,1]))
all_scores_DT <- data.matrix(all_scores_DT)
stopifnot(is.numeric(all_scores_DT[1,2]))
# for ranking the p-val: 
all_scores_DT[, grep("_pval_", colnames(all_scores_DT))] <- 1-all_scores_DT[, grep("_pval_", colnames(all_scores_DT))]
# order by decreasing since ratio is Obs/Perm
all_scores_rank_DT <- apply(-all_scores_DT,2,rank, ties="min")
colnames(all_scores_rank_DT) <- paste0(colnames(all_scores_rank_DT), "_rank")
# reverse the pval after ranking
all_scores_DT[, grep("_pval_", colnames(all_scores_DT))] <- 1-all_scores_DT[, grep("_pval_", colnames(all_scores_DT))]
all_DT <- cbind(all_scores_DT, all_scores_rank_DT)
all_DT <- all_DT[,sort(colnames(all_DT))]      

stopifnot(paste0(ref_to_rank, "_auc", "_permGenes_rank") %in% colnames(all_DT))
all_DT <- all_DT[order(all_DT[,paste0(ref_to_rank, "_auc", "_permGenes_rank")]),]

stopifnot(is.numeric(all_DT[1,1]))
# all_scores_DT <- round(all_DT, 2)
all_scores_DT <- round(all_DT, 4)
all_scores_DT <- as.data.frame(all_scores_DT)
tmpCol <- colnames(all_scores_DT)
all_scores_DT$dataset <- rownames(all_scores_DT)

all_scores_DT$process <- unlist(sapply(as.character(all_scores_DT$dataset), function(x) {
  process <- processDT$process[processDT$dataset == x ]
  process <- gsub("_", " ", process)
  if(length(process) == 0) process <- "undef"
  return(process)
}))

tmpCol <- tmpCol[-which(tmpCol == paste0(ref_to_rank, "_auc_permGenes_rank"))]
tmpCol <- tmpCol[-which(tmpCol == paste0(ref_to_rank, "_auc_permGenes"))]

tmpCol <- tmpCol[-which(tmpCol == paste0(ref_to_rank, "_pval_permGenes_rank"))]
tmpCol <- tmpCol[-which(tmpCol == paste0(ref_to_rank, "_pval_permGenes"))]

all_scores_DT <- all_scores_DT[,c("dataset", "process", 
                                  paste0(ref_to_rank, "_auc_permGenes"), paste0(ref_to_rank, "_auc_permGenes_rank"), 
                                  paste0(ref_to_rank, "_pval_permGenes"), paste0(ref_to_rank, "_pval_permGenes_rank"), tmpCol)]

# dtFile <- paste0(curr_outFold, "/rank_all_scores_DT_withPval_withShuff.txt")
# write.table(all_scores_DT, file = dtFile, row.names=F, col.names=T, quote=F, sep="\t")
# cat(paste0("... written: ", dtFile, "\n"))

all_scores_samples_DT <- all_scores_DT

all_scores_samples_DT$cond1 <- unlist(sapply(as.character(all_scores_samples_DT$dataset), function(x) {
  cond1 <- processDT_withSamp$cond1[processDT_withSamp$dataset == x ]
  cond1 <- gsub("_", " ", cond1)
  if(length(cond1) == 0) cond1 <- "undef"
  return(cond1)
}))

all_scores_samples_DT$cond2 <- unlist(sapply(as.character(all_scores_samples_DT$dataset), function(x) {
  cond2 <- processDT_withSamp$cond2[processDT_withSamp$dataset == x ]
  cond2 <- gsub("_", " ", cond2)
  if(length(cond2) == 0) cond2 <- "undef"
  return(cond2)
}))
all_scores_samples_DT$nSamp1 <- unlist(sapply(as.character(all_scores_samples_DT$dataset), function(x) {
  nSamp1 <- processDT_withSamp$nSamp1[processDT_withSamp$dataset == x ]
  nSamp1 <- gsub("_", " ", nSamp1)
  if(length(nSamp1) == 0) nSamp1 <- "undef"
  return(nSamp1)
}))
all_scores_samples_DT$nSamp2 <- unlist(sapply(as.character(all_scores_samples_DT$dataset), function(x) {
  nSamp2 <- processDT_withSamp$nSamp2[processDT_withSamp$dataset == x ]
  nSamp2 <- gsub("_", " ", nSamp2)
  if(length(nSamp2) == 0) nSamp2 <- "undef"
  return(nSamp2)
}))

all_scores_samples_DT <- all_scores_samples_DT[,c("dataset", "process", "cond1", "cond2", "nSamp1", "nSamp2", 
                                                  paste0(ref_to_rank, "_auc_permGenes"), paste0(ref_to_rank, "_auc_permGenes_rank"),  
                                                  paste0(ref_to_rank, "_pval_permGenes"), paste0(ref_to_rank, "_pval_permGenes_rank"),  
                                                  tmpCol)]
stopifnot(all_scores_samples_DT$dataset == rownames(all_scores_samples_DT))
stopifnot(all_ratio_DT$dataset == rownames(all_ratio_DT))
stopifnot(setequal(rownames(all_scores_samples_DT), rownames(all_ratio_DT)))


all_scores_ratios_samples_DT <- merge(all_scores_samples_DT, all_ratio_DT, by="dataset")


dtFile <- paste0(curr_outFold, "/all_scores_ratios_samples_DT.Rdata")
save(all_scores_ratios_samples_DT, file = dtFile)
cat(paste0("... written: ", dtFile, "\n"))


all_scores_ratios_samples_DT[, paste0("ratio_ratioDown_", ratioDown_thresh)] <- round(all_scores_ratios_samples_DT[, paste0("ratio_ratioDown_", ratioDown_thresh)] , 4)
all_scores_ratios_samples_DT[, paste0("ratio_prodSignedRatio_", prodSignedRatio_thresh)] <- round(all_scores_ratios_samples_DT[, paste0("ratio_prodSignedRatio_", prodSignedRatio_thresh)] , 4)

dtFile <- paste0(curr_outFold, "/rank_all_scores_DT_withSamp_withPval_withShuff_withRatios.txt")
write.table(all_scores_ratios_samples_DT, file = dtFile, row.names=F, col.names=T, quote=F, sep="\t")
cat(paste0("... written: ", dtFile, "\n"))



# plotFile <- file.path(curr_outFold, "rank_all_scores_DT_withPval_withShuff_withRatios.svg")
# plotTable <- ggtexttable(all_scores_DT, 
#                          rows = NULL,
#                          theme = ttheme("mCyanWhite"))
# ggsave(plotTable, file = plotFile, height = 20, width=49)
# cat(paste0("... written: ", plotFile, "\n"))                        
# 
# 
# plotFile <- file.path(curr_outFold, "rank_all_scores_DT_withSamp_withPval_withShuff_withRatios.svg")
# plotTable <- ggtexttable(all_scores_samples_DT, 
#                          rows = NULL,
#                          theme = ttheme("mCyanWhite"))
# ggsave(plotTable, file = plotFile, height = 20, width=49)
# cat(paste0("... written: ", plotFile, "\n"))                        



cat("*** DONE\n")      
cat(paste0(startTime, "\n", Sys.time(), "\n"))

