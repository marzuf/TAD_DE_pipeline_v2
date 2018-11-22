library(foreach)
library(doMC)

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

registerDoMC(ifelse(SSHFS, 2, 20))

# script0_name <- "0_prepGeneData"
# script3_name <- "3_runMeanTADLogFC"
# script4_name <- "4_runMeanTADCorr"
# script6_name <- "6_runPermutationsMeanLogFC"
# script7_name <- "7_runPermutationsMeanTADCorr"
# script8_name <- "8c_runAllDown"

pipScriptDir <- paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2")

script170_name <- "170_score_auc_pval_withShuffle"
script20_name <- "20_getSignifCombined"
script21_name <- "21_getSignifFDR"

cat(paste0("> Rscript cmp_ratioAUC_nbrSignif.R\n"))

caller <- "TopDom"


pipScriptDir <- paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2")


# cat(paste0("setDir = ", setDir, "\n"))
# source("main_settings.R") # setDir is the main_settings not in run_settings
source(paste0(pipScriptDir, "/", "TAD_DE_utils.R"))

caller <- "TopDom"

mainOutFold <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller), "OUTPUT_FOLDER")

outFold <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller), "cmp_ratioAUC_nbrSignif")
system(paste0("mkdir -p ", outFold))

myHeight <- 10
myWidth <- 12

all_datasets <- list.files(mainOutFold, full.names = TRUE)

all_datasets_list_FDR <- foreach(curr_dataset = all_datasets, .combine="rbind") %dopar% {
  
  all_ratiosAUC <- eval(parse(text = load(file.path(curr_dataset, script170_name, "allratio_auc_pval.Rdata"))))
  
  all_empFDR <- eval(parse(text = load(file.path(curr_dataset, script21_name, "all_empFDR_signifThresh.Rdata"))))
  
  data.frame(
  dataset = basename(curr_dataset),
  thresh_adjPval = eval(parse(text = load(file.path(curr_dataset, script20_name, "adjCombinedPval_thresh.Rdata")))),
  thresh_empFDR = all_empFDR["desired_empFDR_thresh"],
  nbrSignif_adjPvalCombined = length(eval(parse(text = load(file.path(curr_dataset, script20_name, "signifTADs.Rdata"))))),
  nbrSignif_empFDR = length(eval(parse(text = load(file.path(curr_dataset, script21_name, "signifTADs.Rdata"))))),
  ratioDown_ratioAUC_permGenes = all_ratiosAUC["ratioDown_auc_permGenes"],
  ratioDown_ratioAUC_shuffTADs = all_ratiosAUC["ratioDown_auc_shuffTADs"],
  prodSignedRatio_ratioAUC_permGenes = all_ratiosAUC["prodSignedRatio_auc_permGenes"],
  prodSignedRatio_ratioAUC_shuffTADs = all_ratiosAUC["prodSignedRatio_auc_shuffTADs"],
  stringsAsFactors = FALSE)
  
}


