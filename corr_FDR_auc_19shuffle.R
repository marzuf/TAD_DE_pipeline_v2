library(foreach)
library(doMC)

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

registerDoMC(ifelse(SSHFS, 2, 20))

pipScriptDir <- paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2")

script170_name <- "170_score_auc_pval_withShuffle"
script19_name <- "19shuffle_SAM_emp_measurement"

cat(paste0("> Rscript cmp_19_SAM.R\n"))

# cat(paste0("setDir = ", setDir, "\n"))
# source("main_settings.R") # setDir is the main_settings not in run_settings
source(paste0(pipScriptDir, "/", "TAD_DE_utils.R"))

caller <- "TopDom"

mainOutFold <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller), "OUTPUT_FOLDER")

outFold <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller), "corr_FDR_auc_19shuffle")
system(paste0("mkdir -p ", outFold))

myHeight <- 10
myWidth <- 12

all_datasets <- list.files(mainOutFold, full.names = TRUE)


all_datasets_list_FDR <- foreach(curr_dataset = all_datasets) %dopar% {
  eval(parse(text = load(file.path(curr_dataset, script19_name, "empFDR_list.Rdata"))))
}
names(all_datasets_list_FDR) <- basename(all_datasets)

all_datasets_list_AUC <- foreach(curr_dataset = all_datasets) %dopar% {
  eval(parse(text = load(file.path(curr_dataset, script170_name, "allratio_auc_pval.Rdata"))))
}
names(all_datasets_list_AUC) <- basename(all_datasets)

all_ratios <- c( "prodSignedRatio", "ratioDown")

# all_ratios <- all_ratios[1]

for(curr_ratio in all_ratios) {
  
  cat("> START", curr_ratio, "\n")
  
  
  # retrieve the AUC permGenes for this ratio
  ratio_auc_permGenes <- sapply(all_datasets_list_AUC, function(x) x[[paste0(curr_ratio, "_auc_permGenes")]])
  
  # retrieve the AUC permGenes for this ratio
  ratio_auc_shuffTADs <- sapply(all_datasets_list_AUC, function(x) x[[paste0(curr_ratio, "_auc_shuffTADs")]])
  
  
  # retrieve the slope empFDR for this ratio
  slope_FDR <- sapply(all_datasets_list_FDR, function(x) x[[paste0("slopeEmpFDR_", curr_ratio)]])
  
  
  intersectReg <- Reduce(intersect, list(names(ratio_auc_permGenes), names(ratio_auc_shuffTADs), names(slope_FDR)))
  
  stopifnot(length(intersectReg) > 0)
  
  ratio_auc_permGenes <- ratio_auc_permGenes[intersectReg]
  ratio_auc_shuffTADs <- ratio_auc_shuffTADs[intersectReg]
  slope_FDR <- slope_FDR[intersectReg]
  
  ######## PLOT SLOPE FDR WITH AUC PERMGENES
  outFile <- paste0(outFold, "/", curr_ratio, "_slopeEmpFDR_ratioAUCpermGenes.svg")
  svg(outFile, height = myHeight, width = myWidth)
  plot(x = slope_FDR,
       y = ratio_auc_permGenes,
       main = paste0("Emp. FDR slope vs. AUC permGenes"),
       xlab = paste0(curr_ratio, " - slope emp. FDR shuffleTADs"),
       ylab = paste0(curr_ratio , " - ratio AUC permGenes"),
       cex = 0.7, pch = 16)
  mtext(paste0(caller, " - ", curr_ratio, " - all datasets"), font=3)
  corText <- paste0("PCC = ", sprintf("%.4f", round(cor(x=ratio_auc_permGenes, y=slope_FDR), 4)), "\nn = ", length(intersectReg))
  legend("topright", legend = corText, bty="n")
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
  ######## PLOT SLOPE FDR WITH AUC SHUFFTADS
  outFile <- paste0(outFold, "/", curr_ratio, "_slopeEmpFDR_ratioAUCshuffTADs.svg")
  svg(outFile, height = myHeight, width = myWidth)
  plot(x = slope_FDR,
       y = ratio_auc_shuffTADs,
       main = paste0("Emp. FDR slope vs. AUC shuffTADs"),
       xlab = paste0(curr_ratio, " - slope emp. FDR shuffleTADs"),
       ylab = paste0(curr_ratio , " - ratio AUC shuffTADs"),
       cex = 0.7, pch = 16)
  mtext(paste0(caller, " - ", curr_ratio, " - all_datasets"), font=3)
  corText <- paste0("PCC = ", sprintf("%.4f", round(cor(x=ratio_auc_shuffTADs, y=slope_FDR), 4)), "\nn = ", length(intersectReg))
  legend("topright", legend = corText, bty="n")
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  
}

