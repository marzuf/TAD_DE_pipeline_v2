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

script17_name <- "17ipreRatios_build_score_table_limited_pval_withShuffle"

cat(paste0("> Rscript cmp_ratioAUC_ratioThresh.R\n"))

# cat(paste0("setDir = ", setDir, "\n"))
# source("main_settings.R") # setDir is the main_settings not in run_settings
# source(paste0(pipScriptDir, "/", "TAD_DE_utils.R"))

caller <- "TopDom"

dtFile <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller), "BUILDDT_OUTPUT_FOLDER", script17_name, "all_scores_ratios_samples_DT.Rdata")
stopifnot(file.exists(dtFile))

outFold <- file.path(setDir, paste0("/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2_", caller), "cmp_ratioAUC_ratioThresh")
system(paste0("mkdir -p ", outFold))

myHeight <- 10
myWidth <- 10

ratio_DT <- eval(parse(text = load(dtFile)))

ratioDownThresh <- gsub("ratio_ratioDown_", "", colnames(ratio_DT)[grep("ratio_ratioDown_", colnames(ratio_DT))])
prodSignedRatioThresh <- gsub("ratio_prodSignedRatio_", "", colnames(ratio_DT)[grep("ratio_prodSignedRatio_", colnames(ratio_DT))])

all_ratios <- c("ratioDown", "prodSignedRatio")
all_perms <- c("permGenes", "shuffTADs")

for(curr_ratio in all_ratios) {
  for(curr_perm in all_perms) {
    curr_thresh <- ifelse(curr_ratio == "ratioDown", ratioDownThresh,
                          ifelse(curr_ratio == "prodSignedRatio", prodSignedRatioThresh, NA))
    stopifnot(!is.na(curr_thresh))
    # y-axis = ratio of the AUC
    curr_y <- ratio_DT[, paste0(curr_ratio, "_auc_", curr_perm)] 
    # x-axis = ratio of the TADs
    curr_x <- ratio_DT[, paste0("ratio_", curr_ratio, "_", curr_thresh)]
    
    outFile <- file.path(outFold, paste0(curr_ratio, "AUCratio_vs_ratioTADs_", curr_perm, ".svg") )
    svg(outFile, height = myHeight, width = myWidth)
    plot(x = curr_x,
         y = curr_y,
         main = paste0("AUC ratio vs. ratio TADs thresh."),
         xlab = paste0("Ratio TADs ", curr_ratio, " <= ", curr_thresh, " || ratioDown >= ", (1-as.numeric(as.character(curr_thresh))) ),
         ylab = paste0("AUC ratio - ", curr_ratio, " - ", curr_perm),
         cex = 0.7, pch=16)
    mtext(side = 3, text = paste0(caller, " - ", curr_perm, " - ", curr_ratio))
    text(x= curr_x, y = curr_y, labels = ratio_DT[, "dataset"], cex = 0.7)
    corText <- paste0("PCC = ", sprintf("%.4f", round(cor(x=curr_x, y=curr_y), 4)))
    legend("topright", bty="n", legend = corText)
    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
  }
}


cat("*** DONE\n")







