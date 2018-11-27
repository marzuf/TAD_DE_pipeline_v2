startTime <- Sys.time()

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1)
settingF <- args[1]
stopifnot(file.exists(settingF))

pipScriptDir <- paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline")

script0_name <- "0_prepGeneData"
script1_name <- "1_runGeneDE"
script8_name <- "8_runRatioDown"
script9_name <- "9_runEmpPvalMeanTADLogFC"
script10_name <- "10_runEmpPvalMeanTADCorr"
script11_name <- "11_runEmpPvalCombined"
script_name <- "17_build_score_table"
stopifnot(file.exists(paste0(pipScriptDir, "/", script_name, ".R")))
cat(paste0("> START ", script_name,  "\n"))

source("main_settings.R")
#source("run_settings.R")
source(settingF)
source(paste0(pipScriptDir, "/", "TAD_DE_utils.R"))
source(paste0(pipScriptDir, "/", "calc_smile_score.R"))
source(paste0(pipScriptDir, "/", "calc_fit_ks.R"))
suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) # error bar
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) # error bar

registerDoMC(ifelse(SSHFS,2, nCpu)) # loaded from main_settings.R

# create the directories
curr_outFold <- paste0(pipOutFold, "/", script_name)
system(paste0("mkdir -p ", curr_outFold))

pipLogFile <- paste0(pipOutFold, "/", script_name, "_logFile.txt")
system(paste0("rm -f ", pipLogFile))


plotType <- "svg"
myHeight <- 7
myWidth <- 10

#*********************************************************
# build table with following columns:

# dataset_name | smile score | slope_obs | slope_meanPerm | mean_slopeAllPerm | slope_obs/slope_meanPerm | KS pval sigFit | KS D stat sigFit | KS pval cumsum | KS D stat cumsum

#*********************************************************

txt <- paste0(toupper(script_name), "> Build score table for following datasets: ", paste0(comp_folders, collapse = ", "), "\n")
printAndLog(txt, pipLogFile)

########################################### INTERSECT FOR THE MEAN TAD CORR


all_scores_DT <- foreach(curr_ds = comp_folders, .combine = 'rbind') %dopar% {

    # 1) smile score
    curr_smile_score <- calc_smile_score(currOutFold = curr_ds)

    # retrieve all the stat in the same func
    all_ks_scores <- calc_all_ks_scores(curr_ds) 
          # c(cumsum_ks_pval = as.numeric(obs_perm_ks_test$p.value), cumsum_ks_D_stat = as.numeric(obs_perm_ks_test$statistic)))
          # c(obs_slope_meanPerm = obs_slope, meanPerm_slope = mean_perm_slope, fitMeanPerm_ks_pval = as.numeric(vect_ks_test$p.value), fitMeanPerm_ks_D_stat = as.numeric(vect_ks_test$statistic)))
          # c(obs_slope_meanAllPerm = obs_slope, meanAllPerm_slope = perm_mean_slope, fit_meanAllPerm_ks_pval = as.numeric(vect_ks_test$p.value), fit_meanAllPerm_ks_D_stat = as.numeric(vect_ks_test$statistic)))

    # 2) slope_obs
    obs_sigFit_slope1 <- all_ks_scores["obs_slope_meanPerm"]
    obs_sigFit_slope2 <- all_ks_scores["obs_slope_meanAllPerm"]
    stopifnot(obs_sigFit_slope1 == obs_sigFit_slope2)


    # 3) slope_meanPerm
    slope_meanPerm <- all_ks_scores["meanPerm_slope"]

    # 4) KS sigFit meanPerm
    ks_pval_meanPerm <- all_ks_scores["fitMeanPerm_ks_pval"]
    ks_dstat_meanPerm <- all_ks_scores["fitMeanPerm_ks_D_stat"]


    # 5) mean_slopeAllPerm
    slope_meanAllPerm <- all_ks_scores["meanAllPerm_slope"]


    # 6) KS sigFit mean allPerm
    ks_pval_meanAllPerm <- all_ks_scores["fit_meanAllPerm_ks_pval"]
    ks_dstat_meanAllPerm <- all_ks_scores["fit_meanAllPerm_ks_D_stat"]


    # 7) KS cumsum
    ks_pval_cumsum_permMax <- all_ks_scores["cumsum_permMax_ks_pval"]
    ks_dstat_cumsum_permMax <- all_ks_scores["cumsum_permMax_ks_D_stat"]

    ks_pval_cumsum_permMean <- all_ks_scores["cumsum_permMean_ks_pval"]
    ks_dstat_cumsum_permMean <- all_ks_scores["cumsum_permMean_ks_D_stat"]


    data.frame(dataset = curr_ds, 
                smile_score = curr_smile_score,
                obsSlope_meanPermSlope = obs_sigFit_slope2/slope_meanPerm,
                oneMinus_obsSlope_meanPermSlope = 1-obs_sigFit_slope2/slope_meanPerm,
                obsSlope_meanAllPermSlope = obs_sigFit_slope2/slope_meanAllPerm,
                oneMinus_obsSlope_meanAllPermSlope = 1-obs_sigFit_slope2/slope_meanAllPerm,               
                obs_sigFit_slope = obs_sigFit_slope2,
                meanPerm_slope = slope_meanPerm,
                meanPerm_ks_pval = ks_pval_meanPerm,
                meanPerm_ks_dstat = ks_dstat_meanPerm,
                meanAllPerm_slope = slope_meanAllPerm,
                meanAllPerm_ks_pval = ks_pval_meanAllPerm,
                meanAllPerm_ks_dstat = ks_dstat_meanAllPerm,
                cumsum_ks_pval_permMean = ks_pval_cumsum_permMean,
                cumsum_ks_dstat_permMean = ks_dstat_cumsum_permMean,
                cumsum_ks_pval_permMax = ks_pval_cumsum_permMax,
                cumsum_ks_dstat_permMax = ks_dstat_cumsum_permMax)
}

outFile <- paste0(curr_outFold, "/", "all_scores_DT.txt" )
write.table(all_scores_DT, file = outFile, col.names = T, row.names=F, sep="\t", quote=F)
cat(paste0("... written: ", outFile, "\n"))

all_scores_DT_round <- all_scores_DT
all_scores_DT_round[,2:ncol(all_scores_DT_round)] <- round(all_scores_DT_round[,2:ncol(all_scores_DT_round)],2)
outFile <- paste0(curr_outFold, "/", "all_scores_DT_rounded.txt" )
write.table(all_scores_DT_round, file = outFile, col.names = T, row.names=F, sep="\t", quote=F)
cat(paste0("... written: ", outFile, "\n"))

txt <- paste0(startTime, "\n", Sys.time(), "\n")
printAndLog(txt, pipLogFile)
cat(paste0("*** DONE: ", script_name, "\n"))

rank_DT <- all_scores_DT[,c("smile_score", "oneMinus_obsSlope_meanPermSlope", "oneMinus_obsSlope_meanAllPermSlope", "cumsum_ks_dstat_permMean", "cumsum_ks_dstat_permMax")]
# minus sign as we want decreasing
rownames(rank_DT) <- all_scores_DT$dataset
rank_DT <- apply(rank_DT , 2, function(x) rank(-x, ties="min"))
rank_DT <- as.data.frame(rank_DT)
rank_DT <- rank_DT[order(rank_DT$smile_score),]
# rank_DT <- rank_DT[order(rank_DT[,c("smile_score")]),]
rank_DT$dataset <- rownames(rank_DT)
rank_DT <- rank_DT[,c("dataset", "smile_score", "oneMinus_obsSlope_meanPermSlope", "oneMinus_obsSlope_meanAllPermSlope", "cumsum_ks_dstat_permMean", "cumsum_ks_dstat_permMax")]
colnames(rank_DT) <- c("dataset", "rank_smile_score", "rank_oneMinus_obsSlope_meanPermSlope", "rank_oneMinus_obsSlope_meanAllPermSlope", "rank_cumsum_ks_dstat_permMean", "rank_cumsum_ks_dstat_permMax")

outFile <- paste0(curr_outFold, "/", "all_score_ranks_DT.txt" )
write.table(rank_DT, file = outFile, col.names = T, row.names=F, sep="\t", quote=F)
cat(paste0("... written: ", outFile, "\n"))


############################################################################### PLOT THE COMPARISON

all_scores_DT <- all_scores_DT[order(all_scores_DT$obs_sigFit_slope),]

outFile <- paste0(curr_outFold, "/", "smile_score_obs_sigFit_slope.", plotType)
do.call(plotType, list(outFile, height = myHeight, width = myWidth))
plot(smile_score ~ obs_sigFit_slope, data=all_scores_DT, 
     xlab = "obs. slope ecdf sigmoid fit",
     ylab = "smile score",
     pch = 16, cex = 0.7,
     type="b", bty="l")
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


all_scores_DT <- all_scores_DT[order(all_scores_DT$oneMinus_obsSlope_meanPermSlope),]

outFile <- paste0(curr_outFold, "/", "smile_score_obs_meanPerm_slope.", plotType)
do.call(plotType, list(outFile, height = myHeight, width = myWidth))
plot(smile_score ~ oneMinus_obsSlope_meanPermSlope, data=all_scores_DT, 
     xlab = "1 - obs. slope/mean perm. slope ecdf sigmoid fit",
     ylab = "smile score",
     pch = 16, cex = 0.7,
     type="b", bty="l")
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- paste0(curr_outFold, "/", "cumsum_permMean_ks_D_obs_meanPerm_slope.", plotType)
do.call(plotType, list(outFile, height = myHeight, width = myWidth))
plot(cumsum_ks_dstat_permMean ~ oneMinus_obsSlope_meanPermSlope, data=all_scores_DT, 
     xlab = "1 - obs. slope/mean perm. slope",
     ylab = "KS D stat cumsum permMean ratioConcordant departure 0.5",
     pch = 16, cex = 0.7,
     type="b", bty="l")
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- paste0(curr_outFold, "/", "cumsum_permMax_ks_D_obs_meanPerm_slope.", plotType)
do.call(plotType, list(outFile, height = myHeight, width = myWidth))
plot(cumsum_ks_dstat_permMax ~ oneMinus_obsSlope_meanPermSlope, data=all_scores_DT, 
     xlab = "1 - obs. slope/mean perm. slope",
     ylab = "KS D stat cumsum permMax ratioConcordant departure 0.5",
     pch = 16, cex = 0.7,
     type="b", bty="l")
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


all_scores_DT <- all_scores_DT[order(all_scores_DT$oneMinus_obsSlope_meanAllPermSlope),]

outFile <- paste0(curr_outFold, "/", "smile_score_obs_allPerm_mean_slope.", plotType)
do.call(plotType, list(outFile, height = myHeight, width = myWidth))
plot(smile_score ~ oneMinus_obsSlope_meanAllPermSlope, data=all_scores_DT, 
     xlab = "1 - obs. slope/all perm. mean slope ecdf sigmoid fit",
     ylab = "smile score",
     pch = 16, cex = 0.7,
     type="b", bty="l")
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- paste0(curr_outFold, "/", "cumsum_permMean_ks_D_obs_allPerm_mean_slope.", plotType)
do.call(plotType, list(outFile, height = myHeight, width = myWidth))
plot(cumsum_ks_dstat_permMean ~ oneMinus_obsSlope_meanAllPermSlope, data=all_scores_DT, 
     xlab = "1 - obs. slope/all perm. mean slope",
     ylab = "KS D stat cumsum permMean ratioConcordant departure 0.5",
     pch = 16, cex = 0.7,
     type="b", bty="l")
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- paste0(curr_outFold, "/", "cumsum_permMax_ks_D_obs_allPerm_mean_slope.", plotType)
do.call(plotType, list(outFile, height = myHeight, width = myWidth))
plot(cumsum_ks_dstat_permMax ~ oneMinus_obsSlope_meanAllPermSlope, data=all_scores_DT, 
     xlab = "1 - obs. slope/all perm. mean slope",
     ylab = "KS D stat cumsum permMax ratioConcordant departure 0.5",
     pch = 16, cex = 0.7,
     type="b", bty="l")
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

