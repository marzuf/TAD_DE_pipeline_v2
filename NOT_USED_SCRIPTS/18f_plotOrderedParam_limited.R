startTime <- Sys.time()

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

# Rscript 18f_plotOrderedParam_limited.R BUILDDT_SETTING_FILES/run_settings_buildDT.R OUTPUT_FOLDER/SCORE_DT/18_plotOrderedParam/rank_scores_to_plot_DT.txt

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)
settingF <- args[1]
scoreF <- args[2]
stopifnot(file.exists(settingF))
stopifnot(file.exists(scoreF))

pipScriptDir <- paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline")

script0_name <- "0_prepGeneData"
script1_name <- "1_runGeneDE"
script8_name <- "8c_runAllDown"
script9_name <- "9_runEmpPvalMeanTADLogFC"
script10_name <- "10_runEmpPvalMeanTADCorr"
script11_name <- "11_runEmpPvalCombined"
script_name <- "18f_plotOrderedParam_limited"
stopifnot(file.exists(paste0(pipScriptDir, "/", script_name, ".R")))
cat(paste0("> START ", script_name,  "\n"))

# settingF <- "BUILDDT_SETTING_FILES/run_settings_buildDT.R"
# scoreF <- "OUTPUT_FOLDER/SCORE_DT/18_plotOrderedParam/rank_scores_to_plot_DT.txt"

source("main_settings.R")
#source("run_settings.R")
source(settingF)
source(paste0(pipScriptDir, "/", "TAD_DE_utils.R"))
source(paste0(pipScriptDir, "/", "calc_smile_score.R"))
source(paste0(pipScriptDir, "/", "calc_fit_ks.R"))
source(paste0(pipScriptDir, "/", "plot_smile_fct.R"))
suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) # error bar
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) # error bar
suppressPackageStartupMessages(library(ggplot2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) # error bar
suppressPackageStartupMessages(library(ggpubr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) # error bar

registerDoMC(ifelse(SSHFS,2, nCpu)) # loaded from main_settings.R

# create the directories
curr_outFold <- paste0(pipOutFold, "/", script_name)
system(paste0("mkdir -p ", curr_outFold))

pipLogFile <- paste0(pipOutFold, "/", script_name, "_logFile.txt")
system(paste0("rm -f ", pipLogFile))

processDT <- read.delim(process_file,header=T, sep="\t", stringsAsFactors=F, col.names=c("dataset", "process"))
processDT_withSamp <- read.delim(process_file_with_samp,header=T, sep="\t", stringsAsFactors=F)

# plotType <- "svg"
# myHeight <- 7
# myWidth <- 10

# plotType is pdf by default => multi-pages plot

nRowPlot <- 4
nColPlot <- 2

# nRowPlot <- 1
# nColPlot <- 2
sizeHeight <- 7
sizeWidth <- 12

myHeight <- sizeHeight*nRowPlot
myWidth <-  sizeWidth * nColPlot


# "permThresh" the quantile of permutations to take is loaded from main_settings.R
# retrieved from main_settings.R toPlotRanked
#*********************************************************
# build table with following columns:

# dataset_name | smile score | slope_obs | slope_meanPerm | mean_slopeAllPerm | slope_obs/slope_meanPerm | KS pval sigFit | KS D stat sigFit | KS pval cumsum | KS D stat cumsum

#*********************************************************

toPlotRanked <- c("ratioDown","rescWeighted", "rescWeightedQQ", "prodSignedRatio")

txt <- paste0(toupper(script_name), "> Build score table for following datasets: ", paste0(comp_folders, collapse = ", "), "\n")
printAndLog(txt, pipLogFile)

txt <- paste0(toupper(script_name), "> Output ranked plots for: ", paste0(toPlotRanked, collapse = ", "), "\n")
printAndLog(txt, pipLogFile)

# for the ouputfile, sort according to the 1st of the list toPlotRanked
ref_to_rank <- toPlotRanked[1]

stopifnot(ref_to_rank %in% toPlotRanked)

########################################### INTERSECT FOR THE MEAN TAD CORR

all_scores_DT <- read.delim(paste0(scoreF), header=T, stringsAsFactors = F)
rownames(all_scores_DT) <- all_scores_DT$dataset

# LOADED FROM MAIN_SETTINGS.R
#nRandomPermut
#useTADonly
gene2tadDT <- read.delim(gene2tadDT_file, header=F, col.names = c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = F)
gene2tadDT$entrezID <- as.character(gene2tadDT$entrezID)

# 
for(curr_ratio_type in toPlotRanked) {
  rank_DT <- all_scores_DT[order(all_scores_DT[,paste0(curr_ratio_type, "_auc_rank")]),]
  outFile <- file.path(curr_outFold, paste0("smilePlot_vs_cumsum_", curr_ratio_type, ".pdf"))
  pdf(outFile, width = myWidth, height = myHeight)
  par(mfrow = c(nRowPlot, nColPlot))
  for(curr_ds in as.character(rownames(rank_DT))) {
    curr_ds_full <- comp_folders[grep(curr_ds, comp_folders)]
    stopifnot(length(curr_ds_full) == 1)

    cat(paste0("... plotting: ", curr_ratio_type, " - ", curr_ds,"\n" ))

    mytit <- paste0(curr_ds, " - ", curr_ratio_type, " rank ", all_scores_DT[curr_ds,paste0(curr_ratio_type, "_auc_rank") ])
    # plot the smile on the left, and the cumsum on the right

    ### DATA FOR THE SMILE SCORE
    file1 <- paste0(curr_ds_full, "/", script8_name, "/all_obs_", curr_ratio_type, ".Rdata")
    if(!file.exists(file1)) {
      stop(paste0("ERROR: ",  file1 , " does not exist !\n"))
    }
    cat(paste0(file1, "\n"))
    file2 <- paste0(curr_ds_full, "/", script8_name, "/", curr_ratio_type, "_permDT.Rdata")
    if(!file.exists(file2)) {
      stop(paste0("ERROR: ",  file2, " does not exist !\n"))
    }
    cat(paste0(file2, "\n"))
    obs_down <- eval(parse(text = load(file1)))
    permutDT_down <- eval(parse(text = load(file2)))

    ### DATA FOR THE CUMSUM
    obs_curr_down <- obs_down
    permut_currDown <- permutDT_down
    # ensure I used the same set of TADs for the permutation and for the calculation
    # (NB: would also be possible to filter the obs_curr_down, but not the permut_currDown)
    stopifnot(all(names(obs_curr_down) %in% rownames(permut_currDown)))
    stopifnot(all(rownames(permut_currDown) %in% names(obs_curr_down)))
    interReg <- intersect(names(obs_curr_down),rownames(permut_currDown) )
    filter_obs_curr_down <- sort(obs_curr_down[interReg], decreasing = T)
    filter_permut_currDown_unsort <- permut_currDown[interReg,]
    stopifnot(length(filter_obs_curr_down) == nrow(filter_permut_currDown_unsort))
    filter_permut_currDown <- apply(filter_permut_currDown_unsort, 2, sort, decreasing=T)
    rownames(filter_permut_currDown) <- NULL
    stopifnot(length(filter_obs_curr_down) == nrow(filter_permut_currDown_unsort))


    if(curr_ratio_type == "ratioDown") {
      curr_filter_obs_curr_down <- abs(filter_obs_curr_down - 0.5) + 0.5
      curr_filter_permut_currDown <- abs(filter_permut_currDown - 0.5) + 0.5
      departFromValue <- 0.5
      curr_ratio_type_concord <- paste0(curr_ratio_type, "_concord")
    } else if(curr_ratio_type == "rescWeightedQQ" | curr_ratio_type == "rescWeighted" ) {
      curr_filter_obs_curr_down <- abs(filter_obs_curr_down - 0.5) + 0.5
      curr_filter_permut_currDown <- abs(filter_permut_currDown - 0.5) + 0.5
      departFromValue <- 0.5
      curr_ratio_type_concord <- paste0(curr_ratio_type, "_concord")
    } else if(curr_ratio_type == "prodSignedRatio") {
      curr_filter_obs_curr_down <- filter_obs_curr_down
      curr_filter_permut_currDown <- filter_permut_currDown
      departFromValue <- 0
      curr_ratio_type_concord <- paste0(curr_ratio_type)
    } else{
      stop("should not happen")
    }

    ### PLOT SMILE PLOT
    plot_smile(obsvect_data_down =obs_down,
               mytit = mytit,
                           permutDT_data_down=permutDT_down,
                           nRandom=nRandomPermut,
                           TADonly=useTADonly,
                           g2tDT = gene2tadDT,
                           histBreakStep = 0.1,
                           cexMain=2,
                           toPlot = "top",
                           obsCol = "bisque",
                           obsColText = "bisque2",
                           permutCol = "mediumseagreen",
                           signedStat = (regexpr("signed", tolower(curr_ratio_type)) > 0))

    ############################### PLOT CUMSUM POLYGON 1) DEPARTURE 0.5, RATIO
    plot_cumsumDiff05(observ_vect = curr_filter_obs_curr_down,
                      permut_DT = curr_filter_permut_currDown,
                      my_stat = curr_ratio_type,
                      cexMain=2, departureValue = departFromValue)
    auc_ratio_stat_dep05 <- calc_aucObsPerm_cumsum(observ_vect = curr_filter_obs_curr_down,
                                                   permut_DT=curr_filter_permut_currDown,
                                            thresh_perm =permThresh,
                                            doPlot=F,
                                            departureValue=departFromValue)

    mtext(paste0("cumsum depart. ", departFromValue, " AUC Obs/Perm = ", round(auc_ratio_stat_dep05,2)))



  } # END ITERATING OVER DATASETS
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
}# END ITERATING OVER RATIO TYPES

for(curr_ratio_type in toPlotRanked) {
  rank_DT <- all_scores_DT[order(all_scores_DT[,paste0(curr_ratio_type, "_auc_rank")]),]
  outFile <- file.path(curr_outFold, paste0("cumsumArea_vs_cumsumThresh_", curr_ratio_type, ".pdf"))
  pdf(outFile, width = myWidth, height = myHeight)
  par(mfrow = c(nRowPlot, nColPlot))
  for(curr_ds in as.character(rownames(rank_DT))) {
    curr_ds_full <- comp_folders[grep(curr_ds, comp_folders)]
    stopifnot(length(curr_ds_full) == 1)

    cat(paste0("... plotting: ", curr_ratio_type, " - ", curr_ds,"\n" ))

    mytit <- paste0(curr_ds, " - ", curr_ratio_type, " rank ", all_scores_DT[curr_ds,paste0(curr_ratio_type, "_auc_rank") ])
    # plot the smile on the left, and the cumsum on the right

    ### DATA FOR THE SMILE SCORE
    file1 <- paste0(curr_ds_full, "/", script8_name, "/all_obs_", curr_ratio_type, ".Rdata")
    if(!file.exists(file1)) {
      stop(paste0("ERROR: ",  file1 , " does not exist !\n"))
    }
    cat(paste0(file1, "\n"))
    file2 <- paste0(curr_ds_full, "/", script8_name, "/", curr_ratio_type, "_permDT.Rdata")
    if(!file.exists(file2)) {
      stop(paste0("ERROR: ",  file2, " does not exist !\n"))
    }
    cat(paste0(file2, "\n"))
    obs_down <- eval(parse(text = load(file1)))
    permutDT_down <- eval(parse(text = load(file2)))

    ### DATA FOR THE CUMSUM
    obs_curr_down <- obs_down
    permut_currDown <- permutDT_down
    # ensure I used the same set of TADs for the permutation and for the calculation
    # (NB: would also be possible to filter the obs_curr_down, but not the permut_currDown)
    stopifnot(all(names(obs_curr_down) %in% rownames(permut_currDown)))
    stopifnot(all(rownames(permut_currDown) %in% names(obs_curr_down)))
    interReg <- intersect(names(obs_curr_down),rownames(permut_currDown) )
    filter_obs_curr_down <- sort(obs_curr_down[interReg], decreasing = T)
    filter_permut_currDown_unsort <- permut_currDown[interReg,]
    stopifnot(length(filter_obs_curr_down) == nrow(filter_permut_currDown_unsort))
    filter_permut_currDown <- apply(filter_permut_currDown_unsort, 2, sort, decreasing=T)
    rownames(filter_permut_currDown) <- NULL
    stopifnot(length(filter_obs_curr_down) == nrow(filter_permut_currDown_unsort))


    if(curr_ratio_type == "ratioDown") {
      curr_filter_obs_curr_down <- abs(filter_obs_curr_down - 0.5) + 0.5
      curr_filter_permut_currDown <- abs(filter_permut_currDown - 0.5) + 0.5
      departFromValue <- 0.5
      curr_ratio_type_concord <- paste0(curr_ratio_type, "_concord")
    } else if(curr_ratio_type == "rescWeightedQQ" | curr_ratio_type == "rescWeighted" ) {
      curr_filter_obs_curr_down <- abs(filter_obs_curr_down - 0.5) + 0.5
      curr_filter_permut_currDown <- abs(filter_permut_currDown - 0.5) + 0.5
      departFromValue <- 0.5
      curr_ratio_type_concord <- paste0(curr_ratio_type, "_concord")
    } else if(curr_ratio_type == "prodSignedRatio") {
      curr_filter_obs_curr_down <- filter_obs_curr_down
      curr_filter_permut_currDown <- filter_permut_currDown
      departFromValue <- 0
      curr_ratio_type_concord <- paste0(curr_ratio_type)
    }

    ############################### PLOT CUMSUM POLYGON
    plot_cumsumDiff05(observ_vect = curr_filter_obs_curr_down,
                      permut_DT = curr_filter_permut_currDown,
                      my_stat = curr_ratio_type_concord,
                      cexMain=1, departureValue = departFromValue)

    ############################### PLOT CUMSUM LINE PERMUT THRESHOLD
    auc_ratio_stat_dep05 <- calc_aucObsPerm_cumsum(observ_vect = curr_filter_obs_curr_down,
                                                   permut_DT=curr_filter_permut_currDown,
                                                   thresh_perm =permThresh,
                                                   my_stat = curr_ratio_type_concord,
                                                   doPlot=T,
                                                   departureValue=departFromValue)

    mtext(paste0("cumsum depart. ", departFromValue, " AUC Obs/Perm = ", round(auc_ratio_stat_dep05,2)))

  } # END ITERATING OVER DATASETS
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
}# END ITERATING OVER RATIO TYPES


cat("*** DONE\n")      
cat(paste0(startTime, "\n", Sys.time(), "\n"))


# 
# all_scores_DT <- read.delim(paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline/OUTPUT_FOLDER/SCORE_DT/18_plotOrderedParam/rank_scores_to_plot_DT.txt"))
# plotTable<-ggtexttable(all_scores_DT, 
#             rows = NULL,
#             theme = ttheme("mCyanWhite"))
# ggsave(plotTable, file = "tmp2.svg", height = 20, width=49)
# ggsave(plotTable, file = "tmp2.png", height = 800, width=2000)


