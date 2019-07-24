startTime <- Sys.time()

options(scipen=100)

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
script_name <- "170prunePermGenes_score_auc_pval_withShuffle"
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

# ADDED 27.11.2018 to check using other files
txt <- paste0("gene2tadDT_file\t=\t", gene2tadDT_file, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0("TADpos_file\t=\t", TADpos_file, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0("settingF\t=\t", settingF, "\n")
printAndLog(txt, pipLogFile)

# "permThresh" the quantile of permutations to take is loaded from main_settings.R

# allDown retrieved from main_settings.R
# toPlotRanked retrieved from run_settings_<>.R

myHeight <- myWidth <- 7

i_increment <- 10

########################################### INTERSECT FOR THE MEAN TAD CORR


#all_scores_DT <- foreach(pipOutFold = comp_folders, .combine = 'rbind') %do% {
  
  cat(paste0("**** START computing scores for: ", basename(pipOutFold), "\n"))
  
  
  # 2) then iterate over the ratio type to get the AUC ratio
  allratio_prune_thresh_permGenes <- foreach(curr_ratio_type = allDown, .combine = 'c') %do% {
    cat(paste0("*** START ", curr_ratio_type, "\n"))
    
    cat(paste0("... load: ", paste0(pipOutFold, "/", script8_name, "/all_obs_", curr_ratio_type, ".Rdata"), "\n"))
    obs_curr_down <- eval(parse(text = load(paste0(pipOutFold, "/", script8_name, "/all_obs_", curr_ratio_type, ".Rdata"))))
    
    cat(paste0("... load: ", paste0(pipOutFold, "/", script8_name, "/", curr_ratio_type, "_permDT.Rdata"), "\n"))
    permut_currDown <- eval(parse(text = load(paste0(pipOutFold, "/", script8_name, "/", curr_ratio_type, "_permDT.Rdata"))))
    
    # cat(paste0("... load: ", paste0(pipOutFold, "/", script8_name, "/", curr_ratio_type, "_randomShuffleList.Rdata"), "\n"))
    # shuff_currDown <-  eval(parse(text = load(paste0(pipOutFold, "/", script8_name, "/", curr_ratio_type, "_randomShuffleList.Rdata"))))
    
    # ensure I used the same set of TADs for the permutation and for the calculation
    # (NB: would also be possible to filter the obs_curr_down, but not the permut_currDown)
    stopifnot(all(names(obs_curr_down) %in% rownames(permut_currDown)))
    stopifnot(all(rownames(permut_currDown) %in% names(obs_curr_down)))
    interReg <- intersect(names(obs_curr_down),rownames(permut_currDown) )
    ############################################################# prepare the vectors and tables
    cat("... prepare observed data\n")
    filter_obs_curr_down <- sort(obs_curr_down[interReg], decreasing = T)
    
    cat("... prepare permut. data\n")
    filter_permut_currDown_unsort <- permut_currDown[interReg,]
    stopifnot(length(filter_obs_curr_down) == nrow(filter_permut_currDown_unsort))
    filter_permut_currDown <- apply(filter_permut_currDown_unsort, 2, sort, decreasing=T)
    rownames(filter_permut_currDown) <- NULL
    stopifnot(length(filter_obs_curr_down) == nrow(filter_permut_currDown_unsort))
    
    # cat("... prepare shuff. data\n")
    # shuff_currDown <- lapply(shuff_currDown, sort, decreasing=T)
    
    # return(NULL)
    
    if(curr_ratio_type == "ratioDown") {
      curr_filter_obs_curr_down <- abs(filter_obs_curr_down - 0.5) + 0.5
      curr_filter_permut_currDown <- abs(filter_permut_currDown - 0.5) + 0.5
      # curr_shuff_currDown <- lapply(shuff_currDown, function(x) abs(x - 0.5) + 0.5)
      departFromValue <- 0.5
    } else if(curr_ratio_type == "rescWeightedQQ" | curr_ratio_type == "rescWeighted") {
      curr_filter_obs_curr_down <- abs(filter_obs_curr_down - 0.5) + 0.5
      curr_filter_permut_currDown <- abs(filter_permut_currDown - 0.5) + 0.5
      # curr_shuff_currDown <- lapply(shuff_currDown, function(x) abs(x - 0.5) + 0.5)
      departFromValue <- 0.5
    } else if(curr_ratio_type == "prodSignedRatio") {
      curr_filter_obs_curr_down <- filter_obs_curr_down
      curr_filter_permut_currDown <- filter_permut_currDown
      # curr_shuff_currDown <- shuff_currDown
      departFromValue <- 0
    } else {
      stop(paste0("curr_ratio_type: ", curr_ratio_type, "\n", "should not happen"))
    }
    
    ############################### GET THE AUC AND THE PVAL - FOR GENE PERMUTATION DATA
    stopifnot(length(curr_filter_obs_curr_down) == nrow(curr_filter_permut_currDown))
    
    curr_filter_obs_curr_down <- sort(curr_filter_obs_curr_down, decreasing = T)
    curr_filter_permut_currDown <- apply(curr_filter_permut_currDown, 2, sort, decreasing=T)
    # curr_shuff_currDown <- lapply(curr_shuff_currDown, function(x) { sort(x, decreasing=TRUE) })
    
    i <- 1
    i_permGenes <- 0
    
    
    while( i_permGenes == 0 & i <= length(curr_filter_obs_curr_down) ) {
      cat(paste0("... ", curr_ratio_type, " - start pruning i=", i, "\n"))
      
      if(i_permGenes == 0) {
        toRowDT <- nrow(curr_filter_permut_currDown) - i + 1
        auc_ratio_stat_dep05 <- calc_aucObsPerm_cumsum(observ_vect = curr_filter_obs_curr_down[i:length(curr_filter_obs_curr_down)],
                                                       permut_DT=curr_filter_permut_currDown[1:toRowDT,],
                                                       # observ_vect = curr_filter_obs_curr_down[i:length(curr_filter_obs_curr_down)],
                                                       # permut_DT=curr_filter_permut_currDown[i:nrow(curr_filter_permut_currDown),],
                                                       # observ_vect = curr_filter_obs_curr_down,
                                                       # permut_DT=curr_filter_permut_currDown,
                                                       thresh_perm =permThresh,
                                                       doPlot=F,
                                                       departureValue=departFromValue)
        
        if(auc_ratio_stat_dep05 < 1) {
          i_permGenes <- i
          ratio_pruneThresh_permGenes <- curr_filter_obs_curr_down[i]
          auc_ratio_pruneThresh_permGenes <- auc_ratio_stat_dep05
          cat(paste0("----- ", curr_ratio_type, " - break permGenes: i=", i_permGenes, " (ratio=", round(ratio_pruneThresh_permGenes,2), ")\t found AUC=", auc_ratio_stat_dep05, "\n"))
        } else{
          cat("...... ", curr_ratio_type,  " - permGenes AUC = ", auc_ratio_stat_dep05, "\t -> go to next iteration\n")
        }
      }

            
      i <- i + i_increment
      # break
    } # end-while increasing i
    
    # i_permGenes <-  100
    
    outFile <- paste0(curr_outFold, "/", curr_ratio_type, "_cumsum_plot_init_and_prune.svg")
    svg(outFile, width = myWidth * 2.2, height = myHeight)
    par(mfrow=c(1,2))
    
    foo1 <- calc_aucObsPerm_cumsum(observ_vect = curr_filter_obs_curr_down,
                                                   permut_DT=curr_filter_permut_currDown,
                                                   thresh_perm =permThresh,
                                                   doPlot=TRUE,
                                                   departureValue=departFromValue)
    
    toRowDT <- nrow(curr_filter_permut_currDown) - i_permGenes + 1
    
    foo1 <- calc_aucObsPerm_cumsum(observ_vect = curr_filter_obs_curr_down[i_permGenes:length(curr_filter_obs_curr_down)],
                                                   #permut_DT=curr_filter_permut_currDown[i_permGenes:nrow(curr_filter_permut_currDown),],
                                   permut_DT=curr_filter_permut_currDown[1:toRowDT,],
                                                   thresh_perm =permThresh,
                                                   doPlot=TRUE,
                                                   departureValue=departFromValue)
    

    foo <- dev.off()
    cat(paste0("... written: ", outFile, "\n"))
    
    c(i_permGenes, ratio_pruneThresh_permGenes, auc_ratio_pruneThresh_permGenes)
      # i_shuffTADs, ratio_pruneThresh_shuffTADs, auc_ratio_pruneThresh_shuffTADs)
    
  } # end iterating over ratios
  
  names(allratio_prune_thresh_permGenes) <- paste0(rep(allDown,each=3), c("_idx_prune_permGenes", "_ratio_prune_permGenes", "_auc_prune_permGenes"))#, 
                                                                # "_idx_prune_shuffTADs", "_ratio_prune_shuffTADs", "_auc_prune_shuffTADs"))

increment_step <- c(increment_step=i_increment)
  
allratio_prune_thresh_permGenes <- c(increment_step,allratio_prune_thresh_permGenes)

outFile <-  paste0(curr_outFold, "/", "allratio_prune_thresh_permGenes.Rdata")
save(allratio_prune_thresh_permGenes, file = outFile)
cat(paste0("... written: ", outFile, "\n"))




################################****************************************************************************************
####################################################### WRITE OUTPUT
################################****************************************************************************************


txt <- paste0(startTime, "\n", Sys.time(), "\n")
printAndLog(txt, pipLogFile)

cat(paste0("*** DONE: ", script_name, "\n"))






