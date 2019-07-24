startTime <- Sys.time()

options(scipen=100)

### => in the folder 170v0: this was with the other version of computing the AUC ratio for the shuffTADs
### i.e. by taking the difference between the 2 curves without cutting at the end of the permutated curve (shorter)
### I change the way to calculate and store in 170v2 only the data with the shuffTADs ratios
### so now in 170_ => save all ratios (as in 170v0), but with the way of computing ratios for shuffTADs as in 170v2)

# here: I compute in this script only permGenes ratio, as in 170v2 

SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1)
settingF <- args[1]
stopifnot(file.exists(settingF))

pipScriptDir <- paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline_v2")

script0_name <- "0_prepGeneData"
script8_name <- "8c_runAllDown"
script_name <- "170revision2EZH2filter_score_auc_pval_permGenes"
stopifnot(file.exists(paste0(pipScriptDir, "/", script_name, ".R")))
cat(paste0("> START ", script_name,  "\n"))


histFC_file <- file.path(setDir, "/mnt/ed4/marie/scripts/EZH2_final_MAPQ/revision2_suppTable2.csv")
histFC_DT <- read.delim(histFC_file, sep=",", header=T, stringsAsFactors = FALSE)



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

cat("settingF = ", settingF, "\n")

# TADs_histFC <- setNames(histFC_DT$`Mean.H3K27me3.logFC`, histFC_DT$`TAD.ID`)

if(grepl("OUTPUT_FOLDER_EZH2_MAPQ_v2_CL", pipOutFold)) {
  column_id <- "Mean.H3K27me3.logFC"
} else {
  stop("error no hist FC\n")
}

TADs_histFC <- setNames(histFC_DT[, column_id], histFC_DT[, "TAD.ID"])

histFC_thresh <- 1
TADs_toKeep_histFC <- names(TADs_histFC[TADs_histFC <=  histFC_thresh])
stopifnot( TADs_histFC[names(TADs_histFC) %in% TADs_toKeep_histFC] <= histFC_thresh)

txt <- paste0(toupper(script_name), "> !!! HARD-CODED PATH TO histFC_file: ", histFC_file, "\n")
printAndLog(txt, pipLogFile)

txt <- paste0(toupper(script_name), "> !!! HARD-CODED COLUMN NAME WITH HIST. FC: ", column_id, "\n")
printAndLog(txt, pipLogFile)

txt <- paste0(toupper(script_name), "> Retain TADs with histFC <= ", histFC_thresh, "\n")
printAndLog(txt, pipLogFile)

txt <- paste0(toupper(script_name), "> # TADs to keep: ", length(TADs_toKeep_histFC), "\n")
printAndLog(txt, pipLogFile)


########################################### INTERSECT FOR THE MEAN TAD CORR


#all_scores_DT <- foreach(pipOutFold = comp_folders, .combine = 'rbind') %do% {
  
  cat(paste0("**** START computing scores for: ", basename(pipOutFold), "\n"))
  
  
  # 2) then iterate over the ratio type to get the AUC ratio
  allratio_auc_pval <- foreach(curr_ratio_type = allDown, .combine = 'c') %do% {
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
    
    
    stopifnot(TADs_toKeep_histFC %in% interReg)   ########################################### 30.08.2018 - added for filterTADs
    
    #****************************************************** CHANGED HERE FOR THE WAVE WITH ONLY SUBSET OF TADs
    txt <- paste0("... # number of TADs before filtering histFC:\t", length(interReg), "\n")
    printAndLog(txt, pipLogFile)
    
    interReg <- interReg[interReg %in% TADs_toKeep_histFC]
    txt <- paste0("... # number of TADs after filtering histFC:\t", length(interReg), "\n")
    printAndLog(txt, pipLogFile)
    
    
    outFile <- file.path(curr_outFold,"TADs_toKeep_histFC.Rdata")
    cat("... save the list of TADs used for the wave plot\n")
    save(TADs_toKeep_histFC, file = outFile)
    cat(paste0("... written: ", outFile, "\n"))
    
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
    
    
    
   auc_ratio_stat_dep05 <- calc_aucObsPerm_cumsum(observ_vect = curr_filter_obs_curr_down,
                                                  permut_DT=curr_filter_permut_currDown,
                                                  thresh_perm =permThresh,
                                                  doPlot=F,
                                                  departureValue=departFromValue)
   
   # use plot_cumsumDiff05_revision2 if permut and observed do not have same number of TADs
   

   auc_ratio_stat_dep05_pval <- plot_cumsumDiff05_AUC2(observ_vect = curr_filter_obs_curr_down,
                                                       permut_DT=curr_filter_permut_currDown,
                                                       toPlot = FALSE,
                                                       departureValue = departFromValue,
                                                       my_stat = curr_ratio_type)
    
                      # auc_ratio_stat_dep05_shuffTADs <- calc_aucObsPerm_cumsum_list_v2(observ_vect = curr_filter_obs_curr_down,
                      #                                                               shuff_List=curr_shuff_currDown,
                      #                                                               thresh_perm =permThresh,
                      #                                                               doPlot=F,
                      #                                                               departureValue=departFromValue)
                  #    
                  #    auc_ratio_stat_dep05_pval_shuffTADs <-   plot_cumsumDiff05_AUC2_list(observ_vect_cumsum=curr_filter_obs_curr_down,
                  #                                                                         permut_List_cumsum =curr_shuff_currDown,
                  #                                                                         toPlot = FALSE,
                  #                                                                         my_stat = curr_ratio_type,
                  #                                                                         departureValue = departFromValue)
 
              #    c(auc_ratio_stat_dep05, auc_ratio_stat_dep05_pval, auc_ratio_stat_dep05_shuffTADs, auc_ratio_stat_dep05_pval_shuffTADs)
                  # c(auc_ratio_stat_dep05_shuffTADs)
    c(auc_ratio_stat_dep05, auc_ratio_stat_dep05_pval)
    
  }
  # ds_vect <- c(basename(pipOutFold),  allratio_auc)
  # names(ds_vect) <- c("dataset", paste0(allDown, "_auc"))
  # names(ds_vect) <- c("dataset", paste0(rep(allDown,each=2), c("_auc", "_pval")))
  # names(ds_vect) <- c("dataset", paste0(rep(allDown,each=4), c("_auc_permGenes", "_pval_permGenes", "_auc_shuffTADs", "_pval_shuffTADs")))
#  names(allratio_auc_pval) <- paste0(rep(allDown,each=4), c("_auc_permGenes", "_pval_permGenes", "_auc_shuffTADs", "_pval_shuffTADs"))
                        # names(allratio_auc_pval) <- paste0(rep(allDown,each=1), c("_auc_shuffTADs"))
  names(allratio_auc_pval) <- paste0(rep(allDown,each=2), c("_auc_permGenes"))


outFile <-  paste0(curr_outFold, "/", "auc_ratios.Rdata")
save(allratio_auc_pval, file = outFile)
cat(paste0("... written: ", outFile, "\n"))




################################****************************************************************************************
####################################################### WRITE OUTPUT
################################****************************************************************************************


txt <- paste0(startTime, "\n", Sys.time(), "\n")
printAndLog(txt, pipLogFile)

cat(paste0("*** DONE: ", script_name, "\n"))






