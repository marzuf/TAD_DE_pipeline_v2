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
script8_name <- "8c_runAllDown"
script9_name <- "9_runEmpPvalMeanTADLogFC"
script10_name <- "10_runEmpPvalMeanTADCorr"
script11_name <- "11_runEmpPvalCombined"
script_name <- "17e_build_score_table_limited"
stopifnot(file.exists(paste0(pipScriptDir, "/", script_name, ".R")))
cat(paste0("> START ", script_name,  "\n"))

# settingF <- "BUILDDT_SETTING_FILES/run_settings_buildDT.R"

source("main_settings.R")
#source("run_settings.R")
source(settingF)
source(paste0(pipScriptDir, "/", "TAD_DE_utils.R"))
source(paste0(pipScriptDir, "/", "calc_smile_score.R"))
source(paste0(pipScriptDir, "/", "calc_fit_ks.R"))
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

stopifnot(all(processDT$dataset == processDT_withSamp$dataset))

plotType <- "svg"

ref_to_rank <- "ratioDown"

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

# allDown <- c("ratioDown","rescWeighted", "rescWeightedQQ", "prodSignedRatio")
# -> loaded from main_settings.R

all_scores_DT <- foreach(curr_ds = comp_folders, .combine = 'rbind') %dopar% {
  
  cat(paste0("**** START computing scores for: ", basename(curr_ds), "\n"))
  

    # 2) then iterate over the ratio type to get the AUC ratio
    allratio_auc <- foreach(curr_ratio_type = allDown, .combine = 'c') %do% {
      cat(paste0("*** START ", curr_ratio_type, "\n"))
      cat(paste0("... load: ", paste0(curr_ds, "/", script8_name, "/all_obs_", curr_ratio_type, ".Rdata"), "\n"))
      cat(paste0("... load: ", paste0(curr_ds, "/", script8_name, "/", curr_ratio_type, "_permDT.Rdata"), "\n"))
      obs_curr_down <- eval(parse(text = load(paste0(curr_ds, "/", script8_name, "/all_obs_", curr_ratio_type, ".Rdata"))))
      permut_currDown <- eval(parse(text = load(paste0(curr_ds, "/", script8_name, "/", curr_ratio_type, "_permDT.Rdata"))))
      # ensure I used the same set of TADs for the permutation and for the calculation
      # (NB: would also be possible to filter the obs_curr_down, but not the permut_currDown)
      stopifnot(all(names(obs_curr_down) %in% rownames(permut_currDown)))
      stopifnot(all(rownames(permut_currDown) %in% names(obs_curr_down)))
      interReg <- intersect(names(obs_curr_down),rownames(permut_currDown) )
      ############################################################# prepare the vectors and tables
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
      } else if(curr_ratio_type == "rescWeightedQQ" | curr_ratio_type == "rescWeighted") {
        curr_filter_obs_curr_down <- abs(filter_obs_curr_down - 0.5) + 0.5
        curr_filter_permut_currDown <- abs(filter_permut_currDown - 0.5) + 0.5
        departFromValue <- 0.5
      } else if(curr_ratio_type == "prodSignedRatio") {
        curr_filter_obs_curr_down <- filter_obs_curr_down
        curr_filter_permut_currDown <- filter_permut_currDown
        departFromValue <- 0
      } else {
        stop(paste0("curr_ratio_type: ", curr_ratio_type, "\n", "should not happen"))
      }

      ############################### PLOT AUC CUMSUM POLYGON 1) DEPARTURE 0.5, RATIO
      auc_ratio_stat_dep05 <- calc_aucObsPerm_cumsum(observ_vect = curr_filter_obs_curr_down, 
                                                     permut_DT=curr_filter_permut_currDown, 
                                                     thresh_perm =permThresh, 
                                                     doPlot=F,
                                                     departureValue=departFromValue)
      auc_ratio_stat_dep05
    }
    ds_vect <- c(basename(curr_ds),  allratio_auc)
    names(ds_vect) <- c("dataset", paste0(allDown, "_auc"))
    ds_vect  
}
cat(paste0("... DT built, assemble and rank. Ranking will be down according to: ", ref_to_rank, "\n"))

all_scores_DT_s <- all_scores_DT
all_scores_DT <- data.frame(all_scores_DT, stringsAsFactors = F)
rownames(all_scores_DT) <- all_scores_DT$dataset
all_scores_DT$dataset <- NULL
# ensure this is not a characer
stopifnot(is.character(all_scores_DT[1,1]))
all_scores_DT <- data.matrix(all_scores_DT)
# order by decreasing since ratio is Obs/Perm
all_scores_rank_DT <- apply(-all_scores_DT,2,rank, ties="min")
colnames(all_scores_rank_DT) <- paste0(colnames(all_scores_rank_DT), "_rank")
all_DT <- cbind(all_scores_DT, all_scores_rank_DT)
all_DT <- all_DT[,sort(colnames(all_DT))]      
      
stopifnot(paste0(ref_to_rank, "_auc", "_rank") %in% colnames(all_DT))
all_DT <- all_DT[order(all_DT[,paste0(ref_to_rank, "_auc", "_rank")]),]

stopifnot(is.numeric(all_DT[1,1]))
all_scores_DT <- round(all_DT, 2)
all_scores_DT <- as.data.frame(all_scores_DT)
tmpCol <- colnames(all_scores_DT)
all_scores_DT$dataset <- rownames(all_scores_DT)

all_scores_DT$process <- unlist(sapply(as.character(all_scores_DT$dataset), function(x) {
                                process <- processDT$process[processDT$dataset == x ]
                                process <- gsub("_", " ", process)
                                if(length(process) == 0) process <- "undef"
                                return(process)
                            }))

tmpCol <- tmpCol[-which(tmpCol == paste0(ref_to_rank, "_auc", "_rank"))]
tmpCol <- tmpCol[-which(tmpCol == paste0(ref_to_rank, "_auc"))]

all_scores_DT <- all_scores_DT[,c("dataset", "process", paste0(ref_to_rank, "_auc"), paste0(ref_to_rank, "_auc_rank"),  tmpCol)]
dtFile <- paste0(curr_outFold, "/rank_all_scores_DT.txt")
write.table(all_scores_DT, file = dtFile, row.names=F, col.names=T, quote=F, sep="\t")
cat(paste0("... written: ", dtFile, "\n"))

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

all_scores_samples_DT <- all_scores_samples_DT[,c("dataset", "process", "cond1", "cond2", "nSamp1", "nSamp2", paste0(ref_to_rank, "_auc"), paste0(ref_to_rank, "_auc_rank"),  tmpCol)]
dtFile <- paste0(curr_outFold, "/rank_all_scores_DT_withSamp.txt")
write.table(all_scores_samples_DT, file = dtFile, row.names=F, col.names=T, quote=F, sep="\t")
cat(paste0("... written: ", dtFile, "\n"))

# stopifnot(paste0(ref_to_plot, "_rank") %in% colnames(all_DT))
# all_DT <- all_DT[order(paste0(ref_to_plot, "_rank")),]

# par(mfrow=c(6,2))
# for(i in 1:nrow(all_DT)) {
#   mainFold <- comp_folders[grep(rownames(all_DT)[i], comp_folders)]
#   stopifnot(length(mainFold == 1))
#   
#   mytit <- paste0(ref_to_plot, " - rank ", all_DT[i,paste0(ref_to_plot, "_rank")], ": ", rownames(all_DT)[i] )
#   
#   # plot the smile score
#   #outFile2_A <- paste0(curr_outFold, "/", curr_ratio_type, "departure05_cumsum_obs_permut_ratio0_1.", plotType)
# 
#   plot_cumsumDiff05(observ_vect = filter_obs_curr_down, permut_DT = filter_permut_currDown,  
#                     my_stat = curr_ratio_type, departureValue = 0.5)
#   mtext(subtitDir, font=3)
#   foo <- dev.off()
#   cat(paste0("... written: ", outFile2_A, "\n"))
#   
#   # and plot the cumsum
# }
#  
plotFile <- file.path(curr_outFold, "rank_all_scores_DT.svg")
plotTable <- ggtexttable(all_scores_DT, 
                        rows = NULL,
                        theme = ttheme("mCyanWhite"))
# title <- textGrob("Scores summary table", #y=unit(0.5,"npc") + 0.5*h, 
#                   vjust=0, gp=gpar(fontsize=20))
# grid.arrange(title, plotTable, ncol = 1, heights=c(0.1,1))
ggsave(plotTable, file = plotFile, height = 20, width=49)
cat(paste0("... written: ", plotFile, "\n"))                        


plotFile <- file.path(curr_outFold, "rank_all_scores_DT_withSamp.svg")
plotTable <- ggtexttable(all_scores_samples_DT, 
                         rows = NULL,
                         theme = ttheme("mCyanWhite"))
# title <- textGrob("Scores summary table", #y=unit(0.5,"npc") + 0.5*h, 
#                   vjust=0, gp=gpar(fontsize=20))
# grid.arrange(title, plotTable, ncol = 1, heights=c(0.1,1))
ggsave(plotTable, file = plotFile, height = 20, width=49)
cat(paste0("... written: ", plotFile, "\n"))                        



cat("*** DONE\n")      
cat(paste0(startTime, "\n", Sys.time(), "\n"))



# all_scores_DT <- read.delim(paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline/OUTPUT_FOLDER/SCORE_DT/17c_build_score_table_short/rank_all_scores_DT.txt"))
# plotTable<-ggtexttable(all_scores_DT, 
#             rows = NULL,
#             theme = ttheme("mCyanWhite"))
# ggsave(plotTable, file = "tmp.svg", height = 20, width=49)
# ggsave(plotTable, file = "tmp.png", height = 800, width=2000)
