#!/usr/bin/Rscript

options(scipen=100)

startTime <- Sys.time()

################  USE THE FOLLOWING FILES FROM PREVIOUS STEPS

# - script0: rna_geneList.Rdata
# - script0: pipeline_geneList.Rdata

# - script1: DE_topTable.Rdata
# - script1: DE_geneList.Rdata
# - script1: DE_rnaseqDT.Rdata

# - script8: all_obs_ratioDown.Rdata
# - script9: emp_pval_meanLogFC.Rdata
# - script10: emp_pval_meanCorr.Rdata
################################################################################

################  OUTPUT
# - emp_pval_combined.Rdata + tables
################################################################################

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
script_name <- "14f5_cumulAllDown_limited_AUC_flexQuant"
stopifnot(file.exists(paste0(pipScriptDir, "/", script_name, ".R")))
cat(paste0("> START ", script_name,  "\n"))

source("main_settings.R")
source(settingF)
source(paste0(pipScriptDir, "/", "TAD_DE_utils.R"))
suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) # error bar
suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) # error bar
suppressPackageStartupMessages(library(flux, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) # error bar

registerDoMC(ifelse(SSHFS,2, nCpu)) # loaded from main_settings.R

# create the directories
curr_outFold <- paste0(pipOutFold, "/", script_name)
system(paste0("mkdir -p ", curr_outFold))

pipLogFile <- paste0(pipOutFold, "/", format(Sys.time(), "%Y%d%m%H%M%S"),"_", script_name, "_logFile.txt")
system(paste0("rm -f ", pipLogFile))

# ADDED 27.11.2018 to check using other files
txt <- paste0("inputDataType\t=\t", inputDataType, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0("gene2tadDT_file\t=\t", gene2tadDT_file, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0("TADpos_file\t=\t", TADpos_file, "\n")
printAndLog(txt, pipLogFile)
txt <- paste0("settingF\t=\t", settingF, "\n")
printAndLog(txt, pipLogFile)

all_aucQuantThresh <- c(0.95)

# 14f:   leg_xpos <- ifelse(observed_auc <= max(density_permut$x), observed_auc,  max(density_permut$x))
# 14f2:  leg_xpos <- observed_auc

################################****************************************************************************************
####################################################### PREPARE INPUT
################################****************************************************************************************

# INPUT DATA
gene2tadDT <- read.delim(gene2tadDT_file, header=F, col.names = c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = F)
gene2tadDT$entrezID <- as.character(gene2tadDT$entrezID)

# UPDATE SELECT THE GENES ACCORDING TO THE SETTINGS PREPARED IN 0_PREPGENEDATA
initList <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/rna_geneList.Rdata"))))
geneList <- eval(parse(text = load(paste0(pipOutFold, "/", script0_name, "/pipeline_geneList.Rdata"))))
txt <- paste0(toupper(script_name), "> Start with # genes: ", length(geneList), "/", length(initList), "\n")
printAndLog(txt, pipLogFile)

stopifnot(!any(duplicated(names(geneList))))

gene2tadDT <- gene2tadDT[gene2tadDT$entrezID %in% geneList,]
geneNbr <- setNames(as.numeric(table(gene2tadDT$region)), names(table(gene2tadDT$region)))

### SET OUTPUT
plotType <- "svg"
myHeight <- ifelse(plotType == "png", 480 , 7)
myWidth <- ifelse(plotType == "png", 600, 10)

# if plotSeparated == TRUE => 1 plot per ratio, otherwise all on the same figure (#x2 plots)
plotSeparated <- F

# "permThresh" the quantile of permutations to take is loaded from main_settings.R

###############
##### retrieve the direction of up/down
###############
# retrieve the direction of up/down
DE_topTable <- eval(parse(text = load(paste0(pipOutFold, "/", script1_name, "/DE_topTable.Rdata"))))
DE_geneList <- eval(parse(text = load(paste0(pipOutFold, "/", script1_name, "/DE_geneList.Rdata"))))
exprDT <- eval(parse(text = load(paste0(pipOutFold, "/", script1_name, "/DE_rnaseqDT.Rdata"))))
# samp1 <- eval(parse(text=load(paste0(setDir, "/", sample1_file))))
# samp2 <- eval(parse(text=load(paste0(setDir, "/", sample2_file))))
samp1 <- eval(parse(text=load(paste0(sample1_file))))
samp2 <- eval(parse(text=load(paste0(sample2_file))))
DE_topTable <- DE_topTable[DE_topTable$genes %in% names(DE_geneList),]
stopifnot(all(DE_topTable$genes %in% names(DE_geneList)))
stopifnot(!any(duplicated(names(DE_geneList))))
stopifnot(all(colnames(exprDT) %in% c(samp1, samp2)))
stopifnot(all(samp1 %in% colnames(exprDT)))
stopifnot(all(samp2 %in% colnames(exprDT)))
maxDownGene <- DE_topTable$genes[which.min(DE_topTable$logFC)]
stopifnot(maxDownGene %in% rownames(exprDT))
mean_expr1 <- mean(unlist(c(exprDT[maxDownGene, samp1])), na.rm=T)
mean_expr2 <- mean(unlist(c(exprDT[maxDownGene, samp2])), na.rm=T)

if(mean_expr1 > mean_expr2) {
  subtitDir <- paste0("down: ", toupper(cond1), " > ", toupper(cond2))
} else{
  subtitDir <- paste0("down: ", toupper(cond2), " > ", toupper(cond1))
}

if(! plotSeparated) {
  nColPlot <- 2
																												#  nRowPlot <- length(allDown_limited)*2/nColPlot
  nRowPlot <- length(allDown_limited)*1/nColPlot
  outFile <- paste0(curr_outFold, "/allRatios_cumsum_obs_permut.", plotType)
  do.call(plotType, list(outFile, height=myHeight*nRowPlot, width=myWidth*nColPlot))
  par(mfrow=c(nRowPlot, nColPlot))
}

################################****************************************************************************************
####################################################### ITERATE OVER RATIOS TO PLOT
################################****************************************************************************************

all_ratios_all_auc_values <- list()

for(curr_ratio_type in allDown_limited) {
    cat(paste0("*** START ", curr_ratio_type, "\n"))
    
    obs_curr_down <- eval(parse(text = load(paste0(pipOutFold, "/", script8_name, "/all_obs_", curr_ratio_type, ".Rdata"))))
    permut_currDown <- eval(parse(text = load(paste0(pipOutFold, "/", script8_name, "/", curr_ratio_type, "_permDT.Rdata"))))
  
    # ensure I used the same set of TADs for the permutation and for the calculation
    # (NB: would also be possible to filter the obs_curr_down, but not the permut_currDown)
    stopifnot(all(names(obs_curr_down) %in% rownames(permut_currDown)))
    stopifnot(all(rownames(permut_currDown) %in% names(obs_curr_down)))

    interReg <- intersect(names(obs_curr_down),rownames(permut_currDown) )

    ############################################################
    ############################################################ # filter the TADs and sort
    ############################################################
    filter_obs_curr_down <- sort(obs_curr_down[interReg], decreasing = T)

    filter_permut_currDown_unsort <- permut_currDown[interReg,]
    stopifnot(length(filter_obs_curr_down) == nrow(filter_permut_currDown_unsort))

    filter_permut_currDown <- apply(filter_permut_currDown_unsort, 2, sort, decreasing=T)
    rownames(filter_permut_currDown) <- NULL
    stopifnot(length(filter_obs_curr_down) == nrow(filter_permut_currDown_unsort))

    # FOR ratioDown => plot ratioConcord, departure from 0.5
    if(curr_ratio_type == "ratioDown") {
      my_stat_curr_ratio <- "ratioDown_Concord"
      departureFromValue <- 0.5
      # => Concord, departure 0.5
      # Transform ratioDown -> ratioConcord
      # change so that the ratioDown ranges between 0.5 and 1 (-> e.g. treats 0.1 as 0.9)
      # transf. ratioDown -> ratioConcord
      filter_obs_curr_down_half <- abs(filter_obs_curr_down - 0.5) + 0.5
      filter_permut_currDown_half <- abs(filter_permut_currDown - 0.5) + 0.5
    } else if(curr_ratio_type == "rescWeightedQQ" | curr_ratio_type == "rescWeighted" ) {
        my_stat_curr_ratio <- paste0(curr_ratio_type, "_Concord")
        departureFromValue <- 0.5
        # => Concord, departure 0.5
        # Transform rescWeightedQQ -> rescWeightedQQConcord
        # change so that the ratioDown ranges between 0.5 and 1 (-> e.g. treats 0.1 as 0.9)
        # transf. ratioDown -> ratioConcord
        filter_obs_curr_down_half <- abs(filter_obs_curr_down - 0.5) + 0.5
        filter_permut_currDown_half <- abs(filter_permut_currDown - 0.5) + 0.5
    } else if(curr_ratio_type == "prodSignedRatio") {
      my_stat_curr_ratio <- "prodSignedRatio"
      departureFromValue <- 0
      # => raw (departure 0)
      # prodSignedRatio -> does not need to be transformed
      filter_obs_curr_down_half <- filter_obs_curr_down
      filter_permut_currDown_half <- filter_permut_currDown
    }  
    # PLOT THE 1ST PLOT
      if(plotSeparated) {
        outFile <- paste0(curr_outFold, "/", curr_ratio_type, "_departure05_cumsum_obs_permut.", plotType)
        do.call(plotType, list(outFile, height=myHeight, width=myWidth))
      }
      all_auc_values <- plot_cumsumDiff05_withLines(observ_vect = filter_obs_curr_down_half, 
                        permut_DT = filter_permut_currDown_half, 
						all_quantThresh = all_aucQuantThresh,
                        my_stat = my_stat_curr_ratio,
                        departureValue = departureFromValue, 
						drawline=TRUE)
      mtext(subtitDir, font=3)
      if(plotSeparated) {
        foo <- dev.off()
        cat(paste0("... written: ", outFile, "\n"))
      }
																									#      # PLOT THE 2ND PLOT
																									#      if(plotSeparated){
																									#        outFile <- paste0(curr_outFold, "/", curr_ratio_type, "_departure05_cumsum_obs_permut_AUC.", plotType)
																									#        do.call(plotType, list(outFile, height=myHeight, width=myWidth))
																									#      }
																									#      plot_cumsumDiff05_AUC2(filter_obs_curr_down_half, 
																									#                            filter_permut_currDown_half, 
																									#                            my_stat = my_stat_curr_ratio,
																									#                            departureValue = departureFromValue)
																									#      mtext(subtitDir, font=3)
																									#      if(plotSeparated){
																									#        foo <- dev.off()
																									#        cat(paste0("... written: ", outFile, "\n"))

  all_ratios_all_auc_values[[curr_ratio_type]] <- all_auc_values
																									#      }
}

if(!plotSeparated){
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
}


outFile <- file.path(curr_outFold, "all_ratios_all_auc_values.Rdata")
save(all_ratios_all_auc_values, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

txt <- paste0(startTime, "\n", Sys.time(), "\n")
printAndLog(txt, pipLogFile)

cat(paste0("*** DONE: ", script_name, "\n"))


