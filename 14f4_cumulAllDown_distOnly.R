#!/usr/bin/Rscript

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
script_name <- "14f4_cumulAllDown_distOnly"
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

# 14f:   leg_xpos <- ifelse(observed_auc <= max(density_permut$x), observed_auc,  max(density_permut$x))
# 14f2:  leg_xpos <- observed_auc

linePermutCol <-   rgb(0/255,76/255,153/255)

lineDistCol <- "orange"

stopifnot(exists("permThresh"))

chgMar <- c(0,0,0,3)

initMar <- par()$mar


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
myWidth <- ifelse(plotType == "png", 800, 12)

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
  nColPlot <- 3
  nRowPlot <- length(allDown_limited)*2/nColPlot
  outFile <- paste0(curr_outFold, "/allRatios_dist_obs_permut.", plotType)
  do.call(plotType, list(outFile, height=myHeight*nRowPlot, width=myWidth*nColPlot))
  par(mfrow=c(nRowPlot, nColPlot))
}


################################****************************************************************************************
####################################################### ITERATE OVER RATIOS TO PLOT
################################****************************************************************************************

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
    
    obs_vect <- sort(filter_obs_curr_down_half, decreasing = T)
    permut_DT <- apply(filter_permut_currDown_half, 2, sort, decreasing=T)
    cumsum_permut_DT <- apply(permut_DT, 2, function(x) cumsum(abs(x-departureFromValue)))
    cumsum_permut_vect <- apply(cumsum_permut_DT, 1 ,function(x) as.numeric(quantile(x, probs=permThresh)))
    
    
    my_main <- paste0(curr_ratio_type, ": cumsum departure from ", departureFromValue)
    my_ylab <- paste0("cumsum(abs(", curr_ratio_type, " - ", departureFromValue,"))")
    my_xlab <- paste0("regions ranked by decreasing ", curr_ratio_type)
    
    ### PLOT 1) FULL CURVES
    full_curve_obs_vect <- cumsum(abs(obs_vect-departureFromValue)) 
    full_curve_perm_vect <- cumsum_permut_vect
    full_dist_vect <- full_curve_obs_vect - full_curve_perm_vect
    
    if(plotSeparated) {
      outFile <- paste0(curr_outFold, "/", curr_ratio_type, "_fullCurve_dist_obs_permut.", plotType)
      do.call(plotType, list(outFile, height=myHeight, width=myWidth))
    }
    par(mar = initMar + chgMar)
    plot_with_departure(vect_obs = full_curve_obs_vect, 
                        vect_perm = full_curve_perm_vect,
                        col_obs = "black",
                        col_perm = linePermutCol,
                        vect_axis4 = full_dist_vect,
                        col_axis4 = lineDistCol, 
                        lab_axis4 = paste0('Distance between the 2 curves'),
                        main = my_main, xlab = my_xlab, ylab = my_ylab)
    
    mod0 <- lm(full_dist_vect ~ c(1:length(full_dist_vect)))
    slope_fullDist <- as.numeric(coef(mod0)[2])
    abline(mod0, col="darkorange", lty=2)
    legTxt <- c("fitted lm", paste0("(slope=", round(slope_fullDist,4), ")"))
    legend("topleft", legTxt, lty=c(2,NA), col = c("darkorange", "black"), bty="n")
    
    if(plotSeparated){
      foo <- dev.off()
      cat(paste0("... written: ", outFile, "\n"))
    }
    par(mar = initMar)
    
    
    ### PLOT 2) CURVES 1st PART = FOR ONLY THE PART WITH THE 1 FOR THE OBSERVED RATIO
    idx_onlyOne <- which(sort(obs_vect, decreasing=T) == 1 )
    if(length(idx_onlyOne) > 0) {
      stopifnot(diff(idx_onlyOne) == 1)
      stopifnot(sort(obs_vect, decreasing=T)[idx_onlyOne] == 1)
      first_curve_obs_vect <- full_curve_obs_vect[idx_onlyOne] 
      first_curve_perm_vect <- full_curve_perm_vect[idx_onlyOne]
      first_dist_vect <- first_curve_obs_vect - first_curve_perm_vect
      if(plotSeparated) {
        outFile <- paste0(curr_outFold, "/", curr_ratio_type, "_1stPartCurve_dist_obs_permut.", plotType)
        do.call(plotType, list(outFile, height=myHeight, width=myWidth))
      }
      par(mar = initMar + chgMar)
      plot_with_departure(vect_obs = first_curve_obs_vect, 
                          vect_perm = first_curve_perm_vect,
                          col_obs = "black",
                          col_perm = linePermutCol,
                          vect_axis4 = first_dist_vect,
                          col_axis4 = lineDistCol, 
                          lab_axis4 = paste0('Distance between the 2 curves'),
                          main = my_main, xlab = my_xlab, ylab = my_ylab)
      
      mod1 <- lm(first_dist_vect ~ c(1:length(first_dist_vect)))
      slope_firstDist <- as.numeric(coef(mod1)[2])
      abline(mod1, col="darkorange", lty=2)
      legTxt <- c("fitted lm", paste0("(slope=", round(slope_firstDist,4), ")"))
      legend("topleft", legTxt, lty=c(2,NA), col = c("darkorange", "black"), bty="n")
      
      if(plotSeparated){
        foo <- dev.off()
        cat(paste0("... written: ", outFile, "\n"))
      }
      par(mar = initMar)
    } else{
      plot(NULL, xlim=c(0,1), ylim=c(0,1), axes = F, xlab="", ylab="")
    }
    
    ### PLOT 3) CURVES 2ND PART = FOR ONLY THE PART WITH THE 1 FOR THE OBSERVED RATIO
    idx_noOne <- which(sort(obs_vect, decreasing=T) < 1 )
    if(length(idx_noOne) > 0) {
      stopifnot(diff(idx_noOne) == 1)
      stopifnot(sort(obs_vect, decreasing=T)[1:(idx_noOne[1]-1)] == 1)
      scd_curve_obs_vect <- full_curve_obs_vect[idx_noOne] 
      scd_curve_perm_vect <- full_curve_perm_vect[idx_noOne]
      scd_dist_vect <- scd_curve_obs_vect - scd_curve_perm_vect

      if(plotSeparated) {
        outFile <- paste0(curr_outFold, "/", curr_ratio_type, "_fullCurve_dist_obs_permut.", plotType)
        do.call(plotType, list(outFile, height=myHeight, width=myWidth))
      }
      par(mar = initMar + chgMar)
      plot_with_departure(vect_obs = scd_curve_obs_vect, 
                          vect_perm = scd_curve_perm_vect,
                          col_obs = "black",
                          col_perm = linePermutCol,
                          vect_axis4 = scd_dist_vect,
                          col_axis4 = lineDistCol, 
                          lab_axis4 = paste0('Distance between the 2 curves'),
                          main = my_main, xlab = my_xlab, ylab = my_ylab)
      
      mod2 <- lm(scd_dist_vect ~ c(1:length(scd_dist_vect)))
      slope_scdDist <- as.numeric(coef(mod2)[2])
      abline(mod2, col="darkorange", lty=2)
      legTxt <- c("fitted lm", paste0("(slope=", round(slope_scdDist,4), ")"))
      legend("topleft", legTxt, lty=c(2,NA), col = c("darkorange", "black"), bty="n")
      
      if(plotSeparated){
        foo <- dev.off()
        cat(paste0("... written: ", outFile, "\n"))
      }
      par(mar = initMar)
    } else {
      plot(NULL, xlim=c(0,1), ylim=c(0,1), axes = F, xlab="", ylab="")
    }    
    
}

if(!plotSeparated){
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
}


slopes_obsPermDist <- c(
  slope_fullCurve = slope_fullDist,
  slope_firstPart = slope_firstDist,
  slope_scdPart = slope_scdDist
)

outFile <-  file.path(curr_outFold, "slopes_obsPermDist.Rdata")
save(slopes_obsPermDist, file = outFile)
cat(paste0("... written: ", outFile, "\n"))

txt <- paste0(startTime, "\n", Sys.time(), "\n")
printAndLog(txt, pipLogFile)

cat(paste0("*** DONE: ", script_name, "\n"))


