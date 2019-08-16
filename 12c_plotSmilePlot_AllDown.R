#!/usr/bin/Rscript

options(scipen=100)

startTime <- Sys.time()

################  USE THE FOLLOWING FILES FROM PREVIOUS STEPS
# - script0: rna_geneList.Rdata
# - script0: pipeline_geneList.Rdata
# - script8: <ratio>_permDT.Rdata
# - script8: all_obs_<ratio>.Rdata
################################################################################

################  OUTPUT
# - emp_pval_meanCorr.Rdata + plots
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
script_name <- "12c_plotSmilePlot_AllDown"
stopifnot(file.exists(paste0(pipScriptDir, "/", script_name, ".R")))
cat(paste0("> START ", script_name,  "\n"))

source("main_settings.R")
#source("run_settings.R")
source(settingF)
source(paste0(pipScriptDir, "/", "TAD_DE_utils.R"))
suppressPackageStartupMessages(library(Hmisc, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) # error bar


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

#allDown <- c("ratioDown", "FCdown", "prodConcord", "prodMeanConcord", "prodLogRatioNbr", "prodRatioSum") LOADED FROM main_settings.R


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

## TO SET
## hard-coded:
obsCol <- "bisque"
obsColText <- "bisque2"
permutCol <- "mediumseagreen"
plotType <- "svg"

myHeight <- ifelse(plotType == "png", 1028 , 15)
myWidth <- ifelse(plotType == "png", 686, 10)

#****************************************************** IMPORTANT HARD-CODED SETTINGS HERE !!!
histBreakStep <- 0.1
stopifnot( (1/histBreakStep) %% 1 == 0 )
#***********************************************************************************


################################****************************************************************************************
####################################################### ITERATE OVER RATIOS AND PLOT
################################****************************************************************************************


for(curr_ratio_type in allDown) {

  curr_ratio_breaks <- switch(1+as.numeric(regexpr("signed", tolower(curr_ratio_type)) > 0), seq(0,1,histBreakStep), seq(-1,1,histBreakStep))

  cat(paste0("*** START smile plot for: ", curr_ratio_type, "\n"))
  
  obs_curr_down <- eval(parse(text = load(paste0(pipOutFold, "/", script8_name, "/all_obs_", curr_ratio_type, ".Rdata"))))
  permut_currDown_DT <- eval(parse(text = load(paste0(pipOutFold, "/", script8_name, "/", curr_ratio_type, "_permDT.Rdata"))))
  
  ### !!! NEED TO CHECK WHY THERE ARE SOME NA
  permut_currDown_DT <- na.omit(permut_currDown_DT)
  
  nRandom <- nRandomPermut
  
  intersectRegions <- intersect(names(obs_curr_down), rownames(permut_currDown_DT))
  gene2tadDT <- gene2tadDT[gene2tadDT$region %in% intersectRegions,]
  # some regions could have been discarded because no gene
  intersectRegions <- intersect(intersectRegions, gene2tadDT$region)
  
  # take only the TADs
  ### filter TAD only regions
  if(useTADonly) {
      if( length(grep("_TAD", rownames(permut_currDown_DT))) > 0 ) {
          txt <- paste0(toupper(script_name), "> !!! WARNING: permutation ", curr_ratio_type, " data contain non-TAD regions as well !!!\n")
          printAndLog(txt, pipLogFile)    
      }
      initLen <- length(intersectRegions)
      intersectRegions <- intersectRegions[grep("_TAD", intersectRegions)]
      txt <- paste0(toupper(script_name), "> Take only TAD regions: ", length(intersectRegions), "/", initLen, "\n")
      printAndLog(txt, pipLogFile)    
  }
  
  gene2tadDT <- gene2tadDT[gene2tadDT$region %in% intersectRegions,]
  nbrGenes <- setNames(as.numeric(table(gene2tadDT$region)), names(table(gene2tadDT$region)))
  
  outFile <- paste0(curr_outFold, "/", curr_ratio_type, "smilePlotObservExpect_by", histBreakStep, "_filter.", plotType)
  outFile2 <- paste0(curr_outFold, "/", curr_ratio_type, "other_smilePlotObservExpect_by", histBreakStep, "_filter.", plotType)
  
  initNrow <- nrow(permut_currDown_DT)
  permut_currDown_DT <- permut_currDown_DT[intersectRegions,]
  txt <- paste0(toupper(script_name), "> Take only the regions form permut for which genes were retained: ", nrow(permut_currDown_DT), "/", initNrow, "\n")
  printAndLog(txt, pipLogFile)
  
  initLen <- length(obs_curr_down)
  obs_curr_down <- obs_curr_down[intersectRegions]
  txt <- paste0(toupper(script_name), "> Take only the regions form observ. for which genes were retained: ", length(obs_curr_down), "/", initLen, "\n")
  printAndLog(txt, pipLogFile)
  
  stopifnot(nrow(permut_currDown_DT) == length(obs_curr_down))
  
  ###################################################################################################
  
  cat("... Warning - considering only TAD regions: \n")
  
  cat("> Prepare simulation data \n")
  nTAD <- sum(regexpr("_TAD",  rownames(permut_currDown_DT)) > 0)
  cat(paste0("... TADs: ", nTAD, "/", nrow(permut_currDown_DT), "\n"))
  permut_currDown_vect <- as.vector(permut_currDown_DT[grep("_TAD", rownames(permut_currDown_DT)),])
  stopifnot(length(permut_currDown_vect) == ncol(permut_currDown_DT) * nTAD)
  
  cat(paste0("length(permut_currDown_vect) = ", length(permut_currDown_vect), "\n"))
  
  cat("> Prepare observed data \n")
  nTAD <- sum(regexpr("_TAD",  names(obs_curr_down)) > 0)
  cat(paste0("... TADs: ", nTAD, "/", length(obs_curr_down), "\n"))
  obs_curr_down <- obs_curr_down[grep("_TAD", names(obs_curr_down))]
  
  cat("> Prepare data for plotting \n")
  
  obs_rel_hist <- try(hist(obs_curr_down, breaks=curr_ratio_breaks, plot=F))
  #plot(obs_rel_hist)
  # some are not bounded 0-1 => go to next directly !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"WILL NEED REFINEMENT!!!
  if(class(obs_rel_hist) == "try-error") {
    warning(paste0("!!! CANNOT PLOT HISTOGRAM FOR: ", curr_ratio_type, " !!! \n" ))
    next
  }
  
  cat(paste0("length(permut_currDown_vect) = ", length(permut_currDown_vect), "\n"))
  cat(paste0("length(permut_currDown_vect) = ", length(na.omit(permut_currDown_vect)), "\n"))
  shuff_rel_hist <-  hist(permut_currDown_vect,breaks=curr_ratio_breaks, plot=F)
  # plot(shuff_rel_hist)
  
  cat(paste0("ncol(permut_currDown_DT) = ", ncol(permut_currDown_DT), "\n"))
  cat(paste0("sum(shuff_rel_hist$counts) = ", sum(shuff_rel_hist$counts), "\n"))
  cat(paste0("sum(shuff_rel_hist$counts/nRandom) = ", sum(shuff_rel_hist$counts/nRandom), "\n"))
  cat(paste0("sum(obs_rel_hist$counts) = ", sum(obs_rel_hist$counts),"\n"))
  stopifnot(sum(shuff_rel_hist$counts/nRandom) == sum(obs_rel_hist$counts))   # 830
  stopifnot(sum(shuff_rel_hist$counts) == length(permut_currDown_vect))
  stopifnot(sum(obs_rel_hist$counts) == length(na.omit(obs_curr_down)))
  
  rel_freqValues <- rbind(obs_rel_hist$counts, shuff_rel_hist$counts/nRandom)
  rownames(rel_freqValues) <- c("observed", "randomized")
  
  # FOR THE ERROR BARS:
  # for each randomization, calculate the observed, then divide expected by observed
  # table 1 column = 1 randomization, rows = breaks of the hist.
  rel_HistValues <- apply(permut_currDown_DT[grep("_TAD", rownames(permut_currDown_DT)),],
                          2, function(x) hist(x,breaks=curr_ratio_breaks, plot=F)$counts)
  stopifnot(nrow(rel_HistValues) == length(obs_rel_hist$counts))
  # for each randomization, divide observ./exp.
  rel_obsOverExp <- apply(rel_HistValues, 2, function(x) log2(obs_rel_hist$counts/x))
  # now get the SEM for each break of the hist (i.e. each row)
  rel_obsOverExp_sem <- apply(rel_obsOverExp, 1, function(x) {
    x <- na.omit(x[is.finite(x)])
    sd(x, na.rm=T)/sqrt(length(x))
  })
  
  # Calculate observed/expected ratio
  rel_logRatio <- log2(rel_freqValues["observed",]/rel_freqValues["randomized",])
  # put y coord at the middle of the breaks
  rel_ycoord <- (curr_ratio_breaks * 100)[-1]
  rel_ycoord <- rel_ycoord-(histBreakStep*100)/2
  stopifnot(length(rel_ycoord) == length(rel_logRatio))
  toKeep <- !is.na(rel_logRatio) & abs(rel_logRatio) != "Inf"
  rel_logRatio <- rel_logRatio[toKeep]
  rel_ycoord <- rel_ycoord[toKeep]
  rel_obsOverExp_sem <- rel_obsOverExp_sem[toKeep] 
  
  smile_score <- round(get_smileScore(ylogRatioVect=rel_logRatio, xBreakVect=rel_ycoord, withPlot=FALSE), 2)
  
  my_ylab <- paste0(rel_ycoord - (histBreakStep*100)/2, "-", rel_ycoord + (histBreakStep*100)/2, "%")
  
  # retrieve the direction of up/down
  DE_topTable <- eval(parse(text = load(paste0(pipOutFold, "/", script1_name, "/DE_topTable.Rdata"))))
  DE_geneList <- eval(parse(text = load(paste0(pipOutFold, "/", script1_name, "/DE_geneList.Rdata"))))
  exprDT <- eval(parse(text = load(paste0(pipOutFold, "/", script1_name, "/DE_rnaseqDT.Rdata"))))
  samp1 <- eval(parse(text=load(paste0(setDir, "/", sample1_file))))
  samp2 <- eval(parse(text=load(paste0(setDir, "/", sample2_file))))
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
  
  # plot the % only if not signed
  if(regexpr("signed", tolower(curr_ratio_type)) > 0) {
    rel_ycoord <- rel_ycoord/100
  }
  
  
  ####### SMILE PLOT HERE
  do.call(plotType, list(outFile, height=myHeight, width=myWidth))
  par(mfrow = c(2,1))
  # top panel: line plot
  myxlab <- ifelse(regexpr("signed", tolower(curr_ratio_type)) > 0, curr_ratio_type, paste0("% of ", curr_ratio_type, " per TAD"))
  plot(rel_logRatio ~ rel_ycoord, type="l",axes=F, cex.main=0.9,
       main = "Smile plot",
       xlab = myxlab, ylab = "log2(Observed/Randomized)")
  mtext(paste0("smile score = ", smile_score, "; ", subtitDir))
  box(bty="l")
  axis(2)

    if(regexpr("signed", tolower(curr_ratio_type)) > 0) {
      axis(1)
    } else{
      axis(1, at=rel_ycoord, labels = my_ylab)
    }

  abline(h=0, lty=2)
  errbar(x=rel_ycoord, y=rel_logRatio, yplus=rel_logRatio+rel_obsOverExp_sem, yminus=rel_logRatio-rel_obsOverExp_sem, add=T)
  # bottom panel: bar plot
  barplot(rel_freqValues, beside=T, col=c(obsCol, permutCol), ylab="Frequency", 
          xlab=myxlab)
  legend("topright", c( "observed Freq.", "randomized Freq."),
         text.col=c(obsColText, permutCol), bty='n', cex=1)
  legend("topleft", legend = subtitDir,
          bty='n', cex=1)
  # xlabpos <- seq(2, 30, by=3)
  rel_ycoord <- (curr_ratio_breaks * 100)[-1]
  rel_ycoord <- rel_ycoord-(histBreakStep*100)/2
  
  # plot the % only if not signed
  if(regexpr("signed", tolower(curr_ratio_type)) > 0) {
    rel_ycoord <- rel_ycoord/100
  }
  
  my_ylab <- switch(1+as.numeric(regexpr("signed", tolower(curr_ratio_type)) > 0), 
                    paste0(rel_ycoord - (histBreakStep*100)/2, "-", rel_ycoord + (histBreakStep*100)/2, "%"),
                    rel_ycoord)
  
  xlabpos <- switch(1+as.numeric(regexpr("signed", tolower(curr_ratio_type)) > 0), 
                    seq(2, 1/histBreakStep * 3, by=3),
                    seq(2, 1/histBreakStep * 6, by=3))
  
  # from 0 to (nGroup+1)*nbreaks
  box(bty="l")
  axis(1, at=xlabpos, labels = my_ylab)
  

  # axis(1, at = c(0,xlabpos, 30), labels=c("", my_ylab, ""))
  axis(2, at=c(0, max(rel_freqValues, na.rm=T)), labels=F)
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
  # go further if can divide by two
  halfStep <- 0.5* 1/histBreakStep
  if(regexpr("signed", tolower(curr_ratio_type)) > 0) {
    halfStep <- halfStep * 2
  }
  stopifnot( halfStep %% 1 == 0 )
  
  par(mfrow=c(1,1))
  ####### merged bar plot
  #barplot(rel_freqValues, beside=T, col=c(obsCol, permutCol), ylab="Frequency", xlab="% of down FC per TAD")
  #legend("topright", c( "observed Freq.", "randomized Freq."),
  #       text.col=c(obsColText, permutCol), bty='n', cex=1)
  #rel_ycoord <- (curr_ratio_breaks * 100)[-1]
  #rel_ycoord <- rel_ycoord-(histBreakStep*100)/2
  #my_ylab <- paste0(rel_ycoord - (histBreakStep*100)/2, "-", rel_ycoord + (histBreakStep*100)/2, "%")
  #xlabpos <- seq(2, 1/histBreakStep * 3, by=3)
  ## from 0 to (nGroup+1)*nbreaks
  #box(bty="l")
  #axis(1, at = c(xlabpos), labels=c(my_ylab))
  ##### after merging
  do.call(plotType, list(outFile2, height=myHeight, width=myWidth))
  rel_freq_1stpart <- rel_freqValues[,1:halfStep]
  if(regexpr("signed", tolower(curr_ratio_type)) > 0) {
    rel_freq_2ndpart <- rel_freqValues[,(halfStep + 1):(2*1/histBreakStep)]
  } else{
    rel_freq_2ndpart <- rel_freqValues[,(halfStep + 1):(1/histBreakStep)] 
  }
  rev_rel_freq_2ndpart <- rel_freq_2ndpart[,rev(1:ncol(rel_freq_2ndpart))]
  merged_freq <- rel_freq_1stpart + rev_rel_freq_2ndpart
  myxlab <- ifelse(regexpr("signed", tolower(curr_ratio_type)) > 0, curr_ratio_type, paste0("% of ", curr_ratio_type, " per TAD"))
  barplot(merged_freq, beside=T, col=c(obsCol, permutCol), 
          ylab="Frequency", xlab=myxlab, border=NA)
  legend("topleft", c( "observed Freq.", "randomized Freq."),
         text.col=c(obsColText, permutCol), bty='n', cex=1)
  rel_ycoord <- (curr_ratio_breaks * 100)[-1]
  rel_ycoord <- rel_ycoord[1:halfStep]
  rel_ycoord <- rel_ycoord-(histBreakStep*100)/2
  
  # plot the % only if not signed
  if(regexpr("signed", tolower(curr_ratio_type)) > 0) {
    rel_ycoord <- rel_ycoord/100
  }
  
  my_ylab <- switch(1+as.numeric(regexpr("signed", tolower(curr_ratio_type)) > 0), 
                    paste0(rel_ycoord - (histBreakStep*100)/2, "-", rel_ycoord + (histBreakStep*100)/2, "%"),
                    rel_ycoord
                )
                    
  xlabpos <- switch(1+as.numeric(regexpr("signed", tolower(curr_ratio_type)) > 0),
                    seq(2, 1/histBreakStep * 3, by=3)[1:halfStep],
                    seq(2, 1/histBreakStep * 3, by=3)[1:halfStep])
  # from 0 to (nGroup+1)*nbreaks
  box(bty="l")

  axis(1, at = c(xlabpos), labels=c(my_ylab))
  
  lines(x = xlabpos-0.5, y=merged_freq["observed", ], col=obsCol)
  lines(x = xlabpos+0.5, y=merged_freq["randomized", ], col=permutCol)
  foo <- dev.off()
  cat(paste0("... written: ", outFile2, "\n"))
  
  cat("... range number of genes by TAD: ", paste0(range(nbrGenes), collapse = "-"), "\n")

} # end iterating over ratio type

txt <- paste0(startTime, "\n", Sys.time(), "\n")
printAndLog(txt, pipLogFile)

cat(paste0("*** DONE: ", script_name, "\n"))




