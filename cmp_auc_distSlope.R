library(foreach)
library(doMC)


SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

registerDoMC(ifelse(SSHFS, 2, 20))

# caller <- "DI"

myHeight <- myWidth <- 7

outFold <- "cmp_auc_distSlope"
system(paste0("mkdir -p ", outFold))

# all_callers <- c("DI", "TopDom")
all_callers <- "TopDom"
all_ratios <- c("ratioDown", "prodSignedRatio") 
all_perms <- c("shuffTADs", "permGenes")

nVectEmt <- 14

for(caller in all_callers) {
    
  main_outFold <- paste0(setDir, "/mnt/ed4/marie/scripts/", "TAD_DE_pipeline_v2_", caller, "/OUTPUT_FOLDER/")
  
  all_outFolders <- list.files(main_outFold, full.names = T)
  stopifnot(length(all_outFolders) > 0)
  
  all_ratio_datasets <- foreach(myfold = all_outFolders) %dopar% {

    ratioAUCfile <- paste0(myfold, "/170_score_auc_pval_withShuffle/allratio_auc_pval.Rdata")
    if(file.exists(ratioAUCfile)) {
      allratio_auc_pval <- eval(parse(text = load(paste0(ratioAUCfile))))  
    }else{
      allratio_auc_pval <- NULL
    }
    
    distSlopeFile_shuffTADs <- paste0(myfold, "/14i4_cumulAllDown_distOnly_randomTADsShuffle/slopes_obsPermDist.Rdata")
    if(file.exists(distSlopeFile_shuffTADs)) {
      all_ratio_split_shuffTADs <- eval(parse(text = load(paste0(distSlopeFile_shuffTADs))))  
    }else{
      all_ratio_split_shuffTADs <- NULL
    }
    names(all_ratio_split_shuffTADs) <- paste0(names(all_ratio_split_shuffTADs), "_shuffTADs")
    
    distSlopeFile_permGenes <- paste0(myfold, "/14f4_cumulAllDown_distOnly/slopes_obsPermDist.Rdata")
    if(file.exists(distSlopeFile_permGenes)) {
      all_ratio_split_permGenes <- eval(parse(text = load(paste0(distSlopeFile_permGenes))))  
    }else{
      all_ratio_split_permGenes <- NULL
    }
    names(all_ratio_split_permGenes) <- paste0(names(all_ratio_split_permGenes), "_permGenes")

    c(allratio_auc_pval, all_ratio_split_permGenes, all_ratio_split_shuffTADs)
  }
  stopifnot(length(all_ratio_datasets) == length(all_outFolders))
  names(all_ratio_datasets) <- basename(all_outFolders)
  
  all_ratio_datasets <- Filter(function(x) length(x) == nVectEmt, all_ratio_datasets)
  stopifnot(length(all_ratio_datasets) > 0)
  
  for(curr_perm in all_perms) {
    for(curr_ratio_type in all_ratios) {
      
      ratioAUC_fullCurve <- sapply(all_ratio_datasets, function(x)  {
        if(is.null(x)) return(NA)
        if(length(x) < nVectEmt) return(NA)
        x[[paste0(curr_ratio_type, "_auc_", curr_perm)]] 
      })
      ratioAUC_fullCurve <- ratioAUC_fullCurve[!is.na(ratioAUC_fullCurve)]
      stopifnot(length(ratioAUC_fullCurve) > 0)
      
      slope_fullCurve <- sapply(all_ratio_datasets, function(x)  {
        if(is.null(x)) return(NA)
        if(length(x) < nVectEmt) return(NA)
        x[[paste0("slope_fullCurve_", curr_perm)]] 
      })
      slope_fullCurve <- slope_fullCurve[!is.na(slope_fullCurve)]
      stopifnot(length(slope_fullCurve) > 0)
      
      slope_firstPart <- sapply(all_ratio_datasets, function(x)  {
        if(is.null(x)) return(NA)
        if(length(x) < nVectEmt) return(NA)
        x[[paste0("slope_firstPart_", curr_perm)]] 
      })
      slope_firstPart <- slope_firstPart[!is.na(slope_firstPart)]
      stopifnot(length(slope_firstPart) > 0)
      
      slope_scdPart <- sapply(all_ratio_datasets, function(x)  {
        if(is.null(x)) return(NA)
        if(length(x) < nVectEmt) return(NA)
        x[[paste0("slope_scdPart_", curr_perm)]] 
      })
      slope_scdPart <- slope_scdPart[!is.na(slope_scdPart)]
      stopifnot(length(slope_scdPart) > 0)
      
      curr_datasets <- Reduce(intersect, list(names(ratioAUC_fullCurve), names(slope_firstPart), names(slope_scdPart), names(slope_fullCurve)))
  
      ratioAUC_fullCurve <- ratioAUC_fullCurve[curr_datasets]
      slope_firstPart <- slope_firstPart[curr_datasets]
      slope_scdPart <- slope_scdPart[curr_datasets]
      slope_fullCurve <- slope_fullCurve[curr_datasets]
      
      outFile <- file.path(outFold, paste0(caller, "_", curr_ratio_type, "_", curr_perm, "_", "ratioAUC_vs_slopePart1", ".svg"))
      svg(outFile, width = myWidth, height = myHeight)
      plot(x = slope_firstPart, y = ratioAUC_fullCurve,
           xlab = "slope obs-perm dist 1st part",
           ylab = "ratio AUC full curve",
           main = "Slope dist obs-perm 1st part vs. ratio AUC full",
           cex=0.7, pch=16)
      mtext(side = 3, paste0(caller, " - ", curr_ratio_type,  " - ", curr_perm))
      text(x=slope_firstPart, y = ratioAUC_fullCurve, labels = curr_datasets, cex=0.8)
      corTxt <- paste0("PCC = ", round(cor(slope_firstPart, ratioAUC_fullCurve), 4), "\nn = ", length(curr_datasets))
      legend("topleft", legend = corTxt, bty="n")
      foo <- dev.off()
      cat(paste0("... written: ", outFile, "\n"))
      
      outFile <- file.path(outFold, paste0(caller, "_", curr_ratio_type, "_", curr_perm, "_", "ratioAUC_vs_slopePart2", ".svg"))
      svg(outFile, width = myWidth, height = myHeight)
      plot(x = slope_scdPart, y = ratioAUC_fullCurve,
           xlab = "slope obs-perm dist 2nd part",
           ylab = "ratio AUC full curve",
           main = "Slope dist obs-perm 2nd part vs. ratio AUC full",
           cex=0.7, pch=16)
      mtext(side = 3, paste0(caller, " - ", curr_ratio_type,  " - ", curr_perm))
      text(x=slope_scdPart, y = ratioAUC_fullCurve, labels = curr_datasets, cex=0.8)
      corTxt <- paste0("PCC = ", round(cor(slope_scdPart, ratioAUC_fullCurve), 4), "\nn = ", length(curr_datasets))
      legend("topleft", legend = corTxt, bty="n")
      foo <- dev.off()
      cat(paste0("... written: ", outFile, "\n"))
      
      outFile <- file.path(outFold, paste0(caller, "_", curr_ratio_type, "_", curr_perm, "_", "ratioAUC_vs_slopeFullCurve", ".svg"))
      svg(outFile, width = myWidth, height = myHeight)
      plot(x = slope_fullCurve, y = ratioAUC_fullCurve,
           xlab = "slope obs-perm dist full curve",
           ylab = "ratio AUC full curve",
           main = "Slope dist obs-perm full curve vs. ratio AUC full",
           cex=0.7, pch=16)
      mtext(side = 3, paste0(caller, " - ", curr_ratio_type,  " - ", curr_perm))
      text(x=slope_fullCurve, y = ratioAUC_fullCurve, labels = curr_datasets, cex=0.8)
      corTxt <- paste0("PCC = ", round(cor(slope_fullCurve, ratioAUC_fullCurve), 4), "\nn = ", length(curr_datasets))
      legend("topleft", legend = corTxt, bty="n")
      foo <- dev.off()
      cat(paste0("... written: ", outFile, "\n"))
      

    }
  }
  
}
  










