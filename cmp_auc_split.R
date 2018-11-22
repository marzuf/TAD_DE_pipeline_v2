library(foreach)
library(doMC)


SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

registerDoMC(ifelse(SSHFS, 2, 20))

# caller <- "DI"

myHeight <- myWidth <- 7

outFold <- "cmp_auc_split"
system(paste0("mkdir -p ", outFold))

# all_callers <- c("DI", "TopDom")
all_callers <- "TopDom"
all_ratios <- c("ratioDown", "prodSignedRatio") 
all_perms <- c("shuffTADs", "permGenes")

for(caller in all_callers) {
    
  main_outFold <- paste0(setDir, "/mnt/ed4/marie/scripts/", "TAD_DE_pipeline_v2_", caller, "/OUTPUT_FOLDER/")
  
  all_outFolders <- list.files(main_outFold, full.names = T)
  stopifnot(length(all_outFolders) > 0)
  
  all_ratio_datasets <- foreach(myfold = all_outFolders) %dopar% {
    
    # all_ratio_norm <- load("OUTPUT_FOLDER/TCGAbrca_lum_bas/14i2splitSlope_cumulAllDown_limited_AUC_randomTADsShuffle/all_ratio_split.Rdata
    
    permGenefile <- paste0(myfold, "/14f2split_cumulAllDown_limited_AUC/all_ratio_split.Rdata")
    
    if(file.exists(permGenefile)) {
      all_ratio_permGenes <- eval(parse(text = load(paste0(permGenefile))))  
    }else {
      all_ratio_permGenes <- NULL
    }
    
    shuffTADfile <- paste0(myfold, "/14i2split_cumulAllDown_limited_AUC_randomTADsShuffle/all_ratio_split.Rdata")
    
    if(file.exists(shuffTADfile)) {
      all_ratio_shuffTADs <- eval(parse(text = load(paste0(shuffTADfile))))  
    }else{
      all_ratio_shuffTADs <- NULL
    }
    
    c(all_ratio_permGenes, all_ratio_shuffTADs)
  }
  stopifnot(length(all_ratio_datasets) == length(all_outFolders))
  names(all_ratio_datasets) <- basename(all_outFolders)
  
  all_ratio_datasets <- Filter(function(x) length(x) == 12, all_ratio_datasets)
  stopifnot(length(all_ratio_datasets) > 0)
  
  for(curr_perm in all_perms) {
    for(curr_ratio_type in all_ratios) {
      
      auc_first_and_second <- sapply(all_ratio_datasets, function(x)  {
        if(is.null(x)) return(NA)
        if(length(x) < 12) return(NA)
        x[[paste0(curr_ratio_type, "_auc_", curr_perm, "_all")]] 
        })
      auc_first_and_second <- auc_first_and_second[!is.na(auc_first_and_second)]
        
      auc_first <- sapply(all_ratio_datasets, function(x)  {
        if(is.null(x)) return(NA)
        if(length(x) < 12) return(NA)
        x[[paste0(curr_ratio_type, "_auc_", curr_perm, "_one")]] 
      })
      auc_first <- auc_first[!is.na(auc_first)]
      
      auc_second <- sapply(all_ratio_datasets, function(x)  {
        if(is.null(x)) return(NA)
        if(length(x) < 12) return(NA)
        x[[paste0(curr_ratio_type, "_auc_", curr_perm, "_2ndPart")]] 
      })
      auc_second <- auc_second[!is.na(auc_second)]
      
      curr_datasets <- Reduce(intersect, list(names(auc_first_and_second), names(auc_first), names(auc_second)))
  
      auc_first_and_second <- auc_first_and_second[curr_datasets]
      auc_first <- auc_first[curr_datasets]
      auc_second <- auc_second[curr_datasets]
      
      outFile <- file.path(outFold, paste0(caller, "_", curr_ratio_type, "_", curr_perm, "_", "part1_vs_all", ".svg"))
      svg(outFile, width = myWidth, height = myHeight)
      plot(x = auc_first_and_second, y = auc_first,
           xlab = "AUC ratio both parts",
           ylab = "AUC ratio 1-part",
           main = "AUC ratio 1 only vs. AUC ratio full",
           cex=0.7, pch=16)
      mtext(side = 3, paste0(caller, " - ", curr_ratio_type,  " - ", curr_perm))
      text(x=auc_first_and_second, y = auc_first, labels = curr_datasets, cex=0.8)
      corTxt <- paste0("PCC = ", round(cor(auc_first_and_second, auc_first), 4))
      legend("topleft", legend = corTxt, bty="n")
      foo <- dev.off()
      cat(paste0("... written: ", outFile, "\n"))
      
      outFile <- file.path(outFold, paste0(caller, "_", curr_ratio_type, "_", curr_perm, "_", "part2_vs_all", ".svg"))
      svg(outFile, width = myWidth, height = myHeight)
      plot(x = auc_first_and_second, y = auc_second,
           xlab = "AUC ratio both parts",
           ylab = "AUC ratio 2nd part",
           main = "AUC ratio 2nd part only vs. AUC ratio full",
           cex=0.7, pch=16)
      mtext(side = 3, paste0(caller, " - ", curr_ratio_type,  " - ", curr_perm))
      text(x=auc_first_and_second, y = auc_second, labels = curr_datasets, cex=0.8)
      corTxt <- paste0("PCC = ", round(cor(auc_first_and_second, auc_second), 4))
      legend("topleft", legend = corTxt, bty="n")
      foo <- dev.off()
      cat(paste0("... written: ", outFile, "\n"))
      
      
      outFile <- file.path(outFold, paste0(caller, "_", curr_ratio_type, "_", curr_perm, "_", "part2_vs_part1", ".svg"))
      svg(outFile, width = myWidth, height = myHeight)
      plot(x = auc_first, y = auc_second,
           xlab = "AUC ratio 1-part",
           ylab = "AUC ratio 2nd part",
           main = "AUC ratio 2nd part only vs. AUC ratio 1 only",
           cex=0.7, pch=16)
      mtext(side = 3, paste0(caller, " - ", curr_ratio_type,  " - ", curr_perm))
      text(x=auc_first, y = auc_second, labels = curr_datasets, cex=0.8)
      corTxt <- paste0("PCC = ", round(cor(auc_first, auc_second), 4))
      legend("topleft", legend = corTxt, bty="n")
      foo <- dev.off()
      cat(paste0("... written: ", outFile, "\n"))
      
      ratio_second_first <- auc_second/auc_first
      outFile <- file.path(outFold, paste0(caller, "_", curr_ratio_type, "_", curr_perm, "_", "part2overpart1_vs_all", ".svg"))
      svg(outFile, width = myWidth, height = myHeight)
      plot(x = auc_first_and_second, y = ratio_second_first,
           xlab = "AUC ratio both parts",
           ylab = "AUC ratio 2nd part/AUC ratio 1-part",
           main = "Ratio AUC ratio 2nd part/1-part vs. AUC ratio full",
           cex=0.7, pch=16)
      mtext(side = 3, paste0(caller, " - ", curr_ratio_type,  " - ", curr_perm))
      text(x=auc_first_and_second, y = ratio_second_first, labels = curr_datasets, cex=0.8)
      corTxt <- paste0("PCC = ", round(cor(auc_first_and_second, ratio_second_first), 4))
      legend("topleft", legend = corTxt, bty="n")
      foo <- dev.off()
      cat(paste0("... written: ", outFile, "\n"))
      
      
      
      
    }
  }
  
}
  










