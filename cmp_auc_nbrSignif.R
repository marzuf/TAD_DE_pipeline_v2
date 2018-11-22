library(foreach)
library(doMC)


SSHFS <- F
setDir <- ifelse(SSHFS, "/media/electron", "")

registerDoMC(ifelse(SSHFS, 2, 20))

# caller <- "DI"

myHeight <- myWidth <- 7

outFold <- "cmp_auc_nbrSignif"
system(paste0("mkdir -p ", outFold))

# all_callers <- c("DI", "TopDom")
all_callers <- "TopDom"
all_ratios <- c("ratioDown", "prodSignedRatio") 
all_perms <- c("shuffTADs", "permGenes")


for(caller in all_callers) {
    
  main_outFold <- paste0(setDir, "/mnt/ed4/marie/scripts/", "TAD_DE_pipeline_v2_", caller, "/OUTPUT_FOLDER/")
  
  all_outFolders <- list.files(main_outFold, full.names = T)
  stopifnot(length(all_outFolders) > 0)
  
  signif_all_datasets_DT <- foreach(myfold = all_outFolders, .combine='rbind') %dopar% {

    ratioAUCfile <- paste0(myfold, "/170_score_auc_pval_withShuffle/allratio_auc_pval.Rdata")
    allratio_auc_pval <- eval(parse(text = load(paste0(ratioAUCfile))))  
    
    ratioDown_auc_permGenes <- allratio_auc_pval["ratioDown_auc_permGenes"]
    ratioDown_auc_shuffTADs <- allratio_auc_pval["ratioDown_auc_shuffTADs"]
    prodSignedRatio_auc_permGenes <- allratio_auc_pval["prodSignedRatio_auc_permGenes"]
    prodSignedRatio_auc_shuffTADs <- allratio_auc_pval["prodSignedRatio_auc_shuffTADs"]
    
    signifPvalFile <- paste0(myfold, "/20_getSignifCombined/signifTADs.Rdata")
    signifTADs_adjPvalCombined <- eval(parse(text = load(paste0(signifPvalFile))))
    nSignifPval <- length(signifTADs_adjPvalCombined)
    
    thresh_signifPvalFile <- paste0(myfold, "/20_getSignifCombined/adjCombinedPval_thresh.Rdata")
    signifPval_thresh <- as.numeric(eval(parse(text = load(paste0(thresh_signifPvalFile)))))
    
    signifEmpFDRfile <- paste0(myfold, "/21_getSignifFDR/signifTADs.Rdata")
    signifTADs_empFDR <- eval(parse(text = load(paste0(signifEmpFDRfile))))
    nSignifFDR <- length(signifTADs_empFDR)
    
    thresh_FDRfile <- paste0(myfold, "/21_getSignifFDR/all_empFDR_signifThresh.Rdata")
    empFDR_thresh <- eval(parse(text = load(paste0(thresh_FDRfile))))
    empFDR_thresh <- as.numeric(empFDR_thresh["desired_empFDR_thresh"])
    
    
    data.frame(
      dataset=basename(myfold),
      ratioDown_auc_permGenes=ratioDown_auc_permGenes,
      ratioDown_auc_shuffTADs=ratioDown_auc_shuffTADs,
      prodSignedRatio_auc_permGenes=prodSignedRatio_auc_permGenes,
      prodSignedRatio_auc_shuffTADs=prodSignedRatio_auc_shuffTADs,
      nSignifPval = nSignifPval,
      signifPval_thresh = signifPval_thresh,
      nSignifFDR = nSignifFDR,
      empFDR_thresh =empFDR_thresh,
      nIntersectSignif = length(intersect(signifTADs_adjPvalCombined,signifTADs_empFDR )),
      stringsAsFactors = F
    )
    
  }
  rownames(signif_all_datasets_DT) <- NULL
  outFile <- paste0(outFold, "/", "signif_all_datasets_DT.Rdata")
  save(signif_all_datasets_DT, file=outFile)
  cat(paste0("... written: ", outFile, "\n"))
  
    
  stopifnot(length(unique(signif_all_datasets_DT$signifPval_thresh)) == 1)
  stopifnot(length(unique(signif_all_datasets_DT$empFDR_thresh)) == 1)
  
  pvalComb_thresh <- unique(signif_all_datasets_DT$signifPval_thresh)
  empFDR_thresh <- unique(signif_all_datasets_DT$empFDR_thresh)
  
  for(curr_ratio_type in all_ratios) {
    
    for(curr_perm in all_perms) {
      
      curr_x <- signif_all_datasets_DT[,paste0(curr_ratio_type, "_auc_", curr_perm)]
      
      # PLOT THE SIGNIF PVAL COMBINED
      outFile <- file.path(outFold, paste0(caller, "_", curr_ratio_type, "_", curr_perm, "_", "nbrSignifPvalcombined_vs_ratioAUC", ".svg"))
      svg(outFile, width = myWidth, height = myHeight)
      plot(x = curr_x, y = signif_all_datasets_DT[,"nSignifPval"],
           xlab = paste0(curr_ratio_type, " - AUC ratio"),
           ylab = "nbr signif TADs adj. combined Pval",
           main = "Nbr signif Pval combined vs. ratio AUC",
           cex=0.7, pch=16)
      mtext(side = 3, paste0(caller, " - ", curr_ratio_type,  " - ", curr_perm, " - signif. thresh=", pvalComb_thresh))
      # text(x=curr_x, y = signif_all_datasets_DT[,"nSignifPval"], labels = signif_all_datasets_DT$dataset, cex=0.8)
      corTxt <- paste0("PCC = ", round(cor(curr_x, signif_all_datasets_DT[,"nSignifPval"]), 4), "\nn = ", nrow(signif_all_datasets_DT))
      legend("topleft", legend = corTxt, bty="n")
      foo <- dev.off()
      cat(paste0("... written: ", outFile, "\n"))
      
      # PLOT THE SIGNIF EMP FDR
      outFile <- file.path(outFold, paste0(caller, "_", curr_ratio_type, "_", curr_perm, "_", "nbrSignifEmpFDR_vs_ratioAUC", ".svg"))
      svg(outFile, width = myWidth, height = myHeight)
      plot(x = curr_x, y = signif_all_datasets_DT[,"nSignifFDR"],
           xlab = paste0(curr_ratio_type, " - AUC ratio"),
           ylab = "nbr signif TADs emp. FDR",
           main = "Nbr signif emp. FDR vs. ratio AUC",
           cex=0.7, pch=16)
      mtext(side = 3, paste0(caller, " - ", curr_ratio_type,  " - ", curr_perm, " - signif. thresh=", empFDR_thresh))
      # text(x=curr_x, y = signif_all_datasets_DT[,"nSignifFDR"], labels = signif_all_datasets_DT$dataset, cex=0.8)
      corTxt <- paste0("PCC = ", round(cor(curr_x, signif_all_datasets_DT[,"nSignifFDR"]), 4), "\nn = ", nrow(signif_all_datasets_DT))
      legend("topleft", legend = corTxt, bty="n")
      foo <- dev.off()
      cat(paste0("... written: ", outFile, "\n"))
      
      
      # PLOT THE INTERSECT
      outFile <- file.path(outFold, paste0(caller, "_", curr_ratio_type, "_", curr_perm, "_", "nbrIntersectSignif_vs_ratioAUC", ".svg"))
      svg(outFile, width = myWidth, height = myHeight)
      plot(x = curr_x, y = signif_all_datasets_DT[,"nIntersectSignif"],
           xlab = paste0(curr_ratio_type, " - AUC ratio"),
           ylab = "nbr signif TADs intersect Pval/FDR",
           main = "Nbr signif intersect vs. ratio AUC",
           cex=0.7, pch=16)
      mtext(side = 3, paste0(caller, " - ", curr_ratio_type,  " - ", curr_perm))
      # text(x=curr_x, y = signif_all_datasets_DT[,"nIntersectSignif"], labels = signif_all_datasets_DT$dataset, cex=0.8)
      corTxt <- paste0("PCC = ", round(cor(curr_x, signif_all_datasets_DT[,"nIntersectSignif"]), 4), "\nn = ", nrow(signif_all_datasets_DT))
      legend("topleft", legend = corTxt, bty="n")
      foo <- dev.off()
      cat(paste0("... written: ", outFile, "\n"))
      
    }
  }
} # end iterating over caller










