inFile <- paste0("OUTPUT_FOLDER/SCORE_DT/17e_build_score_table_limited/rank_all_scores_DT_withSamp.txt")

scoresDT <- read.delim(inFile, stringsAsFactors = F, header=T)

scoresDT_noOut <- scoresDT[-which(scoresDT$dataset == "TCGAbrca_lum_bas"),]

all_values <- c(scoresDT$prodSignedRatio_auc, scoresDT$ratioDown_auc, scoresDT$rescWeightedQQ_auc, scoresDT$rescWeighted_auc)

outFold <- paste0("SCORES_COMPARISON")
system(paste0("mkdir -p ", outFold))
plotType <- "svg"
myHeight <- 7
myWidth <- 7

############################################################## COMPARE PAIRS OF SCORES

outFile <- paste0(outFold, "/prodSignedRatio_vs_ratioDown.", plotType)
do.call(plotType, list(outFile, height = myHeight, width=myWidth))
plot(prodSignedRatio_auc ~ratioDown_auc, data = scoresDT, 
     main = "Comparison AUC prodSignedRatio vs. ratioDown",
     xlim = range(all_values),
     ylim = range(all_values),
     xlab = "AUC obs/perm [cumsum ratioConcord-0.5]",
     ylab = "AUC obs/perm [cumsum prodSignedRatio]",
     pch=16, cex=0.7 )
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- paste0(outFold, "/rescWeightedQQ_vs_ratioDown.", plotType)
do.call(plotType, list(outFile, height = myHeight, width=myWidth))
plot(rescWeightedQQ_auc ~ratioDown_auc, data = scoresDT, 
     main = "Comparison AUC rescWeightedQQ vs. ratioDown",
     xlim = range(all_values),
     ylim = range(all_values),
     xlab = "AUC obs/perm [cumsum ratioConcord-0.5]",
     ylab = "AUC obs/perm [cumsum rescWeightedQQConcord-0.5]",
     pch=16, cex=0.7 )
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- paste0(outFold, "/rescWeighted_vs_ratioDown.", plotType)
do.call(plotType, list(outFile, height = myHeight, width=myWidth))
plot(rescWeighted_auc ~ratioDown_auc, data = scoresDT, 
     main = "Comparison AUC rescWeighted vs. ratioDown",
     xlim = range(all_values),
     ylim = range(all_values),
     xlab = "AUC obs/perm [cumsum ratioConcord-0.5]",
     ylab = "AUC obs/perm [cumsum rescWeightedConcord-0.5]",
     pch=16, cex=0.7 )
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- paste0(outFold, "/rescWeightedQQ_vs_prodSignedRatio.", plotType)
do.call(plotType, list(outFile, height = myHeight, width=myWidth))
plot(rescWeightedQQ_auc ~prodSignedRatio_auc, data = scoresDT, 
     main = "Comparison AUC rescWeightedQQ vs. prodSignedRatio",
     xlim = range(all_values),
     ylim = range(all_values),
     xlab = "AUC obs/perm [cumsum prodSignedRatio]",
     ylab = "AUC obs/perm [cumsum rescWeightedQQConcord-0.5]",
     pch=16, cex=0.7 )
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- paste0(outFold, "/rescWeighted_vs_prodSignedRatio.", plotType)
do.call(plotType, list(outFile, height = myHeight, width=myWidth))
plot(rescWeighted_auc ~prodSignedRatio_auc, data = scoresDT, 
     main = "Comparison AUC rescWeighted vs. prodSignedRatio",
     xlim = range(all_values),
     ylim = range(all_values),
     xlab = "AUC obs/perm [cumsum prodSignedRatio]",
     ylab = "AUC obs/perm [cumsum rescWeightedConcord-0.5]",
     pch=16, cex=0.7 )
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))

outFile <- paste0(outFold, "/rescWeightedQQ_vs_rescWeighted.", plotType)
do.call(plotType, list(outFile, height = myHeight, width=myWidth))
plot(rescWeightedQQ_auc ~rescWeighted_auc, data = scoresDT, 
     main = "Comparison AUC rescWeightedQQ vs. rescWeighted",
     xlim = range(all_values),
     ylim = range(all_values),
     xlab = "AUC obs/perm [cumsum rescWeightedConcord-0.5]",
     ylab = "AUC obs/perm [cumsum rescWeightedQQConcord-0.5]",
     pch=16, cex=0.7 )
foo <- dev.off()
cat(paste0("... written: ", outFile, "\n"))


############################################################## COMPARE SCORES AND # OF SAMPLES

all_scores <- paste0(c("ratioDown", "prodSignedRatio", "rescWeightedQQ", "rescWeighted"), "_auc")

nColPlot <- 3
nRowPlot <- 1

for(curr_score in all_scores){
  outFile <- paste0(outFold, "/", curr_score, "_vs_samplesNbrRatio.", plotType)
  do.call(plotType, list(outFile, width = nColPlot*myWidth, height = nRowPlot*myHeight))
  par(mfrow=c(nRowPlot,nColPlot))
  plot(scoresDT[,curr_score] ~ I(scoresDT$nSamp1+scoresDT$nSamp2),
       main = paste0("AUC ", gsub("_auc", "", curr_score) , " and total # of samples"),
       ylim = range(all_values),
       xlab = "Total # of samples (nSamp1+nSamp2)",
       ylab = paste0("AUC obs/perm ", gsub("_auc", "", curr_score)),
       pch=16, cex=0.7 )
  plot(scoresDT_noOut[,curr_score] ~ I(scoresDT_noOut$nSamp1+scoresDT_noOut$nSamp2),
       main = paste0("AUC ", gsub("_auc", "", curr_score) , " and total # of samples (without TCGAbrca)"),
       ylim = range(all_values),
       xlab = "Total # of samples (nSamp1+nSamp2)",
       ylab = paste0("AUC obs/perm ", gsub("_auc", "", curr_score)),
       pch=16, cex=0.7 )
  plot(scoresDT[,curr_score] ~ I(scoresDT$nSamp1/scoresDT$nSamp2),
       main = paste0("AUC ", gsub("_auc", "", curr_score) , " and ratio # of samples"),
       ylim = range(all_values),
       xlab = "Ratio # of samples (nSamp1/nSamp2)",
       ylab = paste0("AUC obs/perm ", gsub("_auc", "", curr_score)),
       pch=16, cex=0.7 )
  # plot(scoresDT_noOut[,curr_score] ~ I(scoresDT_noOut$nSamp1/scoresDT_noOut$nSamp2),
  #      main = paste0("AUC ", gsub("_auc", "", curr_score) , " and ratio # of samples (without TCGAbrca)"),
  #      ylim = range(all_values),
  #      xlab = "Ratio # of samples (nSamp1/nSamp2)",
  #      ylab = paste0("AUC obs/perm ", gsub("_auc", "", curr_score)),
  #      pch=16, cex=0.7 )
  # 
  foo <- dev.off()
  cat(paste0("... written: ", outFile, "\n"))
  
}

