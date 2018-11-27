suppressPackageStartupMessages(library(Hmisc, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) # error bar

#################################################################3

plot_smile <- function(obsvect_data_down,
                       permutDT_data_down,
                       nRandom,
                       TADonly,
                       g2tDT,
                       mytit = NULL,
                       cexMain=1,
                       # currOutFold,
                       # setDir = "", 
                       # script1_name="1_runGeneDE",
                       histBreakStep = 0.1,
                       toPlot = NULL,
                       obsCol = "bisque",
                       obsColText = "bisque2",
                       permutCol = "mediumseagreen",
                       signedStat = FALSE) {
  
  curr_ratio_breaks <- switch(1+as.numeric(signedStat), seq(0,1,histBreakStep), seq(-1,1,histBreakStep))

  stopifnot( (1/histBreakStep) %% 1 == 0 )

  if(is.null(mytit)) mytit <- "Smile plot"
  
  if(is.null(toPlot)) {
    toPlot <- "both"
  } else{
    stopifnot(toPlot == "top" | toPlot == "bottom")
  }
  cat(paste0("*** START smile plot \n"))  
  
  ### !!! NEED TO CHECK WHY THERE ARE SOME NA
  permutDT_data_down <- na.omit(permutDT_data_down)
  
  intersectRegions <- intersect(names(obsvect_data_down), rownames(permutDT_data_down))
  g2tDT <- g2tDT[g2tDT$region %in% intersectRegions,]
  # some regions could have been discarded because no gene
  intersectRegions <- intersect(intersectRegions, g2tDT$region)
  
  ### if needed, filter TAD only regions
  if(TADonly) {
    initLen <- length(intersectRegions)
    intersectRegions <- intersectRegions[grep("_TAD", intersectRegions)]
    txt <- paste0(toupper(script_name), "> Take only TAD regions: ", length(intersectRegions), "/", initLen, "\n")
    cat(txt)    
  }

  initNrow <- nrow(permutDT_data_down)
  permutDT_data_down <- permutDT_data_down[intersectRegions,]
  txt <- paste0(toupper(script_name), "> Take only the regions form permut for which genes were retained: ", nrow(permutDT_data_down), "/", initNrow, "\n")
  cat(txt)
  
  initLen <- length(obsvect_data_down)
  obsvect_data_down <- obsvect_data_down[intersectRegions]
  txt <- paste0(toupper(script_name), "> Take only the regions form observ. for which genes were retained: ", length(obsvect_data_down), "/", initLen, "\n")
  cat(txt)
  
  stopifnot(nrow(permutDT_data_down) == length(obsvect_data_down))
  
  ###################################################################################################
  
  cat("... Warning - considering only TAD regions: \n")
  
  cat("> Prepare simulation data \n")
  nTAD <- sum(regexpr("_TAD",  rownames(permutDT_data_down)) > 0)
  cat(paste0("... TADs: ", nTAD, "/", nrow(permutDT_data_down), "\n"))
  permut_currDown_vect <- as.vector(permutDT_data_down[grep("_TAD", rownames(permutDT_data_down)),])
  stopifnot(length(permut_currDown_vect) == ncol(permutDT_data_down) * nTAD)
  
  cat(paste0("length(permut_currDown_vect) = ", length(permut_currDown_vect), "\n"))
  
  cat("> Prepare observed data \n")
  nTAD <- sum(regexpr("_TAD",  names(obsvect_data_down)) > 0)
  cat(paste0("... TADs: ", nTAD, "/", length(obsvect_data_down), "\n"))
  obsvect_data_down <- obsvect_data_down[grep("_TAD", names(obsvect_data_down))]
  
  cat("> Prepare data for plotting \n")
  
  obs_rel_hist <- try(hist(obsvect_data_down, breaks=curr_ratio_breaks, plot=F))
  #plot(obs_rel_hist)
  # some are not bounded 0-1 => go to next directly !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"WILL NEED REFINEMENT!!!
  if(class(obs_rel_hist) == "try-error") stop("cannot plot smile plot ! values not in the right range\n")
  
  cat(paste0("length(permut_currDown_vect) = ", length(permut_currDown_vect), "\n"))
  cat(paste0("length(permut_currDown_vect) = ", length(na.omit(permut_currDown_vect)), "\n"))
  shuff_rel_hist <-  hist(permut_currDown_vect,breaks=curr_ratio_breaks, plot=F)
  # plot(shuff_rel_hist)
  
  cat(paste0("ncol(permutDT_data_down) = ", ncol(permutDT_data_down), "\n"))
  cat(paste0("sum(shuff_rel_hist$counts) = ", sum(shuff_rel_hist$counts), "\n"))
  cat(paste0("sum(shuff_rel_hist$counts/nRandom) = ", sum(shuff_rel_hist$counts/nRandom), "\n"))
  cat(paste0("sum(obs_rel_hist$counts) = ", sum(obs_rel_hist$counts),"\n"))
  stopifnot(sum(shuff_rel_hist$counts/nRandom) == sum(obs_rel_hist$counts))   # 830
  stopifnot(sum(shuff_rel_hist$counts) == length(permut_currDown_vect))
  stopifnot(sum(obs_rel_hist$counts) == length(na.omit(obsvect_data_down)))
  
  rel_freqValues <- rbind(obs_rel_hist$counts, shuff_rel_hist$counts/nRandom)
  rownames(rel_freqValues) <- c("observed", "randomized")
  
  # FOR THE ERROR BARS:
  # for each randomization, calculate the observed, then divide expected by observed
  # table 1 column = 1 randomization, rows = breaks of the hist.
  rel_HistValues <- apply(permutDT_data_down[grep("_TAD", rownames(permutDT_data_down)),],
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
  # DE_topTable <- eval(parse(text = load(paste0(currOutFold, "/", script1_name, "/DE_topTable.Rdata"))))
  # DE_geneList <- eval(parse(text = load(paste0(currOutFold, "/", script1_name, "/DE_geneList.Rdata"))))
  # exprDT <- eval(parse(text = load(paste0(currOutFold, "/", script1_name, "/DE_rnaseqDT.Rdata"))))
  # samp1 <- eval(parse(text=load(paste0(setDir, "/", sample1_file))))
  # samp2 <- eval(parse(text=load(paste0(setDir, "/", sample2_file))))
  # stopifnot(all(DE_topTable$genes %in% names(DE_geneList)))
  # stopifnot(!any(duplicated(names(DE_geneList))))
  # stopifnot(all(colnames(exprDT) %in% c(samp1, samp2)))
  # stopifnot(all(samp1 %in% colnames(exprDT)))
  # stopifnot(all(samp2 %in% colnames(exprDT)))
  # maxDownGene <- DE_topTable$genes[which.min(DE_topTable$logFC)]
  # stopifnot(maxDownGene %in% rownames(exprDT))
  # mean_expr1 <- mean(unlist(c(exprDT[maxDownGene, samp1])), na.rm=T)
  # mean_expr2 <- mean(unlist(c(exprDT[maxDownGene, samp2])), na.rm=T)
  
  # if(mean_expr1 > mean_expr2) {
  #   subtitDir <- paste0("down: ", toupper(cond1), " > ", toupper(cond2))
  # } else{
  #   subtitDir <- paste0("down: ", toupper(cond2), " > ", toupper(cond1))
  # }
  
  if(signedStat)
    rel_ycoord <- rel_ycoord/100
  
  ####### START SMILE PLOT HERE
  if(toPlot == "both") {
    par(mfrow = c(2,1))  
  }else{
    # par(mfrow = c(1,1))
  }
  if(toPlot == "both" | toPlot == "top") {
    myxlab <- ifelse(signedStat, curr_ratio_type, paste0("% of ", curr_ratio_type, " per TAD"))
    # top panel: line plot
    plot(rel_logRatio ~ rel_ycoord, type="l",axes=F, cex.main=cexMain,
         main = mytit,
         xlab = myxlab, ylab = "log2(Observed/Randomized)")
    # mtext(paste0("smile score = ", smile_score, "; ", subtitDir))
    mtext(paste0("smile score = ", smile_score))
    box(bty="l")
    axis(2)
    if(signedStat){
      axis(1)
    }else{
      axis(1, at=rel_ycoord, labels = my_ylab) 
    }
    abline(h=0, lty=2)
    errbar(x=rel_ycoord, y=rel_logRatio, yplus=rel_logRatio+rel_obsOverExp_sem, yminus=rel_logRatio-rel_obsOverExp_sem, add=T)
  }
  
  if(toPlot == "both" | toPlot == "bottom") {
   if(toPlot == "both") mytit <- ""
     # bottom panel: bar plot
    barplot(rel_freqValues, beside=T, col=c(obsCol, permutCol), ylab="Frequency", xlab=paste0("% of ", curr_ratio_type, " per TAD"), main =mytit)
    legend("topright", c( "observed Freq.", "randomized Freq."),
           text.col=c(obsColText, permutCol), bty='n', cex=1)
    # legend("topleft", legend = subtitDir,
    #        bty='n', cex=1)
    # xlabpos <- seq(2, 30, by=3)
    rel_ycoord <- (curr_ratio_breaks * 100)[-1]
    rel_ycoord <- rel_ycoord-(histBreakStep*100)/2
    my_ylab <- paste0(rel_ycoord - (histBreakStep*100)/2, "-", rel_ycoord + (histBreakStep*100)/2, "%")
    xlabpos <- seq(2, 1/histBreakStep * 3, by=3)
    # from 0 to (nGroup+1)*nbreaks
    box(bty="l")
    axis(1, at = c(xlabpos), labels=c(my_ylab))
    # axis(1, at = c(0,xlabpos, 30), labels=c("", my_ylab, ""))
    axis(2, at=c(0, max(rel_freqValues, na.rm=T)), labels=F)
  }
}

  

