#############################################################################################################################
############################################################################################################################# 
#############################################################################################################################

get_meanCorr_value <- function(exprMatrix, inside_genes, outside_genes, cormet) {
  stopifnot(inside_genes %in% rownames(exprMatrix))
  stopifnot(outside_genes %in% rownames(exprMatrix))
  stopifnot(setequal(c(inside_genes, outside_genes), rownames(exprMatrix)))
  
  nAllGenes <- length(inside_genes) + length(outside_genes)
  
  coexprMatrix <- cor(t(exprMatrix), method = cormet)
  stopifnot(dim(coexprMatrix) == nAllGenes)
  
  coexprMatrix[lower.tri(coexprMatrix, diag = TRUE)] <- NA   # because after I filter that 1 gene should be inside, and 1 gene should be outside -> can never happen to take the diag. value of coexpression
  coexprMatrix <- na.omit(melt(coexprMatrix))
  colnames(coexprMatrix)[1:2] <- c("Var1", "Var2")
  stopifnot(colnames( coexprMatrix)[3] == "value" )
  coexprMatrix$Var1 <- as.character(coexprMatrix$Var1)
  coexprMatrix$Var2 <- as.character(coexprMatrix$Var2)
  
  stopifnot(coexprMatrix$Var1 %in% outside_genes | coexprMatrix$Var1 %in% inside_genes)
  stopifnot(coexprMatrix$Var2 %in% outside_genes | coexprMatrix$Var2 %in% inside_genes)
  stopifnot(inside_genes %in% coexprMatrix$Var1 | inside_genes %in% coexprMatrix$Var2)
  stopifnot(outside_genes %in% coexprMatrix$Var1 | outside_genes %in% coexprMatrix$Var2)
  
  # take only if one of the two genes outside and the other inside
  coexprMatrix <- coexprMatrix[!  (coexprMatrix$Var1 %in% outside_genes & coexprMatrix$Var2 %in% outside_genes),]  # do not take correlation between pairs of genes in  the same TAD
  coexprMatrix <- coexprMatrix[ ! (coexprMatrix$Var1 %in% inside_genes & coexprMatrix$Var2 %in% inside_genes),]  # do not take correlation between pairs of genes in  the same TAD
  
  stopifnot(     (coexprMatrix$Var1 %in% outside_genes & coexprMatrix$Var2 %in% inside_genes) | 
                   (coexprMatrix$Var2 %in% outside_genes & coexprMatrix$Var1 %in% inside_genes) )
  
  

  stopifnot(nrow(coexprMatrix) == (length(inside_genes) * length(outside_genes) ))
  
  meanCorr_value <- mean(coexprMatrix$value)
  stopifnot(!is.na(meanCorr_value))
  return(meanCorr_value)
}


#############################################################################################################################
############################################################################################################################# ADDED 26.11.2019
#############################################################################################################################

get_meanCorrPartial_value <- function(exprMatrix, inside_genes, outside_genes, purity_vals, cormet, nause) {
  stopifnot(inside_genes %in% rownames(exprMatrix))
  stopifnot(outside_genes %in% rownames(exprMatrix))
  stopifnot(setequal(c(inside_genes, outside_genes), rownames(exprMatrix)))
  
  nAllGenes <- length(inside_genes) + length(outside_genes)
  
  full_coexprMatrix <- cor(t(exprMatrix), method = cormet)
  stopifnot(dim(full_coexprMatrix) == nAllGenes)
  
  #### CORRELATION ===> change here for partial correlation (26.11.2019)
  fpkmdt <- exprMatrix
  stopifnot(colnames(fpkmdt) %in% names(purity_vals))
  fpkmdt_with_pur <- as.data.frame(fpkmdt)
  fpkmdt_with_pur["purity",] <- purity_vals[colnames(fpkmdt_with_pur)]
  stopifnot(rownames(fpkmdt_with_pur)[nrow(fpkmdt_with_pur)] == "purity")
  
  
  
  # partial correlations for a set (x) of variables with set (y) removed.
  partialcorr <- as(partial.r(t(fpkmdt_with_pur),
                              x = 1:(nrow(fpkmdt_with_pur)-1),
                              y=nrow(fpkmdt_with_pur),
                              method=cormet,
                              use=nause), "matrix")
  
  
  stopifnot(isSymmetric(partialcorr))
  stopifnot(rownames(partialcorr) == rownames(full_coexprMatrix))
  stopifnot(colnames(partialcorr) == colnames(full_coexprMatrix))
  
  
  coexprMatrix <- partialcorr
  # => end added here
    
  coexprMatrix[lower.tri(coexprMatrix, diag = TRUE)] <- NA   # because after I filter that 1 gene should be inside, and 1 gene should be outside -> can never happen to take the diag. value of coexpression
  coexprMatrix <- na.omit(melt(coexprMatrix))
  colnames(coexprMatrix)[1:2] <- c("Var1", "Var2")
  stopifnot(colnames( coexprMatrix)[3] == "value" )
  coexprMatrix$Var1 <- as.character(coexprMatrix$Var1)
  coexprMatrix$Var2 <- as.character(coexprMatrix$Var2)
  
  stopifnot(coexprMatrix$Var1 %in% outside_genes | coexprMatrix$Var1 %in% inside_genes)
  stopifnot(coexprMatrix$Var2 %in% outside_genes | coexprMatrix$Var2 %in% inside_genes)
  stopifnot(inside_genes %in% coexprMatrix$Var1 | inside_genes %in% coexprMatrix$Var2)
  stopifnot(outside_genes %in% coexprMatrix$Var1 | outside_genes %in% coexprMatrix$Var2)
  
  # take only if one of the two genes outside and the other inside
  coexprMatrix <- coexprMatrix[!  (coexprMatrix$Var1 %in% outside_genes & coexprMatrix$Var2 %in% outside_genes),]  # do not take correlation between pairs of genes in  the same TAD
  coexprMatrix <- coexprMatrix[ ! (coexprMatrix$Var1 %in% inside_genes & coexprMatrix$Var2 %in% inside_genes),]  # do not take correlation between pairs of genes in  the same TAD
  
  stopifnot(     (coexprMatrix$Var1 %in% outside_genes & coexprMatrix$Var2 %in% inside_genes) | 
                   (coexprMatrix$Var2 %in% outside_genes & coexprMatrix$Var1 %in% inside_genes) )
  
  

  stopifnot(nrow(coexprMatrix) == (length(inside_genes) * length(outside_genes) ))
  
  meanCorr_value <- mean(coexprMatrix$value)
  stopifnot(!is.na(meanCorr_value))
  return(meanCorr_value)
}


#############################################################################################################################
############################################################################################################################# 
#############################################################################################################################


get_SAM_FDR_aroundTADs <- function(obs_vect, permut_values, cut_off, symDir, nPermut = NULL, variableName = "", 
                                   withPlot=T, plotOffsetY = 0, minQuant = 0.05, maxQuant = 0.95, inputList=FALSE, sortToPlot=FALSE) {
  
  stopifnot(symDir %in% c("symmetric", "higher", "lower"))
  # N = number of real solutions higher than cut-off
  # then, for each permutation (column), how many solutions greater than k ?
  # R = average number of random solutions greater than k
  # get the average number of items called signif from the permut
  # (estimated number of false discoveries = average # of genes called signif from the permut.)
  
  # => with the values from sampling across boundaries, no need to take the average
  
  ## NEED TO DIVIDE BY NUMBER OF SAMPLING
  if(is.null(nPermut)) {
    nPseudoPermut <- length(permut_values)/length(obs_vect)  # ????????????????????? 
  } else { 
    nPseudoPermut <- nPermut
    stopifnot(is.numeric(nPermut)) 
  }
  
  
  if(symDir == "symmetric") {
    observ_N <- sum(abs(obs_vect) >= abs(cut_off))  
    random_R <- sum(abs(permut_values) >= abs(cut_off))/nPseudoPermut
  } else if(symDir == "higher"){
    observ_N <- sum(obs_vect >= cut_off)
    random_R <-  sum(permut_values >= cut_off)/nPseudoPermut
  } else if (symDir =="lower"){
    observ_N <- sum(obs_vect <= cut_off)
    random_R <-  sum(permut_values <= cut_off)/nPseudoPermut
  } else{
    stop("should never happen\n")
  }
  
  
  empFDR <- random_R/observ_N
  # then the empirical FDR is R/N
  # = average signif from permut / # observed signif
  
  if(withPlot){
    if(inputList) {
      warning("... cannot draw for Shuffle data\n")
    } else{
      simuColArea <- rgb(0.2, 0.2, 1, 0.3)
      
      sampledRange <- quantile(permut_values, probs = c(minQuant, maxQuant))
      nValues <- length(obs_vect)
      
      if(sortToPlot) {
        plot_vect <- sort(obs_vect, decreasing=TRUE)
      } else {
        plot_vect <- obs_vect
      }
      
      plot(plot_vect, cex=0.7, pch=16, xlab="",
           main = variableName, 
           ylab= paste0(variableName, " observed values"), 
           xlim = c(0, nValues), ylim = c(min(obs_vect)-plotOffsetY, max(obs_vect)+plotOffsetY),
           axes=F, bty="l")
      box(bty="l")
      axis(2)
      mtext(paste0("empirical FDR = ", round(empFDR, 2), " % (observed signif: ", observ_N, ")"), line=-0.5)
      abline(h = cut_off, lty=2, col = "firebrick3")
      text(x = 10, y = (cut_off + 0.1*cut_off) , paste0("cut-off = ", cut_off),  adj = c(0,0), col = "firebrick3", offset=2)
      if(symDir == "symmetric")
        abline(h = -cut_off, lty=2, col = "firebrick3")
      polygon(c(rev(1:nValues), 1:nValues), c(rev( rep(sampledRange[1], nValues)), rep(sampledRange[2], nValues)), 
              col = simuColArea, border = NA)
      legend("topleft", legend = paste0(minQuant, "-", maxQuant, " quantile sampling"), col=simuColArea, bty="n", text.col=simuColArea)
      
    }
  }
  
  
  invisible(empFDR)
}


#############################################################################################################################
#############################################################################################################################
#############################################################################################################################

plot_FDR_with_observedSignif <- function(yaxis_empFDR_vect, xaxis_cutoff, y2_obsSignif, variableName,  feature_name = "genes", rightAxisCol = "steelblue") {
  initmar <- par()$mar
  par(mar = c(5,5,2,5))
  plot((100*yaxis_empFDR_vect) ~ xaxis_cutoff, type="o", pch=16, 
       main=paste0("Empirical FDR with variable cut-off of ", variableName),
       xlab=paste0(variableName, " cut-off"),
       ylab="% FDR", axes=F)
  axis(1, at=xaxis_cutoff, labels = xaxis_cutoff)
  axis(2)
  par(new = T)
  plot(log10(y2_obsSignif) ~ xaxis_cutoff, type="o", 
       col=rightAxisCol, pch=16, ylim=c(0, max(log10(y2_obsSignif))),
       axes=F, xlab=NA, ylab=NA, bty="l")
  box(bty="u")
  axis(side = 4, col=rightAxisCol)
  mtext(side = 4, line = 3, paste0('Observed # of significant ', feature_name, ' (log10)'), col = rightAxisCol)
  on.exit(par(initmar))
}

