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
