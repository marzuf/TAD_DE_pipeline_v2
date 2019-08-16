

################################################################# !!!!!!!!!!!!!!!!!!! warning hard-coded
# - genome assembly used with biomart: GRCh37
##########################


#######################################################################################################################
#######################################################################################################################
#######################################################################################################################

quantNorm <- function(x) {qqnorm(x, plot.it=F)$x}

#######################################################################################################################
#######################################################################################################################
#######################################################################################################################

madNorm <- function(x) {
  rna_rnaseqDT_tmp <- x
  # median center
  all_meds <- apply(rna_rnaseqDT_tmp, 1, median)
  rna_rnaseqDT_tmp <- sweep(rna_rnaseqDT_tmp, 1, all_meds, "-")  # substract from each row their corresponding median
  # mad norm
  # In order to use the MAD as a consistent estimator for the estimation of the standard deviation σ, one takes
  # σ ^ = k ⋅ MAD where k is a constant scale factor, which depends on the distribution.[1]
  # For normally distributed data k is taken to be: ~1.4826
  all_mads <-1.4826 * apply(abs(rna_rnaseqDT_tmp), 1, median)
  rna_madnorm_rnaseqDT <- sweep(rna_rnaseqDT_tmp, 1, all_mads, "/") 
  stopifnot( all ( dim(x) == dim(rna_madnorm_rnaseqDT)))
  rownames(rna_madnorm_rnaseqDT) <- rownames(x)
  colnames(rna_madnorm_rnaseqDT) <- colnames(x)
  return(rna_madnorm_rnaseqDT)
}

#######################################################################################################################
#######################################################################################################################
#######################################################################################################################  
get_geneList_fromEntrez <- function(refList, g2t, histDT_file) {
  # look if find other history mapping
  # replace if some of the input entrezID have a correspondance to the newest one stored in entrezDT 
  historyDT <- read.delim(paste0(histDT_file), header=TRUE, stringsAsFactors = F)
  historyDT$entrezID <- as.character(historyDT$entrezID)
  historyDT$mappingID <- as.character(historyDT$mappingID)
  refList_v2 <- unlist(sapply(refList, function(x) ifelse (x %in% historyDT$mappingID,
                                                                              historyDT$entrezID[historyDT$mappingID == x], x)))
  # the name is now refList, but the value not necessarily geneList
  if(!is.null(g2t))
    refList_v2 <- refList_v2[refList_v2 %in% g2t$entrezID]
  # remove if any duplicated
  dupID <- unique(refList_v2[which(duplicated(refList_v2))])
  refList_v2 <- refList_v2[!refList_v2 %in% dupID]
  stopifnot(all(names(refList_v2) %in% refList))
  stopifnot(!any(duplicated(refList_v2)))
  # return in correct order
  refList <- refList[refList %in% names(refList_v2)]
  returnList <- refList_v2[refList]
  stopifnot(!any(duplicated(returnList)))
  stopifnot(!any(duplicated(names(returnList))))
  returnList != names(returnList)
  return(returnList)
}


#######################################################################################################################
#######################################################################################################################
#######################################################################################################################  

# refList is the list of symbol
get_geneList_fromSymbol <- function(refList, g2t, symbDT_file) {
  symbDT <- read.delim(paste0(symbDT_file), header=TRUE, stringsAsFactors = F)
  symbDT$entrezID <- as.character(symbDT$entrezID)
  symbDT <- symbDT[symbDT$symbol %in% refList,]
  if(!is.null(g2t))
    symbDT <- symbDT[symbDT$entrezID %in% g2t$entrezID, ]
  
  dubSymb <- unique(symbDT$symbol[which(duplicated(symbDT$symbol))])
  symbDT <- symbDT[!symbDT$symbol %in% dubSymb, ]
  dubEntrez <- unique(symbDT$entrezID[which(duplicated(symbDT$entrezID))])
  symbDT <- symbDT[!symbDT$entrezID %in% dubEntrez, ]
  refList_v2 <- setNames(symbDT$entrezID, symbDT$symbol)
  refList <- refList[refList %in% names(refList_v2)]
  returnList <- refList_v2[refList]
  stopifnot(!any(duplicated(returnList)))
  stopifnot(!any(duplicated(names(returnList))))
  return(returnList)
}

#######################################################################################################################
#######################################################################################################################
#######################################################################################################################  

get_geneList_fromEnsembl <- function(refList, g2t, ensDT_file) {

  ensDT <- read.delim(paste0(ensDT_file), header=TRUE, stringsAsFactors = F)
  ensDT$entrezID <- as.character(ensDT$entrezID)
  ensDT <- ensDT[ensDT$ensemblID %in% refList,]
  
  if(!is.null(g2t))
    ensDT <- ensDT[ensDT$entrezID %in% g2t$entrezID, ]
  
  dubEns <- unique(ensDT$ensemblID[which(duplicated(ensDT$ensemblID))])
  ensDT <- ensDT[!ensDT$ensemblID %in% dubEns, ]
  
  dubEntrez <- unique(ensDT$entrezID[which(duplicated(ensDT$entrezID))])
  ensDT <- ensDT[!ensDT$entrezID %in% dubEntrez, ]
  refList_v2 <- setNames(ensDT$entrezID, ensDT$ensemblID)
  refList <- refList[refList %in% names(refList_v2)]
  returnList <- refList_v2[refList]
  stopifnot(!any(duplicated(returnList)))
  stopifnot(!any(duplicated(names(returnList))))
  return(returnList)
  

}

#######################################################################################################################
#######################################################################################################################
#######################################################################################################################

printAndLog <- function(mytext, mylogf) {
  cat(mytext)
  cat(mytext, append=T, file = mylogf)
}

#######################################################################################################################
#######################################################################################################################
#######################################################################################################################

get_multiShuffledPositions_vFunct <- function(g2TADdt, RNAdt, geneIDlist, nClass, withExprClass, TADonly, nSimu, nCpu, aggregFun) {
  suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))  
  suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))  
  if(withExprClass) {
    stopifnot(!is.null(RNAdt) & !is.null(nClass))
  }
  g2TADdt$entrezID <- as.character(g2TADdt$entrezID)
  g2TADdt$region <- as.character(g2TADdt$region)
  registerDoMC(nCpu)
  # need to ensure that I get the same order for the genes
  # do the first one
  allT <- get_ShuffledPositions_vFunct(g2TADdt = g2TADdt, RNAdt = RNAdt, geneIDlist = geneIDlist, 
                                      nClass = nClass, TADonly = TADonly, withExprClass = withExprClass, aggregFun=aggregFun) 
  colnames(allT) <- c(colnames(allT)[1], paste0(colnames(allT)[2], "1")) # region1
  genes1 <- allT$entrezID
  if(nSimu >1){
    tmpDT <- foreach(i=2:nSimu, .combine='cbind') %dopar% {
      if(withExprClass) {
        cat(paste0("... WITH CLASS ", aggregFun, " - shuffle: ", i, "/", nSimu, "\n"))
      } else{
        cat(paste0("... NO CLASS - shuffle: ", i, "/", nSimu, "\n"))
      }
      x <- get_ShuffledPositions_vFunct(g2TADdt = g2TADdt, RNAdt = RNAdt, geneIDlist = geneIDlist, 
                                       nClass = nClass, TADonly = TADonly, withExprClass = withExprClass, aggregFun=aggregFun,saveIdx=i) 
      stopifnot(all(genes1 == x[,1]))
      x[,2]
    }
    colnames(tmpDT) <- paste0("region", 2:nSimu)
    allT <- cbind(allT, tmpDT)
  }
  return(allT)  
}  

#######################################################################################################################
#######################################################################################################################
#######################################################################################################################

### same as _vJune but can pass aggregFun for aggregating the expression values

get_ShuffledPositions_vFunct <- function(g2TADdt, RNAdt, geneIDlist, nClass, TADonly, withExprClass, aggregFun, saveIdx=0) {
  warning("geneIDlist argument should correspond to rownames of RNAdt")
  warning("duplicated - ambiguous - are removed !!! ")
  duplicatedID <- geneIDlist[duplicated(geneIDlist)]
  RNAdt <- RNAdt[which(! geneIDlist %in% duplicatedID),]
  geneIDlist <- geneIDlist[which(! geneIDlist %in% duplicatedID)]
  rownames(RNAdt) <- geneIDlist
  if(withExprClass) {
    stopifnot(!is.null(RNAdt) & !is.null(nClass))
    stopifnot(!is.null(aggregFun))
  }
  g2TADdt$entrezID <- as.character(g2TADdt$entrezID)
  g2TADdt$region <- as.character(g2TADdt$region)
  # take only the genes for which we have their positions
  # and subset the rnaseq data for these genes only
  if(TADonly) {
    # take only the genes that are in TAD
    geneListTAD <- geneIDlist[geneIDlist  %in% g2TADdt$entrezID[grep("_TAD", g2TADdt$region)] ]
    RNAdt <- RNAdt[geneListTAD,]
  } else {
    geneListTAD <- geneIDlist[geneIDlist  %in% g2TADdt$entrezID]
    RNAdt <- RNAdt[geneListTAD,]
  }
  ##########
  ##### DO IT BY SHUFFLING THE LABELS BY CLASS OF EXPRESSION
  ##########
  if(withExprClass) {
    # define expression classes based on median expression
    geneAggregExpression <- data.frame(gene = geneListTAD, expValue = apply(RNAdt, 1, aggregFun))
    rownames(geneAggregExpression) <- NULL
    # rank by expression value
    geneAggregExpression <- geneAggregExpression[order(geneAggregExpression$expValue),]
    # split into 'nClass' groups 
    nR <- nrow(geneAggregExpression)
    nGeneByClass <- rep(nR %/% nClass, nClass)
    # add the remainder (1 to each, to avoid having a lot more in one class)
    if(nR%%nClass > 0)
      nGeneByClass[1:(nR%%nClass)] <- nR %/% nClass + 1
    stopifnot(sum(nGeneByClass) == nR)
    geneClass <- rep(c(1:nClass), nGeneByClass)
    stopifnot(length(geneClass) == nR)
    # add a column to DF with their corresponding class
    geneAggregExpression$class <- geneClass
    #save(geneAggregExpression, file="geneAggregExpression.Rdata")
    stopifnot(length(unique(geneAggregExpression$class)) == nClass)
    # add a column with their initial TAD
    geneAggregExpression$initRegion <- sapply(geneAggregExpression$gene, function(x) g2TADdt$region[g2TADdt$entrezID==x])
    # now, for each class, reshuffle the TAD -> new column with the reshuffled positions
    geneAggregExpression$shuffRegion <- foreach(i_cl = 1:nClass, .combine='c') %do% {
      subDT <- geneAggregExpression[geneAggregExpression$class == i_cl,]
      initPos <- subDT$initRegion
      newPos <- sample(initPos, size=length(initPos), replace=F)  
      stopifnot(all(unique(as.character(initPos)) %in% unique(as.character(newPos))))
      newPos
    }

save(geneAggregExpression, file=paste0("geneAggregExpression", saveIdx, ".Rdata"))

    shuffGenePosDT <- data.frame(entrezID = geneAggregExpression$gene, 
                                 region = as.character(geneAggregExpression$shuffRegion),
                                 stringsAsFactors = F)
  } else {
    ##########
    ##### JUST SHUFFLE THE LABELS, IRRESPECTIVE OF GENE EXPRESSION
    ##########
    newPos <- sample(g2TADdt$region[g2TADdt$entrezID %in% geneListTAD], replace=F)
    shuffGenePosDT <- data.frame(entrezID = geneListTAD, 
                                 region = as.character(newPos),
                                 stringsAsFactors = F)
  }
  # to be compatible with the functions that use gene2tadDT, return a similar DF
  return(shuffGenePosDT)
}




#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
# this version uses left_joining and aggregate, return the ratio for all regions
get_downByRegion_v2 <- function(g2TADdt, DEdt) {
  suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
  g2TADdt$entrezID <- as.character(g2TADdt$entrezID)
  g2TADdt$region <- as.character(g2TADdt$region)
  DEdt$genes <- as.character(DEdt$genes)
  g2TADdt <- g2TADdt[g2TADdt$entrezID %in% DEdt$genes,]
  mergedDT <- left_join(DEdt[,c("genes", "logFC")], g2TADdt[,c("entrezID", "region")], by=c("genes" = "entrezID"))
  mergedDT <- na.omit(mergedDT)
  nGenesDT <- aggregate(genes ~ region, data = mergedDT, FUN = length)
  nGenes <- setNames(nGenesDT$genes, nGenesDT$region)
  ratioDownDT <- aggregate(logFC ~ region, data = mergedDT, FUN = function(x) sum(x<0))
  ratioDown <- setNames(ratioDownDT$logFC, ratioDownDT$region)
  all_ratio <- ratioDown/nGenes[names(ratioDown)]
  stopifnot(all_ratio >= 0 & all_ratio <= 1)
  return(all_ratio)
}  

#######################################################################################################################
# instead of getting the ratio of genes that are down-regulated, returns the ratio of FC that is down 
get_FCdownByRegion_v2 <- function(g2TADdt, DEdt) {
  suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
  g2TADdt$entrezID <- as.character(g2TADdt$entrezID)
  g2TADdt$region <- as.character(g2TADdt$region)
  DEdt$genes <- as.character(DEdt$genes)
  g2TADdt <- g2TADdt[g2TADdt$entrezID %in% DEdt$genes,]
  mergedDT <- left_join(DEdt[,c("genes", "logFC")], g2TADdt[,c("entrezID", "region")], by=c("genes" = "entrezID"))
  mergedDT <- na.omit(mergedDT)
  stopifnot(is.numeric(mergedDT$logFC[1]))
  # get the absolute value of tot. logFC by region
  totLogFC <- aggregate(logFC ~ region, data = mergedDT, FUN = function(x) sum(abs(x)))
  totLogFC <- setNames(totLogFC$logFC, totLogFC$region)
  # get the absolute value of negative logFC  
  negLogFC <- aggregate(logFC ~ region, data = mergedDT, FUN = function(x) sum(abs(x[x<0])))
  negLogFC <- setNames(negLogFC$logFC, negLogFC$region)
  stopifnot(all(names(totLogFC) == names(negLogFC)))
  all_ratio <- negLogFC/totLogFC[names(negLogFC)]
  
  # CHANGED 24.01.2018: PROBLEM IF A TAD CONTAINS ALL GENES WITH 0 LOGFC (HAPPENS TOPDOM GSE71119 DEDIFF MFSM)
  # in this case ratioDown should be zero not NA
  poss_NA <- which(negLogFC==0 & totLogFC==0)
  
  if(length(poss_NA) > 0){
    stopifnot(is.na(all_ratio[poss_NA]))
    all_ratio[poss_NA] <- 0
  }
  
  stopifnot(all_ratio >= 0 & all_ratio <= 1)
  return(all_ratio)
}  


#######################################################################################################################
#######################################################################################################################
#######################################################################################################################

# !!!! WARNING PIPELINE VERSION GENE NAMES ARE ROWNAMES NOT 1ST COL

# in this version also I only take as reference the regions I have in the permut data  (not the regions in the gene2tadDT)
get_statFromShuffle_para  <- function(DEdt, shuffData, stat, geneIDlist=NULL, ncpu=2, TADonly=F) {
  stopifnot(stat %in% c("ratioDown", "FCdown"))
  suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
  suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
  registerDoMC(cores = ncpu)
  # ensure always have the same
  regions <- as.character(unique(shuffData[,1]))
  if(TADonly)
    regions <- regions[grep("_TAD", regions)]
  regions <- sort(regions)

  statDT <- foreach(i_perm = 1:ncol(shuffData), .combine='cbind') %dopar% {
	cat(paste0("...... ratioDown permut: ", i_perm, "/", ncol(shuffData), "\n"))
    shuff_g2TAD <- data.frame(entrezID = rownames(shuffData), region = shuffData[,i_perm], stringsAsFactors =F)
    shuff_g2TAD$entrezID <- as.character(shuff_g2TAD$entrezID)
    shuff_g2TAD$region <- as.character(shuff_g2TAD$region)
    if (stat == "ratioDown") {
      x <- get_downByRegion_v2(shuff_g2TAD, DEdt) 
    } else if(stat == "FCdown") { 
      x <- get_FCdownByRegion_v2(shuff_g2TAD, DEdt)     
    } else  {
      stop("should not happen\n")
    }
    x[regions]
  }
  stopifnot(all(rownames(statDT) == regions))
  return(statDT)
}


#######################################################################################################################
#######################################################################################################################
#######################################################################################################################

stouffer <- function(ps, two.tails=FALSE) {
	if(two.tails) ps = ps / 2
	# transform p-values into quantiles of standard normal
	qps = qnorm(1-ps, lower.tail = TRUE)
	# take the average of the quantiles
	iqps = sum(qps) / sqrt(length(qps))
	# find p-val of the average
	p = 1 - pnorm(iqps, lower.tail = TRUE)
	if(two.tails) p = 2 * p
	return(p)
}

#the individual p-values are transformed into the quantiles of a standard normal
#the p-val of the average of the quantiles is then found

# Fisher's method treates large and small p-values asymmetrically
# e.g. combining 0.999 and 0.001 gives 0.008
# is asymmetrically sensitive to small p-values

# Z-transform test: one-to-one mapping of the standard normal curve of the p-value
# Z = standard deviate (= a number drawn from from normal distribution with mean 0 and sd 1)
# the test converts p-values into standard normal deviates
# sum of the Z divided by square root of the number of tests has a std normal distribution
# no asymmetry problem

# weighted version of the Z-transform test: a weight can be assigned to each test
# if each test is given equal weight, this reduces to Z-transform test
# ideally, each study is weighted proportional to the inverse of its error variance
# (i.e. by the reciprocal of its squared standard error)
# more generally, the weights should be the inverse of the squared standard error fo the effect
# size estimate for each study

######################################################################################################################################################################################################
######################################################################################################################################################################################################
######################################################################################################################################################################################################

# abs(logRatio) * -signFirstDeriv if 1st half
# abs(logRatio) * signFirstDeriv if 2nd half
get_smileScore <- function(ylogRatioVect, xBreakVect, withPlot=FALSE) {
  stopifnot(length(ylogRatioVect) == length(xBreakVect))
  quadFunc_firstDeriv <- function(x, mymodel) {
    my_coef <- coef(mymodel)
    my_intercept <- as.numeric(my_coef["(Intercept)"])
    my_x2 <-  as.numeric(my_coef[grep("^I\\(.+\\^2\\)", names(my_coef))])
    my_x <- as.numeric(my_coef[- c(grep("Intercept", names(my_coef)), grep("^I\\(.+\\^2\\)", names(my_coef)))])
    (2*my_x2) * x + my_x
  }
  # fit 2nd order linear model
  smileMod <- lm(ylogRatioVect ~ xBreakVect + I(xBreakVect^2))
  signFirstDeriv <- sign(quadFunc_firstDeriv(xBreakVect, smileMod))
  stopifnot(length(ylogRatioVect) == length(signFirstDeriv))
  mid <- length(xBreakVect)/2
  if(length(xBreakVect) %% 2 == 0 ) {
    firstH <- -signFirstDeriv[1:mid] * abs(ylogRatioVect[1:mid])
    scdH <- signFirstDeriv[(mid+1):length(signFirstDeriv)] * abs(ylogRatioVect[(mid+1):length(ylogRatioVect)])
    stopifnot(length(firstH) == length(scdH))
    stopifnot(length(firstH) == mid)
  } else {
    warning("not even number of breaks; take half of middle break in both halves\n")
    # e.g. if 5, firstHalf => 1-2, half of 3; secondHalf => 4-5, half of 3
    firstH <- -signFirstDeriv[1:floor(mid)] * abs(ylogRatioVect[1:floor(mid)])
    firstH_last <- (-signFirstDeriv[ceiling(mid)] * abs(ylogRatioVect[ceiling(mid)]))/2
    firstH <- c(firstH, firstH_last)
    stopifnot(length(firstH) == ceiling(mid))
    scdH <- signFirstDeriv[(ceiling(mid) +1):length(signFirstDeriv)] * abs(ylogRatioVect[(ceiling(mid) +1):length(ylogRatioVect)])
    scdH_last <- (signFirstDeriv[ceiling(mid)] * abs(ylogRatioVect[ceiling(mid)]))/2
    scdH <- c(scdH, scdH_last)
    stopifnot(length(scdH) == ceiling(mid))
  }
  smileScore <- sum(c(firstH, scdH))
  if(withPlot) {
    plot(ylogRatioVect ~ xBreakVect, type="l", bty="l", main=paste0("smile score = ", round(smileScore, 2)))
    lines(xBreakVect, predict(smileMod), col="red")
  }
  return(smileScore)
}


######################################################################################################################################################################################################
######################################################################################################################################################################################################
######################################################################################################################################################################################################
# NEWER VERSION WITH INDICATION OF GENE RANK IN X AXIS LABEL
# 
# plot_lolliTAD <- function(TAD_to_plot, 
#                           meanExpr_vect,
#                           DE_table,
#                           g2t_table,
#                           id2name_table,
#                           geneList,
#                           orderBy = "logFC",
#                           upcolor = "limegreen", downcolor = "orangered",
#                           palettefunc_pos = colorRampPalette(c("#e6ffe6", "#00b300")),
#                           palettefunc_neg = colorRampPalette(c("#ffe6e6", "#b30000")), 
#                           graybars = FALSE,
#                           textLeft = FALSE,
# 						  cond1 = "cond1",
# 						  cond2 = "cond2"
#                           ) {
# 
#   suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
#   suppressPackageStartupMessages(library(ggplot2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
#   suppressPackageStartupMessages(library(plotrix, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
#   suppressPackageStartupMessages(library(gridExtra, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
#   suppressPackageStartupMessages(library(ggpubr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
# 
#   stopifnot("entrezID" %in% colnames(id2name_table))
#   stopifnot("geneName" %in% colnames(id2name_table))
#   
#   rownames(id2name_table) <- id2name_table$entrezID
#   rownames(g2t_table) <- g2t_table$entrezID
#   rownames(DE_table) <- DE_table$genes
# 
#   if(!orderBy %in% c("logFC", "meanExpr_vect", "startPos"))
#     stop("ERROR - invalid \"orderBy\" argument !\n")
#   negPalette <- palettefunc_neg(11)
#   posPalette <- palettefunc_pos(11)
#   
#   genes_to_plot <- g2t_table$entrezID[g2t_table$region == TAD_to_plot]
#   genes_to_plot <- genes_to_plot[genes_to_plot %in% geneList]
#   if(length(genes_to_plot) == 0) {
#     warning("ERROR - no genes to plot \n")
# 	return(grob())
#    }
#   TAD_genes_DT <- data.frame(gene = genes_to_plot, 
#                              symbol = id2name_table[genes_to_plot, "geneName"],
#                              start = g2t_table[genes_to_plot, "start"],
#                              mean_expr = meanExpr_vect[genes_to_plot],
#                              log_FC = DE_table[genes_to_plot, "logFC"],
#                              voom_adj_pval = DE_table[genes_to_plot, "adj.P.Val"])
#   
#   stopifnot(!any(is.na(TAD_genes_DT)))
#   
#   TAD_genes_DT$start <- as.numeric(as.character(TAD_genes_DT$start))
#   TAD_genes_DT$gene <- as.character(TAD_genes_DT$gene)
#   
#   if(orderBy == "logFC") {
#     TAD_genes_DT <- TAD_genes_DT[order(TAD_genes_DT$log_FC),]
#   } else if(orderBy == "meanExpr") {
#     TAD_genes_DT <- TAD_genes_DT[order(TAD_genes_DT$mean_expr),]
#   } else if(orderBy == "startPos") {
#     TAD_genes_DT <- TAD_genes_DT[order(TAD_genes_DT$start),]
#   } else {
#     stop("error")
#   }   
#   
#   TAD_genes_DT$gene <- factor(TAD_genes_DT$gene, levels = as.character(TAD_genes_DT$gene))
#   
#   TAD_genes_DT$signif <-  unlist(sapply(TAD_genes_DT$voom_adj_pval, function(x)
#         ifelse(x < 0.001, "***", 
#                ifelse(x < 0.01, "**", 
#                       ifelse(x < 0.05, "*", "")))))
#   
#   TAD_genes_DT$logFC_color <- unlist(sapply(TAD_genes_DT$log_FC, function(x){
#     if(x < 0)
#       return(downcolor)
#     if(x > 0) 
#       return(upcolor)
#     return("black")
#   }))
#   
#   resc_mean_expr <- plotrix::rescale(TAD_genes_DT$mean_expr, newrange=c(0,1))
#   
#   TAD_genes_DT$meanExpr_color <- unlist(foreach(i = 1:nrow(TAD_genes_DT), .combine='c') %do% {
#     exprVal <- resc_mean_expr[i]
#     ifelse(TAD_genes_DT$log_FC[i] < 0, negPalette[round(exprVal*10)+1], 
#            ifelse(TAD_genes_DT$log_FC[i] > 0, posPalette[round(exprVal*10)+1], "gray"))
#   })
#   
#   if(graybars) {
#     barColor <- "lightgray"
#   } else{
#     barColor <- TAD_genes_DT$meanExpr_color
#   }
#   
#   # use gene entrez ID not symbols here (some genes duplicated name with different entrez ID)
#   p <- ggdotchart(TAD_genes_DT, x = "gene", y = "log_FC",
#                   title = TAD_to_plot,
#              color = TAD_genes_DT$logFC_color ,             # color of the dots
#              ylab = "Log2(fold change)",
#              xlab="",
#               add = "segments",                             # Add segments from y = 0 to dots
#              add.params = list(color = barColor, size = 2), # Change segment color and size
#              
#              # for correct sorting
#              order = as.character(TAD_genes_DT$gene),
#              sort.by.groups = FALSE,
#              group = "gene",
#              dot.size = 10,                                 # Large dot size
#              label = TAD_genes_DT$signif,                        # Add mpg values as dot labels
#               font.label = list(color = "white", size = 14, 
#                                vjust = 0.5),               # Adjust label parameters
#              ggtheme = theme_pubr()                        # ggplot2 theme
#   ) + 
#     geom_hline(yintercept = 0, linetype = 2, color = "lightgray") +
#     theme(plot.title = element_text(hjust = 0.5, face=2, size=18))
#   
#   # change the x labels to be the gene names
#  p <- p+scale_x_discrete(labels=TAD_genes_DT$symbol)
#   
#   p <- ggpar(p, ylim=c(-max(abs(TAD_genes_DT$log_FC), na.rm=T), max(abs(TAD_genes_DT$log_FC), na.rm=T)))
#   p <- p + annotate("text", x = ifelse(textLeft, 0 , nrow(TAD_genes_DT))+0.5, y = max(abs(TAD_genes_DT$log_FC), na.rm=T), 
#                     label = paste0("expr. ", cond1, " >\nexpr. ", cond2), color=upcolor, hjust =  ifelse(textLeft,0,1), fontface="bold")
#   p <- p + annotate("text", x = ifelse(textLeft, 0 , nrow(TAD_genes_DT))+0.5, y = -max(abs(TAD_genes_DT$log_FC), na.rm=T), 
#                     label = paste0("expr. ", cond2, " >\nexpr. ", cond1), color=downcolor, hjust = ifelse(textLeft,0,1), fontface="bold")
#   return(p)
# }


######################################################################################################################################################################################################
######################################################################################################################################################################################################
######################################################################################################################################################################################################

plot_lolliTAD <- function(TAD_to_plot, 
                          meanExpr_vect,
                          DE_table,
                          g2t_table,
                          id2name_table,
                          geneList,
                          orderBy = "logFC",
                          upcolor = "limegreen", 
                          downcolor = "orangered",
                          palettefunc_pos = colorRampPalette(c("#e6ffe6", "#00b300")),
                          palettefunc_neg = colorRampPalette(c("#ffe6e6", "#b30000")), 
                          graybars = FALSE,
                          textLeft = FALSE,
						  cond1 = "cond1",
						  cond2 = "cond2",
						  labelWithRank=FALSE,
						  mytitle=NULL
                          ) {

  suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
  suppressPackageStartupMessages(library(ggplot2, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
  suppressPackageStartupMessages(library(plotrix, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
  suppressPackageStartupMessages(library(gridExtra, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
  suppressPackageStartupMessages(library(ggpubr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 


library(rlang)
as_dictionary <- as_data_pronoun

  stopifnot("entrezID" %in% colnames(id2name_table))
  stopifnot("geneName" %in% colnames(id2name_table))
  
  rownames(id2name_table) <- id2name_table$entrezID
  rownames(g2t_table) <- g2t_table$entrezID
  rownames(DE_table) <- DE_table$genes



  if(!orderBy %in% c("logFC", "meanExpr_vect", "startPos"))
    stop("ERROR - invalid \"orderBy\" argument !\n")
  negPalette <- palettefunc_neg(11)
  posPalette <- palettefunc_pos(11)
  
  genes_to_plot <- g2t_table$entrezID[g2t_table$region == TAD_to_plot]
  genes_to_plot <- genes_to_plot[genes_to_plot %in% geneList]
  if(length(genes_to_plot) == 0) {
    warning("ERROR - no genes to plot \n")
	return(grob())
   }
  TAD_genes_DT <- data.frame(gene = genes_to_plot, 
                             symbol = id2name_table[genes_to_plot, "geneName"],
                             start = g2t_table[genes_to_plot, "start"],
                             mean_expr = meanExpr_vect[genes_to_plot],
                             log_FC = DE_table[genes_to_plot, "logFC"],
                             voom_adj_pval = DE_table[genes_to_plot, "adj.P.Val"])
  
  stopifnot(!any(is.na(TAD_genes_DT)))
  
  TAD_genes_DT$start <- as.numeric(as.character(TAD_genes_DT$start))
  TAD_genes_DT$gene <- as.character(TAD_genes_DT$gene)
  
  if(orderBy == "logFC") {
    TAD_genes_DT <- TAD_genes_DT[order(TAD_genes_DT$log_FC),]
  } else if(orderBy == "meanExpr") {
    TAD_genes_DT <- TAD_genes_DT[order(TAD_genes_DT$mean_expr),]
  } else if(orderBy == "startPos") {
    TAD_genes_DT <- TAD_genes_DT[order(TAD_genes_DT$start),]
  } else {
    stop("error")
  }   
  
  TAD_genes_DT$gene <- factor(TAD_genes_DT$gene, levels = as.character(TAD_genes_DT$gene))
  
  TAD_genes_DT$signif <-  unlist(sapply(TAD_genes_DT$voom_adj_pval, function(x)
        ifelse(x < 0.001, "***", 
               ifelse(x < 0.01, "**", 
                      ifelse(x < 0.05, "*", "")))))
  
  TAD_genes_DT$logFC_color <- unlist(sapply(TAD_genes_DT$log_FC, function(x){
    if(x < 0)
      return(downcolor)
    if(x > 0) 
      return(upcolor)
    return("black")
  }))
  
  ### ADDED FOR THE RANK - 11.04.18
  DE_table <- DE_table[order(DE_table$adj.P.Val),]

  TAD_genes_DT$gene_rank <-  unlist(sapply(TAD_genes_DT$gene, function(x) which(as.character(DE_table$genes) == as.character(x)) ))

  resc_mean_expr <- plotrix::rescale(TAD_genes_DT$mean_expr, newrange=c(0,1))
  
  TAD_genes_DT$meanExpr_color <- unlist(foreach(i = 1:nrow(TAD_genes_DT), .combine='c') %do% {
    exprVal <- resc_mean_expr[i]
    ifelse(TAD_genes_DT$log_FC[i] < 0, negPalette[round(exprVal*10)+1], 
           ifelse(TAD_genes_DT$log_FC[i] > 0, posPalette[round(exprVal*10)+1], "gray"))
  })
  
  if(graybars) {
    barColor <- "lightgray"
  } else{
    barColor <- TAD_genes_DT$meanExpr_color
  }
  
  my_xlab <- ifelse(labelWithRank, paste0("# genes =  ",  nrow(DE_table)), "")

#cat("AAA\n")
  
  lolliTitle <- ifelse(is.null(mytitle), TAD_to_plot, mytitle) 
  
  
  # use gene entrez ID not symbols here (some genes duplicated name with different entrez ID)
  p <- ggdotchart(TAD_genes_DT, x = "gene", y = "log_FC",
                  title = lolliTitle,
             color = TAD_genes_DT$logFC_color ,             # color of the dots
             ylab = "Log2(fold change)",
             xlab=my_xlab,
              add = "segments",                             # Add segments from y = 0 to dots
             add.params = list(color = barColor, size = 2), # Change segment color and size
             
             # for correct sorting
             order = as.character(TAD_genes_DT$gene),
             sort.by.groups = FALSE,
             group = "gene",
             dot.size = 10,                                 # Large dot size
             label = TAD_genes_DT$signif,                        # Add mpg values as dot labels
              font.label = list(color = "white", size = 14, 
                               vjust = 0.5),               # Adjust label parameters
             ggtheme = theme_pubr()                        # ggplot2 theme
  ) + 
    geom_hline(yintercept = 0, linetype = 2, color = "lightgray") +
    theme(plot.title = element_text(hjust = 0.5, face=2, size=18),
          axis.title.x = element_text(hjust = 0.5, vjust = 1, face=3, size=10)
    )

#cat("BBB\n")
  
# change the x labels to be the gene names (ok because TAD_genes_DT has been sorted !!!)
 if(labelWithRank){
   TAD_genes_DT$symbol <- paste0(TAD_genes_DT$symbol, "\n(", TAD_genes_DT$gene_rank, ")")
 }
 p <- p+scale_x_discrete(labels=TAD_genes_DT$symbol)
  
  p <- ggpar(p, ylim=c(-max(abs(TAD_genes_DT$log_FC), na.rm=T), max(abs(TAD_genes_DT$log_FC), na.rm=T)))
  p <- p + annotate("text", x = ifelse(textLeft, 0 , nrow(TAD_genes_DT))+0.5, y = max(abs(TAD_genes_DT$log_FC), na.rm=T), 
                    label = paste0("expr. ", cond1, " >\nexpr. ", cond2), color=upcolor, hjust =  ifelse(textLeft,0,1), fontface="bold")
  p <- p + annotate("text", x = ifelse(textLeft, 0 , nrow(TAD_genes_DT))+0.5, y = -max(abs(TAD_genes_DT$log_FC), na.rm=T), 
                    label = paste0("expr. ", cond2, " >\nexpr. ", cond1), color=downcolor, hjust = ifelse(textLeft,0,1), fontface="bold")
  return(p)
}



######################################################################################################################################################################################################
######################################################################################################################################################################################################
######################################################################################################################################################################################################


# draw the diagram
# Dark2 and Set3 are fine for example for the palette
# vd_tit =  in the argument to set the title of the VD plot
draw_venn <- function(...){
  set_to_cmp <- list(...)
  vd_tit <- set_to_cmp[["vd_tit"]]
  set_to_cmp[["vd_tit"]] <- NULL
  vd <- venn.diagram(set_to_cmp, fill = brewer.pal(length(set_to_cmp),"Dark2")[1:length(set_to_cmp)],
                     alpha = 0.3, filename = NULL, height = 3000, width = 3000,
                     main=ifelse(length(vd_tit) == 0, "", vd_tit), main.cex=1.5)
  grid.newpage()
  grid.draw(vd)
  invisible(vd)
}

# UPDATED VERSION WHERE my_font ARGUMENT CAN BE USED TO CHANGE FONT OF THE PLOT
draw_venn_font <- function(...){
  suppressPackageStartupMessages(library(VennDiagram, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
  suppressPackageStartupMessages(library(RColorBrewer, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
  set_to_cmp <- list(...)
  vd_tit <- set_to_cmp[["vd_tit"]]
  set_to_cmp[["vd_tit"]] <- NULL
  
  myfont <- set_to_cmp[["my_font"]]
  set_to_cmp[["my_font"]] <- NULL


  mycex <- set_to_cmp[["my_cex"]]
  set_to_cmp[["my_cex"]] <- NULL    
  
  myfont <- ifelse(is.null(myfont), "serif", myfont)

  mycex <- ifelse(is.null(mycex), 1, as.numeric(mycex))
  
  vd <- venn.diagram(set_to_cmp, fill = brewer.pal(length(set_to_cmp),"Dark2")[1:length(set_to_cmp)],
                     alpha = 0.3, filename = NULL, height = 3000, width = 3000,
                     main=ifelse(length(vd_tit) == 0, "", vd_tit), main.cex=mycex*1.5, 
                     cat.cex = mycex, cex = mycex,
                     main.fontfamily = myfont)
  grid.newpage()
  grid.draw(vd)
  invisible(vd)
}

######################################################################################################################################################################################################
######################################################################################################################################################################################################
######################################################################################################################################################################################################


plot_ecdf <- function(observ_vect, permut_DT, lineObsCol = "black", halfOnly = F,
                      pointObsCol = rgb(160/255,160/255,160/255, 0.9), polygonPermutCol =  rgb(0/255,76/255,153/255, 0.3),
                      my_stat = "ratioDown") {
  x_val <- seq(ifelse(halfOnly, 0.5, 0),1,0.1)
  observ_vect <- sort(observ_vect, decreasing = T)
  permut_DT <- apply(permut_DT, 2, sort, decreasing=T)
  
  ecdf_obs <- ecdf(observ_vect)(observ_vect)
  ecdf_permut_polygon <- apply(permut_DT, 2, function(x) ecdf(x) (x_val))
  
  plot(ecdf(observ_vect), bty="l", 
       xlab = paste0(my_stat, "ratioDown values"),
       ylab = "cumulative density",
       main = paste0(my_stat, ": cumulative distribution comparison"),
       cex=0.5, col = pointObsCol)
  lines(x = x_val, ecdf(observ_vect)(x_val), col =lineObsCol)
  polygon(x = c(x_val, rev(x_val)), 
          y = c( apply(ecdf_permut_polygon, 1, function(x) min(x)), rev(apply(ecdf_permut_polygon, 1, function(x) max(x)))),
          border=NA,
          col = polygonPermutCol)
  legend(x=1, y= 0, 
         xjust=0.5, yjust=0,
         lty=c(1, NA), pch = c(NA, 15), 
         legend = c("observed", "min-max permut"), 
         pt.cex = c(NA, 2),
         col = c(lineObsCol, polygonPermutCol),
         bty="n")
}


######################################################################################################################################################################################################
######################################################################################################################################################################################################
######################################################################################################################################################################################################

# 
# pointObsCol = "black"
# my_stat = "ratioDown"
# my_main = NULL
# my_ylab = NULL
# my_xlab = NULL
# cexMain = 1
# polygonPermutCol =  rgb(0/255,76/255,153/255, 0.3)
# departureValue = 0.5
# observ_vect=filter_obs_curr_down_half 
# permut_DT=filter_permut_currDown_half
# my_stat = my_stat_curr_ratio
# departureValue = departureFromValue

plot_cumsumDiff05 <- function(observ_vect, permut_DT, pointObsCol = "black", 
           my_stat = "ratioDown",
           my_main = NULL,
           my_ylab = NULL,
           my_xlab = NULL,
           cexMain = 1,
           polygonPermutCol =  rgb(0/255,76/255,153/255, 0.3),
           departureValue = 0.5, drawline=FALSE) {
  
  plotTrueType <- ifelse(drawline, "l", "p")

  if(is.null(my_main)) my_main <- paste0(my_stat, ": cumsum departure from ", departureValue)
  if(is.null(my_ylab)) my_ylab <- paste0("cumsum(abs(", my_stat, " - ", departureValue,"))")
  if(is.null(my_xlab)) my_xlab <- paste0("regions ranked by decreasing ", my_stat)

  observ_vect <- sort(observ_vect, decreasing = T)
  permut_DT <- apply(permut_DT, 2, sort, decreasing=T)
  
  x_val <- c(1:length(observ_vect))
  diff_05_permut <- apply(permut_DT, 2, function(x) cumsum(abs(x-departureValue)))
  
  plot(cumsum(abs(observ_vect - departureValue)) ~ x_val,
       main= my_main,
       cex.main = cexMain,
       type = plotTrueType,
       pch = 16, cex = 0.7,
       xlab= my_xlab, 
       ylab= my_ylab,
       bty="l")
  polygon(x = c(x_val, rev(x_val)), 
          y = c( apply(diff_05_permut, 1, function(x) min(x)), rev(apply(diff_05_permut, 1, function(x) max(x)))),
          border=NA,
          col = polygonPermutCol)
  legend("topleft",
         xjust=0.5, yjust=0,
         pch = c(16, 15), 
         legend = c(paste0("observed (n=", length(observ_vect), ")"), "min-max permut"), 
         pt.cex = c(0.7, 2),
         col = c(pointObsCol, polygonPermutCol),
         bty="n")
}

# to handle different number of observed and permut vect.
plot_cumsumDiff05_revision2 <- function(observ_vect, permut_DT, pointObsCol = "black", 
           my_stat = "ratioDown",
           my_main = NULL,
           my_ylab = NULL,
           my_xlab = NULL,
           cexMain = 1,
           polygonPermutCol =  rgb(0/255,76/255,153/255, 0.3),
           departureValue = 0.5, drawline=FALSE) {
  
  plotTrueType <- ifelse(drawline, "l", "p")

  if(is.null(my_main)) my_main <- paste0(my_stat, ": cumsum departure from ", departureValue)
  if(is.null(my_ylab)) my_ylab <- paste0("cumsum(abs(", my_stat, " - ", departureValue,"))")
  if(is.null(my_xlab)) my_xlab <- paste0("regions ranked by decreasing ", my_stat)

  observ_vect <- sort(observ_vect, decreasing = T)
  permut_DT <- apply(permut_DT, 2, sort, decreasing=T)
  
  x_val <- c(1:length(observ_vect))
  diff_05_permut <- apply(permut_DT, 2, function(x) cumsum(abs(x-departureValue)))


  x_val_max <- max(length(observ_vect), nrow(permut_DT))
  x_val_permut <- 1:nrow(permut_DT)

  y_val_permut <-  c( apply(diff_05_permut, 1, function(x) min(x)), rev(apply(diff_05_permut, 1, function(x) max(x))))
  y_val_obs <- cumsum(abs(observ_vect - departureValue))


y_val_max <- max(c(y_val_permut, y_val_obs))


  stopifnot(!is.na(x_val_max))
  
  plot(y_val_obs ~ x_val,
		
	  xlim = c(0, x_val_max),
	  ylim = c(0, y_val_max),
	
       main= my_main,
       cex.main = cexMain,
       type = plotTrueType,
       pch = 16, cex = 0.7,
       xlab= my_xlab, 
       ylab= my_ylab,
       bty="l")
  polygon(x = c(x_val_permut, rev(x_val_permut)), 
          y = y_val_permut,
          border=NA,
          col = polygonPermutCol)
  legend("topleft",
         xjust=0.5, yjust=0,
         pch = c(16, 15), 
         legend = c(paste0("observed (n=", length(observ_vect), ")"), "min-max permut"), 
         pt.cex = c(0.7, 2),
         col = c(pointObsCol, polygonPermutCol),
         bty="n")
}

######################################################################################################################################################################################################
######################################################################################################################################################################################################
######################################################################################################################################################################################################

plot_cumsumDiff05_withLines <- function(observ_vect, permut_DT, 
           all_quantThresh = NULL,
		   pointObsCol = "black", 
           my_stat = "ratioDown",
           my_main = NULL,
           my_ylab = NULL,
           my_xlab = NULL,
           cexMain = 1,
           polygonPermutCol =  rgb(0/255,76/255,153/255, 0.3),			
		   linePermutCol = "navyblue",
           departureValue = 0.5, drawline=FALSE) {
  
  plotTrueType <- ifelse(drawline, "l", "p")

  if(is.null(my_main)) my_main <- paste0(my_stat, ": cumsum departure from ", departureValue)
  if(is.null(my_ylab)) my_ylab <- paste0("cumsum(abs(", my_stat, " - ", departureValue,"))")
  if(is.null(my_xlab)) my_xlab <- paste0("regions ranked by decreasing ", my_stat)

  observ_vect <- sort(observ_vect, decreasing = T)
  permut_DT <- apply(permut_DT, 2, sort, decreasing=T)
  
  x_val <- c(1:length(observ_vect))
  diff_05_permut <- apply(permut_DT, 2, function(x) cumsum(abs(x-departureValue)))


  observed_yval <- cumsum(abs(observ_vect - departureValue))

  permut_yval_min <- apply(diff_05_permut, 1, function(x) min(x))
  permut_yval_max <- apply(diff_05_permut, 1, function(x) max(x))
  permut_yval_mean <- apply(diff_05_permut, 1, function(x) mean(x))
  stopifnot(!any(is.na(permut_yval_min)))
  stopifnot(!any(is.na(permut_yval_max)))
  stopifnot(!any(is.na(permut_yval_mean)))
  
  plot(observed_yval ~ x_val,
       main= my_main,
       cex.main = cexMain,
       type = plotTrueType,
       pch = 16, cex = 0.7,
       xlab= my_xlab, 
       ylab= my_ylab,
       bty="l")
  polygon(x = c(x_val, rev(x_val)), 
          y = c( permut_yval_min, rev(permut_yval_max)),
          border=NA,
          col = polygonPermutCol)

  # add lines:
  lines(x = c(x_val),
      y = c(permut_yval_min),
      col = linePermutCol)

  lines(x = c(x_val),
      y = c(permut_yval_max),
      col = linePermutCol)

  lines(x = c(x_val),
      y = c(permut_yval_mean),
      col = linePermutCol)

  if(!is.null(all_quantThresh)) {
	stopifnot(is.numeric(all_quantThresh))
    stopifnot(all_quantThresh >= 0 & all_quantThresh <= 1)

    all_auc_quantThresh <- list()

	for(quantThresh in all_quantThresh) {
		permut_yval_currThresh <- apply(diff_05_permut, 1, function(x) quantile(x, probs = quantThresh, names=F) )

	    lines(x = c(x_val),
		      y = c(permut_yval_currThresh),

		      col = linePermutCol)

         aucPermut_currThresh <- auc(x = x_val, y = permut_yval_currThresh)
		 all_auc_quantThresh[[as.character(quantThresh)]] <- aucPermut_currThresh

	}

    legend("topleft",
         xjust=0.5, yjust=0,
         pch = c(16, 15, -1),
         lty = c(-1, -1, 1), 
         legend = c(paste0("observed (n=", length(observ_vect), ")"), "min-max permut", paste0("min-mean-max-quant. thresh: ", paste0(all_quantThresh, collapse=", "))), 
         pt.cex = c(0.7, 2, -1),
         col = c(pointObsCol, polygonPermutCol, linePermutCol),
         bty="n")


    all_auc_quantThresh_vect <- setNames(unlist(all_auc_quantThresh), paste0("quant_", names(all_auc_quantThresh)))

  } else {

    legend("topleft",
         xjust=0.5, yjust=0,
         pch = c(16, 15, -1),
         lty = c(-1, -1, 1), 
         legend = c(paste0("observed (n=", length(observ_vect), ")"), "min-max permut", paste0("min-mean-max")), 
         pt.cex = c(0.7, 2, -1),
         col = c(pointObsCol, polygonPermutCol, linePermutCol),
         bty="n")

	all_auc_quantThresh_vect <- NULL

}


  # compute and return the AUC

  auc_values <- c(
	observed_auc = auc(x = x_val, y = observed_yval),
	minPermut_auc = auc(x = x_val, y = permut_yval_min),
	meanPermut_auc = auc(x = x_val, y = permut_yval_mean),
	maxPermut_auc = auc(x = x_val, y = permut_yval_max)
	)
  auc_values <- c(auc_values, all_auc_quantThresh_vect)

	return(auc_values)

}



######################################################################################################################################################################################################
######################################################################################################################################################################################################
######################################################################################################################################################################################################

# The same as plot_cumsumDiff05, but here the input permut_DT is a list !

plot_cumsumDiff05_list <- function(observ_vect, permut_List, pointObsCol = "black", 
                              my_stat = "ratioDown",
                              my_main = NULL,
                              my_ylab = NULL,
                              my_xlab = NULL,
                              cexMain = 1,
                              polygonPermutCol =  rgb(0/255,76/255,153/255, 0.3),
                              departureValue = 0.5) {
  
  if(is.null(my_main)) my_main <- paste0(my_stat, ": cumsum departure from ", departureValue)
  if(is.null(my_ylab)) my_ylab <- paste0("cumsum(abs(", my_stat, " - ", departureValue,"))")
  if(is.null(my_xlab)) my_xlab <- paste0("regions ranked by decreasing ", my_stat)
  
  observ_vect <- sort(observ_vect, decreasing = T)
  permut_List <- apply(permut_List, 2, sort, decreasing=T)
  
  x_val <- c(1:length(observ_vect))
  diff_05_permut <- apply(permut_List, 2, function(x) cumsum(abs(x-departureValue)))
  
  plot(cumsum(abs(observ_vect - departureValue)) ~ x_val,
       main= my_main,
       cex.main = cexMain,
       pch = 16, cex = 0.7,
       xlab= my_xlab, 
       ylab= my_ylab,
       bty="l")
  polygon(x = c(x_val, rev(x_val)), 
          y = c( apply(diff_05_permut, 1, function(x) min(x)), rev(apply(diff_05_permut, 1, function(x) max(x)))),
          border=NA,
          col = polygonPermutCol)
  legend("topleft",
         xjust=0.5, yjust=0,
         pch = c(16, 15), 
         legend = c(paste0("observed (n=", length(observ_vect), ")"), "min-max permut"), 
         pt.cex = c(0.7, 2),
         col = c(pointObsCol, polygonPermutCol),
         bty="n")
}


######################################################################################################################################################################################################
######################################################################################################################################################################################################
######################################################################################################################################################################################################


# => for each permutation, compute the AUC
# => for the true data, compute the AUC
# => get an empirical p-value for the observed AUC

### for debug:
# pointObsCol = "black" 
# my_main = NULL
# my_ylab = NULL
# my_xlab = NULL
# cexMain = 1
# polygonPermutCol =  rgb(0/255,76/255,153/255, 0.3)
# observ_vect=filter_obs_curr_down
# permut_DT=filter_permut_currDown
# my_stat = my_stat_curr_ratio
# departureValue = departureFromValue

plot_cumsumDiff05_AUC <- function(observ_vect, permut_DT, pointObsCol = "black", 
                                  my_stat = "ratioDown",
                                  my_main = NULL,
                                  my_ylab = NULL,
                                  my_xlab = NULL,
                                  cexMain = 1,
                                  polygonPermutCol =  rgb(0/255,76/255,150/255, 0.3),
                                  departureValue = 0.5) {
  def_mar <- par()$mar
  def_xpd <- par()$xpd
  on.exit( {par(xpd=def_xpd); par(mar=def_mar)})
  
  suppressPackageStartupMessages(library(flux, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) # error bar
  
  
  if(is.null(my_main)) my_main <- paste0(my_stat, ": AUC cumsum departure from ", departureValue)
  if(is.null(my_xlab)) my_xlab <- paste0("cumsum(abs(", my_stat, " - ", departureValue,"))")
  if(is.null(my_ylab)) my_ylab <- paste0("density")
  
  observ_vect <- sort(observ_vect, decreasing = T)
  permut_DT <- apply(permut_DT, 2, sort, decreasing=T)
  x_val <- c(1:length(observ_vect))
  # # for the standard wave plot:
  # diff_05_permut <- apply(permut_DT, 2, function(x) cumsum(abs(x-departureValue)))
  # plot(cumsum(abs(observ_vect - departureValue)) ~ x_val,
  #      main= my_main,
  #      cex.main = cexMain,
  #      pch = 16, cex = 0.7,
  #      xlab= my_xlab, 
  #      ylab= my_ylab,
  #      bty="l")
  # polygon(x = c(x_val, rev(x_val)), 
  #         y = c( apply(diff_05_permut, 1, function(x) min(x)), rev(apply(diff_05_permut, 1, function(x) max(x)))),
  #         border=NA,
  #         col = polygonPermutCol)
  # legend("topleft",
  #        xjust=0.5, yjust=0,
  #        pch = c(16, 15), 
  #        legend = c(paste0("observed (n=", length(observ_vect), ")"), "min-max permut"), 
  #        pt.cex = c(0.7, 2),
  #        col = c(pointObsCol, polygonPermutCol),
  #        bty="n")
  
  cat(".... compute observed AUC \n")
  observed_line <- cumsum(abs(observ_vect - departureValue)) 
  observed_auc <- auc(x = x_val, y = observed_line)
  permut_lines <- apply(permut_DT, 2, function(x) cumsum(abs(x-departureValue)))
  cat(paste0(".... compute permut. AUC [start: ", Sys.time(), "]\t"))
  # permut_auc <- apply(permut_DT, 2, function(perm_line) auc(x=x_val, y = perm_line))
  permut_auc <- foreach(i=1:ncol(permut_lines), .combine='c') %dopar% { auc(x=x_val, y = permut_lines[,i])} 
  cat(paste0("[end: ", Sys.time(), "]\n"))
  stopifnot(length(permut_auc) == ncol(permut_DT))  
  emp_p_val <- (sum(permut_auc >= observed_auc)+1)/(length(permut_auc) + 1)
  # draw the distribution of the permutated AUC and observed AUC
  density_permut <- density(permut_auc)
  
  # foo <- dev.off()
  # guarantee enough space in case observed AUC is not in the range of permut AUC
  par(mar=par()$mar+c(0,0,0,6))
  
  plot(density_permut,
       bty="l",
       # xlim = range(c(permut_auc, observed_auc)),
       main = my_main,
       xlab=my_xlab, ylab=my_ylab)
  
  polygon(x = c(density_permut$x, rev(density_permut$x)),
          y = c( rep(0, length(density_permut$y)), rev(density_permut$y)),
          border=NA,
          col = polygonPermutCol)
  # 
  abline(v = observed_auc, lty=2)
  
  leg_xpos <- ifelse(observed_auc <= max(density_permut$x), observed_auc,  max(density_permut$x))

  if( sum(permut_auc >= observed_auc) == 0) {
   labelTxt_pval <- paste0("(emp. p-val. < ", 1/length(permut_auc), ")")
  } else{
   labelTxt_pval <- paste0("(emp. p-val. = ", formatC(emp_p_val, format = "e", digits = 2), ")")
  }

  par(xpd=TRUE)                   
  text(x = leg_xpos, y = max(density_permut$y), 
       adj= c(0,1),
       labels = paste0("observ. AUC: ", sprintf("%.2f", observed_auc), "\n",
                       # "(emp. p-val. = ", sprintf("%.2f", emp_p_val), ")"))
                       labelTxt_pval))

  
}

######################################################################################################################################################################################################
######################################################################################################################################################################################################
######################################################################################################################################################################################################


# => for each permutation, compute the AUC
# => for the true data, compute the AUC
# => get an empirical p-value for the observed AUC

### for debug:
# pointObsCol = "black" 
# my_main = NULL
# my_ylab = NULL
# my_xlab = NULL
# cexMain = 1
# polygonPermutCol =  rgb(0/255,76/255,153/255, 0.3)
# observ_vect=filter_obs_curr_down
# permut_DT=filter_permut_currDown
# my_stat = my_stat_curr_ratio
# departureValue = departureFromValue

#..._AUC():
#  leg_xpos <- ifelse(observed_auc <= max(density_permut$x), observed_auc,  max(density_permut$x))

#..._AUC2():
#  leg_xpos <- observed_auc
plot_cumsumDiff05_AUC2 <- function(observ_vect, permut_DT, pointObsCol = "black", 
                                  my_stat = "ratioDown",
                                  my_main = NULL,
                                  my_ylab = NULL,
                                  my_xlab = NULL,
                                  cexMain = 1,
                                  polygonPermutCol =  rgb(0/255,76/255,150/255, 0.3),
                                  departureValue = 0.5,
                                  toPlot=TRUE) {
  
  if(toPlot){
    def_mar <- par()$mar
    def_xpd <- par()$xpd
    on.exit( {par(xpd=def_xpd); par(mar=def_mar)})
    if(is.null(my_main)) my_main <- paste0(my_stat, ": AUC cumsum departure from ", departureValue)
    if(is.null(my_xlab)) my_xlab <- paste0("cumsum(abs(", my_stat, " - ", departureValue,"))")
    if(is.null(my_ylab)) my_ylab <- paste0("density")
  }
  suppressPackageStartupMessages(library(flux, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) # error bar
  
  observ_vect <- sort(observ_vect, decreasing = T)
  permut_DT <- apply(permut_DT, 2, sort, decreasing=T)
  x_val <- c(1:length(observ_vect))
  # # for the standard wave plot:
  # diff_05_permut <- apply(permut_DT, 2, function(x) cumsum(abs(x-departureValue)))
  # plot(cumsum(abs(observ_vect - departureValue)) ~ x_val,
  #      main= my_main,
  #      cex.main = cexMain,
  #      pch = 16, cex = 0.7,
  #      xlab= my_xlab, 
  #      ylab= my_ylab,
  #      bty="l")
  # polygon(x = c(x_val, rev(x_val)), 
  #         y = c( apply(diff_05_permut, 1, function(x) min(x)), rev(apply(diff_05_permut, 1, function(x) max(x)))),
  #         border=NA,
  #         col = polygonPermutCol)
  # legend("topleft",
  #        xjust=0.5, yjust=0,
  #        pch = c(16, 15), 
  #        legend = c(paste0("observed (n=", length(observ_vect), ")"), "min-max permut"), 
  #        pt.cex = c(0.7, 2),
  #        col = c(pointObsCol, polygonPermutCol),
  #        bty="n")
  
  cat(".... compute observed AUC \n")
  observed_line <- cumsum(abs(observ_vect - departureValue)) 
  observed_auc <- auc(x = x_val, y = observed_line)
  permut_lines <- apply(permut_DT, 2, function(x) cumsum(abs(x-departureValue)))
  cat(paste0(".... compute permut. AUC [start: ", Sys.time(), "]\t"))
  # permut_auc <- apply(permut_DT, 2, function(perm_line) auc(x=x_val, y = perm_line))
  permut_auc <- foreach(i=1:ncol(permut_lines), .combine='c') %dopar% { auc(x=x_val, y = permut_lines[,i])} 
  cat(paste0("[end: ", Sys.time(), "]\n"))
  stopifnot(length(permut_auc) == ncol(permut_DT))  
  emp_p_val <- (sum(permut_auc >= observed_auc)+1)/(length(permut_auc) + 1)
  if(!toPlot) return(emp_p_val)
  # draw the distribution of the permutated AUC and observed AUC
  density_permut <- density(permut_auc)
  
  # foo <- dev.off()
  # guarantee enough space in case observed AUC is not in the range of permut AUC
  par(mar=par()$mar+c(0,0,0,6))
  
  plot(density_permut,
       bty="l",
       xlim = range(c(permut_auc, observed_auc)),
       main = my_main,
       xlab=my_xlab, ylab=my_ylab)
  
  polygon(x = c(density_permut$x, rev(density_permut$x)),
          y = c( rep(0, length(density_permut$y)), rev(density_permut$y)),
          border=NA,
          col = polygonPermutCol)
  # 
  abline(v = observed_auc, lty=2)
  # leg_xpos <- ifelse(observed_auc <= max(density_permut$x), observed_auc,  max(density_permut$x))
  leg_xpos <- observed_auc

  if( sum(permut_auc >= observed_auc) == 0) {
   labelTxt_pval <- paste0("(emp. p-val. < ", 1/length(permut_auc), ")")
  } else{
   labelTxt_pval <- paste0("(emp. p-val. = ", formatC(emp_p_val, format = "e", digits = 2), ")")
  }

  par(xpd=TRUE)                   
  text(x = leg_xpos, y = max(density_permut$y), 
       adj= c(0,1),
       labels = paste0("observ. AUC: ", sprintf("%.2f", observed_auc), "\n",
                       # "(emp. p-val. = ", sprintf("%.2f", emp_p_val), ")"))
                       labelTxt_pval))

  
}

########################################################################################################
########################################################################################################
########################################################################################################

plot_cumsumDiff05_AUC2_list <- function(observ_vect_cumsum, permut_List_cumsum, pointObsCol = "black", 
                                        my_stat = "ratioDown",
                                        my_main = NULL,
                                        my_ylab = NULL,
                                        my_xlab = NULL,
                                        cexMain = 1,
                                        polygonPermutCol =  rgb(0/255,76/255,150/255, 0.3),
                                        departureValue = 0.5,
                                        toPlot=TRUE) {
  if(toPlot){
    def_mar <- par()$mar
    def_xpd <- par()$xpd
    on.exit( {par(xpd=def_xpd); par(mar=def_mar)})
    if(is.null(my_main)) my_main <- paste0(my_stat, ": AUC cumsum departure from ", departureValue)
    if(is.null(my_xlab)) my_xlab <- paste0("cumsum(abs(", my_stat, " - ", departureValue,"))")
    if(is.null(my_ylab)) my_ylab <- paste0("density")
  }
  suppressPackageStartupMessages(library(flux, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) # error bar
  
  cat(".... compute observed AUC \n")
  x_val <- c(1:length(observ_vect_cumsum))
  
  # when working with shuffle TAD, input already cumsum !
  observed_line <- observ_vect_cumsum 
  observed_auc <- auc(x = x_val, y = observed_line)
  
  cat(paste0(".... compute permut. AUC [start: ", Sys.time(), "]\t"))
  
  permut_auc <- foreach(i=1:length(permut_List_cumsum), .combine='c') %dopar% { auc(x=c(1:length( permut_List_cumsum[[i]])), y = permut_List_cumsum[[i]])} 
  
  cat(paste0("[end: ", Sys.time(), "]\n"))
  stopifnot(length(permut_auc) == length(permut_List_cumsum))  
  emp_p_val <- (sum(permut_auc >= observed_auc)+1)/(length(permut_auc) + 1)
  if(!toPlot) return(emp_p_val)
  # draw the distribution of the permutated AUC and observed AUC
  density_permut <- density(permut_auc)
  
  # foo <- dev.off()
  # guarantee enough space in case observed AUC is not in the range of permut AUC
  par(mar=par()$mar+c(0,0,0,6))
  
  plot(density_permut,
       bty="l",
       xlim = range(c(permut_auc, observed_auc)),
       main = my_main,
       xlab=my_xlab, ylab=my_ylab)
  
  polygon(x = c(density_permut$x, rev(density_permut$x)),
          y = c( rep(0, length(density_permut$y)), rev(density_permut$y)),
          border=NA,
          col = polygonPermutCol)
   
  abline(v = observed_auc, lty=2)
  # leg_xpos <- ifelse(observed_auc <= max(density_permut$x), observed_auc,  max(density_permut$x))
  leg_xpos <- observed_auc

  if( sum(permut_auc >= observed_auc) == 0) {
   labelTxt_pval <- paste0("(emp. p-val. < ", 1/length(permut_auc), ")")
  } else{
   labelTxt_pval <- paste0("(emp. p-val. = ", formatC(emp_p_val, format = "e", digits = 2), ")")
  }

  par(xpd=TRUE)                   
  text(x = leg_xpos, y = max(density_permut$y), 
       adj= c(0,1),
       labels = paste0("observ. AUC: ", sprintf("%.2f", observed_auc), "\n",
                       # "(emp. p-val. = ", sprintf("%.2f", emp_p_val), ")"))
                       labelTxt_pval))
  
}
#system.time(  foreach(i=1:length(tmpList), .combine='c') %dopar% { auc(x=c(1:length( tmpList[[i]])), y = tmpList[[i]])} )
## user  system elapsed 
## 0.545   0.506   0.688 
#system.time(  unlist(lapply(tmpList, function(i)  { auc(x=c(1:length(i)), y = i)})))
## user  system elapsed 
## 1.066   1.278   0.635 

########################################################################################################
########################################################################################################
########################################################################################################

plot_cumsumDiff05_logAUC <- function(observ_vect, permut_DT, pointObsCol = "black", 
                                  my_stat = "ratioDown",
                                  my_main = NULL,
                                  my_ylab = NULL,
                                  my_xlab = NULL,
                                  cexMain = 1,
                                  polygonPermutCol =  rgb(0/255,76/255,150/255, 0.3),
                                  departureValue = 0.5) {
  def_mar <- par()$mar
  def_xpd <- par()$xpd
  on.exit( {par(xpd=def_xpd); par(mar=def_mar)})
  
  suppressPackageStartupMessages(library(flux, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) # error bar
  
  
  if(is.null(my_main)) my_main <- paste0(my_stat, ": log10AUC cumsum departure from ", departureValue)
  if(is.null(my_xlab)) my_xlab <- paste0("cumsum(abs(", my_stat, " - ", departureValue,"))")
  if(is.null(my_ylab)) my_ylab <- paste0("density")
  
  observ_vect <- sort(observ_vect, decreasing = T)
  permut_DT <- apply(permut_DT, 2, sort, decreasing=T)
  x_val <- c(1:length(observ_vect))
  
  cat(".... compute observed log10AUC \n")
  observed_line <- cumsum(abs(observ_vect - departureValue)) 
  observed_auc <- log10(auc(x = x_val, y = observed_line))
  permut_lines <- apply(permut_DT, 2, function(x) cumsum(abs(x-departureValue)))
  cat(paste0(".... compute permut. log10AUC [start: ", Sys.time(), "]\t"))
  # permut_auc <- apply(permut_DT, 2, function(perm_line) auc(x=x_val, y = perm_line))
  permut_auc <- foreach(i=1:ncol(permut_lines), .combine='c') %dopar% { log10(auc(x=x_val, y = permut_lines[,i]))} 
  cat(paste0("[end: ", Sys.time(), "]\n"))
  stopifnot(length(permut_auc) == ncol(permut_DT))  
  emp_p_val <- (sum(permut_auc >= observed_auc)+1)/(length(permut_auc) + 1)
  # draw the distribution of the permutated AUC and observed AUC
  density_permut <- density(permut_auc)
  
  # foo <- dev.off()
  # guarantee enough space in case observed AUC is not in the range of permut AUC
  par(mar=par()$mar+c(0,0,0,6))
  
  plot(density_permut,
       bty="l",
       xlim = range(c(permut_auc, observed_auc)),
       main = my_main,
       xlab=my_xlab, ylab=my_ylab)
  
  polygon(x = c(density_permut$x, rev(density_permut$x)),
          y = c( rep(0, length(density_permut$y)), rev(density_permut$y)),
          border=NA,
          col = polygonPermutCol)
  # 
  abline(v = observed_auc, lty=2)
  # leg_xpos <- ifelse(observed_auc <= max(density_permut$x), observed_auc,  max(density_permut$x))  
  leg_xpos <- observed_auc
  par(xpd=TRUE)                   
  text(x = leg_xpos, y = max(density_permut$y), 
       adj= c(0,1),
       labels = paste0("observ. log10AUC: ", sprintf("%.2f", observed_auc), "\n",
                       "(emp. p-val. = ", formatC(emp_p_val, format = "e", digits = 2), ")"))
                       # "(emp. p-val. = ", sprintf("%.2f", emp_p_val), ")"))
  
}


########################################################################################################
########################################################################################################
########################################################################################################
plot_diff_obs_perm_cumsum_curve <- function(obs_vect, perm_vect, 
                                            my_main = "Difference betw. perm. and obs. cumsum curves",
                                            my_xlab = "regions ranked by increasing ratioConcordant values",
                                            my_ylab = "diff. obs. - mean perm. cumsum"){
  stopifnot(length(perm_vect) == length(obs_vect))
  diff_obs_perm <- obs_vect - perm_vect
  plot(diff_obs_perm ~ c(1:length(diff_obs_perm)),
       bty = "l",
       xaxs ="i",yaxs="i",
       type = "b", cex = 0.7, pch=16, main = my_main, xlab=my_xlab, ylab=my_ylab)
}

######################################################################################################################################################################################################
######################################################################################################################################################################################################
######################################################################################################################################################################################################


plot_cumsumDiff05_line_mean_perm_ks <- function(observ_vect, permut_DT, 
           pointObsCol = "black", 
           my_linetype = "l",                 
           my_stat = "ratioDown",
           my_main = NULL,
           my_ylab =  NULL,
           my_xlab = NULL,
           lineObsCol = "black",
           linePermutCol =  "deepskyblue2") {


  if(is.null(my_main)) my_main <- paste0(my_stat, ": cumsum departure from 0.5")
  if(is.null(my_xlab)) my_xlab <- paste0("regions ranked by decreasing ", my_stat)
  if(is.null(my_ylab)) my_ylab <- paste0("cumsum(abs(", my_stat, " - 0.5))")

  stopifnot(length(observ_vect) == nrow(permut_DT))
  
  # observed vect
  observ_vect <- sort(observ_vect, decreasing = T)
  obs_cumsum <- cumsum(abs(observ_vect - 0.5))
  
  # permut vect - v1 first cumsum than means
  permut_DT_sort <- apply(permut_DT, 2, sort, decreasing=T)
  rownames(permut_DT_sort) <- NULL
  permut_DT_cumsum <- apply(permut_DT_sort, 2, function(x) cumsum(abs(x-0.5)))
  permut_cumsum <- rowMeans(permut_DT_cumsum, na.rm=T)
  
  # permut vect - v0 first means than cumsum - WRONG
  # permut_means_vect <- sort( rowMeans(permut_DT, na.rm=T), decreasing = T)
  # permut_cumsum <- cumsum(abs(permut_means_vect - 0.5))
  # stopifnot(all(names(permut_means_vect) %in% names(observ_vect)))
  # stopifnot(all(names(observ_vect) %in% names(permut_means_vect))

  x_val <- c(1:length(observ_vect))
  
  # ks test
  obs_perm_ks_test <- ks.test(obs_cumsum, permut_cumsum)
  
  obs_ecdf <- ecdf(obs_cumsum)
  obs_ecdf_values <- obs_ecdf(obs_cumsum)
  perm_ecdf <- ecdf(permut_cumsum)
  perm_ecdf_values <- perm_ecdf(permut_cumsum)
  
  x <- seq(min(obs_cumsum, permut_cumsum), max(obs_cumsum, permut_cumsum), length.out=length(obs_cumsum))
  # find the point where the difference between the 2 ecdf is the biggest
  x0 <- x[which( abs(obs_ecdf(x) - perm_ecdf(x)) == max(abs(obs_ecdf(x) - perm_ecdf(x))) )][1]
  y0 <- obs_ecdf(x0)
  y1 <- perm_ecdf(x0)
  
  par(mfrow=c(1,3))
  
  ### PLOT 1: CUMSUM(ABS(RATIOCONCORDANT-0.5))
  plot( NULL, xaxs ="i",yaxs="i",
        xlim = c(1, length(x_val)), ylim = range(c(obs_cumsum,permut_cumsum)),
        main= my_main,
        xlab= my_xlab, 
        ylab= my_ylab,
        bty="l")
  lines(x = x_val, y = obs_cumsum, type = my_linetype, pch = 16, cex = 0.7, col = lineObsCol)
  lines(x = x_val, y = permut_cumsum, type = my_linetype, pch = 16, cex = 0.7, col = linePermutCol)
  
  legend("topleft",
         xjust=0.5, yjust=0,
         inset = c(0.05,0),
         pch = c(16, 16), 
         lty = c(1,1),
         legend = c(paste0("observed (n=", length(obs_cumsum), ")"), c(paste0("permut means (n=", length(permut_cumsum), ")"))), 
         pt.cex = c(0.7, 0.7),
         col = c(lineObsCol, linePermutCol),
         bty="n")
  
  ### PLOT 2: THE DIFFERENCE IN CUMSUM CURVES
  plot_diff_obs_perm_cumsum_curve(obs_vect = obs_cumsum, perm_vect = permut_cumsum)
  
  y0 <- x0
  x_obs <- x_val[which(abs(obs_cumsum-y0) == min(abs(obs_cumsum-y0)))]
  x_perm <- x_val[which(abs(permut_cumsum-y0) == min(abs(permut_cumsum-y0)))][1]
  points(c(x_obs, x_perm), c(y0, y0), pch=16, col="red")
  segments(x_obs, y0, x_perm, y0, col="firebrick3", lty="dotted")
  textX <- mean(c(x_obs, x_perm))
  text("max ecdf dist.", x = textX, y = y0, pos =3, col ="red")  
  points(c(1), c(y0), pch=16, col="red")
  segments(0, y0, min(c(x_obs, x_perm)), y0, col="firebrick3", lty="dotted")
  textX <- mean(c(0, min(c(x_obs, x_perm))))
  text(paste0("y = ", round(y0, 2)), x = textX, y = y0, pos =3, col ="red")  
  
  x <- seq(min(obs_cumsum, permut_cumsum), max(obs_cumsum, permut_cumsum), length.out=length(obs_cumsum))
  # find the point where the difference between the 2 ecdf is the biggest
  x0 <- x[which( abs(obs_ecdf(x) - perm_ecdf(x)) == max(abs(obs_ecdf(x) - perm_ecdf(x))) )][1]
  y0 <- obs_ecdf(x0)
  y1 <- perm_ecdf(x0)
  
  #### PLOT 3: DENSITY OF THE CUMSUM WITH KS
  plot(obs_ecdf, 
       # yaxs="i", 
       bty="l",
       col= lineObsCol, cex = 0.2,
       main = paste0("Density function - ", my_main), 
       xlab = my_ylab,
       ylab = "Density")
  lines(perm_ecdf, col= linePermutCol, cex = 0.2)
  points(c(x0, x0), c(y0, y1), pch=16, col="red")
  segments(x0, y0, x0, y1, col="firebrick3", lty="dotted")
  segments(x0, 0, x0, y1, col="firebrick3", lty="dotted")
  text(x=x0, y=y1, labels=paste0("(", round(x0,2), "; ", round(y1,2), ")"), col="red", pos = 4)
  text(x=x0, y=y0, labels=paste0("(", round(x0,2), "; ", round(y0,2), ")"), col="red", pos = 4)
  legTxt <- paste0("Two-sided KS test:\np-val = ", round(obs_perm_ks_test$p.value, 4), "\nD = ", round(as.numeric(obs_perm_ks_test$statistic), 4))
  legend("topleft", 
         legend = c(paste0("Two-sided KS test:"), paste0("p-val = ", round(obs_perm_ks_test$p.value, 4)), 
                    paste0("D = ", round(as.numeric(obs_perm_ks_test$statistic), 4))),
         text.font = c(2,1,1),
         inset = c(0,0.05),
          bty="n")
  legend("bottomright",
         inset = c(0,0.05),
         pch = c(16, 16), 
         lty = c(1,1),
         legend = c(paste0("observed (n=", length(obs_cumsum), ")"), c(paste0("permut means (n=", length(permut_cumsum), ")"))), 
         pt.cex = c(0.7, 0.7),
         col = c(lineObsCol, linePermutCol),
         bty="n")
  
}


######################################################################################################################################################################################################
######################################################################################################################################################################################################
######################################################################################################################################################################################################


plot_fittedLogReg <- function(observ_vect,  permut_DT,
                             my_stat = "ratioDown",
                             pointObsCol = "coral", fittObsCol = "coral3",
                              pointPermutCol = "lightcyan2", fittPermutCol = "lightcyan4", halfOnly=F) {
  logisticFunc <- function(x, m, b){
    exp(m*x +b)/(1+exp(m*x+b))
  }
  logisticDeriv <- function(x, m, b) {
    (m*exp(m*x+b))/((1+exp(m*x+b))^2)
  }
  logisticGetX <- function(y, m, b) {
    (log(y/(1-y)) - b)/m
  }
  # prepare the OBSERVED data
  observ_vect <- sort(observ_vect, decreasing = T)
  observ_vect_ecdf <- ecdf(observ_vect)(observ_vect)
  # fit the logistic model
  mod_obs <- glm(observ_vect_ecdf ~ observ_vect, family=binomial(link=logit))
  obs_ipt <- as.numeric(coef(mod_obs)["(Intercept)"])
  obs_slope <- as.numeric(coef(mod_obs)[2])
  obs_der05 <- logisticDeriv(0.5, m=obs_slope, b=obs_ipt)
  obs_y05 <- logisticGetX(0.5, m=obs_slope, b=obs_ipt)
  obsTxt <- paste0("m_obs = ", round(obs_slope,2),
                   "\nicpt_obs = ",round(obs_ipt,2), 
                   " \n(y05_obs = ", round(obs_y05,2), ")", 
                    "\nderiv05_obs = ", round(obs_der05,2), "\n")
  # prepare the PERMUT data
  permut_DT <- apply(permut_DT, 2, sort, decreasing=T)
  ecdf_permut_DT <- apply(permut_DT, 2, function(x) ecdf(x) (x))
  permut_vect <- rowMeans(permut_DT)
  permut_vect_ecdf <- rowMeans(ecdf_permut_DT)
  mod_permut <- glm(permut_vect_ecdf ~ permut_vect, family=binomial(link=logit))
  permut_ipt <- as.numeric(coef(mod_permut)["(Intercept)"])
  permut_slope <- as.numeric(coef(mod_permut)[2])
  permut_der05 <- logisticDeriv(0.5, m=permut_slope, b=permut_ipt)
  permut_y05 <- logisticGetX(0.5, m=permut_slope, b=permut_ipt)
  permutTxt <- paste0("m_perm = ", round(permut_slope,2), 
                                "\n icpt_perm = ",round(permut_ipt,2), 
                      "\n(y05_perm = ", round(permut_y05,2), ")",          
                      "\nderiv05_perm = ", round(permut_der05,2))
  # plot the lines and points
  plot(NULL, xlim=c(ifelse(halfOnly, 0.5, 0),1), ylim=c(0,1), bty="l",
        xlab = paste0(my_stat, " values"),
        ylab ="cumulative density", 
        main=paste0(my_stat, ": comparison ECDF Sigmoid fit"))
  points(observ_vect_ecdf ~ observ_vect, 
         pch=16, cex =0.7, col=pointObsCol)
  lines(x = observ_vect, y = fitted(mod_obs), col=fittObsCol)
  points(x = permut_vect, y = permut_vect_ecdf, 
         pch=16, cex =0.7, 
         col=pointPermutCol)
  lines(x = permut_vect, y = fitted(mod_permut), col=fittPermutCol)
  
  legend("topleft", legend = c(obsTxt, permutTxt), 
         text.col = c(fittObsCol, fittPermutCol),
         bty="n")
  legend("bottomright", 
            legend=c("observed", "obs. fit", "mean permut.", "permut. fit"), 
            pch = c(16,NA,16,NA),
            lty=c(NA,1,NA,1),
            col = c(pointObsCol, fittObsCol,pointPermutCol, fittPermutCol),
            bty="n")
}

########################################################################################################
########################################################################################################
########################################################################################################

plot_fittedLogReg_meanPerm <- function(observ_vect,  permut_DT,
                                     my_stat ="ratioDown",
                                       pointObsCol = "coral", fittObsCol = "coral3",
                                       pointPermutCol = "lightcyan2", fittPermutCol = "lightcyan4", halfOnly=F) {
  logisticFunc <- function(x, m, b){
    exp(m*x +b)/(1+exp(m*x+b))
  }
  logisticDeriv <- function(x, m, b) {
    (m*exp(m*x+b))/((1+exp(m*x+b))^2)
  }
  logisticGetX <- function(y, m, b) {
    (log(y/(1-y)) - b)/m
  }

  xseq <- seq(0,1,0.01)

  # prepare the OBSERVED data
  observ_vect <- sort(observ_vect, decreasing = T)
  observ_vect_ecdf <- ecdf(observ_vect)(observ_vect)
  # fit the logistic model
  mod_obs <- glm(observ_vect_ecdf ~ observ_vect, family=binomial(link=logit))
  obs_ipt <- as.numeric(coef(mod_obs)["(Intercept)"])
  obs_slope <- as.numeric(coef(mod_obs)[2])
  obs_der05 <- logisticDeriv(0.5, m=obs_slope, b=obs_ipt)
  obs_y05 <- logisticGetX(0.5, m=obs_slope, b=obs_ipt)
  
  # fit logistic for each permut
  mod_permut_list <- apply(permut_DT, 2, function(i_perm_vec){
    # prepare the OBSERVED data
    i_perm_vec <- sort(i_perm_vec, decreasing = T)
    i_perm_vect_ecdf <- ecdf(i_perm_vec)(i_perm_vec)
    # fit the logistic model
    i_perm_obs <- glm(i_perm_vect_ecdf ~ i_perm_vec, family=binomial(link=logit))
    i_perm_ipt <- as.numeric(coef(i_perm_obs)["(Intercept)"])
    i_perm_slope <- as.numeric(coef(i_perm_obs)[2])
    i_perm_der05 <- logisticDeriv(0.5, m=i_perm_slope, b=i_perm_ipt)
    i_perm_y05 <- logisticGetX(0.5, m=i_perm_slope, b=i_perm_ipt)
    c(perm_ipt = i_perm_ipt, perm_slope= i_perm_slope, 
      perm_der05 = i_perm_der05, perm_y05 = i_perm_y05  )
  })
  stopifnot(ncol(mod_permut_list) == ncol(permut_DT))
  stopifnot(nrow(mod_permut_list) == 4 )
  
  # test emp-pval: test if
  # - m     => obs < perm
  # - ipt  => obs > perm
  # - y05   => obs > perm
  # -der05  => obs < perm
  emp_pval_obs_slope <- (sum(mod_permut_list["perm_slope",] <= obs_slope) + 1)/(ncol(mod_permut_list) + 1)
  emp_pval_obs_ipt <- (sum(mod_permut_list["perm_ipt",] >= obs_ipt) + 1)/(ncol(mod_permut_list) + 1)
  emp_pval_obs_y05 <- (sum(mod_permut_list["perm_y05",] <= obs_y05) + 1)/(ncol(mod_permut_list) + 1)
  emp_pval_obs_der05 <- (sum(mod_permut_list["perm_der05",] <= obs_der05) + 1)/(ncol(mod_permut_list) + 1)
  
  
  obsTxt <- paste0("m_obs = ", round(obs_slope,2), " (emp. p-val: ", round(emp_pval_obs_slope, 4), ")", 
                   "\nicpt_obs = ",round(obs_ipt,2),  " (emp. p-val: ", round(emp_pval_obs_ipt, 4), ")", 
                   " \n(y05_obs = ", round(obs_y05,2), "; emp. p-val: ", round(emp_pval_obs_y05, 4), ")", 
                   "\nderiv05_obs = ", round(obs_der05,2), " (emp. p-val:" , round(emp_pval_obs_der05, 4), ")\n")
  
  
  # prepare the PERMUT data
  permutTxt <- paste0("mean_m_perm = ", round(mean(mod_permut_list["perm_slope",],na.rm=T),2), 
                      "\n mean_icpt_perm = ",round(mean(mod_permut_list["perm_ipt",],na.rm=T),2), 
                      "\n(mean_y05_perm = ", round(mean(mod_permut_list["perm_y05",],na.rm=T),2), ")",          
                      "\nmean_deriv05_perm = ", round(mean(mod_permut_list["perm_der05",],na.rm=T),2))
  # plot the lines and points
  plot(NULL, xlim=c(ifelse(halfOnly, 0.5, 0),1), ylim=c(0,1), bty="l",
       xlab = paste0(my_stat, " values"),
       ylab ="cumulative density", 
       main=paste0(my_stat, ": comparison ECDF Sigmoid fit"))
  points(observ_vect_ecdf ~ observ_vect, 
         pch=16, cex =0.7, col=pointObsCol)
  lines(x = observ_vect, y = fitted(mod_obs), col=fittObsCol)
  
  
  # x_vect <- seq(ifelse(halfOnly, 0.5, 0),1,0.1)
  # meanLogFit <- logisticFunc(x_vect, m=mean(mod_permut_list["perm_slope",],na.rm=T), b=mean(mod_permut_list["perm_ipt",],na.rm=T))
  # v2: take the rowMeans to get the fitted values ??????????????????????????????????????????????????????????????????????????????
   permut_DT_sort_tmp <- apply(permut_DT, 2, function(x){sort(x, decreasing=T)})
   permut_DT_vect <- rowMeans(permut_DT_sort_tmp, na.rm=T)  
   meanLogFit <- logisticFunc(permut_DT_vect, m=mean(mod_permut_list["perm_slope",],na.rm=T), b=mean(mod_permut_list["perm_ipt",],na.rm=T))
  # no KS here, so let v2
  # v3: for the KS, take the same sequence from 0 to 1 for both observed and permut ??????????????????????????????????????????????????????????????????????????????  
  #meanLogFit <- logisticFunc(xseq, m=mean(mod_permut_list["perm_slope",],na.rm=T), b=mean(mod_permut_list["perm_ipt",],na.rm=T))
  
  # points(x = permut_vect, y = permut_vect_ecdf, 
  #        pch=16, cex =0.7, 
  #        col=pointPermutCol)
  lines(x = permut_DT_vect, y = meanLogFit, col=fittPermutCol)
  
  legend("topleft", legend = c(obsTxt, permutTxt), 
         text.col = c(fittObsCol, fittPermutCol),
         bty="n")
  legend("bottomright", 
         legend=c("observed", "obs. fit", "permut. mean fit"), 
         pch = c(16,NA,NA),
         lty=c(NA,1,1),
         col = c(pointObsCol, fittObsCol, fittPermutCol),
         bty="n")
}



########################################################################################################
########################################################################################################
########################################################################################################


plot_diff_obs_perm_ecdf_curve <- function(obs_vect, perm_vect, 
                                          my_stat = "ratioDown",
                                          my_main = NULL,
                                          my_xlab = NULL,
                                          my_ylab = NULL){

  if(is.null(my_main)) my_main <- paste0("Difference betw. perm. and obs. ecdf curves (with x = ", my_stat, " obs. vect. + perm. vect.)")
  if(is.null(my_xlab)) my_xlab <- paste0("regions ranked by increasing ", my_stat, " values")
  if(is.null(my_ylab)) my_ylab <- paste0("diff. obs. ecdf - mean perm. ecdf (on obs. + perm. ", my_stat, " vect)")

  stopifnot(length(perm_vect) == length(obs_vect))
  diff_obs_perm <- obs_vect - perm_vect
  plot(rev(diff_obs_perm) ~ c(1:length(diff_obs_perm)),
       xaxs ="i",yaxs="i",
       bty = "l",
       type = "b", cex = 0.7, pch=16, main = my_main, xlab=my_xlab, ylab=my_ylab)
}



########################################################################################################
########################################################################################################
########################################################################################################

plot_fittedLogReg_fitAll_meanPerm_ks <- function(observ_vect,  
                                                 permut_DT,
                                                 my_stat = "ratioDown",
                                                 my_ylab ="cumulative density",
                                                 my_xlab = NULL,
                                                 my_main = NULL,
                                                 pointObsCol = "coral", 
                                                 fittObsCol = "coral3",
                                                 pointPermutCol = "lightcyan2", 
                                                 fittPermutCol = "lightcyan4", 
                                                 halfOnly=F) {

  if(is.null(my_main)) my_main <- paste0(my_stat, ": comparison ECDF Sigmoid fit")
  if(is.null(my_xlab)) my_xlab <- paste0(my_stat, " values")

  logisticFunc <- function(x, m, b){
    exp(m*x +b)/(1+exp(m*x+b))
  }
  logisticDeriv <- function(x, m, b) {
    (m*exp(m*x+b))/((1+exp(m*x+b))^2)
  }
  logisticGetX <- function(y, m, b) {
    (log(y/(1-y)) - b)/m
  }
  #v4: for the KS, take the same sequence from 0 to 1 for both observed and permut ??????????????????????????????????????????????????????????????????????????????
  xseq <- seq(0, 1, 0.01)

  # prepare the OBSERVED data
  observ_vect <- sort(observ_vect, decreasing = T)
  observ_vect_ecdf <- ecdf(observ_vect)(observ_vect)
  # fit the logistic model
  mod_obs <- glm(observ_vect_ecdf ~ observ_vect, family=binomial(link=logit))
  obs_ipt <- as.numeric(coef(mod_obs)["(Intercept)"])
  obs_slope <- as.numeric(coef(mod_obs)[2])
  obs_der05 <- logisticDeriv(0.5, m=obs_slope, b=obs_ipt)
  obs_y05 <- logisticGetX(0.5, m=obs_slope, b=obs_ipt)
  # v4: for the KS, take the same sequence from 0 to 1 for both observed and permut ??????????????????????????????????????????????????????????????????????????????
  # obs_logFit <- fitted(mod_obs)
  # obs_logFit2 <- logisticFunc(observ_vect, m=obs_slope, b=obs_ipt)
  # all(obs_logFit == obs_logFit2)
  # # TRUE
  
  # fit logistic for each permut
  mod_permut_list <- apply(permut_DT, 2, function(i_perm_vec){
    # prepare the OBSERVED data
    i_perm_vec <- sort(i_perm_vec, decreasing = T)
    i_perm_vect_ecdf <- ecdf(i_perm_vec)(i_perm_vec)
    # fit the logistic model
    i_perm_obs <- glm(i_perm_vect_ecdf ~ i_perm_vec, family=binomial(link=logit))
    i_perm_ipt <- as.numeric(coef(i_perm_obs)["(Intercept)"])
    i_perm_slope <- as.numeric(coef(i_perm_obs)[2])
    i_perm_der05 <- logisticDeriv(0.5, m=i_perm_slope, b=i_perm_ipt)
    i_perm_y05 <- logisticGetX(0.5, m=i_perm_slope, b=i_perm_ipt)
    c(perm_ipt = i_perm_ipt, perm_slope= i_perm_slope, 
      perm_der05 = i_perm_der05, perm_y05 = i_perm_y05  )
  })
  stopifnot(ncol(mod_permut_list) == ncol(permut_DT))
  stopifnot(nrow(mod_permut_list) == 4 )
  
  perm_mean_slope <- mean(mod_permut_list["perm_slope",], na.rm =  T)
  
  # test emp-pval: test if
  # - m     => obs < perm
  # - ipt   => obs > perm
  # - y05   => obs > perm
  # -der05  => obs < perm
  emp_pval_obs_slope <- (sum(mod_permut_list["perm_slope",] <= obs_slope) + 1)/(ncol(mod_permut_list) + 1)
  emp_pval_obs_ipt <- (sum(mod_permut_list["perm_ipt",] >= obs_ipt) + 1)/(ncol(mod_permut_list) + 1)
  emp_pval_obs_y05 <- (sum(mod_permut_list["perm_y05",] <= obs_y05) + 1)/(ncol(mod_permut_list) + 1)
  emp_pval_obs_der05 <- (sum(mod_permut_list["perm_der05",] <= obs_der05) + 1)/(ncol(mod_permut_list) + 1)
  
  # obsTxt <- paste0("m_obs = ", round(obs_slope,2), " (emp. p-val: ", round(emp_pval_obs_slope, 4), ")", 
  #                  "\nicpt_obs = ",round(obs_ipt,2),  " (emp. p-val: ", round(emp_pval_obs_ipt, 4), ")", 
  #                  " \n(y05_obs = ", round(obs_y05,2), "; emp. p-val: ", round(emp_pval_obs_y05, 4), ")", 
  #                  "\nderiv05_obs = ", round(obs_der05,2), " (emp. p-val:" , round(emp_pval_obs_der05, 4), ")\n")
  # 
  # # prepare the PERMUT data
  # permutTxt <- paste0("mean_m_perm = ", round(mean(mod_permut_list["perm_slope",],na.rm=T),2), 
  #                     "\n mean_icpt_perm = ",round(mean(mod_permut_list["perm_ipt",],na.rm=T),2), 
  #                     "\n(mean_y05_perm = ", round(mean(mod_permut_list["perm_y05",],na.rm=T),2), ")",          
  #                     "\nmean_deriv05_perm = ", round(mean(mod_permut_list["perm_der05",],na.rm=T),2))

  # x_vect <- seq(ifelse(halfOnly, 0.5, 0),1,0.1)
  # meanLogFit <- logisticFunc(x_vect, m=mean(mod_permut_list["perm_slope",],na.rm=T), b=mean(mod_permut_list["perm_ipt",],na.rm=T))
  # points(x = permut_vect, y = permut_vect_ecdf, 
  #        pch=16, cex =0.7, 
  #        col=pointPermutCol)
  # lines(x = x_vect, y = meanLogFit, col=fittPermutCol)
  # v2: to have the same number of points, fit the observ_vect
  # perm_meanLogFit <- logisticFunc(observ_vect, m=perm_mean_slope, b=mean(mod_permut_list["perm_ipt",],na.rm=T))

  # v3: take the rowMeans to get the fitted values ??????????????????????????????????????????????????????????????????????????????
  # permut_DT_sort_tmp <- apply(permut_DT, 2, function(x){sort(x, decreasing=T)})
  # permut_DT_vect <- rowMeans(permut_DT_sort_tmp, na.rm=T)  
  # perm_meanLogFit <- logisticFunc(permut_DT_vect, m=perm_mean_slope, b=mean(mod_permut_list["perm_ipt",],na.rm=T))

  #v4: for the KS, take the same sequence from 0 to 1 for both observed and permut ??????????????????????????????????????????????????????????????????????????????
  perm_meanLogFit <- logisticFunc(xseq, m=perm_mean_slope, b=mean(mod_permut_list["perm_ipt",],na.rm=T))
  obs_logFit <- logisticFunc(xseq, m=obs_slope, b=obs_ipt)

  stopifnot(length(perm_meanLogFit) == length(obs_logFit))

  # perm_Func <- function(x) logisticFunc(x, m=perm_mean_slope, b=mean(mod_permut_list["perm_ipt",],na.rm=T))
  # obs_Func <- function(x) logisticFunc(x, m=obs_slope, b=obs_ipt)
  # vect_ks_test <- ks.test(obs_Func, perm_Func)
    
  # KS test between the observed and permut
  vect_ks_test <- ks.test(obs_logFit, perm_meanLogFit)
  obs_ecdf <- ecdf(obs_logFit)
  obs_ecdf_values <- obs_ecdf(obs_logFit)
  perm_ecdf <- ecdf(perm_meanLogFit)
  perm_ecdf_values <- perm_ecdf(perm_meanLogFit)
  
  x <- seq(min(obs_logFit, perm_meanLogFit), max(obs_logFit, perm_meanLogFit), length.out=length(obs_logFit))
  # find the point where the difference between the 2 ecdf is the biggest
  x0 <- x[which( abs(obs_ecdf(x) - perm_ecdf(x)) == max(abs(obs_ecdf(x) - perm_ecdf(x))) )][1]
  y0 <- obs_ecdf(x0)
  y1 <- perm_ecdf(x0)
  
  par(mfrow=c(1,3))
  #### PLOT1:  THE SIGMOID FIT FOR ECDF OF RATIO DOWN VALUES
  # plot the lines and points
  plot(NULL, xlim=c(ifelse(halfOnly, 0.5, 0),1), ylim=c(0,1), bty="l",
       xlab = my_xlab,
       ylab = my_ylab, 
       xaxs ="i",yaxs="i",
       main=my_main)
#  mtext("(1 fit/permut, then mean of curve param. for the KS, with x = observ. values)", cex = 0.9, font=3)
  # v4:
  mtext("(1 fit/permut, then mean of curve parameters\nfor the curves, x = seq(0,1,0.01)", cex = 0.8, font=3, line=-1)
  points(observ_vect_ecdf ~ observ_vect, 
         pch=16, cex =0.7, col=pointObsCol)
  # lines(x = observ_vect, y = obs_logFit , col=fittObsCol)
  # lines(x =observ_vect, y = perm_meanLogFit, col=fittPermutCol)
  #v4: for the KS, take the same sequence from 0 to 1 for both observed and permut ??????????????????????????????????????????????????????????????????????????????
  lines(x = xseq, y = obs_logFit , col=fittObsCol)
  lines(x =xseq, y = perm_meanLogFit, col=fittPermutCol)


  # legend("topleft", legend = c(obsTxt, permutTxt), 
  #        text.col = c(fittObsCol, fittPermutCol),
  #        bty="n")
  legTxt <- paste0("obs. slope/mean perm. slopes\n= ", round(obs_slope,2), "/", round(perm_mean_slope,2), " = ", round(obs_slope/perm_mean_slope, 2))
  legend("topleft", 
         legend = legTxt,
         bty="n")
  legend("bottomright", 
         legend=c(paste0("observed (n=", length(obs_logFit), ")"), "obs. fit", "permut. mean fit"), 
         pch = c(16,NA,NA),
         lty=c(NA,1,1),
         col = c(pointObsCol, fittObsCol, fittPermutCol),
         bty="n")
  
#  y0 <- x0
#  x_obs <- observ_vect[which(abs(obs_logFit-y0) == min(abs(obs_logFit-y0)))][1]
#  x_perm <- observ_vect[which(abs(perm_meanLogFit-y0) == min(abs(perm_meanLogFit-y0)))][1]
#  points(c(x_obs, x_perm), c(y0, y0), pch=16, col="red")
#  segments(x_obs, y0, x_perm, y0, col="firebrick3", lty="dotted")
#  textX <- mean(c(x_obs, x_perm))
#  text("max ecdf dist.", x = textX, y = y0, pos =3, col ="red")  
#  points(c(0), c(y0), pch=16, col="red")
#  segments(0, y0, min(c(x_obs, x_perm)), y0, col="firebrick3", lty="dotted")
#  textX <- mean(c(0, min(c(x_obs, x_perm))))
#  text(paste0("y = ", round(y0, 2)), x = textX, y = y0, pos =3, col ="red")  
  
  
  #### PLOT2:  DISTANCE BETWEEN OBSERVED AND PERMUT ECDF SIGMOID FITS 
  # get the distance between the sigmoid fit
  # plot_diff_obs_perm_ecdf_curve(obs_vect = obs_logFit, perm_vect = perm_meanLogFit)
  # v3: they should have been calculated for the same x value !  ! ???????????????????????????????????????????????????????????????
  # permut_obs_vect <- sort(c(observ_vect, permut_DT_vect))
  # perm_meanLogFit_tmp <- logisticFunc(permut_obs_vect, m=perm_mean_slope, b=mean(mod_permut_list["perm_ipt",],na.rm=T))
  # obs_logFit_tmp <- logisticFunc(permut_obs_vect, m=obs_slope, b=obs_ipt)
  # plot_diff_obs_perm_ecdf_curve(obs_vect = obs_logFit_tmp, perm_vect = perm_meanLogFit_tmp)

  # v4: for the KS, take the same sequence from 0 to 1 for both observed and permut ??????????????????????????????????????????????????????????????????????????????
  perm_meanLogFit_tmp <- logisticFunc(xseq, m=perm_mean_slope, b=mean(mod_permut_list["perm_ipt",],na.rm=T))
  obs_logFit_tmp <- logisticFunc(xseq, m=obs_slope, b=obs_ipt)
  plot_diff_obs_perm_ecdf_curve(obs_vect = obs_logFit_tmp, perm_vect = perm_meanLogFit_tmp, my_main = paste0(my_stat, ": difference betw. perm. and obs. ecdf curves (with x = xseq)"))
    
  #### PLOT3:  DENSITY OF THE ECDF SIGMOID FITS WITH KS TEST RESULTS
  plot(obs_ecdf, 
       # yaxs="i", 
       bty="l",
       col= fittObsCol, cex = 0.2,
       # main = paste0(my_main, "  - density function"), 
       main =paste0(my_stat, ": density function of the ECDF sigmoid fit"),
       xlab = my_ylab,
       ylab = "Density")
  lines(perm_ecdf, col= fittPermutCol, cex = 0.2)
  points(c(x0, x0), c(y0, y1), pch=16, col="red")
  segments(x0, y0, x0, y1, col="firebrick3", lty="dotted")
  segments(x0, 0, x0, y1, col="firebrick3", lty="dotted")
  text(x=x0, y=y1, labels=paste0("(", round(x0,2), "; ", round(y1,2), ")"), col="red", pos = 4)
  text(x=x0, y=y0, labels=paste0("(", round(x0,2), "; ", round(y0,2), ")"), col="red", pos = 4)
  legTxt <- paste0("Two-sided KS test:\np-val = ", round(vect_ks_test$p.value, 4), "\nD = ", round(as.numeric(vect_ks_test$statistic), 4))
  legend("topleft", 
         legend = c(paste0("Two-sided KS test:"), paste0("p-val = ", round(vect_ks_test$p.value, 4)), 
                    paste0("D = ", round(as.numeric(vect_ks_test$statistic), 4))),
         text.font = c(2,1,1),
         inset = c(0,0.05),
         bty="n")
  legend("bottomright",
         inset = c(0,0.05),
         pch = c(16, 16), 
         lty = c(1,1),
         legend = c(paste0("observed (n=", length(obs_logFit), ")"), c(paste0("permut means (n=", length(perm_meanLogFit), ")"))), 
         pt.cex = c(0.7, 0.7),
         col = c(fittObsCol, fittPermutCol),
         bty="n")
  
}

########################################################################################################
########################################################################################################
########################################################################################################


plot_fittedLogReg_fitMeanPerm_ks <- function(observ_vect,  
                                             permut_DT,
                                             my_stat = "ratioDown",
                                             my_ylab ="cumulative density",
                                             my_xlab = NULL,
                                             my_main = NULL,
                                             pointObsCol = "coral", 
                                             fittObsCol = "coral3",
                                             pointPermutCol = "lightcyan2", 
                                             fittPermutCol = "lightcyan4", 
                                             halfOnly=F) {
  if(is.null(my_main)) my_main <- paste0(my_stat, ": comparison ECDF Sigmoid fit")
  if(is.null(my_xlab)) my_xlab <- paste0(my_stat, " values")

  logisticFunc <- function(x, m, b){
    exp(m*x +b)/(1+exp(m*x+b))
  }
  logisticDeriv <- function(x, m, b) {
    (m*exp(m*x+b))/((1+exp(m*x+b))^2)
  }
  logisticGetX <- function(y, m, b) {
    (log(y/(1-y)) - b)/m
  }
  #v4: for the KS, take the same sequence from 0 to 1 for both observed and permut ??????????????????????????????????????????????????????????????????????????????
  xseq <- seq(0, 1, 0.01)

  # v4 => to fit the logistic function -> use the available x data, but for the KS, use the sequence

  # prepare the OBSERVED data
  observ_vect <- sort(observ_vect, decreasing = T)
  observ_vect_ecdf <- ecdf(observ_vect)(observ_vect)
  # fit the logistic model
  mod_obs <- glm(observ_vect_ecdf ~ observ_vect, family=binomial(link=logit))
  obs_ipt <- as.numeric(coef(mod_obs)["(Intercept)"])
  obs_slope <- as.numeric(coef(mod_obs)[2])
  obs_der05 <- logisticDeriv(0.5, m=obs_slope, b=obs_ipt)
  obs_y05 <- logisticGetX(0.5, m=obs_slope, b=obs_ipt)
  #v4: for the KS, take the same sequence from 0 to 1 for both observed and permut ??????????????????????????????????????????????????????????????????????????????
  # obs_logFit <- fitted(mod_obs)
  
  # take the means of the permut and fit a single logistic 
  permut_DT_sort <- apply(permut_DT, 2, function(x){sort(x, decreasing=T)})
  permut_DT_vect <- rowMeans(permut_DT_sort, na.rm=T)  
  stopifnot(length(permut_DT_vect) == nrow(permut_DT))
  permut_DT_vect_ecdf <- ecdf(permut_DT_vect)(permut_DT_vect)
  # fit the logistic model
  mean_perm_obs <- glm(permut_DT_vect_ecdf ~ permut_DT_vect, family=binomial(link=logit))
  mean_perm_ipt <- as.numeric(coef(mean_perm_obs)["(Intercept)"])
  mean_perm_slope <- as.numeric(coef(mean_perm_obs)[2])
  mean_perm_der05 <- logisticDeriv(0.5, m=mean_perm_slope, b=mean_perm_ipt)
  mean_perm_y05 <- logisticGetX(0.5, m=mean_perm_slope, b=mean_perm_ipt)

  # v2: to have the same number of points, fit the observ_vect
  # mean_perm_logFit <- logisticFunc(observ_vect, m=mean_perm_slope, b=mean_perm_ipt)
  
  # v3: take the rowMeans to get the fitted values ??????????????????????????????????????????????????????????????????????????????
  #mean_perm_logFit <- logisticFunc(permut_DT_vect, m=mean_perm_slope, b=mean_perm_ipt)

  #v4: for the KS, take the same sequence from 0 to 1 for both observed and permut ??????????????????????????????????????????????????????????????????????????????
  mean_perm_logFit <- logisticFunc(xseq, m=mean_perm_slope, b=mean_perm_ipt)
  obs_logFit <- logisticFunc(xseq, m=obs_slope, b=obs_ipt)


  stopifnot(length(mean_perm_logFit) == length(obs_logFit))
  
  # KS test between the observed and permut
  vect_ks_test <- ks.test(obs_logFit, mean_perm_logFit)
  obs_ecdf <- ecdf(obs_logFit)
  obs_ecdf_values <- obs_ecdf(obs_logFit)
  # plot(obs_ecdf)
  # lines(obs_ecdf_values ~ obs_logFit, col="red")
  perm_ecdf <- ecdf(mean_perm_logFit)  
  perm_ecdf_values <- perm_ecdf(mean_perm_logFit)
  # plot(perm_ecdf)
  # lines(perm_ecdf_values ~ mean_perm_logFit, col="red")
  
  x <- seq(min(obs_logFit, mean_perm_logFit), max(obs_logFit, mean_perm_logFit), length.out=length(obs_logFit))
  # find the point where the difference between the 2 ecdf is the biggest
  x0 <- x[which( abs(obs_ecdf(x) - perm_ecdf(x)) == max(abs(obs_ecdf(x) - perm_ecdf(x))) )][1]
  y0 <- obs_ecdf(x0)
  y1 <- perm_ecdf(x0)
  
  par(mfrow=c(1,3))

  #### PLOT1:  THE SIGMOID FIT FOR ECDF OF RATIO DOWN VALUES
  # plot the lines and points
  plot(NULL, xlim=c(ifelse(halfOnly, 0.5, 0),1), ylim=c(0,1), bty="l",
       xlab = my_xlab,
       ylab = my_ylab, 
       xaxs ="i",yaxs="i",
       main=my_main)
  # mtext("(row means of sorted permut, 1 fit for all permuts\nto get curve param. for the KS, with x = observ. values)", cex = 0.8, font=3, line=-1)
  # v4:
  mtext("(row means of sorted permut, 1 fit for all permuts\nfor the curves, x = seq(0,1,0.01)", cex = 0.8, font=3, line=-1)
  points(observ_vect_ecdf ~ observ_vect, 
         pch=16, cex =0.7, col=pointObsCol)
  # lines(x = observ_vect, y = obs_logFit , col=fittObsCol)
  # lines(x = observ_vect, y = mean_perm_logFit, col=fittPermutCol)
  #v4: for the KS, take the same sequence from 0 to 1 for both observed and permut ??????????????????????????????????????????????????????????????????????????????
  lines(x = xseq, y = obs_logFit , col=fittObsCol)
  lines(x = xseq, y = mean_perm_logFit, col=fittPermutCol)
  # legend("topleft", legend = c(obsTxt, permutTxt), 
  #        text.col = c(fittObsCol, fittPermutCol),
  #        bty="n")
  legTxt <- paste0("obs. slope/mean perm. slopes\n= ", round(obs_slope,2), "/", round(mean_perm_slope,2), " = ", round(obs_slope/mean_perm_slope, 2))
  legend("topleft", 
         legend = legTxt,
         bty="n")
  legend("bottomright", 
         legend=c(paste0("observed (n=", length(obs_logFit), ")"), "obs. fit", "permut. mean fit"), 
         pch = c(16,NA,NA),
         lty=c(NA,1,1),
         col = c(pointObsCol, fittObsCol, fittPermutCol),
         bty="n")
  
  # v4: do not plot this
#  y0 <- x0
#  x_obs <- observ_vect[which(abs(obs_logFit-y0) == min(abs(obs_logFit-y0)))][1]
#  x_perm <- observ_vect[which(abs(mean_perm_logFit-y0) == min(abs(mean_perm_logFit-y0)))][1]
#  points(c(x_obs, x_perm), c(y0, y0), pch=16, col="red")
#  segments(x_obs, y0, x_perm, y0, col="firebrick3", lty="dotted")
#  textX <- mean(c(x_obs, x_perm))
#  text("max ecdf dist.", x = textX, y = y0, pos =3, col ="red")  
#  points(c(0), c(y0), pch=16, col="red")
#  segments(0, y0, min(c(x_obs, x_perm)), y0, col="firebrick3", lty="dotted")
#  textX <- mean(c(0, min(c(x_obs, x_perm))))
#  text(paste0("y = ", round(y0, 2)), x = textX, y = y0, pos =3, col ="red")  
  
  #### PLOT2:  DISTANCE BETWEEN OBSERVED AND PERMUT ECDF SIGMOID FITS 
  # get the distance between the sigmoid fit
  # plot_diff_obs_perm_ecdf_curve(obs_vect = obs_logFit, perm_vect = mean_perm_logFit)
  # change after v3: they should have been calculated for the same x value ! ???????????????????????????????????????????????????????????????
  # permut_obs_vect <- sort(c(observ_vect, permut_DT_vect))
  # mean_perm_logFit_tmp <- logisticFunc(permut_obs_vect, m=mean_perm_slope, b=mean_perm_ipt)
  # obs_logFit_tmp <- logisticFunc(permut_obs_vect, m=obs_slope, b=obs_ipt)
  # plot_diff_obs_perm_ecdf_curve(obs_vect = obs_logFit_tmp, perm_vect = mean_perm_logFit_tmp)

  #v4: for the KS, take the same sequence from 0 to 1 for both observed and permut ??????????????????????????????????????????????????????????????????????????????
  mean_perm_logFit_tmp <- logisticFunc(xseq, m=mean_perm_slope, b=mean_perm_ipt)
  obs_logFit_tmp <- logisticFunc(xseq, m=obs_slope, b=obs_ipt)
  plot_diff_obs_perm_ecdf_curve(obs_vect = obs_logFit_tmp, perm_vect = mean_perm_logFit_tmp, my_main = paste0(my_stat, ": difference betw. perm. and obs. ecdf curves (with x = xseq)"))

    
  #### PLOT3:  DENSITY OF THE ECDF SIGMOID FITS WITH KS TEST RESULTS
  plot(obs_ecdf,
       # yaxs="i",
       bty="l",
       col= fittObsCol, cex = 0.2,
       # main = paste0(my_main, "  - density function"),
       main = paste0("Density function of the ECDF sigmoid fit"),
       xlab = my_ylab,
       ylab = "Density")
  # same as:
  # plot(y=obs_ecdf_values, x= obs_logFit,
  #      # yaxs="i",
  #      bty="l",
  #      col= "blue", cex = 0.2,
  #      main = paste0(my_main, "  - density function"),
  #      xlab = my_ylab,
  #      ylab = "Density")
  # lines(perm_ecdf, col= "red", cex = 0.2)
  # same as
  # plot(perm_ecdf, 
  #      # yaxs="i", 
  #      bty="l",
  #      col= fittObsCol, cex = 0.2,
  #      main = paste0(my_main, "  - density function"), 
  #      xlab = my_ylab,
  #      ylab = "Density")
  # plot(y=perm_ecdf_values, x= perm_meanLogFit,
  #      # yaxs="i",
  #      bty="l",
  #      col= "red", cex = 0.2,
  #      main = paste0(my_main, "  - density function"),
  #      xlab = my_ylab,
  #      ylab = "Density")
  # lines(obs_ecdf, col= "blue", cex = 0.2)
  # xx=perm_ecdf(obs_logFit)
  # xx=perm_ecdf(perm_meanLogFit)
  # yy=obs_ecdf(obs_logFit)
  # xx <- perm_ecdf(mean_perm_logFit)
  lines(perm_ecdf, col= fittPermutCol, cex = 0.2)
  points(c(x0, x0), c(y0, y1), pch=16, col="red")
  segments(x0, y0, x0, y1, col="firebrick3", lty="dotted")
  segments(x0, 0, x0, y1, col="firebrick3", lty="dotted")
  text(x=x0, y=y1, labels=paste0("(", round(x0,2), "; ", round(y1,2), ")"), col="red", pos = 4)
  text(x=x0, y=y0, labels=paste0("(", round(x0,2), "; ", round(y0,2), ")"), col="red", pos = 4)
  legTxt <- paste0("Two-sided KS test:\np-val = ", round(vect_ks_test$p.value, 4), "\nD = ", round(as.numeric(vect_ks_test$statistic), 4))
  legend("topleft", 
         legend = c(paste0("Two-sided KS test:"), paste0("p-val = ", round(vect_ks_test$p.value, 4)), 
                    paste0("D = ", round(as.numeric(vect_ks_test$statistic), 4))),
         text.font = c(2,1,1),
         inset = c(0,0.05),
         bty="n")
  legend("bottomright",
         inset = c(0,0.05),
         pch = c(16, 16), 
         lty = c(1,1),
         legend = c(paste0("observed (n=", length(obs_logFit), ")"), c(paste0("permut means (n=", length(plot_cumsumDiff05_line_mean_perm_ks), ")"))), 
         pt.cex = c(0.7, 0.7),
         col = c(fittObsCol, fittPermutCol),
         bty="n")

}

########################################################################################################
########################################################################################################
########################################################################################################


plot_cumsum_down_and_up <- function(ratio_down_vect, downCol = "darkslateblue", upCol = "sienna1", my_stat = "ratioDown/Up", lgTxt = c("ratioDown", "ratioUp")) {
  ratio_mostDown_vect <- ratio_down_vect[ratio_down_vect > 0.5]
  ratio_mostUp_vect <- ratio_down_vect[ratio_down_vect < 0.5]
  # convert so that 0.2 of the ratioDown -> 0.8 of the ratioUp
  ratio_mostUp_vect <- abs(ratio_mostUp_vect - 1)
  ratio_mostDown_vect_sort <- sort(ratio_mostDown_vect, decreasing=T)
  ratio_mostUp_vect_sort <- sort(ratio_mostUp_vect, decreasing=T)
  
  ratioDown_cumsum <- cumsum(ratio_mostDown_vect_sort-0.5)
  stopifnot(all(ratioDown_cumsum > 0))
  stopifnot(all(diff(ratioDown_cumsum) > 0))
  ratioUp_cumsum <- cumsum(ratio_mostUp_vect_sort-0.5)
  stopifnot(all(ratioUp_cumsum > 0))
  stopifnot(all(diff(ratioUp_cumsum) > 0))
  
  x_down <- c(1:length(ratioDown_cumsum))
  x_up <- c(1:length(ratioUp_cumsum))
  
  plot( NULL, xaxs ="i",yaxs="i",
        main = paste0("Cumsum ", my_stat, " departure from 0.5"), 
        ylab = paste0("cumsum(", my_stat, " - 0.5)"), 
        xlab = paste0("regions ranked by decreasing ", my_stat),
        xlim = c(1, max(c(x_down, x_up))), 
        ylim = range(c(ratioDown_cumsum,ratioUp_cumsum)),
        bty="l")
  lines(x = x_down, y = ratioDown_cumsum, type = "b", pch = 16, cex = 0.7, col = downCol)
  lines(x = x_up, y = ratioUp_cumsum, type = "b", pch = 16, cex = 0.7, col = upCol)
  legend("bottomright", 
         lty = 1,
         col = c(downCol, upCol),
         legend = lgTxt, 
         text.col = c(downCol, upCol),
         bty="n")
}

########################################################################################################
########################################################################################################
########################################################################################################

calc_aucObsPerm_cumsum <- function(observ_vect, permut_DT, thresh_perm =0.95, 
                     obsLineCol = "black",
                     permutLineCol = "dodgerblue",
                     polygonPermutCol =  rgb(255/255,102/255,102/255, 0.3),
                     departureValue=0.5, 
                     my_stat=NULL, 
                     my_mtext_prefix = NULL,
                     doPlot=T,
                     plotDeviation = FALSE){

  suppressPackageStartupMessages(library(flux, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))  

  if(is.null(my_stat)) {
    main_prefix <- ""
    ylab_suffix <- "stat."
  } else{
    main_prefix <- paste0(my_stat, ": ")
    ylab_suffix <- my_stat
  }
    
  stopifnot(thresh_perm >= 0 & thresh_perm <=1)
  cat(paste0("... cumsum departure from ", departureValue, "\n"))
  cat(paste0("... calc AUC. by taking the ", thresh_perm,  "% quantile of permut values at each rank\n"))
  
  observ_vect <- sort(observ_vect, decreasing = T)
  permut_DT <- apply(permut_DT, 2, sort, decreasing=T)
  stopifnot(nrow(permut_DT) == length(observ_vect))
  
  if( ! all( abs(observ_vect - departureValue) <= 1 &  abs(observ_vect - departureValue) >= 0 )) {
    stop("ERROR: input should be ratio not cumsum !\n")
  }

  cumsum_obs <- cumsum(abs(observ_vect - departureValue))


  if( ! all( abs(permut_DT - departureValue) <= 1 &  abs(permut_DT - departureValue) >= 0 )) {
    stop("ERROR: input should be ratio not cumsum !\n")
  }


  cumsum_permDT <- apply(permut_DT, 2, function(x) cumsum(abs(x-departureValue)))

  
  # take the values of the permutations
  cumsum_perm <- apply(cumsum_permDT, 1 ,function(x) as.numeric(quantile(x, probs=thresh_perm)))
  stopifnot(length(cumsum_perm) == length(observ_vect))
  
  x_val <- c(1:length(observ_vect))

  auc_obs <- auc(x=x_val, y=cumsum_obs)
  auc_perm <- auc(x=x_val, y = cumsum_perm)
  auc_ratio <- auc_obs/auc_perm

  if(doPlot) {
      plot(cumsum_obs ~ x_val,type="l",
           main= paste0(main_prefix, "cumsum departure from ", departureValue),
           lwd=2, col = obsLineCol,
           xlab= paste0("increasing rank by decreasing ", ylab_suffix), 
           ylab= paste0("cumsum departure from ", departureValue),
           bty="l")
      lines(x=x_val, y=cumsum_perm, lwd=2, col = permutLineCol)
      # polygon(x = c(x_val, rev(x_val)), 
      #         y = c( cumsum_perm,
      #                rev(cumsum_obs)),
      #         border=NA,
      #         col = polygonPermutCol)
      my_mtext <- paste0("obs./perm. AUC = ", round(auc_ratio, 2))
      if(!is.null(my_mtext_prefix)) my_mtext <- paste0(my_mtext_prefix, "; ", my_mtext) 
      mtext(my_mtext, cex = 1.2, font=3)#, line=-1)
      
      legend("topleft",
             xjust=0.5, yjust=0,
             # pch = c(16, 15), 
             legend = c(paste0("observed (n=", length(observ_vect), "; AUC=", round(auc_obs,2), ")"), 
                        paste0(thresh_perm * 100, "% permut (AUC=", round(auc_perm,2), ")")), 
             lwd=2,
             lty=c(1,1),
             col = c(obsLineCol, permutLineCol),
             bty="n")
  }
  
  
  
  if(plotDeviation) {
    
    stopifnot(length(cumsum_perm) == length(cumsum_obs) )
    
    diff_ObsPerm <- cumsum_obs - cumsum_perm
    
    plot(diff_ObsPerm ~ x_val,type="l",
         main= paste0(main_prefix, " deviation Observed - Random"),
         lwd=2, col = obsLineCol,
         xlab=paste0("increasing rank by decreasing ", ylab_suffix), 
         ylab= paste0("Deviation observed - ", thresh_perm, "random"),
         bty="l")

    
    my_mtext <- paste0("obs./perm. AUC = ", round(auc_ratio, 2))
    if(!is.null(my_mtext_prefix)) my_mtext <- paste0(my_mtext_prefix, "; ", my_mtext) 
    mtext(my_mtext, cex = 1.2, font=3)#, line=-1)
    
    legend("topleft",
           xjust=0.5, yjust=0,
           # pch = c(16, 15), 
           legend = c(paste0("observed (n=", length(observ_vect), "; AUC=", round(auc_obs,2), ")"), 
                      paste0(thresh_perm * 100, "% permut (AUC=", round(auc_perm,2), ")")), 
           lwd=2,
           lty=c(1,1),
           col = c(obsLineCol, permutLineCol),
           bty="n")
    
  }
  
  return(auc_ratio)
}

########################################################################################################
######################################################################################################## same as previous one but for shuffleTADs (with list instead of DT)
########################################################################################################

calc_aucObsPerm_cumsum_list <- function(observ_vect, shuff_List, 
                    thresh_perm =0.95, 
                     obsLineCol = "black",
                     permutLineCol = "dodgerblue",
                     polygonPermutCol =  rgb(255/255,102/255,102/255, 0.3),
                     departureValue=0.5, 
                     my_stat=NULL, 
                     my_mtext_prefix = NULL,
                     doPlot=TRUE,
                     plotDeviation = FALSE ){

  suppressPackageStartupMessages(library(flux, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))  

  if(is.null(my_stat)) {
    main_prefix <- ""
    ylab_suffix <- "stat."
  } else{
    main_prefix <- paste0(my_stat, ": ")
    ylab_suffix <- my_stat
  }
    
  stopifnot(thresh_perm >= 0 & thresh_perm <=1)
  cat(paste0("... cumsum departure from ", departureValue, "\n"))
  cat(paste0("... calc AUC. by taking the ", thresh_perm,  "% quantile of permut values at each rank\n"))
  
  cat("... prepare observed data\n")
   
  if( ! all( abs(observ_vect - departureValue) <= 1 &  abs(observ_vect - departureValue) >= 0 )) {
    stop("ERROR: input should be ratio not cumsum !\n")
  }

  observ_vect <- sort(observ_vect, decreasing = T)
  cumsum_obs <- cumsum(abs(observ_vect - departureValue))
  x_val_obs <- c(1:length(observ_vect))
  
  cat("... prepare shuffled data\n")
  shuff_List <- lapply(shuff_List, function(x) {
    sort(x, decreasing=TRUE)
  })
  cumsum_shuff_List <- lapply(shuff_List, function(x) {
    x <- sort(x, decreasing=TRUE)
      if( ! all( abs(x - departureValue) <= 1 &  abs(x - departureValue) >= 0 )) {
        stop("ERROR: input should be ratio not cumsum !\n")
      }
    cumsum(abs(x - departureValue))
  })
  # extend the permutation data so that they are all of same length
  maxElements <- max(unlist(lapply(cumsum_shuff_List, length)))  
  stopifnot(!is.na(maxElements))
  cumsum_shuff_List_sameLength <- lapply(cumsum_shuff_List, function(x) x[1:maxElements] <- x[1:maxElements])
  stopifnot(length(unique(sapply(cumsum_shuff_List_sameLength, length))) == 1)
  stopifnot(unique(sapply(cumsum_shuff_List_sameLength, length)) == maxElements)
  # must not have NA as element name -> does not accept NA row names
  # cumsum_shuff_DT <- as.data.frame(cumsum_shuff_List_sameLength)
  cumsum_shuff_DT <- as.data.frame(cumsum_shuff_List_sameLength, 
                    row.names=paste0("region", 1:maxElements), col.names=paste0("shuff", 1:length(cumsum_shuff_List_sameLength)))
  stopifnot(nrow(cumsum_shuff_DT) == maxElements)
  stopifnot(ncol(cumsum_shuff_DT) == nRandomPermutShuffle)
  # na.rm = TRUE here because might have NA when not the same number of TADs
  cumsum_perm <- apply(cumsum_shuff_DT, 1 ,function(x) as.numeric(quantile(x, probs=thresh_perm, na.rm=T)))
  stopifnot(length(cumsum_perm) == maxElements)
  stopifnot(!any(is.na(cumsum_perm)))
  x_val_perm <- 1:maxElements
  
  # compute AUCs
  cat("... compute the AUCs\n")
  auc_obs <- auc(x=x_val_obs, y=cumsum_obs)
  auc_perm <- auc(x=x_val_perm, y = cumsum_perm)
  auc_ratio <- auc_obs/auc_perm

  # compute AUCs
  cat("... compute the AUCs - restrict \n")
  max_idx <- min(c(length(x_val_obs), length(x_val_perm)))
  cat(paste0("...... restrict to: ", max_idx, "/", length(x_val_perm), "\n"))
  new_auc_obs <- auc(x=x_val_obs[1:max_idx], y=cumsum_obs[1:max_idx])
  new_auc_perm <- auc(x=x_val_perm[1:max_idx], y = cumsum_perm[1:max_idx])
  new_auc_ratio <- new_auc_obs/new_auc_perm

  cat(paste0("auc_obs=", round(auc_obs,2), "\tnew_auc_obs", round(new_auc_obs,2), "\n"))
  cat(paste0("auc_perm=", round(auc_perm,2), "\tnew_auc_perm", round(new_auc_perm,2), "\n"))
  cat(paste0("auc_ratio=", round(auc_ratio,2), "\tnew_auc_ratio", round(new_auc_ratio,2), "\n"))
  
  # range of the plot:
  maxTADs <- max(c(maxElements, length(cumsum_obs) ))
  x_range <- c(1,maxTADs)
  y_range <- range(c(unlist(cumsum_shuff_List), cumsum_obs))
  stopifnot(!is.na(y_range))
  stopifnot(!is.na(x_range))
  
  if(doPlot) {

#cat(paste0("x_range=", x_range, "\n"))
#cat(paste0("x_range=", x_range, "\n"))
#cat(paste0("range(x_val_obs)=", range(x_val_obs), "\n"))
#cat(paste0("range(x_val_perm)=", range(x_val_perm), "\n"))
#cat(paste0("range(cumsum_obs)=", range(cumsum_obs), "\n"))
#cat(paste0("range(cumsum_perm)=", range(cumsum_perm), "\n"))
    
    plot(NULL,
         xlim = x_range,
         ylim = y_range,
         main= paste0(main_prefix, "cumsum departure from ", departureValue),
         xlab= paste0("increasing rank by decreasing ", ylab_suffix), 
         ylab= paste0("cumsum departure from ", departureValue),
         bty="l")
    
    lines(x=x_val_obs, y=cumsum_obs, lwd=2, col = obsLineCol)
    lines(x=x_val_perm, y=cumsum_perm, lwd=2, col = permutLineCol)

    my_mtext <- paste0("obs./perm. AUC = ", round(auc_ratio, 2))
    if(!is.null(my_mtext_prefix)) my_mtext <- paste0(my_mtext_prefix, "; ", my_mtext) 
    mtext(my_mtext, cex = 1.2, font=3)#, line=-1)
    
    legend("topleft",
           xjust=0.5, yjust=0,
           # pch = c(16, 15), 
           legend = c(paste0("observed (n=", length(observ_vect), "; AUC=", round(auc_obs,2), ")"), 
                      paste0(thresh_perm * 100, "% permut (AUC=", round(auc_perm,2), ")")), 
           lwd=2,
           lty=c(1,1),
           col = c(obsLineCol, permutLineCol),
           bty="n")
  }


  if(plotDeviation) {
    
    x_ObsPerm <- min(length(cumsum_obs), length(cumsum_perm))
    stopifnot(!is.na(x_ObsPerm))
    stopifnot(x_ObsPerm > 1)
    diff_ObsPerm <- cumsum_obs[1:x_ObsPerm] - cumsum_perm[1:x_ObsPerm]
    
    plot(x= c(1:length(diff_ObsPerm)),
         y = diff_ObsPerm,
         type="l",
         main= paste0(main_prefix, " deviation Observed - Random"),
         xlab= paste0("increasing rank by decreasing ", ylab_suffix), 
         ylab= paste0("Deviation observed - ", thresh_perm, "random"),
         bty="l")
    
    my_mtext <- paste0("obs./perm. AUC = ", round(auc_ratio, 2))
    if(!is.null(my_mtext_prefix)) my_mtext <- paste0(my_mtext_prefix, "; ", my_mtext) 
    mtext(my_mtext, cex = 1.2, font=3)#, line=-1)
    
    legend("topleft",
           xjust=0.5, yjust=0,
           # pch = c(16, 15), 
           legend = c(paste0("observed (n=", length(observ_vect), "; AUC=", round(auc_obs,2), ")"), 
                      paste0(thresh_perm * 100, "% permut (AUC=", round(auc_perm,2), ")")), 
           lwd=2,
           lty=c(1,1),
           col = c(obsLineCol, permutLineCol),
           bty="n")
    
    
  }

  return(auc_ratio)
}

#######################################################################################################################
####################################################################################################################### v2 -> limit the AUC to length of shuffTADs
#######################################################################################################################

calc_aucObsPerm_cumsum_list_v2 <- function(observ_vect, shuff_List, 
                    thresh_perm =0.95, 
                     obsLineCol = "black",
                     permutLineCol = "dodgerblue",
                     polygonPermutCol =  rgb(255/255,102/255,102/255, 0.3),
                     departureValue=0.5, 
                     my_stat=NULL, 
                     my_mtext_prefix = NULL,
                     doPlot=TRUE,
                     plotDeviation = FALSE ){

  suppressPackageStartupMessages(library(flux, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))  

  if(is.null(my_stat)) {
    main_prefix <- ""
    ylab_suffix <- "stat."
  } else{
    main_prefix <- paste0(my_stat, ": ")
    ylab_suffix <- my_stat
  }
    
  stopifnot(thresh_perm >= 0 & thresh_perm <=1)
  cat(paste0("... cumsum departure from ", departureValue, "\n"))
  cat(paste0("... calc AUC. by taking the ", thresh_perm,  "% quantile of permut values at each rank\n"))
  
  cat("... prepare observed data\n")
   
  if( ! all( abs(observ_vect - departureValue) <= 1 &  abs(observ_vect - departureValue) >= 0 )) {
    stop("ERROR: input should be ratio not cumsum !\n")
  }

  observ_vect <- sort(observ_vect, decreasing = T)
  cumsum_obs <- cumsum(abs(observ_vect - departureValue))
  x_val_obs <- c(1:length(observ_vect))
  
  cat("... prepare shuffled data\n")
  shuff_List <- lapply(shuff_List, function(x) {
    sort(x, decreasing=TRUE)
  })
  cumsum_shuff_List <- lapply(shuff_List, function(x) {
    x <- sort(x, decreasing=TRUE)
      if( ! all( abs(x - departureValue) <= 1 &  abs(x - departureValue) >= 0 )) {
        stop("ERROR: input should be ratio not cumsum !\n")
      }
    cumsum(abs(x - departureValue))
  })
  # extend the permutation data so that they are all of same length
  maxElements <- max(unlist(lapply(cumsum_shuff_List, length)))  
  stopifnot(!is.na(maxElements))
  cumsum_shuff_List_sameLength <- lapply(cumsum_shuff_List, function(x) x[1:maxElements] <- x[1:maxElements])
  stopifnot(length(unique(sapply(cumsum_shuff_List_sameLength, length))) == 1)
  stopifnot(unique(sapply(cumsum_shuff_List_sameLength, length)) == maxElements)
  # must not have NA as element name -> does not accept NA row names
  # cumsum_shuff_DT <- as.data.frame(cumsum_shuff_List_sameLength)
  cumsum_shuff_DT <- as.data.frame(cumsum_shuff_List_sameLength, 
                    row.names=paste0("region", 1:maxElements), col.names=paste0("shuff", 1:length(cumsum_shuff_List_sameLength)))
  stopifnot(nrow(cumsum_shuff_DT) == maxElements)
  stopifnot(ncol(cumsum_shuff_DT) == nRandomPermutShuffle)
  # na.rm = TRUE here because might have NA when not the same number of TADs
  cumsum_perm <- apply(cumsum_shuff_DT, 1 ,function(x) as.numeric(quantile(x, probs=thresh_perm, na.rm=T)))
  stopifnot(length(cumsum_perm) == maxElements)
  stopifnot(!any(is.na(cumsum_perm)))
  x_val_perm <- 1:maxElements
  
  # compute AUCs
  cat("... compute the AUCs\n")
  auc_obs <- auc(x=x_val_obs, y=cumsum_obs)
  auc_perm <- auc(x=x_val_perm, y = cumsum_perm)
  auc_ratio <- auc_obs/auc_perm

  # compute AUCs
  cat("... compute the AUCs - restrict \n")
  max_idx <- min(c(length(x_val_obs), length(x_val_perm)))
  cat(paste0("...... restrict to: ", max_idx, "/", length(x_val_perm), "\n"))
  new_auc_obs <- auc(x=x_val_obs[1:max_idx], y=cumsum_obs[1:max_idx])
  new_auc_perm <- auc(x=x_val_perm[1:max_idx], y = cumsum_perm[1:max_idx])
  new_auc_ratio <- new_auc_obs/new_auc_perm

  cat(paste0("auc_obs=", round(auc_obs,2), "\tnew_auc_obs", round(new_auc_obs,2), "\n"))
  cat(paste0("auc_perm=", round(auc_perm,2), "\tnew_auc_perm", round(new_auc_perm,2), "\n"))
  cat(paste0("auc_ratio=", round(auc_ratio,2), "\tnew_auc_ratio", round(new_auc_ratio,2), "\n"))
  

  auc_obs <- new_auc_obs
  auc_perm <- new_auc_perm
  auc_ratio <- new_auc_ratio

  # range of the plot:
  maxTADs <- max(c(maxElements, length(cumsum_obs) ))
  x_range <- c(1,maxTADs)
  y_range <- range(c(unlist(cumsum_shuff_List), cumsum_obs))
  stopifnot(!is.na(y_range))
  stopifnot(!is.na(x_range))
  
  if(doPlot) {

#cat(paste0("x_range=", x_range, "\n"))
#cat(paste0("x_range=", x_range, "\n"))
#cat(paste0("range(x_val_obs)=", range(x_val_obs), "\n"))
#cat(paste0("range(x_val_perm)=", range(x_val_perm), "\n"))
#cat(paste0("range(cumsum_obs)=", range(cumsum_obs), "\n"))
#cat(paste0("range(cumsum_perm)=", range(cumsum_perm), "\n"))
    
    plot(NULL,
         xlim = x_range,
         ylim = y_range,
         main= paste0(main_prefix, "cumsum departure from ", departureValue),
         xlab= paste0("increasing rank by decreasing ", ylab_suffix), 
         ylab= paste0("cumsum departure from ", departureValue),
         bty="l")
    
    lines(x=x_val_obs, y=cumsum_obs, lwd=2, col = obsLineCol)
    lines(x=x_val_perm, y=cumsum_perm, lwd=2, col = permutLineCol)

    my_mtext <- paste0("obs./perm. AUC = ", round(auc_ratio, 2))
    if(!is.null(my_mtext_prefix)) my_mtext <- paste0(my_mtext_prefix, "; ", my_mtext) 
    mtext(my_mtext, cex = 1.2, font=3)#, line=-1)
    
    legend("topleft",
           xjust=0.5, yjust=0,
           # pch = c(16, 15), 
           legend = c(paste0("observed (n=", length(observ_vect), "; AUC=", round(auc_obs,2), ")"), 
                      paste0(thresh_perm * 100, "% permut (AUC=", round(auc_perm,2), ")")), 
           lwd=2,
           lty=c(1,1),
           col = c(obsLineCol, permutLineCol),
           bty="n")
  }


  if(plotDeviation) {
    
    x_ObsPerm <- min(length(cumsum_obs), length(cumsum_perm))
    stopifnot(!is.na(x_ObsPerm))
    stopifnot(x_ObsPerm > 1)
    diff_ObsPerm <- cumsum_obs[1:x_ObsPerm] - cumsum_perm[1:x_ObsPerm]
    
    plot(x= c(1:length(diff_ObsPerm)),
         y = diff_ObsPerm,
         type="l",
         main= paste0(main_prefix, " deviation Observed - Random"),
         xlab= paste0("increasing rank by decreasing ", ylab_suffix), 
         ylab= paste0("Deviation observed - ", thresh_perm, "random"),
         bty="l")
    
    my_mtext <- paste0("obs./perm. AUC = ", round(auc_ratio, 2))
    if(!is.null(my_mtext_prefix)) my_mtext <- paste0(my_mtext_prefix, "; ", my_mtext) 
    mtext(my_mtext, cex = 1.2, font=3)#, line=-1)
    
    legend("topleft",
           xjust=0.5, yjust=0,
           # pch = c(16, 15), 
           legend = c(paste0("observed (n=", length(observ_vect), "; AUC=", round(auc_obs,2), ")"), 
                      paste0(thresh_perm * 100, "% permut (AUC=", round(auc_perm,2), ")")), 
           lwd=2,
           lty=c(1,1),
           col = c(obsLineCol, permutLineCol),
           bty="n")
    
    
  }

  return(auc_ratio)
}

#######################################################################################################################
####################################################################################################################### split table for ezh2
#######################################################################################################################

split_by_rownames <- function(inputDT, splitpattern) {
	singleGene_exprDT <- inputDT[ ! grepl(splitpattern, rownames(inputDT)),]
	manyGenes_exprDT <- inputDT[  grepl(splitpattern, rownames(inputDT)),]

	#  FOR EZH2: processed upstream: no duplicates (but still need to remove then that would induce duplicate rows for a given TAD)
	all_genes_init <- unlist(sapply(rownames(inputDT), function(x) unlist(strsplit(x,split=splitpattern))))
	stopifnot(!any(duplicated(all_genes_init)))

	full_many_exprDT <- foreach(i = 1:dim(manyGenes_exprDT)[1], .combine='rbind') %dopar% {
	  gene_list <- rownames(manyGenes_exprDT)[i]
	  all_genes <-  unlist(strsplit(gene_list,split= splitpattern))
	  tmpDT <- data.frame(lapply(manyGenes_exprDT[i,], rep, length(all_genes)))
	  tmpDT$genes <- all_genes
	  tmpDT
	}
	rownames(full_many_exprDT) <- full_many_exprDT$genes
	full_many_exprDT$genes <- NULL

	full_rnaseqDT <- rbind(singleGene_exprDT, full_many_exprDT)
	stopifnot(nrow(full_rnaseqDT) == length(all_genes_init))

	return(full_rnaseqDT)
}


#######################################################################################################################
####################################################################################################################### assign genes2tad for permutation (step 5b,5c,5d)
#######################################################################################################################



assignGene2TADs <- function(regionDT, geneDT, assignMethod) {
  stopifnot(!any(duplicated(regionDT$region)))
  stopifnot(!any(duplicated(geneDT$entrezID)))
  stopifnot(nrow(regionDT) > 0, nrow(geneDT) > 0)
  stopifnot(assignMethod %in% c("maxOverlap","startPos"))
  stopifnot(c("entrezID", "start","end", "chromo", "strand") %in% colnames(geneDT) )
  stopifnot(c("chromo", "start","end", "region") %in% colnames(regionDT) )

  if(! "package:foreach" %in% search()) suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
  if(! "package:dplyr" %in% search()) suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
  if(! "package:IRanges" %in% search()) suppressPackageStartupMessages(library(IRanges, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))

  # if by startPos -> a bit easier, because I can iterate over TADs (faster)
  if(assignMethod == "startPos") {
    # if on - strand -> start is in fact the end
    geneDT$startPos <- ifelse(geneDT$strand == "-", geneDT$end, geneDT$start)
    gene2tad_dt <- foreach(i=1:nrow(regionDT), .combine="rbind") %do%{
      tad_start <- regionDT$start[i]
      tad_end <- regionDT$end[i]
      tad_chromo <- regionDT$chromo[i]
      # select all the genes that fall within this TAD
      # I CAN DO IN THIS WAY BECAUSE HERE I ASSIGN GENES ACCORDING TO THEIR START POSITION !!!
      tmpDT <- geneDT[geneDT$startPos >= tad_start & geneDT$startPos <= tad_end & geneDT$chromo == tad_chromo,]    
      stopifnot(!any(duplicated(tmpDT$entrezID)))
      if(nrow(tmpDT) == 0) return(NULL)
      tmpDT$region <- regionDT$region[i]
      tmpDT
    }
  } else if(assignMethod == "maxOverlap") {
    all_chr <- intersect(as.character(regionDT$chromo), as.character(geneDT$chromo) )

    
    gene2tad_dt <- foreach(chromo = all_chr, .combine="rbind") %do% {
      
      regionDT_chr <- regionDT[as.character(regionDT$chromo) == chromo,]
      geneDT_chr <- geneDT[as.character(geneDT$chromo) == chromo,]
      
      tad_ranges_chr <- IRanges(start=regionDT_chr$start, end = regionDT_chr$end)
      metadata(tad_ranges_chr)$region <- regionDT_chr$region
      metadata(tad_ranges_chr)$chromo <- regionDT_chr$chromo
      
      gene_ranges_chr <- IRanges(start=geneDT_chr$start, end = geneDT_chr$end)
      metadata(gene_ranges_chr)$entrezID <- as.character(geneDT_chr$entrezID)
      metadata(gene_ranges_chr)$chromo <- geneDT_chr$chromo
      
      stopifnot(all(regionDT[regionDT$chromo == chromo,"start"] == start(tad_ranges_chr)))
      stopifnot(all(regionDT[regionDT$chromo == chromo,"end"] == end(tad_ranges_chr)))
      stopifnot(all(regionDT[regionDT$chromo == chromo,"region"] == metadata(tad_ranges_chr)$region))
      
      stopifnot(all(geneDT[geneDT$chromo == chromo,"start"] == start(gene_ranges_chr)))
      stopifnot(all(geneDT[geneDT$chromo == chromo,"end"] == end(gene_ranges_chr)))
      stopifnot(all(geneDT[geneDT$chromo == chromo,"region"] == metadata(gene_ranges_chr)$region))
      
      # search for hit of genes in TADs
      tad_gene_hits <- findOverlaps(query=gene_ranges_chr, subject= tad_ranges_chr)
      # Compute percent overlap and filter the hits:
      sizeOverlaps <- pintersect(gene_ranges_chr[queryHits(tad_gene_hits)], tad_ranges_chr[subjectHits(tad_gene_hits)])
      sizeOverlapDT <- tad_gene_hits
      # for each hit, get the bp of the overlap
      metadata(sizeOverlapDT)$overlapBp <- width(sizeOverlaps)
      metadata(sizeOverlapDT)$entrezID <- metadata(gene_ranges_chr)$entrezID[queryHits(tad_gene_hits)]
      stopifnot(metadata(gene_ranges_chr)$chromo[queryHits(tad_gene_hits)] == chromo)
      metadata(sizeOverlapDT)$region <- metadata(tad_ranges_chr)$region[subjectHits(tad_gene_hits)]
      
      overlapDT <- data.frame(metadata(sizeOverlapDT), stringsAsFactors = FALSE)
      # overlapBp entrezID    region
      # 43276     2782 chr1_TAD2
      # 62552     2782 chr1_TAD3
      # 314617    23261 chr1_TAD7
      # 669766    23261 chr1_TAD8
      
      
      ## ADDED 21.02 -> not sure !!!
      if(nrow(overlapDT) == 0) {
        cat("!!! WARNING: no gene assigned for", chromo, "!!! \n")
        return(data.frame(entrezID=character(0), region = character(0), stringsAsFactors = FALSE))
      }
      
      # select the hit with the highest overlap
      assignDT <- do.call(rbind,by(overlapDT, overlapDT$entrezID, function(x) x[which.max(as.numeric(as.character(x$overlapBp))), c("entrezID", "region")]))
      rownames(assignDT) <- NULL
      assignDT$entrezID <- as.character(assignDT$entrezID)
      assignDT$region <- as.character(assignDT$region)
      stopifnot(!any(duplicated(assignDT$entrezID)))
      stopifnot(grepl(chromo, assignDT$region))
      stopifnot(length(unique(assignDT$entrezID)) == length(unique(overlapDT$entrezID)))
      
      assignDT
    }
  }
  
  if(nrow(gene2tad_dt) == 0) {
    cat("!!! WARNING: no gene assigned for ALL chromo !!! \n")
    return(data.frame(entrezID=character(0), region = character(0), stringsAsFactors = FALSE))
  }
  
  # gene2tad_dt[gene2tad_dt$entrezID %in% c("2782", "23261"),]
  # entrezID    region
  # 23261 chr1_TAD8
  # 2782 chr1_TAD3
  stopifnot(!any(duplicated(gene2tad_dt$entrezID)))
  # check that there is an overlap  
  tmpDT <- left_join(gene2tad_dt, regionDT,by="region")
  tmpDT <- left_join(tmpDT, geneDT, by="entrezID") 
  stopifnot(!any(duplicated(tmpDT$entrezID)))
  stopifnot(tmpDT$chromo.x == tmpDT$chromo.y)
  stopifnot((tmpDT$end.y <= tmpDT$end.x & tmpDT$end.y >= tmpDT$start.x) |
    (tmpDT$end.x <= tmpDT$end.y & tmpDT$end.x >= tmpDT$start.y))
  
  return(gene2tad_dt[,c("entrezID", "region")])
} 


#######################################################################################################################
####################################################################################################################### empirical FDR
#######################################################################################################################



get_SAM_FDR <- function(obs_vect, permutDT, cut_off, symDir, variableName = "", withPlot=T, plotOffsetY = 0, minQuant = 0.05, maxQuant = 0.95, inputList=FALSE) {
  # permutDT should be the data frame with the rows as being the features and columns the permutations
  if(!inputList)
    permutDT <- permutDT[names(obs_vect),]
  stopifnot(symDir %in% c("symmetric", "higher", "lower"))
  # N = number of real solutions higher than cut-off
  # then, for each permutation (column), how many solutions greater than k ?
  # R = average number of random solutions greater than k
  # get the average number of items called signif from the permut
  # (estimated number of false discoveries = average # of genes called signif from the permut.)
  if(symDir == "symmetric") {
    observ_N <- sum(abs(obs_vect) >= abs(cut_off))  
    if(inputList)
      random_R <- mean(sapply(permutDT, function(x) sum(abs(x) >= abs(cut_off))))
    else 
      random_R <- mean(apply(permutDT, 2, function(x) sum(abs(x) >= abs(cut_off))))
  } else if(symDir == "higher"){
    observ_N <- sum(obs_vect >= cut_off)
    if(inputList)
      random_R <- mean(sapply(permutDT, function(x) sum(x >= cut_off)))
    else
      random_R <- mean(apply(permutDT, 2, function(x) sum(x >= cut_off)))
  } else if (symDir =="lower"){
    observ_N <- sum(obs_vect <= cut_off)
    if(inputList)
      random_R <- mean(sapply(permutDT, function(x) sum(x <= cut_off)))
    else
      random_R <- mean(apply(permutDT, 2, function(x) sum(x <= cut_off)))
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
      quantileShuff <- apply(permutDT, 1,  function(x) quantile(x, probs = c(minQuant, maxQuant)))
      plot(obs_vect, cex=0.7, pch=16, xlab="",
           main = variableName, 
           ylab= paste0(variableName, " observed values"), 
           xlim = c(0, ncol(quantileShuff)), ylim = c(min(obs_vect)-plotOffsetY, max(obs_vect)+plotOffsetY),
           axes=F, bty="l")
      box(bty="l")
      #axis(2, xaxs="i", yaxs="i")
      axis(2)
      #axis(2, at=c(min(obs_vect), max(obs_vect)), labels=F, lwd.ticks=-1, xaxs="r", yaxs="r")
      #axis(2,  labels=F, lwd.ticks=-1, xaxs="r", yaxs="r")
      #axis(1, at=c(0, ncol(quantileShuff)), labels=F, lwd.ticks=-1, xaxs="r", yaxs="r")
      mtext(paste0("empirical FDR = ", round(empFDR, 2), " % (observed signif: ", observ_N, ")"), line=-0.5)
      abline(h = cut_off, lty=2, col = "firebrick3")
      text(x = 10, y = (cut_off + 0.1*cut_off) , paste0("cut-off = ", cut_off),  adj = c(0,0), col = "firebrick3", offset=2)
      if(symDir == "symmetric")
        abline(h = -cut_off, lty=2, col = "firebrick3")
      # polygon(c(x, rev(x)), c(y.high, rev(y.low)),
      polygon(c(rev(1:ncol(quantileShuff)), 1:ncol(quantileShuff)), c(rev(quantileShuff[1,]), quantileShuff[2,]), 
              col = simuColArea, border = NA)
      legend("topleft", legend = paste0(minQuant, "-", maxQuant, " quantile permutations"), col=simuColArea, bty="n", text.col=simuColArea)
      
    }
  }
  invisible(empFDR)
}

#######################################################################################################################
####################################################################################################################### PLOT 2 VECTOR OF POINTS + DISTANCE BETWEEN THEM
#######################################################################################################################
    
    plot_with_departure <- function(vect_obs, vect_perm,col_obs, col_perm,  vect_axis4=NULL, col_axis4 = NULL, lab_axis4="", ...) {
      stopifnot(length(vect_obs) == length(vect_perm))
      vect_x <- 1:length(vect_obs)
      plot(NULL, xlim=range(vect_x, na.rm=T), ylim=range(c(vect_obs, vect_perm), na.rm=T), ...)
      lines(x = vect_x, y = vect_obs, col = col_obs)
      lines(x = vect_x, y = vect_perm, col = col_perm)
      
      if(!is.null(vect_axis4)) {
        stopifnot(length(vect_axis4) == length(vect_x))
        par(new = T)
        plot(NULL, 
             xlim=range(vect_x, na.rm=T),
             ylim = range(vect_axis4, na.rm=T),
             col=col_axis4, 
             pch=16, 
             axes=F, xlab=NA, ylab=NA, bty="l")
        lines(x = vect_x, y = vect_axis4, col=col_axis4)
        box(bty="u")
        axis(side = 4, col=col_axis4)
        mtext(side = 4, line = 3, lab_axis4, col = col_axis4, cex=0.9)
              # paste0('Distance between the 2 curves'), col = "red")
      }
    }

#######################################################################################################################
#######################################################################################################################
#######################################################################################################################

plot_multiDens <- function(size_list, plotTit="", legTxt=NULL, legPos="topright", my_ylab="density", my_xlab="") {
  
  dens <- lapply(size_list, function(x) density(na.omit(x)))
  names(dens) <- names(size_list)
  
  lengthDens <- unlist(lapply(size_list, function(x) length(na.omit(x))))
  
  plot(NA, xlim=range(sapply(dens, "[", "x")), ylim=range(sapply(dens, "[", "y")), 
       main=plotTit, xlab=my_xlab, ylab=my_ylab)
  foo <- mapply(lines, dens, col=1:length(dens))
  if(is.null(legTxt)){
    # legTxt <- names(dens)
    legTxt <- paste0(names(dens), " (n=", lengthDens, ")")
  }
  legend(legPos, legend=legTxt, fill=1:length(dens), bty='n')
}



#######################################################################################################################
#######################################################################################################################
#######################################################################################################################


barplot_funcEnrich_results <- function(resultDT, nTop, xVar, g2t_DT, DE_genes, barcol="slategray", plotHoriz=FALSE, wrap_width = 40){
  
  stopifnot(xVar %in% colnames(resultDT))
  stopifnot(nTop > 0)
  stopifnot(!any(duplicated(resultDT$ID)))
  
  # if table from gseGO  
  if("core_enrichment" %in% colnames(resultDT)) {
    colnames(resultDT)[which(colnames(resultDT) == "core_enrichment")] <- "geneID"
  }
  
  if(xVar == "GeneRatio") {
    resultDT$GeneRatio_nbr <- sapply(resultDT$GeneRatio, function(x) eval(parse(text=x)))  
  } else if(xVar == "BgRatio") {
    resultDT$BgRatio_nbr <- sapply(resultDT$BgRatio, function(x) eval(parse(text=x)))  
  }
  
  if(xVar == "p.adjust" | xVar == "pvalue") {
    resultDT[,xVar] <- -log10(resultDT[,xVar])
    my_ylab <- paste0(xVar, " [-log10]")
  } else{
    my_ylab <- xVar
  }
  
  plotDT <- resultDT[1:nTop, c(xVar, "ID", "Description", "geneID")]
  stopifnot(nrow(plotDT) > 0)
  
  if(plotHoriz) {
    plotDT <- plotDT[nrow(plotDT):1,]
  } 
  
  # ALL THE GENES FROM A GO: retrieve the TADs 
  plotDT$regionID_tot <-  unlist(sapply(1:nrow(plotDT), function(i) {
    all_GO_genes <- strsplit(x= plotDT$geneID[i] , split="/")[[1]]
    stopifnot(all_GO_genes %in% g2t_DT$entrezID)
    paste0(sort(unique(unlist(sapply(all_GO_genes, function(x) as.character(g2t_DT$region[g2t_DT$entrezID == x]))))), collapse = "/")
  }))
  # ALL THE GENES FROM A GO: count the TADs
  plotDT$nbr_TADs_tot <- unlist(sapply(1:nrow(plotDT), function(i) {
    length(strsplit(x= plotDT$regionID_tot[i] , split="/")[[1]])
  }))
  # ALL THE GENES FROM A GO: count the genes
  plotDT$nbr_genes_tot <- unlist(sapply(1:nrow(plotDT), function(i) {
    length(strsplit(x= plotDT$geneID[i] , split="/")[[1]])
  }))
  
  # ONLY DE GENES FROM A GO: retrieve the TADs 
  plotDT$regionID_DE <-  unlist(sapply(1:nrow(plotDT), function(i) {
    all_GO_genes <- strsplit(x= plotDT$geneID[i] , split="/")[[1]]
    all_GO_genes <- all_GO_genes[all_GO_genes %in% DE_genes] 
    stopifnot(all_GO_genes %in% g2t_DT$entrezID)
    paste0(sort(unique(unlist(sapply(all_GO_genes, function(x) as.character(g2t_DT$region[g2t_DT$entrezID == x]))))), collapse = "/")
  }))
  # ONLY DE GENES FROM A GO: count the TADs
  plotDT$nbr_TADs_DE <- unlist(sapply(1:nrow(plotDT), function(i) {
    length(strsplit(x= plotDT$regionID_DE[i] , split="/")[[1]])
  }))
  
  # ONLY DE GENES FROM A GO: COUNT THE GENES
  plotDT$nbr_genes_DE <-  unlist(sapply(1:nrow(plotDT), function(i) {
    all_genes <- strsplit(x= plotDT$geneID[i] , split="/")[[1]]
    sum(all_genes %in% DE_genes)
  }))
  
  
  txtBars <- paste0("# DE genes: ", plotDT$nbr_genes_DE, 
                    "\n(from ", plotDT$nbr_TADs_DE, " TADs)")
  
  plotDT$Description_wrap <- sapply(strwrap(plotDT$Description, width = wrap_width, simplify = FALSE), paste, collapse="\n")
  
  # try(dev.off())
  if(plotHoriz) {
    # par(oma = par()$oma + c(0,8,0,0))
    # par(mfg = par()$mfg + c(0,0,0,20))
    par(mar = par()$mar + c(0,10,0,2))
    my_xlab <- my_ylab
    my_ylab <- ""
    
    barpos <- barplot(plotDT[,xVar], names.arg = plotDT$Description_wrap, 
                      xlab = my_xlab,
                      ylab = my_ylab,
                      cex.names = 0.8,
                      # xlim = c(0, max(plotDT[,xVar])) + c(0,1),
                      xlim = c(0, max(plotDT[,xVar])) + c(0,max(plotDT[,xVar])/10),
                      las=ifelse(plotHoriz, 1, 2),
                      col = barcol,
                      border=F,
                      space=0.1,
                      horiz = plotHoriz,
                      axes=T)
  par(xpd=T)  
  text(y = barpos, x = plotDT[,xVar], pos = 4, offset=0.9, labels = txtBars, cex = 0.6)
  # axis(1, at = c(0, ceiling(max(plotDT[,xVar]))), labels = c("", ceiling(max(plotDT[,xVar]))), tick = T)
  axis(1, at = c(0, max(plotDT[,xVar])), labels = c("", ""), tick = T)
    
  } else{
    # par(oma = par()$oma + c(10,0,0,0))  
    par(mar = par()$mar + c(10,0,0,0))  
    my_xlab <- ""
    
    barpos <- barplot(plotDT[,xVar], names.arg = plotDT$Description_wrap, 
                      xlab = my_xlab,
                      ylab = my_ylab,
                      cex.names = 0.8,
                      # ylim = c(0, max(plotDT[,xVar])) + c(0,1),
                      ylim = c(0, max(plotDT[,xVar])) + c(0,max(plotDT[,xVar])/10),
                      las=ifelse(plotHoriz, 1, 2),
                      col = barcol,
                      border=F,
                      space=0.1,
                      horiz = plotHoriz,
                      axes=T)
    par(xpd=T) 
    text(x = barpos, y = plotDT[,xVar], pos = 3, offset=0.9, labels = txtBars, cex = 0.6)
    # axis(2, at = c(0, ceiling(max(plotDT[,xVar]))), labels = c("", ceiling(max(plotDT[,xVar]))), tick=T)
    axis(2, at = c(0, max(plotDT[,xVar])), labels = c("", ""), tick=T)
  }
    
  
  
}




