
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


# update 16.08.2019: this was done in the child function before
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






if(withExprClass) {

# update 16.08.2019: this was done in the child function before
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


} else {

geneAggregExpression <- NULL
}


# save(geneAggregExpression, file="geneAggregExpression.Rdata", version=2) # used to compare with slower version


  # need to ensure that I get the same order for the genes
  # do the first one
  allT <- get_ShuffledPositions_vFunct(g2TADdt = g2TADdt, geneIDlist = geneIDlist,  #  aggregFun=aggregFun, RNAdt = RNAdt, 
                                      nClass = nClass, TADonly = TADonly, withExprClass = withExprClass, geneAggregExpressionDT = geneAggregExpression )  
  colnames(allT) <- c(colnames(allT)[1], paste0(colnames(allT)[2], "1")) # region1
  genes1 <- allT$entrezID
  if(nSimu >1){
    tmpDT <- foreach(i=2:nSimu, .combine='cbind') %dopar% {
      if(withExprClass) {
        cat(paste0("... WITH CLASS ", aggregFun, " - shuffle: ", i, "/", nSimu, "\n"))
      } else{
        cat(paste0("... NO CLASS - shuffle: ", i, "/", nSimu, "\n"))
      }
    x <- get_ShuffledPositions_vFunct(g2TADdt = g2TADdt,  geneIDlist = geneIDlist,  #  aggregFun=aggregFun, RNAdt = RNAdt, 
                                         nClass = nClass, TADonly = TADonly, withExprClass = withExprClass, geneAggregExpressionDT = geneAggregExpression ) 


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

# UPDATE 16.08.2019 => geneAggregExpression is built only once in the multiShuffled function ! and passed here
# => aggreg expression, computed only once, RNAdt and aggregFun no need to be passed anymore

get_ShuffledPositions_vFunct <- function(g2TADdt, geneIDlist, nClass, TADonly, withExprClass, geneAggregExpressionDT=NULL) { # removed aggregFun and rnaDT


  warning("geneIDlist argument should correspond to rownames of RNAdt")
  warning("duplicated - ambiguous - are removed !!! ")
  duplicatedID <- geneIDlist[duplicated(geneIDlist)]

  if(withExprClass) {
    stopifnot( !is.null(nClass))
  }
  g2TADdt$entrezID <- as.character(g2TADdt$entrezID)
  g2TADdt$region <- as.character(g2TADdt$region)
  # take only the genes for which we have their positions
  # and subset the rnaseq data for these genes only
  if(TADonly) {
    # take only the genes that are in TAD
    geneListTAD <- geneIDlist[geneIDlist  %in% g2TADdt$entrezID[grep("_TAD", g2TADdt$region)] ]
  } else {
    geneListTAD <- geneIDlist[geneIDlist  %in% g2TADdt$entrezID]
  }
  ##########
  ##### DO IT BY SHUFFLING THE LABELS BY CLASS OF EXPRESSION
  ##########
  if(withExprClass) {

    stopifnot(!is.null(geneAggregExpressionDT))     # 16.08.2019 now computed once and passed from parent function
    geneAggregExpression <- geneAggregExpressionDT

    # now, for each class, reshuffle the TAD -> new column with the reshuffled positions
    geneAggregExpression$shuffRegion <- foreach(i_cl = 1:nClass, .combine='c') %do% {
      subDT <- geneAggregExpression[geneAggregExpression$class == i_cl,]
      initPos <- subDT$initRegion
      newPos <- sample(initPos, size=length(initPos), replace=F)  
      stopifnot(all(unique(as.character(initPos)) %in% unique(as.character(newPos))))
      newPos
    }
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

