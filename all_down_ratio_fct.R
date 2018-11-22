## 
# source first, then overwrite some function
source(paste0(setDir, "/mnt/ed4/marie/scripts/TAD_DE_pipeline/TAD_DE_utils.R"))


### DEFINE THE FOLLOWING FUNCTIONS:
#get_statFromShuffle_para
#get_ratioDownByRegion_v2
#get_FCdownByRegion_v2
#get_prodConcordByRegion_v2
#get_prodMeanConcordByRegion_v2
#get_prodLogNbrLogFCByRegion_v2
#get_prodLogRatioNbrByRegion_v2
#get_prodRatioSumByRegion_v2

get_statFromShuffle_para  <- function(DEdt, shuffData, stat_fct, geneIDlist=NULL, ncpu=2, TADonly=F) {
  stopifnot(exists(as.character(stat_fct)))
  suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
  suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
  registerDoMC(cores = ncpu)
  stat_name <- sub("ByRegion_v2", "",  stat_fct)
  stat_name <- sub("get_", "", stat_fct)  
  # ensure always have the same
  regions <- as.character(unique(shuffData[,1]))
  if(TADonly)
    regions <- regions[grep("_TAD", regions)]
  regions <- sort(regions)
  statDT <- foreach(i_perm = 1:ncol(shuffData), .combine='cbind') %dopar% {
	cat(paste0("...... ", stat_name, " permut: ", i_perm, "/", ncol(shuffData), "\n"))
    shuff_g2TAD <- data.frame(entrezID = rownames(shuffData), region = shuffData[,i_perm], stringsAsFactors =F)
    shuff_g2TAD$entrezID <- as.character(shuff_g2TAD$entrezID)
    shuff_g2TAD$region <- as.character(shuff_g2TAD$region)
    x <- do.call(stat_fct, list(g2TADdt=shuff_g2TAD, DEdt = DEdt))
    x[regions]
  }
  stopifnot(all(rownames(statDT) == regions))
  return(statDT)
}

#######################################################################################################################
#######################################################################################################################
#######################################################################################################################

# !!!! version for the random TADs -> I cannot store in a data frame because not always the same set of TADs !!!
# _para_list => the input is a data frame (i.e. always the same set of genes)

get_statFromShuffle_para_list  <- function(DEdt, shuffData, stat_fct, geneIDlist=NULL, ncpu=2, TADonly=F) {
  stopifnot(exists(as.character(stat_fct)))
  suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
  suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
  registerDoMC(cores = ncpu)
  stat_name <- sub("ByRegion_v2", "",  stat_fct)
  stat_name <- sub("get_", "", stat_fct)  
  # ensure always have the same  # => this is not true for the random TADs !!!
  # regions <- as.character(unique(shuffData[,1]))
  # if(TADonly)
  #   regions <- regions[grep("_TAD", regions)]
  # regions <- sort(regions)
  statList <- foreach(i_perm = 1:ncol(shuffData)) %dopar% {
    cat(paste0("...... ", stat_name, " permut: ", i_perm, "/", ncol(shuffData), "\n"))
    shuff_g2TAD <- data.frame(entrezID = rownames(shuffData), region = shuffData[,i_perm], stringsAsFactors =F)
    shuff_g2TAD$entrezID <- as.character(shuff_g2TAD$entrezID)
    shuff_g2TAD$region <- as.character(shuff_g2TAD$region)
    regions <- sort(unique(shuff_g2TAD$region))
    x <- do.call(stat_fct, list(g2TADdt=shuff_g2TAD, DEdt = DEdt))
    x[regions]
  }
  stopifnot(length(statList) == ncol(shuffData))
  # stopifnot(all(rownames(statDT) == regions))# => not true for the random TADs
  return(statList)
}

# _list_para_list => the input is a list (i.e. possibly different set of genes at each permutation)
# ! here the filtering of the DEdt is done internally
get_statFromShuffle_list_para_list  <- function(DEdt, shuffData, stat_fct, geneIDlist=NULL, ncpu=2, TADonly=F) {
  stopifnot(exists(as.character(stat_fct)))
  suppressPackageStartupMessages(library(foreach, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
  suppressPackageStartupMessages(library(doMC, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
  registerDoMC(cores = ncpu)
  stat_name <- sub("ByRegion_v2", "",  stat_fct)
  stat_name <- sub("get_", "", stat_fct)  
  # ensure always have the same  # => this is not true for the random TADs !!!
  # regions <- as.character(unique(shuffData[,1]))
  # if(TADonly)
  #   regions <- regions[grep("_TAD", regions)]
  # regions <- sort(regions)
  statList <- foreach(i_perm = 1:length(shuffData)) %dopar% {
    cat(paste0("...... ", stat_name, " permut: ", i_perm, "/", length(shuffData), "\n"))
    # shuff_g2TAD <- data.frame(entrezID = rownames(shuffData), region = shuffData[,i_perm], stringsAsFactors =F)
    shuffList <- shuffData[[i_perm]]
    shuff_g2TAD <- data.frame(entrezID = names(shuffList), region = shuffList, stringsAsFactors =F)
    shuff_g2TAD$entrezID <- as.character(shuff_g2TAD$entrezID)
    shuff_g2TAD$region <- as.character(shuff_g2TAD$region)
    regions <- sort(unique(shuff_g2TAD$region))
    DEdt_shuff <- DEdt[DEdt$genes %in% shuff_g2TAD$entrezID,]
    stopifnot(!any(duplicated(DEdt$genes)))
    stopifnot(nrow(DEdt_shuff) == nrow(shuff_g2TAD))
    x <- do.call(stat_fct, list(g2TADdt=shuff_g2TAD, DEdt = DEdt))
    x[regions]
  }
  stopifnot(length(statList) == ncol(shuffData))
  # stopifnot(all(rownames(statDT) == regions))# => not true for the random TADs
  return(statList)
}

########################################################################################################################
########################################################################################################################
########################################################################################################################
# rename the function defined in TAD_DE_utils.R
get_ratioDownByRegion_v2 <- get_downByRegion_v2
get_FCdownByRegion_v2 <- get_FCdownByRegion_v2

get_ratioConcordByRegion_v2 <- function(g2TADdt, DEdt) {
   ratioDown <- get_ratioDownByRegion_v2(g2TADdt, DEdt)
   ratioConcord <- abs(ratioDown - 0.5) + 0.5
   stopifnot(all(ratioConcord >= 0 & ratioConcord <= 1))
   return(ratioConcord)
}

get_FCconcordByRegion_v2 <- function(g2TADdt, DEdt) {
   FCdown <- get_FCdownByRegion_v2(g2TADdt, DEdt)
   FCconcord <- abs(FCdown - 0.5) + 0.5
   stopifnot(all(FCconcord >= 0 & FCconcord <= 1))
   return(FCconcord)
}

get_meanRatioFCByRegion_v2 <- function(g2TADdt, DEdt) {
   FCdown <- get_FCdownByRegion_v2(g2TADdt, DEdt)
   ratioDown <- get_ratioDownByRegion_v2(g2TADdt,DEdt)
   stopifnot(all(FCdown >= 0 & FCdown <= 1))
   stopifnot(all(ratioDown >= 0 & ratioDown <= 1))
  stopifnot(all(names(FCdown) == names(ratioDown)))
   meanRatio <- (FCdown+ratioDown)/2
   return(meanRatio)
}

get_meanConcordRatioFCByRegion_v2 <- function(g2TADdt, DEdt) {
   ratioDown <- get_ratioDownByRegion_v2(g2TADdt,DEdt)
   FCdown <- get_FCdownByRegion_v2(g2TADdt, DEdt)
  stopifnot(all(names(FCdown) == names(ratioDown)))
   meanConcord <- unlist(sapply(seq_along(ratioDown), function(i) {
    if(ratioDown[i] > 0.5) {
#      resc_rD <- ratioDown[i] * 2 - 1
      return(0.5*(ratioDown[i] + FCdown[i]))
    } else{
#      resc_rD <- (1-ratioDown[i]) * 2 - 1
      return(0.5*(1-ratioDown[i] + 1-FCdown[i]))
    } 
   }))
   stopifnot(length(ratioDown) == length(FCdown))
   stopifnot(length(ratioDown) == length(meanConcord))
  return(meanConcord)
}


########################################################################################################################
########################################################################################################################
########################################################################################################################
#### COMBINE THE RATIO DOWN IN TERMS OF NUMBER OF GENES AND IN TERMS OF FOLD CHANGE
# if ratioDown (# genes) > 0.5 => ratioDown * FCdown
# else  => (1-ratioDown) * (1-FCdown)
get_prodConcordByRegion_v2 <- function(g2TADdt, DEdt) {
  ratioDown <- get_downByRegion_v2(g2TADdt, DEdt)
  FCdown <- get_FCdownByRegion_v2(g2TADdt, DEdt)
  stopifnot(all(names(ratioDown) == names(FCdown)))
  prodConcord <- foreach(reg = names(ratioDown), .combine='c') %dopar% {
    if(ratioDown[reg] > 0.5) {
      return(ratioDown[reg] * FCdown[reg])
    }else{
      return((1-ratioDown[reg]) * (1-FCdown[reg]))
    }
  }
  stopifnot(all(prodConcord >= 0 & prodConcord <= 1 ))
  return(prodConcord)
}

########################################################################################################################
########################################################################################################################
########################################################################################################################
#### like prod down by region by take the mean
# for each TAD, rescale the FC to the max FC of the TAD (to ensure [0-1] bounding)
# if ratioDown > 0.5 => ratioDown * abs(meanDownFC)
# else => (1-ratioDown) * meanUpFC)
get_prodMeanConcordByRegion_v2 <- function(g2TADdt, DEdt) {
  suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
  g2TADdt$entrezID <- as.character(g2TADdt$entrezID)
  DEdt$genes <- as.character(DEdt$genes)
  stopifnot(all(c("genes", "logFC") %in% colnames(DEdt)))
  stopifnot(all(c("entrezID", "region") %in% colnames(g2TADdt)))
  
  # for each region, rescaled the FC by dividing by the max FC of the region
  # => ensures that the meanFC will be [0-1]
  DE_region <- left_join(DEdt[,c("genes", "logFC")], g2TADdt[,c("entrezID", "region")], by=c("genes"= "entrezID"))
  
  prodMean <- do.call('c', list(by(DE_region, DE_region$region, function(x) {
    rescaledFC <- x$logFC / max(abs(x$logFC), na.rm=T)
    ratioDown <- sum(x$logFC < 0)/length(x$logFC)
    if(ratioDown > 0.5) {
      return(ratioDown * abs(mean(rescaledFC[rescaledFC < 0], na.rm=T)))
    } else {
      return( (1-ratioDown) * mean(rescaledFC[rescaledFC > 0], na.rm=T))
    }
  } )))
  stopifnot(length(prodMean) == length(unique(DE_region$region)))
  stopifnot(all(prodMean >= 0 & prodMean <= 1))
  # stopifnot(all(names(prodMean) == as.character(unique(DE_region$region)))) prodMean returned as sorted
  return(prodMean)
}

########################################################################################################################
########################################################################################################################
########################################################################################################################
# ratioFC  => FC- / FC+
# ratioNbr => nbr- / nbr+
# => log2(ratioFC) * log2(ratioNbr)
get_prodLogNbrLogFCByRegion_v2 <- function(g2TADdt, DEdt) {
  suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
  g2TADdt$entrezID <- as.character(g2TADdt$entrezID)
  DEdt$genes <- as.character(DEdt$genes)
  stopifnot(all(c("genes", "logFC") %in% colnames(DEdt)))
  stopifnot(all(c("entrezID", "region") %in% colnames(g2TADdt)))
  
  # for each region, rescaled the FC by dividing by the max FC of the region
  # => ensures that the meanFC will be [0-1]
  DE_region <- left_join(DEdt[,c("genes", "logFC")], g2TADdt[,c("entrezID", "region")], by=c("genes"= "entrezID"))
  
  prodLogRatio <- do.call('c', list(by(DE_region, DE_region$region, function(x) {
    ratioFC <- sum(abs(x$logFC[x$logFC < 0]))/sum(abs(x$logFC[x$logFC > 0]))
    stopifnot(ratioFC >= 0)
    stopifnot(!is.na(ratioFC))
    ratioNbr <- sum(x$logFC < 0)/sum(x$logFC > 0)
    stopifnot(ratioNbr >= 0)
    return(log2(ratioFC) * log2(ratioNbr))
  } )))
  
  stopifnot(length(prodLogRatio) == length(unique(DE_region$region)))
  
  # stopifnot(all(names(prodMean) == as.character(unique(DE_region$region)))) prodMean returned as sorted
  return(prodLogRatio)
}

########################################################################################################################
########################################################################################################################
########################################################################################################################
# ratioFC  => FC- / (FC+ + FC-) OR 
# ratioNbr => nbr- / (nbr+ + nbr-)
# => if ratioDown > 0.5  1-(-log2(ratioFC) * -log2(ratioNbr))
# => else 1-(-log2(1-ratioFC) * -log2(1-ratioNbr))
## => cannot be bounded with the log2 => remove it
#get_prodLogRatioNbrByRegion_v2 <- function(g2TADdt, DEdt) {
#  suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
#  g2TADdt$entrezID <- as.character(g2TADdt$entrezID)
#  DEdt$genes <- as.character(DEdt$genes)
#  stopifnot(all(c("genes", "logFC") %in% colnames(DEdt)))
#  stopifnot(all(c("entrezID", "region") %in% colnames(g2TADdt)))
#  
#  # for each region, rescaled the FC by dividing by the max FC of the region
#  # => ensures that the meanFC will be [0-1]
#  DE_region <- left_join(DEdt[,c("genes", "logFC")], g2TADdt[,c("entrezID", "region")], by=c("genes"= "entrezID"))
#  
#  prodLogRatio <- do.call('c', list(by(DE_region, DE_region$region, function(x) {
#    ratioFC <- sum(abs(x$logFC[x$logFC < 0]))/(sum(abs(x$logFC[x$logFC > 0])) + sum(abs(x$logFC[x$logFC < 0])))
#    stopifnot(ratioFC >= 0 & ratioFC <= 1)
#    stopifnot(!is.na(ratioFC))
#    ratioNbr <- sum(x$logFC < 0)/(sum(x$logFC > 0) + sum(x$logFC < 0))
#    stopifnot(ratioNbr >= 0 & ratioNbr <= 1)
#    if(ratioNbr > 0.5) {
#      # 02.10 => change log2 to -log2 to bound between 0 and 1, then take 1- ... => so that 0 worst score, 1 best score
#      #return(log2(ratioFC) * log2(ratioNbr))        
#      return( 1-(-log2(ratioFC) * -log2(ratioNbr))  )   
#    } else {
#      #return(log2(1-ratioFC) * log2(1-ratioNbr))      
#      return( 1-(-log2(1-ratioFC) * -log2(1-ratioNbr)) )
#    }
#  } )))
#  
#  stopifnot(length(prodLogRatio) == length(unique(DE_region$region)))
#  
#  # stopifnot(all(names(prodMean) == as.character(unique(DE_region$region)))) prodMean returned as sorted
#  return(prodLogRatio)
#}

########################################################################################################################
########################################################################################################################
########################################################################################################################
 # (ratioDown-0.5) * (sum(FC-) - sum(FC+))
# ?? should I divide the diffSum by tot sum FC to bound it ?

get_prodRatioSumByRegion_v2 <- function(g2TADdt, DEdt) {
  suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
  g2TADdt$entrezID <- as.character(g2TADdt$entrezID)
  DEdt$genes <- as.character(DEdt$genes)
  stopifnot(all(c("genes", "logFC") %in% colnames(DEdt)))
  stopifnot(all(c("entrezID", "region") %in% colnames(g2TADdt)))
  DE_region <- left_join(DEdt[,c("genes", "logFC")], g2TADdt[,c("entrezID", "region")], by=c("genes"= "entrezID"))
  prodRatioSum <- do.call('c', list(by(DE_region, DE_region$region, function(x) {
    diffSum <- sum(abs(x$logFC[x$logFC < 0])) - sum(abs(x$logFC[x$logFC > 0])) 
    ratioDown <- sum(x$logFC < 0)/length(x$logFC)
    # 02.10: to get something bounded, divide diffSum by tot FC and + 0.5
    diffSum <- diffSum/sum(abs(x$logFC))
    return( (ratioDown-0.5) * diffSum + 0.5)
  })))
  stopifnot(length(prodRatioSum) == length(unique(DE_region$region)))
  return(prodRatioSum)
}

########################################################################################################################
########################################################################################################################
########################################################################################################################
#### 
# signed: ranges from -1 to 1
# will be 0 if either FCdown=FCup or ratioDown=ratioUp
# 2* added so that the left part ranges -1 to 1
# negative is FC and ratio uncoordinated
# 2*(ratioDown-0.5) * (FCdown-FCup)/FCtot
# => this is the same as get_prodSignedConcordByRegion_v2()
get_prodSignedRatioByRegion_v2 <- function(g2TADdt, DEdt) {
  suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
  g2TADdt$entrezID <- as.character(g2TADdt$entrezID)
  DEdt$genes <- as.character(DEdt$genes)
  stopifnot(all(c("genes", "logFC") %in% colnames(DEdt)))
  stopifnot(all(c("entrezID", "region") %in% colnames(g2TADdt)))
  
  # for each region, rescaled the FC by dividing by the max FC of the region
  # => ensures that the meanFC will be [0-1]
  DE_region <- left_join(DEdt[,c("genes", "logFC")], g2TADdt[,c("entrezID", "region")], by=c("genes"= "entrezID"))
  
  all_prodRatio <- do.call('c', list(by(DE_region, DE_region$region, function(x) {
    diffRatioFC <- (sum(abs(x$logFC[x$logFC < 0])) - sum(abs(x$logFC[x$logFC > 0])) ) /sum(abs(x$logFC))
    ratioDown <- sum(x$logFC < 0)/length(x$logFC)
    stopifnot(ratioDown >=0 & ratioDown <=1)
    stopifnot(diffRatioFC >=-1 & diffRatioFC <=1)
    prodRatio <- 2*(ratioDown - 0.5) * diffRatioFC
    stopifnot(all(prodRatio >= -1 & prodRatio <= 1))
    return(prodRatio)
    })))
  stopifnot(length(all_prodRatio) == length(unique(DE_region$region)))
  return(all_prodRatio)
}

### since (FC- - FC+)/FC = (FC- - (FC - FC-))/FC = 2FCdown -1, can be rewritten as:
########################################################################################################################
########################################################################################################################
########################################################################################################################
#### 
# signed: ranges from -1 to 1
# will be 0 if either FCdown=FCup or ratioDown=ratioUp

get_prodSignedConcordByRegion_v2 <- function(g2TADdt, DEdt) {
  suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE)) 
  g2TADdt$entrezID <- as.character(g2TADdt$entrezID)
  DEdt$genes <- as.character(DEdt$genes)
  stopifnot(all(c("genes", "logFC") %in% colnames(DEdt)))
  stopifnot(all(c("entrezID", "region") %in% colnames(g2TADdt)))
  ratioDown <- get_ratioDownByRegion_v2(g2TADdt,DEdt)
  ratioDown <- ratioDown * 2 - 1
  FCdown <- get_FCdownByRegion_v2(g2TADdt, DEdt)
  FCdown <- FCdown * 2 - 1
  stopifnot(all(names(FCdown) == names(ratioDown)))
  signedProdConcord <- ratioDown * FCdown
  stopifnot(all(signedProdConcord >= -1 & signedProdConcord <= 1))
  return(signedProdConcord)
}

get_prodSignedRatioByRegion_v2 <- get_prodSignedConcordByRegion_v2

########################################################################################################################
########################################################################################################################
########################################################################################################################
#### 


#DEdt <- DE_topTable
#ggdensity(DEdt, x = "logFC", rug=T)

#DEdt$resc1 <- sign(DEdt$logFC) * rescFunc(abs(DEdt$logFC), newMin = 0, newMax = 1)
#ggdensity(DEdt, x = "resc1", rug=T)
#plot(DEdt$logFC ~ DEdt$resc1, pch=16, cex=0.7)

# DEdt$resc1b <- sign(DEdt$logFC) * abs(rescFunc(DEdt$logFC, newMin = -1, newMax = 1))
# ggdensity(DEdt, x = "resc1b", rug=T) # => does not respect the shape of the initial distribution [bimodal], 2 modes near 0.5
# plot(DEdt$logFC ~ DEdt$resc1b, pch=16, cex=0.7)


#DEdt$resc3 <- sign(DEdt$logFC) * rescFunc(abs(quantNorm(DEdt$logFC)), newMin = 0, newMax = 1)
#ggdensity(DEdt, x = "resc3", rug=T)
#plot(DEdt$logFC ~ DEdt$resc3, pch=16, cex=0.7)

# DEdt$resc5 <- sign(DEdt$logFC) * rescFunc(quantNorm(abs(DEdt$logFC)), newMin = 0, newMax = 1)
# ggdensity(DEdt, x = "resc5", rug=T) # => does not respect the shape of the initial distribution [bimodal], 2 modes near 0.5
# plot(DEdt$logFC ~ DEdt$resc5, pch=16, cex=0.7)



rescFunc <- function(x, newMin, newMax, oldMin=min(x,na.rm=T), oldMax=max(x, na.rm=T)) {
  (newMax-newMin)/(oldMax-oldMin) * (x-oldMin) + newMin
}



get_rescWeightedByRegion_v2 <- function(g2TADdt, DEdt) {
  suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
  g2TADdt$entrezID <- as.character(g2TADdt$entrezID)
  g2TADdt$region <- as.character(g2TADdt$region)
  DEdt$genes <- as.character(DEdt$genes)
  g2TADdt <- g2TADdt[g2TADdt$entrezID %in% DEdt$genes,]
  mergedDT <- left_join(DEdt[,c("genes", "logFC")], g2TADdt[,c("entrezID", "region")], by=c("genes" = "entrezID"))
  mergedDT <- na.omit(mergedDT)
  stopifnot(is.numeric(mergedDT$logFC[1]))
  # rescaled logFC genome-wide
  mergedDT$rescFC <- sign(mergedDT$logFC) * rescFunc(abs(mergedDT$logFC), newMin = 0, newMax = 1)
  weightDownDT <- aggregate(rescFC ~ region, data=mergedDT, function(x) mean(x[x<0], na.rm=T))
  weightDown <- setNames(weightDownDT$rescFC, weightDownDT$region)
  stopifnot(all(na.omit(weightDown) < 0))
  weightUpDT <- aggregate(rescFC ~ region, data=mergedDT, function(x) mean(x[x>0], na.rm=T))
  weightUp <- setNames(weightUpDT$rescFC, weightUpDT$region)
  stopifnot(all(na.omit(weightUp) > 0))
  ratioDownDT <- aggregate(rescFC ~ region, data=mergedDT, function(x) sum(x < 0)/length(x[x!=0]))
  ratioDown <- setNames(ratioDownDT$rescFC, ratioDownDT$region)
  # ratioDown <- get_ratioDownByRegion_v2(g2TADdt,DEdt)
  # there is a problem with this function that does not exclude the ==0,
  # e.g. if one is down and 2 other are 0
  # => the ratioDown will be < 0.5 but weightUp will be NaN
  stopifnot(!any(is.na(ratioDown)))
  stopifnot(all(names(ratioDown) %in% names(weightDown)))
  stopifnot(all(names(ratioDown) %in% names(weightUp)))
  whichNaDown <- names(weightDown)[is.na(weightDown)]
  whichNaUp <- names(weightUp)[is.na(weightUp)]
  stopifnot(all(mergedDT[mergedDT$region %in% whichNaUp,"rescFC"] <= 0))
  stopifnot(all(mergedDT[mergedDT$region %in% whichNaDown,"rescFC"] >= 0))
  # might no be true because ratioDown takes into account when FC = 0
  # stopifnot(all(ratioDown[which(names(ratioDown) %in% whichNaDown)] == 0))
  # stopifnot(all(ratioDown[which(names(ratioDown) %in% whichNaUp)] == 1))
  final_ratio <- unlist(sapply(names(ratioDown), function(x) {
    if(ratioDown[x] > 0.5) {
      ratioWeight <- 0.5+abs(weightDown[x])*0.5
    }else if(ratioDown[x] < 0.5){
      ratioWeight <- 0.5-weightUp[x]*0.5
    }else{
      ratioWeight <- 0.5
    }
    return(ratioWeight)
  }))
  names(final_ratio) <- names(ratioDown)
  stopifnot(all(final_ratio >= 0 & final_ratio <= 1))
  return(final_ratio)
}

get_rescWeightedQQByRegion_v2 <- function(g2TADdt, DEdt) {
  suppressPackageStartupMessages(library(dplyr, warn.conflicts = FALSE, quietly = TRUE, verbose = FALSE))
  g2TADdt$entrezID <- as.character(g2TADdt$entrezID)
  g2TADdt$region <- as.character(g2TADdt$region)
  DEdt$genes <- as.character(DEdt$genes)
  g2TADdt <- g2TADdt[g2TADdt$entrezID %in% DEdt$genes,]
  mergedDT <- left_join(DEdt[,c("genes", "logFC")], g2TADdt[,c("entrezID", "region")], by=c("genes" = "entrezID"))
  mergedDT <- na.omit(mergedDT)
  stopifnot(is.numeric(mergedDT$logFC[1]))
  # rescaled logFC genome-wide
  mergedDT$rescFC <- sign(mergedDT$logFC) * rescFunc(abs(quantNorm(mergedDT$logFC)), newMin = 0, newMax = 1)
  weightDownDT <- aggregate(rescFC ~ region, data=mergedDT, function(x) mean(x[x<0], na.rm=T))
  weightDown <- setNames(weightDownDT$rescFC, weightDownDT$region)
  stopifnot(all(na.omit(weightDown) < 0))
  weightUpDT <- aggregate(rescFC ~ region, data=mergedDT, function(x) mean(x[x>0], na.rm=T))
  weightUp <- setNames(weightUpDT$rescFC, weightUpDT$region)
  stopifnot(all(na.omit(weightUp) > 0))
  ratioDownDT <- aggregate(rescFC ~ region, data=mergedDT, function(x) sum(x < 0)/length(x[x!=0]))
  ratioDown <- setNames(ratioDownDT$rescFC, ratioDownDT$region)
  # ratioDown <- get_ratioDownByRegion_v2(g2TADdt,DEdt)
  # there is a problem with this function that does not exclude the ==0,
  # e.g. if one is down and 2 other are 0
  # => the ratioDown will be < 0.5 but weightUp will be NaN
  stopifnot(!any(is.na(ratioDown)))
  stopifnot(all(names(ratioDown) %in% names(weightDown)))
  stopifnot(all(names(ratioDown) %in% names(weightUp)))
  whichNaDown <- names(weightDown)[is.na(weightDown)]
  whichNaUp <- names(weightUp)[is.na(weightUp)]
  stopifnot(all(mergedDT[mergedDT$region %in% whichNaUp,"rescFC"] <= 0))
  stopifnot(all(mergedDT[mergedDT$region %in% whichNaDown,"rescFC"] >= 0))
  # might no be true because ratioDown takes into account when FC = 0
  # stopifnot(all(ratioDown[which(names(ratioDown) %in% whichNaDown)] == 0))
  # stopifnot(all(ratioDown[which(names(ratioDown) %in% whichNaUp)] == 1))
  final_ratio <- unlist(sapply(names(ratioDown), function(x) {
    if(ratioDown[x] > 0.5) {
      ratioWeight <- 0.5+abs(weightDown[x])*0.5
    }else if(ratioDown[x] < 0.5){
      ratioWeight <- 0.5-weightUp[x]*0.5
    }else{
      ratioWeight <- 0.5
    }
    return(ratioWeight)
  }))
  names(final_ratio) <- names(ratioDown)
  stopifnot(all(final_ratio >= 0 & final_ratio <= 1))
  return(final_ratio)
}

