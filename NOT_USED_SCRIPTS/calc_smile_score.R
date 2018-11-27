

calc_smile_score <- function(currOutFold, mystat = "ratioDown") {

    #****************************************************************************************************************************
    # INPUT DATA
    gene2tadDT <- read.delim(gene2tadDT_file, header=F, col.names = c("entrezID", "chromo", "start", "end", "region"), stringsAsFactors = F)
    gene2tadDT$entrezID <- as.character(gene2tadDT$entrezID)

    # UPDATE SELECT THE GENES ACCORDING TO THE SETTINGS PREPARED IN 0_PREPGENEDATA
    initList <- eval(parse(text = load(paste0(currOutFold, "/", script0_name, "/rna_geneList.Rdata"))))
    geneList <- eval(parse(text = load(paste0(currOutFold, "/", script0_name, "/pipeline_geneList.Rdata"))))
    txt <- paste0(toupper(script_name), "> Start with # genes: ", length(geneList), "/", length(initList), "\n")
    cat(txt)

    stopifnot(!any(duplicated(names(geneList))))

    gene2tadDT <- gene2tadDT[gene2tadDT$entrezID %in% geneList,]

    #****************************************************** IMPORTANT SETTINGS HERE !!!
    histBreakStep <- 0.1
    stopifnot( (1/histBreakStep) %% 1 == 0 )
    #***********************************************************************************

    ###################################################################################################

    obs_statVect <- eval(parse(text = load(paste0(currOutFold, "/", script8_name, "/all_obs_", mystat, ".Rdata"))))
    permut_stat_DT <- eval(parse(text = load(paste0(currOutFold, "/", script8_name, "/", mystat, "_permDT.Rdata"))))

    ### !!! NEED TO CHECK WHY THERE ARE SOME NA
    permut_stat_DT <- na.omit(permut_stat_DT)

    nRandom <- nRandomPermut

    intersectRegions <- intersect(names(obs_statVect), rownames(permut_stat_DT))
    gene2tadDT <- gene2tadDT[gene2tadDT$region %in% intersectRegions,]
    # some regions could have been discarded because no gene
    intersectRegions <- intersect(intersectRegions, gene2tadDT$region)

    # take only the TADs

    ### filter TAD only regions
    if(useTADonly) {
        if( length(grep("_TAD", rownames(permut_stat_DT))) > 0 ) {
            txt <- paste0(toupper(script_name), "> !!! WARNING: permutation ", mystat, " data contain non-TAD regions as well !!!\n")
            cat(txt)
        }
        initLen <- length(intersectRegions)
        intersectRegions <- intersectRegions[grep("_TAD", intersectRegions)]
        txt <- paste0(toupper(script_name), "> Take only TAD regions: ", length(intersectRegions), "/", initLen, "\n")
        cat(txt)
    }

    gene2tadDT <- gene2tadDT[gene2tadDT$region %in% intersectRegions,]
    nbrGenes <- setNames(as.numeric(table(gene2tadDT$region)), names(table(gene2tadDT$region)))

    initNrow <- nrow(permut_stat_DT)
    permut_stat_DT <- permut_stat_DT[intersectRegions,]
    txt <- paste0(toupper(script_name), "> Take only the regions form permut for which genes were retained: ", nrow(permut_stat_DT), "/", initNrow, "\n")
    cat(txt)

    initLen <- length(obs_statVect)
    obs_statVect <- obs_statVect[intersectRegions]
    txt <- paste0(toupper(script_name), "> Take only the regions form observ. for which genes were retained: ", length(obs_statVect), "/", initLen, "\n")
    cat(txt)

    stopifnot(nrow(permut_stat_DT) == length(obs_statVect))

    ###################################################################################################

    cat("... Warning - considering only TAD regions: \n")

    cat("> Prepare simulation data \n")
    nTAD <- sum(regexpr("_TAD",  rownames(permut_stat_DT)) > 0)
    cat(paste0("... TADs: ", nTAD, "/", nrow(permut_stat_DT), "\n"))
    permut_stat_vect <- as.vector(permut_stat_DT[grep("_TAD", rownames(permut_stat_DT)),])
    stopifnot(length(permut_stat_vect) == ncol(permut_stat_DT) * nTAD)

    cat(paste0("length(permut_stat_vect) = ", length(permut_stat_vect), "\n"))

    cat("> Prepare observed data \n")
    nTAD <- sum(regexpr("_TAD",  names(obs_statVect)) > 0)
    cat(paste0("... TADs: ", nTAD, "/", length(obs_statVect), "\n"))
    obs_statVect <- obs_statVect[grep("_TAD", names(obs_statVect))]

    cat("> Prepare data for plotting \n")

    obs_rel_hist <- hist(obs_statVect, breaks=seq(0,1,histBreakStep), plot=F)
    #plot(obs_rel_hist)

    cat(paste0("length(permut_stat_vect) = ", length(permut_stat_vect), "\n"))
    cat(paste0("length(permut_stat_vect) = ", length(na.omit(permut_stat_vect)), "\n"))
    shuff_rel_hist <-  hist(permut_stat_vect,breaks=seq(0,1,histBreakStep), plot=F)
    # plot(shuff_rel_hist)

    cat(paste0("ncol(permut_stat_DT) = ", ncol(permut_stat_DT), "\n"))
    cat(paste0("sum(shuff_rel_hist$counts) = ", sum(shuff_rel_hist$counts), "\n"))
    cat(paste0("sum(shuff_rel_hist$counts/nRandom) = ", sum(shuff_rel_hist$counts/nRandom), "\n"))
    cat(paste0("sum(obs_rel_hist$counts) = ", sum(obs_rel_hist$counts),"\n"))
    stopifnot(sum(shuff_rel_hist$counts/nRandom) == sum(obs_rel_hist$counts))   # 830
    stopifnot(sum(shuff_rel_hist$counts) == length(permut_stat_vect))
    stopifnot(sum(obs_rel_hist$counts) == length(na.omit(obs_statVect)))

    rel_freqValues <- rbind(obs_rel_hist$counts, shuff_rel_hist$counts/nRandom)
    rownames(rel_freqValues) <- c("observed", "randomized")

    # Calculate observed/expected ratio
    rel_logRatio <- log2(rel_freqValues["observed",]/rel_freqValues["randomized",])
    # put y coord at the middle of the breaks
    rel_ycoord <- (seq(0,1,histBreakStep) * 100)[-1]
    rel_ycoord <- rel_ycoord-(histBreakStep*100)/2
    stopifnot(length(rel_ycoord) == length(rel_logRatio))
    toKeep <- !is.na(rel_logRatio) & abs(rel_logRatio) != "Inf"
    rel_logRatio <- rel_logRatio[toKeep]
    rel_ycoord <- rel_ycoord[toKeep]

    smile_score <- get_smileScore(ylogRatioVect=rel_logRatio, xBreakVect=rel_ycoord, withPlot=FALSE)

    cat(paste0("> RETURN smile_score = ", round(smile_score,2), "\n"))

    return(smile_score)

}
